/**
 * Scans a result file for admixture. Admixture is defined as members of ingroup
 * being in a clade with members of admixer that does not include members of
 * outgroup.
 */
#include <getopt.h>
#include <argp.h>
#include <zlib.h>
#include <math.h>
#include <time.h>
#include <string>
#include <algorithm>
#include <vector>
#include <deque>
#include <set>
#include <iterator>
#include <map>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <zlib.h>
#include <set>
#include "sets.h"
#include "treeNode.h"
#include "serialize.h"
#include "common.h"

using std::cout;
using std::endl;
using namespace std;

/**
 * Print a help message to the terminal and exit.
 */
void help(int code){
   fprintf(stderr, "admix_scan [OPTIONS]\n");
   fprintf(stderr, "Scans a SARGE output file for admixture and outputs statistics about it.\
Admixture is defined as a clade containing members of ingroup and admixer, but \
not outgroup.\n");
   fprintf(stderr, "[OPTIONS]:\n");
   fprintf(stderr, "    --ingroup -i One or more populations that are admixed (can \
specify more than once)\n");
   fprintf(stderr, "    --outgroup -o One or more populations that are not admixed (can \
specify more than once)\n");
   fprintf(stderr, "    --admixer -a The admixing population of interest (optionally specify multiple times)\n");
fprintf(stderr, "    --indvs -v (OPTIONAL) The file containing names of individuals \n");
   fprintf(stderr, "    --ploidy -n (OPTIONAL) The number of haplotypes per individual (\
default 2) \n");
    fprintf(stderr, "    --require_mutations -M (OPTIONAL) if set, a clade will only count \
as indicative of admixture if it is supported by actual mutations (and not just inferred \
by ancestral recombination events)\n");
    fprintf(stderr, "    --pops -p A file mapping individual names to population names, \
tab-separated\n");
    fprintf(stderr, "    --map_prefix -m (OPTIONAL) A prefix to create output files of the format \
<PREFIX>.hapname.bed. Each will be a bed file of regions where the given haplotype is \
admixed.\n");
    fprintf(stderr, "    --sites_prefix -s (OPTIONAL) A prefix to create output files of the format \
<PREFIX>.hapname.sites. Each line in each file will consist of chromosome and position \
(tab separated) at a site that informed a clade indicating the given haplotype was admixed.\n");
    fprintf(stderr, "    --bufsize -b The number of characters to read from the \
input file at a time\n");
    fprintf(stderr, "    --polytomy -P the maximum size of a polytomy to allow as an admixed clade\n");
exit(code);
}

// Define a structure to store summary data about admixture within a population.
struct pop_stats {
    string popname;
    float perc_admixed;
    float perc_admixed_sites;
    float tract_len;
    long int tract_len_min;
    long int tract_len_max;
    float allele_freq;
};

bool pop_stats_comp(pop_stats s1, pop_stats s2){
    if (s1.perc_admixed == s2.perc_admixed){
        return (s1.popname.compare(s2.popname) > 0);
    }
    else{
        return s1.perc_admixed > s2.perc_admixed;
    }
}


float tmrca(treeNode* node1, treeNode* node2, treeNode* root, int& steps){
    long int rootlen = root->mutations.size();
    if (rootlen == 0){
        rootlen = 1;
    }
    steps = 0;
    treeNode* parent = NULL;
    treeNode* child = NULL;
    float dist_to_parent = 0.0;
    if (issuperset_bitset(node1->subtree_leaves(), node2->subtree_leaves())){
        parent = node1;
        child = node2;
    }
    else if (issuperset_bitset(node2->subtree_leaves(), node1->subtree_leaves())){
        parent = node2;
        child = node1;
    }
    else{
        parent = node1;
        child = node2;
        while(parent != NULL && !issuperset_bitset(parent->subtree_leaves(), child->subtree_leaves())){
            dist_to_parent += parent->dist_norm;
            steps += 1;
            //dist_to_parent += ((float)parent->mutations.size()/(float)rootlen);
            //dist_to_parent += parent->dist;
            parent = parent->parent;
        }
    }
    treeNode* origparent = parent;
    treeNode* origchild = child;
    float dist_to_child = 0.0;
    while(child != NULL && child != parent){
        dist_to_child += child->dist_norm;
        steps += 1;
        //dist_to_child += ((float)child->mutations.size()/(float)rootlen);
        //dist_to_child += child->dist;
        child = child->parent;
    }
    return (dist_to_child + dist_to_parent)/2.0;
}

float dist_down_aux(treeNode* parent, cladeset& hapset, float dist){
    for (vector<treeNode*>::iterator child = parent->children.begin();
        child != parent->children.end(); ++child){
        cladeset st = (*child)->subtree_leaves();
        if ((st & hapset).count() > 0){
            float dist_new = dist_down_aux(parent, hapset, dist + (*child)->dist);
            dist = dist_new;
        }
    }
    return dist;
}

float dist_down(treeNode* parent, unsigned int hap, int num_haplotypes){
    cladeset hap_set;
    hap_set.set(num_haplotypes-1-hap);
    return dist_down_aux(parent, hap_set, 0.0);
}

struct admixinfo{
    bool admixed;
    // Lowest TMRCA to an admixing haplotype
    float tmrca_admixer;
    // Lowest TMRCA to an outgroup haplotype
    float tmrca_outgroup;
    int num_ingroup;
    float dist_below;
    admixinfo(bool a, float t1, float t2, int ni, float db){
        this->admixed = a;
        this->tmrca_admixer = t1;
        this->tmrca_outgroup = t2;
        this->num_ingroup = ni;
        this->dist_below = db;
    };
};

bool admixed_clades(treeNode* n, 
    cladeset& admixers_all,
    vector<cladeset>& admixers_single, 
    cladeset& ingroup_all, 
    cladeset& outgroup_all,
    unordered_set<cladeset>& clades,
    bool require_mutations,
    int polytomy_size){
    
    float outgroup_thresh = 0.1;
    
    for (vector<cladeset>::iterator adm = admixers_single.begin(); adm != admixers_single.end();
        ++adm){
        treeNode* adm_n = n->get_smallest_containing(*adm);
        bool stop = false;
        treeNode* parent = adm_n;
        while(!stop){
            cladeset st = parent->subtree_leaves();
            cladeset st_ig = (st & ingroup_all);
            if (st_ig.count() > 0 && st_ig.count() < ingroup_all.count()){
                stop = true;
                if (!require_mutations || parent->mutations.size() > 0){
                    if (outgroup_all.count() == 0 || 
                        ((float)(st & outgroup_all).count()/(float)outgroup_all.count()) < outgroup_thresh){
                        if (polytomy_size == -1 || st.count() <= polytomy_size){
                            clades.insert(st);
                        }
                    }
                }
            }
            else{
                if (parent->parent == NULL){
                    stop = true;
                }
                else{
                    parent = parent->parent;
                }
            }
        }
    }
}

void end_runs(string& admixer_pop,
    set<unsigned int>& runs_end, 
    map<unsigned int, map<long int, admixinfo> >& admix_map, 
    string& prevchrom,
    long int cursite,
    vector<string>& indvs,
    map<string, string>& indv2pop,
    bool write_bed,
    map<string, FILE*>& bed_out,
    bool eof,
    map<string, float>& pop_sums,
    map<string, int>& pop_counts,
    map<string, float>& pop_perc_sums,
    map<string, int>& pop_perc_counts){

    long int propdist = 25000;
    
    for (set<unsigned int>::iterator hap = runs_end.begin(); hap != runs_end.end();){
        
        // Deal with any run of admixture that exists
        if (admix_map.count(*hap) > 0 && admix_map[*hap].size() > 1){
            long int prevsite = -1;
            long int runstart = -1;
            string popkey = indv2pop[indvs[*hap]];
            long int adm_sum = 0;
            
            float tmrca_sum = 0;
            int tmrca_count = 0;
            float tmrca_outside_sum = 0;
            float tmrca_min = -1;
            float tmrca_outside_min = -1;
            
            float dist_below_sum = 0;
            float dist_below_min = -1;
            
            float size_sum = 0;
            bool prev_admixed = false;
            
            for (map<long int, admixinfo>::iterator m = admix_map[*hap].begin();
                m != admix_map[*hap].end();){
                if (m->second.admixed && !prev_admixed){
                    if (prevsite != -1){
                        runstart = prevsite + (long int)round((float)(m->first-prevsite)/2.0);
                        if (runstart < m->first - propdist){
                            runstart = m->first - propdist;
                        }
                    }
                    else{
                        runstart = m->first;
                    }
                    tmrca_sum = m->second.tmrca_admixer;
                    tmrca_count = 1;
                    tmrca_outside_sum = m->second.tmrca_outgroup;
                    tmrca_min = m->second.tmrca_admixer;
                    tmrca_outside_min = m->second.tmrca_outgroup;
                    size_sum = m->second.num_ingroup;
                    dist_below_sum = m->second.dist_below;
                    dist_below_min = m->second.dist_below;
                }
                else if (m->second.admixed && prev_admixed){
                    tmrca_sum += m->second.tmrca_admixer;
                    tmrca_outside_sum += m->second.tmrca_outgroup;
                    tmrca_count++;
                    size_sum += m->second.num_ingroup;
                    dist_below_sum += m->second.dist_below;
                    if (tmrca_min == -1 || m->second.tmrca_admixer < tmrca_min){
                        tmrca_min = m->second.tmrca_admixer;
                    }
                    if (tmrca_outside_min == -1 || m->second.tmrca_outgroup < tmrca_outside_min){
                        tmrca_outside_min = m->second.tmrca_outgroup;
                    }
                    if (dist_below_min == -1 || m->second.dist_below < dist_below_min){
                        dist_below_min = m->second.dist_below;
                    }
                }
                else if (!m->second.admixed && prev_admixed){
                    long int runend;
                    if (prevsite == -1){
                        runend = m->first;
                    }
                    else{
                        runend = prevsite + (long int)round((float)(m->first-prevsite)/2.0);
                    }
                    if (tmrca_count == 0){
                        tmrca_count = 1;
                    }
                    float tmrca_mean = tmrca_sum / (float)tmrca_count;
                    float tmrca_outside_mean = tmrca_outside_sum / (float)tmrca_count;
                    float dist_below_mean = dist_below_sum / (float)tmrca_count;
                    if (write_bed){
                        FILE* outf = bed_out[indv2pop[indvs[*hap]]];
                        fprintf(outf, "%s\t%ld\t%ld\t%ld\t%f\t%f\t%f\t%f\t%s\t%f\t%f\t%f\t%s\n", prevchrom.c_str(), 
                            runstart-1, runend, runend-runstart+1, tmrca_min,
                            tmrca_mean, tmrca_outside_min, tmrca_outside_mean, 
                            indvs[*hap].c_str(), size_sum / (float)tmrca_count,
                            dist_below_mean, dist_below_min, admixer_pop.c_str());
                    }
                    tmrca_sum = 0;
                    tmrca_count = 0;
                    size_sum = 0;
                    adm_sum += (runend-runstart+1);
                    if (pop_sums.count(popkey) == 0){
                        pop_sums.insert(make_pair(popkey, 0.0));
                        pop_counts.insert(make_pair(popkey, 0));
                    }
                    pop_sums[popkey] += (runend-runstart+1);
                    pop_counts[popkey]++;
                }
        
                prev_admixed = m->second.admixed;
                prevsite = m->first;
                admix_map[*hap].erase(m++);
            }
            
            if (eof && prev_admixed && runstart != -1){
                if (tmrca_count == 0){
                    tmrca_count = 1;
                }
                if (write_bed){
                    float tmrca_mean = tmrca_sum / (float)tmrca_count;
                    float tmrca_outside_mean = tmrca_outside_sum / (float)tmrca_count;
                    float dist_below_mean = dist_below_sum / (float)tmrca_count;
                    FILE* outf = bed_out[indv2pop[indvs[*hap]]];
                    fprintf(outf, "%s\t%ld\t%ld\t%ld\t%f\t%f\t%f\t%f\t%s\t%f\t%f\t%f\t%s\n", prevchrom.c_str(), 
                        runstart-1, cursite, cursite-runstart+1, tmrca_min, tmrca_mean,
                        tmrca_outside_min, tmrca_outside_mean, indvs[*hap].c_str(),
                        size_sum / (float)tmrca_count, dist_below_mean, dist_below_min,
                        admixer_pop.c_str());
                }
                tmrca_sum = 0;
                tmrca_count = 0;
                adm_sum += (cursite-runstart+1);
                if (pop_sums.count(popkey) == 0){
                    pop_sums.insert(make_pair(popkey, 0.0));
                    pop_counts.insert(make_pair(popkey, 0));
                }
                pop_sums[popkey] += (cursite-runstart+1);
                pop_counts[popkey]++;
            }
            if (eof){
                fprintf(stdout, "%s\t%f\n", indvs[*hap].c_str(), (float)adm_sum/(float)cursite);
            }
            if (pop_perc_sums.count(popkey) == 0){
                pop_perc_sums.insert(make_pair(popkey, 0.0));
                pop_perc_counts.insert(make_pair(popkey, 0));
            }
            pop_perc_sums[popkey] += (float)adm_sum/(float)cursite;
            pop_perc_counts[popkey] += 1;
        }
        
        if (!eof){
            // Mark the current site un-admixed.
            admix_map[*hap].insert(make_pair(cursite, admixinfo(false, -1, -1, -1, -1)));
        }
        runs_end.erase(hap++);
    }
}

int main(int argc, char *argv[]) {    
    // Define arguments 
    // http://www.gnu.org/software/libc/manual/html_node/Getopt-Long-Options.html#Getopt-Long-Options
    // http://www.gnu.org/software/libc/manual/html_node/Getopt-Long-Option-Example.html#Getopt-Long-Option-Example
    // Fields for each argument: name, has_arg (values: no_argument, required_argument,
    //     optional_argument)
    // flag = int value to store flag for the option, or NULL if option is string
    // val = short name for string option, or NULL
    static struct option long_options[] = {
       {"ingroup", required_argument, 0, 'i'},
       {"outgroup", required_argument, 0, 'o'},
       {"admixer", required_argument, 0, 'a'},
       {"indvs", required_argument, 0, 'v'},
       {"ploidy", optional_argument, 0, 'n'},
       {"pops", required_argument, 0, 'p'},
       {"bufsize", optional_argument, 0, 'b'},
       {"map_prefix", required_argument, 0, 'm'},
       {"require_mutations", no_argument, 0, 'M'},
       {"sites_prefix", required_argument, 0, 's'},
       {"polytomy", required_argument, 0, 'P'},
       {"help", optional_argument, 0, 'h'},
       {0, 0, 0, 0} 
    };
    
    // Set default values
    int bufsize = 1048576;
    int ploidy = 2;
    
    vector<string> ingroup_pops;
    vector<string> outgroup_pops;
    vector<string> admixers;
    string indvfilename;
    string popfilename;
    string map_prefix;
    string sites_prefix;
    bool require_mutations = false;
    int polytomy_size = -1;
    
    int option_index = 0;
    int ch;
    
    if (argc == 1){
        help(0);
    }
    while((ch = getopt_long(argc, argv, "i:o:a:v:n:p:b:P:m:s:Mh", long_options, &option_index )) != -1){
        switch(ch){
            case 0:
                // This option set a flag. No need to do anything here.
                break;
            case 'i':
                ingroup_pops.push_back(string(optarg));
                break;
            case 'o':
                outgroup_pops.push_back(string(optarg));
                break;
            case 'a':
                admixers.push_back(string(optarg));
                break;
            case 'n':
                ploidy = atoi(optarg);
                break;
            case 'v':
                indvfilename = string(optarg);
                break;
            case 'p':
                popfilename = string(optarg);
                break;
            case 'b':
                bufsize = atoi(optarg);
                break;
            case 'P':
                polytomy_size = atoi(optarg);
                break;
            case 'm':
                map_prefix = optarg;
                break;
            case 's':
                sites_prefix = optarg;
                break;
            case 'M':
                require_mutations = true;
                break;
            case '?':
                //help(0);
                break;
            case 'h':
                help(0);
                break;
            default:
                help(0);
        }    
    }
    
    if (ingroup_pops.size() == 0 || admixers.size() == 0){
        fprintf(stderr, "ERROR: you must provide at least one ingroup and at least one \
admixer.\n");
        exit(1);
    }
    
    // Get file reading stuff ready
    FILE *instream = stdin;
    if (instream == NULL){
        fprintf(stderr, "Error opening input stream\n");
        exit(1);
    }
    gzFile fp = gzdopen(fileno(instream), "rb");
    if (!fp){
        fprintf(stderr, "ERROR: unable to read from stdin.\n");
        exit(1);
    }
    bool write_bed = false;
    if (map_prefix.length() > 0){
        write_bed = true;
        if (indvfilename.length() == 0){
            fprintf(stderr, "ERROR: you must provide a file of haplotype names to write \
out ancestry maps in BED format.\n");
            exit(1);
        }
    }
    bool write_sites = false;
    if (sites_prefix.length() > 0){
        write_sites = true;
        if (indvfilename.length() == 0){
            fprintf(stderr, "ERROR: you must provide a file of haplotype names to write \
out files of admixed sites.\n");
        }
    }
    instream_info is;
    instream_init(is, &fp, bufsize);
    
    // Read file header.
    int num_haplotypes;
    cladeset mask;
    read_header(is, num_haplotypes, mask);
    
    vector<string> indvs;
    parse_indvs(indvs, indvfilename);
    
    map<string, int> indv2hap;
    parse_indvs_map(indv2hap, indvfilename);
    map<string, string> indv2pop;
    map<string, vector<string> > pop2indv;
    parse_pops(indv2pop, pop2indv, popfilename, ploidy);
    
    vector<string> leafname_pops;
    for (vector<string>::iterator indvit = indvs.begin(); indvit != indvs.end();
        ++indvit){
        //string indv_sub = (*indvit).substr(0, (*indvit).find('-'));
        leafname_pops.push_back(indv2pop[*indvit]);
    }
    
    // Build ingroup and outgroup clades
    
    set<unsigned int> ingroup_set;
    for (vector<string>::iterator ingroup_pop = ingroup_pops.begin();
        ingroup_pop != ingroup_pops.end(); ++ingroup_pop){
        for (vector<string>::iterator ingroup_name = pop2indv[*ingroup_pop].begin();
            ingroup_name != pop2indv[*ingroup_pop].end(); ++ingroup_name){
            if (indv2hap.count(*ingroup_name) > 0){
                ingroup_set.insert(indv2hap[*ingroup_name]);
                //fprintf(stderr, "ingroup %d %s\n", indv2hap[*ingroup_name], ingroup_name->c_str());
            }
        }
    }
    
    // Create cladesets for each individual ingroup haplotype
    vector<cladeset> ingroups_split;
    for (set<unsigned int>::iterator ig = ingroup_set.begin(); ig != ingroup_set.end();
        ++ig){
        cladeset igcs;
        igcs.set(num_haplotypes-1-*ig);
        ingroups_split.push_back(igcs);
    }
    
    cladeset ingroup_clade = set2bitset(ingroup_set, num_haplotypes);
    
    set<unsigned int> outgroup_set;
    for (vector<string>::iterator outgroup_pop = outgroup_pops.begin();
        outgroup_pop != outgroup_pops.end(); ++outgroup_pop){
        for (vector<string>::iterator outgroup_name = pop2indv[*outgroup_pop].begin();
            outgroup_name != pop2indv[*outgroup_pop].end(); ++outgroup_name){
            if (indv2hap.count(*outgroup_name) > 0){
                outgroup_set.insert(indv2hap[*outgroup_name]);
                //fprintf(stderr, "outgroup %d %s\n", indv2hap[*outgroup_name], outgroup_name->c_str());
            }
        }
    }

    cladeset outgroup_clade = set2bitset(outgroup_set, num_haplotypes);
    
    vector<cladeset > admix_clades;
    cladeset admixers_all;
    unordered_map<cladeset, string> admixer2pop;
    map<string, int> admixer2ind;
    for (int i = 0; i < admixers.size(); ++i){
        admixer2ind.insert(make_pair(admixers[i], i));
    }
    map<int, cladeset> pops_admclades;
    for (int i = 0; i < admixers.size(); ++i){
        cladeset cl;
        pops_admclades.insert(make_pair(i, cl));
    }
    int popind = 0;
    for (vector<string>::iterator admixer = admixers.begin(); admixer != admixers.end(); ++admixer){
        for (vector<string>::iterator admix_hap = pop2indv[*admixer].begin();
            admix_hap != pop2indv[*admixer].end(); ++admix_hap){
            if (indv2hap.count(*admix_hap) > 0){
                set<unsigned int> admix_hap_set;
                admix_hap_set.insert(indv2hap[*admix_hap]);
                cladeset admix_clade = set2bitset(admix_hap_set, num_haplotypes);
                admixers_all |= admix_clade;
                admix_clades.push_back(admix_clade);
                admixer2pop.insert(make_pair(admix_clade, *admixer));
                pops_admclades[popind] |= admix_clade;
                //fprintf(stderr, "admix %d %s\n", indv2hap[*admix_hap], admix_hap->c_str());
            }
        }
        ++popind;
    }
    if (admixers.size() > 1){
        // Add "unknown" population
        admixers.push_back("UNK");
        admixer2ind.insert(make_pair("UNK", admixers.size()-1));
    }
    
        // Create file handles for writing BED and sites files
    vector<map<string, FILE*> > ig_pops_bed_out;
    vector<map<string, FILE*> > ig_pops_sites_out;
    
    if (write_bed || write_sites){
        for (int adm_ind = 0; adm_ind < admixers.size(); ++adm_ind){
            map<string, FILE*> ig_pops_bed_row;
            map<string, FILE*> ig_pops_sites_row;
            for(vector<string>::iterator ig = ingroup_pops.begin(); ig != ingroup_pops.end();
                ++ig){
                // Create output files.
                if (write_bed){
                    string outname = map_prefix + "_" + admixers[adm_ind] + "_" + *ig + ".bed";
                    FILE* outf_bed = fopen(outname.c_str(), "w");
                    if (!outf_bed){
                        fprintf(stderr, "ERROR opening file %s for writing.\n", outname.c_str());
                        exit(1);
                    }
                    ig_pops_bed_row.insert(make_pair(*ig, outf_bed));
                }
                if (write_sites){
                    string outname = sites_prefix + "_" + admixers[adm_ind] + "_" + *ig + ".sites";
                    FILE* outf_sites = fopen(outname.c_str(), "w");
                    if (!outf_sites){
                        fprintf(stderr, "ERROR opening file %s for writing.\n", outname.c_str());
                        exit(1);
                    }
                    ig_pops_sites_row.insert(make_pair(*ig, outf_sites));
                }
            }
            ig_pops_bed_out.push_back(ig_pops_bed_row);
            ig_pops_sites_out.push_back(ig_pops_sites_row);
        }
    }
    
    fprintf(stderr, "Ingroup haps: %ld Outgroup haps: %ld Admixing haps: %ld\n", 
        ingroup_set.size(), outgroup_set.size(), admix_clades.size());
    
    /*
    if (ingroup_set.size() == 0 || outgroup_set.size() == 0 || admix_clades.size() == 0){
        fprintf(stderr, "ERROR: the haplotypes you have subsampled do not include any individuals \
from either the ingroup, outgroup or admixing group.\n");
        exit(1);
    }
    */
    
    // Prep data structures for storing stats
    // For each individual, we want to store the sum of Neanderthal bases,
    // the number of runs of Neanderthal bases, and the sum of lengths of runs.
    
    // For each haplotype, we're either currently NOT in a run of Neanderthal
    // ancestry (-1), or we'll store where our current run of Neanderthal ancestry
    // started (long int)
    map<int, long int> cur_run_start;
    
    long int prevpos;
    long int firstpos = -1;
    
    long int progress = 0;
    
    string prevchrom;
    
    map<unsigned int, long int> hap_bases;
    
    vector<map<unsigned int, map<long int, admixinfo> > > admix_map;
    for (int i = 0; i < admixers.size(); ++i){
        map<unsigned int, map<long int, admixinfo> > m;
        admix_map.push_back(m);
    }
    
    unordered_set<cladeset> clades_adm_prev;
    
    vector<set<unsigned int> > runs_end;
    for (int i = 0; i < admixers.size(); ++i){
        set<unsigned int> s;
        runs_end.push_back(s);
    }
    
    map<string, float> pop_sums;
    map<string, int> pop_counts;
    
    map<string, float> pop_perc_sums;
    map<string, int> pop_perc_counts;
    
    // Read from stdin.
    while(!is.finished()){

        string chrom;
        long int pos;
        
        treeNode tree;

        read_sitedata(num_haplotypes, is, chrom, pos, tree);

        // Remove all ancestral stuff
        //tree.leaves &= mask;
        
        if (pos - progress > 100000){
            fprintf(stderr, "Read %s %ld\r", chrom.c_str(), pos);
            progress = pos;
        }
        if (firstpos == -1){
            firstpos = pos;
        }
        
        // Find all admixer haplotypes
        set<treeNode*> admixers_nodes;
        map<treeNode*, string> admixer_node_pops;
        for (vector<cladeset>::iterator ac = admix_clades.begin(); ac != admix_clades.end();
            ++ac){
            treeNode* an = tree.get_smallest_containing(*ac);
            admixers_nodes.insert(an);
            admixer_node_pops.insert(make_pair(an, admixer2pop[*ac]));
        }
        
        if (true){
            unordered_set<cladeset> clades_adm;
            admixed_clades(&tree, admixers_all, admix_clades, ingroup_clade, outgroup_clade, 
                clades_adm, require_mutations, polytomy_size);
            /*
            for (unordered_set<cladeset>::iterator cl = clades_adm_prev.begin(); cl != clades_adm_prev.end();){
                bool still_admixed = true;
                if (clades_adm.find(*cl) == clades_adm.end()){
                    cladeset clcpy = *cl;
                    //if (tree.get_clade_match(clcpy) == NULL){
                    if (tree.has_clade_invalidates(clcpy)){
                        // Not admixed.
                        still_admixed = false;
                    }
                    else{
                        
                        treeNode* n = tree.get_smallest_containing(clcpy);
                        // If this is closer to other ingroup haplotypes outside the set
                        // than to admixers, then it's not admixed.
                        bool stop = n->parent == NULL;
                        cladeset ig_other = set_diff_bitset(ingroup_clade, clcpy);
                        while (!stop){
                            n = n->parent;
                            cladeset st = n->subtree_leaves();
                            st = set_diff_bitset(st, clcpy);
                            if ((st & ig_other).count() > 0 && (st & admixers_all).count() == 0){
                                still_admixed = false;
                                stop = true;
                            }
                            else if (n->parent == NULL){
                                stop = true;
                            }
                        }
                        
                    }
                }
                if (still_admixed){
                    clades_adm.insert(*cl);
                }
                clades_adm_prev.erase(cl++);
            }
            */
            // Make one of these per admixer
            vector<set<unsigned int> > haps_adm;
            for (int i = 0; i < admixers.size(); ++i){
                set<unsigned int> s;
                haps_adm.push_back(s);
            }
            for (unordered_set<cladeset>::iterator cl = clades_adm.begin(); cl != clades_adm.end(); ++cl){
                treeNode* cln = tree.get_smallest_containing(*cl);
                if (cln->subtree_leaves() == *cl){
                    set<unsigned int> haps_adm_clade = bitset2set(*cl, num_haplotypes);
                
                    float min_tmrca_outside = -1;
                    /*
                    treeNode* parent = cln->parent;
                    bool stop_looking = false;
                    while (!stop_looking){
                        cladeset ig_outclade = set_diff_bitset(parent->subtree_leaves(), *cl);
                        ig_outclade = (ig_outclade & ingroup_clade);
                        if (ig_outclade.count() > 0){
                            stop_looking = true;
                            for (vector<treeNode*>::iterator child = parent->children.begin();
                                child != parent->children.end(); ++child){
                                if (*child != cln){
                                    float this_tmrca = tmrca(cln, *child, &tree);
                                    if (min_tmrca_outside == -1 || this_tmrca < min_tmrca_outside){
                                        min_tmrca_outside = this_tmrca;
                                    }
                                }
                            }
                        }
                        if (parent->parent == NULL){
                            stop_looking = true;
                        }
                        if (!stop_looking){
                            parent = parent->parent;
                        }
                    }
                    */
                    
                    for (set<unsigned int>::iterator h = haps_adm_clade.begin(); h != haps_adm_clade.end(); ++h){
                        
                        float min_tmrca = -1;
                        int min_tmrca_steps = -1;
                        set<string> min_tmrca_with;
                        set<string> min_tmrca_steps_with;
                        cladeset hk;
                        hk.set(num_haplotypes-1-*h);
                        treeNode* hn = tree.get_smallest_containing(hk);
                        for (set<treeNode*>::iterator an = admixers_nodes.begin(); an != admixers_nodes.end(); ++an){
                            int steps = 0;
                            float this_tmrca = tmrca(hn, *an, &tree, steps);
                            if (min_tmrca == -1 || this_tmrca < min_tmrca){
                                min_tmrca = this_tmrca;
                                min_tmrca_with.clear();
                                min_tmrca_with.insert(admixer_node_pops[*an]);
                            }
                            else if (min_tmrca == this_tmrca){
                                min_tmrca_with.insert(admixer_node_pops[*an]);
                            }
                            if (min_tmrca_steps == -1 || steps < min_tmrca_steps){
                                min_tmrca_steps = steps;
                                min_tmrca_steps_with.clear();
                                min_tmrca_steps_with.insert(admixer_node_pops[*an]);
                            }
                            else if (min_tmrca_steps == steps){
                                min_tmrca_steps_with.insert(admixer_node_pops[*an]);
                            }
                        }
                    
                        float denom = cln->dist_below + cln->dist_above;
                        if (denom == 0){
                            denom = 1.0;
                        }
                        
                        // Determine what population it belongs to.
                        int pop_ind = -1;
                        if (min_tmrca_steps_with.size() > 1){
                        //if (min_tmrca_with.size() > 1){
                            pop_ind = admixer2ind["UNK"];
                        }
                        else{
                            pop_ind = admixer2ind[*min_tmrca_steps_with.begin()];
                            //pop_ind = admixer2ind[*min_tmrca_with.begin()];
                        }
                        haps_adm[pop_ind].insert(*h);
                        if (admix_map[pop_ind].count(*h) == 0){
                            map<long int, admixinfo> m;
                            admix_map[pop_ind].insert(make_pair(*h, m));
                        }
                        admix_map[pop_ind][*h].insert(make_pair(pos, admixinfo(true, min_tmrca, min_tmrca_outside, cl->count(), cln->dist_below / denom)));
                    
                    }
                }
                else{
                    /*
                    // Clade not invalidated but doesn't exist in current tree. Do nothing
                    // until we see evidence that stuff recombined out of it.
                    set<unsigned int> haps_adm_clade = bitset2set(*cl, num_haplotypes);
                    for (set<unsigned int>::iterator h = haps_adm_clade.begin();
                        h != haps_adm_clade.end(); ++h){
                        haps_adm.insert(*h);
                    }
                    */
                }
            }
            
            for (int pop_ind = 0; pop_ind < admixers.size(); ++pop_ind){
                for (set<unsigned int>::iterator h = ingroup_set.begin(); h != ingroup_set.end(); ++h){
                    if (haps_adm[pop_ind].find(*h) == haps_adm[pop_ind].end()){
                    
                        if (admix_map[pop_ind].count(*h) == 0){
                            map<long int, admixinfo> m;
                            admix_map[pop_ind].insert(make_pair(*h, m));
                        }
                        admix_map[pop_ind][*h].insert(make_pair(pos, admixinfo(false, -1, -1, -1, -1)));
                        runs_end[pop_ind].insert(*h);
                    }
                }
                 end_runs(admixers[pop_ind], runs_end[pop_ind], admix_map[pop_ind], chrom, pos, indvs, 
                    indv2pop, write_bed, ig_pops_bed_out[pop_ind], false, pop_sums, 
                    pop_counts, pop_perc_sums, pop_perc_counts);
            }
            
            //clades_adm_prev = clades_adm;
           
        }

        prevchrom = chrom;
        prevpos = pos;
        
    }
    
    long int len = prevpos;
    
    for (int pop_ind = 0; pop_ind < admixers.size(); ++pop_ind){
        end_runs(admixers[pop_ind], runs_end[pop_ind], admix_map[pop_ind], prevchrom, 
            prevpos, indvs, indv2pop, write_bed,
            ig_pops_bed_out[pop_ind], true, pop_sums, pop_counts, pop_perc_sums, pop_perc_counts);
    }
    fprintf(stderr, "\n");
    
    
    if (write_bed){
        for (vector<map<string, FILE*> >::iterator bo_p = ig_pops_bed_out.begin();
            bo_p != ig_pops_bed_out.end(); ++bo_p){
            for (map<string, FILE*>::iterator bo = bo_p->begin(); bo != bo_p->end(); ++bo){
                if (bo->second != NULL){
                    fclose(bo->second);
                }
            }
        }
    }
    if (write_sites){
        for (vector<map<string, FILE*> >::iterator so_p = ig_pops_sites_out.begin();
            so_p != ig_pops_sites_out.end(); ++so_p){
            for (map<string, FILE*>::iterator so = so_p->begin(); so != so_p->end(); ++so){
                if (so->second != NULL){
                    fclose(so->second);
                }
            }
        }
    }
    return 0;
}
