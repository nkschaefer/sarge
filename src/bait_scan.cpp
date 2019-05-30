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
   fprintf(stderr, "bait_scan [OPTIONS]\n");
   fprintf(stderr, "Scans a SARGE output file for candidate admixed clades, using \
admixer haplotypes as \"bait.\" Whereas admix_scan requires the presence of an \
admixer haplotype in every appearance of an admixed clade, this method uses admixer \
haplotypes to find clades, then follows them upstream and downstream until \
ingroup haplotypes recombine out of them.\n");
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
    fprintf(stderr, "    --mapfile -m The file to which (BED format) admixed tracks \
for individual ingroup haplotypes will be written\n");
    fprintf(stderr, "    --cladefile -c The file to which summary information about each \
admixed clade (regardless of which haplotypes it contains) will be written\n");
    fprintf(stderr, "    --bufsize -b The number of characters to read from the \
input file at a time\n");
    fprintf(stderr, "    --polytomy -P the maximum size of a polytomy to allow as an admixed clade\n");
    fprintf(stderr, "    --dist -d The distance upstream to store trees when determining haplotypes. \
If this number is too big, it will use extra memory and slow down. If too small, complete \
haplotypes will not be recovered. (Default = 75000)\n");
    fprintf(stderr, "    --minsize -s the minimum number of haplotypes for a clade \
to contain to be considered\n");
    fprintf(stderr, "    --maxsize -S the maximum number of haplotypes for a clade \
to contain to be considered\n");
exit(code);
}

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
       {"mapfile", required_argument, 0, 'm'},
       {"cladefile", required_argument, 0, 'c'},
       {"require_mutations", no_argument, 0, 'M'},
       {"polytomy", required_argument, 0, 'P'},
       {"dist", required_argument, 0, 'd'},
       {"minsize", required_argument, 0, 's'},
       {"maxsize", required_argument, 0, 'S'},
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
    string mapfile;
    string cladefile;
    int minsize = 0;
    int maxsize = -1;
    long int dist = 75000;
    
    bool require_mutations = false;
    int polytomy_size = -1;
    
    int option_index = 0;
    int ch;
    
    if (argc == 1){
        help(0);
    }
    while((ch = getopt_long(argc, argv, "i:o:a:v:n:p:b:P:m:c:d:s:S:Mh", long_options, &option_index )) != -1){
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
                mapfile = optarg;
                break;
            case 'c':
                cladefile = optarg;
                break;
            case 'd':
                dist = atol(optarg);
                break;
            case 's':
                minsize = atoi(optarg);
                break;
            case 'S':
                maxsize = atoi(optarg);
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
    
    if (dist < 0){
        fprintf(stderr, "ERROR: invalid distance value of %ld provided\n", dist);
        exit(1);
    }
    
    if (ingroup_pops.size() == 0 || admixers.size() == 0){
        fprintf(stderr, "ERROR: you must provide at least one ingroup and at least one \
admixer.\n");
        exit(1);
    }
    if (mapfile == "" || cladefile == ""){
        fprintf(stderr, "ERROR: you must provide a map file name and clade file name.\n");
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
   
    instream_info is;
    instream_init(is, &fp, bufsize);
    
    // Read file header.
    int num_haplotypes;
    cladeset mask;
    read_header(is, num_haplotypes, mask);
    
    if (minsize < 0 || minsize > num_haplotypes){
        fprintf(stderr, "ERROR: invalid minimum number of haplotypes: %d\n", minsize);
        exit(1);
    }
    if (maxsize != -1 && (maxsize <= minsize || maxsize > num_haplotypes)){
        fprintf(stderr, "ERROR: invalid maximum number of haplotypes: %d\n", maxsize);
        exit(1);
    }
    
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
    unordered_map<cladeset, string> admixer2name;
    
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
                admixer2name.insert(make_pair(admix_clade, *admix_hap));
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
    
    FILE* mapf = fopen(mapfile.c_str(), "w");
    if (mapf == NULL){
        fprintf(stderr, "ERROR: unable to open file %s for writing.\n", mapfile.c_str());
        exit(1);
    }
    FILE* cladef = fopen(cladefile.c_str(), "w");
    if (cladef == NULL){
        fprintf(stderr, "ERROR: unable to open file %s for writing.\n", cladefile.c_str());
    }
    
    fprintf(stderr, "Ingroup haps: %ld Outgroup haps: %ld Admixing haps: %ld\n", 
        ingroup_set.size(), outgroup_set.size(), admix_clades.size());
    
    // Store all trees (and clades) within a set distance of the current site
    // This is needed because when we find a "tagged" clade (containing admixer
    // haplotype(s)), we look back to see how far upstream it extended, possibly
    // minus the admixer(s). We then look how far downstream it extends, possibly
    // minus the admixer(s).
    vector<treeNode*> prevtrees;
    vector<long int> prevpos;
    
    // Map clade to start & end coordinates
    unordered_map<cladeset, long int> admclades_begin;
    unordered_map<cladeset, long int> admclades_end;
    
    // Map clade to TMRCA mean sum & count
    unordered_map<cladeset, float> admclades_tmrca_sum;
    unordered_map<cladeset, float> admclades_adm_tmrca_sum;
    unordered_map<cladeset, float> admclades_tmrca_count;
    unordered_map<cladeset, float> admclades_adm_tmrca_count;
    
    // Map clade to admixer haplotypes contained within
    unordered_map<cladeset, set<string> > admclades_admixer_pop;
    unordered_map<cladeset, set<string> > admclades_admixer_hap;
    
    long int progress = 10000;
    long int prev_printed = 0;
    
    string prevchrom;
    
    // Read from stdin.
    while(!is.finished()){

        string chrom;
        long int pos;
        
        treeNode* tree = new treeNode();
        tree->set_haps(num_haplotypes);
        
        read_sitedata(num_haplotypes, is, chrom, pos, *tree);

        if (pos - prev_printed > progress){
            fprintf(stderr, "Read %s %ld\r", chrom.c_str(), pos);
            prev_printed = pos;
        }

        // Find all admixer haplotypes
        set<treeNode*> admixers_nodes;
        map<treeNode*, string> admixer_node_pops;
        for (vector<cladeset>::iterator ac = admix_clades.begin(); ac != admix_clades.end();
            ++ac){
            treeNode* an = tree->get_smallest_containing(*ac);
            admixers_nodes.insert(an);
            admixer_node_pops.insert(make_pair(an, admixer2pop[*ac]));
        }
        
        unordered_set<cladeset> clades_adm;
        admixed_clades(tree, admixers_all, admix_clades, ingroup_clade, outgroup_clade, 
            clades_adm, require_mutations, polytomy_size);
        
        
        // Which of the clades we're currently tracking was found in the 
        // current site's tree?
        unordered_set<cladeset> admclades_tracking_found;
        
        for (unordered_set<cladeset>::iterator ca = clades_adm.begin(); ca != clades_adm.end();
            ++ca){
            cladeset without_adm = set_diff_bitset(*ca, admixers_all);
            
            if (without_adm.count() < minsize || (maxsize != -1 && without_adm.count() > maxsize)){
                // Don't track this clade.
                continue;
            }
            else{
                // Are we already tracking this clade?
                if (admclades_begin.count(without_adm) > 0){
                    admclades_tracking_found.insert(without_adm);
                    // Get TMRCA -- if the admixed haplotypes exist as a clade
                    // without including admixers, take the TMRCA of that clade.
                    // If not, take the TMRCA of the smallest clade containing
                    // all of them.
                    admclades_end[without_adm] = pos;
                    treeNode* tn = tree->get_smallest_containing(without_adm);
                    admclades_tmrca_sum[without_adm] += (tn->dist_below/(tn->dist_below+tn->dist_above));
                    admclades_tmrca_count[without_adm]++;
                    treeNode* tn2 = tn;
                    while ((tn2->subtree_leaves() & admixers_all).count() == 0 && tn2->parent != NULL){
                        tn2 = tn2->parent;
                    }
                    admclades_adm_tmrca_sum[without_adm] += (tn2->dist_below/(tn2->dist_below+tn2->dist_above));
                    admclades_adm_tmrca_count[without_adm]++;
                    
                    // Determine populations of all admixers in this clade.
                    for (unordered_map<cladeset, string>::iterator a2p = admixer2pop.begin();
                        a2p != admixer2pop.end(); ++a2p){
                        if ((a2p->first & *ca).count() > 0){
                            admclades_admixer_pop[without_adm].insert(a2p->second);
                            admclades_admixer_hap[without_adm].insert(admixer2name[a2p->first]);
                        }
                    }
                }
                else{
                    // Start tracking this clade.
                    admclades_tracking_found.insert(without_adm);
                    
                    // Add TMRCA, etc. information from current site.
                    treeNode* tn = tree->get_smallest_containing(without_adm);
                    if (tn == NULL){
                        fprintf(stderr, "NULL?1\n");
                        exit(1);
                    }
                    admclades_tmrca_sum.insert(make_pair(without_adm, 
                        (tn->dist_below/(tn->dist_below+tn->dist_above))));
                    
                    admclades_tmrca_count.insert(make_pair(without_adm, 1));
                    
                    treeNode* tn2 = tn;
                    while((tn2->subtree_leaves() & admixers_all).count() == 0 && tn2->parent != NULL){
                        tn2 = tn2->parent;
                    }
                    admclades_adm_tmrca_sum.insert(make_pair(without_adm,
                        (tn2->dist_below/(tn2->dist_below+tn2->dist_above))));
                        
                    admclades_adm_tmrca_count.insert(make_pair(without_adm, 1));
                    
                    admclades_end.insert(make_pair(without_adm, pos));
                    
                    admclades_begin.insert(make_pair(without_adm, pos));
                    
                    set<string> s;
                    set<string> s2;
                    admclades_admixer_pop.insert(make_pair(without_adm, s));
                    admclades_admixer_hap.insert(make_pair(without_adm, s2));
                    
                    for (unordered_map<cladeset, string>::iterator a2p = admixer2pop.begin();
                        a2p != admixer2pop.end(); ++a2p){
                        if ((a2p->first & *ca).count() > 0){
                            admclades_admixer_pop[without_adm].insert(a2p->second);
                            admclades_admixer_hap[without_adm].insert(admixer2name[a2p->first]);
                        }
                    }
                    
                    // Look up where clade first appeared in previous trees.
                    for (int i = prevtrees.size()-1; i >= 0; i--){
                        
                        treeNode* parent = prevtrees[i]->get_smallest_containing(without_adm);
                        cladeset others = set_diff_bitset(parent->subtree_leaves(), without_adm);
                        
                        if (others.count() == 0 || (others.count() == (others & admixers_all).count())){
                            // Either no other haplotypes included, or all other
                            // haplotypes included are admixers.
                            admclades_begin[without_adm] = prevpos[i];
                            
                            treeNode* tn = prevtrees[i]->get_smallest_containing(without_adm);
                            if (tn == NULL){
                                fprintf(stderr, "NULL?2\n");
                                exit(1);
                            }
                            admclades_tmrca_sum[without_adm] += (tn->dist_below/(tn->dist_below+tn->dist_above));
                            admclades_tmrca_count[without_adm]++;
                            
                            if (others.count() > 0){
                                // Look up populations of admixers.
                                for (unordered_map<cladeset, string>::iterator a2p = admixer2pop.begin();
                                    a2p != admixer2pop.end(); ++a2p){
                                    if ((a2p->first & others).count() > 0){
                                        admclades_admixer_pop[without_adm].insert(a2p->second);
                                        admclades_admixer_hap[without_adm].insert(admixer2name[a2p->first]);
                                    }
                                }
                            }
                        }
                        else{
                            // Clade contains other haplotypes that aren't admixers.
                            // This doesn't count.
                            break;
                        }
                    }
                }
            }
        }
        // Check to see whether any currently-tracked clades exist without
        // admixers in the current tree.
        for (unordered_map<cladeset, long int>::iterator ab = admclades_begin.begin();
            ab != admclades_begin.end(); ){
            bool rm = false;
            if (admclades_tracking_found.find(ab->first) != admclades_tracking_found.end()){
                // Already found (with admixer haplotypes included).
                // Skip it.
            }
            else if (tree->get_clade_match(ab->first) != NULL){
                // Found without admixer haplotypes.
                admclades_end[ab->first] = pos;
                // Add information about TMRCA, etc.
                treeNode* tn = tree->get_smallest_containing(ab->first);
                if (tn == NULL){
                    fprintf(stderr, "NULL?3\n");
                    exit(1);
                }
                admclades_tmrca_sum[ab->first] += (tn->dist_below/(tn->dist_below+tn->dist_above));
                admclades_tmrca_count[ab->first]++;
                // No admixers found in this clade - leave that alone.
            }
            else{
                // Not found. This one is finished.
                float denom = 1.0;
                if (admclades_tmrca_count[ab->first] > 0){
                    denom = admclades_tmrca_count[ab->first];
                }
                float denom2 = 1.0;
                if (admclades_adm_tmrca_count[ab->first] > 0){
                    denom2 = admclades_adm_tmrca_count[ab->first];
                }
                
                // Compile names of all admixers
                string admhaps = "";
                string admpops = "";
                for (set<string>::iterator h = admclades_admixer_hap[ab->first].begin();
                    h != admclades_admixer_hap[ab->first].end(); ++h){
                    admhaps += *h + ",";
                }
                for (set<string>::iterator p = admclades_admixer_pop[ab->first].begin();
                    p != admclades_admixer_pop[ab->first].end(); ++p){
                    admpops += *p + ",";
                }
                admhaps = admhaps.substr(0, admhaps.length()-1);
                admpops = admpops.substr(0, admpops.length()-1);
                
                // Compile all ingroup populations in this clade
                string igpops_clade = "";
                set<unsigned int> hapinds_clade = bitset2set(ab->first, num_haplotypes);
                set<string> igpops_cl_set;
                for (set<unsigned int>::iterator hic = hapinds_clade.begin();
                    hic != hapinds_clade.end(); ++hic){
                    igpops_cl_set.insert(indv2pop[indvs[*hic]]);
                }
                for (set<string>::iterator ics = igpops_cl_set.begin();
                    ics != igpops_cl_set.end(); ++ics){
                    igpops_clade += *ics + ",";
                }
                igpops_clade = igpops_clade.substr(0, igpops_clade.length()-1);
                
                set<unsigned int> ighaps_clade = bitset2set(ab->first, num_haplotypes);
                for (set<unsigned int>::iterator igc = ighaps_clade.begin();
                    igc != ighaps_clade.end(); ++igc){
                    fprintf(mapf, "%s\t%ld\t%ld\t%ld\t%f\t%f\t%s\t%s\t%s\t%s\n", 
                        chrom.c_str(),
                        ab->second, 
                        admclades_end[ab->first], 
                        admclades_end[ab->first] - ab->second + 1,
                        admclades_tmrca_sum[ab->first]/denom,
                        admclades_adm_tmrca_sum[ab->first]/denom2,
                        indvs[*igc].c_str(),
                        indv2pop[indvs[*igc]].c_str(),
                        admhaps.c_str(),
                        admpops.c_str());
                }
                fprintf(cladef, "%s\t%ld\t%ld\t%ld\t%f\t%f\t%ld\t%s\t%s\n",
                    chrom.c_str(),
                    ab->second,
                    admclades_end[ab->first],
                    admclades_end[ab->first] - ab->second + 1,
                    admclades_tmrca_sum[ab->first]/denom,
                    admclades_adm_tmrca_sum[ab->first]/denom2,
                    ab->first.count(),
                    igpops_clade.c_str(),
                    admpops.c_str());
                
                admclades_end.erase(ab->first);
                admclades_tmrca_sum.erase(ab->first);
                admclades_tmrca_count.erase(ab->first);
                admclades_adm_tmrca_sum.erase(ab->first);
                admclades_adm_tmrca_count.erase(ab->first);
                admclades_admixer_pop[ab->first].clear();
                admclades_admixer_pop.erase(ab->first);
                admclades_admixer_hap[ab->first].clear();
                admclades_admixer_hap.erase(ab->first);
                
                admclades_begin.erase(ab++);
                rm = true;
            }
            if (!rm){
                ++ab;
            }
        }

        // Delete all trees/sites beyond propagation distance.
        bool deleting = prevpos.size() > 0;
        while (deleting){
            if (prevpos.size() == 0){
                deleting = false;
            }
            else{
                if (*prevpos.begin() < pos - dist){
                    delete *prevtrees.begin();
                    prevtrees.erase(prevtrees.begin());
                    prevpos.erase(prevpos.begin());
                }
                else{
                    deleting = false;
                }
            }
        }
        
        prevpos.push_back(pos);
        prevtrees.push_back(tree);
        
        prevchrom = chrom;
    }
    
    // Finish all current runs
    for (unordered_map<cladeset, long int>::iterator ab = admclades_begin.begin();
        ab != admclades_begin.end(); ){
        // Not found. This one is finished.
        float denom = 1.0;
        if (admclades_tmrca_count[ab->first] > 0){
            denom = admclades_tmrca_count[ab->first];
        }
        float denom2 = 1.0;
        if (admclades_adm_tmrca_count[ab->first] > 0){
            denom2 = admclades_adm_tmrca_count[ab->first];
        }
        
        // Compile names of all admixers
        string admhaps = "";
        string admpops = "";
        for (set<string>::iterator h = admclades_admixer_hap[ab->first].begin();
            h != admclades_admixer_hap[ab->first].end(); ++h){
            admhaps += *h + ",";
        }
        for (set<string>::iterator p = admclades_admixer_pop[ab->first].begin();
            p != admclades_admixer_pop[ab->first].end(); ++p){
            admpops += *p + ",";
        }
        admhaps = admhaps.substr(0, admhaps.length()-1);
        admpops = admpops.substr(0, admpops.length()-1);
        
        // Compile all ingroup populations in this clade
        string igpops_clade = "";
        set<unsigned int> hapinds_clade = bitset2set(ab->first, num_haplotypes);
        set<string> igpops_cl_set;
        for (set<unsigned int>::iterator hic = hapinds_clade.begin();
            hic != hapinds_clade.end(); ++hic){
            igpops_cl_set.insert(indv2pop[indvs[*hic]]);
        }
        for (set<string>::iterator ics = igpops_cl_set.begin();
            ics != igpops_cl_set.end(); ++ics){
            igpops_clade += *ics + ",";
        }
        igpops_clade = igpops_clade.substr(0, igpops_clade.length()-1);
        
        set<unsigned int> ighaps_clade = bitset2set(ab->first, num_haplotypes);
        for (set<unsigned int>::iterator igc = ighaps_clade.begin();
            igc != ighaps_clade.end(); ++igc){
                fprintf(mapf, "%s\t%ld\t%ld\t%ld\t%f\t%f\t%s\t%s\t%s\t%s\n", 
                    prevchrom.c_str(),
                    ab->second, 
                    admclades_end[ab->first], 
                    admclades_end[ab->first] - ab->second + 1,
                    admclades_tmrca_sum[ab->first]/denom,
                    admclades_adm_tmrca_sum[ab->first]/denom2,
                    indvs[*igc].c_str(),
                    indv2pop[indvs[*igc]].c_str(),
                    admhaps.c_str(),
                    admpops.c_str());
        }
        fprintf(cladef, "%s\t%ld\t%ld\t%ld\t%f\t%f\t%ld\t%s\t%s\n",
            prevchrom.c_str(),
            ab->second,
            admclades_end[ab->first],
            admclades_end[ab->first] - ab->second + 1,
            admclades_tmrca_sum[ab->first]/denom,
            admclades_adm_tmrca_sum[ab->first]/denom2,
            ab->first.count(),
            igpops_clade.c_str(),
            admpops.c_str());
        
        admclades_end.erase(ab->first);
        admclades_tmrca_sum.erase(ab->first);
        admclades_tmrca_count.erase(ab->first);
        admclades_admixer_pop[ab->first].clear();
        admclades_admixer_pop.erase(ab->first);
        admclades_admixer_hap[ab->first].clear();
        admclades_admixer_hap.erase(ab->first);
        
        admclades_begin.erase(ab++);
    }
    
    fclose(mapf);
    fclose(cladef);
    
    return 0;
}
