#include <getopt.h>
#include <argp.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include "sets.h"
#include "treeNode.h"
#include <math.h>
#include <zlib.h>
#include "serialize.h"
#include "common.h"
#include <unordered_set>

using std::cout;
using std::endl;
using namespace std;

// Store information about an individual clade
struct haplen_dat{
    int traversal_index;
    float size;
    long int persistence;
    float tmrca;
    float tmrca_p;
    bool selected_haps;
    float tmrca_avgtree;
    set<long int> sites;
    
    haplen_dat(int ti, float s, float pers, float t, float p, bool sh, float ta){
        this->traversal_index = ti;
        this->size = s;
        this->persistence = pers;
        this->tmrca = t;
        this->tmrca_p = p;
        this->selected_haps = sh;
        this->tmrca_avgtree = ta;
        this->sites.clear();
    };
};

bool operator<(const haplen_dat& h1, const haplen_dat& h2){
    if (h1.persistence > h2.persistence){
        return true;
    }
    else{
        return false;
    }
}

bool p_comp(const haplen_dat& h1, const haplen_dat& h2){
    if (h1.tmrca_p < h2.tmrca_p){
        return true;
    }
    else{
        return false;
    }
}

void compile_haplens(treeNode* tree, 
    int& traversal_index, 
    bool freqs_given,
    map<pair<long double, long double>, vector<long double> >& ages_freqs,  
    bool use_selected_haps, 
    const cladeset& selected_haps,
    bool use_excl_haps,
    const cladeset& excl_haps,
    vector<haplen_dat>& haps,
    bool use_avg_tree,
    treeNode* avg_tree,
    int num_haplotypes){
    
    if (tree->parent != NULL && tree->persistence > 0){
        cladeset st = tree->subtree_leaves();
        
        if (st.count() > 1 && st.count() < tree->num_haps){
        //if (st.count() >= minsize && tree->persistence >= minlen && (maxsize == -1 ||
        //    st.count() <= maxsize) && (!use_selected_haps || (st & selected_haps).count() > 0)){
            long double tmrca = tree->dist_below/(tree->dist_below+tree->dist_above);
            if (tree->dist_below + tree->dist_above == 0){
                tmrca = 0;
            }
            
            if (true){
            //if (tmrca <= tmrca_cutoff){
                if (freqs_given){
                    // Calculate probability of clade size based on age
                    if (tmrca == 0 || tree->dist_below + tree->dist_above == 0 || tmrca < ages_freqs.begin()->first.first){
                        tmrca = ages_freqs.begin()->first.first;
                    }
                    bool found = false;
                    for (map<pair<long double, long double>, vector<long double> >::iterator af = 
                        ages_freqs.begin(); af != ages_freqs.end(); ++af){
                        if (tmrca >= af->first.first && tmrca < af->first.second){
                            found = true;
                            double probsum = 0;
                            for (int i = st.count()-1; i < af->second.size(); ++i){
                                probsum += af->second[i];
                            }
                            if (true){
                            //if (p_cutoff == 1 || probsum < p_cutoff){
                                bool has_selected_haps = true;
                                if (use_selected_haps){
                                    if ((st & selected_haps).count() > 0){
                                        has_selected_haps = true;
                                    }
                                    else{
                                        has_selected_haps = false;
                                    }
                                }
                                // Apply exclude haplotypes filter if applicable
                                if (has_selected_haps && use_excl_haps){
                                    if ((st & excl_haps).count() > 0){
                                        has_selected_haps = false;
                                    }
                                }
                                
                                // Find TMRCA of clade in average tree
                                float tmrca_avg = -1;
                                if (use_avg_tree){
                                    cladeset st = tree->subtree_leaves();
                                    treeNode* n_avg = avg_tree->get_smallest_containing(st);
                                    tmrca_avg = n_avg->dist_below/(n_avg->dist_below+n_avg->dist_above);
                                }
                                
                                haps.push_back(haplen_dat(traversal_index, (float)st.count()/(float)num_haplotypes, 
                                    tree->persistence, tmrca, probsum, has_selected_haps, tmrca_avg));
                                
                                for (set<long int>::iterator site = 
                                    tree->mutations.begin(); site != tree->mutations.end();
                                    ++site){
                                    haps.rbegin()->sites.insert(*site);
                                }
                                //fprintf(stdout, "%s\t%ld\t%d\t%ld\t%ld\t%Lf\t%f\n", chrom.c_str(), pos,
                                //    traversal_index, st.count(), tree->persistence, tmrca, probsum);
                            }
                        }
                    }
                    if (!found){
                        fprintf(stderr, "\n?? !found\n");
                        fprintf(stderr, "%.8Lf\n", tmrca);
                        fprintf(stderr, "%.8f %.8f\n", tree->dist_below, tree->dist_above);
                        for (map<pair<long double, long double>, vector<long double> >::iterator af = 
                            ages_freqs.begin(); af != ages_freqs.end(); ++af){
                            fprintf(stderr, "%.8Lf %.8Lf\n", af->first.first, af->first.second);
                        }
                        exit(1);
                    }
                }
                else{
                    //fprintf(stdout, "%s\t%ld\t%d\t%ld\t%ld\t%Lf\n", chrom.c_str(), pos,
                    //    traversal_index, st.count(), tree->persistence, tmrca);
                    bool has_selected_haps = true;
                    if (use_selected_haps){
                        if ((st & selected_haps).count() > 0){
                            has_selected_haps = true;
                        }
                        else{
                            has_selected_haps = false;
                        }
                    }
                    // Find TMRCA of clade in average tree
                    float tmrca_avg = -1;
                    if (use_avg_tree){
                        cladeset st = tree->subtree_leaves();
                        treeNode* n_avg = avg_tree->get_smallest_containing(st);
                        tmrca_avg = n_avg->dist_below/(n_avg->dist_below+n_avg->dist_above);
                    }
                    
                    haps.push_back(haplen_dat(traversal_index, (float)st.count()/(float)num_haplotypes,
                        tree->persistence, tmrca, -1, has_selected_haps, tmrca_avg));
                    
                    for (set<long int>::iterator site = tree->mutations.begin();
                        site != tree->mutations.end(); ++site){
                        haps.rbegin()->sites.insert(*site);
                    }
                }
            }
        }
    }
    traversal_index++;
    for (vector<treeNode*>::iterator child = tree->children.begin(); child !=
        tree->children.end(); ++child){
        compile_haplens(*child, traversal_index, freqs_given, 
            ages_freqs, use_selected_haps, selected_haps, use_excl_haps, excl_haps,
            haps, use_avg_tree, avg_tree, num_haplotypes);
    }
}

void print_haplens(vector<haplen_dat>& haps, 
    string& chrom, 
    long int pos, 
    bool freqs_given,
    float minsize, 
    float maxsize, 
    long int minlen, 
    float p_cutoff,
    float tmrca_cutoff,
    bool use_avg_tree){
    
    // Map traversal index to rank percentile (two different ways)
    map<int, float> rank_persistence;
    map<int, float> rank_tmrca_p;
    
    float pers_sum = 0;
    float tmrca_p_sum = 0;
    for (vector<haplen_dat>::iterator h = haps.begin(); h != haps.end(); ++h){
        pers_sum += h->persistence;
        if (freqs_given){
            tmrca_p_sum += h->tmrca_p;
        }
    }
    float tmrca_p_mean = 0;
    float pers_mean = 0;
    if (haps.size() > 0){
        pers_mean = (float)pers_sum/(float)haps.size();
        if (freqs_given){
            tmrca_p_mean = (float)tmrca_p_sum/(float)haps.size();
        }
    }
    
    for (vector<haplen_dat>::iterator hap = haps.begin(); hap != haps.end(); ++hap){
        if (hap->selected_haps && hap->tmrca <= tmrca_cutoff && 
            (minsize == -1 || hap->size >= minsize) &&
            (maxsize == -1 || hap->size <= maxsize) && 
            hap->persistence >= minlen &&
            (!freqs_given || hap->tmrca_p <= p_cutoff)){
            
            string sitestr = "";
            for (set<long int>::iterator site = hap->sites.begin();
                site != hap->sites.end(); ++site){
                char buf[10];
                sprintf(buf, "%ld", *site);
                string bufstr = buf;
                sitestr += bufstr + ",";
            }
            if (sitestr.length() == 0){
                sitestr = "NA";
            }
            else{
                sitestr = sitestr.substr(0, sitestr.length()-1);
            }
            
            fprintf(stdout, "%s\t%ld\t%s\t%d\t%f\t%ld\t%f", chrom.c_str(), pos, sitestr.c_str(),
                hap->traversal_index, hap->size, hap->persistence, hap->tmrca);
            if (freqs_given){
                fprintf(stdout, "\t%f", hap->tmrca_p);
            } 
            float p_persistence = exp(-1.0/pers_mean*hap->persistence);
            fprintf(stdout, "\t%f", p_persistence);
            if (freqs_given){
                float p_tmrca = 1-exp(-1.0/tmrca_p_mean*hap->tmrca_p);
                fprintf(stdout, "\t%f\t%f", p_tmrca,
                    p_persistence*p_tmrca);
            }
            if (use_avg_tree){
                fprintf(stdout, "\t%f", hap->tmrca_avgtree);
            }
            fprintf(stdout, "\n");
        }
    }
    
    
    return;
    
    
    static int count = 0;
    count++;
    if (count == 1000){
    for (vector<haplen_dat>::iterator h = haps.begin(); h != haps.end(); ++h){
        fprintf(stderr, "%f\t%f\t%ld\n", h->tmrca, h->size, h->persistence);   
    }
    exit(1);
    }
    if (freqs_given){
        // Sort in increasing order of TMRCA p-value
        sort(haps.begin(), haps.end(), p_comp);
        int i = 0;
        for (vector<haplen_dat>::iterator h = haps.begin(); h != haps.end(); ++h){
            int count_geq = i;
            vector<haplen_dat>::iterator h2 = h;
            while (h2 != haps.end()){
                if (h2->tmrca_p == h->tmrca_p){
                    ++count_geq;
                }
                else{
                    break;
                }
                ++h2;
            }
            rank_tmrca_p.insert(make_pair(h->traversal_index, (float)count_geq/(float)haps.size()));
        }
    }
    
    // Sort in decreasing order of persistence
    sort(haps.begin(), haps.end());
    int i = 0;
    for (vector<haplen_dat>::iterator h = haps.begin(); h != haps.end(); ++h){
        int count_geq = i;
        vector<haplen_dat>::iterator h2 = h;
        while (h2 != haps.end()){
            if (h2->persistence == h->persistence){
                ++count_geq;
            }
            else{
                break;
            }
            ++h2;
        }
        rank_persistence.insert(make_pair(h->traversal_index, (float)count_geq/(float)haps.size()));
        ++i;
    }

    for (vector<haplen_dat>::iterator hap = haps.begin(); hap != haps.end(); ++hap){
        if (hap->selected_haps && hap->tmrca <= tmrca_cutoff && 
            (minsize == -1 || hap->size >= minsize) &&
            (maxsize == -1 || hap->size <= maxsize) && 
            hap->persistence >= minlen &&
            (!freqs_given || hap->tmrca_p <= p_cutoff)){
            
            string sitestr = "";
            for (set<long int>::iterator site = hap->sites.begin();
                site != hap->sites.end(); ++site){
                char buf[10];
                sprintf(buf, "%ld", *site);
                string bufstr = buf;
                sitestr += bufstr + ",";
            }
            if (sitestr.length() == 0){
                sitestr = "NA";
            }
            else{
                sitestr = sitestr.substr(0, sitestr.length()-1);
            }
            fprintf(stdout, "%s\t%ld\t%s\t%d\t%f\t%ld\t%f", chrom.c_str(), pos, 
                sitestr.c_str(), hap->traversal_index,
                hap->size, hap->persistence, hap->tmrca);
            if (freqs_given){
                fprintf(stdout, "\t%f", hap->tmrca_p);
            } 
            fprintf(stdout, "\t%f", rank_persistence[hap->traversal_index]);
            if (freqs_given){
                fprintf(stdout, "\t%f\t%f", rank_tmrca_p[hap->traversal_index],
                    rank_persistence[hap->traversal_index]*rank_tmrca_p[hap->traversal_index]);
            }
            fprintf(stdout, "\n");
        }
    }
}

void parse_freqs_file(string freqsfile, map<pair<long double, long double>, vector<long double> >& ages_freqs){
    std::ifstream infile;
    infile.open(freqsfile.c_str(), std::ifstream::in);
    if (!infile.good()){
        fprintf(stderr, "ERROR: could not open file %s for reading\n", freqsfile.c_str());
        exit(1);
    }
    bool headerline = true;
    for (string line; getline(infile, line);){
        if (headerline){
            headerline = false;
            continue;
        }
        else{
            // Split by tab
            istringstream tabsplitter(line);
            string field;
            int field_index = 0;
        
            long double lower;
            long double upper;
            vector<long double> vals;
            
            while(std::getline(tabsplitter, field, '\t')){
                if (field_index == 0){
                    // Bin ranges
                    size_t space = field.find(" ");
                    if (space == string::npos){
                        fprintf(stderr, "ERROR: could not parse bin ranges in input file\n");
                        exit(1);
                    }
                    string lowerstr = field.substr(0, space);
                    string upperstr = field.substr(space+1, field.length()-space);
                    lower = strtold(lowerstr.c_str(), NULL);
                    upper = strtold(upperstr.c_str(), NULL);
                }
                else{
                    // Value.
                    float val = strtold(field.c_str(), NULL);
                    vals.push_back(val);
                }
                field_index++;
            }
            
            ages_freqs.insert(make_pair(make_pair(lower, upper), vals));
        }   
    }
}

/**
 * Print a help message to the terminal and exit.
 */
void help(int code){
   fprintf(stderr, "dump_haplens [OPTIONS]\n");
   fprintf(stderr, "Once haplotype lengths have been computed and stored (using the \
haplens program), this prints clade sizes and persistences in tab-delimited format.\n");
    fprintf(stderr, "[OPTIONS]:\n");
    fprintf(stderr, "   --bufsize -b (OPTIONAL) The number of characters to read from the \
input file at a time\n");
    fprintf(stderr, "   --minsize -s (OPTIONAL) minimum number of haplotypes for a clade \
to contain in order to be printed (must be > 1)\n");
    fprintf(stderr, "   --maxsize -S (OPTIONAL) maximum number of haplotypes for a clade \
to contain in order to be printed (must be < num hapotypes -1)\n");
    fprintf(stderr, "   --length -l (OPTIONAL) minimum number of bases for a clade to \
persist in order to be printed (must be >= 1)\n");
    fprintf(stderr, "   --freqs -f (OPTIONAL) to calculate probabilities of clade sizes \
based on TMRCA, provide table computed by utilities/calc_ages_freqs.py\n");
    fprintf(stderr, "   --p_cutoff -P (OPTIONAL) if providing a frequency table (-f argument), \
the maximum p-value allowed to print a clade\n");
    fprintf(stderr, "   --tmrca_cutoff, -t (OPTIONAL) a TMRCA value (as percent of total divergence \
to outgroup) above which not to print clades\n");
    fprintf(stderr, "   --indvs -v (OPTIONAL) A file listing haplotype names - required \
if using -p and -i options\n");
    fprintf(stderr, "   --pops -p (OPTIONAL) file mapping haplotype IDs to populations - \
required if using -i and -v options\n");
    fprintf(stderr, "   --include -i (OPTIONAL) require clades to include haplotypes \
from a population given in the --pops file. Requires -p and -v options; may specify more than \
once (use -i before each population name).\n");
    fprintf(stderr, "   --exclude -e (OPTIONAL) require clades to exclude haplotypes \
from a population given in the --pops file. Requires -p and -v options; may specify more than \
once (use -i before each population name).\n");
    fprintf(stderr, "   --tree -T (OPTIONAL) if you run upgma on the entire chromosome, \
you can provide the resulting tree as a parameter here, and each printed clade will \
list its TMRCA in the average tree. Higher TMRCA means more unusual grouping of haplotypes.\n");
    exit(code);
}

int main(int argc, char *argv[]) {
    static struct option long_options[] = {
       {"bufsize", required_argument, 0, 'b'},
       {"help", optional_argument, 0, 'h'},
       {"minsize", required_argument, 0, 's'},
       {"maxsize", required_argument, 0, 'S'},
       {"length", required_argument, 0, 'l'},
       {"freqs", required_argument, 0, 'f'},
       {"p_cutoff", required_argument, 0, 'P'},
       {"tmrca_cutoff", required_argument, 0, 't'},
       {"indvs", required_argument, 0, 'v'},
       {"pops", required_argument, 0, 'p'},
       {"include", required_argument, 0, 'i'},
       {"exclude", required_argument, 0, 'e'},
       {"tree", required_argument, 0, 'T'},
       {0, 0, 0, 0} 
    };
    
    int bufsize = 1048576;
    long int minsize = 2;
    long int maxsize = -1;
    long int minlen = 1;
    int option_index = 0;
    string freqsfile;
    bool freqs_given = false;
    float p_cutoff = 1.0;
    float tmrca_cutoff = 1.0;
    string hapsfile;
    bool hapsfile_given = false;
    string popsfile;
    bool popsfile_given = false;
    vector<string> include_pops;
    vector<string> exclude_pops;
    string treefilename;
    bool treefile_given = false;
    
    int ch;
    /*
    if (argc == 1){
        help(0);
    }
    */
    while((ch = getopt_long(argc, argv, "b:s:S:l:f:P:t:v:p:i:e:T:h", long_options, &option_index )) != -1){
        switch(ch){
            case 0:
                // This option set a flag. No need to do anything here.
                break;
            case 'b':
                bufsize = atoi(optarg);
                break;
            case 's':
                minsize = atol(optarg);
                break;
            case 'S':
                maxsize = atol(optarg);
                break;
            case 'l':
                minlen = atol(optarg);
                break;
            case 'f':
                freqsfile = optarg;
                freqs_given = true;
                break;
            case 'P':
                p_cutoff = atof(optarg);
                break;
            case 't':
                tmrca_cutoff = atof(optarg);
                break;
            case 'v':
                hapsfile = optarg;
                hapsfile_given = true;
                break;
            case 'p':
                popsfile = optarg;
                popsfile_given = true;
                break;
            case 'i':
                include_pops.push_back((string)optarg);
                break;
            case 'e':
                exclude_pops.push_back((string)optarg);
                break;
            case 'T':
                treefilename = optarg;
                treefile_given = true;
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
    
    if (minsize < 2){
        fprintf(stderr, "ERROR: minimum size must be greater than 1\n");
        exit(1);
    }
    if (minlen < 1){
        fprintf(stderr, "ERROR: minimum length must be greater than or equal \
to 1\n");
        exit(1);
    }
    if (p_cutoff < 1.0 && !freqs_given){
        fprintf(stderr, "ERROR: you cannot specify a p-value cutoff without a frequency table\n");
        exit(1);
    }
    if (tmrca_cutoff <= 0.0 || tmrca_cutoff > 1.0){
        fprintf(stderr, "ERROR: invalid TMRCA cutoff given.\n");
        exit(1);
    }
    if (include_pops.size() > 0 && (!hapsfile_given || !popsfile_given)){
        fprintf(stderr, "ERROR: you have specified one or more populations to require \
in target clades, but have not specified the haplotype names file -v and/or have \
not specified the haplotype/population mapping file -p.\n");
        exit(1);
    }
    if (exclude_pops.size() > 0 && (!hapsfile_given || !popsfile_given)){
        fprintf(stderr, "ERROR: you have specified one or more populations to exclude \
in target clades, but have not specified the haplotype names file -v and/or have \
not specified the haplotype/population mapping file -p.\n");
        exit(1);
    }
    FILE *instream = stdin;
    if (instream == NULL){
        fprintf(stderr, "Error opening input stream\n");
        exit(1);
    }
    gzFile fp = gzdopen(fileno(instream), "rb");
    if (fp == NULL){
        fprintf(stderr, "ERROR: unable to read from stdin.\n");
        exit(1);
    }
    
    map<pair<long double, long double>, vector<long double> > ages_freqs;
    if (freqs_given){
        parse_freqs_file(freqsfile, ages_freqs);
    }
    
    instream_info is;
    instream_init(is, &fp, bufsize);
    
    // Read file header.
    int num_haplotypes;
    cladeset mask;
    read_header(is, num_haplotypes, mask);
    
    
    fprintf(stderr, "Read %d haplotypes\n", num_haplotypes);
    
    if (maxsize != -1 && maxsize > num_haplotypes-1){
        fprintf(stderr, "ERROR: maximum clade size must be <= %d\n", num_haplotypes-1);
        exit(1);
    }
    
    float maxsize_float = -1;
    float minsize_float = -1;
    if (minsize != -1){
        minsize_float = (float)minsize/(float)num_haplotypes;
    }
    if (maxsize != -1){
        maxsize_float = (float)maxsize/(float)num_haplotypes;
    }
    
    treeNode* avg_tree;
    
    if (treefile_given){
        // Load tree from file.
        gzFile treefile;
        treefile = gzopen(treefilename.c_str(), "rb");
        if (!treefile){
            fprintf(stderr, "ERROR opening file %s for reading.\n", treefilename.c_str());
            exit(1);
        }
        instream_info is;
        instream_init(is, &treefile, bufsize);
        read_header(is, num_haplotypes, mask);
        avg_tree = new treeNode();
        
        unserialize_treeNode(*avg_tree, num_haplotypes, is);
        avg_tree->setRoot();
        unnorm_branchlens(avg_tree);
    }
    
    vector<string> haps;
    map<string, int> str2hap;
    if (hapsfile_given){
        parse_indvs(haps, hapsfile);
        parse_indvs_map(str2hap, hapsfile);
    }
    
    // Detect ploidy from haplotypes file
    string firsthap = "";
    int ploidy = 1;
    for (int i = 0; i < haps.size(); ++i){
        int last = haps[i].rfind("-");
        string hap = haps[i];
        if (last != string::npos){
            hap = hap.substr(0, last);
        }
        if (firsthap == ""){
            firsthap = hap;
        }
        else if (hap == firsthap){
            ploidy++;
        }
        else{
            break;
        }
    }
    
    map<string, string> indv2pop; 
    map<string, vector<string> > pop2indv;
    set<unsigned int> selected_haps_set;
    cladeset selected_haps;
    bool use_selected_haps = false;
    set<unsigned int> excl_haps_set;
    cladeset excl_haps;
    bool use_excl_haps = false;
    
    if (popsfile_given){
        parse_pops(indv2pop, pop2indv, popsfile, ploidy);
    }
    
    if (popsfile_given && hapsfile_given && include_pops.size() > 0){
        use_selected_haps = true;
        // Now translate the selected populations into a bitset of 
        // eligible haplotypes.
        for (vector<string>::iterator pop = include_pops.begin(); 
            pop != include_pops.end(); ++pop){
            for (vector<string>::iterator hapname = pop2indv[*pop].begin();
                hapname != pop2indv[*pop].end(); ++hapname){
                int hapind = str2hap[*hapname];
                //selected_haps.set(num_haplotypes-1-hapind);
                selected_haps_set.insert((unsigned int)hapind);
            }   
        }
        selected_haps = set2bitset(selected_haps_set, num_haplotypes);
    }
    if (popsfile_given && hapsfile_given && exclude_pops.size() > 0){
        use_excl_haps = true;
        for (vector<string>::iterator pop = exclude_pops.begin();
            pop != exclude_pops.end(); ++pop){
            for (vector<string>::iterator hapname = pop2indv[*pop].begin();
                hapname != pop2indv[*pop].end(); ++hapname){
                int hapind = str2hap[*hapname];
                excl_haps_set.insert((unsigned int)hapind);
            }
        }
        excl_haps = set2bitset(excl_haps_set, num_haplotypes);
    }
    
    
    long int last_printed = 0;
    long int progress = 5000;
    
    fprintf(stdout, "#chrom\tpos\tsites\ttrav_index\thaps\tpersistence\ttmrca");
    if (freqs_given){
        fprintf(stdout, "\tp");
    }
    fprintf(stdout, "\trank_pers");
    if (freqs_given){
        fprintf(stdout, "\trank_tmrca_p\tcombined_p");
    }
    if (treefile_given){
        fprintf(stdout, "\tavg_tree_tmrca");
    }
    fprintf(stdout, "\n");
    
    while(!is.finished()){
        string chrom;
        long int pos;
        treeNode* tree = new treeNode();
        tree->set_haps(num_haplotypes);
        read_sitedata(num_haplotypes, is, chrom, pos, *tree);
        int ti = 0;
        vector<haplen_dat> haps;
        compile_haplens(tree, ti, freqs_given, ages_freqs, use_selected_haps,
            selected_haps, use_excl_haps, excl_haps, haps, treefile_given, avg_tree, num_haplotypes);
        print_haplens(haps, chrom, pos, freqs_given, minsize_float, 
            maxsize_float, minlen,
            p_cutoff, tmrca_cutoff, treefile_given);
        delete tree;
    }
    
}
