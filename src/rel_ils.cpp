/**
 * Computes the extent of ILS between groups in ARG output.
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
    fprintf(stderr, "rel_ils [OPTIONS]\n");
    fprintf(stderr, "Given SARGE output, and populations of interest, counts the number \
of sites at which the populations undergo ILS with other populations.\n");
    fprintf(stderr, "--indvs -v A file listing haplotype IDs (required)\n");
    fprintf(stderr, "--pops -p A file mapping haplotype IDs to population names (required)\n");
    fprintf(stderr, "--ingroup -i A population of interest (can specify more than once)\n");
    fprintf(stderr, "--exclude_bed -b A BED file of regions to exclude from analysis (OPTIONAL)\n");
    fprintf(stderr, "--out_hap -H The name of the output file containing per-haplotype ILS information \
(REQUIRED)\n");
    fprintf(stderr, "--out_pop -P The name of the output file containing per-population ILS information \
(REQUIRED)\n");
    fprintf(stderr, "--bins -B The number of bins to use for histograms (default = 1000)\n");
    fprintf(stderr, "--shorten -s A file listing branch shortening values (OPTIONAL). \
Each line should list a haplotype, then tab, then branch shortening value (as fraction of \
total height of the tree)\n");
exit(code);
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
       {"help", optional_argument, 0, 'h'},
       {"indvs", required_argument, 0, 'v'},
       {"pops", required_argument, 0, 'p'},
       {"ingroup", required_argument, 0, 'i'},
       {"exclude_bed", required_argument, 0, 'b'},
       {"out_hap", required_argument, 0, 'H'},
       {"out_pop", required_argument, 0, 'P'},
       {"bins", required_argument, 0, 'B'},
       {"shorten", required_argument, 0, 's'},
       {0, 0, 0, 0} 
    };
    
    // Set default values
    int bufsize = 1048576;
    
    string indvfilename;
    string popfilename;
    vector<string> ingroup;
    string exclbedfilename;
    bool exclbed_given = false;
    string outpopname;
    string outhapname;
    string brshortenfilename;
    bool brshorten_given = false;
    int bins = 1000;
    
    int option_index = 0;
    int ch;
    
    if (argc == 1){
        help(0);
    }
    while((ch = getopt_long(argc, argv, "v:p:i:b:H:P:B:s:h", long_options, &option_index )) != -1){
        switch(ch){
            case 0:
                // This option set a flag. No need to do anything here.
                break;
            case '?':
                //help(0);
                break;
            case 'h':
                help(0);
                break;
            case 'v':
                indvfilename = optarg;
                break;
            case 'p':
                popfilename = optarg;
                break;
            case 'i':
                ingroup.push_back((string)optarg);
                break;
            case 'b':
                exclbedfilename = optarg;
                exclbed_given = true;
                break;
            case 'H':
                outhapname = optarg;
                break;
            case 'P':
                outpopname = optarg;
                break;
            case 'B':
                bins = atoi(optarg);
                break;
            case 's':
                brshortenfilename = optarg;
                brshorten_given = true;
                break;
            default:
                help(0);
        }    
    }
    
    if (indvfilename.length() == 0){
        fprintf(stderr, "ERROR: you must provide a file listing haplotype names\n");
        exit(1);
    }
    if (popfilename.length() == 0){
        fprintf(stderr, "ERROR: you must provide a file mapping haplotype names to \
population IDs.\n");
        exit(1);
    }
    if (ingroup.size() == 0){
        fprintf(stderr, "ERROR: you must provide one or more ingroup populations.\n");
        exit(1);
    }
    if (outhapname.length() == 0){
        fprintf(stderr, "ERROR: you must provide an output file name for haplotype information\n");
        exit(1);
    }
    if (outpopname.length() == 0){
        fprintf(stderr, "ERROR: you must provide an output file name for population information\n");
        exit(1);
    }
    
    FILE* outhapf = fopen(outhapname.c_str(), "w");
    if (outhapf == NULL){
        fprintf(stderr, "ERROR opening %s for writing.\n", outhapname.c_str());
        exit(1);
    }
    FILE* outpopf = fopen(outpopname.c_str(), "w");
    if (outpopf == NULL){
        fprintf(stderr, "ERROR opening %s for writing.\n", outpopname.c_str());
        exit(1);
    }
    string outhaphistname = outhapname + ".hist";
    FILE* outhaphistf = fopen(outhaphistname.c_str(), "w");
    if (outhaphistf == NULL){
        fprintf(stderr, "ERROR opening %s for writing.\n", outhaphistname.c_str());
        exit(1);
    }    
    string outpophistname = outpopname + ".hist";
    FILE* outpophistf = fopen(outpophistname.c_str(), "w");
    if (outpophistf == NULL){
        fprintf(stderr, "ERROR opening %s for writing.\n", outpophistname.c_str());
        exit(1);
    }
    
    // Parse input files
    vector<string> indvs;
    parse_indvs(indvs, indvfilename);
    map<string, int> indv2ind;
    for (int i = 0; i < indvs.size(); ++i){
        indv2ind.insert(make_pair(indvs[i], i));
    }
    
    map<string, string> indv2pop; 
    map<string, vector<string> > pop2indv;
    parse_pops(indv2pop, pop2indv, popfilename, 2);
    
    map<string, set<pair<long int, long int> > > bed;
    if (exclbed_given){
        parse_exclude_bed(exclbedfilename, bed);
    }
    set<pair<long int, long int> >::iterator cur_bed;
    string cur_chrom = "";
    
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
    
    // Read branch shortening values
    unordered_map<cladeset, float> brshorten_vals; 
    if (brshorten_given){
        parse_brshorten(brshorten_vals, indvs, brshortenfilename, num_haplotypes);
    }
    
    // Create bitsets for lookup of ingroup haplotypes
    // Map population name to all haplotypes within
    map<string, cladeset> pop_clades;
    // Map population name to bitsets for individual haplotypes
    map<string, vector<cladeset> > ingroup_clades;
    unordered_map<cladeset, string> hapnames_cl;
    
    cladeset others;
    for (int i = 0; i < num_haplotypes; ++i){
        others.set(num_haplotypes-1-i);
    }
    
    for (vector<string>::iterator pop = ingroup.begin(); pop != ingroup.end(); ++pop){
        cladeset pop_clade;
        vector<cladeset> pop_indvs;
        for (vector<string>::iterator hap = pop2indv[*pop].begin();
            hap != pop2indv[*pop].end(); ++hap){
            pop_clade.set(num_haplotypes-1-indv2ind[*hap]);
            cladeset cl;
            cl.set(num_haplotypes-1-indv2ind[*hap]);
            others.reset(num_haplotypes-1-indv2ind[*hap]);
            pop_indvs.push_back(cl);
            hapnames_cl.insert(make_pair(cl, *hap));
        }
        pop_clades.insert(make_pair(*pop, pop_clade));
        ingroup_clades.insert(make_pair(*pop, pop_indvs));
    }
    pop_clades.insert(make_pair("OTHER", others));
    
    // Store results
    // For each haplotype, what population does its closest neighbor(s) belong
    // to? 
    map<string, map<string, long int> > haps_ils_count;
    
    // For computing mean TMRCAs
    map<string, map<string, float> > haps_ils_tmrca_sum;

    // For each population (if sorted), what population do its closest neighbor(s)
    // belong to?
    map<string, map<string, long int> > pops_ils_count;
    
    // For computing mean TMRCAs
    map<string, map<string, float> > pops_ils_tmrca_sum;
    
    // For creating histogram of TMRCAs
    map<string, map<string, map<float, float> > > haps_ils_tmrca_hist;
    map<string, map<string, map<float, float> > > pops_ils_tmrca_hist;
    
    long int progress = 10000;
    long int last_printed = 0;
    
    // Read from stdin.
    while(!is.finished()){

        string chrom;
        long int pos;
        
        treeNode* tree = new treeNode();
        tree->set_haps(num_haplotypes);
        
        read_sitedata(num_haplotypes, is, chrom, pos, *tree);
        
        if (pos > last_printed + progress){
            fprintf(stderr, "Read site %ld\r", pos);
            last_printed = pos;
        }
        
        bool skip_site = false;
        if (exclbed_given){
            if (cur_chrom == ""){
                if (bed.count(chrom) > 0){
                    cur_chrom = chrom;
                    cur_bed = bed[cur_chrom].begin();
                }
            }
            else if (chrom != cur_chrom){
                if (bed.count(chrom) > 0){
                    cur_bed = bed[cur_chrom].begin();
                    cur_chrom = chrom;
                }
            }
            if (cur_chrom != ""){
                while(cur_bed->second < pos && cur_bed != bed[cur_chrom].end()){
                    cur_bed++;
                }
                if (cur_bed->first+1 <= pos && cur_bed->second >= pos){
                    skip_site = true;
                }
            }
        }
        
        if (!skip_site){
            // Look for individual haplotypes
            for (map<string, vector<cladeset> >::iterator igc = ingroup_clades.begin();
                igc != ingroup_clades.end(); ++igc){
                
                for (vector<cladeset>::iterator hap = igc->second.begin(); hap != igc->second.end();
                    ++hap){
                    
                    if (haps_ils_count.count(hapnames_cl[*hap]) == 0){
                        map<string, long int> m;
                        haps_ils_count.insert(make_pair(hapnames_cl[*hap], m));
                        map<string, float> m2;
                        haps_ils_tmrca_sum.insert(make_pair(hapnames_cl[*hap], m2));
                        map<string, map<float, float> > m3;
                        haps_ils_tmrca_hist.insert(make_pair(hapnames_cl[*hap], m3));
                    }
                    
                    // Find haplotype in tree.
                    treeNode* tn = tree->get_smallest_containing(*hap);
                    if (tn->subtree_leaves().count() == 1){
                        tn = tn->parent;
                    }
                    
                    cladeset st = tn->subtree_leaves();
                    
                    // ILS if contains other pop(s) before members of own pop
                    if ((st & pop_clades[igc->first]).count() < st.count()){
                        for (map<string, cladeset>::iterator pop_clade = pop_clades.begin();
                            pop_clade != pop_clades.end(); ++pop_clade){
                            if (pop_clade->first != igc->first){
                                if ((st & pop_clade->second).count() > 0){
                                    // Has members of this pop.
                                    float tmrca = tn->dist_below/(tn->dist_below+tn->dist_above);
                                    // Correct for branch shortening
                                    for (unordered_map<cladeset, float>::iterator bs = brshorten_vals.begin();
                                        bs != brshorten_vals.end(); ++bs){
                                        if ((st & bs->first).count() > 0){
                                            tmrca += (1/(float)st.count()) * bs->second;
                                        }
                                    }
                                    if (tmrca > 1.0){
                                        tmrca = 1.0;
                                    }
                                    
                                    if (haps_ils_count[hapnames_cl[*hap]][pop_clade->first] == 0){
                                        haps_ils_count[hapnames_cl[*hap]].insert(make_pair(pop_clade->first, (long int)0));
                                    }
                                    haps_ils_count[hapnames_cl[*hap]][pop_clade->first]++;
                                    if (haps_ils_tmrca_sum[hapnames_cl[*hap]][pop_clade->first] == 0){
                                        haps_ils_tmrca_sum[hapnames_cl[*hap]].insert(make_pair(pop_clade->first, 0.0));
                                    }
                                    haps_ils_tmrca_sum[hapnames_cl[*hap]][pop_clade->first] += tmrca;
                                    
                                    if (haps_ils_tmrca_hist[hapnames_cl[*hap]].count(pop_clade->first) == 0){
                                        map<float, float> m;
                                        haps_ils_tmrca_hist[hapnames_cl[*hap]].insert(make_pair(pop_clade->first, m));
                                    }
                                    //float tmrca = tn->dist_below/(tn->dist_below+tn->dist_above);
                                    tmrca = round(tmrca/(1/(float)bins))*(1/(float)bins);
                                    if (haps_ils_tmrca_hist[hapnames_cl[*hap]][pop_clade->first].count(tmrca) == 0){
                                        haps_ils_tmrca_hist[hapnames_cl[*hap]][pop_clade->first].insert(make_pair(tmrca, 0.0));
                                    }
                                    haps_ils_tmrca_hist[hapnames_cl[*hap]][pop_clade->first][tmrca]++;
                                }
                            }
                        }
                    }
                }
            }
            // Look for sorted populations
            for (map<string, cladeset>::iterator pc = pop_clades.begin(); pc != pop_clades.end();
                ++pc){
                if (pops_ils_count.count(pc->first) == 0){
                    map<string, long int> m;
                    pops_ils_count.insert(make_pair(pc->first, m));
                    map<string, float> m2;
                    pops_ils_tmrca_sum.insert(make_pair(pc->first, m2));
                    map<string, map<float, float> > m3;
                    pops_ils_tmrca_hist.insert(make_pair(pc->first, m3));
                }
                treeNode* tn = tree->get_clade_match(pc->second);
                if (tn != NULL && tn->parent != NULL){
                    // This population exists as a clade.
                    cladeset st = tn->parent->subtree_leaves();
                    for (map<string, cladeset>::iterator pc2 = pop_clades.begin(); pc2 != pop_clades.end();
                        ++pc2){
                        if (pc2->first != pc->first){
                            if ((st & pc2->second).count() > 0){
                                // Closest relative includes this pop.
                                float tmrca = (tn->parent->dist_below/(tn->parent->dist_below+tn->parent->dist_above));
                                // Correct for branch shortening
                                for (unordered_map<cladeset, float>::iterator bs = brshorten_vals.begin();
                                    bs != brshorten_vals.end(); ++bs){
                                    if ((st & bs->first).count() > 0){
                                        tmrca += (1/(float)st.count()) * bs->second;
                                    }
                                }
                                if (tmrca > 1.0){
                                    tmrca = 1.0;
                                }
                                if (pops_ils_count[pc->first].count(pc2->first) == 0){
                                    pops_ils_count[pc->first].insert(make_pair(pc2->first, (long int)0));
                                }
                                pops_ils_count[pc->first][pc2->first]++;
                                if (pops_ils_tmrca_sum[pc->first].count(pc2->first) == 0){
                                    pops_ils_tmrca_sum[pc->first].insert(make_pair(pc2->first, 0.0));
                                }
                                pops_ils_tmrca_sum[pc->first][pc2->first] += tmrca;
                                
                                if (pops_ils_tmrca_hist[pc->first].count(pc2->first) == 0){
                                    map<float, float> m;
                                    pops_ils_tmrca_hist[pc->first].insert(make_pair(pc2->first, m));
                                }
                                //float tmrca = tn->parent->dist_below/(tn->parent->dist_below+tn->parent->dist_above);
                                tmrca = round(tmrca/(1/(float)bins))*(1/(float)bins);
                                if (pops_ils_tmrca_hist[pc->first][pc2->first].count(tmrca) == 0){
                                    pops_ils_tmrca_hist[pc->first][pc2->first].insert(make_pair(tmrca, 0.0));
                                }
                                pops_ils_tmrca_hist[pc->first][pc2->first][tmrca]++;
                            }
                        }
                    }
                }
            }
        }
        
        delete tree;
    }
    
    fprintf(stderr, "\n");
    
    for (map<string, map<string, long int> >::iterator hic = haps_ils_count.begin(); hic != haps_ils_count.end(); ++hic){
        for (map<string, long int>::iterator hic2 = hic->second.begin(); hic2 != hic->second.end();
            ++hic2){
            if (haps_ils_tmrca_sum.count(hic->first) > 0 && haps_ils_tmrca_sum[hic->first].count(hic2->first) > 0){
                fprintf(outhapf, "%s\t%s\t%ld\t%f\n", hic->first.c_str(), hic2->first.c_str(), 
                    hic2->second, haps_ils_tmrca_sum[hic->first][hic2->first]/(float)hic2->second);
            }
        }
    }
    fclose(outhapf);
    
    for (map<string, map<string, long int> >::iterator pic = pops_ils_count.begin(); pic != pops_ils_count.end(); ++pic){
        for (map<string, long int>::iterator pic2 = pic->second.begin(); pic2 != pic->second.end();
            ++pic2){
            if (pops_ils_tmrca_sum.count(pic->first) > 0 && pops_ils_tmrca_sum[pic->first].count(pic2->first) > 0){
                fprintf(outpopf, "%s\t%s\t%ld\t%f\n", pic->first.c_str(), pic2->first.c_str(), 
                    pic2->second, pops_ils_tmrca_sum[pic->first][pic2->first]/(float)pic2->second);
            }
        }
    }
    
    for (map<string, map<string, map<float, float> > >::iterator hit = haps_ils_tmrca_hist.begin();
        hit != haps_ils_tmrca_hist.end(); ++hit){
        for (map<string, map<float, float> >::iterator hit2 = hit->second.begin(); hit2 != hit->second.end();
            ++hit2){
            for (map<float, float>::iterator hit3 = hit2->second.begin(); hit3 != hit2->second.end(); ++hit3){
                fprintf(outhaphistf, "%s\t%s\t%f\t%f\n", hit->first.c_str(),
                    hit2->first.c_str(), hit3->first, hit3->second);
            }
        }
    }
    
    fclose(outhaphistf);
    
    for (map<string, map<string, map<float, float> > >::iterator pit = pops_ils_tmrca_hist.begin();
        pit != pops_ils_tmrca_hist.end(); ++pit){
        for (map<string, map<float, float> >::iterator pit2 = pit->second.begin();
            pit2 != pit->second.end(); ++pit2){
            for (map<float, float>::iterator pit3 = pit2->second.begin();
                pit3 != pit2->second.end(); ++pit3){
                fprintf(outpophistf, "%s\t%s\t%f\t%f\n", pit->first.c_str(),
                    pit2->first.c_str(), pit3->first, pit3->second);
            }
        }
    } 
    fclose(outpophistf);
    
    
    return 0;
}
