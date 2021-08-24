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
#include <cmath>
#include <zlib.h>
#include "sets.h"
#include "treeNode.h"
#include "serialize.h"
#include "common.h"

using std::cout;
using std::endl;
using namespace std;

cladeset kc_aux(treeNode* t, 
    int dist_down, 
    float len_down,
    map<pair<int, int>, int>& kc_top,
    map<pair<int, int>, float>& kc_br,
    int num_haplotypes){
    cladeset st = t->subtree_leaves();
    for (vector<treeNode*>::iterator child = t->children.begin(); child != t->children.end();
        ++child){
        // Ignore singleton nodes
        if (((*child)->subtree_leaves()).count() > 1){
            cladeset cst = kc_aux(*child, dist_down + 1, 
                len_down + (*child)->dist_norm, kc_top, kc_br, num_haplotypes);
            //st = set_diff_bitset(st, cst);
        }
    }
    set<unsigned int> remaining = bitset2set(st, num_haplotypes);
    for (set<unsigned int>::iterator hap1 = remaining.begin(); hap1 != remaining.end(); ++hap1){
        set<unsigned int>::iterator hap2 = hap1;
        hap2++;
        while(hap2 != remaining.end()){
            pair<int, int> key = make_pair(*hap1, *hap2);
            if (kc_top.count(key) == 0){
                kc_top.insert(make_pair(key, dist_down));
                kc_br.insert(make_pair(key, len_down));
            }
            ++hap2; 
        }
    }
    //return t->subtree_leaves();
    return st;
}

/** 
 * Computes Kendall-Cojin vector (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5026250/)
 * for a provided tree.
 */
void kendall_cojin(treeNode* t, map<pair<int, int>, int>& kc_top, 
    map<pair<int, int>, float>& kc_br,
    int num_haplotypes){
    kc_aux(t, 0, 0.0, kc_top, kc_br, num_haplotypes);
}

void combine_vecs(map<pair<int, int>, int>& kc_top,
    map<pair<int, int>, float>& kc_br,
    map<pair<int, int>, float>& result,
    float lambda){
    
    for (map<pair<int, int>, int>::iterator m = kc_top.begin();
        m != kc_top.end(); ++m){
        result.insert(make_pair(m->first, (1-lambda)*m->second + lambda*kc_br[m->first]));
    }
}

/**
 * lambda = percentage that branch lengths are allowed to contribute to 
 * metric
 */
float kc_vec_compare(map<pair<int, int>, int>& kc_top1,
    map<pair<int, int>, float>& kc_br1,
    map<pair<int, int>, int>& kc_top2,
    map<pair<int, int>, float>& kc_br2,
    float lambda){
    
    map<pair<int, int>, float> v_lambda1;
    combine_vecs(kc_top1, kc_br1, v_lambda1, lambda);
    map<pair<int, int>, float> v_lambda2;
    combine_vecs(kc_top2, kc_br2, v_lambda2, lambda);
    
    float sum = 0;
    for (map<pair<int, int>, float>::iterator item1 = v_lambda1.begin();
        item1 != v_lambda1.end(); ++item1){
        sum += (pow(item1->second - v_lambda2[item1->first], 2));
    }
    return sqrt(sum);
    
}

int main(int argc, char *argv[]) {    
    int bufsize = 1048576;
    
    string msfilename;
    if (argc < 2){
        fprintf(stderr, "USAGE: score_trees_ms [MSFILE] <lambda> <mean>\n");
        fprintf(stderr, "<lambda> is optional parameter - a float between 0 and 1 \
that determines the extent to which branch lengths contribute to Kendall-Cojin distance\n");
        fprintf(stderr, "<mean> is mean (OPTIONAL), to summarize output instead of \
printing stats for every tree\n");
        exit(1);
    }
    
    float lambda = 0.0;
    if (argc > 2){
        lambda = atof(argv[2]);
    }
    
    string msfile = string(argv[1]);
    
    std::ifstream msfilep;
    msfilep.open(msfile.c_str(), std::ifstream::in);
    if (!msfilep.good()){
        fprintf(stderr, "ERROR: could not open MS file %s for reading\n", msfile.c_str());
        exit(1);
    }
    
    // Determine whether to just average stuff or output for every site
    bool mean = false;
    if (argc > 3 && strcmp(argv[3], "mean") == 0){
        mean = true;
    }
    
    // Parse trees from MS run
    //map<pair<long int, long int>, treeNode> ms_trees;
    //fprintf(stderr, "Parsing trees from MS run...\n");
    //int nhaps = parse_ms_file(ms_trees, msfile);
    //fprintf(stderr, "Finished\n");
    
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
    //if (mask.size() == 0){
    //    mask.set();
    //}
    //map<pair<long int, long int>, treeNode>::iterator ms_tree = ms_trees.begin();
    //if (mask.count() != mask.size()){
    //    ms_tree->second.mask_haps(mask);
    //}
    map<pair<int, int>, int> ms_kc_top;
    map<pair<int, int>, float> ms_kc_br;
    //kendall_cojin(&ms_tree->second, ms_kc_top, ms_kc_br, num_haplotypes);
            
    map<pair<int, int>, int> kc_top;
    map<pair<int, int>, float> kc_br;
    
    float nodes_corr_sum = 0;
    float nodes_tot_sum = 0;
    float msnodes_corr_sum = 0;
    float msnodes_tot_sum = 0;
    float sites_count = 0;
    float kc_sum = 0;
    float kc_sum2 = 0;
    
    printf("nodes_correct\tnodes_tot\tnodes_correct_ms\tnodes_tot_ms\tk-c\tk-cl0.5\n");
    
    long int msfilestart = -1;
    long int msfileend = -1;
    treeNode mstree;
    
    long int last_printed = 1;
    long int progress = 1000;
    
    // Read from stdin.
    while(!is.finished()){

        string chrom;
        long int pos;
        treeNode tree;
        tree.delete_children();
        read_sitedata(num_haplotypes, is, chrom, pos, tree);
        
        if (pos - last_printed >= progress){
            fprintf(stderr, "Read site %ld\r", pos);
            last_printed = pos;
        }
        // Catch up MS file.
        bool newtree = false;
        while (msfileend < pos){
            mstree.delete_children();
            bool success = ms_next_tree(msfilep, mstree, msfilestart, msfileend);
            if (!success){
                break;
            }
            newtree = true;
        }
        
        if (pos >= msfilestart && pos <= msfileend){
            vector<string> x;
        
            if (newtree){
                ms_kc_top.clear();
                ms_kc_br.clear();
                kendall_cojin(&mstree, ms_kc_top, ms_kc_br, num_haplotypes);    
            }
        
            kc_top.clear();
            kc_br.clear();
            kendall_cojin(&tree, kc_top, kc_br, num_haplotypes);
        
            float kcdist = kc_vec_compare(kc_top, kc_br, ms_kc_top, ms_kc_br, lambda);
            float kcdist2 = kc_vec_compare(kc_top, kc_br, ms_kc_top, ms_kc_br, 0.5);
        
            kc_sum += kcdist;
            kc_sum2 += kcdist2;

            int nodes_correct = tree.clades_correct(&mstree);
            int nodes_tot = tree.clades_tot();
        
            int msnodes_correct = mstree.clades_correct(&tree);
            int msnodes_tot = mstree.clades_tot();
        
            if (nodes_tot == 0){
                nodes_tot = 1;
            }
            if (msnodes_tot == 0){
                msnodes_tot = 1;
            }
            
            if (!mean){
                printf("%f\t%d\t%f\t%d\t%f\t%f\n",
                    (float)nodes_correct/(float)nodes_tot, 
                    nodes_tot, 
                    (float)msnodes_correct/(float)msnodes_tot,
                    msnodes_tot,
                    kcdist,
                    kcdist2);
            }
            else{
                nodes_corr_sum += ((float)nodes_correct/(float)nodes_tot);
                nodes_tot_sum += nodes_tot;
            
                msnodes_corr_sum += ((float)msnodes_correct/(float)msnodes_tot);
                msnodes_tot_sum += msnodes_tot;
                sites_count += 1;
            }
        }
    }
    
    fprintf(stderr, "\n");
    
    if (mean){
        if (sites_count == 0){
            sites_count = 1;
        }
        sites_count = (float) sites_count;
        printf("%f\t%f\t%f\t%f\t%f\t%f\n", 
            nodes_corr_sum/sites_count,
            nodes_tot_sum/sites_count,
            msnodes_corr_sum/sites_count,
            msnodes_tot_sum/sites_count,
            kc_sum/sites_count,
            kc_sum2/sites_count);
    }
}
