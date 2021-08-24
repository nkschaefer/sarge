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

int main(int argc, char *argv[]) {    
    int bufsize = 1048576;
    
    string msfilename;
    if (argc < 2){
        fprintf(stderr, "USAGE: ms_corr_dist [MSFILE]\n");
        exit(1);
    }
    
    string msfile = string(argv[1]);
    
    // Parse trees from MS run
    map<pair<long int, long int>, treeNode> ms_trees;
    fprintf(stderr, "Parsing trees from MS run...\n");
    int nhaps = parse_ms_file(ms_trees, msfile);
    fprintf(stderr, "Finished\n");
    
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
    map<pair<long int, long int>, treeNode>::iterator ms_tree = ms_trees.begin();
    //if (mask.count() != mask.size()){
    //    ms_tree->second.mask_haps(mask);
    //}
    
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
        while (ms_tree->first.second < pos){
            ms_tree++;
        }
        if (ms_tree == ms_trees.end()){
            break;
        }
        
        if (pos >= ms_tree->first.first && pos <= ms_tree->first.second){
            
            unordered_set<cladeset> wrong_clades;
            tree.get_wrong_clades(&ms_tree->second, wrong_clades);
            
            // Look how far away we need to go to see the correct clade.
            for (unordered_set<cladeset>::iterator wc = wrong_clades.begin(); wc != wrong_clades.end();
                ++wc){
                long int dist_f = -1;
                map<pair<long int, long int>, treeNode>::iterator ms_tree_f(ms_tree);
                ms_tree_f++;
                bool found_f = false;
                if (ms_tree_f == ms_trees.end()){
                    found_f = true;
                }
                while (!found_f){
                    if (ms_tree_f->second.get_clade_match(*wc) != NULL){
                        dist_f = ms_tree_f->first.first - pos;
                        found_f = true;
                    }
                    if (!found_f){
                        ++ms_tree_f;
                        if (ms_tree_f == ms_trees.end()){
                            found_f = true;
                        }
                    }
                }
                
                long int dist_r = -1;
                map<pair<long int, long int>, treeNode>::reverse_iterator ms_tree_r(ms_tree);
                bool found_r = false;
                if (ms_tree_r == ms_trees.rend()){
                    found_r = true;
                }
                while (!found_r){
                    if (ms_tree_r->second.get_clade_match(*wc) != NULL){
                        dist_r = pos - ms_tree_r->first.second;
                        found_r = true;
                    }
                    if (!found_r){
                        if (dist_f != -1 && pos - ms_tree_r->first.second > dist_f){
                            found_r = true;
                        }
                        else{
                            ++ms_tree_r;
                            if (ms_tree_r == ms_trees.rend()){
                                found_r = true;
                            }
                        }
                    }
                }
                if (dist_f  == -1 && dist_r == -1){
                    printf("-1\n");
                }
                else if (dist_f == -1 && dist_r != -1){
                    printf("%ld\n", dist_r);
                }
                else if (dist_r == -1 && dist_f != -1){
                    printf("%ld\n", dist_f);
                }
                else if (dist_r < dist_f){
                    printf("%ld\n", dist_r);
                }
                else{
                    printf("%ld\n", dist_f);
                }
            }
        }
    }
    
    fprintf(stderr, "\n");
    
}
