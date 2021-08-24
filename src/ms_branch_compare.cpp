#include <zlib.h>
#include <math.h>
#include <time.h>
#include <string>
#include <algorithm>
#include <getopt.h>
#include <argp.h>
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

/**
 * Print a help message to the terminal and exit.
 */
void help(int code){
   fprintf(stderr, "ms_branch_compare [OPTIONS]\n");
   fprintf(stderr, "For clades in a file that also exist in the true MS trees, compare branch lengths.\n");
    fprintf(stderr, "[OPTIONS]:\n");
    fprintf(stderr, "   --ms -m The MS output file (REQUIRED)\n");
    fprintf(stderr, "   --unnorm -N do not normalize branch lengths for file being compared to MS\n");
    fprintf(stderr, "   --scale1 -1 Divide branch lengths in test data by this number\n");
    fprintf(stderr, "   --scale2 -2 Divide branch lengths in MS data by this number\n");
exit(code);
}

int main(int argc, char *argv[]) {    
    int bufsize = 1048576;
    
    string msfilename;
    if (argc < 2){
        help(1);
    }
    
    static struct option long_options[] = {
       {"ms", required_argument, 0, 'm'},
       {"unnorm", no_argument, 0, 'N'},
       {"scale1", required_argument, 0, '1'},
       {"scale2", required_argument, 0, '2'},
       {"help", optional_argument, 0, 'h'},
       {0, 0, 0, 0} 
    };
    
    if (argc == 1){
        help(0);
    }
    
    float scale1 = 1.0;
    float scale2 = 1.0;
    bool unnorm = false;
    
    int option_index = 0;
    int ch;
    
    while((ch = getopt_long(argc, argv, "m:1:2:Nh", long_options, &option_index )) != -1){
        switch(ch){
            case 0:
                // This option set a flag. No need to do anything here.
                break;
            case 'm':
                msfilename = optarg;
                break;
            case '1':
                scale1 = atof(optarg);
                break;
            case '2':
                scale2 = atof(optarg);
                break;
            case 'N':
                unnorm = true;
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
    
    if (msfilename.length() == 0){
        fprintf(stderr, "ERROR: ms file required\n");
        exit(1);
    }
    
    std::ifstream msfilep;
    msfilep.open(msfilename.c_str(), std::ifstream::in);
    if (!msfilep.good()){
        fprintf(stderr, "ERROR: could not open MS file %s for reading\n", msfilename.c_str());
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
    
    printf("brlen_true\tbrlen_est\n");
    
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
            vector<pair<treeNode*, treeNode*> > cladepairs;
            tree.get_matching_clades(&mstree, cladepairs);
            
            for (vector<pair<treeNode*, treeNode*> >::iterator cp = cladepairs.begin();
                cp != cladepairs.end(); ++cp){
                // First is MS tree - do not normalize
                float d1 = cp->first->dist;
                float d2 = cp->second->dist_norm;
                if (unnorm){
                    d2 = cp->second->dist;
                }
                // Scale by correct numbers
                d1 = d1 * scale1;
                d2 = d2 * scale2;
                fprintf(stdout, "%f %f\n", d1, d2);
            }
            
        }
    }
    
    fprintf(stderr, "\n");
    
}
