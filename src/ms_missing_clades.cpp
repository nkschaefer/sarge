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

/**
 * When an MS clade is NOT identified by SARGE, outputs the sizes of all such clades. 
 * Also outputs the sizes of ALL MS clades (as background).
 */
 
using std::cout;
using std::endl;
using namespace std;

int main(int argc, char *argv[]) {    
    int bufsize = 1048576;
    
    string msfilename;
    if (argc < 2){
        fprintf(stderr, "USAGE: ms_missing_clades [MSFILE]\n");
        exit(1);
    }
    
    string msfile = string(argv[1]);
    
    std::ifstream msfilep;
    msfilep.open(msfile.c_str(), std::ifstream::in);
    if (!msfilep.good()){
        fprintf(stderr, "ERROR: could not open MS file %s for reading\n", msfile.c_str());
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
    //if (mask.size() == 0){
    //    mask.set();
    //}
    //map<pair<long int, long int>, treeNode>::iterator ms_tree = ms_trees.begin();
    //if (mask.count() != mask.size()){
    //    ms_tree->second.mask_haps(mask);
    //}
    
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
            mstree.print_missing_sizes(&tree);
        }
    }
    
    fprintf(stderr, "\n");
    
}
