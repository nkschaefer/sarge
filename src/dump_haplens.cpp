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

void print_haplens(treeNode* tree, int& traversal_index, string& chrom, 
    long int pos, long int minsize, long int minlen){
    if (tree->parent != NULL && tree->persistence > 0){
        cladeset st = tree->subtree_leaves();
        if (st.count() >= minsize && tree->persistence >= minlen){
            fprintf(stdout, "%s\t%ld\t%d\t%ld\t%ld\n", chrom.c_str(), pos,
                traversal_index, st.count(), tree->persistence);
        }
    }
    traversal_index++;
    for (vector<treeNode*>::iterator child = tree->children.begin(); child !=
        tree->children.end(); ++child){
        print_haplens(*child, traversal_index, chrom, pos, minsize, minlen);
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
    fprintf(stderr, "   --size -s (OPTIONAL) minimum number of haplotypes for a clade \
to contain in order to be printed (must be > 1)\n");
    fprintf(stderr, "   --length -l (OPTIONAL) minimum number of bases for a clade to \
persist in order to be printed (must be >= 1)\n");
    exit(code);
}

int main(int argc, char *argv[]) {
    static struct option long_options[] = {
       {"bufsize", required_argument, 0, 'b'},
       {"help", optional_argument, 0, 'h'},
       {"size", required_argument, 0, 's'},
       {"length", required_argument, 0, 'l'},
       {0, 0, 0, 0} 
    };
    
    int bufsize = 1048576;
    long int minsize = 2;
    long int minlen = 1;
    int option_index = 0;
    int ch;
    /*
    if (argc == 1){
        help(0);
    }
    */
    while((ch = getopt_long(argc, argv, "b:s:l:h", long_options, &option_index )) != -1){
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
            case 'l':
                minlen = atol(optarg);
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
    
    instream_info is;
    instream_init(is, &fp, bufsize);
    
    // Read file header.
    int num_haplotypes;
    cladeset mask;
    read_header(is, num_haplotypes, mask);
    
    fprintf(stderr, "Read %d haplotypes\n", num_haplotypes);
    
    long int last_printed = 0;
    long int progress = 5000;
    
    while(!is.finished()){
 
        string chrom;
        long int pos;
        treeNode* tree = new treeNode();
        tree->set_haps(num_haplotypes);
        read_sitedata(num_haplotypes, is, chrom, pos, *tree);
        int ti = 0;
        print_haplens(tree, ti, chrom, pos, minsize, minlen);
        delete tree;
    }
    
}
