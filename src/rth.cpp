/**
 * Given ARG output, computes the "relative TMRCA halflife" (RTH or TMRCA50)
 * statistic used by the ArgWeaver paper (Rasmussen et al 2014).
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
   fprintf(stderr, "rth [OPTIONS]\n");
   fprintf(stderr, "Prints the RTH statistic for every tree encountered \
in SARGE output (pipe from stdin).\n");
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
       {0, 0, 0, 0} 
    };
    
    // Set default values
    int bufsize = 1048576;
    
    int option_index = 0;
    int ch;
    
    while((ch = getopt_long(argc, argv, "h", long_options, &option_index )) != -1){
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
            default:
                help(0);
        }    
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
        
        map<float, vector<treeNode*> > dist_map;
        tree->to_dist_map(0, dist_map);
        
        set<treeNode*> visited;
        treeNode* min_node = NULL;
        float min_tmrca = 2.0;
        
        // In this data structure, leaves should come early.
        for (map<float, vector<treeNode*> >::iterator dm = dist_map.begin();
            dm != dist_map.end(); ++dm){
            if (dm->first > min_tmrca){
                break;
            }
            for (vector<treeNode*>::iterator n = dm->second.begin();
                n != dm->second.end(); ++n){
                treeNode* curnode = *n;
                bool stop = false;
                while (!stop){
                    if ((float)curnode->subtree_leaves().count() >= 0.5*(float)num_haplotypes){
                        // This node is a candidate.
                        if ((curnode->dist_below / (curnode->dist_below + curnode->dist_above)) < min_tmrca){
                            min_tmrca = curnode->dist_below / (curnode->dist_below + curnode->dist_above);
                            min_node = curnode;
                        }
                        stop = true;
                    }
                    else{
                        if (curnode->parent == NULL){
                            stop = true;
                        }
                        else if (visited.find(curnode->parent) != visited.end()){
                            stop = true;
                        }
                        else{
                            curnode = curnode->parent;   
                        }
                    }
                }
            }
        }
        
        fprintf(stdout, "%s\t%ld\t%f\n", chrom.c_str(), pos, min_tmrca);
        delete tree;
    }
    
    fprintf(stderr, "\n");
    
    return 0;
}
