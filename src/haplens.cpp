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

/**
 * Print a help message to the terminal and exit.
 */
void help(int code){
   fprintf(stderr, "haplens [OPTIONS]\n");
   fprintf(stderr, "Labels each clade of each tree with the physical distance along the chromosome \
for which that clade persists. Outputs (gzipped) SARGE format, with the persistence \
field of each clade populated. This can then be accessed by other programs.\n");
    fprintf(stderr, "[OPTIONS]:\n");
    fprintf(stderr, "   --bufsize -b (OPTIONAL) The number of characters to read from the \
input file at a time\n");
    fprintf(stderr, "   --out -o The name of the output file to create. Output will be \
in gzipped binary format.\n");
    fprintf(stderr, "   --prop_dist -p (OPTIONAL) The length each mutation was allowed to persist \
when the ARG was inferred. If two sites are more than 2* this distance apart, clades \
will be not be inferred to span across the distance between.\n");
    exit(code);
}

int main(int argc, char *argv[]) {
    static struct option long_options[] = {
       {"bufsize", required_argument, 0, 'b'},
       {"out", required_argument, 0, 'o'},
       {"prop_dist", required_argument, 0, 'p'},
       {"help", optional_argument, 0, 'h'},
       {0, 0, 0, 0} 
    };
    
    int bufsize = 1048576;
    string outname;
    int prop_dist = -1;
    
    int option_index = 0;
    int ch;
    
    if (argc == 1){
        help(0);
    }
    while((ch = getopt_long(argc, argv, "b:o:p:h", long_options, &option_index )) != -1){
        switch(ch){
            case 0:
                // This option set a flag. No need to do anything here.
                break;
            case 'b':
                bufsize = atoi(optarg);
                break;
            case 'o':
                outname = optarg;
                break;
            case 'p':
                prop_dist = atoi(optarg);
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
    
    if (outname.length() == 0){
        fprintf(stderr, "ERROR: you must provide an output file name.\n");
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
    
    gzFile outp = NULL;
    if (outname == "-"){
        outp = gzdopen(fileno(stdout), "wb");
        if (outp == NULL){
            fprintf(stderr, "ERROR: unable to open stdout for writing.\n");
            exit(1);
        }
    }
    else{
        outp = gzopen(outname.c_str(), "wb");
        if (outp == NULL){
            fprintf(stderr, "ERROR: unable to open file %s for writing.\n", outname.c_str());
            exit(1);
        }
    }
    
    instream_info is;
    instream_init(is, &fp, bufsize);
    
    // Read file header.
    int num_haplotypes;
    cladeset mask;
    read_header(is, num_haplotypes, mask);
    
    fprintf(stderr, "Read %d haplotypes\n", num_haplotypes);
    
    // Write file header.
    serialize_header(outp, num_haplotypes, mask);
    
    // Store all already-seen trees that haven't been written to disk yet
    map<long int, treeNode*> prevtrees;
    
    // Map trees to sets of clades within them that still exist in the most-recently-read
    // tree
    multimap<treeNode*, cladeset> tree_clades;
    
    string prevchrom = "";
    long int prevpos = -1;
    
    long int last_printed = 0;
    long int progress = 5000;
    
    while(!is.finished()){
 
        string chrom;
        long int pos;
        treeNode* tree = new treeNode();
        tree->set_haps(num_haplotypes);
        read_sitedata(num_haplotypes, is, chrom, pos, *tree);
        
        unordered_set<cladeset> curtree_flat;
        tree->flatten_set(curtree_flat);
        
        if (prop_dist != -1 && pos - prevpos > 2*prop_dist){
            // Too much space in between sites. We can't really infer clades
            // to span this distance, since it might be across a centromere,
            // etc.
            for (map<long int, treeNode*>::iterator pt = prevtrees.begin();
                pt != prevtrees.end(); ++pt){
                // Delete all stored clades tied to this tree.
                for (multimap<treeNode*, cladeset>::iterator pc = 
                    tree_clades.equal_range(pt->second).first; 
                    pc != tree_clades.equal_range(pt->second).second;){
                    // End clades at previous site (pre-gap)
                    long int dist = prevpos - pt->first + 1;
                    treeNode* clade = pt->second->get_clade_match(pc->second);
                    if (clade == NULL){
                        fprintf(stderr, "ERROR: clade not found in tree\n");
                        print_bitset_set(pc->second, num_haplotypes);
                        exit(1);
                    }
                    else{
                        clade->persistence = dist;
                    }
                    tree_clades.erase(pc++);
                }
            }
            for (map<long int, treeNode*>::iterator pt = prevtrees.begin();
                pt != prevtrees.end();){
                serialize_sitedata(outp, chrom, pt->first, *pt->second);
                delete pt->second;
                prevtrees.erase(pt++);
            }
        }
        else{
            for (map<long int, treeNode*>::iterator pt = prevtrees.begin(); 
                pt != prevtrees.end(); ++pt){
                for (multimap<treeNode*, cladeset>::iterator pc = 
                    tree_clades.equal_range(pt->second).first; pc != 
                    tree_clades.equal_range(pt->second).second;){
                    if (curtree_flat.find(pc->second) == curtree_flat.end()){
                        // Not found in current tree.
                        // Determine distance that clade exists.
                        long int dist = prevpos - pt->first + 1;
                        treeNode* clade = pt->second->get_clade_match(pc->second);
                        if (clade == NULL){
                            fprintf(stderr, "ERROR: clade not found in tree\n");
                            print_bitset_set(pc->second, num_haplotypes);
                            //print_node_lite(pt->second, num_haplotypes);
                            exit(1);
                        }
                        else{
                            clade->persistence = dist;
                        }
                        tree_clades.erase(pc++);
                    }
                    else{
                        // Found in current tree.
                        ++pc;
                    }
                }
            }
            // Erase and print out anything already fully processed.
            for (map<long int, treeNode*>::iterator pt = prevtrees.begin();
                pt != prevtrees.end();){
                if (tree_clades.find(pt->second) == tree_clades.end()){
                    // No clades still stored for this tree. We're done with it.
                    serialize_sitedata(outp, chrom, pt->first, *pt->second);
                    delete pt->second;
                    prevtrees.erase(pt++);
                }
                else{
                    // Not done with this tree, so don't try to print anything that
                    // comes after it yet.
                    break;
                }
            }
        }
        
        // Store stuff from current tree in the data structure.
        for (unordered_set<cladeset>::iterator cl = curtree_flat.begin();
            cl != curtree_flat.end();){
            tree_clades.insert(make_pair(tree, *cl));
            curtree_flat.erase(cl++);
        }
        prevtrees.insert(make_pair(pos, tree));
        if (pos > last_printed + progress){
            fprintf(stderr, "Processed %s\t%ld\r", chrom.c_str(), pos);
            last_printed = pos;      
        }
        prevchrom = chrom;
        prevpos = pos;
    }
    
    // Finalize all remaining ranges.
    for (map<long int, treeNode*>::iterator pt = prevtrees.begin(); pt != 
        prevtrees.end(); ++pt){
        for (multimap<treeNode*, cladeset>::iterator pc = 
            tree_clades.equal_range(pt->second).first; pc != 
            tree_clades.equal_range(pt->second).second;){
            long int dist = prevpos - pt->first + 1;
            treeNode* clade = pt->second->get_clade_match(pc->second);
            if (clade == NULL){
                fprintf(stderr, "ERROR: clade not found in tree\n");
                print_bitset_set(pc->second, num_haplotypes);
                //print_node_lite(pt->second, num_haplotypes);
                exit(1);
            }
            else{
                clade->persistence = dist;
            }
            tree_clades.erase(pc++);
        }
    }
    // Print out last trees.
    for (map<long int, treeNode*>::iterator pt = prevtrees.begin();
        pt != prevtrees.end();){
        serialize_sitedata(outp, prevchrom, pt->first, *pt->second);
        delete pt->second;
        prevtrees.erase(pt++);
    }
    
    gzclose(outp);
    
}
