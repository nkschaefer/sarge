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

using std::cout;
using std::endl;
using namespace std;

/**
 * Print a help message to the terminal and exit.
 */
void help(int code){
   fprintf(stderr, "newick2sarge [OPTIONS]\n");
   fprintf(stderr, "Given tab-separated newick format data (chrom\tpos\tnewick) - i.e. from utilities/sarge_rentplus.py, \
outputs in SARGE format to the given file.\n");
    fprintf(stderr, "[OPTIONS]:\n");
    fprintf(stderr, "   --input the input .geno.gz file for computing number of haplotypes\n");
    fprintf(stderr, "   --outfile -o The file to create in SARGE format.\n");
    fprintf(stderr, "   --no_branchlengths -b Specify if there are no branch lengths \
in the trees\n");
    fprintf(stderr, "   --decrement -d Specify to reduce each haplotype index by 1\n");
exit(code);
}


int main(int argc, char *argv[]) {
    static struct option long_options[] = {
       {"input", required_argument, 0, 'i'},
       {"outfile", required_argument, 0, 'o'},
       {"no_branchlengths", no_argument, 0, 'b'},
       {"decrement", no_argument, 0, 'd'},
       {"help", optional_argument, 0, 'h'},
       {0, 0, 0, 0} 
    };
    
    string infile;
    string outfile;
    
    bool no_brlen = false;
    bool decrement = false;
    
    int bufsize = 1048576;
 
    int option_index = 0;
    int ch;
    
    if (argc == 1){
        help(0);
    }
    while((ch = getopt_long(argc, argv, "i:o:bdh", long_options, &option_index )) != -1){
        switch(ch){
            case 0:
                // This option set a flag. No need to do anything here.
                break;
            case 'i':
                infile = optarg;
                break;
            case 'o':
                outfile = optarg;
                break;
            case 'b':
                no_brlen = true;
                break;
            case 'd':
                decrement = true;
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
    
    if (infile.length() == 0 || outfile.length() == 0){
        fprintf(stderr, "Missing one ore more required arguments.\n");
        exit(1);
    }
    
    gzFile in_gz = gzopen(infile.c_str(), "r");
    if (!in_gz){
        fprintf(stderr, "Unable to open input file %s\n", infile.c_str());
        exit(1);
    }
    int num_haplotypes = get_num_haplotypes(in_gz, bufsize);
    gzclose(in_gz);
    
    gzFile out_gz = gzopen(outfile.c_str(), "w");
    if (!out_gz){
        fprintf(stderr, "Unable to open output file %s\n", outfile.c_str());
        exit(1);
    }
    
    cladeset mask;
    for (int i = 0; i < num_haplotypes; ++i){
        mask.set(num_haplotypes-1-i);
    }
    serialize_header(out_gz, num_haplotypes, mask);
    
    string line;
    while(getline(cin, line)){
        istringstream splitter(line);
        string field;
        short fld_index = 0;
        
        string chrom;
        long int pos;
        string treedat;
        
        while(getline(splitter, field, '\t')){
            if (fld_index == 0){
                chrom = field;
            }
            else if (fld_index == 1){
                pos = atol(field.c_str());
            }
            else if (fld_index == 2){
                treedat = field;
            }
            fld_index++;
        }
        
        
        treeNode tree;
        if (no_brlen){
            treedat = add_brlens(treedat);
        }
        if (decrement){
            treedat = decrement_haps_newick(treedat);
        }

        parse_newick(tree, treedat, num_haplotypes);
        
        // Only print data about the chromosome if we haven't already seen it.
        // This will save a ton of space.
        static string prevchrom = "";
    
        if (chrom.compare(prevchrom) != 0){
            // Need to store string data.
            serialize_bool(out_gz, true);
            serialize_str(out_gz, chrom);
        }
        else{
            // No need to store string data.
            serialize_bool(out_gz, false);
        }
        prevchrom = chrom;
    
        serialize_longint(out_gz, pos);
        
        // Parse and serialize tree
        tree.serialize(out_gz);
    }
    
    gzclose(out_gz);
    
}
