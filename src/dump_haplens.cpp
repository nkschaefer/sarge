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
    long int pos, long int minsize, long int minlen, bool freqs_given,
    map<pair<long double, long double>, vector<long double> >& ages_freqs, float p_cutoff,
    float tmrca_cutoff){
    if (tree->parent != NULL && tree->persistence > 0){
        cladeset st = tree->subtree_leaves();
        if (st.count() >= minsize && tree->persistence >= minlen){
            long double tmrca = tree->dist_below/(tree->dist_below+tree->dist_above);
            if (tree->dist_below + tree->dist_above == 0){
                tmrca = 0;
            }
            if (tmrca <= tmrca_cutoff){
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
                            if (p_cutoff == 1 || probsum < p_cutoff){
                                fprintf(stdout, "%s\t%ld\t%d\t%ld\t%ld\t%Lf\t%f\n", chrom.c_str(), pos,
                                    traversal_index, st.count(), tree->persistence, tmrca, probsum);
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
                    fprintf(stdout, "%s\t%ld\t%d\t%ld\t%ld\t%Lf\n", chrom.c_str(), pos,
                        traversal_index, st.count(), tree->persistence, tmrca);
                }
            }
        }
    }
    traversal_index++;
    for (vector<treeNode*>::iterator child = tree->children.begin(); child !=
        tree->children.end(); ++child){
        print_haplens(*child, traversal_index, chrom, pos, minsize, minlen, freqs_given, 
            ages_freqs, p_cutoff, tmrca_cutoff);
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
    fprintf(stderr, "   --size -s (OPTIONAL) minimum number of haplotypes for a clade \
to contain in order to be printed (must be > 1)\n");
    fprintf(stderr, "   --length -l (OPTIONAL) minimum number of bases for a clade to \
persist in order to be printed (must be >= 1)\n");
    fprintf(stderr, "   --freqs -f (OPTIONAL) to calculate probabilities of clade sizes \
based on TMRCA, provide table computed by utilities/calc_ages_freqs.py\n");
    fprintf(stderr, "   --p_cutoff -p (OPTIONAL) if providing a frequency table (-f argument), \
the maximum p-value allowed to print a clade\n");
    fprintf(stderr, "   --tmrca_cutoff, -t (OPTIONAL) a TMRCA value (as percent of total divergence \
to outgroup) above which not to print clades\n");
    exit(code);
}

int main(int argc, char *argv[]) {
    static struct option long_options[] = {
       {"bufsize", required_argument, 0, 'b'},
       {"help", optional_argument, 0, 'h'},
       {"size", required_argument, 0, 's'},
       {"length", required_argument, 0, 'l'},
       {"freqs", required_argument, 0, 'f'},
       {"p_cutoff", required_argument, 0, 'p'},
       {"tmrca_cutoff", required_argument, 0, 't'},
       {0, 0, 0, 0} 
    };
    
    int bufsize = 1048576;
    long int minsize = 2;
    long int minlen = 1;
    int option_index = 0;
    string freqsfile;
    bool freqs_given = false;
    float p_cutoff = 1.0;
    float tmrca_cutoff = 1.0;
    
    int ch;
    /*
    if (argc == 1){
        help(0);
    }
    */
    while((ch = getopt_long(argc, argv, "b:s:l:f:p:t:h", long_options, &option_index )) != -1){
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
            case 'f':
                freqsfile = optarg;
                freqs_given = true;
                break;
            case 'p':
                p_cutoff = atof(optarg);
                break;
            case 't':
                tmrca_cutoff = atof(optarg);
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
    
    long int last_printed = 0;
    long int progress = 5000;
    
    fprintf(stdout, "#chrom\tpos\ttrav_index\thaps\tpersistence\ttmrca");
    if (freqs_given){
        fprintf(stdout, "\tp");
    }
    fprintf(stdout, "\n");
    
    while(!is.finished()){
 
        string chrom;
        long int pos;
        treeNode* tree = new treeNode();
        tree->set_haps(num_haplotypes);
        read_sitedata(num_haplotypes, is, chrom, pos, *tree);
        int ti = 0;
        print_haplens(tree, ti, chrom, pos, minsize, minlen, freqs_given, ages_freqs, p_cutoff, tmrca_cutoff);
        delete tree;
    }
    
}
