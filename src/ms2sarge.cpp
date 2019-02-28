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
#include <regex>
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
    if (argc < 4){
        fprintf(stderr, "Usage:\n");
        fprintf(stderr, "ms2sarge <msfile> <sitesfile> <outfile>\n");
        exit(1);
    }
    
    string msfile = string(argv[1]);
    string sitesfile = string(argv[2]);
    string outfile = string(argv[3]);
    
    // Parse sites
    fprintf(stderr, "Parsing sites file\n");
    map<string, set<pair<long int, long int> > > excl_bed;
    vector<loc> sites = parse_locs(sitesfile, excl_bed);
    fprintf(stderr, "Finished\n");
    
    // Open output file
    gzFile fp_out = gzopen(outfile.c_str(), "w");
    if (!fp_out){
        fprintf(stderr, "ERROR: could not open file %s for writing.\n", outfile.c_str());
        exit(1);
    }
    
    long int progress = 10000;
    long int prevpos = 0;
    
    vector<loc>::iterator cursite = sites.begin();
    
    std::ifstream infile;
    infile.open(msfile.c_str(), std::ifstream::in);
    if (!infile.good()){
        fprintf(stderr, "ERROR: could not open file %s for reading\n", msfile.c_str());
        exit(1);
    }
    
    regex treematch(R"(^\[([0-9]+)](.+)$)", regex::extended);
    
    long int start = 1;
    
    int num_haplotypes = -1;
    
    bool reading_trees = false;
    
    for (string line; getline(infile, line);){
        // Need to get the number of haplotypes -- this can be found in 
        // the first line, which is a copy of the command used to run MS
        if (num_haplotypes == -1){
            istringstream space_splitter(line);
            string field;
            int field_index = 0;
            while(std::getline(space_splitter, field, ' ')){
                if (field_index == 0){
                    // This is the path to MS.
                }
                else if (field_index == 1){
                    num_haplotypes = atoi(field.c_str());
                    break;
                }
                field_index++;
            }
            if (num_haplotypes == -1){
                fprintf(stderr, "ERROR: unable to read the number of haplotypes from the MS output.\n");
                exit(1);
            }
            else{
                fprintf(stderr, "Read %d haplotypes\n", num_haplotypes);
            }
            cladeset mask;
            serialize_header(fp_out, num_haplotypes, mask);
        }
        else{
            // Check that this is a tree line and extract useful information
            cmatch match;

            if (regex_match(line.c_str(), match, treematch)){
                reading_trees = true;
            
                string len = extract_regex_match(line.c_str(), match[0], match[1]);
                string newick = extract_regex_match(line.c_str(), match[0], match[2]);
                
                long int len_num = strtol(len.c_str(), 0, 10);
            
                long int end = start + len_num-1;
            
                treeNode root;
                
                string newick_decremented = decrement_haps_newick(newick);

                string newick_finished = replace_scinot_newick(newick_decremented);
                
                parse_newick(root, newick_finished, num_haplotypes);
            
                while(cursite->pos < start){
                    cursite++;
                    if (cursite == sites.end()){
                        fprintf(stderr, "\nFound end of sites after loading tree at position (%ld %ld)\n", start, end);
                        gzclose(fp_out);
                        return 0;
                    }
                }
                
                if (cursite->pos >= start && cursite->pos <= end){
                    serialize_sitedata(fp_out, cursite->chrom, cursite->pos, root);
                }
                
                if (cursite->pos - prevpos > progress){
                    fprintf(stderr, "Processed %s %ld\r", cursite->chrom.c_str(), cursite->pos);
                    prevpos = cursite->pos;
                }
                
                start += len_num;
            }
            else if (reading_trees){
            
                // We've read past the tree part of the file.
                fprintf(stderr, "\nFinished reading trees at position %ld\n", start);
                gzclose(fp_out);
                return 0;
            }
        }
    }
    fprintf(stderr, "\n");
    gzclose(fp_out);
    return 0;
}
