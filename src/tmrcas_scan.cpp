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
   fprintf(stderr, "tmrcas_scan [OPTIONS]\n");
   fprintf(stderr, "Prints out TMRCAs of given groups at individual sites.\n");
    fprintf(stderr, "[OPTIONS]:\n");
    fprintf(stderr, "    --indvs -v The file containing names of individuals \n");
    fprintf(stderr, "   --bufsize -b (OPTIONAL) The number of characters to read from the \
input file at a time\n");
    fprintf(stderr, "    --ingroup -i One or more populations that are admixed (can \
specify more than once)\n");
    fprintf(stderr, "    --indvs -v The file containing names of individuals \n");
    fprintf(stderr, "    --pops -p A file mapping individual names to population names, \
tab-separated\n");
    fprintf(stderr, "   --superpops -s A file mapping population IDs to superpop IDs (optional)\n");
exit(code);
}

float tmrca_group_aux(treeNode* parent, cladeset& group_whole){
    cladeset st = parent->subtree_leaves();
    if (issubset_bitset(st, group_whole)){
        // The dist below this is a mean of branch lengths to members of this group.
        return (parent->dist + parent->dist_below) / (parent->dist_below + parent->dist_above);
    }
    float children_sum = 0;
    float children_count = 0;
    for (vector<treeNode*>::iterator child = parent->children.begin(); 
        child != parent->children.end(); ++child){
        if (((*child)->subtree_leaves() & group_whole).count() > 0){
            children_sum += tmrca_group_aux(*child, group_whole);
            children_count += 1;
        }
    }
    if (children_count == 0){
        children_count = 1;
    }
    return parent->dist_norm + (children_sum/children_count);
}

float tmrca_group(treeNode* tree, cladeset& group_whole){
    treeNode* parent = tree->get_smallest_containing(group_whole);
    if (parent->subtree_leaves() == group_whole){
        // Stop here.
        return parent->dist_below / (parent->dist_below + parent->dist_above);
    }
    else{
        // This is harder
        float children_sum = 0;
        float children_count = 0;
        for (vector<treeNode*>::iterator child = parent->children.begin();
            child != parent->children.end(); ++child){
            if (((*child)->subtree_leaves() & group_whole).count() > 0){
                children_sum += tmrca_group_aux(*child, group_whole);
                children_count += 1;
            }
        }
        if (children_count == 0){
            children_count = 1;
        }
        return children_sum / children_count;
    }
    
}

int main(int argc, char *argv[]) {
    static struct option long_options[] = {
       {"indvs", required_argument, 0, 'v'},
       {"bufsize", required_argument, 0, 'b'},
       {"ingroup", required_argument, 0, 'i'},
       {"indvs", required_argument, 0, 'v'},
       {"pops", required_argument, 0, 'p'},
       {"superpops", required_argument, 0, 's'},
       {"help", optional_argument, 0, 'h'},
       {0, 0, 0, 0} 
    };
    
    vector<string> hapnames;
    int num_haplotypes;
    string indvfilename;
    vector<string> ingroup_pops;
    string popfilename;
    string superpopfilename;
    bool superpops_given = false;
    bool tmrcafile_given = false;
    
    int ploidy = 2;
    
    int bufsize = 1048576;
    
    int option_index = 0;
    int ch;
    
    if (argc == 1){
        help(0);
    }
    while((ch = getopt_long(argc, argv, "v:b:i:v:p:s:h", long_options, &option_index )) != -1){
        switch(ch){
            case 0:
                // This option set a flag. No need to do anything here.
                break;
            case 'v':                
                indvfilename = string(optarg);
                break;
            case 'b':
                bufsize = atoi(optarg);
                break;
            case 'i':
                ingroup_pops.push_back(optarg);
                break;
            case 'p':
                popfilename = optarg;
                break;
            case 's':
                superpopfilename = optarg;
                superpops_given = true;
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

    if (indvfilename.length() == 0){
        fprintf(stderr, "ERROR: you must provide a file listing haplotype names.\n");
        exit(1);
    }
    parse_indvs(hapnames, indvfilename);
    
    // Get file reading stuff ready
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
    cladeset mask;
    read_header(is, num_haplotypes, mask);
    
    fprintf(stderr, "Read %d haplotypes\n", num_haplotypes);
    
    map<string, int> indv2hap;
    map<int, string> hap2pop;
    parse_indvs_map(indv2hap, indvfilename);
    
    map<string, string> indv2pop;
    map<string, vector<string> > pop2indv;
    parse_pops(indv2pop, pop2indv, popfilename, ploidy);
    
    map<string, string> pop2superpop;
    map<string, vector<string> > superpop2pop;
    if (superpops_given){
        parse_pops(pop2superpop, superpop2pop, superpopfilename, 0);
    }
    
    // Transform populations into clades
    map<string, cladeset> ingroup_clades;
    for (vector<string>::iterator ingroup_pop = ingroup_pops.begin();
        ingroup_pop != ingroup_pops.end(); ++ingroup_pop){
        cladeset popclade;
        for (vector<string>::iterator indvname = pop2indv[*ingroup_pop].begin();
            indvname != pop2indv[*ingroup_pop].end(); ++indvname){
            if (indv2hap.count(*indvname) > 0){
                popclade.set(num_haplotypes-1-indv2hap[*indvname]);
            }
        }  
        ingroup_clades.insert(make_pair(*ingroup_pop, popclade));
    }
    
    // Transform superpopulations into clades
    for (map<string, vector<string> >::iterator sp = superpop2pop.begin(); sp != superpop2pop.end(); ++sp){
        cladeset spclade;
        for (vector<string>::iterator pop = sp->second.begin(); pop != sp->second.end();
            ++pop){
            if (ingroup_clades.count(*pop) > 0){
                spclade |= ingroup_clades[*pop];
            }
        }
        ingroup_clades.insert(make_pair(sp->first, spclade));
    }
    
    int denom = 0;
    int num = 0;
    
    int progress = 10000;
    int last_printed = 0;
    
    long int run_start = -1;
    long int run_end = -1;
    
    string prevchrom = "";
    
    while(!is.finished()){
 
        string chrom;
        long int pos;
        treeNode tree;
        read_sitedata(num_haplotypes, is, chrom, pos, tree);

        if (pos > last_printed + progress){
            fprintf(stderr, "Processed %s\t%ld\r", chrom.c_str(), pos);
            last_printed = pos;      
        }
        
        fprintf(stdout, "%s\t%ld\t%f\tROOT\n", chrom.c_str(), pos, tree.dist_below / (tree.dist_below + tree.dist_above));
        
        for (map<string, cladeset>::iterator pop_clade = ingroup_clades.begin();
            pop_clade != ingroup_clades.end(); ++pop_clade){
            float tmrca_clade = tmrca_group(&tree, pop_clade->second);
            fprintf(stdout, "%s\t%ld\t%f\t%s\n", chrom.c_str(), pos, tmrca_clade,
                pop_clade->first.c_str());
        }
        
        prevchrom = chrom;
        
    }
    
}
