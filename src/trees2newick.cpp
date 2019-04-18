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

void print_persistence(treeNode* tree, bool use_hapnames, 
    vector<string>& hapnames, vector<string>& text){
    if (tree->persistence > 0 && tree->parent != NULL){
        string clade_text = "";
        set<unsigned int> leaves_set = bitset2set(tree->subtree_leaves(), tree->num_haps);
        int i = 0; 
        for (set<unsigned int>::iterator li = leaves_set.begin(); li != leaves_set.end();
            ++li){
            if (use_hapnames){
                clade_text += hapnames[*li];
            }
            else{
                char intstr[10];
                sprintf(intstr, "%d", *li);
                string intstr2 = intstr;
                clade_text += intstr2;
            }
            if (i < leaves_set.size()-1){
                clade_text += ",";
            }
            ++i;
        }
        clade_text += ":";
        char pers_str[11];
        sprintf(pers_str, "%ld", tree->persistence);
        string pers_str2 = pers_str;
        clade_text += pers_str2;
        text.push_back(clade_text);
    }
    for (vector<treeNode*>::iterator child = tree->children.begin();
        child != tree->children.end(); ++child){
        print_persistence(*child, use_hapnames, hapnames, text);
    }
}

void print_site(treeNode& tree, string chrom, long int pos, bool use_hapnames, 
    vector<string>& hapnames, bool recomb, bool metadata){
    //float treedepth = tree.branchlens2percent();
    //float treedepth = 0;
    fprintf(stdout, "%s\t%ld\t%s", chrom.c_str(), pos, tree.newick(recomb, use_hapnames, hapnames).c_str());    
    if (metadata){
        // Print leaf names, separated by commas, then colon and persistence time
        vector<string> meta_text;
        print_persistence(&tree, use_hapnames, hapnames, meta_text);
        int i = 0;
        fprintf(stdout, "\t");
        for (vector<string>::iterator mt = meta_text.begin(); mt != meta_text.end(); ++mt){
            fprintf(stdout, "%s", mt->c_str());
            if (i < meta_text.size()-1){
                fprintf(stdout, ";");
            }
            ++i;
        } 
    }
    fprintf(stdout, "\n");
}

void divide_branches(treeNode* node, float rootdist){
    if (rootdist == 0.0){
        rootdist = 1.0;
    }
    if (rootdist == 1.0){
        return;
    }
    if (node->parent != NULL){
        node->dist /= rootdist;
    }
    for (vector<treeNode*>::iterator child = node->children.begin(); child != node->children.end(); ++child){
        divide_branches(*child, rootdist);
    }
}

void add_pseudocount(treeNode* node, float ps){
    node->dist_norm += ps;
    for (vector<treeNode*>::iterator child = node->children.begin(); child != node->children.end(); ++child){
        add_pseudocount(*child, ps);
    }
}

treeNode* get_site_clade(treeNode* node, int site){
    if (node->mutations.find(site) != node->mutations.end()){
        return node;
    }
    for (vector<treeNode*>::iterator child = node->children.begin(); child != node->children.end(); ++child){
        treeNode* result = get_site_clade(*child, site);
        if (result != NULL){
            return result;
        }
    }
    return NULL;
}

void insert_alleles(treeNode* n, map<long int, string>& alleles){
    if (n->parent != NULL){
        string namestr;
        for (set<long int>::iterator site = n->mutations.begin(); site != n->mutations.end();
            ++site){
            char buf[15];
            sprintf(buf, "%ld_", *site);
            namestr += string(buf) + alleles[*site] + "|";
        }
        if (n->mutations.size() > 0){
            n->name = namestr;
        }
    }
    for (vector<treeNode*>::iterator child = n->children.begin(); child != n->children.end();
        ++child){
        insert_alleles(*child, alleles);
    }
}

/**
 * Print a help message to the terminal and exit.
 */
void help(int code){
   fprintf(stderr, "trees2newick [OPTIONS]\n");
   fprintf(stderr, "Given output of SARGE, \
and a mapping of haplotype indices to individual names, converts all trees to \
human-readable Newick format with text haplotype labels inserted.\n");
    fprintf(stderr, "[OPTIONS]:\n");
    fprintf(stderr, "    --indvs -v (OPTIONAL) The file containing names of individuals \n");
    fprintf(stderr, "    --derived -d (OPTIONAL) Only show derived members of each tree -- \
remove any haplotypes that are direct leaves of the root\n");
    fprintf(stderr, "   --bufsize -b (OPTIONAL) The number of characters to read from the \
input file at a time\n");
    fprintf(stderr, "   --flatten_root -f (OPTIONAL) specify to change root branch length to \
0, in case there is a long root that prevents you from being able to see the other clades\n");
    fprintf(stderr, "   --site -S (OPTIONAL) a site position to create a FigTree annotation file for\n");
    fprintf(stderr, "   --outfile -o (OPTIONAL) with -s option, the name of the file to print \
the FigTree annotation information in.\n");
    fprintf(stderr, "   --remove_branchlengths -B take all branch lengths out of output trees\n");
    fprintf(stderr, "   --alleles -a (OPTIONAL) the .alleles file; will be used to annotate \
nodes in trees (SITES FILE ALSO REQUIRED).\n");
    fprintf(stderr, "   --sites -s (OPTIONAL) the .sites file; will be used to annotate \
nodes in trees (ALLELES FILE ALSO REQUIRED).\n");
    fprintf(stderr, "   --pseudocount -p (OPTIONAL) set to add a small value to each branch \
length, to make clades with 0 branch length more easy to see\n");
    fprintf(stderr, "   --metadata -m (OPTIONAL) output metadata about nodes (i.e. persistence \
along the chromosome) as an additional column. This can be read by view_trees.py.\n");
exit(code);
}

int main(int argc, char *argv[]) {
    static struct option long_options[] = {
       {"indvs", optional_argument, 0, 'v'},
       {"derived", no_argument, 0, 'd'},
       {"bufsize", optional_argument, 0, 'b'},
       {"site", required_argument, 0, 'S'},
       {"outfile", required_argument, 0, 'o'},
       {"alleles", required_argument, 0, 'a'},
       {"sites", required_argument, 0, 's'},
       {"pseudocount", required_argument, 0, 'p'},
       {"metadata", no_argument, 0, 'm'},
       {"help", optional_argument, 0, 'h'},
       {0, 0, 0, 0} 
    };
    
    vector<string> hapnames;
    int num_haplotypes;
    bool numhaps_read = false;
    string indvfilename;
    bool indvs_given = false;
    string allelesfile;
    string locsfile;
    float pseudocount = 0.0;
    bool use_pseudocount = false;
    bool metadata = false;
    
    int site = -1;
    string outfile;
    
    int bufsize = 1048576;
    
    bool derived_only = false;
    
    int option_index = 0;
    int ch;
    
    if (argc == 1){
        help(0);
    }
    while((ch = getopt_long(argc, argv, "v:n:b:S:o:a:s:p:mdh", long_options, &option_index )) != -1){
        switch(ch){
            case 0:
                // This option set a flag. No need to do anything here.
                break;
            case 'v':
                if (strcmp("n", optarg) == 0){
                    // No individuals to use.
                }
                else{
                    indvfilename = string(optarg);
                    indvs_given = true;
                }
                break;
            case 'd':
                derived_only = true;
                break;
            case 'm':
                metadata = true;
                break;
            case 'b':
                bufsize = atoi(optarg);
                break;
            case 'S':
                site = atoi(optarg);
                break;
            case 'o':
                outfile = optarg;
                break;
            case 'a':
                allelesfile = optarg;
                break;
            case 's':
                locsfile = optarg;
                break;
            case 'p':
                pseudocount = atof(optarg);
                use_pseudocount = true;
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

    FILE* outsitefile = NULL;
    if (site != -1){
        if (outfile.length() == 0){
            fprintf(stderr, "ERROR: you have given a site but not specified an output file for \
the FigTree annotation information.\n");
            exit(1);
        }
        else{
            outsitefile = fopen(outfile.c_str(), "w");
            if (!outsitefile){
                fprintf(stderr, "ERROR opening %s for writing.\n", outfile.c_str());
                exit(1);
            }
        }
    }
    
    if (indvs_given){
        parse_indvs(hapnames, indvfilename);
    }
    if ((locsfile.length() > 0 && allelesfile.length() == 0) ||
        (allelesfile.length() > 0 && locsfile.length() == 0)){
        fprintf(stderr, "ERROR: both the .sites and .alleles file are required \
if you want to annotate nodes.\n");
        exit(1);
    }
    
    map<string, set<pair<long int, long int> > > excl_bed;
    map<long int, string> alleles;
    bool has_alleles = false;
    if (allelesfile.length() > 0 && locsfile.length() > 0){
        vector<loc> locs = parse_locs(locsfile, excl_bed);
        ifstream infile(allelesfile);
        string anc;
        string der;
        int index = 0;
        has_alleles = true;
        while (infile >> anc >> der){
            alleles.insert(make_pair(locs[index].pos, der));
            ++index;
        }
    }
    
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
    
    while(!is.finished()){
 
        string chrom;
        long int pos;
        treeNode tree;
        read_sitedata(num_haplotypes, is, chrom, pos, tree);

        if (derived_only){
            tree.leaves.reset();
            tree.clear_cache();
        }
        
        if (use_pseudocount){
            add_pseudocount(&tree, pseudocount);
            tree.dist = pseudocount;
            tree.dist_norm = pseudocount;
        }
        else{
            tree.dist = 0.0;
            tree.dist_norm = 0.0;
        }

        if (has_alleles){
            insert_alleles(&tree, alleles);
        }
        
        if (site != -1 && pos == site){
            // Create FigTree annotation data for this site.
            fprintf(outsitefile, "taxa\tsite_%d\n", site);
            treeNode* der_clade = get_site_clade(&tree, site);
            cladeset leaves = der_clade->subtree_leaves();
            for (int i = 0; i < num_haplotypes; ++i){
                if (leaves.test(num_haplotypes-i-1)){
                    if (indvs_given){
                        fprintf(outsitefile, "%s\tderived\n", hapnames[i].c_str());
                    }
                    else{
                        fprintf(outsitefile, "%d\tderived\n", i);
                    }
                }
                else{
                    if (indvs_given){
                        fprintf(outsitefile, "%s\tancestral\n", hapnames[i].c_str());
                    }
                    else{
                        fprintf(outsitefile, "%d\tancestral\n", i);
                    }
                }
            }
        }
        print_site(tree, chrom, pos, indvs_given, hapnames, false, metadata);

    }
    
    if (outsitefile != NULL){
        fclose(outsitefile);
    }
}
