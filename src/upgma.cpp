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
#include <bitset>
#include <fstream>
#include <thread>
#include <ctime>
#include <cstdlib>
#include <sys/stat.h>
#include "sets.h"
#include "treeNode.h"
#include "serialize.h"
#include "common.h"

using std::cout;
using std::endl;
using namespace std;

// How often should we skip computing distances between haplotypes within
// a given "read" block of the input file?
float read_skip_prob = 0;

/**
 * Print a help message to the terminal and exit.
 */
void help(int code){
   fprintf(stderr, "upgma [OPTIONS]\n");
   fprintf(stderr, "Given an input haplotype matrix, builds a UPGMA tree. This can be used to choose haplotypes \
that will sample the deepest-diverging clades.\n");
   fprintf(stderr, "[OPTIONS]:\n");
   fprintf(stderr, "    --bufsize -b The number of characters to read from the \
input file at a time\n");
   fprintf(stderr, "    --input -i The main input file; a matrix of genotype data \
where rows are sites and columns are individuals. It must contain 1s and 0s only, \
with no spaces\n");
    fprintf(stderr, "   --threads -t The number of threads to use for parallel \
processing (default 4)\n");
    fprintf(stderr, "   --loci -l The file listing all sites in the haplotype matrix\n");
    fprintf(stderr, "   --numsites -n The maximum number of sites to consider when \
building the tree\n");
    fprintf(stderr, "   --sample -s The number of haplotypes to sample (optional)\n");
    fprintf(stderr, "   --prev -p A file of previously sampled haplotypes to exclude (optional, can specify multiple times)\n");
    fprintf(stderr, "   --regenerate -r If a tree file from a previous run is found and this \
option is specified, then a new tree will be built to replace the old one.\n");
    fprintf(stderr, "   --indv -v A file listing sample IDs (in order of index)\n");
    fprintf(stderr, "   --view -V Output the tree in Newick format\n");
    fprintf(stderr, "   --distmat -d Write distance matrix to given file name\n");
    fprintf(stderr, "   --mask -m A haplotype index to remove from analysis (may specify multiple times)\n");
    fprintf(stderr, "   --quote_names -q Add single quotes around leaf names. This is needed \
if you want to load the file into R using the phylogram package's read.dendrogram method.\n");
    exit(code);
}

/**
 * Function to read in a chunk from the input file.
 * Returns the number of new sites added; modifies hapList, infile, and eof by reference.
 * Unlike the similar function in sparge1, this stores data in columns rather
 * than rows, since we're interested more in haplotypes than sites for this step.
 */
long int expand_matrix(deque<cladeset > &hapList, 
    gzFile &infile, 
    long int bufsize,
    const unsigned int num_haplotypes, 
    bool &eof,
    set<int>& mask_haps){
    
    long unsigned int new_sites = 0;
    char chunk[bufsize+1];
    
    // Avoid letting garbage into this buffer next time
    memset(chunk, 0, bufsize+1);
    
    int chunk_size;
    chunk_size = gzread(infile, chunk, bufsize);
    chunk[chunk_size] = '\0';
    
    eof = gzeof(infile);
    
    if (!eof){
        // Find the last newline in the chunk and rewind to just before that position.
        int steps_back = 0;
        // Should not have to do this if bufsize is properly set.
        if (chunk[chunk_size-1] != '\n'){
            for (int i=chunk_size-1; i >= 0; i--){
                if (chunk[i] == '\n'){
                    // Truncate the string here.
                    chunk[i] = '\0';
                    gzseek(infile, -steps_back, SEEK_CUR);
                    break;
                }
                steps_back++;
            }
            chunk_size -= steps_back;
        }
        else{
            chunk[chunk_size-1] = '\0';
        }
    }
    
    // If we can skip this block, do that now.
    double randsamp = ((double) rand() / (RAND_MAX));
    if (randsamp < read_skip_prob){
        return 0;
    }
    
    // Turn the text into haplotypes.
    // First need to intialize a char array for each haplotype, of the appropriate
    // size.
    int numsites = (int) ceil((float) chunk_size/ ((float)num_haplotypes+1));
    char hap_strs[num_haplotypes][numsites + 1];
    
    for (int i = 0; i < num_haplotypes; ++i){
        memset(hap_strs[i], '\0', numsites+1);
    }
    
    int hap_index = 0;
    int site_index = 0;
    
    for (int i = 0; i < chunk_size; ++i){
        if (chunk[i] != '0' && chunk[i] != '1' && chunk[i] != '\n'){
            break;
        }
        if (chunk[i] == '\n'){
            hap_index = 0;
            site_index++;
        }
        else{
            hap_strs[hap_index][site_index] = chunk[i];
            hap_index++;
        }
    }
    
    new_sites += site_index+1;
    
    // Convert these into bitsets representing haplotypes across all sites read
    for (int i = 0; i < num_haplotypes; ++i){
        //fprintf(stderr, "%s\n", hap_strs[i]);
        string line = string(hap_strs[i]);
        cladeset bitstr(line);
        hapList.push_back(bitstr);
    }
    
    // Remove the unwanted haplotypes from analysis.
    for (set<int>::iterator mh = mask_haps.begin(); mh != mask_haps.end(); ++mh){
        hapList[*mh].reset();
    }
    
    return new_sites;
}

/**
 * Initialize a pre-declared distance matrix by setting all elements that will
 * be used and accessed to 0 and all other elements to -1.
 */
void init_distmat(vector<vector<float> >&dist_mat, int dim){
    for (int i = 0; i < dim-1; ++i){
        vector<float> row;
        for (int j = 0; j <= i; ++j){
            row.push_back(-1);
        }
        for (int j = i+1; j < dim; ++j){
            row.push_back(0);
        }
        dist_mat.push_back(row);
    }
}

/**
 * Print out a distance matrix to stdout.
 */
void print_distmat(vector<vector<float> >& dist_mat){
    for (vector<vector<float> >::iterator rowit = dist_mat.begin(); rowit != dist_mat.end();
        ++rowit){
        for (vector<float>::iterator colit = rowit->begin(); colit != rowit->end();
            ++colit){
            printf("%0.2f ", *colit);
        }
        printf("\n");
    }
}

/**
 * In a given window, compute pairwise distances between all haplotypes and
 * insert the results into the provided distance matrix. This will be 
 * done within a thread.
 */
void compute_diffs_window(vector<vector<float> >& dist_mat, 
    deque<cladeset >::iterator hapList, 
    const int num_sites,
    const int num_haplotypes){
    for (int i = 0; i < num_haplotypes - 1; ++i){
        for (int j = i+1; j < num_haplotypes; ++j){
            dist_mat[i][j] += ((float) ((*(hapList + i)) ^ (*(hapList + j))).count() / (float) num_sites);
            //dist_mat[i][j] += ((float) (hapList[i] ^ hapList[j]).count() / (float) num_sites);
        }
    }
    
}

/** 
 * Given a distance matrix from a single thread run and the master distance matrix,
 * dumps all values from the individual run into the master distance matrix.
 */
void combine_distmats(vector<vector<float> >& thread_dist_mat,
    vector<vector<float> >& master_dist_mat, 
    const int num_haplotypes){
    for (int i = 0; i < num_haplotypes - 1; ++i){
        for (int j = i+1; j < num_haplotypes; ++j){
            master_dist_mat[i][j] += thread_dist_mat[i][j];
        }
    }
}

/**
 * Counts the number of lines in a file. This is for the purpose of determining
 * how many sites are in the input file (by scanning the "loci" file instead of
 * the much-larger haplotype matrix).
 */
long int count_lines(const char* filename, long int bufsize){
    long int numlines = 0;
    FILE* fp = fopen(filename, "r");
    char buffer[bufsize];
    
    while(!feof(fp)){
        long int readsize;
        readsize = fread(&buffer, 1, bufsize, fp);
        for (int i = 0; i < readsize; ++i){
            if (buffer[i] == '\n'){
                numlines++;
            }
        }
        if (feof(fp) && buffer[readsize-1] != '\n'){
            numlines++;
        }
    }
    fclose(fp);
    return numlines;
}

void unnorm_branchlens(treeNode* node){
    node->dist_norm = node->dist;
    for (vector<treeNode*>::iterator c = node->children.begin(); c != node->children.end(); ++c){
        unnorm_branchlens((*c));
    }
}

void norm_brlens(treeNode* tree, float denom){
    tree->dist_norm = tree->dist_norm / denom;
    for (vector<treeNode*>::iterator c = tree->children.begin(); c != tree->children.end(); ++c){
        norm_brlens(*c, denom);
    }
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
       {"bufsize", optional_argument, 0, 'b'},
       {"input", required_argument, 0, 'i'},
       {"threads", required_argument, 0, 't'},
       {"loci", required_argument, 0, 'l'},
       {"numsites", required_argument, 0, 'n'},
       {"sample", required_argument, 0, 's'},
       {"prev", required_argument, 0, 'p'},
       {"regenerate", no_argument, 0, 'r'},
       {"indv", required_argument, 0, 'v'},
       {"view", no_argument, 0, 'V'},
       {"distmat", required_argument, 0, 'd'},
       {"quote_names", no_argument, 0, 'q'},
       {0, 0, 0, 0} 
    };
    
    // Set default values
    long int bufsize = 1048576;
    std::string infilename;
    int num_threads = 4;
    
    long int max_sites = -1;
    //long int max_sites = 1000000;
    
    const char* locifile;
    bool locifile_given = false;
    
    bool regenerate = false;
    bool view = false;
    
    bool prevsamples = false;
    vector<string> prevsampfiles;
    
    bool indvfile_given = false;
    string indvfile;
    
    struct stat stats;
    
    int sample = 0;
    
    bool write_distmat = false;
    string distmatfile;
    
    bool quote_names = false;
    
    int option_index = 0;
    int ch;
    
    
    set<int> mask_ints;
    
    if (argc == 1){
        help(0);
    }
    while((ch = getopt_long(argc, argv, "b:i:t:l:n:s:p:v:d:m:rqV", long_options, &option_index )) != -1){
        switch(ch){
            case 0:
                // This option set a flag. No need to do anything here.
                break;
            case 'b':
                bufsize = atol(optarg);
                //bufsize = strtol(optarg, NULL, 0);
                break;
            case 'i':
                infilename = string(optarg);
                break;
            case 't':
                num_threads = atoi(optarg);
                break;
            case 'l':
                locifile = (const char*) optarg;
                locifile_given = true;
                break;
            case 'n':
                max_sites = atol(optarg);
                break;
            case 's':
                sample = atoi(optarg);
                break;
            case 'p':
                if (stat(optarg, &stats) == 0){
                    prevsamples = true;
                    prevsampfiles.push_back(string(optarg));
                }
                else{
                    fprintf(stderr, "ERROR: previous sample file %s does not exist.\n", optarg);
                    exit(1);
                }
                break;
            case 'r':
                regenerate = true;
                break;
            case 'v':
                indvfile_given = true;
                indvfile = optarg;
                break;
            case 'V':
                view = true;
                break;
            case 'd':
                write_distmat = true;
                distmatfile = optarg;
                break;
            case 'm':
                mask_ints.insert(atoi(optarg));
                break;
            case 'q':
                quote_names = true;
                break;
            case '?':
                //help(0);
                break;
            default:
                help(0);
        }    
    }
    
    // Check arguments.
    if (sample > 0 && view){
        fprintf(stderr, "You have specified a number of haplotypes to sample AND \
to output the tree in Newick format. Only one of these options is permitted per run.\n");
        exit(1);
    }
    if (!locifile_given){
        fprintf(stderr, "ERROR: you must provide the file containing the list of input sites.\n");
        exit(1);
    }
    
    
    // Initialize random number seed
    srand(time(NULL));
    
    // Open main input file for reading
    gzFile infile;
    infile = gzopen(infilename.c_str(), "r");
    if (!infile){
        fprintf(stderr, "ERROR: the input file %s does not exist. Exiting.\n", infilename.c_str());
        exit(1);
    }
    
    // Load indv names, if given
    vector<string> indvs;
    if (indvfile_given){
        parse_indvs(indvs, indvfile);
    }
    else{
        fprintf(stderr, "ERROR: You must provide an individual ID input file\n");
        exit(1);
    }
    
    // Determine name of output/tree file
    string treefilename;
    int suffixstart = infilename.find_last_of('.');
    if (suffixstart > 0){
        treefilename = infilename.substr(0, suffixstart);
        // Check to see if this was just ".gz" -- if so, need to strip off
        // the next extension too.
        if (infilename.substr(suffixstart, infilename.length()).compare(".gz") == 0){
            suffixstart = treefilename.find_last_of('.');
            if (suffixstart > 0){
                treefilename = treefilename.substr(0, suffixstart);
            }
        }
    }
    treefilename += ".tree";
    
    int num_haplotypes;
    
    // This will store the tree (whether it's loaded from a file or built during
    // this run)    
    treeNode* tree;
    
    // Determine whether a tree was already created from this input file on
    // a previous run.
    bool treefile_exists = (stat(treefilename.c_str(), &stats) == 0);
    
    cladeset mask;
    
    if (treefile_exists && !regenerate){
        // Load tree from file.
        gzFile treefile;
        treefile = gzopen(treefilename.c_str(), "rb");
        if (!treefile){
            fprintf(stderr, "ERROR opening file %s for reading.\n", treefilename.c_str());
            exit(1);
        }
        instream_info is;
        instream_init(is, &treefile, bufsize);
        read_header(is, num_haplotypes, mask);
        tree = new treeNode();
        
        unserialize_treeNode(*tree, num_haplotypes, is);
        tree->setRoot();
        unnorm_branchlens(tree);
        //gzclose(treefile);
    }
    else{
        
        long int numsites = count_lines(locifile, bufsize);
        fprintf(stderr, "Counted %ld sites\n", numsites);
        bool eof = false;
        
        // Find out how many haplotypes there are by reading the first line of the input
        // file
        
        // Ensure that the buffer is big enough to read in at least one whole row of 
        // haplotype data from the input file at a time
        num_haplotypes = get_num_haplotypes(infile, bufsize);
        fprintf(stderr, "Read %d haplotypes\n", num_haplotypes);
        
        fprintf(stderr, "Computing distances\n");
        
        // Declare a distance matrix to use for building the UPGMA tree
        // Initialize to zero
        // This will be the "master" distance matrix that will store the sum of all
        // distances computed by all threads
        vector<vector<float> > dist_mat;
        init_distmat(dist_mat, num_haplotypes);

        // Now create one subordinate distance matrix per thread.
        vector<vector<vector<float> > > thread_dist_mats;
        thread threads[num_threads];
        bool launched[num_threads];
        int thread_index = 0;
        for (int i = 0; i < num_threads; ++i){
            // Initialize subordinate distance matrix
            vector<vector<float> > thread_dist_mat;
            init_distmat(thread_dist_mat, num_haplotypes);
            thread_dist_mats.push_back(thread_dist_mat);
            launched[i] = false;
        }

        // We'll now make sure that bufsize is a multiple of num_haplotypes+1 so 
        // we never have to gzseek().
        if (bufsize < num_haplotypes + 1){
            bufsize = num_haplotypes + 1;
        }
        else{
            bufsize = (long int) round( (float) (num_haplotypes+1) * round((float) bufsize / (float) (num_haplotypes+1)));
        }
        
        // Figure out the probability of skipping a given read block.
        float sites_per_read = round((float)bufsize / (float)(num_haplotypes+1));
        if (sites_per_read > numsites){
            sites_per_read = (float)numsites;
        }
        if (max_sites == -1 || max_sites > numsites){
            max_sites = numsites;
        }

        long int reads_needed = (long int) ceil((float)max_sites/sites_per_read);
        long int reads_tot = (long int) floor((float)numsites/sites_per_read);
        
        read_skip_prob = ((float)reads_tot - (float)reads_needed)/((float)reads_tot);
        if (max_sites == numsites || read_skip_prob < 0){
            read_skip_prob = 0;
        }
        
        // Create haplotype matrix (which can be resized)
        // This is a vector of bitsets, since each haplotype is just a bit string (0s and 1s)
        // The std::bitset implementation requires size to be declared at compile time,
        // however, and since the length of these should be equal to num_haplotypes,
        // we use the Boost dynamic bitset implementation instead.
        
        long unsigned int sites_tot = 0;

        vector<deque<cladeset > > hapLists;
        for (int i = 0; i < num_threads; ++i){
            deque<cladeset > hapList;
            hapLists.push_back(hapList);
        }
            
        while (!eof){

             // Check to see if a thread is already running in this slot.
            if (launched[thread_index]){
                threads[thread_index].join();
                combine_distmats(thread_dist_mats[thread_index], dist_mat, num_haplotypes);
                launched[thread_index] = false;
            }
            
            hapLists[thread_index].clear();
            long unsigned int num_sites = expand_matrix(hapLists[thread_index], infile, bufsize, 
                num_haplotypes, eof, mask_ints);
            
            // Only do stuff if we didn't skip this block.
            if (num_sites > 0){
                sites_tot += num_sites;
                
                fprintf(stderr, "Read %ld sites\n", sites_tot);
                
                // Process and remove
                
                // Create a new thread to run the computation.
                threads[thread_index] = thread(compute_diffs_window,
                    ref(thread_dist_mats[thread_index]),
                    hapLists[thread_index].begin(),
                    num_sites,
                    num_haplotypes);
                launched[thread_index] = true;
                
                thread_index++;
                if (thread_index >= num_threads){
                    thread_index = 0;
                }
            }
            
            /*
            // Bail out if we've seen enough sites already.
            if (sites_tot >= max_sites){
                break;
            }
            */
        }
        gzclose(infile);

        // Join all threads; merge results.
        for (int i = 0; i < num_threads; ++i){
            if (!launched[i]){
                // Just clear out distance matrix.
                for (vector<vector<float> >::iterator it1 = thread_dist_mats[i].begin();
                    it1 != thread_dist_mats[i].end(); ++it1){
                    it1->clear();
                }
                thread_dist_mats[i].clear();
            }
            else{
                fprintf(stderr, "(Final) joining thread %d\n", i);
                threads[i].join();
                combine_distmats(thread_dist_mats[i], dist_mat, num_haplotypes);
                // Throw out stored stuff.
                for (vector<vector<float> >::iterator it1 = thread_dist_mats[i].begin();
                    it1 != thread_dist_mats[i].end(); ++it1){
                    it1->clear();
                }
                thread_dist_mats[i].clear();
                launched[i] = false;
            }
        }
        thread_dist_mats.clear();
        
        // Now we can use UPGMA to build a tree.
        
        // For each entry in the current distance matrix, store a pointer to a treeNode
        // that that distance matrix entry represents.
        
        fprintf(stderr, "Building tree\n");
        
        vector<treeNode*> dist_mat_nodes;
        vector<float> node_heights;
        for (int i = 0 ; i < num_haplotypes; ++i){
            cladeset leaves;
            leaves.set(num_haplotypes-i-1);
            treeNode* node = new treeNode(num_haplotypes, leaves);
            dist_mat_nodes.push_back(node);
            node_heights.push_back(0);
        }
        
        int dist_mat_dim = num_haplotypes;
        
        if (write_distmat){
            FILE* distmatf = fopen(distmatfile.c_str(), "w");
            if (indvs.size() == num_haplotypes){
                for (int i = 0; i < num_haplotypes; ++i){
                    if (mask_ints.find(i) == mask_ints.end()){
                        fprintf(distmatf, "\t%s", indvs[i].c_str());
                    }
                }
            }
            else{
                for (int i = 0; i < num_haplotypes; ++i){
                    if (mask_ints.find(i) == mask_ints.end()){
                        fprintf(distmatf, "\thap%d", i+1);
                    }
                }
            }
            fprintf(distmatf, "\n");
            for (int i = 0; i < num_haplotypes; ++i){
                if (mask_ints.find(i) == mask_ints.end()){
                    if (indvs.size() == num_haplotypes){
                        fprintf(distmatf, "%s\t", indvs[i].c_str());
                    }
                    else{
                        fprintf(distmatf, "hap%d\t", i+1);
                    }
                    // Lower values
                    for (int j = 0; j < i; ++j){
                        if (mask_ints.find(j) == mask_ints.end()){
                            fprintf(distmatf, "\t%f", dist_mat[j][i]);
                        }
                    }
                    // Matching value
                    fprintf(distmatf, "\t0");
                    // Higher values
                    for (int j = i+1; j < num_haplotypes; ++j){
                        if (mask_ints.find(j) == mask_ints.end()){
                            fprintf(distmatf, "\t%f", dist_mat[i][j]);
                        }
                    }
                    fprintf(distmatf, "\n");
                }
            }
            fclose(distmatf);
        }
        
        while(dist_mat_dim > 1){
            fprintf(stderr, "Collapsing %d nodes\n", dist_mat_dim);
            
            // Get the smallest distance in the matrix.
            float smallest_dist = dist_mat[0][1];
            int smallest_i = 0;
            int smallest_j = 1;
            for (int i = 0; i < dist_mat_dim-1; ++i){
                for (int j = i+1; j < dist_mat_dim; ++j){
                    // Use -1 to represent things removed from the matrix.
                    if (dist_mat[i][j] != -1 && dist_mat[i][j] < smallest_dist){
                        smallest_dist = dist_mat[i][j];
                        smallest_i = i;
                        smallest_j = j;
                    }
                }
            }
            
            // Build a node connecting the closest things.
            treeNode* a = dist_mat_nodes[smallest_i];
            treeNode* b = dist_mat_nodes[smallest_j];
            cladeset a_subtree = a->subtree_leaves();
            cladeset b_subtree = b->subtree_leaves();
            cladeset parent_leaves;
            treeNode* ab = new treeNode(num_haplotypes, parent_leaves);
            a->parent = ab;
            b->parent = ab;
            
            a->dist = a->dist_norm = smallest_dist/2 + node_heights[smallest_j]/2;
            b->dist = b->dist_norm = smallest_dist/2 + node_heights[smallest_i]/2;
            
            ab->children.push_back(a);
            ab->children.push_back(b);
            
            int a_weight = a_subtree.count();
            int b_weight = b_subtree.count();
            
            float ab_height = (node_heights[smallest_i] + node_heights[smallest_j] + smallest_dist);
            
            // Rebuild distance matrix.
            dist_mat_dim--;
            
            // Create new distance matrix and matrix of tree pointers.
            // Allocate the same amount of space, so assignment continues to work,
            // but ignore the last elements.
            vector<vector<float> > dist_mat_new;
            init_distmat(dist_mat_new, dist_mat_dim);
            vector<treeNode*> dist_mat_nodes_new;
            vector<float> node_heights_new;
            
            // Replace previous i index with new node. Delete previous j index.
            for (int i = 0; i < dist_mat_dim; ++i){
                if (i != smallest_j){
                    int new_i = i;
                    if (i > smallest_j){
                        new_i--;
                    }
                    for (int j = i+1; j < dist_mat_dim+1; ++j){
                        if (j != smallest_j){
                            int new_j = j;
                            if (j > smallest_j){
                                new_j--;
                            }
                            if (i == smallest_i){
                                // Distance will be from the new node to the other stuff.
                                if (j > smallest_j){
                                    dist_mat_new[new_i][new_j] = (dist_mat[smallest_i][j]*a_weight + dist_mat[smallest_j][j]*b_weight)/(a_weight+b_weight);
                                }
                                else{
                                    dist_mat_new[new_i][new_j] = (dist_mat[smallest_i][j]*a_weight + dist_mat[j][smallest_j]*b_weight)/(a_weight+b_weight);
                                }
                            }
                            else{
                                // Copy over previous value.
                                dist_mat_new[new_i][new_j] = dist_mat[i][j];
                            }
                        }
                    }
                }
            }
            
            // Update other stuff
            for (int i = 0; i < dist_mat_dim+1; ++i){
                if (i == smallest_i){
                    // Add info for new node
                    dist_mat_nodes_new.push_back(ab);
                    node_heights_new.push_back(ab_height);
                }
                else if (i == smallest_j){
                    // Skip index.
                }
                else{
                    dist_mat_nodes_new.push_back(dist_mat_nodes[i]);
                    node_heights_new.push_back(node_heights[i]);
                }
            }
            
            dist_mat = dist_mat_new;
            dist_mat_nodes = dist_mat_nodes_new;
            node_heights = node_heights_new;
            
        }
        
        gzFile treefile = gzopen(treefilename.c_str(), "w");
        serialize_header(treefile, num_haplotypes, mask);
        dist_mat_nodes[0]->serialize(treefile);
        gzclose(treefile);
        
        tree = dist_mat_nodes[0];
    }
    
    // We now have a tree loaded in memory. 
    
    // Load previously sampled haplotypes, if provided.
    if (prevsamples){
        cladeset prevsamp;
        // Dump in all previous samples 
        for (vector<string>::iterator prevsampfile = prevsampfiles.begin();
            prevsampfile != prevsampfiles.end(); ++prevsampfile){
            load_hap_sample(*prevsampfile, prevsamp, num_haplotypes);
        }
        prevsamp.flip();
        // Set all previously-sampled indices to ancestral in the tree.
        tree->mask_haps(prevsamp);
    }
    
    if (sample > 0){
        // Sample the requested number of deepest-diverging haplotypes from the
        // tree.
        set<unsigned int> sampled = sample_haplotypes(*tree, sample, num_haplotypes);
        for (set<unsigned int>::iterator it = sampled.begin(); it != sampled.end();
            ++it){
            if (indvfile_given){
                // Print sample ID
                printf("%s\n", indvs[*it].c_str());
            }
            else{
                // Just print haplotype index
                printf("%d\n", *it);
            }
        }
    }
    if (view){
        
        norm_brlens(tree, tree->dist_below);
        
        if (quote_names){
            for (int i = 0; i < indvs.size(); ++i){
                indvs[i] = "'" + indvs[i] + "'";
            }
        }
        
        // Print out the tree in Newick format.
        if (indvfile_given){
            printf("%s\n", tree->newick(false, true, indvs).c_str());
        }
        else{
            printf("%s\n", tree->newick(false, false, indvs).c_str());
        }
    }
    
    //gzclose(infile);
}
