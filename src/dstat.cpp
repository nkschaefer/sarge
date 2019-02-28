/**
 * Computes admixture using the D/f-hat statistic from SARGE input files.
 * Uses the same parameters as admix_scan. Can be used to compare the results
 * of the two programs.
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
 * Function to read in a chunk from the input file.
 * Returns the number of new sites added; modifies hapList, infile, and eof by reference.
 */
void read_genotypes(map<string, cladeset>& ingroup_clades,
    cladeset& outgroup,
    cladeset& admix_in,
    cladeset& admix_out,
    cladeset& admixers,
    bool f_hat,
    float& admix_num_sum,
    map<string, float>& sums_numerator,
    map<string, float>& sums_denominator,
    int& count_sites,
    gzFile &infile, 
    long int bufsize,
    const unsigned int num_haplotypes){
    
    long unsigned int new_sites = 0;
    char* chunk = (char*) malloc(bufsize+1);
    
    count_sites = 0;
    
    bool eof = false;
    
    int line_index = 0;
    
    while(!eof){
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
            }
            else{
                chunk[chunk_size-1] = '\0';
            }
        }
        
        // Split lines   
        istringstream linesplitter(chunk);
        string line;
        int token_index = 0;

        while(std::getline(linesplitter, line, '\n')){
            if (line.size() == num_haplotypes){
                
                // Skip lines that were in excluded regions (from BED file).
                if (del_sites.size() > 0 && *del_sites.begin() == line_index){
                    del_sites.erase(del_sites.begin());
                }
                else{
                    new_sites++;
                        
                    // Convert to bit string and store
                    cladeset bitstr(line);
                    
                    if (bitstr.count() < num_haplotypes && bitstr.count() > 0){
                        // Convert statistic according to Durand et al 2011
                        // p1 = outgroup
                        // p2 = ingroup
                        // p3 = admixer
                        // p4 = ancestral (always 0)
                        float p1 = (float) (bitstr & outgroup).count() / (float) outgroup.count();
                        float p3 = (float) (bitstr & admixers).count() / (float) admixers.count();
                        float p4 = 0.0;
                        
                        for (map<string, cladeset>::iterator ig = ingroup_clades.begin();
                            ig != ingroup_clades.end(); ++ig){
                            
                            float p2 = (float) (bitstr & ig->second).count() / (float) ig->second.count();
                            sums_numerator[ig->first] += ((1-p1)*p2*p3*(1-p4) - p1*(1-p2)*p3*(1-p4));
                            sums_denominator[ig->first] += ((1-p1)*p2*p3*(1-p4) + p1*(1-p2)*p3*(1-p4));
                        }
                        
                        if (f_hat){
                            float p3_1 = (float) (bitstr & admix_in).count() / (float) admix_in.count();
                            float p3_2 = (float) (bitstr & admix_out).count() / (float) admix_out.count();
                            admix_num_sum += ((1-p1)*p3_1*p3_2*(1-p4) - p1*(1-p3_1)*p3_2*(1-p4));
                        }
                        
                        count_sites++;
                    }
                }
                line_index++;
            }
        }
        // Avoid letting garbage into this buffer next time
        memset(chunk, 0, bufsize);
    }
    free(chunk);
}

/**
 * Print a help message to the terminal and exit.
 */
void help(int code){
   fprintf(stderr, "dstat [OPTIONS]\n");
   fprintf(stderr, "Computes admixture (using same parameters as admix_scan) \
but uses the D/f-hat statistics on the SARGE input file rather than looking at \
output trees.\n");
   fprintf(stderr, "[OPTIONS]:\n");
    fprintf(stderr, "   --geno -g Input genotype file (can be gzipped)\n");
   fprintf(stderr, "    --ingroup -i One or more populations that are admixed (can \
specify more than once)\n");
   fprintf(stderr, "    --outgroup -o One or more populations that are not admixed (can \
specify more than once)\n");
   fprintf(stderr, "    --admixer -a The admixing population of interest\n");
fprintf(stderr, "    --indvs -v (OPTIONAL) The file containing names of individuals \n");
   fprintf(stderr, "    --ploidy -n (OPTIONAL) The number of haplotypes per individual (\
default 2) \n");
    fprintf(stderr, "    --pops -p A file mapping individual names to population names, \
tab-separated\n");
    fprintf(stderr, "    --proportion -f Specify to calculate admixture proportion per \
population instead of D statistic. Will randomly hold out 50%% of outgroup haplotypes\n");
    fprintf(stderr, "    --exclude_bed -B Specify a bed file of sites to EXCLUDE from analysis, \
if desired. If you have a BED file of sites to INCLUDE instead, please use bedtools complement \
to convert it to excludable sites.\n");
    fprintf(stderr, "    --sites -s The sites file to be used as SARGE input.\n");
    fprintf(stderr, "    --bufsize -b The number of characters to read from the \
input file at a time\n");
exit(code);
}

// Define a structure to store summary data about admixture within a population.
struct pop_stats {
    string popname;
    float perc_admixed;
    float perc_admixed_sites;
    float tract_len;
    long int tract_len_min;
    long int tract_len_max;
    float allele_freq;
};

bool pop_stats_comp(pop_stats s1, pop_stats s2){
    if (s1.perc_admixed == s2.perc_admixed){
        return (s1.popname.compare(s2.popname) > 0);
    }
    else{
        return s1.perc_admixed > s2.perc_admixed;
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
       {"ingroup", required_argument, 0, 'i'},
       {"outgroup", required_argument, 0, 'o'},
       {"admixer", required_argument, 0, 'a'},
       {"indvs", required_argument, 0, 'v'},
       {"ploidy", optional_argument, 0, 'n'},
       {"pops", required_argument, 0, 'p'},
       {"bufsize", optional_argument, 0, 'b'},
       {"geno", required_argument, 0, 'g'},
       {"proportion", no_argument, 0, 'f'},
       {"exclude_bed", required_argument, 0, 'B'},
       {"sites", required_argument, 0, 's'},
       {"help", optional_argument, 0, 'h'},
       {0, 0, 0, 0} 
    };
    
    // Set default values
    int bufsize = 1048576;
    int ploidy = 2;
    
    vector<string> ingroup_pops;
    vector<string> outgroup_pops;
    string admixer;
    string indvfilename;
    string popfilename;
    string genofilename;
    string bedfilename;
    string sitesfilename;
    bool f_hat = false;
    
    int option_index = 0;
    int ch;
    
    if (argc == 1){
        help(0);
    }
    while((ch = getopt_long(argc, argv, "i:o:a:v:n:p:b:g:B:s:fh", long_options, &option_index )) != -1){
        switch(ch){
            case 0:
                // This option set a flag. No need to do anything here.
                break;
            case 'i':
                ingroup_pops.push_back(string(optarg));
                break;
            case 'o':
                outgroup_pops.push_back(string(optarg));
                break;
            case 'a':
                admixer = string(optarg);
                break;
            case 'n':
                ploidy = atoi(optarg);
                break;
            case 'v':
                indvfilename = string(optarg);
                break;
            case 'p':
                popfilename = string(optarg);
                break;
            case 'b':
                bufsize = atoi(optarg);
                break;
            case 'g':
                genofilename = optarg;
                break;
            case 'f':
                f_hat = true;
                break;
            case 'B':
                bedfilename = optarg;
                break;
            case 's':
                sitesfilename = optarg;
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
    
    if (ingroup_pops.size() == 0 || outgroup_pops.size() == 0 || admixer.length() == 0){
        fprintf(stderr, "ERROR: you must provide at least one ingroup, outgroup, and admixing population.\n");
        exit(1);
    }
    if (genofilename.length() == 0){
        fprintf(stderr, "ERROR: input genotype file name not provided.\n");
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
    if (sitesfilename.length() == 0){
        fprintf(stderr, "ERROR: you must provide a file listing sites.\n");
        exit(1);
    }
    
    map<string, set<pair<long int, long int> > > exclude_bed;
    if (bedfilename.length() > 0){
        parse_exclude_bed(bedfilename, exclude_bed);
    }
    vector<loc> locs = parse_locs(sitesfilename, exclude_bed);
    
    vector<string> indvs;
    parse_indvs(indvs, indvfilename);
    
    map<string, int> indv2hap;
    parse_indvs_map(indv2hap, indvfilename);
    map<string, string> indv2pop;
    map<string, vector<string> > pop2indv;
    parse_pops(indv2pop, pop2indv, popfilename, ploidy);
    
    vector<string> leafname_pops;
    for (vector<string>::iterator indvit = indvs.begin(); indvit != indvs.end();
        ++indvit){
        //string indv_sub = (*indvit).substr(0, (*indvit).find('-'));
        leafname_pops.push_back(indv2pop[*indvit]);
    }
    
    gzFile geno_in = gzopen(genofilename.c_str(), "r");
    if (!geno_in){
        fprintf(stderr, "ERROR opening %s for reading.\n", genofilename.c_str());
        exit(1);
    }
    // Ensure that the buffer is big enough to read in at least one whole row of 
    // haplotype data from the input file at a time
    const unsigned int num_haplotypes = get_num_haplotypes(geno_in, bufsize);
    
    // We'll now make sure that bufsize is a multiple of num_haplotypes+1 so 
    // we never have to gzseek().
    if (bufsize < num_haplotypes + 1){
        bufsize = num_haplotypes + 1;
    }
    else{
        bufsize = (long int) round( (float) (num_haplotypes+1) * round((float) bufsize / (float) (num_haplotypes+1)));
    }
    
    // Initialize random number seed    
    srand(time(NULL));
    
    // Store ingroup clades separately per population
    map<string, cladeset> ingroup_clades;
    map<string, float> sums_numerator;
    map<string, float> sums_denominator;
    // Store numerator sum for outgroup subsample, so we can calculate f hat
    float admix_num_sum = 0.0;
    
    long int count_ig_haps = 0;
    
    for (vector<string>::iterator ingroup_pop = ingroup_pops.begin();
        ingroup_pop != ingroup_pops.end(); ++ingroup_pop){
        
        sums_numerator.insert(make_pair(*ingroup_pop, 0.0));
        sums_denominator.insert(make_pair(*ingroup_pop, 0.0));
        
        set<unsigned int> ingroup_set;
            
        for (vector<string>::iterator ingroup_name = pop2indv[*ingroup_pop].begin();
            ingroup_name != pop2indv[*ingroup_pop].end(); ++ingroup_name){
            if (indv2hap.count(*ingroup_name) > 0){
                ingroup_set.insert(indv2hap[*ingroup_name]);
                //fprintf(stderr, "ingroup %d %s\n", indv2hap[*ingroup_name], ingroup_name->c_str());
            }
        }
        
        cladeset ingroup_clade = set2bitset(ingroup_set, num_haplotypes);
        count_ig_haps += ingroup_clade.count();
        ingroup_clades.insert(make_pair(*ingroup_pop, ingroup_clade));
    }
    
    set<unsigned int> outgroup_set;
    
    for (vector<string>::iterator outgroup_pop = outgroup_pops.begin();
        outgroup_pop != outgroup_pops.end(); ++outgroup_pop){
        for (vector<string>::iterator outgroup_name = pop2indv[*outgroup_pop].begin();
            outgroup_name != pop2indv[*outgroup_pop].end(); ++outgroup_name){
            if (indv2hap.count(*outgroup_name) > 0){
                outgroup_set.insert(indv2hap[*outgroup_name]);
                //fprintf(stderr, "outgroup %d %s\n", indv2hap[*outgroup_name], outgroup_name->c_str());
            }
        }
    }

    cladeset outgroup_clade = set2bitset(outgroup_set, num_haplotypes);
    
    vector<cladeset > admix_clades;
    set<unsigned int> admix_haps_all;
    for (vector<string>::iterator admix_hap = pop2indv[admixer].begin();
        admix_hap != pop2indv[admixer].end(); ++admix_hap){
        if (indv2hap.count(*admix_hap) > 0){
            admix_haps_all.insert(indv2hap[*admix_hap]);
            set<unsigned int> admix_hap_set;
            admix_hap_set.insert(indv2hap[*admix_hap]);
            cladeset admix_clade = set2bitset(admix_hap_set, num_haplotypes);
            admix_clades.push_back(admix_clade);
            //fprintf(stderr, "admix %d %s\n", indv2hap[*admix_hap], admix_hap->c_str());
        }
    }
    
    cladeset admix_clade;
    for (vector<cladeset>::iterator ac = admix_clades.begin(); ac != admix_clades.end(); ++ac){
        admix_clade |= *ac;
    }
    
    // Clades for calculating f hat
    set<unsigned int> admix_in_set;
    set<unsigned int> admix_out_set;
    cladeset admix_in;
    cladeset admix_out;
    
    if (f_hat){
        // Sort randomly
        vector<unsigned int> admix_unsorted;
        for (set<unsigned int>::iterator a = admix_haps_all.begin(); a != admix_haps_all.end();
            ++a){
            admix_unsorted.push_back(*a);
        }
        random_shuffle(admix_unsorted.begin(), admix_unsorted.end());
        for (int i = 0; i < (int)floor((float)admix_haps_all.size() / 2.0); ++i){
            admix_in_set.insert(admix_unsorted[i]);
        }
        for (int i = admix_in_set.size(); i < admix_unsorted.size(); ++i){
            admix_out_set.insert(admix_unsorted[i]);
        }
        
        if (admix_in_set.size() < 1 || admix_out_set.size() < 1){
            fprintf(stderr, "ERROR: unable to split admixers into two subpopulations for \
calculating f-hat. Please re-run with more haplotypes.\n");
            exit(1);
        }
        admix_in = set2bitset(admix_in_set, num_haplotypes);
        admix_out = set2bitset(admix_out_set, num_haplotypes);
        fprintf(stderr, "Split admixers into %ld ingroup and %ld outgroup haplotypes for f-hat calculation\n",
            admix_in.count(), admix_out.count());
    }
    
    int num_sites = 0;
    
    

    fprintf(stderr, "Ingroup haps: %ld Outgroup haps: %ld Admixing haps: %ld\n", 
        count_ig_haps, outgroup_set.size(), admix_clades.size());
    
    if (count_ig_haps == 0 || outgroup_set.size() == 0 || admix_clades.size() == 0){
        fprintf(stderr, "ERROR: the haplotypes you have subsampled do not include any individuals \
from either the ingroup, outgroup or admixing group.\n");
        exit(1);
    }

    read_genotypes(ingroup_clades, outgroup_clade, admix_in, admix_out, admix_clade, 
        f_hat, admix_num_sum, sums_numerator, sums_denominator, num_sites, geno_in, bufsize, num_haplotypes);
    gzclose(geno_in);
    
    // Print results.
    fprintf(stdout, "== Population-level D-statistics by ingroup population ==\n");
    for (map<string, cladeset>::iterator pop = ingroup_clades.begin(); pop != ingroup_clades.end(); ++pop){
        float d = sums_numerator[pop->first] / sums_denominator[pop->first];
        if (f_hat){
            float f = sums_numerator[pop->first] / admix_num_sum;
            fprintf(stdout, "%s\t%f\t%f\n", pop->first.c_str(), d, f);
        }
        else{
            fprintf(stdout, "%s\t%f\n", pop->first.c_str(), d);
        }
    }
    
}
