#include <getopt.h>
#include <argp.h>
#include <zlib.h>
#include <string>
#include <algorithm>
#include <vector>
#include <iterator>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <cstdlib>
#include <utility>
#include <math.h>
#ifdef HTSLIB
    #include <htslib/vcf.h>
    #include <htslib/kseq.h>
#endif

using std::cout;
using std::endl;
using namespace std;

// Initialize FASTA/FASTQ parser
#ifdef HTSLIB
    KSEQ_INIT(gzFile, gzread);
#endif

/**
 * Trim whitespace and newlines and tabs off the right end of a string.
 */
void rstrip(string& str){
    long int rpos = str.length()-1;
    while (rpos >= 0){
        if (str[rpos] == ' ' || str[rpos] == '\n' || str[rpos] == '\t'){
            rpos--;
        }
        else{
            break;
        }
    }
    if (rpos < str.length()-1){
        str = str.substr(0, rpos+1);
    }
}

/**
 * Trim whitespace and newlines and tabs off the left end of a string.
 */
void lstrip(string& str){
    long int lpos = 0;
    while (lpos < str.length()){
        if (str[lpos] == ' ' || str[lpos] == '\n' || str[lpos] == '\t'){
            lpos++;
        }
        else{
            break;
        }
    }
    if (lpos > 0){
        str = str.substr(lpos, str.length() - lpos);
    }
}

/**
 * Trim whitespace and newlines and tabs off both ends of a string.
 */
void strip(string& str){
    rstrip(str);
    lstrip(str);
}

/**
 * Function to split a string by whitespace
 */
void splitstr(const string& str, vector<string>& output) { 
    istringstream buffer(str);
    copy(istream_iterator<string>(buffer), 
        istream_iterator<string>(),
        back_inserter(output));
}

char upper(char base){
    switch(base){
        case 'a':
        case 'A':
            return 'A';
            break;
        case 'c':
        case 'C':
            return 'C';
            break;
        case 'g':
        case 'G':
            return 'G';
            break;
        case 't':
        case 'T':
            return 'T';
            break;
        default:
            return base;
            break;
    }
    return base;
}

bool isbase(char base){
    if (base == 'A' || base == 'C' || base == 'G' || base == 'T'){
        return true;
    }
    return false;
}

void parse_depthfile(string& depthfile, map<string, pair<int, int> >& depths){
    ifstream infile(depthfile.c_str());
    string id;
    int mindepth;
    int maxdepth;
    while(infile >> id >> mindepth >> maxdepth){
        if (id.length() > 0){
            depths.insert(make_pair(id, make_pair(mindepth, maxdepth)));
        }
    }
}

/**
 * Print a help message to the terminal and exit.
 */
void help(int code){
    fprintf(stderr, "vcf2sarge [OPTIONS]\n");
    fprintf(stderr, "Given a VCF file (from stdin), as well as a reference \
and ancestral FASTA sequence for the same chromosome, converts it to the format \
expected by sarge.\n");
    fprintf(stderr, "[OPTIONS]:\n");
    fprintf(stderr, "    --ancestral -a A FASTA file containing the ancestral sequence \
for the given chromosome (REQUIRED)\n");
    fprintf(stderr, "    --reference -r A FASTA file containing the reference genome \
sequence for the given chromosome (REQUIRED)\n");
    fprintf(stderr, "    --output_prefix -o The base name of all output files (REQUIRED). \
Created files will have the names <prefix>.geno.gz, <prefix>.sites, and <prefix>.haps.\n");
    fprintf(stderr, "    --chrom -c (REQUIRED) The chromosome being considered\n");
    fprintf(stderr, "    --start -s (OPTIONAL) if you are analyzing a region that does not \
include the beginning of the chromosome (i.e. piping from tabix), give the start coordinate \
here, or else fixed differences between all haplotypes and the ancestral sequence upstream \
of this position will be included.\n");
    fprintf(stderr, "    --end -e (OPTIONAL) if you are analyzing a region that does not \
include the end of the chromosome (i.e. piping from tabix), give the end coordinate here, \
or else fixed differences between all haplotypes and the ancestral sequence downstream \
of this position will be included.\n");
    fprintf(stderr, "    --ignore_phasing -p (OPTIONAL) set this flag to not check whether \
sites are phased (it will assume you know these are phased correctly).\n");
    fprintf(stderr, "    --quality -q (OPTIONAL) the minimum genotype quality required for \
a genotype to be included. If any genotype falls below this at a site, the site will be \
skipped.\n");
    fprintf(stderr, "    --depth -d (OPTIONAL) the minimum depth required for a genotype \
to be included. If any genotype falls below this at a site, the site will be excluded.\n");
    fprintf(stderr, "    --depthfile -D (OPTIONAL) A file where each line is a sample name, \
minimum coverage, and maximum coverage for that sample, space separated. This will supersede \
depth (-d) if also given.\n");
    fprintf(stderr, "    --help -h Display this message and exit.\n");
    exit(code);
}

int main(int argc, char *argv[]) {    
    
    /** Define arguments 
     * http://www.gnu.org/software/libc/manual/html_node/Getopt-Long-Options.html#Getopt-Long-Options
     * http://www.gnu.org/software/libc/manual/html_node/Getopt-Long-Option-Example.html#Getopt-Long-Option-Example
     * Fields for each argument: name, has_arg (values: no_argument, required_argument,
     *     optional_argument)
     * flag = int value to store flag for the option, or NULL if option is string
     * val = short name for string option, or NULL
     */
     
    static struct option long_options[] = {
       {"ancestral", required_argument, 0, 'a'},
       {"reference", required_argument, 0, 'r'},
       {"output_prefix", required_argument, 0, 'o'},
       {"chromosome", required_argument, 0, 'c'},
       {"start", required_argument, 0, 's'},
       {"end", required_argument, 0, 'e'},
       {"ignore_phasing", no_argument, 0, 'p'},
       {"quality", required_argument, 0, 'q'},
       {"depth", required_argument, 0, 'd'},
       {"depthfile", required_argument, 0, 'D'},
       {0, 0, 0, 0} 
    };
    
    // Set default values
    string ancestralfile;
    string referencefile;
    string out_prefix;
    string chromosome;
    string depthfile;
    bool depthfile_given = false;
    bool bcf = false;
    bool ignore_phasing = false;
    int start = -1;
    int end = -1;
    int quality = -1;
    int depth = -1;
    
    int option_index = 0;
    int ch;
    
    if (argc == 1){
        help(0);
    }
    while((ch = getopt_long(argc, argv, "a:r:o:c:s:e:q:d:D:ph", long_options, &option_index )) != -1){
        switch(ch){
            case 0:
                // This option set a flag. No need to do anything here.
                break;
            case 'h':
                help(0);
                break;
            case 'a':
                ancestralfile = optarg;
                break;
            case 'r':
                referencefile = optarg;
                break;
            case 'o':
                out_prefix = optarg;
                break;
            case 'c':
                chromosome = optarg;
                break;
            case 's':
                start = atoi(optarg);
                break;
            case 'e':
                end = atoi(optarg);
                break;
            case 'p':
                ignore_phasing = true;
                break;
            case 'q':
                quality = atoi(optarg);
                break;
            case 'd':
                depth = atoi(optarg);
                break;
            case 'D':
                depthfile = optarg;
                depthfile_given = true;
                break;
            default:
                help(0);
                break;
        }    
    }
    
    // Error check arguments.
    if (ancestralfile.length() == 0 || referencefile.length() == 0){
        fprintf(stderr, "ERROR: you must provide an ancestral and reference genome \
file for the given chromosome.\n");
        exit(1);
    }
    
    if (start == 0 || (start != -1 && end != -1 && start > end)){
        fprintf(stderr, "Invalid start & end coordinates: %d %d\n", start, end);
        fprintf(stderr, "Coordinates should be 1-based, not 0-based\n");
        exit(1);
    }
    
    #ifndef HTSLIB
        fprintf(stderr, "ERROR: HTSLib not installed on this system. Conversion of VCF to SARGE \
format requires that you have HTSLib installed and available in the LIBRARY_PATH and LD_LIBRARY_PATH \
environment variables. Then recompile SARGE with HTSLIB=1.\n");
        exit(1);
    #else
    
        // Check that output prefix was provided
        if (out_prefix.length() == 0){
            fprintf(stderr, "ERROR: you must provide a valid output prefix.\n");
            exit(1);
        }
        
        map<string, pair<int, int> > depths;
        if (depthfile_given){
            parse_depthfile(depthfile, depths);
        }
        
        // Create output files
        string out_sites_name = out_prefix + ".sites";
        string out_haps_name = out_prefix + ".haps";
        string out_geno_name = out_prefix + ".geno.gz";
        string out_alleles_name = out_prefix + ".alleles";
        
        FILE* out_sites_f = fopen(out_sites_name.c_str(), "w");
        if (out_sites_f == NULL){
            fprintf(stderr, "ERROR opening file %s for writing.\n", out_sites_name.c_str());
            exit(1);
        }
        FILE* out_haps_f = fopen(out_haps_name.c_str(), "w");
        if (out_haps_f == NULL){
            fprintf(stderr, "ERROR opening file %s for writing.\n", out_haps_name.c_str());
            exit(1);
        }
        gzFile out_geno_f = gzopen(out_geno_name.c_str(), "w");
        if (!out_geno_f){
            fprintf(stderr, "ERROR opening file %s for writing.\n", out_geno_name.c_str());
            exit(1);
        }
        FILE* out_alleles_f = fopen(out_alleles_name.c_str(), "w");
        
        // Initialize fasta parser.
        int anc_progress;
        int ref_progress;
        gzFile anc_fp = gzopen(ancestralfile.c_str(), "r");
        if (!anc_fp){
            fprintf(stderr, "ERROR opening file %s for reading.\n", ancestralfile.c_str());
            exit(1);
        }
        gzFile ref_fp = gzopen(referencefile.c_str(), "r");
        if (!ref_fp){
            fprintf(stderr, "ERROR opening file %s for reading.\n", referencefile.c_str());
            exit(1);
        }
        kseq_t* anc_seq = kseq_init(anc_fp);
        kseq_t* ref_seq = kseq_init(ref_fp);
        
        fprintf(stderr, "Loading ancestral & reference sequences from FASTA files...\n");

        // Load ancestral sequence
        bool anc_found = false;
        while((anc_progress = kseq_read(anc_seq) >= 0)){
            if (strcmp(anc_seq->name.s, chromosome.c_str()) == 0){
                anc_found = true;
                break;
            }
        }
        if (!anc_found){
            fprintf(stderr, "ERROR: chromosome %s not found in ancestral FASTA file\n", chromosome.c_str());
            exit(1);
        }
        
        // Load reference sequence
        bool ref_found = false;
        while((ref_progress = kseq_read(ref_seq) >= 0)){
            if (strcmp(ref_seq->name.s, chromosome.c_str()) == 0){
                ref_found = true;
                break;
            }
        }
        
        if (!ref_found){
            fprintf(stderr, "ERROR: chromosome %s not found in reference FASTA file\n", chromosome.c_str());
            exit(1);
        }
        
        // Find the first and last position in the ancestral & reference sequences for which there is 
        // an allele.
        int prevpos = -1;
        int lastpos = -1;
        
        if (start != -1){
            // Set to one before actual position, so this position gets checked.
            prevpos = start-1;
        }
        else{
            for (int ind = 0; ind < anc_seq->seq.l; ++ind){
                if (isbase(upper(anc_seq->seq.s[ind])) && isbase(upper(ref_seq->seq.s[ind]))){
                    prevpos = ind;
                    break;
                }
            }
        }
        
        int minseqlen = ref_seq->seq.l;
        if (anc_seq->seq.l < minseqlen){
            minseqlen = anc_seq->seq.l;
        }
        
        if (end != -1){
            lastpos = end;
        }
        else{
            for (int ind = minseqlen-1; ind >= 0; --ind){
                if (isbase(upper(anc_seq->seq.s[ind])) && isbase(upper(ref_seq->seq.s[ind]))){
                    lastpos = ind + 1;
                    break;
                }
            }
        }
        
        string prevchrom;
        
        int num_haps = 0;
        
        // How often should progress messages be printed?
        int progress = 50000;
        int last_printed = prevpos;
        
        // Convert depths to look up by index
        vector<pair<int, int> > depths_index;
        
        if (bcf){
            // Read BCF from stdin.
            bcf_hdr_t* bcf_header;
            bcf1_t* bcf_record = bcf_init();
            htsFile* bcf_reader = bcf_open("-", "r");
            if (bcf_reader == NULL){
                fprintf(stderr, "ERROR interpreting stdin as BCF format.\n");
                exit(1);
            }
            bcf_header = bcf_hdr_read(bcf_reader);
            int num_samples = bcf_hdr_nsamples(bcf_header);
            num_haps = num_samples * 2;
            for (int i = 0; i < num_samples; ++i){
                fprintf(out_haps_f, "%s-1\n", bcf_header->samples[i]);
                fprintf(out_haps_f, "%s-2\n", bcf_header->samples[i]);
                if (depthfile_given){
                    depths_index.push_back(depths[bcf_header->samples[i]]);
                }
            }
            fclose(out_haps_f);
            fprintf(stderr, "Read %d haplotypes\n", num_haps);
            
            // Prepare buffer to store genotypes at all sites
            char geno_str[num_haps+2];
            geno_str[num_haps+1] = '\0';
            geno_str[num_haps] = '\n';
            
            // Array to store genotype data
            int32_t* gts = NULL;
            int n_gts = 0;
            
            const char* gq_key = "GQ";
            
            while(bcf_read(bcf_reader, bcf_header, bcf_record) == 0){
                // Make sure we're on the right chromosome.
                string chrom = bcf_hdr_id2name(bcf_header, bcf_record->rid);
                if (chrom != chromosome){
                    fprintf(stderr, "ERROR: unknown chromosome %s encountered.\n",
                        bcf_hdr_id2name(bcf_header, bcf_record->rid));
                    exit(1);
                }
                
                // REALLY IMPORTANT NOTE: pos here is 0-based, whereas in VCF it's 1-based.
                if (start != -1 && bcf_record->pos + 1 < start){
                    fprintf(stderr, "skipping %d\r", bcf_record->pos + 1);
                    continue;
                }
                if (end != -1 && bcf_record->pos + 1 > end){
                    break;
                }
                
                // Display processing message
                if (bcf_record->pos + 1 - last_printed >= progress){
                    fprintf(stderr, "Processed %s\t%d\r", chrom.c_str(), bcf_record->pos);
                    last_printed = bcf_record->pos;
                }
                
                // Check all sites between the previous and current one where
                // ancestral does not match reference.
                if (bcf_record->pos + 1 > prevpos + 1){
                    for (int ind = prevpos + 1; ind < bcf_record->pos + 1; ++ind){
                        char anc_base = upper(anc_seq->seq.s[ind-1]);
                        char ref_base = upper(ref_seq->seq.s[ind-1]);
                        if (isbase(anc_base) && isbase(ref_base) && anc_base != ref_base){
                            // All haplotypes in this file were homozygous reference
                            // at this site, since it's missing from the file (and it
                            // wasn't filtered out, since neither ref_base or anc_base
                            // is masked to N). Therefore, all haplotypes share a
                            // derived allele at this site.
                            for (int i = 0; i < num_haps; ++i){
                                geno_str[i] = '1';
                            }
                            gzwrite(out_geno_f, geno_str, num_haps+1);
                            fprintf(out_sites_f, "%s\t%d\n", chrom.c_str(), ind);
                            fprintf(out_alleles_f, "%c\t%c\n", anc_base, ref_base);
                        }
                    }
                }
                
                prevpos = bcf_record->pos + 1;
                
                if (bcf_record->n_allele > 2){
                    // Multiallelic site; nothing to do.
                    continue;
                }
                
                /*
                // Filter on variant quality?
                if (bcf_record->qual < quality){
                    continue;
                }
                */
                
                // Load ref/alt alleles and other stuff
                // This puts alleles in bcf_record->d.allele[index]
                // Options for parameter 2:
                /*
                BCF_UN_STR  1       // up to ALT inclusive
                BCF_UN_FLT  2       // up to FILTER
                BCF_UN_INFO 4       // up to INFO
                BCF_UN_SHR  (BCF_UN_STR|BCF_UN_FLT|BCF_UN_INFO) // all shared information
                BCF_UN_FMT  8                           // unpack format and each sample
                BCF_UN_IND  BCF_UN_FMT                  // a synonymo of BCF_UN_FMT
                BCF_UN_ALL (BCF_UN_SHR|BCF_UN_FMT) // everything
                */
                
                bcf_unpack(bcf_record, BCF_UN_STR);
                
                bool indel = false;
                for (int i = 0; i < bcf_record->n_allele; ++i){
                    if (strlen(bcf_record->d.allele[i]) > 1){
                        indel = true;
                        break;
                    }
                }
                if (indel){
                    continue;
                }
                else if (!isbase(bcf_record->d.allele[0][0])){
                    // Reference genome has N
                    continue;
                }
                
                bool flip_haps = false;
                char anc_base = upper(anc_seq->seq.s[bcf_record->pos]);
                
                
                if (!isbase(anc_base)){
                    // Unknown ancestral allele - can't do anything.
                    continue;
                }
                else if (bcf_record->n_allele > 1 && anc_base != bcf_record->d.allele[0][0] && 
                    anc_base != bcf_record->d.allele[1][0]){
                    // Ancestral doesn't match either allele at this site
                    continue;
                }
                
                bool all_ref = false;
                bool all_ancestral = true;
                
                if (bcf_record->n_allele == 1){
                    // All genotypes are homozygous reference. Only print out something
                    // at this site if reference is derived (non-ancestral)
                    if (bcf_record->d.allele[0][0] != anc_base){
                        all_ref = true;
                        all_ancestral = false;
                    }
                    else{
                        // Ancestral base matches reference and all are homozygous
                        // reference.
                        continue;
                    }
                }
                
                if (!all_ref){
                    if (anc_base == bcf_record->d.allele[0][0]){
                        flip_haps = false;
                    }
                    else if (anc_base == bcf_record->d.allele[1][0]){
                        flip_haps = true;
                    }
                }
                
                
                // Determine if all genotypes are valid at this site.
                bool geno_pass = true;
                
                // Look through all genotypes.
                
                int num_loaded = bcf_get_genotypes(bcf_header, bcf_record, &gts, &n_gts);
                if (num_loaded <= 0){
                    fprintf(stderr, "ERROR loading genotypes at %s %d\n", 
                        chrom.c_str(), bcf_record->pos);
                    exit(1);
                }
                
                // Filter out any sites where there is a genotype quality below
                // threshold (OPTIONAL)
                if (quality > 0){
                    float* gqs = NULL;
                    int n_gqs = 0;
                    int num_gq_loaded = bcf_get_format_float(bcf_header, bcf_record, "GQ",
                        &gqs, &n_gqs);
                    if (num_gq_loaded > 0){
                        for (int i = 0; i < num_samples; ++i){
                            if (!isnan(gqs[i]) && gqs[i] != bcf_float_missing &&
                                gqs[i] < quality){
                                geno_pass = false;
                                break;
                            }
                        }
                    }
                    free(gqs);
                }
                // Filter out any sites where there is a depth out of threshold
                // (OPTIONAL)
                if ((depth > 0 || depthfile_given) && geno_pass){
                    int32_t* dps = NULL;
                    int n_dps = 0;
                    int num_dp_loaded = bcf_get_format_int32(bcf_header, bcf_record, "DP",
                        &dps, &n_dps);
                    if (num_dp_loaded > 0){
                        for (int i = 0; i < num_samples; ++i){
                            if (dps[i] != bcf_int32_missing){
                                if (depthfile_given){
                                    if (dps[i] < depths_index[i].first || dps[i] > depths_index[i].second){
                                        geno_pass = false;
                                        break;
                                    }
                                }
                                else if (dps[i] < depth){
                                    geno_pass = false;
                                    break;
                                }
                            }
                        }
                    }
                    free(dps);
                }
                
                if (geno_pass){
                    int ploidy = n_gts / num_samples;
                    
                    for (int i = 0; i < num_samples; ++i){
                        int32_t* gtptr = gts + i*ploidy;
                        for (int j = 0; j < ploidy; ++j){
                            if (gtptr[j] == bcf_int32_vector_end){
                                // Lower ploidy?
                                geno_pass = false;
                                break;
                            }
                            else if (bcf_gt_is_missing(gtptr[j])){
                                // Missing genotype.
                                geno_pass = false;
                                break;
                            }
                            else if (!ignore_phasing && !bcf_gt_is_phased(gtptr[j])){
                                // Unphased genotype?
                                // Check whether it's homozygous (which is phased by definition).
                                bool homozygous = true;
                                for (int k = 0; k < ploidy; ++k){
                                    if (k != j && bcf_gt_allele(gtptr[j]) != bcf_gt_allele(gtptr[k])){
                                        homozygous = false;
                                        break;
                                    }
                                }
                                if (!homozygous){
                                    geno_pass = false;
                                    break;
                                }
                            }
                        
                            if (!all_ref){
                                // Don't bother looking at allele indices if all match
                                // the reference allele; we're just checking if the site
                                // should be skipped.
                                
                                // Retrieve allele index
                                int allele_index = bcf_gt_allele(gtptr[j]);
                                if (flip_haps){
                                    if (allele_index == 0){
                                        geno_str[ploidy*i + j] = '1';
                                        all_ancestral = false;
                                    }
                                    else{
                                        geno_str[ploidy*i + j] = '0';
                                    }
                                }
                                else{
                                    if (allele_index == 0){
                                        geno_str[ploidy*i + j] = '0';
                                    }
                                    else{
                                        geno_str[ploidy*i + j] = '1';
                                        all_ancestral = false;
                                    }
                                }
                            }
                        }
                        if (!geno_pass){
                            break;
                        }
                    }
                }
                
                if (!geno_pass){
                    continue;
                } 
                else{
                    if (all_ref){
                        for (int i = 0; i < num_haps; ++i){
                            geno_str[i] = '1';
                        }
                    }
                    else if (all_ancestral){
                        // Don't bother.
                        continue;
                    }
                    
                    // Write site to disk
                    gzwrite(out_geno_f, geno_str, num_haps+1);
                    fprintf(out_sites_f, "%s\t%d\n", chrom.c_str(), bcf_record->pos + 1);
                    
                    // Write alleles
                    if (all_ref){
                        fprintf(out_alleles_f, "%c\t%c\n", anc_base, bcf_record->d.allele[0][0]);
                    }
                    else{
                        if (flip_haps){
                            fprintf(out_alleles_f, "%c\t%c\n", anc_base, bcf_record->d.allele[0][0]);
                        }
                        else{
                            fprintf(out_alleles_f, "%c\t%c\n", anc_base, bcf_record->d.allele[1][0]);
                        }
                    }
                }
            }
            
            bcf_hdr_destroy(bcf_header);
            bcf_destroy(bcf_record);
            bcf_close(bcf_reader);
            free(gts);
        }
        
        // Break carriage return line
        fprintf(stderr, "\n");
        
        // Now handle all sites between the end of the VCF file and the end of the
        // chromosome.
        if (prevpos < lastpos -1){
            for (int ind = prevpos+1; ind <= lastpos; ++ind){
                char anc_base = upper(anc_seq->seq.s[ind-1]);
                char ref_base = upper(ref_seq->seq.s[ind-1]);
                if (isbase(anc_base) && isbase(ref_base) && anc_base != ref_base){
                    char geno_str[num_haps+2];
                    geno_str[num_haps] = '\n';
                    geno_str[num_haps+1] = '\0';
                    for (int i = 0; i < num_haps; ++i){
                        geno_str[i] = '1';
                    }
                    gzwrite(out_geno_f, geno_str, num_haps+1);
                    fprintf(out_sites_f, "%s\t%d\n", prevchrom.c_str(), ind);
                }
            }
        }
        
        
        // Clean up.
        fclose(out_sites_f);
        fclose(out_alleles_f);
        gzclose(out_geno_f);
        kseq_destroy(anc_seq);
        kseq_destroy(ref_seq);
    
    #endif
    
    return 0;
}
