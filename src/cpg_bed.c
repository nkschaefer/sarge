/**
  * Outputs a BED file of CpG sites in a FASTA file.
  *
  * This, run on the reference genome, can be given to sarge as the --exclude_bed option.
  * This can be a good idea because cytosines at CpG sites are more mutable than other
  * bases, which can lead to back mutations, which violate the assumptions of the model.
  */

#include <zlib.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <dirent.h>
#include <ctype.h>
#include <stdint.h>

#ifdef HTSLIB
    #include <htslib/kseq.h>

    KSEQ_INIT(gzFile, gzread);
    /**
     * Function to parse FASTA, including gzipped files.
     */
    inline kseq_t* parse_fasta(char *filename){
        gzFile fp = gzopen(filename, "r");
        if (!fp){
            fprintf(stderr, "ERROR: unable to open file %s\n", filename);
            exit(1);
        } else{
            kseq_t* seq;
            seq = kseq_init(fp);
            return seq;
        }
    }

    /**
     * Function to parse FASTA from stdin.
     */
    inline kseq_t* parse_fasta_stdin(){
        // Useful blog post:
        // http://biowize.wordpress.com/2013/03/05/using-kseq-h-with-stdin/
        FILE *instream = stdin;
        if (instream == NULL){
            fprintf(stderr, "Error opening input stream\n");
            exit(1);
        }
        gzFile fp = gzdopen(fileno(instream), "r");
        if (!fp){
            fprintf(stderr, "ERROR: unable to read FASTA file from stdin.\n");
            exit(1);
        } else{
            kseq_t* seq;
            seq = kseq_init(fp);
            return(seq);
        }
    }
#endif

const char capitalize(char base){
    if (base == 'a'){
        return 'A';
    } else if (base == 'c'){
        return 'C';
    } else if (base == 'g'){
        return 'G';
    } else if (base == 't'){
        return 'T';
    } else if (base == 'n'){
        return 'N';
    }
    return base;
}

/**
 * Finds and prints (in BED format) runs of CpG sites.
 */
void print_cpg(char* id, char* seq, int seqlen){
    int seq_index = 0;
    int prev_char = '\0';
    int cur_char = '\0';
    
    int in_cpg = 0;
    int cpg_start = -1;
    
    for (seq_index; seq_index < seqlen; seq_index++){
        cur_char = capitalize(seq[seq_index]);
        if ((cur_char == 'G' && prev_char == 'C') || (cur_char == 'C' && prev_char == 'G')){
            if (!in_cpg){
                cpg_start = seq_index-1;
                in_cpg = 1;
            }
        } else if (in_cpg){
            // The run of CpG bases ended one base ago. Since BED format requires the
            // "end" index to be the true 0-based end index + 1, we can use the current
            // index as the stopping point.
            
            // Print whatever we had been tracking.
            printf("%s\t%d\t%d\n", id, cpg_start, seq_index);
            in_cpg = 0;
            cpg_start = -1;
        }
        prev_char = cur_char;
    }
    // Handle any CpG sites at the end of the sequence.
    if (in_cpg){
        printf("%s\t%d\t%d\n", id, cpg_start, seqlen);
    }
}

/**
 * Print a help message to the terminal and exit.
 */
void help(int code){
   fprintf(stderr, "cpg_bed [filename/-]\n");
   fprintf(stderr, "Given a FASTA file, finds all CpG sites in all sequences in the file \
and outputs them in BED format. You must either provide the filename of a FASTA file \
(possibly gzipped) as an argument, or - to specify stdin.\n");
   exit(code); 
}

int main(int argc, char *argv[]) {    
    if (argc < 2){
        help(1);
    }
    
    #ifdef HTSLIB
        kseq_t* fasta_in;
        
        if (strcmp(argv[1], "-") == 0){
            // Read from stdin.
            fasta_in = parse_fasta_stdin();    
        } else{
            // Parse file.
            fasta_in = parse_fasta(argv[1]);
        }
        
        int progress = 0;
        
        while(progress = kseq_read(fasta_in) >= 0){
            print_cpg(fasta_in->name.s, fasta_in->seq.s, fasta_in->seq.l);
            
        }
        
        kseq_destroy(fasta_in);
    #else
        fprintf(stderr, "ERROR: HTSLib is required for FASTA parsing. Please install \
HTSLib, ensure it is in your LIBRARY_PATH and LD_LIBRARY_PATH environment variables, \
and recompile SARGE with HTSLIB=1.\n");
        exit(1);
    #endif
    return 0;
}

