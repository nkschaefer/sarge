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
#include "kseq.h"
#ifdef HTSLIB
    #include <htslib/vcf.h>
#endif

using std::cout;
using std::endl;
using namespace std;


int main(int argc, char *argv[]) {    
    
    
    #ifndef HTSLIB
        fprintf(stderr, "ERROR: bcf format not supported. Please recompile with HTSLIB (HTSLIB=1) \
and make sure HTSLib is installed on your system and available in the LIBRARY_PATH and LD_LIBRARY_PATH \
environment variables.\n");
        exit(1);
    #endif
    
    // Map sample IDs -> chrom -> histogram of depths
    map<int, map<string, map<int, long int> > > depths;
    
    int num_samples;
    
    // How often should progress messages be printed?
    int progress = 50000;
    int last_printed = 0;
    
    // What are the names of the samples in the file?
    vector<string> samples;
   
    #ifdef HTSLIB
        // Read BCF from stdin.
        bcf_hdr_t* bcf_header;
        bcf1_t* bcf_record = bcf_init();
        htsFile* bcf_reader = bcf_open("-", "r");
        if (bcf_reader == NULL){
            fprintf(stderr, "ERROR interpreting stdin as BCF format.\n");
            exit(1);
        }
        bcf_header = bcf_hdr_read(bcf_reader);
        num_samples = bcf_hdr_nsamples(bcf_header);
        for (int i = 0; i < num_samples; ++i){
            samples.push_back(bcf_header->samples[i]);
            
        }
        fprintf(stderr, "Read %d samples\n", num_samples);
        
        
        while(bcf_read(bcf_reader, bcf_header, bcf_record) == 0){
            
            string chrom = bcf_hdr_id2name(bcf_header, bcf_record->rid);
            
            // Display processing message
            if (bcf_record->pos + 1 - last_printed >= progress){
                fprintf(stderr, "Processed %s\t%d\r", chrom.c_str(), bcf_record->pos);
                last_printed = bcf_record->pos;
            }
            
            bcf_unpack(bcf_record, BCF_UN_STR);
            
            // Check depths
            int32_t* dps = NULL;
            int n_dps = 0;
            int num_dp_loaded = bcf_get_format_int32(bcf_header, bcf_record, "DP",
                &dps, &n_dps);
            if (num_dp_loaded > 0){
                for (int i = 0; i < num_samples; ++i){
                    if (dps[i] != bcf_int32_missing){
                        if (depths.count(i) == 0){
                            map<string, map<int, long int> > m;
                            depths.insert(make_pair(i, m));
                        }
                        if (depths[i].count(chrom) == 0){
                            map<int, long int> m;
                            depths[i].insert(make_pair(chrom, m));
                        }
                        if (depths[i][chrom].count(dps[i]) == 0){
                            depths[i][chrom].insert(make_pair(dps[i], 0));
                        }
                        depths[i][chrom][dps[i]]++;
                    }
                }
            }
            free(dps);
            
        }
        
        bcf_hdr_destroy(bcf_header);
        bcf_destroy(bcf_record);
        bcf_close(bcf_reader);
        
    #endif
    
    // Break carriage return line
    fprintf(stderr, "\n");
    
    // Print all depth information to stdout.
    for (map<int, map<string, map<int, long int> > >::iterator dp_s = depths.begin();
        dp_s != depths.end(); ++dp_s){
        string samp = samples[dp_s->first];
        for (map<string, map<int, long int> >::iterator dp_c = dp_s->second.begin();
            dp_c != dp_s->second.end(); ++dp_c){
            for (map<int, long int>::iterator dp = dp_c->second.begin();
                dp != dp_c->second.end(); ++dp){
                fprintf(stdout, "%s\t%s\t%d\t%ld\n", samp.c_str(),
                    dp_c->first.c_str(), dp->first, dp->second);
            }   
        }
    }
    return 0;
}
