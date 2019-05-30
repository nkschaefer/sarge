/**
 * Creates an index file for a given file created by SARGE. This index can
 * then be used to retrieve results within a given genomic range.
 */
#include <getopt.h>
#include <argp.h>
#include <zlib.h>
#include <string>
#include <vector>
#include <iterator>
#include <map>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include "treeNode.h"
#include "serialize.h"
#include "common.h"
#include <sys/stat.h>
#include <zlib.h>
#include <regex>

using std::cout;
using std::endl;
using namespace std;

// How many bases do we need to go before writing an index entry?
int pos_step = 5000;
    
inline void write_indexdata(FILE* out, 
    bool haschrom, string& chrom, long int& pos, std::streampos strpos,
    int buffer_index){
    long long int strpos_int = (long long int) strpos;
    fprintf(out, "%s\t%ld\t%lld\t%d\n", chrom.c_str(), pos, strpos_int, buffer_index);
}

inline void create_index(string indexfilename, string infilename){
    
    fprintf(stderr, "building index file for %s\n", infilename.c_str());
    
    int bufsize = 1048576;
    
    FILE* indexfile;
    indexfile = fopen(indexfilename.c_str(), "w");
    if (!indexfile){
        fprintf(stderr, "ERROR: could not open file %s for writing\n", indexfilename.c_str());
        exit(1);
    }
    
    // Open the file and write an index file format that tells us file positions
    // every 100 KB
    gzFile infile = gzopen(infilename.c_str(), "rb");
    if (!infile){
        fprintf(stderr, "ERROR: unable to open file %s\n", infilename.c_str());
        exit(1);
    }
    
    fprintf(indexfile, "%d\n", bufsize);
    
    instream_info is;
    instream_init(is, &infile, bufsize);
    
    int num_haplotypes;
    
    // Read file header.
    cladeset mask;
    read_header(is, num_haplotypes, mask);
    
    string prevchrom = "";
    long int prevpos = 0;
    //long long int prevstreampos = (long long int) bgzf_tell(infile);
    long long int prevstreampos = (long long int) gztell(infile);
    //long long int prevstreampos = infile.tellg();
    int prevbufindex = is.buffer_index;
    
    while(!is.finished()){
        string chrom;
        long int pos;
        treeNode tree;
        read_sitedata(num_haplotypes, is, chrom, pos, tree);
        
        int bytesback = is.bufsize;
        if (is.finished()){
            bytesback = is.last_bytes_read;
        }

        if (chrom.compare(prevchrom) != 0){
            // Write an index entry.
            // Index entries will consist of chrom (if different from before),
            // pos, file position
            prevpos = 0;
            write_indexdata(indexfile, true, chrom, prevpos, prevstreampos - bytesback, prevbufindex);
            prevchrom = chrom;
        }
        else if (pos-prevpos > pos_step){
            // Write an index entry.
            write_indexdata(indexfile, false, chrom, pos, prevstreampos - bytesback, prevbufindex);
            prevpos = pos;
        }
        prevstreampos = gztell(infile);
        prevbufindex = is.buffer_index;
    }
    
    //gzclose(infile);
    fclose(indexfile);
}

void get_region(string indexfilename, string infilename, string chrom, long int start, long int end){
    
    std::ifstream indexfile;
    indexfile.open(indexfilename.c_str(), std::ifstream::in);
    if (!indexfile.good()){
        fprintf(stderr, "ERROR: could not open file %s for reading\n", indexfilename.c_str());
        exit(1);
    }
    
    // First, look through the index file and find the first and last offsets
    // of interest.
    long long int stream_offset_prev = 0;
    int buffer_prev = 0;
    
    long long int stream_offset_start = 0;
    long long int stream_offset_end = 0;
    long long int stream_offset_lastvalid = 0;
    int buffer_start = 0;
    int buffer_end = 0;
    int buffer_lastvalid = 0;
    bool in_interval = false;
    
    int bufsize = -1;
    
    int stream_offset_index = 0;
    
    for (string line; getline(indexfile, line);){
        if (bufsize == -1){
            // We are reading the first line, which gives the buffer size.
            bufsize = atoi(line.c_str());
        }
        else{
            // Split by tab
            istringstream tabsplitter(line);
            string field;
            int field_index = 0;
            
            string curchrom;
            long int pos;
            long long int stream_offset;
            int buffer_index;
            
            while(std::getline(tabsplitter, field, '\t')){
                if (field_index == 0){
                    curchrom = field;
                }
                else if (field_index == 1){
                    pos = strtoul(field.c_str(), 0, 10);
                }
                else if (field_index == 2){
                    stream_offset = strtol(field.c_str(), 0, 10);
                }
                else if (field_index == 3){
                    buffer_index = atoi(field.c_str());
                }
                field_index++;
            }
            
            if (stream_offset != -1){
                stream_offset_lastvalid = stream_offset;
                buffer_lastvalid = buffer_index;
            }
            
            if (stream_offset != stream_offset_prev){
                stream_offset_index++;
            }
            
            if (chrom.compare(curchrom) == 0){
                if (in_interval){
                    if (pos + pos_step >= end){
                        // Finished the interval.
                        in_interval = false;
                        stream_offset_end = stream_offset_prev;
                        buffer_end = buffer_prev;
                        break;
                    }
                }
                else{
                    if (start >= pos && start <= (pos + pos_step)){
                        // Start an interval.
                        stream_offset_start = stream_offset;
                        buffer_start = buffer_index;
                        in_interval = true;
                    }
                }
                stream_offset_prev = stream_offset;
                buffer_prev = buffer_index;
            }
        }
    }
    indexfile.close();
    
    // Open the file and write an index file format that tells us file positions
    // every 100 kb
    FILE *outstream = stdout;
    
    if (outstream == NULL){
        fprintf(stderr, "Error opening stdout for writing\n");
        exit(1);
    }
 
    gzFile out_fp = gzdopen(fileno(outstream), "wb");
    if (!out_fp){
        fprintf(stderr, "ERROR: unable to write gzipped output to stdout.\n");
        exit(1);
    }
    
    gzFile infile = gzopen(infilename.c_str(), "rb");
    if (!infile){
        fprintf(stderr, "ERROR: unable to open file %s\n", infilename.c_str());
        exit(1);
    }
    
    instream_info is;
    instream_init(is, &infile, bufsize);

    // Need to get and serialize num haplotypes first, then rewind the file.
    int num_haplotypes;
    cladeset mask;
    read_header(is, num_haplotypes, mask);
    
    serialize_header(out_fp, num_haplotypes, mask);
    
    if (gztell(infile) != stream_offset_start){
        if (stream_offset_start == -1){
            // Go back to last valid chunk we can read and then read on past it.
            gzseek(infile, stream_offset_lastvalid, SEEK_SET);
            //gzseek(infile, stream_offset_lastvalid-bufsize, SEEK_SET);
            is.buf_stored = gzread(infile, &is.bytes[0], is.bufsize);
            buffer_start = buffer_lastvalid;
        }
        else{
            gzseek(infile, stream_offset_start, SEEK_SET);
            //gzseek(infile, stream_offset_start-bufsize, SEEK_SET);
            is.buf_stored = gzread(infile, &is.bytes[0], is.bufsize);
        }
    }
    is.buffer_index = buffer_start;

    
    // Step 0: get to the right position within the first chunk of data.
    bool start_reached = false;
    bool end_reached = false;
    
    while (!is.finished() && !end_reached){
        string curchrom;
        long int pos;
        treeNode tree;
        read_sitedata(num_haplotypes, is, curchrom, pos, tree);
        if (!start_reached && pos >= start){
            start_reached = true;
        }
        if (pos >= end){
            end_reached = true;
        }
        if (start_reached && (!end_reached || start == end)){
            serialize_sitedata(out_fp, chrom, pos, tree);
        }
    }

    //gzclose(infile);
    gzclose(out_fp);
}

int main(int argc, char *argv[]) {    
    // We require one argument: the name of the file to index.
    if (argc == 1){
        fprintf(stderr, "ERROR: you must provide the name of a file to index.\n");
        exit(1);
    }
    
    bool build_index = true;
    bool retrieve = false;
    // We can have an optional second argument, which is the region to retrieve
    // from the input file. Check to see if this is given.
    string chrom;
    long int start;
    long int end;
    
    // Determine name of index file to create.
    string infilename = string(argv[1]);
    string indexfilename = infilename + ".index";
    
    if (argc > 2){
        retrieve = true;

        // The last argument should be of the form chrom:pos1-pos2
        regex searchmatch(R"((.+):([0-9]+)-([0-9]+))", regex::extended);
        cmatch match;
        if (regex_match(argv[2], match, searchmatch)){
            
            chrom = extract_regex_match(argv[2], match[0], match[1]);
            string startpos = extract_regex_match(argv[2], match[0], match[2]);
            string endpos = extract_regex_match(argv[2], match[0], match[3]);
            start = strtol(startpos.c_str(), 0, 10);
            end = strtol(endpos.c_str(), 0, 10);
            
        }
        else{
            fprintf(stderr, "ERROR: the search string you provided could not be parsed.\n");
            exit(1);
        }
        // Check if index already exists.
        struct stat stats;
        if (stat(indexfilename.c_str(), &stats) == 0){
            build_index = false;
        }
        else{
            build_index = true;
        }
        
    }
    
    if (build_index){
        create_index(indexfilename, infilename);
    }
    if (retrieve){
        // Now we need to load the index into memory.
        get_region(indexfilename, infilename, chrom, start, end);
    }

}
