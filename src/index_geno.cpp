/**
 * Creates an index file and ability to retrieve regions of a .geno.gz file.
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
#include "common.h"
#include <sys/stat.h>
#include <zlib.h>
#include <regex>

using std::cout;
using std::endl;
using namespace std;

// How many bases do we need to go before writing an index entry?
int pos_step = 100000;

/**
 * Function to read in a chunk from the input file.
 * Returns the number of new sites added; modifies hapList, infile, and eof by reference.
 */
long int expand_matrix(gzFile &infile, 
    long int bufsize,
    int& bytes_read,
    const unsigned int num_haplotypes,
    bool &eof,
    bool retrieve_data,
    vector<string>& dat,
    int buf_index){
    
    static long int site_index = 0;
    
    long unsigned int new_sites = 0;
    char* chunk = (char*) malloc(bufsize+1);
    //char chunk[bufsize+1];
    
    int chunk_size;
    chunk_size = gzread(infile, chunk, bufsize);
    bytes_read = chunk_size;
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
    
    if (buf_index != -1){
        // Get rid of stuff before buf_index.
        memmove(&chunk[0], &chunk[buf_index], chunk_size-buf_index);
        chunk[chunk_size-buf_index-1] = '\0';
        //if (chunk_size-buf_index < bufsize+1){
        //    memset(&chunk[chunk_size-buf_index], 0, bufsize-chunk_size-buf_index);
        //}
    }
    
    // Split lines   
    istringstream linesplitter(chunk);
    string line;
    int token_index = 0;

    while(std::getline(linesplitter, line, '\n')){
        if (line.size() == num_haplotypes){
            if (retrieve_data){
                dat.push_back(line);
            }
            //new_sites++;
        }
        new_sites++;
        ++site_index;
    }
    
    // Avoid letting garbage into this buffer next time
    memset(chunk, 0, bufsize);
    
    free(chunk);
    
    //return site_index;
    return new_sites;
}


inline void write_indexdata(FILE* out, 
    string& chrom, 
    long int& pos, 
    int site_index,
    std::streampos strpos,
    int buffer_index){
    long long int strpos_int = (long long int) strpos;
    fprintf(out, "%s\t%ld\t%d\t%lld\t%d\n", chrom.c_str(), pos, site_index, strpos_int, buffer_index);
}

inline void create_index(string indexfilename, string infilename, vector<loc>& sites){
    
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
    
    unsigned int num_haplotypes = get_num_haplotypes(infile, bufsize);
    
    bool eof = false;
    long int site_index = 0;
    string prevchrom = "";
    long int prevpos = 0;
    long long int prevstreampos = (long long int) gztell(infile);
    int prevbufindex = 0;
    
    fprintf(indexfile, "%d\n", bufsize);
    
    vector<string> dat;
    
    while(!eof){
        int bytes_read = 0;
        long int newsites = expand_matrix(infile, bufsize, bytes_read, num_haplotypes, eof, false, dat, -1);
        for (int i = site_index; i < site_index + newsites; ++i){
            if (sites[i].chrom != prevchrom){
                prevpos = 0;
                write_indexdata(indexfile, sites[i].chrom, prevpos, i, prevstreampos, prevbufindex);
                prevchrom = sites[i].chrom;
            }
            else if (sites[i].pos - prevpos > pos_step){
                //write_indexdata(indexfile, sites[i].chrom, sites[i].pos, i, prevstreampos, prevbufindex);
                write_indexdata(indexfile, sites[i].chrom, sites[i].pos, i, prevstreampos-bytes_read, 
                    (i - site_index)*(num_haplotypes+1));
                prevpos = sites[i].pos;
            }
            
            prevstreampos = gztell(infile);

            // Each row takes up num_haplotypes+1 positions in buffer
            //prevbufindex = (i - site_index) * (num_haplotypes+1);
            
        }
        site_index += newsites;
    }
    
    /*
    write_indexdata(indexfile, sites[sites.size()-1].chrom, sites[sites.size()-1].pos,
        sites.size()-2, prevstreampos, prevbufindex);
    */
    fclose(indexfile);
}

void get_region(string indexfilename, 
    string infilename, 
    string chrom, 
    long int start, 
    long int end,
    vector<loc>& sites){
    
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
    int buffer_end = -1;
    int buffer_lastvalid = 0;
    bool in_interval = false;
    int site_index = 0;
    
    int bufsize = -1;
    
    int stream_offset_index = 0;
    
    bool off_the_end = false;
    
    long int pos_start = -1;
    
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
            int cur_site_index = -1;
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
                    cur_site_index = atoi(field.c_str());
                }
                else if (field_index == 3){
                    stream_offset = strtol(field.c_str(), 0, 10);
                }
                else if (field_index == 4){
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
                        site_index = cur_site_index;
                        pos_start = pos;
                    }
                }
                stream_offset_prev = stream_offset;
                buffer_prev = buffer_index;
            }
        }
    }
    
    if (buffer_end == -1){
        off_the_end = true;
    }
    
    indexfile.close();
    
    gzFile infile = gzopen(infilename.c_str(), "rb");
    if (!infile){
        fprintf(stderr, "ERROR: unable to open file %s\n", infilename.c_str());
        exit(1);
    }
    
    unsigned int num_haplotypes = get_num_haplotypes(infile, bufsize);
    
    if (gztell(infile) != stream_offset_start){
        if (stream_offset_start == -1){
            // Go back to last valid chunk we can read and then read on past it.
            gzseek(infile, stream_offset_lastvalid, SEEK_SET);
        }
        else{
            gzseek(infile, stream_offset_start, SEEK_SET);
        }
    }

    bool eof = gzeof(infile);
    bool end_found = false;
    long int base_index = site_index;
    bool start_found = false;
    bool first = true;
    
    while(!end_found){
        vector<string> dat;
        int bytes_read = 0;
        int numsites = expand_matrix(infile, bufsize, bytes_read, num_haplotypes, eof, true, dat, buffer_start);
        
        for (int i = 0; i < dat.size(); ++i){
            if (sites[base_index + i].chrom != chrom || sites[base_index + i].pos > end){
                end_found = true;
                break;
            }
            else if (sites[base_index + i].chrom == chrom && 
                sites[base_index + i].pos >= start && sites[base_index + i].pos <= end){
                fprintf(stdout, "%s\n", dat[i].c_str());
                start_found = true;
            }
            else{
                // Keep going.
            }
        }
        base_index += dat.size();
        if (sites[base_index].chrom == chrom && sites[base_index].pos <= end){
            // Look at next chunk.
            buffer_start = -1;
        }
        if (eof){
            end_found = true;
        }
        first = false;
    }
}

/**
 * Print a help message to the terminal and exit.
 */
void help(int code){
   fprintf(stderr, "index_geno [OPTIONS]\n");
   fprintf(stderr, "Builds an index for quick retrieval of raw genotype data from a .geno.gz\
 file.\n");
    fprintf(stderr, "[OPTIONS]:\n");
    fprintf(stderr, "    --geno -g The .geno.gz file\n");
    fprintf(stderr, "    --sites -s The .sites file\n");
    fprintf(stderr, "    --region -r A region to retrieve in the form chr:start-end\n");
exit(code);
}


int main(int argc, char *argv[]) {   
 
    static struct option long_options[] = {
       {"geno", required_argument, 0, 'g'},
       {"sites", required_argument, 0, 's'},
       {"region", required_argument, 0, 'r'},
       {"help", optional_argument, 0, 'h'},
       {0, 0, 0, 0} 
    };
    
    string genofile;
    string sitesfile;
    string region;
    
    int bufsize = 1048576;
    
    int option_index = 0;
    int ch;
    
    if (argc == 1){
        help(0);
    }
    while((ch = getopt_long(argc, argv, "g:s:r:h", long_options, &option_index )) != -1){
        switch(ch){
            case 0:
                // This option set a flag. No need to do anything here.
                break;
            case 'g':
                genofile = optarg;
                break;
            case 's':
                sitesfile = optarg;
                break;
            case 'r':
                region = optarg;
                break;
            case 'h':
                help(0);
                break;
            default:
                help(0);
        }    
    }

    if (genofile.length() == 0 || sitesfile.length() == 0){
        help(0); 
    }
    
    bool build_index = true;
    bool retrieve = false;
    
    // Determine name of index file to create.
    string indexfilename = genofile + ".index";
    
    // Check if index already exists.
    struct stat stats;
    if (stat(indexfilename.c_str(), &stats) == 0){
        build_index = false;
    }
    else{
        build_index = true;
    }
    
    // We can have an optional second argument, which is the region to retrieve
    // from the input file. Check to see if this is given.
    string chrom;
    long int start;
    long int end;
    
    if (region.length() > 0){
        retrieve = true;

        // The last argument should be of the form chrom:pos1-pos2
        regex searchmatch(R"((.+):([0-9]+)-([0-9]+))", regex::extended);
        cmatch match;
        if (regex_match(region.c_str(), match, searchmatch)){
            chrom = extract_regex_match(region.c_str(), match[0], match[1]);
            string startpos = extract_regex_match(region.c_str(), match[0], match[2]);
            string endpos = extract_regex_match(region.c_str(), match[0], match[3]);
            start = strtol(startpos.c_str(), 0, 10);
            end = strtol(endpos.c_str(), 0, 10);
        }
        else{
            fprintf(stderr, "ERROR: the search string you provided could not be parsed.\n");
            exit(1);
        }
    }
    
    // Parse sites file
    map<string, set<pair<long int, long int> > > excl;
    vector<loc> sites = parse_locs(sitesfile, excl);
    
    if (build_index){
        create_index(indexfilename, genofile, sites);
    }
    if (retrieve){
        // Now we need to load the index into memory.
        get_region(indexfilename, genofile, chrom, start, end, sites);
    }
    if (!build_index && !retrieve){
        fprintf(stderr, "Index already built and no region to retrieve. If you want to rebuild the index, please delete the current file first.\n");
        exit(1);
    }

}
