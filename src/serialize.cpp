// Includes functions used for serializing data / printing in binary format to stdout.
#include <stdlib.h>
#include <stdio.h>
#include <climits>
#include <iterator>
#include <algorithm>
#include <set>
#include <iostream>
#include <string>
#include <math.h>
#include <cstring>
#include <bitset>
#include <zlib.h>
#include "serialize.h"

using std::cout;
using std::endl;
using namespace std;

// Is the system little endian?
bool little_endian = true;

// How big is a char?
int bytesize = CHAR_BIT * sizeof(unsigned char);

// For portability, store default values we'll use for everything.

// We can get away with two bytes here because we'll just be storing numbers of children
// of nodes & node indices. This gives us up to 32,767 (or 65,535 if unsigned).
int int_bytes = 2;
int longint_bytes = 4;
int float_bytes = 4;

// Thanks to http://stackoverflow.com/questions/4181951/how-to-check-whether-a-system-is-big-endian-or-little-endian
// Store whether the current system is little-endian. We'll represent everything as if
// it's little endian and only flip it if required.
int n = 1;
bool is_little_endian = (*(char *)&n == 1);


void instream_init(instream_info& is, gzFile* instream, int bufsize){
    if (bufsize % CHAR_BIT != 0){
        // Make bufsize a round number of bytes.
        bufsize += (bufsize % CHAR_BIT);
    }
    //is.input = &cin;
    is.input = instream;
    is.bufsize = bufsize;
    is.bytes = new char[bufsize+1];
    is.buffer_index = 0;
    is.eof = false;
    is.initialized = true;
    
    // Do a first read.
    //is.input->read(is.bytes, is.bufsize);
    //is.buf_stored = is.input->gcount();
    int bytes_read = gzread(*is.input, &is.bytes[0], is.bufsize);
    if (bytes_read < bufsize){
        is.eof = true;
    }
    
    is.buf_stored = bytes_read;
}

/**
 * Needed by many other functions. If we've reached too close to the end of the buffer
 * to take in the data we need, we have to copy whatever's left at the end of the buffer
 * to the beginning and fill the rest of the buffer from the input stream. If there's a
 * problem doing this, we'll have to bail out.
 */
void check_bytes_left(int bytes_needed, instream_info& is){
    if ((bytes_needed + is.buffer_index) > is.buf_stored){
        if (gzeof(*is.input)){
        //if (is.input->eof()){
            // We've hit the end of the input stream and the buffer is not full.
            // Unfortunately, asking for the amount of data we need will throw us off the
            // end of not-garbage memory. We need to bail out here.
            //fprintf(stderr, "ERROR: request for data beyond what is in input stream. Aborting.\n");
            //exit(1);
            is.eof = true;
            is.buf_stored = -1;
            return;
        }
        char* read_ptr = &is.bytes[0];
        int readsize = is.bufsize;

        // Copy whatever is left of what we need to the beginning of the input buffer.
        
        if (is.buffer_index < (is.buf_stored)){
            int bytestocopy = is.buf_stored-is.buffer_index;
            memmove(&is.bytes[0], &is.bytes[is.buffer_index], bytestocopy);
            //memcpy(&is.bytes[0], &is.bytes[is.buffer_index], bytestocopy);
            // Get ready for the next file read
            read_ptr += bytestocopy;
            readsize -= bytestocopy;
            is.buf_stored = bytestocopy;
        }
        else{
            is.buf_stored = 0;
        }
        // No matter what we just did, we want to look at the first byte in the array
        // first
        is.buffer_index = 0;
        
        // Read in more from the input buffer.
        //is.input->read(read_ptr, readsize);
        //is.buf_stored += is.input->gcount();
        int bytes_read = gzread(*is.input, read_ptr, readsize);
        if (bytes_read < readsize){
            is.eof = true;
            is.last_bytes_read = bytes_read;
        }
        else{
            is.last_bytes_read = is.bufsize;
        }
        is.buf_stored += bytes_read;
    }    
}

/**
 * Print out a one-byte representation of a boolean variable.
 */
void serialize_bool(gzFile& out, bool value){
     char bytes[1];
    
    if (value){
        //fprintf(out, "%c", 0x10);
        bytes[0] = 0x10;
    }
    else{
        //fprintf(out, "%c", 0x00);
        bytes[0] = 0x00;
    }
    //fwrite(&bytes[0], sizeof(char), 1, out);
    gzwrite(out, &bytes[0], 1);
}

void unserialize_bool(instream_info& is, bool& value){
    check_bytes_left(1, is);
    if (is.bytes[is.buffer_index] == 0x10){
        value = true;
    }
    else if (is.bytes[is.buffer_index] == 0x00){
        value = false;
    }
    else{
        fprintf(stderr, "ERROR: input file corrupt; unable to read boolean value.\n");
        exit(1);
    }
    is.buffer_index++;
}

/**
 * numbytes is where we'll store the number of bytes required to represent our haplotype
 * bitsets in binary format. We modify by reference so we only need to calculate it once
 * (otherwise this would be a waste)
 */

void serialize_bitset(gzFile& out, std::bitset<MAXHAPS> set, long int num_haplotypes){
    
    // How many bytes do we need to represent a bitset? we only want to have to calculate
    // this once.
    
    static int numbytes = (int) ceil((float) num_haplotypes / (float) bytesize);
    char bytes[numbytes + 1];
    bytes[numbytes] = '\0';
    memset(&bytes, 0x00, numbytes);
    
    int byte_index = 0;
    
    static bitset<MAXHAPS> mask;
    for (int i = 0; i < bytesize; ++i){
        mask.set(i);
    }
    
    unsigned char buf[numbytes + 1];
    buf[numbytes] = '\0';
    if (is_little_endian){
        for (int i = 0; i < numbytes; ++i){
            buf[i] = static_cast<unsigned char>((set & mask).to_ulong());
            //buf[i] = (char) set & mask;
            set >>= bytesize;
        }
    }
    else{
        // Go in reverse.
        for (int i = numbytes-1; i >= 0; --i){
            buf[i] = static_cast<unsigned char>((set & mask).to_ulong());
            set >>= bytesize;
        }
    }    
    gzwrite(out, &buf[0], numbytes);
}

void unserialize_bitset(instream_info& is,
    std::bitset<MAXHAPS>& bs, 
    long int num_haplotypes){
    
    bs.reset();
    
    // How many bytes are there in a bitset? We only need to calculate once.
    static int numbytes = (int) ceil((float) num_haplotypes / (float) CHAR_BIT);
    
    // Make sure we are able to read this many bytes from the instream.
    check_bytes_left(numbytes, is);

    int bit_index = 0;
    
    int offset;
    if (!is_little_endian){
        offset = num_haplotypes;
    }
    else{
        offset = -CHAR_BIT;
    }
    
    for (int blockstart = 0; blockstart < numbytes; ++blockstart){
        if (!is_little_endian){
            offset -= CHAR_BIT;
        }
        else{
            offset += CHAR_BIT;
        }
        
        for (int bitnum = 0; bitnum < CHAR_BIT; bitnum++){
            int bitnum_add = CHAR_BIT-1-bitnum;
            if (offset + bitnum_add < num_haplotypes && offset + bitnum_add >= 0){
                if (is.bytes[is.buffer_index + blockstart] << bitnum & 0x80){
                    bs.set(offset + (CHAR_BIT-1-bitnum));
                }
            }
        }
    }
    is.buffer_index += numbytes;
    
}

void serialize_int(gzFile& out, int num){
    const char* buf = reinterpret_cast<const char *>(&num);
    
    char bytes[int_bytes + 1];
    memset(&bytes, 0x00, int_bytes);
    bytes[int_bytes] = '\0';
    int byte_index = 0;
    
    // Store in little-endian order, so flip if the machine is big-endian.
    if (is_little_endian){
        for (int i = 0; i < int_bytes; i++){
            //fprintf(out, "%c", buf[i]);
            bytes[byte_index] = buf[i];
            byte_index++;
        }
    }
    else{
        for (int i = int_bytes-1; i >= 0; i--){
            //fprintf(out, "%c", buf[i]);
            bytes[byte_index] = buf[i];
            byte_index++;
        }
    }
    //fwrite(&bytes[0], sizeof(char), int_bytes, out);
    gzwrite(out, &bytes[0], int_bytes);
}

void unserialize_int(instream_info& is, int& num){
    check_bytes_left(int_bytes, is);
    
    // Make sure there's not extra garbage there in case the system's representation
    // of int is bigger than ours
    memset(&num, 0x00, sizeof(int));
    
    // To do: flip byte order if the system is big-endian.
    memcpy(&num, &is.bytes[is.buffer_index], int_bytes);
    //bytes += int_bytes;
    is.buffer_index += int_bytes;
    
}

void serialize_longint(gzFile& out, long int num){
     const char* buf = reinterpret_cast<const char *>(&num);
    
    char bytes[longint_bytes];
    memset(&bytes, 0x00, longint_bytes);
    
    int byte_index = 0;
    
    // Store in little-endian order, so flip if the machine is big-endian.
    if (is_little_endian){
        for (int i = 0; i < longint_bytes; i++){
            //fprintf(out, "%c", buf[i]);
            bytes[byte_index] = buf[i];
            byte_index++;
        }
    }
    else{
        for (int i = int_bytes-1; i >= 0; i--){
            //fprintf(out, "%c", buf[i]);
            bytes[byte_index] = buf[i];
            byte_index++;
        }
    }
    //fwrite(&bytes[0], sizeof(char), longint_bytes, out);
    gzwrite(out, &bytes[0], longint_bytes);
}

void unserialize_longint(instream_info& is, long int& num){
    check_bytes_left(longint_bytes, is);
    
    memset(&num, 0x00, sizeof(long int));
    
    // To do: flip byte order if the system is big-endian.
    memcpy(&num, &is.bytes[is.buffer_index], longint_bytes);
    //bytes += longint_bytes;
    is.buffer_index += longint_bytes;
    
}

void serialize_float(gzFile& out, float num){
    //unsigned char* buf = reinterpret_cast<unsigned char *>(&num);
    char buf[float_bytes];
    memset(&buf, 0, float_bytes);
    memcpy(&buf, &num, float_bytes);
    
    char bytes[float_bytes];
    memset(&bytes, 0x00, float_bytes);
    int byte_index = 0;
    
    // Store in little-endian order, so flip if the machine is big-endian.
    if (is_little_endian){
        for (int i = 0; i < float_bytes; i++){
            //fprintf(out, "%c", buf[i]);
            bytes[byte_index] = buf[i];
            byte_index++;
        }
    }
    else{
        for (int i = float_bytes-1; i >= 0; i--){
            //fprintf(out, "%c", buf[i]);
            bytes[byte_index] = buf[i];
            byte_index++;
        }
    }
    //fwrite(&bytes[0], sizeof(char), float_bytes, out);
    gzwrite(out, &bytes[0], float_bytes);
}

void unserialize_float(instream_info& is, float& num){
    check_bytes_left(float_bytes, is);
    
    memset(&num, 0x00, sizeof(float));
    
    // To do: flip byte order if the system is big-endian
    memcpy(&num, &is.bytes[is.buffer_index], float_bytes);
    //bytes += float_bytes;
    is.buffer_index += float_bytes;
    
}

void serialize_str(gzFile& out, string str){
    serialize_int(out, (int)str.length());
    //for (int i = 0; i < str.length(); i++){
    //    printf("%c", str[i]);
    //}
    //fprintf(out, "%s", str.c_str());
    if (str.length() > 0){
        const char* cstr = str.c_str();
        //fwrite(cstr, sizeof(char), str.length(), out);
        gzwrite(out, cstr, str.length());
    }
}

void unserialize_str(instream_info& is, string& str){
    // We should have stored the length of the string just before the string itself.
    check_bytes_left(int_bytes, is);
    
    int length;
    //memset(&length, 0x00, sizeof(int));
    //memcpy(&length, &is.bytes[is.buffer_index], int_bytes);
    
    unserialize_int(is, length);
    // Why isn't unserialize_int() doing this for us? This is weird.
    //is.buffer_index += int_bytes;
    
    if (length > 0){
        // Now we can get the string
        check_bytes_left(length, is);
        char cstr[length+1];
        //memset(&cstr, 0x00, length+1);
        memcpy(&cstr, &is.bytes[is.buffer_index], length);
        cstr[length] = '\0';
        
        str = string(cstr);
        
        is.buffer_index += length;
    }
}

/**
 * Prints the file header that this program will require to recognize its own output
 */
void serialize_header(gzFile& out, int num_haplotypes, bitset<MAXHAPS> mask){
    char bytes[11];
    memset(&bytes[0], 0x00, 11);
    sprintf(&bytes[0], "SARGEv0.01");
    gzwrite(out, &bytes[0], 10);
    
    serialize_int(out, num_haplotypes);
    serialize_bitset(out, mask, num_haplotypes);
}

void read_header(instream_info& is, int& num_haplotypes, bitset<MAXHAPS>& mask){
    check_bytes_left(10, is);
    char cstr[11];
    memset(cstr, 0, 11);
    memcpy(&cstr, &is.bytes[is.buffer_index], 10);
    cstr[10] = '\0';
    string header = string(cstr);
    
    is.buffer_index += 10;
    if (header.compare("SPARGEv1.0") == 0){
        unserialize_int(is, num_haplotypes);
    }
    else if (header.compare("SARGEv0.01") == 0 || header.compare("SPARGEv1.1") == 0){
        unserialize_int(is, num_haplotypes);
        unserialize_bitset(is, mask, num_haplotypes);
    }
    else{
        fprintf(stderr, "ERROR: file header missing. Are you sure this was created by SARGE?\n");
        exit(1);
    }
    if (num_haplotypes > MAXHAPS){
        fprintf(stderr, "ERROR: you have compiled sparge to allow for at most %d haplotypes, \
but the input file contains %d haplotypes.\n", MAXHAPS, num_haplotypes);
        fprintf(stderr, "Please recompile as follows:\n");
        fprintf(stderr, "make clean\n");
        fprintf(stderr, "make MAXHAPS=[new value]\n");
        fprintf(stderr, "where [new value] is a number greater than or equal to %d\n", num_haplotypes);
        fprintf(stderr, "Note that large propagation distances may require a greater number.\n");
        exit(1);
    }
}
