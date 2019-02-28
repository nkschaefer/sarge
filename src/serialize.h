// header file for functions related to serialization/binary printing in ARG programs

#ifndef SARGESERIALIZE_H
#define SARGESERIALIZE_H

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <zlib.h>
#include <bitset>

extern int bytesize;
extern int int_bytes;
extern int longint_bytes;
extern int float_bytes;
extern bool is_little_endian;

union intunion{
    int x;
    char buf[sizeof(int)];
};

union longintunion{
    long int x;
    char buf[sizeof(long int)];
};

union floatunion{
    float x;
    char buf[sizeof(float)];
};

struct instream_info{
    //std::istream* input;
    gzFile* input;
    int buffer_index;
    int bufsize;
    int buf_stored;
    int last_bytes_read;
    char* bytes;
    bool eof;
    bool initialized;
    bool finished(){ return ((gzeof(*this->input) || this->eof) && (this->buffer_index == this->buf_stored || this->buf_stored == -1)); };
    void print_bufdata(){ fprintf(stderr, "==========\n"); fprintf(stderr, "bufsize %d\n", this->bufsize); fprintf(stderr, "bytes stored %d\n", this->buf_stored); fprintf(stderr, "buffer index %d\n", this->buffer_index); fprintf(stderr, "==========\n");};
    instream_info(){ this->initialized = false; };
    ~instream_info(){ if (this->initialized){ gzclose(*input); delete this->bytes;} };
    //~instream_info(){ if (this->initialized){ gzclose(*input);} };
};

void instream_init(instream_info&, gzFile*, int);
//void instream_init(instream_info&, std::istream*, int);
void check_bytes_left(int, instream_info&);
void serialize_bool(gzFile&, bool);
void unserialize_bool(instream_info&, bool&);
void serialize_bitset(gzFile&, std::bitset<MAXHAPS>, long int);
void unserialize_bitset(instream_info&, std::bitset<MAXHAPS>&, long int);
void serialize_int(gzFile&, int);
void unserialize_int(instream_info&, int&);
void serialize_longint(gzFile&, long int);
void unserialize_longint(instream_info&, long int&);
void serialize_float(gzFile&, float);
void unserialize_float(instream_info&, float&);
void serialize_str(gzFile&, std::string);
void unserialize_str(instream_info&, std::string&);
void serialize_header(gzFile&, int, std::bitset<MAXHAPS>);
void read_header(instream_info&, int&, std::bitset<MAXHAPS>&);

#endif
