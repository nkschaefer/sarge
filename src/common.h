// header file for functions common to multiple SPARGE steps

#ifndef SARGECOMMON_H
#define SARGECOMMON_H

#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <string.h>
#include <vector>
#include <zlib.h>
#include <bitset>
#include <set>
#include <unordered_set>
#include <unordered_map>
#include "treeNode.h"
#include "serialize.h"
#include "sets.h"


/*
// http://stackoverflow.com/questions/12314763/efficient-hashing-of-stdbitset-or-boostdynamic-bitset-for-boosts-unor
// http://stackoverflow.com/questions/3896357/unordered-hash-map-from-bitset-to-bitset-on-boost/3897217#3897217

namespace boost {
    template <typename B, typename A>
    std::size_t hash_value(const boost::dynamic_bitset<B, A>& bs) {            
        return boost::hash_value(bs.m_bits);
    }
}

// Store mappings of cladesets to the number of times they appear
typedef boost::unordered_map<boost::dynamic_bitset<unsigned char>, long int> cladesetcounter;
// Store mappings of cladesets to cladesets that are parents of them, and the number of
// times those appear.
typedef boost::unordered_map<boost::dynamic_bitset<unsigned char>, cladesetcounter> cladesetparentcounter;
*/

typedef std::unordered_map<cladeset, long int> cladecounter;

// Data structure to represent a genomic position (chromosome/scaffold and position)
struct loc {
    std::string chrom;
    long int pos;
};
// Allow sorting of locs
bool operator<(const loc&, const loc&);
bool operator==(const loc&, const loc&);

bool cladesetcomp_func(const std::bitset<MAXHAPS>& bs1, const std::bitset<MAXHAPS>& bs2);

// Allow sorting of bitsets (for std::map)
struct cladesetcomp{
    bool operator()(const std::bitset<MAXHAPS>& bs1, const std::bitset<MAXHAPS>& bs2) const{
        return cladesetcomp_func(bs1, bs2);
    }
};

//typedef std::set<cladeset, cladesetcomp> cladeset_set;
typedef std::unordered_set<cladeset> cladeset_set;

typedef std::unordered_set<nodeset> nodeset_set;

// Store information about recombination events from the step1 recomb file
struct recomb_summary{
    long int left;
    long int right;
    
    bool left_known;
    bool joined_known;
    bool moved_known;
    cladeset left_clade;
    cladeset joined_clade;
    cladeset moved_clade;
    
    bool dd_invalid;
    bool modified;
    
    ~recomb_summary(){
        left_clade.reset();
        joined_clade.reset();
        moved_clade.reset();
    }
};


// store any indices of duplicate sites so they can be ignored.
extern std::vector<long int> del_sites;

// Function to sort a vector of bitsets in increasing order of size
bool bitsets_sort_increasing(const cladeset&, const cladeset&);

// Function to sort a vector of bitsets in decreasing order of size
bool bitsets_sort_decreasing(const cladeset&, const cladeset&);

long int load_clade_counts(instream_info&, int, cladecounter&);

unsigned int get_num_haplotypes(gzFile&, int);

int get_num_sites(gzFile&, int, int);

void parse_exclude_bed(std::string, std::map<std::string, std::set<std::pair<long int, long int> > >&);

// Parse the file of loci corresponding to site indices
std::vector<loc> parse_locs(std::string, std::map<std::string, std::set<std::pair<long int, long int> > >&);
    
// Parse a file of individual IDs, map them to haplotype indices.
void parse_indvs(std::vector<std::string>&, std::string);

// Parse a file of individual IDs but insert them into a map of name -> hap index
void parse_indvs_map(std::map<std::string, int>&, std::string);

// Parse a file mapping individual IDs to population IDs
void parse_pops(std::map<std::string, std::string>&, std::map<std::string, std::vector<std::string> >&, std::string, int ploidy);

// Parse a file containing indices of sampled haplotypes and store as a 
// bitset mask
void load_hap_sample(std::string, cladeset&, int);

// Writes out a given site in binary format.
void serialize_sitedata(gzFile&, const std::string&, const long int&, treeNode&);

// Writes out a given recombination event in binary format.
void serialize_recombdata(gzFile&, const std::string&, const long int&, const long int&, bool, bool, bool, 
    cladeset&,
    cladeset&,
    cladeset&,
    int);

      
// Reads in binary-format data for a given site.
void read_sitedata(const int, instream_info&, std::string&, long int&, treeNode&);

// Reads in binary-format data for a given site (ingoring the tree)
void read_posdata(const int, instream_info&, std::string&, long int&);

void unserialize_recombdata(instream_info&, int, std::string&, 
    long int&, long int&, recomb_summary&);

std::string decrement_haps_newick(std::string);
std::string add_brlens(std::string);
std::string replace_scinot_newick(std::string);

// Tree evaluation functions
int parse_ms_file(std::map<std::pair<long int, long int>, treeNode>&, std::string);
void print_ms_tree_site(std::map<std::pair<long int, long int>, treeNode>&, long int pos);
treeNode* get_ms_tree_site(std::map<std::pair<long int, long int>, treeNode>&, long int);

int count_nodes(treeNode*);
int count_nodes_mutations(treeNode*);
int count_haplotypes(treeNode*);
bool has_empty_nodes(treeNode*);

short compare_bitsets(const cladeset& bs1, const cladeset& bs2);

void print_bitset_set(const cladeset&, int);
void print_bitset_set2(const nodeset&, int);

#endif
