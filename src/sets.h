// header file for functions related to set logic used by ARG programs

#ifndef SARGESET_H
#define SARGESET_H

#include <set>
#include <bitset>

typedef std::bitset<MAXHAPS> cladeset;
typedef std::bitset<NODESETSIZE> nodeset;

/**
 * Prints sets for debugging
 */
void print_set_num(std::set<unsigned int>);
void print_set_num_stdout(std::set<unsigned int>);
void print_set_str(std::set<std::string>);

/**
 * Converts a bitset to a set of ints
 */
std::set<unsigned int> bitset2set(const cladeset&, const unsigned int);

/**
 * Converts a set of ints to a bitset
 */
cladeset set2bitset(std::set<unsigned int>&, const unsigned int);

/**
 * Converts a bitset to a comma-separated string of "on" indices
 */
std::string bitset2str(cladeset, const unsigned int);

/**
 * Randomly subsamples a set number of items from a bitset and returns it as
 * a set of integers representing indices within the bitset.
 */
std::set<unsigned int> subsample_bitset(cladeset&, int, int);
    
/**
 * Performs set intersection
 */
cladeset set_int_bitset(const cladeset&, const cladeset&);
std::set<unsigned int> set_int_num(std::set<unsigned int>&, std::set<unsigned int>&);
std::set<std::string> set_int_str(std::set<std::string>&, std::set<std::string>&);
/**
 * Performs set difference
 */
cladeset set_diff_bitset(const cladeset&, const cladeset&);
nodeset set_diff_nodeset(const nodeset&, const nodeset&);
std::set<unsigned int> set_diff_num(std::set<unsigned int>&, std::set<unsigned int>&);
std::set<std::string> set_diff_str(std::set<std::string>&, std::set<std::string>&);

/**
 * Computes the union of two sets
 */
cladeset set_union_bitset(const cladeset&, const cladeset&);
std::set<unsigned int> set_union_num(std::set<unsigned int>&, std::set<unsigned int>&);
std::set<std::string> set_union_str(std::set<std::string>&, std::set<std::string>&);

/**
 * Returns true if set1 is a superset of (or equal to) set2
 */
bool issuperset_bitset(const cladeset&, const cladeset&);
bool issuperset_nodeset(const nodeset&, const nodeset&);
bool issuperset_num(std::set<unsigned int>&, std::set<unsigned int>&);
bool issuperset_str(std::set<std::string>&, std::set<std::string>&);

/**
 * Returns true if set1 is a subset of (or equal to) set2
 */
bool issubset_bitset(const cladeset&, const cladeset&);
bool issubset_nodeset(const nodeset&, const nodeset&);
bool issubset_num(std::set<unsigned int>&, std::set<unsigned int>&);
bool issubset_str(std::set<std::string>&, std::set<std::string>&);

/**
 * Returns true if set1 and set2 contain the same items
 */
bool seteq_bitset(const cladeset&, const cladeset&);
bool seteq_num(std::set<unsigned int>&, std::set<unsigned int>&);
bool seteq_str(std::set<std::string>&, std::set<std::string>&);

/**
 * Returns true if set contains item.
 */
bool contains_num(std::set<unsigned int>&, unsigned int);
bool contains_str(std::set<std::string>&, std::string&);

/**
 * Performs the four haplotype test on bitsets. Returns true if the test fails, or
 * false if it passes.
 */
bool four_hap_test_bitsets(cladeset&, cladeset&, const unsigned int);
 
/**
 * Performs the four haplotype test on sets rather than bitsets. Returns true
 * if the test fails, or false if it passes.
 */
bool four_hap_test_sets(std::set<unsigned int>&, std::set<unsigned int>&);

/**
 * Compares two bitsets for the purpose of sorting
 */
bool set_comp_bitset(cladeset, cladeset);

/**
 * Compares two sets for the purpose of sorting (numeric)
 */
bool set_comp_num(std::set<unsigned int>, std::set<unsigned int>);

/**
 * Compares two sets for the purpose of sorting (strings)
 */
bool set_comp_str(std::set<std::string>, std::set<std::string>);


#endif
