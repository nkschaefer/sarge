// Includes functions related to set logic used by other pieces of SARGE.

#include <stdlib.h>
#include <stdio.h>
#include <iterator>
#include <algorithm>
#include <set>
#include <iostream>
#include <string>
#include <ctime>
#include <cstdlib>
#include <bitset>
#include "sets.h"

using std::cout;
using std::endl;
using namespace std;

/**
 * Converts a bitset to a set of indices with value == 1
 */
set<unsigned int> bitset2set(const cladeset& b, const unsigned int num_haplotypes){
    set<unsigned int> s;
    // This is important: bitsets are stored in reverse order from the string
    // that created them. This means that haplotype indices are actually 
    // num_haplotypes - 1 - bit index.
    for (int i=0; i < num_haplotypes; i++){
        if (b[i]){
            s.insert(num_haplotypes-1-i);
        }
    }
    return s;
}

/**
 * Converts a set of indices to a bitset with those indices set to 1
 */
cladeset set2bitset(set<unsigned int>& b, const unsigned int num_haplotypes){
    cladeset s;
    for (set<unsigned int>::iterator it = b.begin(); it != b.end(); ++it){
        s.set(num_haplotypes-(*it)-1);
    }
    return s;
}

std::string bitset2str(cladeset b, const unsigned int num_haplotypes){
    std::string s;
    for (int i=num_haplotypes-1; i >= 0; i--){
        if (b[i]){
            char ind_str[12];
            sprintf(ind_str, "%d", num_haplotypes-1-i);
            s += string(ind_str) + ',';
        }
    }
    // Remove the last comma
    s = s.substr(0, s.size()-1);
    return s;
}



/**
 * Given a dynamic bitset and a number of items to sample from it, randomly chooses
 * that number of samples and returns them as a set of indices..
 */
set<unsigned int> subsample_bitset(cladeset& s, int num_haplotypes,
    int num_to_sample){
    static bool initialized = false;
    if (!initialized){
        // Initialize random number seed
        srand(time(NULL));
    }
    int num_sampled = 0;
    set<unsigned int> sampled;
    int total = s.count();
    while(num_sampled < num_to_sample){
        // This will give us the "rank" in the bitset of the next item to sample.
        int sample_index = rand() % total;
        int rank = 0;
        for (int i = num_haplotypes-1; i >=0; i--){
            if (s[i]){
                // This index is set.
                if (rank == sample_index){
                    // Choose it.
                    sampled.insert(num_haplotypes-1-i);
                    num_sampled++;
                    break;
                }
                rank++;
            }
        }
    }
    return sampled;
}

/**
 * Prints a set to stdout. (int)
 */
void print_set_num(set<unsigned int> s){
    fprintf(stderr, "set([ ");
    for (set<unsigned int>::iterator i = s.begin(); i != s.end(); ++i){
        fprintf(stderr, "%d ", (unsigned int) *i); 
    }
    fprintf(stderr, "])\n");
}
void print_set_num_stdout(set<unsigned int> s){
    fprintf(stdout, "set([ ");
    for (set<unsigned int>::iterator i = s.begin(); i != s.end(); ++i){
        fprintf(stdout, "%d ", (unsigned int) *i); 
    }
    fprintf(stdout, "])\n");
}

/**
 * Prints a set to stdout. (str)
 */
void print_set_str(set<string> s){
    fprintf(stderr, "set([ ");
    for (set<string>::iterator i = s.begin(); i != s.end(); ++i){
        fprintf(stderr, "%s ", (*i).c_str()); 
    }
    fprintf(stderr, "])\n");
}

/**
 * Performs set intersection (bitset)
 */
cladeset set_int_bitset(const cladeset& set1, const cladeset& set2){
    return set1 & set2;
}

/**
 * Performs set intersection (int)
 */
set<unsigned int> set_int_num(set<unsigned int>& set1, set<unsigned int>& set2){
    set<unsigned int> intersection;
    set_intersection(set1.begin(), set1.end(), set2.begin(), set2.end(), 
        inserter(intersection, intersection.begin()));
    return intersection;
}

/**
 * Performs set intersection (str)
 */
set<string> set_int_str(set<string>& set1, set<string>& set2){
    set<string> intersection;
    set_intersection(set1.begin(), set1.end(), set2.begin(), set2.end(), 
        inserter(intersection, intersection.begin()));
    return intersection;
}

/**
 * Performs set difference (bitsets)
 */
cladeset set_diff_bitset(const cladeset& set1, const cladeset& set2){
    return set1 & ~set2;
}

nodeset set_diff_nodeset(const nodeset& set1, const nodeset& set2){
    return set1 & ~set2;
}

/**
 * Performs set difference (int)
 */
set<unsigned int> set_diff_num(set<unsigned int>& set1, set<unsigned int>& set2){
    set<unsigned int> difference;
    set_difference(set1.begin(), set1.end(), set2.begin(), set2.end(),
        inserter(difference, difference.begin()));
    return difference;
}

/**
 * Performs set difference (str)
 */
set<string> set_diff_str(set<string>& set1, set<string>& set2){
    set<string> difference;
    set_difference(set1.begin(), set1.end(), set2.begin(), set2.end(),
        inserter(difference, difference.begin()));
    return difference;
}

/**
 * Computes the union of two sets (bitsets)
 */
cladeset set_union_bitset(const cladeset& set1, const cladeset& set2){
    return set1 | set2;
}

/**
 * Computes the union of two sets (int)
 */
set<unsigned int> set_union_num(set<unsigned int>& set1, set<unsigned int>& set2){
    set<unsigned int> setunion;
    set_union(set1.begin(), set1.end(), set2.begin(), set2.end(),
        inserter(setunion, setunion.begin()));
    return setunion;
}

/**
 * Computes the union of two sets (str)
 */
set<string> set_union_str(set<string>& set1, set<string>& set2){
    set<string> setunion;
    set_union(set1.begin(), set1.end(), set2.begin(), set2.end(),
        inserter(setunion, setunion.begin()));
    return setunion;
}

/**
 * Returns true of set1 is a superset of (or equal to) set2 (bitsets)
 */
bool issuperset_bitset(const cladeset& set1, const cladeset& set2){
    if (set1.count() < set2.count()){
        return false;
    }
    else{
        return set_int_bitset(set1, set2).count() == set2.count();
    }
}

bool issuperset_nodeset(const nodeset& set1, const nodeset& set2){
    if (set1.count() < set2.count()){
        return false;
    }
    else{
        return (set1 & set2).count() == set2.count();
        //return set_int_bitset(set1, set2).count() == set2.count();
    }
}

/**
 * Returns true if set1 is a superset of (or equal to) set2 (num)
 */
bool issuperset_num(set<unsigned int>& set1, set<unsigned int>& set2){
    if (set1.size() < set2.size()){
        return false;
    }
    else{
        return set_int_num(set1, set2).size() == set2.size();
    }
}

/**
 * Returns true if set1 is a superset of (or equal to) set2 (str)
 */
bool issuperset_str(set<string>& set1, set<string>& set2){
    if (set1.size() < set2.size()){
        return false;
    }
    else{
        return set_int_str(set1, set2).size() == set2.size();
    }
}

/**
 * Returns true if set1 is a subset of (or equal to) set2 (bitsets)
 */
bool issubset_bitset(const cladeset& set1, const cladeset& set2){
    if (set1.count() > set2.count()){
        return false;
    }
    else{
        return set_int_bitset(set1, set2).count() == set1.count();
    }
}

bool issubset_nodeset(const nodeset& set1, const nodeset& set2){
    if (set1.count() > set2.count()){
        return false;
    }
    else{
        return (set1 & set2).count() == set1.count();
        //return set_int_bitset(set1, set2).count() == set1.count();
    }
}

/**
 * Returns true if set1 is a subset of (or equal to) set2 (int)
 */
bool issubset_num(set<unsigned int>& set1, set<unsigned int>& set2){
    if (set1.size() > set2.size()){
        return false;
    }
    else{
        return set_int_num(set1, set2).size() == set1.size();
    }
}

/**
 * Returns true if set1 is a subset of (or equal to) set2 (str)
 */
bool issubset_str(set<string>& set1, set<string>& set2){

    if (set1.size() > set2.size()){
        return false;
        
    }
    else{
        return set_int_str(set1, set2).size() == set1.size();
    }
}

/** 
 * Returns true if set1 and set2 contain the same items (bitsets)
 */
bool seteq_bitset(const cladeset& set1, const cladeset& set2){
    return set1 == set2;
}

/**
 * Returns true if set1 and set2 contain the same items (int)
 */
bool seteq_num(set<unsigned int>& set1, set<unsigned int>& set2){
    if (set1.size() != set2.size()){
        return false;
    }
    else{
        return set_int_num(set1, set2).size() == set1.size();
    }
}

/**
 * Returns true if set1 and set2 contain the same items (str)
 */
bool seteq_str(set<string>& set1, set<string>& set2){
    if (set1.size() != set2.size()){
        return false;
    }
    else{
        return set_int_str(set1, set2).size() == set1.size();
    }
}

/**
 * Returns true if set contains item. (int)
 */
bool contains_num(set<unsigned int>& s, unsigned int member){
    if (s.find(member) != s.end()){
        return true;
    }
    return false;
}

/**
 * Returns true if set contains item. (str)
 */
bool contains_str(set<string>& s, string& member){
    if (s.find(member) != s.end()){
        return true;
    }
    return false;
}

/**
 * Performs the four haplotype test on two bitsets. Returns true if the test
 * fails, or false if it passes.
 *
 * NO LONGER USED; use compare_bitsets() in sparge_common.cpp because it returns
 * all information about the relationship of the two (subset/superset/no intersection/
 * four hap test failure).
 *
 */
bool four_hap_test_bitsets(cladeset& a, 
    cladeset& b, const unsigned int num_haplotypes){
    bool het1seen = false;
    bool het2seen = false;
    bool homseen = false;
    for (unsigned int i=0; i < num_haplotypes; i++){
        if (!het1seen && a[i] && !b[i]){
            het1seen = true;
        }
        else if (!het2seen && !a[i] && b[i]){
            het2seen = true;
        }
        else if (!homseen && a[i] && b[i]){
            homseen = true;
        }
        if (het1seen && het2seen && homseen){
            return true;
        }
    }
    return false;
}

/**
 * Performs the four haplotype test on sets rather than bitsets. Returns true
 * if the test fails, or false if it passes.
 */
bool four_hap_test_sets(set<unsigned int>& set1, set<unsigned int>& set2){
    if (set_int_num(set1, set2).size() == 0){
        return false;
    }
    else if (issubset_num(set1, set2) || issubset_num(set2, set1)){
        return false;
    }
    else {
        return true;
    }  
}

/**
 * Compares two bitsets for the purpose of sorting
 */
bool set_comp_bitset(cladeset set1, cladeset set2){
    return set1.count() > set2.count();
}

/**
 * Compares two sets for the purpose of sorting (numeric)
 */
bool set_comp_num(set<unsigned int> set1, set<unsigned int> set2){
    return set1.size() > set2.size();
}

/**
 * Compares two sets for the purpose of sorting (strings)
 */
bool set_comp_str(set<string> set1, set<string> set2){
    return set1.size() > set2.size();
}

