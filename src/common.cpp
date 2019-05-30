/**
 * Contains functions used by multiple pieces of the SPARGE pipeline.
 */
#include <zlib.h>
#include <math.h>
#include <string>
#include <algorithm>
#include <vector>
#include <deque>
#include <iterator>
#include <map>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <math.h>
#include <zlib.h>
#include <regex>
#include <bitset>
#include "sets.h"
#include "treeNode.h"
#include "serialize.h"
#include "common.h"

// This code allows us to use a callback function with C++11 std::regex_replace(),
// instead of needing to use Boost. Thanks to John Martin's suggestion on
// https://stackoverflow.com/questions/22617209/regex-replace-with-callback-in-c11

namespace std {
    template<class BidirIt, class Traits, class CharT, class UnaryFunction>
    std::basic_string<CharT> regex_replace(BidirIt first, BidirIt last,
        const std::basic_regex<CharT,Traits>& re, UnaryFunction f){
        std::basic_string<CharT> s;

        typename std::match_results<BidirIt>::difference_type
            positionOfLastMatch = 0;
        auto endOfLastMatch = first;

        auto callback = [&](const std::match_results<BidirIt>& match)
        {
            auto positionOfThisMatch = match.position(0);
            auto diff = positionOfThisMatch - positionOfLastMatch;

            auto startOfThisMatch = endOfLastMatch;
            std::advance(startOfThisMatch, diff);

            s.append(endOfLastMatch, startOfThisMatch);
            s.append(f(match));

            auto lengthOfMatch = match.length(0);

            positionOfLastMatch = positionOfThisMatch + lengthOfMatch;

            endOfLastMatch = startOfThisMatch;
            std::advance(endOfLastMatch, lengthOfMatch);
        };

        std::regex_iterator<BidirIt> begin(first, last, re), end;
        std::for_each(begin, end, callback);

        s.append(endOfLastMatch, last);

        return s;
    }

    template<class Traits, class CharT, class UnaryFunction>
    std::string regex_replace(const std::string& s,
        const std::basic_regex<CharT,Traits>& re, UnaryFunction f)
    {
        return regex_replace(s.cbegin(), s.cend(), re, f);
    }

} // namespace std

using namespace std;

// Comparison operator for locs so we can push them into a std::map.

bool operator<(const loc& loc1, const loc& loc2){
    int comp = loc1.chrom.compare(loc2.chrom);
    //string::size_type strcmp = loc1.chrom.compare(loc2.chrom);
    if (comp < 0){
        return true;
    }
    else if (comp > 0){
        return false;
    }
    else{
        if (loc1.pos < loc2.pos){
            return true;
        }
        else{
            return false;
        }
    }
}

bool operator==(const loc& loc1, const loc& loc2){
    if (loc1.chrom.compare(loc2.chrom) == 0 && loc1.pos == loc2.pos){
        return true;
    }
    else{
        return false;
    }
}


/**
 * NOTE: TRY TO NOT TO USE THIS! 
 */
bool cladesetcomp_func(const bitset<MAXHAPS>& bs1, const bitset<MAXHAPS>& bs2){

    long int shiftbits = std::min(sizeof(unsigned long)*CHAR_BIT, bs1.size());
    
    //static boost::dynamic_bitset<unsigned char> mask(bs1.size());
    static cladeset mask;
    if (mask.count() == 0){
        // initialize one time.
        for (int i = 0; i < shiftbits; ++i){
            mask.set(i);
        }
    }
    
    long unsigned int num_shifted = 0;
    
    // Copy both bitsets so we can destroy them.
    std::bitset<MAXHAPS> bs1cpy = bs1;
    std::bitset<MAXHAPS> bs2cpy = bs2;
    
    while (num_shifted < bs1.size()){
        // Get and compare the next ulong.
        unsigned long l1 = (bs1cpy & mask).to_ulong();
        unsigned long l2 = (bs2cpy & mask).to_ulong();
        if (l1 < l2){
            return true;
        }
        else if (l2 < l1){
            return false;
        }
        // Move on to the next if they were equal.
        bs1cpy <<= shiftbits;
        bs2cpy <<= shiftbits;
        num_shifted += shiftbits;
    }
    
    // They must be equal if we've come this far.
    return false;
}

vector<long int> del_sites;

bool bitsets_sort_increasing(const cladeset& b1,
    const cladeset& b2){
    if (b1.count() < b2.count()){
        return true;
    }
    else{
        return false;
    }
}

bool bitsets_sort_decreasing(const cladeset& b1,
    const cladeset& b2){
    if (b1.count() > b2.count()){
        return true;
    }
    else{
        return false;
    }
}

/**
 * Given a file containing pre-computed counts of cladesets, loads the file
 * into the given data structure.
 *
 * Returns the total for all cladeset counts in the file.
 *
 */
long int load_clade_counts(instream_info& is, 
    int num_haplotypes,
    cladecounter& clade_counts){
    
    long int tot = 0;
    
    while (!is.finished()){
        cladeset bitstr;
        long int count;
        unserialize_bitset(is, bitstr, num_haplotypes);
        unserialize_longint(is, count);
        
        clade_counts[bitstr] = count;
        tot += count;
    }
    
    return tot;
}

/**
 * In the beginning, read in a chunk of the input file using the specified buffer
 * size. Find the first line in it and from it, determine the number of haplotypes
 * in the input file. Then rewind the gzFile so it can be re-read from the start.
 *
 * Having this number also tells us the number of bytes in a line in the file
 * (num_haplotypes + 1), so we can then resize our buffer to read in appropriate
 * numbers of bytes so we always end on a line break. This way, we'll never
 * have to gzseek(), which can be a very slow operation.
 */
unsigned int get_num_haplotypes(gzFile &infile, int bufsize){
    bool eof = false;
    
    while(!eof){
        char chunk[bufsize];
        int chunk_size;
        chunk_size = gzread(infile, chunk, bufsize-1);
        chunk[chunk_size] = '\0';
        eof = gzeof(infile);
        for (int i=0; i < chunk_size; i++){
            if (chunk[i] == '\n'){
                gzrewind(infile);
                if (i > MAXHAPS){
                    fprintf(stderr, "ERROR: you have compiled sparge to allow for at most %d haplotypes, \
but the input file contains %d haplotypes.\n", MAXHAPS, i);
                    fprintf(stderr, "Please recompile as follows:\n");
                    fprintf(stderr, "make clean\n");
                    fprintf(stderr, "make MAXHAPS=[new value]\n");
                    fprintf(stderr, "where [new value] is a number greater than or equal to %d\n", i);
                    fprintf(stderr, "Note that large propagation distances may require a greater number.\n");
                    exit(1);
                }
                return i;
            }
        }
        // Line break not found. Try again, reading in twice as much data this
        // time.
        gzrewind(infile);
        bufsize *= 2;
    }
    
    // Not found. Abort.
    fprintf(stderr, "ERROR: the number of haplotypes could not be determined from the input file.\n");
    exit(1);
}

/**
 * Count the number of lines (sites) in the input file.
 */
int get_num_sites(gzFile& infile, int bufsize, int num_haplotypes){
    
    int numsites = 0;
    
    bool eof = false;
    
    char* chunk = (char*) malloc(bufsize+1);
    int chunk_size;
    
    eof = gzeof(infile);
    
    int intervening_chars = 0;
    
    while (!eof){
        
        chunk_size = gzread(infile, chunk, bufsize);
        chunk[chunk_size] = '\0';
        
        eof = gzeof(infile);
        
        for (int i = 0; i < chunk_size; ++i){
            if (chunk[i] == '\n'){
                if (intervening_chars == num_haplotypes){
                    numsites++;
                }
                intervening_chars = 0;
            }
            else{
                intervening_chars++;
            }
        }
        
        // Avoid letting garbage into this buffer next time
        memset(chunk, 0, bufsize);
    }
    
    gzrewind(infile);
    free(chunk);
    
    return numsites;
}

void parse_exclude_bed(string filename, map<string, set<pair<long int, long int> > >& exclude_bed){
    std::ifstream bedfile;
    bedfile.open(filename.c_str(), std::ifstream::in);
    if (!bedfile.good()){
        fprintf(stderr, "ERROR: could not open file %s for reading\n", filename.c_str());
        exit(1);
    }
    else{
         for (string line; getline(bedfile, line);){
             // Split by tab
            istringstream tabsplitter(line);
            string field;
            int field_index = 0;
            
            string chrom;
            long int start;
            long int end;
            
            while(std::getline(tabsplitter, field, '\t')){
                if (field_index == 0){
                    chrom = field;
                }
                else if (field_index == 1){
                    start = strtol(field.c_str(), 0, 10);
                }
                else if (field_index == 2){
                    end = strtol(field.c_str(), 0, 10);
                }
                field_index++;
            }
            if (exclude_bed.count(chrom) == 0){
                set<pair<long int, long int> > s;
                exclude_bed.insert(make_pair(chrom, s));
            }
            exclude_bed[chrom].insert(make_pair(start, end));
         }
    }
}

/**
 * Parse the file of position indices.
 */
vector<loc> parse_locs(string filename, map<string, set<pair<long int, long int> > >& exclude_bed){
    vector<loc> locs_parsed;
   
    fstream fin;
    fin.open(filename.c_str(), fstream::in); // open a file
    if (!fin.good()){
        fprintf(stderr, "ERROR opening file %s\n", filename.c_str());
        exit(1);
    }
    string line;
    
    string prevchrom = "";
    long int prevpos = 0;
    
    long int line_index = 0;
    
    set<pair<long int, long int> >::iterator bed_pos;
    bool bed_exists = false;
    
    while (!fin.eof()){
        // Read line
        std::getline(fin, line);
        
        if (line.compare("") != 0){
            // Split by tab
            istringstream tokenizer(line);
            string field;
            loc data;
            int token_index = 0;
        
            while(std::getline(tokenizer, field, '\t')){
                if (token_index == 0){
                    data.chrom = field;
                }
                else{
                    data.pos = atol(field.c_str());
                }
                token_index++;
            }
            
            if (data.chrom.compare(prevchrom) != 0 && exclude_bed.count(data.chrom) > 0){
                // Get an iterator to the appropriate place in the BED file.
                bed_pos = exclude_bed[data.chrom].begin();
                bed_exists = true;
            }
            
            // Ensure all positions are unique.
            // If a duplicate site is found, it will be stored here so that it can be
            // set to all-ancestral when it is read.
            if (data.chrom.compare(prevchrom) == 0 && data.pos == prevpos){
                //fprintf(stderr, "DEL SITE %s\t%ld\n", data.chrom.c_str(), data.pos);
                //exit(1);
                del_sites.push_back(line_index);
            }
            else if (bed_exists){
                // Advance to correct place in BED file.
                while(bed_pos->second < data.pos && bed_pos != exclude_bed[data.chrom].end()){
                    ++bed_pos;
                }
                if (bed_pos != exclude_bed[data.chrom].end()){
                    if (data.pos > bed_pos->first && data.pos <= bed_pos->second){
                        del_sites.push_back(line_index);
                    }
                }
            }
            
            // Store in vector
            locs_parsed.push_back(data);
            
            prevchrom = data.chrom;
            prevpos = data.pos;
            line_index++;
        }
    }

    return locs_parsed;
}

/**
 * Parse the list of individuals.
 **/
void parse_indvs(vector<string>& indvs, string filename){
    fstream fin;
    fin.open(filename.c_str(), fstream::in); // open a file
    if (!fin.good()){
        fprintf(stderr, "ERROR opening file %s\n", filename.c_str());
        exit(1);
    }
    string line;
    
    while (!fin.eof()){
        // Read line
        std::getline(fin, line);
        if (line.compare("") != 0){
            // Store in vector
            indvs.push_back(line);
        }
    }
}

/**
 * Parse the list of individuals, but return a mapping of name -> haplotype index,
 * rather than the other way around
 */
void parse_indvs_map(map<string, int>& indvs, string filename){
    fstream fin;
    fin.open(filename.c_str(), fstream::in); // open a file
    if (!fin.good()){
        fprintf(stderr, "ERROR opening file %s\n", filename.c_str());
        exit(1);
    }
    string line;
    
    int hap_index = 0;
    
    while (!fin.eof()){
        // Read line
        std::getline(fin, line);
        if (line.compare("") != 0){
            indvs.insert(make_pair(line, hap_index));
            hap_index++;
        }
    }
}

/**
 * Parses a file mapping individual to population, tab-separated.
 */
void parse_pops(map<string, string>& indv2pop, map<string, vector<string> >& pop2indv, string filename, int ploidy){
    fstream fin;
    fin.open(filename.c_str(), fstream::in);
    if (!fin.good()){
        fprintf(stderr, "ERROR opening file %s\n", filename.c_str());
        exit(1);
    }
    string line;
    while (!fin.eof()){
        std::getline(fin, line);
        if (line.compare("") != 0){
            istringstream tabsplitter(line);
            string field;
            int field_index = 0;
            
            string indvname;
            string popname;
            
            while(std::getline(tabsplitter, field, '\t')){
                if (field_index == 0){
                    indvname = field;
                }
                else if (field_index == 1){
                    popname = field;
                }
                field_index++;
            }
            if (ploidy == 0){
                indv2pop.insert(make_pair(indvname, popname));
                if (pop2indv.count(popname) == 0){
                    vector<string> indvs;
                    pop2indv.insert(make_pair(popname, indvs));
                }
                pop2indv[popname].push_back(indvname);
            }
            else{
                for (int i = 0; i < ploidy; i++){
                    char hapname[line.length() + 2 + 4];
                    sprintf(hapname, "%s-%d", indvname.c_str(), i+1);
                    string hapname_str = string(hapname);
                    indv2pop.insert(make_pair(hapname_str, popname));
                    if (pop2indv.count(popname) == 0){
                        vector<string> indvs;
                        pop2indv.insert(make_pair(popname, indvs));
                    }
                    pop2indv[popname].push_back(hapname_str);
                    //pop2indv.insert(make_pair(popname, hapname_str));
                }
            }
        }
    }
}

/**
 * Parses a file mapping haplotype names to branch shortening values
 */
void parse_brshorten(unordered_map<cladeset, float>& brshorten_vals,
    vector<string>& indvs, string& filename, int num_haplotypes){
    fstream fin;
    fin.open(filename.c_str(), fstream::in);
    if (!fin.good()){
        fprintf(stderr, "ERROR opening file %s\n", filename.c_str());
        exit(1);
    }
    string line;
    while (!fin.eof()){
        std::getline(fin, line);
        if (line.compare("") != 0){
            istringstream tabsplitter(line);
            string field;
            int field_index = 0;
            
            string hapname;
            float brshorten;
            
            while(std::getline(tabsplitter, field, '\t')){
                if (field_index == 0){
                    hapname = field;
                }
                else if (field_index == 1){
                    brshorten = atof(field.c_str());
                }
                field_index++;
            }
            // Translate haplotype name into a cladeset
            int hap_ind = -1;
            for (int i = 0; i < indvs.size(); ++i){
                if (indvs[i] == hapname){
                    hap_ind = i;
                    break;
                }
            }
            if (hap_ind != -1){
                set<unsigned int> bs;
                bs.insert(hap_ind);
                cladeset cs = set2bitset(bs, num_haplotypes);
                brshorten_vals.insert(make_pair(cs, brshorten));
            }
        }
    }
}

void unnorm_branchlens(treeNode* node){
    node->dist_norm = node->dist;
    for (vector<treeNode*>::iterator c = node->children.begin(); c != node->children.end(); ++c){
        unnorm_branchlens((*c));
    }
}

void norm_brlens(treeNode* tree, float denom){
    tree->dist_norm = tree->dist_norm / denom;
    for (vector<treeNode*>::iterator c = tree->children.begin(); c != tree->children.end(); ++c){
        norm_brlens(*c, denom);
    }
}

/**
 * Given a file containing indices of samples haplotypes (one per line), loads it
 * as a boost::dynamic_bitset with sampled indices set to 1.
 */
void load_hap_sample(string filename, cladeset& hapsample_mask, int num_haplotypes){
    
    std::ifstream hapsampfile;
    hapsampfile.open(filename.c_str(), std::ifstream::in);
    if (!hapsampfile.good()){
        fprintf(stderr, "ERROR: could not open file %s for reading\n", filename.c_str());
        exit(1);
    }
    for (string line; getline(hapsampfile, line);){
        // Line should consist of a single haplotype index.
        int hap_index = atoi(line.c_str());
        hapsample_mask.set(num_haplotypes-hap_index-1);
    }
    hapsampfile.close();
}

/**
 * Given data about a particular site, writes it out in binary format.
 */
void serialize_sitedata(gzFile& out, const string& chrom, const long int& pos, treeNode& root){

    // Only print data about the chromosome if we haven't already seen it.
    // This will save a ton of space.
    static string prevchrom = "";
    
    if (chrom.compare(prevchrom) != 0){
        // Need to store string data.
        serialize_bool(out, true);
        serialize_str(out, chrom);
    }
    else{
        // No need to store string data.
        serialize_bool(out, false);
    }
    
    serialize_longint(out, pos);
    root.serialize(out);
    
    prevchrom = chrom;
}

/**
 * Loads all data for a particular site from stdin in binary format.
 */
void read_sitedata(const int num_haplotypes, instream_info& is, std::string& chrom, long int& pos, treeNode& tree){
    // Keep track of previous chromosome read, since we only store chromosome names
    // when they change.
    static string prevchrom = "";
    
    bool has_chrom;
    unserialize_bool(is, has_chrom);
    if (has_chrom){
        unserialize_str(is, chrom);
        prevchrom = chrom;
    }
    else{
        chrom = prevchrom;
    }
    
    unserialize_longint(is, pos);

    unserialize_treeNode(tree, num_haplotypes, is);
    
    tree.setRoot();
    
}

/**
 * Loads just chromosome and position data for a particular site from stdin in binary format.
 */
void read_posdata(const int num_haplotypes, instream_info& is, std::string& chrom, long int&pos){
    // Keep track of previous chromosome read, since we only store chromosome names
    // when they change.
    static string prevchrom = "";
    
    bool has_chrom;
    unserialize_bool(is, has_chrom);
    if (has_chrom){
        unserialize_str(is, chrom);
        prevchrom = chrom;
    }
    else{
        chrom = prevchrom;
    }
    unserialize_longint(is, pos);
}

/**
 * Given data about a particular recombination event, writes it to the specified file
 * in binary format.
 */
void serialize_recombdata(gzFile& out, const string& chrom, const long int& pos1,
    const long int& pos2,
    bool alpha_known, bool beta_known, bool leaving_known, 
    cladeset& alpha_clade,
    cladeset& beta_clade,
    cladeset& leaving_clade,
    int num_haplotypes){
    
    // Only print data about the chromosome if we haven't already seen it.
    static string prevchrom = "";
    if (chrom.compare(prevchrom) != 0){
        // Need to store string data.
        serialize_bool(out, true);
        serialize_int(out, chrom.length());
        serialize_str(out, chrom);
    }
    else{
        // No need to store string data.
        serialize_bool(out, false);
    }
    
    serialize_longint(out, pos1);
    serialize_longint(out, pos2);
    
    serialize_bool(out, alpha_known);
    if (alpha_known){
        serialize_bitset(out, alpha_clade, num_haplotypes);
    }
    serialize_bool(out, beta_known);
    if (beta_known){
        serialize_bitset(out, beta_clade, num_haplotypes);
    }
    serialize_bool(out, leaving_known);
    if (leaving_known){
        serialize_bitset(out, leaving_clade, num_haplotypes);
    }
    prevchrom = chrom;
}

void unserialize_recombdata(instream_info& is, int num_haplotypes, string& chrom, 
    long int& pos1, long int& pos2, recomb_summary& rs){
    // Keep track of previous chromosome read, since we only store chromosome names
    // when they change.
    static string prevchrom = "";
    
    bool has_chrom;
    unserialize_bool(is, has_chrom);
    if (has_chrom){
        unserialize_str(is, chrom);
        prevchrom = chrom;
    }
    else{
        chrom = prevchrom;
    }
    unserialize_longint(is, pos1);
    unserialize_longint(is, pos2);
    
    rs.dd_invalid = false;
    rs.modified = false;
    
    
    unserialize_bool(is, rs.left_known);
    
    if (rs.left_known){
        unserialize_bitset(is, rs.left_clade, num_haplotypes);
    }

    unserialize_bool(is, rs.joined_known);
    if (rs.joined_known){
        unserialize_bitset(is, rs.joined_clade, num_haplotypes);
    }

    unserialize_bool(is, rs.moved_known);
    if (rs.moved_known){
        unserialize_bitset(is, rs.moved_clade, num_haplotypes);
    }
}

string decrement_num(const smatch& match){
    string strmatch(match[1]);
    int num = atoi(strmatch.c_str());
    char buf[15];
    
    string brlen(match[2]);
    
    // Note: we are also chopping off branch lengths to individual leaf nodes.
    // In doing so, we avoid creating singleton cladesets when we load trees.
    //sprintf(buf, "%d:%s", num-1, brlen.c_str());
    sprintf(buf, "%d:", num-1);

    return string(buf);
}

/**
 * Given a Newick-format string, finds all haplotype indices within it and decrements
 * them by one, returning a new Newick-format string. This is done so we can 
 * parse MS output and not have indices be too big (MS starts from index 1)
 */
string decrement_haps_newick(string newick){
    //regex hapmatch(R"(([0-9]+):([0-9\.e\-]+))", regex::extended);
    regex hapmatch(R"(([0-9]+):)", regex::extended);
    string out = regex_replace(newick, hapmatch, decrement_num); 
    return out;
}

string add_brlens_hap(const smatch& match){
    string hapmatch(match[1]);
    string sep = match[2];
    return hapmatch + ":0" + sep;
}

string add_brlens(string newick){
    
    // First add branch lengths to parenthesis
    
    regex hapmatch(R"(([0-9]+)(,|\)))", regex::extended);
    string out = regex_replace(newick, hapmatch, add_brlens_hap);
    
    
    int found = out.find("),");
    while(found != string::npos){
        out = out.substr(0, found+1) + ":0" + out.substr(found+1, out.length()-found-1);
        found = out.find("),");
    }
    found = out.find("))");
    while(found != string::npos){
        out = out.substr(0, found+1) + ":0" + out.substr(found+1, out.length()-found-1);
        found = out.find("))");
    }
    
    return out;
}

string replace_scinot(const smatch& match){
    string before_e(match[1]);
    string after_e(match[2]);
    
    size_t decpos = before_e.find(".");
    int digits_before = 0;
    int digits_after = 0;
    if (decpos != string::npos){
        digits_before = decpos;
        digits_after = before_e.length()-decpos-1;
        before_e = before_e.substr(0, decpos) + before_e.substr(decpos+1, before_e.length()-decpos-1);
    }
    
    string decimal;
    int exp = atoi(after_e.c_str());
    if (exp < 0){
        string prefix = ":0.";
        for (int i = 0; i < (-exp)-digits_before; ++i){
            prefix += "0";
        }
        decimal = prefix + before_e;
    }
    else{
        string suffix = "";
        for (int i = 0; i < digits_after; ++i){
            suffix += "0";
        }
        decimal = ":" + before_e + suffix;
        
    }
    return decimal;
}

string replace_scinot_newick(string newick){
    regex scimatch(R"(:([0-9\.]+)e([0-9\-]+))", regex::extended);
    string out = regex_replace(newick, scimatch, replace_scinot);
    return out;
}

float norm_brlens(treeNode* t){
    float child_dist_sum = 0;
    int num_children = 0;
    for (vector<treeNode*>::iterator c = t->children.begin(); c != t->children.end(); ++c){
        child_dist_sum += norm_brlens(*c);
        num_children++;
    }
    if (num_children == 0){
        num_children = 1;
    }
    child_dist_sum /= (float)num_children;
    t->dist_below = child_dist_sum;
    child_dist_sum += t->dist;
    return child_dist_sum;
}

/**
 * Parses an output file from MS and loads all the trees into the provided 
 * data structure.
 */
int parse_ms_file(map<pair<long int, long int>, treeNode>& trees, string filename){
    std::ifstream infile;
    infile.open(filename.c_str(), std::ifstream::in);
    if (!infile.good()){
        fprintf(stderr, "ERROR: could not open file %s for reading\n", filename.c_str());
        exit(1);
    }
    
    // NOTE: This seems counterintuitive, but when escaping square brackets ([]) and
    // braces ({}), you have to escape the first one, but should not escape the second
    // or it will cause an error (reference: https://stackoverflow.com/questions/36596563/regex-match-for-square-brackets)
    //regex treematch(R"(\[([0-9]+)]([\(\),\.:a-zA-Z0-9\-_;]+))", regex::extended);
    //regex treematch("(^\[([0-9]+)]([.]+)$)", regex::extended);
    regex treematch(R"(^\[([0-9]+)](.+)$)", regex::extended);
    
    long int start = 1;
    
    int num_haplotypes = -1;
    
    bool reading_trees = false;
    
    for (string line; getline(infile, line);){
        // Need to get the number of haplotypes -- this can be found in 
        // the first line, which is a copy of the command used to run MS
        if (num_haplotypes == -1){
            istringstream space_splitter(line);
            string field;
            int field_index = 0;
            while(std::getline(space_splitter, field, ' ')){
                if (field_index == 0){
                    // This is the path to MS.
                }
                else if (field_index == 1){
                    num_haplotypes = atoi(field.c_str());
                    break;
                }
                field_index++;
            }
            if (num_haplotypes == -1){
                fprintf(stderr, "ERROR: unable to read the number of haplotypes from the MS output.\n");
                exit(1);
            }
        }
        else{
            // Check that this is a tree line and extract useful information
            cmatch match;
        
            try{

                if (regex_match(line.c_str(), match, treematch)){
                    reading_trees = true;
                
                    string len = extract_regex_match(line.c_str(), match[0], match[1]);
                    string newick = extract_regex_match(line.c_str(), match[0], match[2]);
                    
                    //fprintf(stderr, "%s\n", len.c_str());
                    //fprintf(stderr, "%s\n", newick.c_str());
                    
                    long int len_num = strtol(len.c_str(), 0, 10);
                
                    long int end = start + len_num-1;
                
                    pair<long int, long int> key = make_pair(start, end);
                    treeNode root;
                    trees.insert(make_pair(key, root));
                
                    string newick_decremented = decrement_haps_newick(newick);
                    string newick_finished = replace_scinot_newick(newick_decremented);
                    //fprintf(stderr, "\n");
                    //fprintf(stderr, "%s\n", newick_finished.c_str());
                    
                    parse_newick(trees[key], newick_finished, num_haplotypes);
                    norm_brlens(&trees[key]);
                    trees[key].setRoot();
                    //vector<string> dummy;
                    //printf("%ld\t%ld\t%s\n", key.first, key.second, trees[key].newick(false, false, dummy).c_str());
                
                    start += len_num;
                }
                else if (reading_trees){
                
                    // We've read past the tree part of the file.
                    return num_haplotypes;
                }
            }
            catch (regex_error& e){
                if (e.code() == regex_constants::error_collate){
                    fprintf(stderr, "regex error_collate\n");
                }
                else if (e.code() == regex_constants::error_ctype){
                    fprintf(stderr, "regex error_ctype\n");
                }
                else if (e.code() == regex_constants::error_escape){
                    fprintf(stderr, "regex error_escape\n");
                }
                else if (e.code() == regex_constants::error_backref){
                    fprintf(stderr, "regex error_backref\n");
                }
                else if (e.code() == regex_constants::error_brack){
                    fprintf(stderr, "regex error_bracket\n");
                }
                else if (e.code() == regex_constants::error_paren){
                    fprintf(stderr, "regex error_paren\n");
                }
                else if (e.code() == regex_constants::error_brace){
                    fprintf(stderr, "regex error_brace\n");
                }
                else if (e.code() == regex_constants::error_badbrace){
                    fprintf(stderr, "regex error_badbrace\n");
                }
                else if (e.code() == regex_constants::error_range){
                    fprintf(stderr, "regex error_range\n");
                }
                else if (e.code() == regex_constants::error_space){
                    fprintf(stderr, "regex error_space\n");
                }
                else if (e.code() == regex_constants::error_badrepeat){
                    fprintf(stderr, "regex error_badrepeat\n");
                }
                else if (e.code() == regex_constants::error_complexity){
                    fprintf(stderr, "regex error_complexity\n");
                }
                else if (e.code() == regex_constants::error_stack){
                    fprintf(stderr, "regex error_stack\n");
                }
                exit(1);
            }
        }
   }
   return num_haplotypes;
}

void print_ms_tree_site(map<pair<long int, long int>, treeNode>& ms_trees, long int pos){
    for (map<pair<long int, long int>, treeNode>::iterator ms_it = ms_trees.begin();
        ms_it != ms_trees.end(); ++ms_it){
        
        if (ms_it->first.first > pos){
            break;
        }
        if (ms_it->first.first <= pos && ms_it->first.second >= pos){
            
            vector<string> dummy;
            fprintf(stderr, "%s\n", ms_it->second.newick(false, false, dummy).c_str());
        }
    }
}

treeNode* get_ms_tree_site(map<pair<long int, long int>, treeNode>& ms_trees, long int pos){
    for (map<pair<long int, long int>, treeNode>::iterator ms_it = ms_trees.begin();
        ms_it != ms_trees.end(); ++ms_it){
        
        if (ms_it->first.first > pos){
            break;
        }
        if (ms_it->first.first <= pos && ms_it->first.second >= pos){
            return &(ms_it->second);
        }
    }
    return NULL;
}

/**
 * Recursive function to score a given tree. Scores reflect how similar the tree is to
 * a well-balanced binary tree.
 */
float score_tree(treeNode* n, int level){
    float scoresum = 0;
    float num_splits = n->leaves.count() + n->children.size();
    for (vector<treeNode*>::iterator child = n->children.begin(); child != n->children.end();
        ++child){
        scoresum += score_tree((*child), level+1);
    }
    scoresum += (1/(float)level)*(1/num_splits);
    return scoresum;
}

/**
 * Recursive function to count the number of nodes in a given node's subtree.
 */
int count_nodes(treeNode* n){
    
    /*
    int count = 1;
    if (n->parent == NULL && n->children.size() == 1){
        count = 0;
    }
    */
    int count = 1;
    if (n->subtree_leaves().count() == n->num_haps){
        count = 0;
    }
    
    for (vector<treeNode*>::iterator child = n->children.begin(); child != n->children.end();
        ++child){
        cladeset child_subtree = (*child)->subtree_leaves();
        if (child_subtree.count() > 1){
            count += count_nodes(*child);
        }
    }
    return count;
}

/**
 * Recursive function to count the number of nodes in a given node's subtree, counting only
 * nodes that have mutations tied to them.
 */
int count_nodes_mutations(treeNode* n){
    
    /*
    int count = 1;
    if (n->parent == NULL && n->children.size() == 1){
        count = 0;
    }
    */
    int count = 1;
    if (n->subtree_leaves().count() == n->num_haps || n->mutations.size() == 0){
        count = 0;
    }
    
    for (vector<treeNode*>::iterator child = n->children.begin(); child != n->children.end();
        ++child){
        cladeset child_subtree = (*child)->subtree_leaves();
        if (child_subtree.count() > 1){
            count += count_nodes_mutations(*child);
        }
    }
    return count;
}
/**
 * Recursive function to count the number of haplotypes that exist as leaves anywhere
 * in the given node's subtree.
 */
int count_haplotypes(treeNode* n){
    int count = (int) n->leaves.count();
    for (vector<treeNode*>::iterator child = n->children.begin(); child != n->children.end();
        ++child){
        count += count_haplotypes((*child));
    } 
    return count;
}

/**
 * Recursive function to determine whether there are any "empty" nodes down the
 * tree (nodes with no children and no leaves)
 */
bool has_empty_nodes(treeNode* n){
    if (n->children.size() == 0 && n->leaves.count() == 0){
        return true;
    }
    for (vector<treeNode*>::iterator child = n->children.begin(); child != n->children.end();
        ++child){
        if (has_empty_nodes((*child))){
            return true;
        }
    }
    return false;
}

/**
 * Compares two bitsets. Tells if one is a parent of the other, or if they fail the
 * four haplotype test.
 *
 * Returns 
 *   0 separate cladesets
 *   1 cladeset 1 superset cladeset 2
 *   2 cladeset 2 superset cladeset 1
 *   3 cladesets are equal
 *   -1 four haplotype test failure
 */
short compare_bitsets(const cladeset& bs1, const cladeset& bs2){
    if (bs1.size() == 0 || bs2.size() == 0){
        fprintf(stderr, "ERROR: one or both bitsets are unitialized.\n");
        exit(1);
    }
    long int int_len = set_int_bitset(bs1, bs2).count();
    if (bs1.count() > bs2.count() && int_len == bs2.count()){
        return 1;
    }
    else if (bs2.count() > bs1.count() && int_len == bs1.count()){
        return 2;
    }
    else if (bs1.count() == bs2.count() && int_len == bs1.count()){
        return 3;
    }
    else if (int_len > 0){
        return -1;
    }
    return 0;
}

void print_bitset_set(const cladeset& bs, int numhaps){
    fprintf(stderr, "set([ ");
    for (int i = numhaps-1; i >= 0; --i){
        if (bs.test(i)){
            fprintf(stderr, "%d ", numhaps-1-i);
        }
    }
    fprintf(stderr, "])\n");
}

void print_bitset_set2(const nodeset& bs, int numhaps){
    fprintf(stderr, "set([ ");
    for (int i = 0; i < bs.size(); ++i){
        if (bs.test(i)){
            fprintf(stderr, "%d ", i);
        }
    }
    fprintf(stderr, "])\n");
}

