// defines treeNode class to represent parts of trees for ARG inference.

#ifndef SARGE_TREENODE
#define SARGE_TREENODE

#include <string>
#include <set>
#include <vector>
#include "sets.h"
#include <fstream>
#include <zlib.h>
#include "serialize.h"
#include <regex>
#include <map>
#include <bitset>
#include <unordered_set>

class treeNode{
    private:
        int support;
        bool visited;
        bool subtree_cache;
        cladeset leaves_below;

    public:
        // Variables
        
        // Branch length (from parent clade)
        float dist;
        
        // Branch length (over entire node range within propagation distance)
        float dist2;
    
        
        // Branch length measured in leaving group counts (over entire node range within propagation distance)
        float dist_lv;
        
        float dist_below;
        float dist_above;
        
        float dist_norm;
        
        treeNode* parent;
        cladeset leaves;
        
        long int persistence;
        
        bool prop;
        
        bool keep_children;
        
        std::string name;
        
        bool recomb_left;
        bool recomb_right;
        bool leaving_clade_left;
        bool leaving_clade_right;
        
        bool present_when_loaded;
        bool top_level_prop;
        
        bool recomb_left_subtree;
        bool recomb_right_subtree;
        
        std::vector<treeNode*> children;
        int num_haps;
        
        // If this clade is marked by recombination (either the alpha or beta clade)
        // and its "partner" (the other of the alpha/beta clades in the other tree)
        // is known, this stores its partner.
        std::vector<cladeset > partner_clades_left;
        std::vector<cladeset > partner_clades_right;
        
        std::set<long int> mutations;
        long int num_mutations;
        
        // Sites learned from recombination events, not SNPs
        //std::set<long int> sites;
        
        // Constructor
        treeNode();
        treeNode(const unsigned int, const cladeset&);
        
        // Destructor
        ~treeNode();
        // Copy constructor
        treeNode(const treeNode&);
        // Assignment
        treeNode& operator=(const treeNode&);
        
        // Other functions
        void set_haps(const int);
        void free_recursive();
        void delete_children();
        void set_dist_above(float);
        void setRoot();
        
        
        bool recomb_blocked_subtree(bool, bool, cladeset&);
        bool get_recombs_blocking(bool, bool, cladeset&,
            std::vector<treeNode*>&);
            
        void clear_cache();
        void clear_cache_above();
        
        void addChild(cladeset&, const float);
        void addLeaf(const unsigned int);
        //void addChild(std::set<std::string>&, const float);
        //void addLeaf(const std::string);
        bool remove_leaves(cladeset&);
        bool remove_leaves_prop(cladeset&);
        bool remove_leaves_prop2(cladeset&);
        
        float sum_branchlens_subtree();
        void divide_branchlens_subtree(float);
        float branchlens2percent();
        
        void serialize(gzFile&);
        
        std::string newick(bool, bool, std::vector<std::string>);
        
        std::set<treeNode*> get_lowest();
        
        void rm_leaving_clades(cladeset&, bool, cladeset&);
        void rm_leaving_clades_aux(cladeset&, bool, cladeset&);
        void adjust_for_recomb(cladeset&, bool);
        void get_prop_clades(std::vector<std::pair<cladeset, float> >&, bool);
        
        cladeset subtree_leaves();
        
        treeNode* get_smallest_containing(const cladeset&);
        //void get_smallest_containing(treeNode*&, cladeset&);
        treeNode* get_clade_match(const cladeset&);
        bool has_clade_invalidates(const cladeset&);
        
        bool get_clades_invalidating(const cladeset&,
            std::vector<cladeset >&);
            
        bool get_clade_invalidates(const cladeset&, 
            std::vector<cladeset >&);
        
        void to_dist_map(float, std::map<float, std::vector<treeNode*> >&);
        void flatten(std::vector<cladeset >&);
        void flatten_set(std::unordered_set<cladeset>&);
        
        int clades_correct(treeNode* correct);
        int clades_tot();
        
        std::pair<float, float> mean_children();
        
        void mask_haps(cladeset&);
};

std::set<unsigned int> sample_haplotypes(treeNode&, int, int);

float unserialize_treeNode(treeNode&, const int, instream_info&);
void add_parents(treeNode&);

void parse_subtree(std::vector<std::string>&, std::string&);

bool parse_leaf(std::string, treeNode&, std::string&);

void parse_newick(treeNode&, std::string, const unsigned int);

std::string extract_regex_match(std::string, const std::sub_match<const char*>&, const std::sub_match<const char*>&);


#endif
