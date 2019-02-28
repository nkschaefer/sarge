// Defines treeNode class to represent pieces of trees for ARG inference.

#include <string>
#include <set>
#include <vector>
#include <stdio.h>
#include "treeNode.h"
#include <iostream>
#include <sstream>
#include "sets.h"
#include "serialize.h"
#include <fstream>
#include <regex>
#include <cmath>
#include <bitset>
#include <zlib.h>

using namespace std;

/**
 * Initialize a treeNode object.
 */
 
/**
 * Empty constructor
 */
treeNode::treeNode(){
    persistence = -1;
    support = 0;
    prop = true;
    visited = false;
    subtree_cache = false;
    
    recomb_left = false;
    recomb_right = false;
    leaving_clade_left = false;
    leaving_clade_right = false;
    
    top_level_prop = false;
    present_when_loaded = false;
    
    recomb_left_subtree = false;
    recomb_right_subtree = false;
    
    parent = NULL;
    cladeset leaves;
    cladeset leaves_below;
    num_haps = 0;
    dist = 1;
    
    keep_children = false;
    
    children.reserve(20);
}

/**
 * Initialize with leaves
 */
treeNode::treeNode(const unsigned int num_haplotypes, const cladeset& add_leaves){
    
    persistence = -1;
    support = 0;
    prop = true;
    visited = false;
    subtree_cache = false;
    recomb_left = false;
    recomb_right = false;
    
    leaving_clade_left = false;
    leaving_clade_right = false;
    
    top_level_prop = false;
    present_when_loaded = false;
    
    recomb_left_subtree = false;
    recomb_right_subtree = false;
    
    keep_children = false;
    
    parent = NULL;
 
    this->leaves = add_leaves;
    
    dist = 1;
    num_haps = num_haplotypes;
}

/**
 * Build bitsets once size is known
 */
void treeNode::set_haps(const int num_haplotypes){
    this->num_haps = num_haplotypes;
}

void treeNode::set_dist_above(float dist_down){
    dist_down += this->dist;
    
    this->dist_above = dist_down;
    
    float denom = this->dist_above + this->dist_below;
    if (denom == 0){
        denom = 1;
    }
    this->dist_norm = this->dist / denom;
    
    for (vector<treeNode*>::iterator child = this->children.begin();
        child != this->children.end(); ++child){
        (*child)->set_dist_above(dist_down);
    }
}

/**
 * Make current node the root of the tree
 */
void treeNode::setRoot(){
    
    this->parent = NULL;
    
    this->set_dist_above(this->dist);
}

/**
 * delete a treeNode object.
 */
treeNode::~treeNode(){
    // Remove all leaves
    /**
    this->leaves.clear();
    
    this->leaves_below.clear();
    
    this->children.clear();
    **/
    
    //if (!this->keep_children){
        this->leaves.reset();
        this->leaves_below.reset();
    
        for (vector<cladeset >::iterator pc = this->partner_clades_left.begin();
            pc != this->partner_clades_left.end();){
            pc->reset();
            this->partner_clades_left.erase(pc);
        }
        for (vector<cladeset >::iterator pc = this->partner_clades_right.begin();
            pc != this->partner_clades_right.end();){
            pc->reset();
            this->partner_clades_right.erase(pc);
        }
    
        this->delete_children();
        this->children.clear();
    //}
}

/**
 * Deletes all treeNodes in the tree below a given node, including this node.
 */
void treeNode::free_recursive(){
    for (vector<treeNode*>::iterator child = this->children.begin();
        child != this->children.end(); ){
        (*child)->free_recursive();
        this->children.erase(child);
    }
    delete (this);
}

/**
 * Deletes all treeNodes in the tree below a given node, not including this
 * node.
 */
void treeNode::delete_children(){
    for (vector<treeNode*>::iterator child = this->children.begin();
        child != this->children.end();){
        (*child)->free_recursive();
        //(*child)->delete_children();
        this->children.erase(child);
        //delete *child;
    }
}

/**
 * Copy constructor.
 */
treeNode::treeNode(const treeNode &t){
    
    this->support = t.support;
    this->visited = t.visited;    
    this->num_haps = t.num_haps;

    // Force regeneration of cached leaves in subtree.
    this->subtree_cache = false;
    //set<string> leaves_below;
    
    if (t.num_haps > 0){
        
        this->leaves = t.leaves;
    }
    
    else{
       
    }
    this->persistence = t.persistence;
    this->keep_children = false;
    this->top_level_prop = t.top_level_prop;
    this->prop = t.prop;
    this->dist = t.dist;
    this->parent = t.parent;
    this->recomb_left = t.recomb_left;
    this->recomb_right = t.recomb_right;
    
    this->leaving_clade_left = t.leaving_clade_left;
    this->leaving_clade_right = t.leaving_clade_right;
    
    this->present_when_loaded = t.present_when_loaded;
    
    this->recomb_left_subtree = t.recomb_left_subtree;
    this->recomb_right_subtree = t.recomb_right_subtree;
    
    this->partner_clades_left = t.partner_clades_left;
    this->partner_clades_right = t.partner_clades_right;
    
    this->name = t.name;
    
    // Copy all children as well.
    vector<treeNode*> children;
    for (int i=0; i < t.children.size(); ++i){
        this->children.push_back(new treeNode((*t.children[i])));
    }
}

/**
 * operator= function
 * Given a treeNode and another treeNode, modifies this treeNode so it exactly
 * resembles the given treeNode
 */
treeNode& treeNode::operator=(const treeNode &t){
    // Check for self assignment
    if (this == &t){
        return *this;
    }
    else{
        this->persistence = t.persistence;
        this->keep_children = false;
        this->support = t.support;
        this->visited = t.visited;
        this->dist = t.dist;
        this->prop = t.prop;
        this->top_level_prop = t.top_level_prop;
        
        this->recomb_left_subtree = t.recomb_left_subtree;
        this->recomb_right_subtree = t.recomb_right_subtree;
        
        if (t.num_haps == 0){
            // If this happens, we'll need to set num_haps at a later step. Otherwise,
            // the bitsets will throw errors.
            
            this->num_haps = 0;
        }
        else if (this->num_haps == t.num_haps){
            // No need to re-declare; just re-set.
            this->leaves.reset();
            this->leaves = t.leaves;
            this->leaves_below.reset();
        }
        else{
            
            this->num_haps = t.num_haps;
        }
        
        // Force regeneration of cached leaves in subtree.
        this->subtree_cache = false;
        
        this->recomb_left = t.recomb_left;
        this->recomb_right = t.recomb_right;
        this->leaving_clade_left = t.leaving_clade_left;
        this->leaving_clade_right = t.leaving_clade_right;
        
        this->present_when_loaded = t.present_when_loaded;
        
        this->partner_clades_left = t.partner_clades_left;
        this->partner_clades_right = t.partner_clades_right;
        
        this->name = t.name;
        
        // Erase children
        this->children.clear();
        
        // Copy all children over
        vector<treeNode*> children;
        for (int i=0; i < t.children.size(); ++i){
            this->children.push_back(new treeNode((*t.children[i])));
        }
        
        return *this;
    }
}

/**
 * Given a node representing the root of a tree, recursively visits all that node's 
 * children and sets parent pointers.
 */
void add_parents(treeNode& t){
    for (vector<treeNode*>::iterator child = t.children.begin(); child != t.children.end();
        ++child){
        (*child)->parent = &t;
        add_parents(*(*child));
    }
}

/**
 * Removes cached information about all leaves in this node's subtree.
 */
void treeNode::clear_cache(){
    //this->leaves_below = set<string>();
    this->leaves_below.reset();
    this->subtree_cache = false;
}

/**
 * Adds a child to the tree below this node; notifies that child that this
 * node is now its parent.
 *
 * http://codereview.stackexchange.com/questions/47395/tree-template-class-implementation-for-c
 */
void treeNode::addChild(cladeset& leaves, float brlen){
    treeNode* child = new treeNode(this->num_haps, leaves);
    child->dist = brlen;
    children.push_back(child);
    child->parent = this;
    
    this->leaves = set_diff_bitset(this->leaves, leaves);
    this->clear_cache();
}

void treeNode::addLeaf(const unsigned int leaf){
    this->leaves.set(this->num_haps-1-leaf);
    this->clear_cache();
}

/**
 * Given a set of leaves to eliminate, finds all such leaves in the subtree
 * of this node and removes them. If this removes all leaves from a given
 * node, that node is also removed.
 *
 * returns true if the current node is to be deleted; false otherwise.
 */
bool treeNode::remove_leaves(cladeset& remove){
    
    if (remove.count() > 0){
        for (vector<treeNode*>::iterator child = this->children.begin();
            child != this->children.end();){
            bool remove_child = (*child)->remove_leaves(remove);
            if (remove_child){
                this->children.erase(child);
            }
            else{
                ++child;
            }
        }
    }
    if (remove.count() > 0){
        cladeset overlap = remove & this->leaves;
        if (overlap.count() > 0){
            this->leaves = set_diff_bitset(this->leaves, overlap);
            remove = set_diff_bitset(remove, overlap);
            
            if (this->leaves.count() == 0 && this->children.size() == 0){
                return true;
            }
        }
    }
    this->clear_cache();
    return false;
}

bool treeNode::remove_leaves_prop(cladeset& remove){
    
    cladeset this_subtree = this->subtree_leaves();
    if (issubset_bitset(this_subtree, remove)){
        return true;
    }
    this->leaves = set_diff_bitset(this->leaves, remove);
    for (vector<treeNode*>::iterator child = this->children.begin(); child != this->children.end();){
        if ((*child)->remove_leaves_prop(remove)){
            this->children.erase(child);
            //this->clear_cache();
        }
        else{
            ++child;
        }
    }
    this->clear_cache();
    this_subtree = this->subtree_leaves();
    if (this_subtree.count() < 2){
        return true;
    }
    else{
        return false;
    }
}

bool treeNode::remove_leaves_prop2(cladeset& remove){
    if (remove.count() > 0){
        for (vector<treeNode*>::iterator child = this->children.begin();
            child != this->children.end();){
            bool remove_child = (*child)->remove_leaves(remove);
            if (remove_child){
                this->children.erase(child);
            }
            else{
                ++child;
            }
        }
    }
    if (remove.count() > 0){
        cladeset overlap = remove & this->leaves;
        if (overlap.count() > 0){
            this->leaves = set_diff_bitset(this->leaves, overlap);
            remove = set_diff_bitset(remove, overlap);
            
            if (this->leaves.count() <= 1 && this->children.size() == 0){
                return true;
            }
        }
    }
    this->clear_cache();
    return false;
}

/**
 * Tells whether recombination is blocked in the given direction anywhere
 * in this node's subtree.
 */
bool treeNode::recomb_blocked_subtree(bool fromleft, bool children_only, 
    cladeset& recomb_clade){
    
    if (!children_only && ((fromleft && this->recomb_left) || (!fromleft && this->recomb_right))){
        recomb_clade = this->subtree_leaves();
        
        return true;
    }
    else{
        for (vector<treeNode*>::iterator child = this->children.begin(); child != this->children.end();
            ++child){
            return (*child)->recomb_blocked_subtree(fromleft, false, recomb_clade);
        }
    }
    return false;
}

bool treeNode::get_recombs_blocking(bool fromleft, bool children_only,
    cladeset& prop_clade,
    vector<treeNode*>& blocked_clades){
    
    cladeset this_subtree = this->subtree_leaves();
    
    if (!children_only && 
        (this_subtree & prop_clade).count() > 0 && (
        (fromleft && this->recomb_left) ||
        (!fromleft && this->recomb_right))){ 
        // Blocked.
        blocked_clades.push_back(this);
        return true;
    }
    else{
        for (vector<treeNode*>::iterator child = this->children.begin(); child != this->children.end();
            ++child){
            return (*child)->get_recombs_blocking(fromleft, false, prop_clade, blocked_clades);
        }
    }
    return false;
}

/**
 * Add up all the branch lengths in this node's subtree and return them
 * (as a float).
 */
float treeNode::sum_branchlens_subtree(){
    float brlen_sum = this->dist;
    for (vector<treeNode*>::iterator child = this->children.begin(); child != this->children.end(); ++child){
        brlen_sum += (*child)->sum_branchlens_subtree();
    }
    return brlen_sum;
}

/**
 * Divides the length of every branch in this node's subtree by the given
 * number.
 */
void treeNode::divide_branchlens_subtree(float num){
    if (num > 0){
        this->dist /= num;
        for (vector<treeNode*>::iterator child = children.begin(); child != children.end(); ++child){
            (*child)->divide_branchlens_subtree(num);
        }
    }
    return;
}

/** 
 * This should be called on the root node.
 * Finds every node in this node's subtree and adds up all of their branch
 * lengths. Returns this to the root. Then makes a second pass of the entire tree
 * and divides every node's branch length by this total value. This results in 
 * branch lengths being percentages of total tree depth, rather than raw numbers,
 * to account for varying mutation rates across the genome.
 *  
 *  Returns:
 *      the total depth of the tree (in mutations; float)
 */
float treeNode::branchlens2percent(){
    float depth_tot = this->sum_branchlens_subtree();
    this->divide_branchlens_subtree(depth_tot);
    return depth_tot;
}

/**
 * Return a Newick-format string representation of this treeNode.
 * 
 * Arguments:
 *  recomb -- true or false, whether or not this node's status as marked by
 *  recombination (or not) in either direction should be marked in the Newick string.
 *  use_hapnames -- true or false -- should haplotype indices be converted into names
 *      using the provided (possibly empty) haplotype index -> name mapping?
 *  hapnames -- a mapping of haplotype index -> string name to use in the tree
 */
string treeNode::newick(bool recomb, bool use_hapnames, vector<std::string> hapnames){
    string recomb_left_str = "";
    
    if (recomb && this->leaving_clade_left){
        recomb_left_str = "{";
    }
    
    if (recomb && this->recomb_left){
        recomb_left_str += "[";
    }
    string recomb_right_str = "";
    
    if (recomb && this->leaving_clade_right){
        recomb_right_str = "}";
    }
    
    if (recomb && this->recomb_right){
        recomb_right_str += "]";
    }
    string wrapperL = "";
    string wrapperR = "";
    string root = "";
    if (this->parent == NULL){
        wrapperL = "(";
        wrapperR = ")";
        
        char rootbuf[15];
        sprintf(rootbuf, ":%5.5f;", this->dist_norm);
        root += string(rootbuf);
        
        //root = ";";
    }

    if (this->children.size() == 0){
        // Base case: this is just leaves.
        string newickstr = "";
        //if (false){
        if (leaves.count() == 1){
            set<unsigned int> leaves_set = bitset2set(leaves, this->num_haps);
            if (use_hapnames){
                newickstr = hapnames[*leaves_set.begin()];
            }
            else{
                char leaf_str[20];
                sprintf(leaf_str, "%d", *leaves_set.begin());
                newickstr = string(leaf_str);
            }
        }
        else if (leaves.count() > 0){
            set<unsigned int> leaves_set = bitset2set(leaves, this->num_haps);
            set<unsigned int>::iterator lastit = leaves_set.end();
            lastit--;
            for (set<unsigned int>::iterator leaf = leaves_set.begin(); leaf != leaves_set.end(); ++leaf){
                if (use_hapnames){
                    string leaf_str = hapnames[*leaf] + ":0";
                    newickstr += leaf_str;
                }
                else{
                    char leaf_str[20];
                    sprintf(leaf_str, "%d:0", (*leaf));
                    newickstr += string(leaf_str);
                }
                //newickstr += (*leaf) + ":0";
                if (leaf != lastit){
                    newickstr += ",";
                }
            }
        }
        return (wrapperL + recomb_left_str + newickstr + recomb_right_str + wrapperR + root);
    }
    else{
        
        string children_rendered = "";
        vector<treeNode*>::iterator lastit = this->children.end();
        lastit--;
        for (vector<treeNode*>::iterator child = this->children.begin(); child != this->children.end(); ++child){
            // Note: adding distances here (only when the method is called on a 
            // child node) means that no distance will be given to the root node.
            // That's okay, since the concept of distance doesn't make sense
            // for it.
            string child_newick = (*child)->newick(recomb, use_hapnames, hapnames);
            
            if (this->name.length() == 0){
                char buf[child_newick.size() + 10 + 3];
                
                //if (false){
                if ((*child)->leaves.count() == 1 && (*child)->children.size() == 0){
                    // Child is a leaf node.
                    sprintf(buf, "%s:%5.5f", child_newick.c_str(), abs((*child)->dist_norm));
                }
                else{
                    sprintf(buf, "(%s):%5.5f", child_newick.c_str(), abs((*child)->dist_norm));
                }
                children_rendered += (string) buf;
            }
            else{
                char buf[child_newick.size() + this->name.length() + 10 + 3 + 1];
                
                sprintf(buf, "(%s)%s:%5.5f", child_newick.c_str(), this->name.c_str(), abs((*child)->dist_norm));
            
                children_rendered += (string) buf;
            }
            if (child != lastit || this->leaves.count() > 0){
                children_rendered += ",";
            }
        }
        
        string leaves_rendered = "";
        if (this->leaves.count() > 0){
            set<unsigned int> leaves_set = bitset2set(this->leaves, this->num_haps);
            set<unsigned int>::iterator lastleaf = leaves_set.end();
            //set<string>::iterator lastleaf = this->leaves.end();
            lastleaf--;
            for (set<unsigned int>::iterator leaf = leaves_set.begin(); leaf != leaves_set.end(); ++leaf){
                if (use_hapnames){
                    string leaf_str = hapnames[*leaf];
                    leaves_rendered += leaf_str + ":0";
                }
                else{
                    char leaf_str[12];
                    sprintf(leaf_str, "%d:0", (*leaf));
                    leaves_rendered += string(leaf_str);
                }
                //leaves_rendered += (*leaf) + ":0";
                if (leaf != lastleaf){
                    leaves_rendered += ",";
                }
            }
        }
          
        return (wrapperL + recomb_left_str + children_rendered + leaves_rendered + 
            recomb_right_str + wrapperR + root);
    }
}

/**
 * Return a vector of pointers to all treeNodes representing lowest-level nodes
 * (no children) below this node in the tree.
 */
set<treeNode*> treeNode::get_lowest(){
    set<treeNode*> lowest_nodes;
    if (children.size() == 0){
        lowest_nodes.insert(this);
    }
    else{
        for (vector<treeNode*>::iterator child = this->children.begin(); child != this->children.end(); ++child){
            set<treeNode*> child_lowest_nodes = (*child)->get_lowest();
            for (set<treeNode*>::iterator child_child = child_lowest_nodes.begin();
                child_child != child_lowest_nodes.end(); ++child_child){
                lowest_nodes.insert(*child_child);
            }
        }
    }
    return lowest_nodes;
}

void treeNode::clear_cache_above(){
    this->clear_cache();
    if (this->parent == NULL){
        return;
    }
    else{
        this->parent->clear_cache_above();
    }
}

void treeNode::rm_leaving_clades_aux(cladeset& prop_subtree,
    bool left, 
    cladeset& partner){
    
    vector<string> dummy;
    
    cladeset this_subtree;
    
    if ((left && this->leaving_clade_left) || (!left && this->leaving_clade_right)){
        this_subtree = this->subtree_leaves();
        if (issubset_bitset(this_subtree, partner)){
            prop_subtree = set_diff_bitset(prop_subtree, this_subtree);
            return;
        }
    }
    
    for (vector<treeNode*>::iterator child = this->children.begin(); 
        child != this->children.end(); ++child){
        (*child)->rm_leaving_clades_aux(prop_subtree, left, partner);
        cladeset child_subtree = (*child)->subtree_leaves();
    }
}

/**
 * This is called on the clade affected by recombination. rm_leaving_clades_aux() will
 * be called on its children that might be marked as leaving clades.
 */
void treeNode::rm_leaving_clades(cladeset& prop_subtree,
    bool left, 
    cladeset& partner){
    
    for (vector<treeNode*>::iterator child = this->children.begin(); 
        child != this->children.end(); ++child){
        cladeset prop_subtree_cpy(prop_subtree);
        (*child)->rm_leaving_clades_aux(prop_subtree_cpy, left, partner);
        prop_subtree = set_int_bitset(prop_subtree, prop_subtree_cpy);
    }
}

void treeNode::adjust_for_recomb(cladeset& prop_subtree, 
    bool left){
    
    cladeset this_subtree = this->subtree_leaves();
    
    if (left && this->recomb_left){
        //fprintf(stderr, "TRUE %ld\n", this->partner_clades_left.size());
        for (vector<cladeset >::iterator pc = this->partner_clades_left.begin();
            pc != this->partner_clades_left.end(); ++pc){
            // We already know it's a superset of the current clade (because this method 
            // is only called on children)
            if (issuperset_bitset(prop_subtree, *pc) && prop_subtree.count() > this_subtree.count()){
                // Do nothing for now; propagating clade is safe.
                //fprintf(stderr, "propsafe\n");
            }
            else{
                //fprintf(stderr, "not safe %ld\n", this->children.size());
                // Remove from the propagating node all leaving clades associated with this
                // recombination event.
                this->rm_leaving_clades(prop_subtree, left, *pc); 
            }
        }
    }
    else if (!left && this->recomb_right){
        for (vector<cladeset >::iterator pc = this->partner_clades_right.begin();
            pc != this->partner_clades_right.end(); ++pc){
            if (issuperset_bitset(prop_subtree, *pc) && prop_subtree.count() > this_subtree.count()){
                // Do nothing for now; propagating clade is safe.
            }
            else{
                // Remove from the propagating node all leaving clades associated with this
                // recombination event.
                this->rm_leaving_clades(prop_subtree, left, *pc);     
            }
        }
    }
    
    // Recursively call on children.
    for (vector<treeNode*>::iterator child = this->children.begin(); 
        child != this->children.end(); ++child){
        
        if ((left && (*child)->recomb_left_subtree) || (!left && (*child)->recomb_right_subtree)){
            cladeset prop_subtree_cpy(prop_subtree);
            (*child)->adjust_for_recomb(prop_subtree_cpy, left);
            prop_subtree = set_int_bitset(prop_subtree, prop_subtree_cpy);
        }
    }
    
}

/**
 * Just propagate top-level clades. If one is invalidated, it'll be replaced with its children.
 */
void treeNode::get_prop_clades(vector<pair<cladeset, float> >& nodes, bool left){

    for (vector<treeNode*>::iterator child = this->children.begin();
        child != this->children.end(); ++child){
        
        cladeset child_subtree = (*child)->subtree_leaves();
        
        // Don't propagate singletons.
        if (child_subtree.count() > 1){
            if ((left && (*child)->recomb_left_subtree) || (!left && (*child)->recomb_right_subtree)){
                (*child)->adjust_for_recomb(child_subtree, left);
            }
            
            if (child_subtree.count() > 1){
                //if (true){
                if ((*child)->dist > 0){
                    // Don't propagate 0 branch length clades (the clades that they're 
                    // built on will exist at other indices)
                    nodes.push_back(make_pair(child_subtree, (*child)->dist));
                }
                // Recursively add all children as well
                vector<pair<cladeset, float> > nodes_children;
                (*child)->get_prop_clades(nodes_children, left);
                for (vector<pair<cladeset, float> >::iterator node_children = 
                    nodes_children.begin(); node_children != nodes_children.end(); ++node_children){
                    nodes.push_back(*node_children);
                }
            }
        }
    }
}

/**
 * Return a set of all leaf names (strings) in this node's entire subtree.
 */
cladeset treeNode::subtree_leaves(){
    //if (this->leaves_below.size() == 0){
    //    this->leaves_below = cladeset(this->num_haps);
    //}
    if (subtree_cache){
        // We have a stored version of this set; no need to re-generate.
        return this->leaves_below;
    }
    else{
        // Need to re-generate.
        this->leaves_below = this->leaves;
        for (vector<treeNode*>::iterator child = this->children.begin();
            child != this->children.end(); ++child){
            cladeset child_subtree = (*child)->subtree_leaves();
            this->leaves_below = set_union_bitset(this->leaves_below, child_subtree);
        }
        this->subtree_cache = true;
        return this->leaves_below;
    }
}

/**
 * Given a set of leaves, goes through this treeNode's subtree until the lowest-possible
 * level is found that contains all of the given leaves.
 */
treeNode* treeNode::get_smallest_containing(const cladeset& clade){
    for (vector<treeNode*>::iterator child = this->children.begin(); child != this->children.end();
        ++child){
        cladeset child_subtree = (*child)->subtree_leaves();
        if (issubset_bitset(clade, child_subtree)){
            return (*child)->get_smallest_containing(clade);
            //(*child)->get_smallest_containing(ptr, clade);
            //return;
        }   
    }
    
    // No child has the clade in its subtree. Is it in this clade's subtree?
    cladeset this_subtree = this->subtree_leaves();
    if (issubset_bitset(clade, this_subtree)){
        //ptr = this;
        //return;
        //this->keep_children = true;
        return this;
    }
    //return;
    return NULL;
}

/**
 * Given a set of leaves, goes through this treeNode's subtree to look for a clade
 * whose subtree exactly matches the set of leaves. Returns NULL on failure.
 */
treeNode* treeNode::get_clade_match(cladeset& clade){
    for (vector<treeNode*>::iterator child = this->children.begin(); child != this->children.end();
        ++child){
        cladeset child_subtree = (*child)->subtree_leaves();
        if (issubset_bitset(clade, child_subtree)){
            return (*child)->get_clade_match(clade);
        }
    }
    
    // No child has the clade in its subtree. Does it exactly match this clade's
    // subtree?
    cladeset this_subtree = this->subtree_leaves();
    if (seteq_bitset(clade, this_subtree)){
        return this;
    }
    
    return NULL;
}

/**
 * Given a clade, looks through this tree's subtree to see if it has any clade
 * that groups some of the clade of interest's members with nonmembers 
 * (invalidates the clade). Returns true if so, or false if not.
 */
bool treeNode::has_clade_invalidates(cladeset& clade){
    cladeset this_subtree = this->subtree_leaves();
    if (!issubset_bitset(this_subtree, clade) && !issuperset_bitset(this_subtree, clade) &&
        set_int_bitset(this_subtree, clade).count() > 0){
        return true;
    }
    for (vector<treeNode*>::iterator child = this->children.begin(); child != this->children.end();
        ++child){
        if ((*child)->has_clade_invalidates(clade)){
            return true;
        }
    }
    return false;
}

/**
 * Gets the lowest-level clades in the tree that invalidate a propagating clade.
 */
bool treeNode::get_clade_invalidates(cladeset& clade,
    vector<cladeset >& invalidates){
    
    bool child_invalidates = false;
    vector<cladeset > invalidates_tmp;
    
    for (vector<treeNode*>::iterator child = this->children.begin();
        child != this->children.end(); ++child){
        if ((*child)->get_clade_invalidates(clade, invalidates_tmp)){
            child_invalidates = true;
            
            for (vector<cladeset >::iterator it = 
                invalidates_tmp.begin(); it != invalidates_tmp.end(); ++it){
                if (find(invalidates.begin(), invalidates.begin(), (*it)) == invalidates.end()){
                    invalidates.push_back((*it));
                }
            }
            
        }
    }
    
    if (!child_invalidates){
        // Check to see if this clade invalidates the clade.
        cladeset this_subtree = this->subtree_leaves();
        if (!issubset_bitset(this_subtree, clade) && !issuperset_bitset(this_subtree, clade) &&
            set_int_bitset(this_subtree, clade).count() > 0){
            invalidates.push_back(this_subtree);    
            return true;
        }
    }
    
    return child_invalidates;
}

bool treeNode::get_clades_invalidating(cladeset& clade,
    vector<cladeset >& invalidates){
    // Goal: get every highest-level invalidating clade in the tree
    cladeset this_subtree = this->subtree_leaves();
    if (set_int_bitset(this_subtree, clade).count() > 0 && !issubset_bitset(this_subtree, clade) &&
        !issubset_bitset(clade, this_subtree)){
        invalidates.push_back(this_subtree);
        return true;
    }
    // Check in children.
    bool child_invalidates = false;
    for (vector<treeNode*>::iterator child = this->children.begin(); child != this->children.end();
        ++child){
        if ((*child)->get_clade_invalidates(clade, invalidates)){
            child_invalidates = true;
        }
    }
    return child_invalidates;

}

/**
 * Given a data structure, which will map total distances from root to nodes
 * to actual treeNode objects, inserts this node and all of its children
 * recursively into the data structure.
 */
void treeNode::to_dist_map(float dist_from_root, std::map<float, vector<treeNode*> >& nodes){
    float this_dist_from_root = this->dist + dist_from_root;
   
        if (nodes.count(this_dist_from_root) == 0){
            // Need to create a new empty vector for this distance from the root.
            vector<treeNode*> nodelist;
            nodes.insert(make_pair(this_dist_from_root, nodelist));
            
        }

        nodes[this_dist_from_root].push_back(this);
        
        for (vector<treeNode*>::iterator child = this->children.begin();
            child != this->children.end(); ++child){
            (*child)->to_dist_map(this_dist_from_root, nodes);
        }

    return;
}

void treeNode::flatten(vector<cladeset >& flat){
    for (vector<treeNode*>::iterator child = this->children.begin();
        child != this->children.end(); ++child){
        (*child)->flatten(flat);
    }
    cladeset subtree = this->subtree_leaves();
    if (subtree.count() < this->num_haps && subtree.count() > 1){
        flat.push_back(subtree);
    }
}

/**
 * Given a "correct" tree (and assuming this is the "test" tree), compute
 * how many of this node's subtree's clades exist in the "correct" tree (and
 * how many don't). Return the total number of nodes in this tree's subtree
 * that are "correct," excluding singletons and roots of both trees.
 */
int treeNode::clades_correct(treeNode* correct){
    int correct_sum = 0;
    if (this->parent != NULL && (this->subtree_leaves()).count() > 1){
        // See if this clade exists in the other tree.
        cladeset subtree = this->subtree_leaves();
        if (correct->get_clade_match(subtree) != NULL){
        //if (!correct->has_clade_invalidates(subtree)){
            correct_sum += 1;
        }
        
        else{
            /*
            vector<string> dummy;
            printf("(%s);\n", this->newick(false, false, dummy).c_str());
            printf("\n");
            treeNode* thisroot = this->parent;
            while(thisroot->parent != NULL){
                thisroot = thisroot->parent;
            }
            printf("%s\n", thisroot->newick(false, false, dummy).c_str());
            printf("\n");
            printf("%s\n", correct->newick(false, false, dummy).c_str());
            //exit(1);
            printf("**********\n");
            */
        }
        
        
    }
    for (vector<treeNode*>::iterator child = this->children.begin();
        child != this->children.end(); ++child){
        correct_sum += (*child)->clades_correct(correct);
    }
    return correct_sum;
}

/**
 * Returns the total number of clades in this clade's subtree, excluding singletons
 * and the root.
 */
int treeNode::clades_tot(){
    int node_sum = 0;
    if (this->parent != NULL && (this->children.size() > 0 || this->leaves.count() > 1)){
        node_sum += 1;
    }
    for (vector<treeNode*>::iterator child = this->children.begin();
        child != this->children.end(); ++child){
        node_sum += (*child)->clades_tot();
    }
    return node_sum;
}

pair<float, float> treeNode::mean_children(){
    float sum = 0;
    float count = 0;
    sum += (float) this->children.size();
    count += 1;
    for (vector<treeNode*>::iterator child = this->children.begin();
        child != this->children.end(); ++child){
        pair<float, float> child_data = (*child)->mean_children();
        // Only count nodes with children.
        if (child_data.first > 0){
            sum += child_data.first;
            count += child_data.second;
        }
    }
    return make_pair(sum, count);
}

set<unsigned int> sample_haplotypes(treeNode& tree, int num_to_sample, int num_haplotypes){
    
    // Start from the root. Need to sample all children.
    vector<cladeset > clades_sampled;
    map<float, vector<treeNode*> > nodes_sorted;
    tree.to_dist_map(0, nodes_sorted);
    int num_sampled = 0;

    map<float, vector<treeNode*> >::iterator it = nodes_sorted.begin();

    while (num_sampled < num_to_sample){
        // Add everything at this level.
        for (vector<treeNode*>::iterator level_it = it->second.begin();
            level_it != it->second.end(); ++level_it){
            cladeset clade_subtree = (*level_it)->subtree_leaves();
            // Remove this node's subtree from the subtrees of every other node 
            // in the list.
            for (vector<cladeset >::iterator otherclade = 
                clades_sampled.begin(); otherclade != clades_sampled.end();){
                *otherclade = set_diff_bitset(*otherclade, clade_subtree);
                if (otherclade->count() == 0){
                    // Remove it.
                    clades_sampled.erase(otherclade);
                    num_sampled--;
                }
                else{
                    // Should have already updated it; just increment the 
                    // iterator.
                    ++otherclade;
                }
            }
            clades_sampled.push_back(clade_subtree);
            num_sampled++;   
        }
        ++it;
    }
    set<unsigned int> haps_sampled;
    for (vector<cladeset >::iterator clade_it = clades_sampled.begin();
        clade_it != clades_sampled.end(); ++clade_it){
        set<unsigned int> samp = subsample_bitset(*clade_it, num_haplotypes, 1);
        haps_sampled.insert((* samp.begin()));
    }
    return haps_sampled;
}

/**
 * Given a bitset representing a chosen set of subsampled haplotypes, goes
 * through this node's entire subtree and removes all haplotypes not in the mask
 * from every node's leaves except for the root (where they are added).
 */
void treeNode::mask_haps(cladeset& mask){
    for (vector<treeNode*>::iterator child = this->children.begin();
        child != this->children.end();){
        (*child)->mask_haps(mask);
        if ((*child)->leaves.count() == 0 && (*child)->children.size() == 0){
            // Remove child.
            this->children.erase(child);
            //delete (*child);
        }
        else{
            ++child;
        }
    }
    
    if (this->parent == NULL){
        // This is the root; add all nodes not in the mask.
        this->leaves = this->leaves | ~mask;
    }
    else{
        this->leaves &= mask;
        // If this was just leaves that were masked, remove the entire node.
        /*
        if (this->leaves.count() == 0 && this->children.size() == 0){
            printf("removing 1\n");
            this->parent->children.erase(std::remove(this->parent->children.begin(),
                this->parent->children.end(), this), this->parent->children.end());
            printf("removing 2\n");
            delete this;
            printf("removed\n");
            return;
        }
        */
    }
    
    return;
}

/**
 * Prints a binary representation of self to stdout.
 *
 * Uses integers to keep track of names of treeNodes. Returns the highest currently-known
 * integer in the tree, so we can keep incrementing.
 */
void treeNode::serialize(gzFile& out){
    // Write in binary:
    // Recombination info (4 bits)
    // Dist (float)
    // Leaves (bitset)
    // Number of children (int)
    // Call this on all children recursively.
    
    // Create one byte worth of recombination data
    // bit 1: recomb left
    // bit 2: recomb right
    // bit 3: leaving clade left
    // bit 4: leaving clade right
    // bit 5: (unused)
    // bit 6: (unused)
    // bit 7: (unused)
    // bit 8: (unused)
    
    /*
    unsigned char recombdata;
    if (this->recomb_left){
        recombdata |= 0x80;
    }
    if (this->recomb_right){
        recombdata |= 0x40;
    }
    if (this->leaving_clade_left){
        recombdata |= 0x20;
    }
    if (this->leaving_clade_right){
        recombdata |= 0x10;
    }
    
    char bytes[1];
    bytes[0] = recombdata;
    fwrite(&bytes[0], sizeof(char), 1, out);
    */
    
    serialize_str(out, this->name);
    
    serialize_longint(out, this->persistence);
    
    serialize_float(out, this->dist);
    
    serialize_float(out, this->dist2);
    
    serialize_float(out, this->dist_lv);
    
    // Serialize mutations (sites)
    serialize_longint(out, this->mutations.size());
    for (set<long int>::iterator site = this->mutations.begin(); site != this->mutations.end(); ++site){
        serialize_longint(out, *site);
    }
    
    /*
    // Serialize non-mutation sites
    serialize_longint(out, this->sites.size());
    for (set<long int>::iterator site = this->sites.begin(); site != this->sites.end(); ++site){
        serialize_longint(out, *site);
    }
    */
    
    serialize_bitset(out, this->leaves, this->num_haps);
    
    /*
    serialize_int(out, (int) this->partner_clades_left.size());
    for (vector<cladeset >::iterator pc = this->partner_clades_left.begin();
        pc != this->partner_clades_left.end(); ++pc){
        serialize_bitset(out, *pc, this->num_haps);
    }
    serialize_int(out, (int) this->partner_clades_right.size());
    for (vector<cladeset >::iterator pc = this->partner_clades_right.begin();
        pc != this->partner_clades_right.end(); ++pc){
        serialize_bitset(out, *pc, this->num_haps);
    }
    */
    
    serialize_int(out, (int) this->children.size());

    for (vector<treeNode*>::iterator it = this->children.begin(); it != this->children.end(); ++it){
        (*it)->serialize(out);
    }
    
}

float unserialize_treeNode(treeNode& node, const int num_haplotypes, instream_info& is){
    // Initialize bitsets to the right size
    node.set_haps(num_haplotypes);
    
    // Since we are loading this from a file, it's safe to say that this clade
    // was present at the time of loading
    node.present_when_loaded = true;
    
    // Next is name
    unserialize_str(is, node.name);
    //fprintf(stderr, "name: %s\n", node.name.c_str());
    
    // Next is persistence
    unserialize_longint(is, node.persistence);
    
    // Next is branch length 
    unserialize_float(is, node.dist);
    if (isnan(node.dist)){
        node.dist = 1.0;
    }
    
    unserialize_float(is, node.dist2);
    if (isnan(node.dist2)){
        node.dist2 = 1.0;
    }
    
    unserialize_float(is, node.dist_lv);
    if (isnan(node.dist_lv)){
        node.dist_lv = 1.0;
    }
    
    // Next is number of mutations
    long int num_muts;
    unserialize_longint(is, num_muts);
    // If this number is less than zero, assume there are no actual mutation
    // sites stored in this treeNode (and this is just for branch lengths)
    if (num_muts < 0){
        node.num_mutations = abs(num_muts);
        num_muts = abs(num_muts);
    }
    else{
        node.num_mutations = num_muts;
        for (long int i = 0; i < num_muts; ++i){
            long int site;
            unserialize_longint(is, site);
            node.mutations.insert(site);
            /*
            char buf[20];
            sprintf(buf, "%ld|", site);
            node.name += string(buf);
            */
        }
    }
    
    
    // Next is number of non-mutation sites.
    // comment this out for backward compatibility
    /*
    long int num_sites;
    unserialize_longint(is, num_sites);
    for (long int i = 0; i < num_sites; ++i){
        long int site;
        unserialize_longint(is, site);
        //node.sites.insert(site);
    }
    */
    
    /*
    if (node.name.length() > 0){
        node.name = "S" + node.name.substr(0, node.name.length()-1);
    }
    */
    // Next is leaves
    unserialize_bitset(is, node.leaves, num_haplotypes);
    
    // Next is number of children.
    int num_children;
    unserialize_int(is, num_children);
    
    
    float child_dist_sum = 0;
    for (int i = 0; i < num_children; i++){
        treeNode* child = new treeNode();
        child_dist_sum = child_dist_sum + unserialize_treeNode(*child, num_haplotypes, is);
        if (child->recomb_left_subtree){
            node.recomb_left_subtree = true;
        }
        if (child->recomb_right_subtree){
            node.recomb_right_subtree = true;
        }
        child->parent = &node;
        node.children.push_back(child);
    }
    if (num_children == 0){
        num_children = 1;
    }
    child_dist_sum /= (float)num_children;
    node.dist_below = child_dist_sum;
    child_dist_sum += node.dist;
    return child_dist_sum;
    
}

/**
 * Called by parse_newick().
 * Finds all Newick strings representing subtrees at the current depth of the
 * Newick string and returns them as a vector of strings. Each can then
 * be parsed independently by parse_newick().
 */
void parse_subtree(vector<string>& subtrees, string& newick){
    //vector<string> subtrees;
    if (newick.size() == 0){
        return;
        //return subtrees;
    }
    
    // Store indices of subtrees.
    vector<pair<int, int> > subtree_inds;
    
    int parenCount = 0;
    int last_subtree_start = 0;
    int last_subtree_end = 0;
    
    for (int chr_index=0; chr_index < newick.size(); ++chr_index){
        if (newick[chr_index] == '('){
            parenCount++;
        }
        else if (newick[chr_index] == ')'){
            parenCount--;
            
            if (parenCount == 0){
                // Two things can happen here: either the string ends,
                // or we find a comma to indicate that this is the
                // start of another subtree at the same level.
                int subtree_end = chr_index+1;
                
                for (int chr_index2 = chr_index+1; chr_index2 < newick.size(); ++chr_index2){
                    if (newick[chr_index2] == ','){
                        subtree_end = chr_index2;
                        break;
                    }
                    else if (chr_index2 == newick.size()-1){
                        // We've hit the end of the tree without finding
                        // a comma.
                        subtree_end = newick.size();
                    }
                }
                
                pair<int, int> inds = make_pair(last_subtree_start, subtree_end);
                subtree_inds.push_back(inds);
                
                // Start next subtree at this index so we can avoid hitting
                // the comma that separates subtrees
                last_subtree_start = subtree_end + 1;
            }
        }
    }
    // Note: string must end with closed parenthesis, so the final subtree
    // indices should already be added to the list.
    for (vector<pair<int, int> >::iterator subtree_ind = subtree_inds.begin();
        subtree_ind != subtree_inds.end(); ++subtree_ind){
        char subtree[subtree_ind->second - subtree_ind->first + 1];
        memcpy(subtree, newick.c_str() + subtree_ind->first, subtree_ind->second - subtree_ind->first);
        subtree[subtree_ind->second - subtree_ind->first] = '\0';
        string subtree_str = string(subtree);
        subtrees.push_back(subtree_str);   
    }
        
    //return subtrees;
}



string extract_regex_match(string subject, const std::sub_match<const char*>& fullmatch, 
    const std::sub_match<const char*>& groupmatch){
    
    long int startpos = groupmatch.first - fullmatch.first;
    long int cpylen = groupmatch.second - groupmatch.first;
    
    char extracted[cpylen+1];
    memcpy(extracted, subject.c_str() + startpos, cpylen);
    extracted[cpylen] = '\0';
    return (string) extracted;
    
}


/**
 * Called by parse_newick().
 * Given a string representing a leaf, parses the name and branch length of the leaf.
 * Determines whether or not the leaf deserves its own node and either creates a node
 * and adds it as a child of the given node, or adds the leaf directly to the node's
 * set of leaves (modifies parent by reference)
 */
bool parse_leaf(string& leaf, string& leaf_add, float& brlen){
    // Regular expression for parsing leaves (and associated branch lengths)
    static regex leafmatch(R"((.+):([0-9\.\-]+))", regex::extended);
    
    // Return true = add leaf as child
    // or false = add leaf as leaf
    cmatch match;
    if (regex_match(leaf.c_str(), match, leafmatch)){
        string leafname = extract_regex_match(leaf, match[0], match[1]);
        string leafdist = extract_regex_match(leaf, match[0], match[2]);
        
        float dist = abs(atof(leafdist.c_str()));
       
        if (dist > 0){
            if (leafname.length() > 0){
                // Create a new node.
                leaf_add = leafname;
                brlen = dist;
                return true;
            }
        }
        else{
            leaf_add = leafname;
            return false;
        }
    }
    else{
        // No branch length provided.
        leaf_add = leaf;
        return false;
    }
    return false;
}

/**
 * (Recursive method)
 * Given a string in Newick format, parses it as a tree. Note that
 * trees nested within this string are also trees in Newick format.
 *
 * leaf_branches = allow branch lengths for leaves (default behavior is to
 *    create a separate node for each leaf that has a nonzero branch
 *    length. If False, all leaves will be treated as branch length zero
 *    and will not be given their own nodes.
**/
void parse_newick(treeNode& root, string newick, const unsigned int num_haplotypes){
    
    root.set_haps(num_haplotypes);
    root.dist = 0.0;
    
    bool is_root = false;
    if (newick[newick.size()-1] == ';'){
        is_root = true;
        newick.erase(newick.end()-1);
    }

    // Precompile all regular expressions
    // NOTE: we use Boost library regexes here, since they're Perl-style.
    try{
        // Regular expression to match whole Newick-format tree nodes (with final 
        // semicolon removed)
        //static regex nodematch(R"(\(([a-zA-Z0-9:\.,\(\)\[\]\{\}\-_]+)\)([a-zA-Z0-9\-_]+)?(:[0-9\.]+)?)", regex::extended);
        static regex nodematch(R"(\(([a-zA-Z0-9:.,()\-_]+)\)([a-zA-Z0-9\-_]+)?(:[0-9\.e\-]+)?)", regex::extended);

        // Get all root-level information out of the tree.
        // Note: we have to look for brackets inside strings, on the left
        // or right side. These indicate clades being blocked off by
        // recombination.
    
        cmatch match;
        if (regex_match(newick.c_str(), match, nodematch)){
    
            // Innertree will be whatever occurs inside the outermost set of 
            // parenthesis.
            string innertree = extract_regex_match(newick.c_str(), match[0], match[1]);
        
            // Look for recombination marks.
        
            // Look for node name
            if (match[2].second != match[2].first){
                //root.name = extract_regex_match(newick.c_str(), match[0], match[2]);
            }
        
            // Look for branch length
            if (match[3].second != match[3].first){
                string dist_str = extract_regex_match(newick.c_str(), match[0], match[3]);
                if (dist_str[0] == ':'){
                    dist_str = dist_str.substr(1, dist_str.length()-1);
                }
                root.dist = abs(atof(dist_str.c_str()));
            }
        
            // Find all leaves attached to this node.
        
            // First, check for leaves on the left.
        
            const char *left_paren_ptr = strchr(innertree.c_str(), '(');
            if (left_paren_ptr && left_paren_ptr != innertree.c_str()){
                int left_paren_index = left_paren_ptr - innertree.c_str();
                char leafStr[left_paren_index+1];
                strncpy(leafStr, innertree.c_str(), left_paren_index);
                leafStr[left_paren_index] = '\0';
                if (leafStr[0] == ','){
                    int trunc_index = strlen(leafStr)-1;
                    memmove(leafStr, leafStr+1, strlen(leafStr)-1);
                    leafStr[trunc_index] = '\0';
                }
                if (leafStr[strlen(leafStr)-1] == ','){
                    leafStr[strlen(leafStr)-1] = '\0';
                }
            
                string leafStr_str(leafStr);
                istringstream tokenizer(leafStr_str);
                string leaf;
            
                vector<string> leaves_str;
                while(std::getline(tokenizer, leaf, ',')){
                    leaves_str.push_back(leaf);
                }
                
                for (vector<string>::iterator ls = leaves_str.begin(); ls != leaves_str.end(); ++ls){
                    string leafstr = "";
                    float brlen = 0;
                    bool add_leaf = parse_leaf(*ls, leafstr, brlen);
                    unsigned int leaf_index = atoi(leafstr.c_str());
                    if (add_leaf){
                        if (false){
                        //if (leaves_str.size() == 1){
                            root.leaves.set(root.num_haps-1-leaf_index);
                            root.dist = brlen;
                        }
                        else{
                            cladeset leaves;
                            leaves.set(root.num_haps-1-leaf_index);
                            root.addChild(leaves, brlen);
                        }
                    }
                    else{
                        root.leaves.set(root.num_haps-1-leaf_index);
                    }
                }
                
                // Pare tree down to part after leaves
                innertree.erase(0, left_paren_index);
            
            }
        
            // Now check for leaves on the right.
            const char *right_paren_ptr = strrchr(innertree.c_str(), ')');
            if (right_paren_ptr && (right_paren_ptr - innertree.c_str()) != strlen(innertree.c_str())){
                int right_paren_index = right_paren_ptr - innertree.c_str();
                char leafStr[strlen(innertree.c_str())-right_paren_index+1];
                memcpy(leafStr, innertree.c_str()+right_paren_index+1, strlen(innertree.c_str())-right_paren_index);
                leafStr[strlen(innertree.c_str())-right_paren_index] = '\0';
            
                // Figure out whether the first "leaf" is actually a branch length.
                if (leafStr[0] == ':'){
                    int brlen_end_index;
                    for (int i = 0; i < strlen(leafStr); ++i){
                        brlen_end_index = i;
                        if (leafStr[i] == ','){
                            break;
                        }
                    }
                    // Include branch length in inner tree.
                    innertree.erase(right_paren_index + brlen_end_index, innertree.length() - right_paren_index - brlen_end_index - 1);
                
                    if (brlen_end_index < strlen(leafStr)-1){
                        // Chop off the branch length portion from leafStr.
                        int trunc_index = strlen(leafStr) - brlen_end_index - 1;
                        memmove(leafStr, leafStr + brlen_end_index + 1, strlen(leafStr) - brlen_end_index-1);
                        leafStr[trunc_index] = '\0';
                    }
                    else{
                        // leafStr is just a branch length 
                        // do nothing -- no leaves to process.
                        leafStr[0] = '\0';
                    }
                }
                else{
                    innertree.erase(right_paren_index+1, innertree.length()-right_paren_index+1);
                }
                // Only proceed if the leaf string is not empty
                if (strlen(leafStr) > 0){
                    string leafStr_str(leafStr);
            
                    istringstream tokenizer(leafStr_str);
                    string leaf;
                
                    while(std::getline(tokenizer, leaf, ',')){
                        string leafstr = "";
                        float brlen = 0;
                        bool add_leaf = parse_leaf(leaf, leafstr, brlen);
                    
                        unsigned int leaf_index = atoi(leafstr.c_str());
                        if (add_leaf){
                            cladeset leaves;
                            leaves.set(root.num_haps-1-leaf_index);
                            root.addChild(leaves, brlen);
                        }
                        else{
                            root.leaves.set(root.num_haps-1-leaf_index);
                        }
                    }
                }
            
                // We now have a leafless tree. Look for all clades contained within
                // and parse each.
            
                // We now have a smaller Newick-format tree or set of trees.
                // parse_subtree will look to see if it's one or multiple trees
                // at this level and parse each.
                vector<string> subtrees;
                parse_subtree(subtrees, innertree);
                
                for (vector<string>::iterator subtree_it = subtrees.begin();
                    subtree_it != subtrees.end(); ++subtree_it){
                    treeNode* child = new treeNode;
                    child->parent = &root;
                    child->set_haps(root.num_haps);
                    parse_newick(*child, (*subtree_it), num_haplotypes);
                    // NOTE: don't touch this. You have to make a copy or else
                    // you'll store a reference to something that gets deleted later on.
                    //root.children.push_back(new treeNode(child));
                
                    root.children.push_back(child);
                }
            
            }
            else{
                // Base case: this tree was just leaves.
                istringstream tokenizer(innertree);
                string leaf;
                while(std::getline(tokenizer, leaf, ',')){
                    string leafstr = "";
                    float brlen = 0;
                    bool add_leaf = parse_leaf(leaf, leafstr, brlen);
                
                    unsigned int leaf_index = atoi(leafstr.c_str());
                    if (add_leaf){
                        cladeset leaves;
                        leaves.set(root.num_haps-1-leaf_index);
                        root.addChild(leaves, brlen);
                    }
                    else{
                        root.leaves.set(root.num_haps-1-leaf_index);
                    }
                }
            }
        }
        else{
            // Not valid Newick format.
            fprintf(stderr, "ERROR parsing tree:\n");
            fprintf(stderr, "%s\n", newick.c_str());
            exit(1);
        }
    } catch (regex_error& e){
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
            fprintf(stderr, "regex error_brack\n");
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
    

    //return root;    
}

