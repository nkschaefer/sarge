#ifndef SARGE_ARGNODE_H
#define SARGE_ARGNODE_H

#include <zlib.h>
#include <math.h>
#include <time.h>
#include <string>
#include <algorithm>
#include <vector>
#include <deque>
#include <set>
#include <iterator>
#include <map>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <cstdlib>
#include <utility>
#include <sys/stat.h>
#include "common.h"
#include <fstream>
#include <bitset>
#include <unordered_map>
#include "sets.h"
#include "treeNode.h"
#include "serialize.h"

struct arg_node;
struct rc_edge;
struct vert_edge;

typedef std::multimap<long int, arg_node*> arg_sitemap;
typedef std::unordered_multimap<cladeset, arg_node*> arg_clademap;

/**
 * Edges between argnodes representing recombination events (lateral edges)
 */
struct rc_edge{
    arg_node* node;
    int type;
    arg_node* final_dest;
    
    ~rc_edge(){
      
    }
    rc_edge(int type, arg_node* node, arg_node* final_dest){
        this->type = type;
        this->node = node;
        this->final_dest = final_dest;
    }
    rc_edge(int type, arg_node* node){
        this->type = type;
        this->node = node;
        this->final_dest = NULL;
    }
    rc_edge(){
        this->type = -1;
        this->node = NULL;
        this->final_dest = NULL;
    }
    rc_edge(const rc_edge& e){
        this->type = e.type;
        this->node = e.node;
        this->final_dest = e.final_dest;
    }
};

/**
 * Edges between argnodes representing vertical (parent/child) relationships.
 */
struct vert_edge{
    long int start;
    long int end;
    arg_node* to;
    
    long int context_start;
    long int context_end;
    
    vert_edge(long int ind1, long int ind2, arg_node* dest){
        this->start = ind1;
        this->end = ind2;
        this->to = dest;
        this->context_start = -1;
        this->context_end = -1;
    };
    vert_edge(const vert_edge& e){
        this->start = e.start;
        this->end = e.end;
        this->to = e.to;
        this->context_start = e.context_start;
        this->context_end = e.context_end;
    }
};

bool operator<(const vert_edge&, const vert_edge&);

struct arg_node;

void del_solved_connections_between(arg_node* main_left, arg_node* main_right);

/**
 * Represents a node in the ancestral recombination graph.
 */
struct arg_node{

    std::vector<rc_edge> edges_from;
    std::vector<rc_edge> edges_to;
    
    std::vector<rc_edge> edges_from_solved;
    std::vector<rc_edge> edges_to_solved;
    
    std::vector<rc_edge> edges_from_unsolvable;
    std::vector<rc_edge> edges_to_unsolvable;
    
    bool newnode;
    
    // Track the range of this node
    long int start;
    long int end;
    
    // Track all individual sites supporting this node within its range
    std::set<long int> sites;
    std::set<long int> mutations;
    // Sites for which this node is a leaving group
    std::set<long int> sites_lvgrp;
    
    // Is this a candidate leaving clade in a recombination event, or a full node?
    bool is_leaving_grp;
    
    // Is this clade marked by an observed site, or inferred by a recombination event?
    bool site_support;
    
    cladeset clade;
    
    // Have we tried to solve recombination events at this site yet?
    short solve_attempts;
    
    // Only store lowest-level parents (that fail 4 hap test with each other)
    //  and highest-level children
    std::set<vert_edge> parents;
    std::set<vert_edge> children;
    
    long int recomb_left;
    long int recomb_right;
    
    long int closest_mut_left;
    long int closest_mut_right;
    
    // Mark if it's meant to be skipped (in the list of removable sites, or all ancestral,
    // etc.)
    bool skip;
    
    // Delete these later
    bool is_new_left;
    bool is_new_right;
    
    // Store all other nodes with pointers to this one
    std::set<arg_node*> nodes_referencing;
    
    bool fixed_left;
    bool fixed_right;
    
    arg_node(){
        this->parents.clear();
        this->children.clear();
        this->is_new_left = false;
        this->is_new_right = false;
        this->recomb_left = -1;
        this->recomb_right = -1;
        this->newnode = false;
        this->site_support = true;
        this->skip = false;
        this->is_leaving_grp = false;
        this->fixed_left = false;
        this->fixed_right = false;
        this->solve_attempts = 0;
        this->closest_mut_left = -1;
        this->closest_mut_right = -1;
    }
    
    long int closest_mut_l(){
        //return *this->sites.rbegin();
        
        if (this->mutations.size() > 0){
            return *this->mutations.rbegin();
        }
        else{
            return this->closest_mut_left;
        }

    };
    
    long int closest_mut_r(){
        //return *this->sites.begin();
        
        if (this->mutations.size() > 0){
            return *this->mutations.begin();
        }
        else{
            return this->closest_mut_right;
        }
    };
    
    ~arg_node(){
        
        // Note: edges to should already be empty.
        // Delete connections to this node from other nodes.
        for (std::set<vert_edge>::iterator child = this->children.begin(); child != this->children.end();){
            for (std::set<vert_edge>::iterator e2 = child->to->parents.begin(); e2 != child->to->parents.end();){
                if (e2->to == this){
                    e2->to->parents.erase(e2++);
                    break;
                }
                else{
                    ++e2;
                }
            }
            this->children.erase(child++);
        }

        for (std::set<vert_edge>::iterator parent = this->parents.begin(); parent != this->parents.end();){
            for (std::set<vert_edge>::iterator e2 = parent->to->children.begin(); e2 != parent->to->children.end();){
                if (e2->to == this){
                    e2->to->children.erase(e2++);
                    break;
                }
                else{
                    ++e2;
                }
            }
            this->parents.erase(parent++);
        }
        
        for (std::vector<rc_edge>::iterator edge_from = this->edges_from.begin(); edge_from != this->edges_from.end();){
            if (edge_from->final_dest != NULL){
                for (std::vector<rc_edge>::iterator et = edge_from->final_dest->edges_to.begin();
                    et != edge_from->final_dest->edges_to.end();){
                    if (et->final_dest == this){
                        for (std::vector<rc_edge>::iterator lvfrom = et->node->edges_from.begin();
                            lvfrom != et->node->edges_from.end();){
                            if (lvfrom->node == edge_from->final_dest && lvfrom->type == et->type){
                                et->node->edges_from.erase(lvfrom);
                                break;
                            }
                            else{
                                ++lvfrom;
                            }
                        }
                        edge_from->final_dest->edges_to.erase(et);
                    }
                    else{
                        ++et;
                    }
                }
            }
            
            for (std::vector<rc_edge>::iterator edge_to = edge_from->node->edges_to.begin();
                edge_to != edge_from->node->edges_to.end();){
                if (edge_to->node == this){
                    edge_from->node->edges_to.erase(edge_to);
                }
                else{
                    ++edge_to;
                }
            }
            this->edges_from.erase(edge_from);
        }
        
        std::set<std::pair<arg_node*, arg_node*> > recpairs;
        for (std::vector<rc_edge>::iterator edge_from = this->edges_from_solved.begin();
            edge_from != this->edges_from_solved.end(); ++edge_from){
            if (edge_from->final_dest != NULL){
                recpairs.insert(std::make_pair(this, edge_from->final_dest));
            }
            else{
                for (std::vector<rc_edge>::iterator to = edge_from->node->edges_to_solved.begin();
                    to != edge_from->node->edges_to_solved.end(); ++to){
                    if (to->final_dest != NULL && to->node == this){
                        recpairs.insert(std::make_pair(to->final_dest, edge_from->node));
                    }
                }
            }
        }
        for (std::vector<rc_edge>::iterator edge_to = this->edges_to_solved.begin();
            edge_to != this->edges_to_solved.end(); ++edge_to){
            if (edge_to->final_dest != NULL){
                recpairs.insert(std::make_pair(edge_to->final_dest, this));
            }
            else{
                for (std::vector<rc_edge>::iterator from = edge_to->node->edges_from_solved.begin();
                    from != edge_to->node->edges_from_solved.end(); ++from){
                    if (from->final_dest != NULL && from->node == this){
                        recpairs.insert(std::make_pair(edge_to->node, from->final_dest));
                    }
                }
            }
        }
        
        for (std::set<std::pair<arg_node*, arg_node*> >::iterator rp = recpairs.begin();
            rp != recpairs.end(); ++rp){
            del_solved_connections_between(rp->first, rp->second);
        }
        
        for (std::vector<rc_edge>::iterator edge_to = this->edges_to.begin(); edge_to != this->edges_to.end();){
            
            for (std::vector<rc_edge>::iterator ef = edge_to->final_dest->edges_from.begin();
                ef != edge_to->final_dest->edges_from.end();){
                if (ef->final_dest == this){
                    for (std::vector<rc_edge>::iterator lvto = ef->node->edges_to.begin();
                        lvto != ef->node->edges_to.end();){
                        if (lvto->node == edge_to->final_dest && lvto->type == ef->type){
                            ef->node->edges_to.erase(lvto);
                            break;
                        }
                        else{
                            ++lvto;
                        }
                    }
                    edge_to->final_dest->edges_from.erase(ef);
                }
                else{
                    ++ef;
                }
            }
            
            for (std::vector<rc_edge>::iterator edge_from = edge_to->node->edges_from.begin();
                edge_from != edge_to->node->edges_from.end();){
                if (edge_from->node == this){
                    edge_to->node->edges_from.erase(edge_from);
                }
                else{
                    ++edge_from;
                }
            }
            this->edges_to.erase(edge_to);
        }
    
        // Erase unsolvable edges.
        for (std::vector<rc_edge>::iterator edge_from = this->edges_from_unsolvable.begin();
            edge_from != this->edges_from_unsolvable.end();){
            for (std::vector<rc_edge>::iterator et = edge_from->final_dest->edges_to_unsolvable.begin();
                et != edge_from->final_dest->edges_to_unsolvable.end();){
                if (et->final_dest == this){
                    edge_from->final_dest->edges_to_unsolvable.erase(et);
                    //break;
                }
                else{
                    ++et;
                }
            }
            this->edges_from_unsolvable.erase(edge_from);
        }
        for (std::vector<rc_edge>::iterator edge_to = this->edges_to_unsolvable.begin();
            edge_to != this->edges_to_unsolvable.end();){
            for (std::vector<rc_edge>::iterator ef = edge_to->final_dest->edges_from_unsolvable.begin();
                ef != edge_to->final_dest->edges_from_unsolvable.end();){
                if (ef->final_dest == this){
                    edge_to->final_dest->edges_from_unsolvable.erase(ef);
                    //break;
                }  
                else{
                    ++ef;
                }
            }
            this->edges_to_unsolvable.erase(edge_to);
        }
        
        this->sites.clear();
        this->mutations.clear();
        this->clade.reset();
        
    }
};

bool operator<(const arg_node&, const arg_node&);
   
bool sort_nodes_size(const arg_node*, const arg_node*);

struct nodeptr_parent_sort{
    bool operator() (const arg_node* n1, const arg_node* n2) const {
        if (n1->clade.count() == n2->clade.count()){
            return n1 < n2;
        }
        else{
            return n1->clade.count() < n2->clade.count();
        }
    }
};

struct nodeptr_child_sort{
    bool operator() (const arg_node* n1, const arg_node* n2) const {
        if (n1->clade.count() == n2->clade.count()){
            return n1 < n2;
        }
        else{
            return n1->clade.count() > n2->clade.count();  
        }
    }
};

void transfer_edges_from(arg_node*, arg_node*, int);

void transfer_edges_to(arg_node*, arg_node*, int);

void add_solved_edges(arg_node*, arg_node*, arg_node*, int, bool, std::vector<arg_node*>&, long int, int);

void add_unsolved_edges(arg_node*, arg_node*, int);

bool ranges_overlap(long int, long int, long int, long int);

bool node_ranges_overlap(arg_node*, arg_node*);

bool site_ranges_overlap(long int, long int, arg_node*);

bool has_node_invalidates(arg_node* parent, long int pos, cladeset& clade);

bool has_node(arg_node* parent, long int pos, cladeset& clade);

arg_node* get_smallest_containing(arg_node* parent, long int pos, cladeset& clade);

void del_leaving_node(arg_node* lv, int num_haplotypes);

long int node_dist(arg_node* node1, arg_node* node2);

void compile_failures(arg_node* child, 
    arg_node* parent, 
    std::set<arg_node*>& failures,
    long int start, 
    long int end, 
    long int prop_dist,
    arg_sitemap& sites_pos,
    int num_haplotypes,
    bool exp_l,
    bool exp_r);

void gather_parents(arg_node* child, 
    arg_node* parent, 
    long int start, 
    long int end, 
    long int prop_dist,
    arg_sitemap& sites_pos,
    arg_clademap& sites_clade,
    arg_clademap& lv_grps,
    std::set<arg_node*, nodeptr_parent_sort>& nodes,
    int num_haplotypes);

arg_node* handle_failure(arg_node* node1, arg_node* node2, 
    arg_sitemap& sites_pos,
    arg_clademap& sites_clade,
    arg_clademap& lv_grps,
    long int prop_dist,
    bool finalize,
    arg_node* root);

void handle_all_failures(arg_node* child,
    arg_node* parent,
    arg_sitemap& sites_pos,
    arg_clademap& sites_clade,
    arg_clademap& lv_grps,
    long int prop_dist,
    long int start,
    long int end,
    arg_node* root,
    std::vector<arg_node*>& newnodes,
    int num_haplotypes);

void adjust_all_failure_indices(arg_node* child, 
    std::set<arg_node*>& failures,
    arg_sitemap& sites_pos,
    arg_clademap& sites_clade,
    arg_clademap& lv_grps,
    long int prop_dist,
    arg_node* root,
    std::vector<arg_node*>& newnodes,
    int num_haplotypes);

void handle_failures(arg_node* child,
    std::set<arg_node*>& failures,
    arg_sitemap& sites_pos,
    arg_clademap& sites_clade,
    arg_clademap& lv_grps,
    long int prop_dist,
    arg_node* root,
    std::vector<arg_node*>& newnodes,
    int num_haplotypes);

void insert_below_parent(arg_node* child, 
    arg_node* parent, 
    arg_sitemap& sites_pos,
    arg_clademap& sites_clade,
    arg_clademap& lv_grps,
    long int prop_dist,
    long int start, 
    long int end,
    arg_node* root,
    std::vector<arg_node*>& newnodes,
    int num_haplotypes);

void elim_impossible_lv_grps(arg_clademap& lv_grps,
    arg_node* node, int num_haplotypes);

bool safe_to_create_node(cladeset& clade,
    arg_sitemap& sites_pos,
    arg_clademap& sites_clade,
    long int pos,
    long int prop_dist,
    int num_haplotypes);

arg_node* create_node(arg_node* node, 
    arg_sitemap& sites_pos,
    arg_clademap& sites_clade,
    arg_clademap& lv_grps,
    long int prop_dist,
    arg_node* root,
    std::vector<arg_node*>& newnodes,
    int num_haplotypes,
    std::set<arg_node*>& nodes_external);

bool sp_aux(arg_node* child, 
    arg_node* parent, 
    std::vector<vert_edge>& edges_in,
    std::vector<vert_edge>& edges_new, 
    bool downward,
    bool make_changes,
    long int prop_dist,
    std::set<arg_node*>& nodes_blocking,
    arg_sitemap& sites_pos,
    arg_clademap& sites_clade,
    arg_clademap& lv_grps,
    std::set<arg_node*, nodeptr_child_sort>& discarded,
    int num_haplotypes);

void set_parents(arg_node* child, 
    arg_node* parent, 
    bool downward_first, 
    long int prop_dist,
    arg_sitemap& sites_pos,
    arg_clademap& sites_clade,
    arg_clademap& lv_grps,
    std::set<arg_node*, nodeptr_child_sort>& discarded,
    int num_haplotypes);

void expand_range_fixed(arg_node* node, long int site,
    arg_sitemap& sites_pos,
    arg_clademap& sites_clade,
    arg_clademap& lv_grps,
    long int prop_dist,
    arg_node* root,
    std::vector<arg_node*>& newnodes,
    bool fix_end,
    int num_haplotypes,
    bool no_failures,
    std::set<arg_node*>&);

void update_recinds_l(arg_node* l);

void update_recinds_r(arg_node* r);

void unfix_edges_right(arg_node* n);
void unfix_edges_left(arg_node* n);

void expand_range(arg_node* node, long int amount, 
    arg_sitemap& sites_pos,
    arg_clademap& sites_clade,
    arg_clademap& lv_grps,
    long int prop_dist,
    arg_node* root,
    std::vector<arg_node*>& newnodes,
    int num_haplotypes,
    bool newsites);

void find_new_parents(arg_node* node, 
    arg_sitemap& sites_pos,
    arg_clademap& sites_clade,
    arg_clademap& lv_grps,
    long int prop_dist,
    long int start,
    long int end,
    arg_node* root,
    std::vector<arg_node*>& newnodes,
    int num_haplotypes);

void find_new_children(arg_node* node,
    arg_sitemap& sites_pos,
    arg_clademap& sites_clade,
    arg_clademap& lv_grps,
    long int prop_dist,
    long int start, 
    long int end,
    arg_node* root,
    std::vector<arg_node*>& newnodes,
    int num_haplotypes);

void shorten_range(arg_node* node, long int amount, 
    arg_sitemap& sites_pos,
    arg_clademap& sites_clade,
    arg_clademap& lv_grps,
    long int prop_dist,
    int num_haplotypes,
    arg_node* root,
    std::vector<arg_node*>& newnodes,
    long int del_lim);

bool lv_grp_compatible(arg_node* left_node, 
    arg_node* right_node,
    int type,
    arg_node* lv_grp);

void merge_lv_nodes(arg_node* lv1, 
    arg_node* lv2, 
    std::map<cladeset, 
    std::set<arg_node*> >& lv_grps);

void merge_lv_nodes_recomb(std::set<arg_node*>& lv_grps_all,
    arg_clademap& lv_grps);
    
arg_node* add_lv_grp_graph(cladeset& clade,
    int type,
    arg_clademap& lv_grps,
    arg_node* left_node,
    arg_node* right_node,
    arg_sitemap& sites_pos,
    arg_clademap& sites_clade,
    long int prop_dist,
    int num_haplotypes);

void del_connections_from_node(arg_node*, arg_node*, bool, std::set<arg_node*>&);
void del_connections_to_node(arg_node*, arg_node*, bool, std::set<arg_node*>&);

void erase_lv_grps(std::set<arg_node*>& lv_grps_delete,
    arg_clademap& lv_grps, int num_haplotypes);

void del_connections_between(arg_node* main_left, arg_node* main_right,
    arg_clademap& lv_grps, std::set<arg_node*>& lv_grps_set, int num_haplotypes);
   

bool path_exists_between(arg_node* left, arg_node* right);

bool l_compat_r(arg_node* left, arg_node* right);

bool l_compat_r_conservative(arg_node* left, arg_node* right, bool recomb_left);

void set_4haptest_failure(arg_node* left_node, 
    arg_node* right_node,
    long int left_site,
    long int right_site,
    arg_sitemap& sites_pos,
    arg_clademap& sites_clade,
    arg_clademap& lv_grps,
    long int prop_dist,
    arg_node* root,
    std::vector<arg_node*>& newnodes,
    int num_haplotypes);

void join_nodes(arg_node* node1,
    arg_node* node2,
    arg_sitemap& sites_pos,
    arg_clademap& sites_clade,
    arg_clademap& lv_grps,
    long int prop_dist,
    bool exp_range,
    arg_node* root,
    std::vector<arg_node*>& newnodes,
    int num_haplotypes);

arg_node* split_node(arg_node* node, 
    long int split_pos, 
    arg_sitemap& sites_pos,
    arg_clademap& sites_clade,
    arg_clademap& lv_grps,
    long int prop_dist,
    arg_node* root,
    std::vector<arg_node*>& newnodes,
    int num_haplotypes);

arg_node* handle_failure(arg_node* node1, arg_node* node2, 
    arg_sitemap& sites_pos,
    arg_clademap& sites_clade,
    arg_clademap& lv_grps,
    long int prop_dist,
    bool finalize,
    arg_node* root,
    std::vector<arg_node*>& newnodes,
    int num_haplotypes);

void expand_range_newsites(arg_node* node, 
    long int amount, 
    arg_sitemap& sites_pos,
    arg_clademap& sites_clade,
    arg_clademap& lv_grps,
    long int prop_dist,
    bool newnode,
    arg_node* root,
    std::vector<arg_node*>& newnodes,
    int num_haplotypes);

void update_recinds_left(arg_node*);
void update_recinds_right(arg_node*);

/**
 * Printing related functions
 */
void print_node_recursive(arg_node* n, int indent_level, int num_haplotypes);
void print_node_recursive_upward(arg_node* n, int indent_level, int num_haplotypes);
void print_node_recursive_pos(arg_node* n, long int pos, int indent_level, int num_haplotypes);
void print_node_recursive_upward_pos(arg_node* n, long int pos, int indent_level, int num_haplotypes);
void print_node_upward_onelevel(arg_node* n, int indent_level, int num_haplotypes);
void print_node_downward_onelevel(arg_node* n, int indent_level, int num_haplotypes);
void print_node_lite(arg_node* n, int num_haplotypes);
void print_graph_sif(std::set<arg_node*>& lv_grps_all, int num_haplotypes);
void print_graph_sif_lim(std::set<arg_node*>& lv_grps_all, std::set<arg_node*>& lnodes, std::set<arg_node*>& rnodes, int num_haplotypes);
void print_recombs_right(arg_node* n, int num_haplotypes);
void print_recombs_right_solved(arg_node* n, int num_haplotypes);
void print_recombs_left(arg_node* n, int num_haplotypes);
void print_recombs_left_solved(arg_node* n, int num_haplotypes);


std::string node2name(arg_node* node, int num_haplotypes);

#endif

