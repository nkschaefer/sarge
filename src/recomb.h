// header file for functions related to determining recombination events
// from four haplotype test failures

#ifndef SARGERECOMB_H
#define SARGERECOMB_H

#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <string.h>
#include <vector>
#include <set>
#include <unordered_set>
#include <map>
#include <zlib.h>
#include <deque>
#include "common.h"
#include "argnode.h"


// Data structure to represent a recombination event.
// Recombination events are built from comparing multiple four-haplotype
// test failures with one another.
struct recomb_event {
    // The rightmost left index of a 4-hap test failure
    long int left;
    // The leftmost right index of a 4-hap test failure
    long int right;
    
    // Is alpha (left) clade known?
    bool alpha_known;
    // Is beta (joined) clade known?
    bool beta_known;
    // Is leaving (moved) clade known?
    bool leaving_known;
    
    // Sets to store actual clades
    cladeset alpha;
    cladeset beta;
    cladeset leaving;
    
    bool up;
    bool down;
    
    ~recomb_event(){
        alpha.reset();
        beta.reset();
        leaving.reset();
    }
    
    recomb_event(){
        alpha_known = false;
        beta_known = false;
        leaving_known = false;
        left = -1;
        right = -1;
        
        up = false;
        down = false;
    }
    
    recomb_event(const recomb_event& r){
        this->left = r.left;
        this->right = r.right;
        this->alpha_known = r.alpha_known;
        this->beta_known = r.beta_known;
        this->leaving_known = r.leaving_known;
        if (this->alpha_known){
            this->alpha = r.alpha;
        }
        if (this->beta_known){
            this->beta = r.beta;
        }
        if (this->leaving_known){
            this->leaving = r.leaving;
        }
        this->up = r.up;
        this->down = r.down;
        
    }
    
    recomb_event& operator=(const recomb_event& r){
        // Check for self assignment
        if (&r == this){
            return *this;
        }
        else{
            this->left = r.left;
            this->right = r.right;
            this->alpha_known = r.alpha_known;
            this->beta_known = r.beta_known;
            this->leaving_known = r.leaving_known;
            if (this->alpha_known){
                this->alpha = r.alpha;
            }
            if (this->beta_known){
                this->beta = r.beta;
            }
            if (this->leaving_known){
                this->leaving = r.leaving;
            }
            
            this->up = r.up;
            this->down = r.down;
            return *this;
        }
    }
};

extern std::vector<recomb_event> recomb_catalog;

void print_recomb(recomb_event&, int);
void print_recomb_stdout(recomb_event&, int);

cladeset adjust_clade_recomb(arg_node* node,
    arg_node* limit,
    bool fromleft,
    std::vector<recomb_event>& recomb_catalog,
    long int& recomb_limit,
    cladeset_set& lv_grps_recomb,
    std::vector<int>& failures_used);
    
void gather_recomb_neighborhood(arg_node* node,
    std::set<arg_node*>& allnodes_l,
    std::set<arg_node*>& allnodes_r,
    std::set<arg_node*>& lv_grps_all,
    arg_clademap& lv_grps,
    int num_haplotypes);

void gather_recomb_neighborhood_solved(arg_node* node,
    std::set<arg_node*>& allnodes_l,
    std::set<arg_node*>& allnodes_r,
    std::set<arg_node*>& lv_grps_all);
    
bool resolve_recomb(arg_node* node,
    arg_clademap& lv_grps,
    arg_sitemap& sites_pos,
    arg_clademap& sites_clade,
    long int& prop_dist,
    arg_node* root,
    bool speculative,
    arg_node*& candidate_other,
    std::vector<recomb_event>&,
    std::vector<arg_node*>&,
    int,
    FILE*,
    long int,
    long int,
    std::set<arg_node*>&);
    
#endif
