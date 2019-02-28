/**
 * Contains functions related to compiling four-haplotype test failures and
 * combining them into recombination events.
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
#include <vector>
#include <set>
#include <unordered_set>
#include <map>
#include <deque>
#include <random>
#include <math.h>
#include <zlib.h>
#include "sets.h"
#include "common.h"
#include "argnode.h"
#include "recomb.h"
#include "rc_trees.h"

using namespace std;

void print_recomb(recomb_event& recomb, int num_haplotypes){
    fprintf(stderr, "RECOMB %ld %ld", recomb.left, recomb.right);
    if (recomb.up){
        fprintf(stderr, " UP");
    }
    else if (recomb.down){
        fprintf(stderr, " DOWN");
    }
    fprintf(stderr, "\n");

    fprintf(stderr, "alpha\n");
    if (recomb.alpha_known){
        print_bitset_set(recomb.alpha, num_haplotypes);
    }
    else{
        fprintf(stderr, "-----\n");
    }
    fprintf(stderr, "beta\n");
    if (recomb.beta_known){
        print_bitset_set(recomb.beta, num_haplotypes);
    }
    else{
        fprintf(stderr, "-----\n");
    }
    fprintf(stderr, "leaving\n");
    if (recomb.leaving_known){
        print_bitset_set(recomb.leaving, num_haplotypes);
    }
    else{
        fprintf(stderr, "-----\n");
    }
}

void print_recomb_stdout(recomb_event& recomb, int num_haplotypes){
    fprintf(stdout, "RECOMB %ld %ld", recomb.left, recomb.right);
    if (recomb.up){
        fprintf(stdout, " UP");
    }
    else if (recomb.down){
        fprintf(stdout, " DOWN");
    }
    fprintf(stdout, "\n");
    
    fprintf(stdout, "alpha\n");
    if (recomb.alpha_known){
        set<unsigned int> s = bitset2set(recomb.alpha, num_haplotypes);
        print_set_num_stdout(s);
    }
    else{
        fprintf(stdout, "-----\n");
    }
    fprintf(stdout, "beta\n");
    if (recomb.beta_known){
        set<unsigned int> s = bitset2set(recomb.beta, num_haplotypes);
        print_set_num_stdout(s);
    }
    else{
        fprintf(stdout, "-----\n");
    }
    fprintf(stdout, "leaving\n");
    if (recomb.leaving_known){
        set<unsigned int> s = bitset2set(recomb.leaving, num_haplotypes);
        print_set_num_stdout(s);
    }
    else{
        fprintf(stdout, "-----\n");
    }
    fprintf(stdout, "\n");
}

cladeset adjust_clade_recomb(arg_node* node,
    arg_node* limit,
    bool fromleft,
    vector<recomb_event>& recomb_catalog,
    long int& recomb_limit,
    cladeset_set& lv_grps_recomb,
    vector<int>& recombs_used){

    cladeset clade = node->clade;
    
    if (fromleft){    
            
        long int start = *node->sites.rbegin();
        long int end = *limit->sites.begin();
        
        recomb_limit = start;
        
        int rec_ind = 0;
        for (vector<recomb_event>::iterator r = recomb_catalog.begin(); r != recomb_catalog.end();
            ++r){
            
            bool has_rec = false;
            
            if (ranges_overlap(r->left, r->right, start, end) && r->right <= end){
                
                if (r->left > recomb_limit){
                    recomb_limit = r->left;
                }
                
                cladeset alphaprime = set_diff_bitset(r->alpha, r->leaving);
                cladeset betaprime = set_diff_bitset(r->beta, r->leaving);
                if (r->up){
                    if (issuperset_bitset(clade, r->alpha) && !issuperset_bitset(clade, betaprime)){
                        clade = set_diff_bitset(clade, r->leaving);
                        lv_grps_recomb.insert(r->leaving);
                        has_rec = true;
                    }
                }
                else if (r->down){
                    if (issuperset_bitset(clade, betaprime) && !issuperset_bitset(clade, r->alpha)){
                        if (clade != betaprime){
                            clade |= r->leaving;
                            lv_grps_recomb.insert(r->leaving);
                            has_rec = true;
                        }
                    }
                }
                else{
                    if (issuperset_bitset(clade, r->alpha) && !issuperset_bitset(clade, r->beta)){
                        clade = set_diff_bitset(clade, r->leaving);
                        lv_grps_recomb.insert(r->leaving);
                        has_rec = true;
                    }
                    else if (issuperset_bitset(clade, betaprime) && !issuperset_bitset(clade, r->leaving)){
                        if (clade != betaprime){
                            clade |= r->leaving;
                            lv_grps_recomb.insert(r->leaving);
                            has_rec = true;
                        }
                    }
                }
            }
            if (has_rec){
                recombs_used.push_back(rec_ind);
            }
            ++rec_ind;
        }
        
    }
    else{
        long int start = *limit->sites.rbegin();
        long int end = *node->sites.begin();
        
        recomb_limit = end;
        
        int rec_ind = recomb_catalog.size()-1;
        for (vector<recomb_event>::reverse_iterator r = recomb_catalog.rbegin(); r != recomb_catalog.rend(); ++r){
            
            bool has_rec = false;
            
            if (ranges_overlap(r->left, r->right, start, end) && r->left >= start){
                
                if (r->right < end){
                    recomb_limit = r->right;
                }
                
                cladeset alphaprime = set_diff_bitset(r->alpha, r->leaving);
                cladeset betaprime = set_diff_bitset(r->beta, r->leaving);
                if (r->up){
                    if (issuperset_bitset(clade, alphaprime) && !issuperset_bitset(clade, r->beta)){
                        if (clade != alphaprime){
                            clade |= r->leaving;
                            lv_grps_recomb.insert(r->leaving);
                            has_rec = true;
                        }
                    }
                }
                else if (r->down){
                    if (issuperset_bitset(clade, r->beta) && !issuperset_bitset(clade, alphaprime)){
                        clade = set_diff_bitset(clade, r->leaving);
                        lv_grps_recomb.insert(r->leaving);
                        has_rec = true;
                    }
                }
                else{
                    if (issuperset_bitset(clade, r->beta) && !issuperset_bitset(clade, r->alpha)){
                        clade = set_diff_bitset(clade, r->leaving);
                        lv_grps_recomb.insert(r->leaving);
                        has_rec = true;
                    }
                    else if (issuperset_bitset(clade, alphaprime) && !issuperset_bitset(clade, r->leaving)){
                        if (clade != alphaprime){
                            clade |= r->leaving;
                            lv_grps_recomb.insert(r->leaving);
                            has_rec = true;
                        }
                    }
                }
            } 
            if (has_rec){
                recombs_used.push_back(rec_ind);
            }
            --rec_ind;
        }
        
    }

    return clade;
}


bool recomb_check_prelim(arg_node* node, 
    arg_clademap& lv_grps, int num_haplotypes){

    set<arg_node*> dummy;
    
    // Why isn't this working when the recomb is solved?
    long int closest_rec_ind = -1;
    //arg_node* closest_rec_node = NULL;
    for (vector<rc_edge>::iterator from = node->edges_from_solved.begin(); from != node->edges_from_solved.end(); ++from){
        if (from->final_dest != NULL){

        }
        if (from->final_dest != NULL && !issuperset_bitset(node->clade, from->final_dest->clade)){
            if (closest_rec_ind == -1 || *from->final_dest->sites.begin() < closest_rec_ind){
                closest_rec_ind = *from->final_dest->sites.begin();
            }
        }
    }
    if (closest_rec_ind != -1){
        set<arg_node*> partners_rm;
        for (vector<rc_edge>::iterator from = node->edges_from.begin(); from != node->edges_from.end(); ++from){
            if (from->final_dest != NULL){
            }
            if (from->final_dest != NULL && *from->final_dest->sites.begin() >= closest_rec_ind){
                partners_rm.insert(from->final_dest);
            }
        }
        
    }
    
    if (node->recomb_right != -1){
        bool has_close_recomb = false;
        set<arg_node*> partners;
        for (vector<rc_edge>::iterator from = node->edges_from.begin(); from != node->edges_from.end(); ++from){
            if (from->final_dest != NULL){
                partners.insert(from->final_dest);
            }
            if (from->final_dest != NULL && *from->final_dest->sites.begin() <= node->recomb_right){
                has_close_recomb = true;
                break;
            }
        }
        if (!has_close_recomb){
            for (set<arg_node*>::iterator p = partners.begin(); p != partners.end(); ++p){
                
            }
            return false;
        }
    }
    return true;
}

/**
 * Given a left node involved in recombination, gathers all possible L, R, and leaving
 * nodes that could be involved in recombination with it.
 */
void gather_recomb_neighborhood(arg_node* node,
    set<arg_node*>& allnodes_l,
    set<arg_node*>& allnodes_r,
    set<arg_node*>& lv_grps_all,
    arg_clademap& lv_grps,
    int num_haplotypes){

    // Store the furthest L/R positions that could be involved in this recombination.
    long int l_boundary = -1;
    long int r_boundary = -1;
    // Gather all possible L nodes (must be compatible with given node)
    for (vector<rc_edge>::iterator from = node->edges_from.begin(); from != node->edges_from.end(); ++from){
        if (from->final_dest != NULL){
            for (vector<rc_edge>::iterator to = from->node->edges_to.begin();
                to != from->node->edges_to.end(); ++to){
                if (*to->node->sites.rbegin() > *node->sites.rbegin()){
                    if (l_compat_r_conservative(node, to->node, true)){
                        allnodes_l.insert(to->node);
                        if (l_boundary == -1 || *to->node->sites.begin() < l_boundary){
                            l_boundary = *to->node->sites.begin();
                        }
                    }
                }
                else{
                    if (l_compat_r_conservative(to->node, node, true) &&
                        (to->node->recomb_right == -1 || to->node->recomb_right >= *node->sites.rbegin())){
                        allnodes_l.insert(to->node);
                        
                        if (l_boundary == -1 || *to->node->sites.begin() < l_boundary){
                            l_boundary = *to->node->sites.begin();
                        }
                    }
                }
            }
        }
    }
    // Gather all possible lv grps
    for (set<arg_node*>::iterator l = allnodes_l.begin(); l != allnodes_l.end(); ++l){
        for (vector<rc_edge>::iterator from = (*l)->edges_from.begin();
            from != (*l)->edges_from.end(); ++from){
            if (from->final_dest != NULL && *from->final_dest->sites.begin() > *node->sites.rbegin()){
                lv_grps_all.insert(from->node);  
                allnodes_r.insert(from->final_dest);
                if (r_boundary == -1 || *from->final_dest->sites.rbegin() > r_boundary){
                    r_boundary = *from->final_dest->sites.rbegin();
                }
            }
        }
    }
    for (set<arg_node*>::iterator l = allnodes_l.begin(); l != allnodes_l.end();){
        bool connected = false;
        for (vector<rc_edge>::iterator from = (*l)->edges_from.begin(); from != (*l)->edges_from.end(); ++from){
            if (from->final_dest != NULL && allnodes_r.find(from->final_dest) != allnodes_r.end()){
                connected = true;
                break;
            }
        }
        if (connected){
            ++l;
        }
        else{
            allnodes_l.erase(l++);
        }
    }

    if (DEBUG_MODE){
        fprintf(stderr, "OUTER BOUNDARIES %ld %ld\n", l_boundary, r_boundary);
    }
    
    // Because of how ranges of leaving nodes are determined, it's possible that
    // a leaving node that should be in this set didn't make it into the set because
    // its range is adjacent to, but not overlapping, the range of another identical
    // leaving node in the set. In other words, one or more leaving nodes in the
    // set may need to be merged with other leaving nodes that have the same clade.
    // Check for this.
    set<arg_node*> lv_grp_add;
    for (set<arg_node*>::iterator lv = lv_grps_all.begin(); lv != lv_grps_all.end();
        ++lv){
        pair<arg_clademap::iterator, arg_clademap::iterator> er = lv_grps.equal_range((*lv)->clade);
        for (arg_clademap::iterator otherlv = er.first; otherlv != er.second; ++otherlv){
            if (otherlv->second != *lv){

                bool overlap_set = false;
                if (site_ranges_overlap(l_boundary, r_boundary, otherlv->second)){
                    overlap_set = true;
                }
                if (overlap_set){
                    if (DEBUG_MODE){
                        fprintf(stderr, "OTHERLV OVERLAP:\n");
                        print_node_lite(otherlv->second, num_haplotypes);
                    }
                    set<arg_node*> lnodes_add;
                    
                    for (vector<rc_edge>::iterator to = otherlv->second->edges_to.begin();
                        to != otherlv->second->edges_to.end(); ++to){
                        if (*to->node->sites.rbegin() > *node->sites.rbegin()){
                            if (l_compat_r_conservative(node, to->node, true)){
                                lnodes_add.insert(to->node);
                                if (DEBUG_MODE){
                                    fprintf(stderr, "COMPAT\n");
                                    print_node_lite(to->node, num_haplotypes);
                                }
                            }
                            else if (DEBUG_MODE){
                                fprintf(stderr, "INCOMPAT\n");
                                print_node_lite(to->node, num_haplotypes);
                            }
                        }
                        else{
                            if (l_compat_r_conservative(to->node, node, true)){
                                lnodes_add.insert(to->node);
                                if (DEBUG_MODE){
                                    fprintf(stderr, "COMPAT\n");
                                    print_node_lite(to->node, num_haplotypes);
                                }
                            }
                            else if (DEBUG_MODE){
                                fprintf(stderr, "INCOMPAT\n");
                                print_node_lite(to->node, num_haplotypes);
                            }
                        }
                    }
                    if (lnodes_add.size() > 0){
                        lv_grp_add.insert(otherlv->second);
                        for (set<arg_node*>::iterator lnode = lnodes_add.begin();
                            lnode != lnodes_add.end(); ++lnode){
                            allnodes_l.insert(*lnode);
                            for (vector<rc_edge>::iterator from = (*lnode)->edges_from.begin();
                                from != (*lnode)->edges_from.end(); ++from){
                                if (from->node == otherlv->second && from->final_dest != NULL){
                                    allnodes_r.insert(from->final_dest);
                                }
                            }
  
                        }
                    }
                }
                else if (DEBUG_MODE){
                    fprintf(stderr, "no overlap\n");
                    print_node_lite(otherlv->second, num_haplotypes);
                }
            }
        }
    }

    for (set<arg_node*>::iterator lvadd = lv_grp_add.begin(); lvadd != lv_grp_add.end(); ++lvadd){
        if (DEBUG_MODE){
            fprintf(stderr, "ADDING LV GRP\n");
            print_node_lite(*lvadd, num_haplotypes);
        }
        lv_grps_all.insert(*lvadd);
    }
    
    if (DEBUG_MODE){
        fprintf(stderr, "neighborhood:\n");
        print_graph_sif_lim(lv_grps_all, allnodes_l, allnodes_r, num_haplotypes);
        fprintf(stderr, "\n");
    }

}

void gather_recomb_neighborhood_solved(arg_node* node,
    set<arg_node*>& allnodes_l,
    set<arg_node*>& allnodes_r,
    set<arg_node*>& lv_grps_all){
    
    // Gather all possible L nodes (must be compatible with given node)
    for (vector<rc_edge>::iterator from = node->edges_from_solved.begin(); from != node->edges_from_solved.end(); ++from){
        for (vector<rc_edge>::iterator to = from->node->edges_to_solved.begin();
            to != from->node->edges_to_solved.end(); ++to){
            if (*to->node->sites.rbegin() > *node->sites.rbegin()){
                if (l_compat_r_conservative(node, to->node, true)){
                    allnodes_l.insert(to->node);
                }
            }
            else{
                if (l_compat_r_conservative(to->node, node, true)){
                    allnodes_l.insert(to->node);
                }
            }
        }
    }
    // Gather all possible lv grps
    for (set<arg_node*>::iterator l = allnodes_l.begin(); l != allnodes_l.end(); ++l){
        for (vector<rc_edge>::iterator from = (*l)->edges_from_solved.begin();
            from != (*l)->edges_from_solved.end(); ++from){
            lv_grps_all.insert(from->node);  
            if (from->final_dest != NULL){
                allnodes_r.insert(from->final_dest);
            } 
        }
    }
}

bool remove_disconnected(set<arg_node*>& allnodes_l, 
    set<arg_node*>& allnodes_r, 
    set<arg_node*>& lv_grps_all,
    bool check_l, 
    bool check_r,
    vector<pair<long int, arg_node*> >& lnodes_sorted,
    vector<pair<long int, arg_node*> >& rnodes_sorted,
    arg_node*& main_left,
    arg_node*& main_right,
    arg_node* node){

    // Compile lv grps from main node
    set<arg_node*> main_lv;
    for (vector<rc_edge>::iterator from = node->edges_from.begin(); from != node->edges_from.end(); ++from){
        if (from->final_dest != NULL && allnodes_r.find(from->final_dest) != allnodes_r.end()){
            main_lv.insert(from->node);
        }
        
    }
    
    if (check_l){
        for (set<arg_node*>::iterator l = allnodes_l.begin(); l != allnodes_l.end();){
            bool connected = false;
            for (vector<rc_edge>::iterator from = (*l)->edges_from.begin();
                from != (*l)->edges_from.end(); ++from){
                if (main_lv.find(from->node) != main_lv.end() && 
                    from->final_dest != NULL && allnodes_r.find(from->final_dest) != allnodes_r.end()){
                    connected = true;
                    break;
                }
            }
            if (connected){
                ++l;
            }
            else{
                for (vector<pair<long int, arg_node*> >::iterator ls = lnodes_sorted.begin();
                    ls != lnodes_sorted.end();){
                    if (ls->second == *l){
                        lnodes_sorted.erase(ls);
                        break;
                    }
                    else{
                        ++ls;
                    }   

                }
                allnodes_l.erase(l++);
            }
        }
        if (allnodes_l.size() == 0){
            return false;
        }
    }
    if (check_r){
        for (set<arg_node*>::iterator r = allnodes_r.begin(); r != allnodes_r.end();){
            bool connected = false;
            for (vector<rc_edge>::iterator to = (*r)->edges_to.begin(); to != (*r)->edges_to.end();
                ++to){
                if (main_lv.find(to->node) != main_lv.end() && 
                    to->final_dest != NULL && allnodes_l.find(to->final_dest) != allnodes_l.end()){
                    connected = true;
                    break;
                }
            }
            if (connected){
                ++r;
            }
            else{
                for (vector<pair<long int, arg_node*> >::iterator rs = rnodes_sorted.begin();
                    rs != rnodes_sorted.end();){
                    if (rs->second == *r){
                        rnodes_sorted.erase(rs);
                        break;
                    }
                    else{
                        ++rs;
                    }
                }
                allnodes_r.erase(r++);
            }
        }
        if (allnodes_r.size() == 0){
            return false;
        }
    }
    // Check leaving groups.
    for (set<arg_node*>::iterator lv = lv_grps_all.begin(); lv != lv_grps_all.end();){
        bool connected_l = false;
        bool connected_r = false;
        for (vector<rc_edge>::iterator to = (*lv)->edges_to.begin();
            to != (*lv)->edges_to.end(); ++to){
            if (allnodes_l.find(to->node) != allnodes_l.end()){
                connected_l = true;
                break;
            }
        }
        if (connected_l){
            for (vector<rc_edge>::iterator from = (*lv)->edges_from.begin();
                from != (*lv)->edges_from.end(); ++from){
                if (allnodes_r.find(from->node) != allnodes_r.end()){
                    connected_r = true;
                    break;
                }
            }
        }
        if (!connected_l || !connected_r){
            lv_grps_all.erase(lv++);
        }
        else{
            ++lv;
        }
    }
    if (lv_grps_all.size() == 0){
        return false;
    }
    
    main_left = lnodes_sorted.rbegin()->second;
    main_right = rnodes_sorted.begin()->second;
    
    return true;
}

bool filter_far_edges(set<arg_node*>& allnodes_l, 
    set<arg_node*>& allnodes_r,
    set<arg_node*>& lv_grps_all,
    vector<pair<long int, arg_node*> >& lnodes_sorted, 
    vector<pair<long int, arg_node*> >& rnodes_sorted,
    arg_node*& main_left,
    arg_node*& main_right,
    arg_node* node, 
    bool force_solve,
    arg_sitemap& sites_pos,
    arg_clademap& lv_grps,
    int num_haplotypes){
    //fprintf(stderr, "Rffe\n");
    if (allnodes_l.size() == 0 || allnodes_r.size() == 0 || lnodes_sorted.size() == 0 || rnodes_sorted.size() == 0){
        return false;
    }
    // Were any nodes deleted from the L or R sets? If so, we need to check later
    // and delete anything that no longer has a recombination partner in the
    // opposite set.
    bool l_removed = false;
    bool r_removed = false;
    
  
    // Remove anything that's blocked by recombination before getting to the other
    // side of the recombination event.
    while(lnodes_sorted.size() > 0 && lnodes_sorted[0].second->recomb_right != -1 && 
        //lnodes_sorted[0].second->recomb_right < rnodes_sorted[0].second->start){
        lnodes_sorted[0].second->recomb_right < *rnodes_sorted[0].second->sites.begin()){
        allnodes_l.erase(allnodes_l.find(lnodes_sorted[0].second));
        l_removed = true;
        
        if (DEBUG_MODE){
            fprintf(stderr, "removing Lmost\n");
            print_node_lite(lnodes_sorted.begin()->second, num_haplotypes);
        }
        lnodes_sorted.erase(lnodes_sorted.begin());
    }
    if (lnodes_sorted.size() == 0){
        return false;
    }
    
    
    while(rnodes_sorted.size() > 0 &&
        rnodes_sorted.rbegin()->second->recomb_left != -1 && 
        //rnodes_sorted.rbegin()->second->recomb_left > lnodes_sorted.rbegin()->second->end){
        rnodes_sorted.rbegin()->second->recomb_left > *lnodes_sorted.rbegin()->second->sites.rbegin()){
        allnodes_r.erase(allnodes_r.find(rnodes_sorted.rbegin()->second));
        r_removed = true;
        vector<pair<long int, arg_node*> >::reverse_iterator lastitem = rnodes_sorted.rbegin();
        vector<pair<long int, arg_node*> >::iterator lastitem_f = (lastitem + 1).base();
        
        if (DEBUG_MODE){
            fprintf(stderr, "removing Rmost\n");
            print_node_lite(rnodes_sorted[rnodes_sorted.size()-1].second, num_haplotypes);
        }
        rnodes_sorted.erase(lastitem_f);
    }
    if (rnodes_sorted.size() == 0){
        return false;
    }
    
    main_left = lnodes_sorted.rbegin()->second;
    main_right = rnodes_sorted.begin()->second;
    
    bool success = true;
    
    if (l_removed || r_removed){
        success = remove_disconnected(allnodes_l, allnodes_r, lv_grps_all, 
            r_removed, l_removed, 
            lnodes_sorted, rnodes_sorted, 
            main_left, main_right,
            node);
        if (DEBUG_MODE){
            fprintf(stderr, "after remove disconnected:\n");
            fprintf(stderr, "L\n");
            for (set<arg_node*>::iterator l = allnodes_l.begin(); l != allnodes_l.end();
                ++l){
                print_node_lite(*l, num_haplotypes);
            }
            fprintf(stderr, "R\n");
            for (set<arg_node*>::iterator r = allnodes_r.begin(); r != allnodes_r.end();
                ++r){
                print_node_lite(*r, num_haplotypes);
            }
        }
    }
    
    return success;
}

void clade_maxdists(map<arg_node*, float>& dists,
    set<arg_node*>& allnodes_l,
    set<arg_node*>& allnodes_r,
    vector<pair<long int, arg_node*> >& lnodes_sorted,
    vector<pair<long int, arg_node*> >& rnodes_sorted,
    vector<pair<int, arg_node*> >& lv_grps_degree,
    arg_sitemap& sites_pos,
    long int prop_dist,
    int num_haplotypes){

    long int left = *lnodes_sorted.rbegin()->second->sites.rbegin();
    long int right = *rnodes_sorted.begin()->second->sites.begin();
    
    int maxdeg = -lv_grps_degree.begin()->first;
    for (vector<pair<int, arg_node*> >::iterator lvd = lv_grps_degree.begin();
        lvd != lv_grps_degree.end(); ++lvd){
        if (-lvd->first < maxdeg){
            break;
        }
        else{
            float dist_sum = 0;
            float dist_tot = 0;
            unordered_set<cladeset> new_l;
            unordered_set<cladeset> new_r;
            bool alpha_l = false;
            bool beta_r = false;
            for (vector<rc_edge>::iterator to = lvd->second->edges_to.begin();
                to != lvd->second->edges_to.end(); ++to){
                if (allnodes_l.find(to->node) != allnodes_l.end()){
                    if (to->type == 1 || to->type == 2){
                        // alpha
                        cladeset newclade = set_diff_bitset(to->node->clade, lvd->second->clade);
                        new_r.insert(newclade);
                        alpha_l = true;
                    }
                    else if (to->type == 3){
                        // beta
                        cladeset newclade = to->node->clade | lvd->second->clade;
                        new_r.insert(newclade);
                    }
                }
            }
            for (vector<rc_edge>::iterator from = lvd->second->edges_from.begin();
                from != lvd->second->edges_from.end(); ++from){
                if (allnodes_r.find(from->node) != allnodes_r.end()){
                    if (from->type == 2 || from->type == 3){
                        // beta
                        cladeset newclade = set_diff_bitset(from->node->clade, lvd->second->clade);
                        new_l.insert(newclade);
                        beta_r = true;
                    }
                    else if (from->type == 1){
                        // alpha
                        cladeset newclade = from->node->clade | lvd->second->clade;
                        new_l.insert(newclade);
                    }
                }
            }
            new_l.insert(lvd->second->clade);
            new_r.insert(lvd->second->clade);
            if (!alpha_l || !beta_r){
                // Get union of all failures
                cladeset un;
                for (set<arg_node*>::iterator l = allnodes_l.begin(); l != allnodes_l.end(); ++l){
                    un |= (*l)->clade;
                }
                for (set<arg_node*>::iterator r = allnodes_r.begin(); r != allnodes_r.end(); ++r){
                    un |= (*r)->clade;
                }
                if (un.count() < num_haplotypes){
                    new_l.insert(un);
                    new_r.insert(un);
                }
            }
            for (unordered_set<cladeset>::iterator nl = new_l.begin(); nl != new_l.end(); ++nl){
                // See how far it can go.
                long int fail_ind = left - prop_dist;
                for (arg_sitemap::reverse_iterator it = arg_sitemap::reverse_iterator(sites_pos.equal_range(left).first);
                    it != sites_pos.rend(); ++it){
                    if (it->first < left - prop_dist){
                        break;
                    }
                    else{
                        short comp = compare_bitsets(it->second->clade, (*nl));
                        if (comp == -1){
                            fail_ind = it->first;
                            break;
                        }
                    }
                }
                dist_tot += 1;
                dist_sum += (left - fail_ind);
            }
            for (unordered_set<cladeset>::iterator nr = new_r.begin(); nr != new_r.end(); ++nr){
                // See how far it can go.
                long int fail_ind = right + prop_dist;
                for (arg_sitemap::iterator it = sites_pos.equal_range(right).second;
                    it != sites_pos.end(); ++it){
                    if (it->first > right + prop_dist){
                        break;
                    }
                    else{
                        short comp = compare_bitsets(it->second->clade, (*nr));
                        if (comp == -1){
                            fail_ind = it->first;
                            break;
                        }
                    }
                }
                dist_tot += 1;
                dist_sum += (fail_ind - right);
            }
            if (dist_tot == 0){
                dist_tot = 1;
            }
            dists.insert(make_pair(lvd->second, dist_sum/dist_tot));
        }
    }
    
}

/**
 * Given a full set of all possible L, R, and leaving nodes, finds all possible
 * groups of these nodes that could coexist with each other in a recombination 
 * event. If there is too much ambiguity, returns false to signify the recombination
 * event might be solved incorrectly. Otherwise, properly adjusts the sets of 
 * nodes to reflect the choice made.
 */
bool choose_nodeset(set<arg_node*>& allnodes_l,
    set<arg_node*>& allnodes_r,
    set<arg_node*>& lv_grps_all,
    vector<pair<long int, arg_node*> >& lnodes_sorted,
    vector<pair<long int, arg_node*> >& rnodes_sorted,
    arg_node* node, 
    arg_node*& main_left,
    arg_node*& main_right,
    bool force_solve,
    bool speculative,
    int num_haplotypes){

    // Remove R nodes that go past the "pivot" node, since they can't be involved in 
    // recombination with it
    for (set<arg_node*>::iterator rn = allnodes_r.begin(); rn != allnodes_r.end();){
        if (*(*rn)->sites.begin() <= *node->sites.rbegin()){
            for (vector<pair<long int, arg_node*> >::iterator rs = rnodes_sorted.begin();
                rs != rnodes_sorted.end();){
                if (rs->second == *rn){
                    rnodes_sorted.erase(rs);
                    break;
                }
                else{
                    ++rs;
                }
            }
            allnodes_r.erase(rn++);
        }
        else{
            ++rn;
        }
    }
    // Remove L nodes that now no longer have a recombination partner in the set of 
    // R nodes.
    for (set<arg_node*>::iterator ln = allnodes_l.begin(); ln != allnodes_l.end();){
        bool connected = false;
        for (vector<rc_edge>::iterator from = (*ln)->edges_from.begin(); from != (*ln)->edges_from.end(); ++from){
            if (from->final_dest != NULL && allnodes_r.find(from->final_dest) != allnodes_r.end()){
                connected = true;
                break;
            }
        }
        if (connected){
            ++ln;
        }
        else{
            for (vector<pair<long int, arg_node*> >::iterator ls = lnodes_sorted.begin();
                ls != lnodes_sorted.end();){
                if (ls->second == *ln){
                    lnodes_sorted.erase(ls);
                    break;
                }
                else{
                    ++ls;
                }
            }
            allnodes_l.erase(ln++);
        }
    }
    
    // If we have more nodes to choose from here than allowable haplotypes, we have 
    // a problem. If compiled with reasonably high MAXHAPS value, we should be okay.
    if (allnodes_l.size() + allnodes_r.size() > NODESETSIZE){
        fprintf(stderr, "ERROR: when solving a recombination event, there were more \
nodes to choose from than available space in bitsets. Please re-compile with a \
higher value of MAXHAPS (>= %.0f) as follows:\n", ceil((float)(allnodes_l.size() + allnodes_r.size())/5.0));
        fprintf(stderr, "make clean\n");
        fprintf(stderr, "make MAXHAPS=[new value]\n");
        exit(1);
    }
    
    // The opposite of choices -- if one node is kept, the other must be too.
    vector<nodeset > recpairs;
    
    // For every l or r node (besides node itself), determine what nodes it's incompatible
    // with. For each node, there is a choice between excluding it or excluding all
    // nodes it's incompatible with.
    int node_index = 0;
    vector<pair<nodeset, nodeset > > choices;
    for (vector<pair<long int, arg_node*> >::iterator l = lnodes_sorted.begin();
        l != lnodes_sorted.end(); ++l){
        // Choice to remove current node
        nodeset choice1;
        choice1.set(node_index);
        // Choice to remove all incompatible nodes
        nodeset choice2;
        // Look left
        bool incompat_found = false;
        for (int l_ind = node_index -1; l_ind >= 0; --l_ind){
            if (incompat_found){
                // Once an incompatibility is found, we must remove all nodes to the left of it
                choice2.set(l_ind);
            }
            else{
                if (!l_compat_r_conservative(lnodes_sorted[l_ind].second, l->second, true)){
                    choice2.set(l_ind);
                    incompat_found = true;
                }
            }
        }
        // Look right
        incompat_found = false;
        for (unsigned int l_ind = node_index + 1; l_ind < lnodes_sorted.size(); ++l_ind){
            if (incompat_found){
                // Once an incompatibility is found, we must remove all nodes to the right of it
                choice2.set(l_ind);
            }
            else{
                if (!l_compat_r_conservative(l->second, lnodes_sorted[l_ind].second, true)){
                    choice2.set(l_ind);
                    incompat_found = true;
                }
            }
        }
        
        if (DEBUG_MODE){
            fprintf(stderr, "--prelim L choice:\n");
            print_bitset_set2(choice1, lnodes_sorted.size() + rnodes_sorted.size());
            print_bitset_set2(choice2, lnodes_sorted.size() + rnodes_sorted.size());
        }
        
        nodeset rec_pair;
        rec_pair.set(node_index);
        
        int best_recpair_ind = -1;
        // Look to R nodes
        for (unsigned int r_ind = 0; r_ind < rnodes_sorted.size(); ++r_ind){
            if (*rnodes_sorted[r_ind].second->sites.begin() < *l->second->sites.rbegin()){
                choice2.set(r_ind+lnodes_sorted.size());
            }
            else if (l->second->recomb_right != -1 && *rnodes_sorted[r_ind].second->sites.begin() ==
                l->second->recomb_right){
                rec_pair.set(r_ind+lnodes_sorted.size());
            }
            else if (l->second->recomb_right != -1 && *rnodes_sorted[r_ind].second->sites.begin() < 
                l->second->recomb_right){
                best_recpair_ind = r_ind + lnodes_sorted.size();
            }
            else if (l->second->recomb_right == -1 || *rnodes_sorted[r_ind].second->sites.begin() > l->second->recomb_right){
                break;
            }
        }
        if (rec_pair.count() == 1 && best_recpair_ind != -1){
            //rec_pair.set(best_recpair_ind);
        }
        bool add = true;
        
        if (choice2.count() == 0){
            // One of the choices is to remove nothing. No need to store the choice then.
            add = false;
        }
        else{
            
            // Determine if there's a preexisting choice that's a subset of this one
            // (delete it) or this is a subset of a preexisting choice (don't add it)
            for (vector<pair<nodeset, nodeset > >::iterator prevchoice = 
                choices.begin(); prevchoice != choices.end(); ++prevchoice){
                
                if (prevchoice->second == choice2){
                    add = false;
                    prevchoice->first |= choice1;
                    break;
                }
                else if (prevchoice->second == choice1 && prevchoice->first == choice2){
                    add = false;
                    break;
                }
                else if (issubset_nodeset(prevchoice->first, choice1) &&
                    issubset_nodeset(prevchoice->second, choice2)){
                    // Prev is subset of current - make it current
                    prevchoice->first |= choice1;
                    prevchoice->second |= choice1;
                    add = false;
                    break;
                }
                else if (issubset_nodeset(prevchoice->second, choice1) &&
                    issubset_nodeset(prevchoice->first, choice2)){
                    // Prev is subset of current - make it current
                    prevchoice->second |= choice1;
                    prevchoice->first |= choice2;
                    add = false;
                    break;
                }
                
                else if (issubset_nodeset(choice1, prevchoice->first) &&
                    issubset_nodeset(choice2, prevchoice->second)){
                    // Current is subset of prev - don't add
                    add = false;
                    break;
                }
                else if (issubset_nodeset(choice1, prevchoice->second) &&
                    issubset_nodeset(choice2, prevchoice->first)){
                    // Current is subset of prev - don't add
                    add = false;
                    break;
                }
            }
            
            if (add){
                choices.push_back(make_pair(choice1, choice2));
            }
        }
        if (rec_pair.count() > 1){
            bool add = true;
            for (vector<nodeset >::iterator recpair_prev = recpairs.begin();
                recpair_prev != recpairs.end(); ++recpair_prev){
                if ((*recpair_prev & rec_pair).count() > 0){
                    *recpair_prev |= rec_pair;
                    add = false;
                    break;
                }
            }
            if (add){
                recpairs.push_back(rec_pair);
            }
        }
        ++node_index;
    }
    for (vector<pair<long int, arg_node*> >::iterator r = rnodes_sorted.begin();
        r != rnodes_sorted.end(); ++r){
        // Choice to remove current node
        nodeset choice1;
        choice1.set(node_index);
        // Choice to remove all incompatible nodes
        nodeset choice2;
        // Look left
        bool incompat_found = false;
        for (unsigned int r_ind = node_index - 1; r_ind >= lnodes_sorted.size(); --r_ind){
            if (incompat_found){
                choice2.set(r_ind);
            }
            else{
                if (!l_compat_r_conservative(rnodes_sorted[r_ind-lnodes_sorted.size()].second, r->second, false)){
                    choice2.set(r_ind);
                    incompat_found = true;
                }
            }
        }
        // Look right
        incompat_found = false;
        for (unsigned int r_ind = node_index + 1; r_ind < lnodes_sorted.size() + rnodes_sorted.size(); ++r_ind){
            if (incompat_found){
                choice2.set(r_ind);
            }
            else{
                if (!l_compat_r_conservative(r->second, rnodes_sorted[r_ind - lnodes_sorted.size()].second, false)){
                    choice2.set(r_ind);
                    incompat_found = true;
                }
            }
        }
        
        nodeset rec_pair;
        rec_pair.set(node_index);
        int best_recpair_ind = -1;
        // Look to L nodes
        for (int l_ind = lnodes_sorted.size()-1; l_ind >= 0; --l_ind){
            if (*lnodes_sorted[l_ind].second->sites.rbegin() > *r->second->sites.begin()){
                choice2.set(l_ind);
            }
            else if (r->second->recomb_left != -1 && r->second->recomb_left == *lnodes_sorted[l_ind].second->sites.rbegin()){
                rec_pair.set(l_ind);
            }
            else if (r->second->recomb_left != -1 && r->second->recomb_left < *lnodes_sorted[l_ind].second->sites.rbegin()){
                best_recpair_ind = l_ind;
            }
            else if (r->second->recomb_left == -1 || *lnodes_sorted[l_ind].second->sites.rbegin() < r->second->recomb_left){
                break;
            }
        }
        if (rec_pair.count() == 1 && best_recpair_ind != -1){
            //rec_pair.set(best_recpair_ind);
        }
        
        if (DEBUG_MODE){
            fprintf(stderr, "--prelim R choice:\n");
            print_bitset_set2(choice1, lnodes_sorted.size() + rnodes_sorted.size());
            print_bitset_set2(choice2, lnodes_sorted.size() + rnodes_sorted.size());
        }
        
        bool add = true;
        
        if (choice2.count() == 0){
            // One of the choices is to remove nothing. No need to store the choice then.
            add = false;
        }
        else{
            
            // Determine if there's a preexisting choice that's a subset of this one
            // (delete it) or this is a subset of a preexisting choice (don't add it)
            for (vector<pair<nodeset, nodeset > >::iterator prevchoice = 
                choices.begin(); prevchoice != choices.end(); ++prevchoice){
                
                
                if (prevchoice->second == choice2){
                    add = false;
                    prevchoice->first |= choice1;
                    break;
                }
                else if (prevchoice->first == choice2 && prevchoice->second == choice1){
                    add = false;
                    break;
                }
            
                else if (issubset_nodeset(prevchoice->first, choice1) &&
                    issubset_nodeset(prevchoice->second, choice2)){
                    // Prev is subset of current - make it current
                    prevchoice->first |= choice1;
                    prevchoice->second |= choice1;
                    add = false;
                    break;
                }
                else if (issubset_nodeset(prevchoice->second, choice1) &&
                    issubset_nodeset(prevchoice->first, choice2)){
                    // Prev is subset of current - make it current
                    prevchoice->second |= choice1;
                    prevchoice->first |= choice2;
                    add = false;
                    break;
                }
                
                else if (issubset_nodeset(choice1, prevchoice->first) &&
                    issubset_nodeset(choice2, prevchoice->second)){
                    // Current is subset of prev - don't add
                    add = false;
                    break;
                }
                else if (issubset_nodeset(choice1, prevchoice->second) &&
                    issubset_nodeset(choice2, prevchoice->first)){
                    // Current is subset of prev - don't add
                    add = false;
                    break;
                }
                
            }
            
            if (add){
                choices.push_back(make_pair(choice1, choice2));
            }
        }
        if (rec_pair.count() > 1){
            bool add = true;
            for (vector<nodeset >::iterator recpair_prev = recpairs.begin();
                recpair_prev != recpairs.end(); ++recpair_prev){
                if ((*recpair_prev & rec_pair).count() > 0){
                    *recpair_prev |= rec_pair;
                    add = false;
                    break;
                }
            }
            if (add){
                recpairs.push_back(rec_pair);
            }
        }
        ++node_index;
    }
    
    if (DEBUG_MODE){
        fprintf(stderr, "\n");
    }
    
    nodeset choice_mask;
    for (int i = 0; i < lnodes_sorted.size() + rnodes_sorted.size(); ++i){
        choice_mask.set(i);
    }
    
    if (choices.size() > 0){
        // Track every possible final combination of nodes.
        vector<nodeset > nodesets;
        
        if (DEBUG_MODE){
            fprintf(stderr, "L nodes:\n");
            for (unsigned int i = 0; i < lnodes_sorted.size(); ++i){
                fprintf(stderr, "[%d]\n", i);
                print_node_lite(lnodes_sorted[i].second, num_haplotypes);
            }
            fprintf(stderr, "\nR nodes:\n");
            for (unsigned int i = 0; i < rnodes_sorted.size(); ++i){
                fprintf(stderr, "[%d]\n", i + (int)lnodes_sorted.size());
                print_node_lite(rnodes_sorted[i].second, num_haplotypes);
            }
            fprintf(stderr, "\n");
            
            // Print everything.
            //fprintf(stderr, "UNFILTERED:\n");
            //print_graph_sif_lim(lv_grps_all, allnodes_l, allnodes_r);
            fprintf(stderr, "node:\n");
            print_node_lite(node, num_haplotypes);
            fprintf(stderr, "\nChoices (%ld):\n", choices.size());
        }
        
        nodeset_set nodesets_prelim;
        
        for (vector<pair<nodeset, nodeset > >::iterator choice = 
            choices.begin(); choice != choices.end(); ++choice){
            
            if (DEBUG_MODE){
                print_bitset_set2(choice->first, lnodes_sorted.size() + rnodes_sorted.size());
                print_bitset_set2(choice->second, lnodes_sorted.size() + rnodes_sorted.size());
            }
            
            if (nodesets_prelim.size() == 0){
                nodeset choice1 = choice->first;
                choice1.flip();
                choice1 &= choice_mask;
                nodeset choice2 = choice->second;
                choice2.flip();
                choice2 &= choice_mask;
                if (choice1.count() > 0){
                    nodesets_prelim.insert(choice1);
                }
                if (choice2.count() > 0){
                    nodesets_prelim.insert(choice2);
                }
            }
            else{
                nodeset_set newnodesets;
                for (nodeset_set::iterator ns = nodesets_prelim.begin();
                    ns != nodesets_prelim.end(); ++ns){
                    
                    if ((*ns & choice->first).count() > 0 && (*ns & choice->second).count() > 0){
                        nodeset nscpy = *ns;
                        
                        nodeset new1 = set_diff_nodeset(nscpy, choice->first);
                        nodeset new2 = set_diff_nodeset(nscpy, choice->second);
                        
                        if (new1.count() > 0){
                            newnodesets.insert(new1);
                        }
                        if (new2.count() > 0){
                            newnodesets.insert(new2);
                        }
                    
                    }
                    else{
                        newnodesets.insert(*ns);
                    }
                }
                
                nodesets_prelim = newnodesets;
                
            }
            
            if (DEBUG_MODE){
                for (unsigned int i = 0; i < lnodes_sorted.size() + rnodes_sorted.size(); ++i){
               
                    if (choice->first.test(i)){
                        if (i < lnodes_sorted.size()){
                            print_node_lite(lnodes_sorted[i].second, num_haplotypes);
                        }
                        else{
                            print_node_lite(rnodes_sorted[i-lnodes_sorted.size()].second, num_haplotypes);
                        }
                    }
                }
                fprintf(stderr, " OR \n");
                
                for (unsigned int i = 0; i < lnodes_sorted.size() + rnodes_sorted.size(); ++i){
                    if (choice->second.test(i)){
                        if (i < lnodes_sorted.size()){
                            print_node_lite(lnodes_sorted[i].second, num_haplotypes);
                        }
                        else{
                            print_node_lite(rnodes_sorted[i-lnodes_sorted.size()].second, num_haplotypes);
                        }
                    }
                }
                fprintf(stderr, "-\n");
            }
        }
        
        if (DEBUG_MODE){
            fprintf(stderr, "RECPAIRS:\n");
        }
        
        bool rp_ov = false;
        for (vector<nodeset >::iterator rp = recpairs.begin();
            rp != recpairs.end(); ++rp){
            
            if (DEBUG_MODE){
                print_bitset_set2(*rp, lnodes_sorted.size() + rnodes_sorted.size());
            }
            if (!rp_ov){
                for (vector<nodeset >::iterator rp2 = recpairs.begin();
                    rp2 != recpairs.end(); ++rp2){
                    if (rp != rp2 && (*rp & *rp2).count() > 0){
                        rp_ov = true;
                        break;
                    }
                }
            }
        }
        if (DEBUG_MODE){
            fprintf(stderr, "\n");
        }
       
        
        for (nodeset_set::iterator ns = nodesets_prelim.begin();
            ns != nodesets_prelim.end(); ++ns){
            
            nodesets.push_back(*ns);
        }
        
        
        long int nodeset_ind = 0;
        for (vector<nodeset >::iterator ns = nodesets.begin();
            ns != nodesets.end(); ){
            bool rm = false;
            
            if (DEBUG_MODE){
                fprintf(stderr, "PRE RECPAIRS:\n");
                print_bitset_set2(*ns, lnodes_sorted.size() + rnodes_sorted.size());
            }
            int match_count = 0;
            for (vector<nodeset >::iterator rp = recpairs.begin();
                rp != recpairs.end(); ++rp){
                if ((*ns & *rp).count() > 0 && !issuperset_nodeset(*ns, *rp)){
                    *ns = set_diff_nodeset(*ns, *rp);
                    //rm = true;
                    //break;
                    if (ns->count() == 0){
                        rm = true;
                        break;
                    }
                }
                else if (issuperset_nodeset(*ns, *rp)){
                    match_count++;
                }
            }
            
            if (rm){
                nodesets.erase(ns);
            }
            else{
                ++ns;
            }
            ++nodeset_ind;
        }
        
        // Filter to unique nodesets.
        for (vector<nodeset >::iterator ns1 = nodesets.begin();
            ns1 != nodesets.end();){
            bool unique = true;
            for (vector<nodeset >::iterator ns2 = nodesets.begin();
                ns2 != ns1; ++ns2){
                if (*ns2 == *ns1){
                    unique = false;
                    break;
                }
            }
            if (unique){
                ++ns1;
            }
            else{
                nodesets.erase(ns1);
            }
        }
        
        if (DEBUG_MODE){
            fprintf(stderr, "POST RECPAIRS:\n");
            for (vector<nodeset >::iterator ns = nodesets.begin();
                ns != nodesets.end(); ++ns){
                print_bitset_set2(*ns, lnodes_sorted.size() + rnodes_sorted.size());
            }
            fprintf(stderr, "\n");
        }
        
        // Track whether a nodeset covers L and R nodes in recombination
        nodeset lmask;
        nodeset rmask;
        for (unsigned int i = 0; i < lnodes_sorted.size(); ++i){
            lmask.set(i);
        }
        for (unsigned int i = lnodes_sorted.size(); i < lnodes_sorted.size() + rnodes_sorted.size(); ++i){
            rmask.set(i);
        }
        
        // Find out which nodes are excluded from all nodesets
        nodeset included;
        
        // Compute left & right indices for recombination according to each nodeset.
        map<int, pair<long int, long int> > nodesets_indices;
        map<int, pair<long int, long int> > nodesets_indices_max;
        
        
        //sort(nodesets.begin(), nodesets.end(), cladesetcomp_func);
        int ns_ind = 0;
        for (vector<nodeset >::iterator ns = nodesets.begin();
            ns != nodesets.end();){
            if (DEBUG_MODE){
                print_bitset_set2(*ns, lnodes_sorted.size() + rnodes_sorted.size());
            }
            
            long int left = -1;
            long int right = -1;
            long int leftmax = -1;
            long int rightmax = -1;
            
            arg_node* main_left_tmp = NULL;
            arg_node* main_right_tmp = NULL;
            
            long int fails_l = -1;
            long int fails_r = -1;
            
            bool erase = false;
            if ((*ns & lmask).count() == 0 || (*ns & rmask).count() == 0){
                erase = true;
            }
            else{

                for (unsigned int i = 0; i < lnodes_sorted.size() + rnodes_sorted.size(); ++i){
                    if (ns->test(i)){
                        if (i < lnodes_sorted.size()){
                            if (left == -1 || *lnodes_sorted[i].second->sites.rbegin() > left){
                                left = *lnodes_sorted[i].second->sites.rbegin();
                            }
                            if (leftmax == -1 || *lnodes_sorted[i].second->sites.rbegin() < leftmax){
                                leftmax = *lnodes_sorted[i].second->sites.rbegin();
                            }
                            
                        }
                        else{
                            if (right == -1 || *rnodes_sorted[i-lnodes_sorted.size()].second->sites.begin() < right){
                                right = *rnodes_sorted[i-lnodes_sorted.size()].second->sites.begin();
                            }
                            if (rightmax == -1 || *rnodes_sorted[i-lnodes_sorted.size()].second->sites.begin() > rightmax){
                                rightmax = *rnodes_sorted[i-lnodes_sorted.size()].second->sites.begin();
                            }
                        }
                    }
                }

                for (unsigned int i = 0; i < lnodes_sorted.size() + rnodes_sorted.size(); ++i){
                    if (ns->test(i)){
                        if (i < lnodes_sorted.size()){
                            if (lnodes_sorted[i].second->recomb_right != -1 &&
                                lnodes_sorted[i].second->recomb_right < right){
                                ns->reset(i);
                            }
                            else{
                                // Make the leftmost one main_left.
                                main_left_tmp = lnodes_sorted[i].second;
                            }
                        }
                        else{
                            if (rnodes_sorted[i-lnodes_sorted.size()].second->recomb_left != -1 &&
                                rnodes_sorted[i-lnodes_sorted.size()].second->recomb_left > left){
                                ns->reset(i);
                            }
                            else{
                                if (main_right_tmp == NULL){
                                    main_right_tmp = rnodes_sorted[i-lnodes_sorted.size()].second;
                                }
                                if (main_left_tmp != NULL && fails_r == -1){
                                    
                                    short comp = compare_bitsets(main_left_tmp->clade,
                                        rnodes_sorted[i-lnodes_sorted.size()].second->clade);
                                    
                                    if (comp == -1){
                                        fails_r = *rnodes_sorted[i-lnodes_sorted.size()].second->sites.begin();
                                    }
                                }
                            }
                        }
                    }
                }
                if ((*ns & lmask).count() == 0 || (*ns & rmask).count() == 0){
                    erase = true;
                }
                else{
                    if (main_left_tmp == NULL || main_right_tmp == NULL){
                        erase = true;
                    }
                    else{
                        if (DEBUG_MODE){
                            fprintf(stderr, "main_left null main_right null %d %d\n",
                                main_left_tmp == NULL, main_right_tmp == NULL);
                            print_node_lite(main_left_tmp, num_haplotypes);
                            print_node_lite(main_right_tmp, num_haplotypes);
                            fprintf(stderr, "fails_r %ld\n", fails_r);
                        }
                        // Make sure "fail_l" and "fail_r" can be obtained.
                        for (int i = (int)lnodes_sorted.size()-1; i >= 0; --i){
                            if (DEBUG_MODE){
                                fprintf(stderr, "testing node %d:\n", i);
                                print_node_lite(lnodes_sorted[i].second, num_haplotypes);
                            }
                            short comp = compare_bitsets(lnodes_sorted[i].second->clade,
                                main_right_tmp->clade);
                            if (comp == -1){
                                fails_l = *lnodes_sorted[i].second->sites.rbegin();
                                break;
                            }
                        }
                        if (fails_l == -1 || fails_r == -1){
                            if (DEBUG_MODE){
                                print_bitset_set2(*ns, lnodes_sorted.size() + rnodes_sorted.size());
                                fprintf(stderr, "fails_l %ld fails_r %ld\n", fails_l, fails_r);
                            }
                            if (!speculative){
                                erase = true;
                            }
                        }
                    }
                }
            }
            if (erase){
                nodesets.erase(ns);
            }
            else{
                included |= *ns;
 
                
                if (left != -1 && right != -1){
                    nodesets_indices.insert(make_pair(ns_ind, make_pair(left, right)));
                    nodesets_indices_max.insert(make_pair(ns_ind, make_pair(leftmax, rightmax)));
                }
                ++ns_ind;
                ++ns;
            }
            
        }
        
        if (nodesets.size() == 0){
            
            if (DEBUG_MODE){
                fprintf(stderr, "NO NODESETS\n");
            }
            
            return false;
        }
        
        if (DEBUG_MODE){
            // Print nodesets.
            fprintf(stderr, "NODESETS:\n");
            for (unsigned int i = 0; i < nodesets.size(); ++i){
                print_bitset_set2(nodesets[i], lnodes_sorted.size() + rnodes_sorted.size());
                fprintf(stderr, "== (%ld, %ld) / (%ld %ld): %ld ==\n", nodesets_indices[i].first, 
                    nodesets_indices[i].second, nodesets_indices_max[i].first, nodesets_indices_max[i].second, 
                    nodesets[i].count());
                fprintf(stderr, "\n");
            }
        }
        
        // Remove anything that's a subset of anything else.
        for (vector<nodeset >::iterator ns = nodesets.begin(); ns != nodesets.end();){
            
            bool rm = false;
            for (vector<nodeset >::iterator nsprev = nodesets.begin(); nsprev != ns; ++nsprev){
                if (issubset_nodeset(*ns, *nsprev)){
                    rm = true;
                    break;
                }
                else if (issubset_nodeset(*nsprev, *ns)){
                    *nsprev = *ns;
                    rm = true;
                    break;
                }
                else if (!force_solve){
                    if ((*ns & *nsprev).count() > 0 && !issubset_nodeset(*ns, *nsprev) &&
                        !issubset_nodeset(*nsprev, *ns)){

                        if (DEBUG_MODE){
                            fprintf(stderr, "Bailing out because of overlapping node sets:\n");
                            print_bitset_set2(*ns, lnodes_sorted.size() + rnodes_sorted.size());
                            print_bitset_set2(*nsprev, lnodes_sorted.size() + rnodes_sorted.size());
                        }
      
                        return false;
                    }
                }
            }
            if (rm){
                nodesets.erase(ns);
            }
            else{
                ++ns;
            }
        }
        
        int orig_index = -1;
        for (unsigned int i = 0; i < lnodes_sorted.size(); ++i){
            if (lnodes_sorted[i].second == node){
                orig_index = i;
                break;
            }
        }

        vector<nodeset > erased;
        for (vector<nodeset >::iterator ns = nodesets.begin(); ns != nodesets.end();){
            if (!ns->test(orig_index)){
                erased.push_back(*ns);
                nodesets.erase(ns);
            }
            else{
                ++ns;
            }
        }
        if (nodesets.size() == 0){
            
            return false;
        }
        if (nodesets.size() > 1 && erased.size() > 0){
            vector<nodeset > ns_filtered;
            for (vector<nodeset >::iterator ns = nodesets.begin();
                ns != nodesets.end(); ++ns){
                bool keep = true;
                for (vector<nodeset >::iterator e = erased.begin();
                    e != erased.end(); ++e){
                    if ((*e & *ns).count() > 0){
                        keep = false;
                        break;
                    }
                }
                if (keep){
                    ns_filtered.push_back(*ns);
                }
            }
            if (ns_filtered.size() > 0){
                nodesets = ns_filtered;
            }
        }
        
        set<int> ns_maxsize;
        long int max_size = -1;
        int chosen;
        for (unsigned int i = 0; i < nodesets.size(); ++i){
            if (DEBUG_MODE){
                fprintf(stderr, "%ld %ld\n", nodesets[i].count(), nodesets[i].size());
            }
            
            if (max_size == -1 || nodesets[i].count() > max_size){
                ns_maxsize.clear();
                ns_maxsize.insert(i);
                max_size = nodesets[i].count();
                
            }
            else if (nodesets[i].count() == (long unsigned int)max_size){
                ns_maxsize.insert(i);
            }
        }
        
        if (DEBUG_MODE){
            fprintf(stderr, "maxsize %ld; number %ld\n", max_size, ns_maxsize.size());
        }
        if (ns_maxsize.size() > 1){
            long int min_interval = -1;
            for (set<int>::iterator nsind = ns_maxsize.begin(); nsind != ns_maxsize.end(); ++nsind){
                if (min_interval == -1 || 
                    nodesets_indices_max[*nsind].second - nodesets_indices_max[*nsind].first < min_interval){
                    min_interval = nodesets_indices_max[*nsind].second - nodesets_indices_max[*nsind].first;
                    chosen = *nsind;
                }
            }
        }
        else{
            chosen = *ns_maxsize.begin();
        }
        
        if (DEBUG_MODE){
            fprintf(stderr, "CHOSEN:\n");
            print_bitset_set2(nodesets[chosen], lnodes_sorted.size() + rnodes_sorted.size());
        }
    
        // Translate chosen nodeset into nodes.
        allnodes_l.clear();
        allnodes_r.clear();
        
        arg_node* lmost = NULL;
        arg_node* rmost = NULL;
        
        for (unsigned int i = 0; i < lnodes_sorted.size() + rnodes_sorted.size(); ++i){
            if (nodesets[chosen].test(i)){
                if (i < lnodes_sorted.size()){
                    allnodes_l.insert(lnodes_sorted[i].second);
                    if (DEBUG_MODE){
                        fprintf(stderr, "L\n");
                        print_node_lite(lnodes_sorted[i].second, num_haplotypes);
                    }
                    if (main_left == NULL || *lnodes_sorted[i].second->sites.rbegin() > 
                        *main_left->sites.rbegin()){
                        main_left = lnodes_sorted[i].second;
                    }
                    if (lmost == NULL || *lnodes_sorted[i].second->sites.rbegin() < 
                        *lmost->sites.rbegin()){
                        lmost = lnodes_sorted[i].second;
                    }
                }
                else{
                    if (DEBUG_MODE){
                        fprintf(stderr, "R\n");
                        print_node_lite(rnodes_sorted[i-lnodes_sorted.size()].second, num_haplotypes);
                    }
                    allnodes_r.insert(rnodes_sorted[i-lnodes_sorted.size()].second);
                    if (main_right == NULL || *rnodes_sorted[i-lnodes_sorted.size()].second->sites.begin() < 
                        *main_right->sites.begin()){
                        main_right = rnodes_sorted[i-lnodes_sorted.size()].second;
                    }
                    if (rmost == NULL || *rnodes_sorted[i-lnodes_sorted.size()].second->sites.begin() >
                        *rmost->sites.begin()){
                        rmost = rnodes_sorted[i-lnodes_sorted.size()].second;
                    }
                }
            }
        }
        
        // Remove sorted L and R nodes that didn't make it in.
        for (vector<pair<long int, arg_node*> >::iterator l = lnodes_sorted.begin();
            l != lnodes_sorted.end();){
            if (allnodes_l.find(l->second) == allnodes_l.end()){
                lnodes_sorted.erase(l);
            }
            else{
                ++l;
            }  
        }
        for (vector<pair<long int, arg_node*> >::iterator r = rnodes_sorted.begin();
            r != rnodes_sorted.end();){
            if (allnodes_r.find(r->second) == allnodes_r.end()){
                rnodes_sorted.erase(r);
            }   
            else{
                ++r;
            }
        }
        
    }
    else{
        if (DEBUG_MODE){
            fprintf(stderr, "no choices to make\n");
        }
        
        main_left = lnodes_sorted.rbegin()->second;
        main_right = rnodes_sorted.begin()->second;
        
        for (set<arg_node*>::iterator l = allnodes_l.begin(); l != allnodes_l.end();){
            if ((*l)->recomb_right != -1 && (*l)->recomb_right < *main_right->sites.begin()){
                if (DEBUG_MODE){
                    fprintf(stderr, "removing (L):\n");
                    print_node_lite(*l, num_haplotypes);
                }
                for (vector<pair<long int, arg_node*> >::iterator ls = lnodes_sorted.begin();
                    ls != lnodes_sorted.end();){
                    if (ls->second == *l){
                        lnodes_sorted.erase(ls);
                        break;
                    }
                    else{
                        ++ls;
                    }  
                }
                allnodes_l.erase(l++);
            }
            else{
                ++l;
            }
        } 
        for (set<arg_node*>::iterator r = allnodes_r.begin(); r != allnodes_r.end();){
            if ((*r)->recomb_left != -1 && (*r)->recomb_left > *main_left->sites.rbegin()){
                if (DEBUG_MODE){
                    fprintf(stderr, "removing (R):\n");
                    print_node_lite(*r, num_haplotypes);
                }
                for (vector<pair<long int, arg_node*> >::iterator rs = rnodes_sorted.begin();
                    rs != rnodes_sorted.end();){
                    if (rs->second == *r){
                        rnodes_sorted.erase(rs);
                        break;
                    }
                    else{
                        ++rs;
                    }  
                }
                allnodes_r.erase(r++);
            }
            else{
                ++r;
            }
        }
        
        main_left = lnodes_sorted.rbegin()->second;
        main_right = rnodes_sorted.begin()->second;
        
        if (allnodes_l.size() == 0 || allnodes_r.size() == 0){
            
            //exit(1);
            return false;
        }
        
        if (DEBUG_MODE){
            fprintf(stderr, "L:\n");
            long int i = 0;
            for (vector<pair<long int, arg_node*> >::iterator l = lnodes_sorted.begin();
                l != lnodes_sorted.end(); ++l){
                fprintf(stderr, "[%ld]\n", i);
                print_node_lite(l->second, num_haplotypes);
                ++i;
            }
            fprintf(stderr, "R:\n");
            for (vector<pair<long int, arg_node*> >::iterator r = rnodes_sorted.begin();
                r != rnodes_sorted.end(); ++r){
                fprintf(stderr, "[%ld]\n", i);
                print_node_lite(r->second, num_haplotypes);
                ++i;
            }
        }
    }
    return true;
}

bool check_multi_recombs(set<arg_node*>& allnodes_l,
    set<arg_node*>& allnodes_r,
    set<arg_node*>& lv_grps_all,
    vector<pair<long int, arg_node*> >& lnodes_sorted,
    vector<pair<long int, arg_node*> >& rnodes_sorted,
    arg_node*& main_left,
    arg_node*& main_right,
    bool force_solve,
    arg_node*& candidate_other,
    int num_haplotypes){

    candidate_other = NULL;
    
    if (force_solve){
        return true;
    }
    
    arg_node* rmost = rnodes_sorted.rbegin()->second;
    arg_node* lmost = lnodes_sorted[0].second;
    
    short comp = compare_bitsets(main_left->clade, main_right->clade);
    if (comp != -1){
        if (DEBUG_MODE){
            fprintf(stderr, "MIDDLE COMPATIBLE\n");
        }
        // See if we might be dealing with multiple recombs.
        bool r_suspect = true;
        if (rmost == main_right || lmost == main_left){
            r_suspect = false;
        }
        else{
            for (vector<rc_edge>::iterator from = main_left->edges_from.begin();
                from != main_left->edges_from.end(); ++from){
                if (from->final_dest != NULL && allnodes_r.find(from->final_dest) != allnodes_r.end() &&
                    from->final_dest != rmost){
                    r_suspect = false;
                    break;
                }
            }
        }
        
        bool l_suspect = true;
        if (lmost == main_left || rmost == main_right){
            l_suspect = false;
        }
        else{
            for (vector<rc_edge>::iterator to = main_right->edges_to.begin();
                to != main_right->edges_to.end(); ++to){
                if (to->final_dest != NULL && allnodes_l.find(to->final_dest) != allnodes_l.end() &&
                    to->final_dest != lmost){
                    l_suspect = false;
                    break;
                }
            }
        
        }
        
        if (DEBUG_MODE){
            fprintf(stderr, "r suspect %d\n", r_suspect);
            fprintf(stderr, "l suspect %d\n", l_suspect);
            fprintf(stderr, "sizes %ld %ld\n", allnodes_l.size(), allnodes_r.size());
            fprintf(stderr, "%d %d\n", allnodes_l.find(main_left) != allnodes_l.end(),
                allnodes_r.find(rmost) != allnodes_r.end());
            fprintf(stderr, "%ld %ld\n", lnodes_sorted.size(), rnodes_sorted.size());
        }
        
        bool r_remove = false;
        bool l_remove = false;
        
        if (r_suspect && allnodes_r.size() > 1 && allnodes_l.size() > 0){
            
            long int alt_r_index = -1;
            for (vector<rc_edge>::iterator from = main_left->edges_from.begin(); from != main_left->edges_from.end(); ++from){
                for (vector<rc_edge>::iterator to = from->node->edges_to.begin(); to != from->node->edges_to.end(); ++to){
                    if (to->node == main_left ||
                        *to->node->sites.begin() > *main_left->sites.rbegin()){
                        for (vector<rc_edge>::iterator from2 = to->node->edges_from.begin();
                            from2 != to->node->edges_from.end(); ++from2){
                            if (from2->final_dest != NULL && *from2->final_dest->sites.begin() > *rmost->sites.begin()
                                ){
                                if (alt_r_index == -1 || *from2->final_dest->sites.begin() < alt_r_index){
                                    if (DEBUG_MODE){
                                        fprintf(stderr, "update alt_r_index\n");
                                        print_node_lite(from2->final_dest, num_haplotypes);
                                    }
                                    alt_r_index = *from2->final_dest->sites.begin();
                                }
                            }
                        }
                    }
                }
            }
            if (alt_r_index != -1){
                // Splitting the recombination event will create two new
                // intervals. If either is narrower than the current one,
                // make the split.
                
                vector<pair<long int, arg_node*> > l_sorted;
                for (set<arg_node*>::iterator l = allnodes_l.begin(); l != allnodes_l.end(); ++l){
                    l_sorted.push_back(make_pair(*(*l)->sites.rbegin(), *l));
                }
                sort(l_sorted.begin(), l_sorted.end());
                vector<pair<long int, arg_node*> > r_sorted;
                for (set<arg_node*>::iterator r = allnodes_r.begin(); r != allnodes_r.end(); ++r){
                    r_sorted.push_back(make_pair(*(*r)->sites.begin(), *r));
                }
                sort(r_sorted.begin(), r_sorted.end());
                
                long int interval1 = alt_r_index - *main_left->sites.rbegin();
                long int interval2 = *main_right->sites.begin() - *l_sorted[l_sorted.size()-2].second->sites.rbegin();
                if (DEBUG_MODE){
                    fprintf(stderr, "interval1 %ld %ld\n", *main_left->sites.rbegin(), alt_r_index);
                    fprintf(stderr, "interval2 %ld %ld\n", *l_sorted[l_sorted.size()-2].second->sites.rbegin(), *main_right->sites.begin());
                    fprintf(stderr, "cur interval %ld %ld\n", *main_left->sites.rbegin(), *main_right->sites.begin());
                }
                long int fail_l = -1;
                for (vector<pair<long int, arg_node*> >::reverse_iterator l = l_sorted.rbegin();
                    l != l_sorted.rend(); ++l){
                    short comp = compare_bitsets(l->second->clade, main_right->clade);
                    if (comp == -1){
                        fail_l = l->first;
                        break;
                    }
                }
                long int fail_r = -1;
                for (vector<pair<long int, arg_node*> >::iterator r = r_sorted.begin();
                    r != r_sorted.end(); ++r){
                    short comp = compare_bitsets(r->second->clade, main_left->clade);
                    if (comp == -1){
                        fail_r = r->first;
                        break;
                    }
                }
                long int interval3 = fail_r - fail_l;
                if (DEBUG_MODE){
                    fprintf(stderr, "conservative interval %ld %ld\n", fail_l, fail_r);
                    fprintf(stderr, "%ld & %ld vs %ld or %ld\n", interval1, interval2, interval3, *main_right->sites.begin()-*main_left->sites.rbegin());
                }
                //if (((float)interval1 + (float)interval2)/2 < interval3){
                if (!force_solve){
                    if (interval1 < interval3 || interval2 < interval3){
                        r_remove = true;
                    }
                }
                else if (interval1 < interval3 && interval2 < interval3){
                    r_remove = true;
                }
                
            }
        }
        if (!r_remove && l_suspect && allnodes_l.size() > 1 && allnodes_r.size() > 0){
            long int alt_l_index = -1;
            for (vector<rc_edge>::iterator to = main_right->edges_to.begin();
                to != main_right->edges_to.end(); ++to){
                for (vector<rc_edge>::iterator from = to->node->edges_from.begin();
                    from != to->node->edges_from.end(); ++from){
                    if (from->node == main_right ||
                        *from->node->sites.rbegin() < *main_right->sites.begin()){
                        for (vector<rc_edge>::iterator to2 = from->node->edges_to.begin();
                            to2 != from->node->edges_to.end(); ++to2){
                            
                            if (to2->final_dest != NULL && *to2->final_dest->sites.rbegin() < *lmost->sites.rbegin()
                                ){
                                if (alt_l_index == -1 || *to2->final_dest->sites.rbegin() > alt_l_index){
                                    if (DEBUG_MODE){
                                        fprintf(stderr, "update alt_l_index\n");
                                        print_node_lite(to2->final_dest, num_haplotypes);
                                    }
                                    alt_l_index = *to2->final_dest->sites.rbegin();
                                }
                            }
                        }
                    }
                }   
            }
            if (alt_l_index != -1){
                vector<pair<long int, arg_node*> > l_sorted;
                for (set<arg_node*>::iterator l = allnodes_l.begin(); l != allnodes_l.end(); ++l){
                    l_sorted.push_back(make_pair(*(*l)->sites.rbegin(), *l));
                }
                sort(l_sorted.begin(), l_sorted.end());
                vector<pair<long int, arg_node*> > r_sorted;
                for (set<arg_node*>::iterator r = allnodes_r.begin(); r != allnodes_r.end(); ++r){
                    r_sorted.push_back(make_pair(*(*r)->sites.begin(), *r));
                }
                sort(r_sorted.begin(), r_sorted.end());
                long int interval1 = *main_right->sites.begin() - alt_l_index;
                long int interval2 = *r_sorted[1].second->sites.begin() - *main_left->sites.rbegin();
                if (DEBUG_MODE){
                    fprintf(stderr, "interval1 %ld %ld\n", alt_l_index, *main_right->sites.begin());
                    fprintf(stderr, "interval2 %ld %ld\n", *main_left->sites.rbegin(), *r_sorted[1].second->sites.begin());
                    fprintf(stderr, "cur interval %ld %ld\n", *main_left->sites.rbegin(), *main_right->sites.begin());
                }
                long int fail_l = -1;
                for (vector<pair<long int, arg_node*> >::reverse_iterator l = l_sorted.rbegin();
                    l != l_sorted.rend(); ++l){
                    short comp = compare_bitsets(l->second->clade, main_right->clade);
                    if (comp == -1){
                        fail_l = l->first;
                        break;
                    }
                }
                long int fail_r = -1;
                for (vector<pair<long int, arg_node*> >::iterator r = r_sorted.begin();
                    r != r_sorted.end(); ++r){
                    short comp = compare_bitsets(r->second->clade, main_left->clade);
                    if (comp == -1){
                        fail_r = r->first;
                        break;
                    }
                }
                long int interval3 = fail_r - fail_l;
                if (DEBUG_MODE){
                    fprintf(stderr, "conservative interval %ld %ld\n", fail_l, fail_r);
                    fprintf(stderr, "%ld & %ld vs %ld or %ld\n", interval1, interval2, interval3, *main_right->sites.begin()-*main_left->sites.rbegin());
                }
                if (!force_solve){
                    if (interval1 < interval3 || interval2 < interval3){
                        l_remove = true;
                    }
                }
                else if (interval1 < interval3 && interval2 < interval3){
                    l_remove = true;
                }
            }
        }
        if (r_remove){
            if (!force_solve){
                candidate_other = main_left;
                return false;
            }
        }
        else if (l_remove){
            
            if (!force_solve){
                candidate_other = lmost;
                
                return false;
            }
        }
    }
    return true;
}

bool check_preexisting_recombs(set<arg_node*>& allnodes_l,
    set<arg_node*>& allnodes_r,
    vector<pair<long int, arg_node*> >& lnodes_sorted,
    vector<pair<long int, arg_node*> >& rnodes_sorted,
    set<arg_node*>& lv_grps_all,
    arg_node*& main_left,
    arg_node*& main_right,
    arg_clademap& lv_grps,
    vector<recomb_event>& recomb_catalog,
    arg_sitemap& sites_pos,
    arg_clademap& sites_clade,
    long int prop_dist,
    arg_node* root,
    vector<arg_node*>& newnodes,
    int num_haplotypes,
    set<arg_node*>& nodes_external){

    vector<pair<arg_node*, arg_node*> > failures;
    for (set<arg_node*>::iterator l = allnodes_l.begin(); l != allnodes_l.end(); ++l){
        for (vector<rc_edge>::iterator from = (*l)->edges_from.begin(); 
            from != (*l)->edges_from.end(); ++from){
            if (from->final_dest != NULL && allnodes_r.find(from->final_dest) != allnodes_r.end()){
                bool found = false;
                for (vector<pair<arg_node*, arg_node*> >::iterator p = failures.begin();
                    p != failures.end(); ++p){
                    if (p->first == *l && p->second == from->final_dest){
                        found = true;
                        break;
                    }
                }
                if (!found){
                    failures.push_back(make_pair(*l, from->final_dest));
                }
            }
        }
    }
    
    // Now see which are explained by preexisting recombination events.
    vector<pair<arg_node*, arg_node*> > failures_explained;
    vector<vector<int> > failures_rec_used;
    long int left = -1;
    long int right = -1;
    cladeset_set lv_grps_recomb;
    
    for (vector<pair<arg_node*, arg_node*> >::iterator f = failures.begin();
        f != failures.end(); ++f){
        
        long int l_limit = *f->first->sites.rbegin();
        long int r_limit = *f->second->sites.begin();
        
        // Adjust L clade
                    
        cladeset lclade = f->first->clade;
        for (vector<rc_edge>::iterator from = f->first->edges_from_solved.begin();
            from != f->first->edges_from_solved.end(); ++from){
            if (from->final_dest != NULL && *from->final_dest->sites.begin() < *f->second->sites.begin()){
                if ((from->type == 1 || from->type == 2) && issuperset_bitset(lclade, from->node->clade)){
                    lclade = set_diff_bitset(lclade, from->node->clade);
                    lv_grps_recomb.insert(from->node->clade);
                }
                else if (from->type == 3){
                    lclade |= from->node->clade;
                    lv_grps_recomb.insert(from->node->clade);
                }
            }
        }
        
        // Adjust R clade
        cladeset rclade = f->second->clade;
        for (vector<rc_edge>::iterator to = f->second->edges_to_solved.begin();
            to != f->second->edges_to_solved.end(); ++to){
            if (to->final_dest != NULL && *to->final_dest->sites.rbegin() >= *f->first->sites.rbegin()){
                if ((to->type == 2 || to->type == 3) && issuperset_bitset(rclade, to->node->clade)){
                    rclade = set_diff_bitset(rclade, to->node->clade);
                    lv_grps_recomb.insert(to->node->clade);
                }
                else if (to->type == 1){
                    rclade |= to->node->clade;
                    lv_grps_recomb.insert(to->node->clade);
                }
            }
        }
        
        short newcomp = compare_bitsets(lclade, rclade);

        vector<int> recombs_used_l;
        vector<int> recombs_used_r;
        
        if (newcomp == -1){

            lclade = adjust_clade_recomb(f->first, f->second, true, recomb_catalog, l_limit, lv_grps_recomb, recombs_used_l);
            
            rclade = adjust_clade_recomb(f->second, f->first, false, recomb_catalog, r_limit, lv_grps_recomb, recombs_used_r);
            
            newcomp = compare_bitsets(lclade, rclade);
        }
        
        if (newcomp != -1 && (lclade != f->first->clade || rclade != f->second->clade)){
            // This pair doesn't belong in the recombination.
            
            failures_explained.push_back(*f);
            if (left == -1 || l_limit > left){
                left = l_limit;
            }
            if (right == -1 || r_limit < right){
                right = r_limit;
            }
            
            // Combine recombs used
            vector<int> rec_used_both = recombs_used_l;
            for (vector<int>::iterator rr = recombs_used_r.begin(); rr != recombs_used_r.end(); ++rr){
                bool found = false;
                for (vector<int>::iterator rb = rec_used_both.begin(); 
                    rb != rec_used_both.end(); ++rb){
                    if (*rb == *rr){
                        found = true;
                        break;
                    }
                }
                if (!found){
                    rec_used_both.push_back(*rr);
                }
            }
            failures_rec_used.push_back(rec_used_both);
        }
    }
    
    if (failures_explained.size() > 0){
        
        
        // Remove anything that was already explained and add finished recomb edges.
        
        int fe_index = 0;
        
        for (vector<pair<arg_node*, arg_node*> >::iterator f = failures_explained.begin();
            f != failures_explained.end(); ++f){
            
            cladeset lv1 = set_diff_bitset(f->first->clade, f->second->clade);
            cladeset lv2 = set_int_bitset(f->first->clade, f->second->clade);
            cladeset lv3 = set_diff_bitset(f->second->clade, f->first->clade);
            
            set<pair<int, arg_node*> > solved_edges;
            
            bool got_new_edges = false;
            
            // Check previously solved recombs from edges
            for (vector<rc_edge>::iterator from = f->first->edges_from_solved.begin();
                from != f->first->edges_from_solved.end(); ++from){
            
                bool applies = false;
                int type = -1;
                arg_node* movednode = NULL;
                if (from->final_dest != NULL && *from->final_dest->sites.begin() <= *f->second->sites.begin()){
                    
                    if (from->node->clade == lv1){
                        type = 1;
                        applies = true;
                        movednode = from->node;
                    }
                    else if (from->node->clade == lv2){
                        type = 2;
                        applies = true;
                        movednode = from->node;
                    }
                    else if (from->node->clade == lv3){
                        type = 3;
                        applies = true;
                        movednode = from->node;
                    }
                }
                if (applies){
                    solved_edges.insert(make_pair(type, movednode));
                }
            }
            for (vector<rc_edge>::iterator to = f->second->edges_to_solved.begin();
                to != f->second->edges_to_solved.end(); ++to){
                
                bool applies = false;
                int type = -1;
                arg_node* movednode = NULL;
                if (to->final_dest != NULL && *to->final_dest->sites.rbegin() >= *f->first->sites.rbegin()){
                    
                    if (to->node->clade == lv1){
                        type = 1;
                        applies = true;
                        movednode = to->node;
                    }
                    else if (to->node->clade == lv2){
                        type = 2;
                        applies = true;
                        movednode = to->node;
                    }
                    else if (to->node->clade == lv3){
                        type = 3;
                        applies = true;
                        movednode = to->node;
                    }
                }
                
                if (applies){
                    solved_edges.insert(make_pair(type, movednode));
                }
            }
            
            for (set<pair<int, arg_node*> >::iterator solveddat = solved_edges.begin();
                solveddat != solved_edges.end(); ++solveddat){
                
                got_new_edges = true;
                add_solved_edges(f->first, f->second, solveddat->second, 
                    solveddat->first, true, newnodes, prop_dist, num_haplotypes);
            }
            
            
            // Check recomb catalog-derived previously solved recombs
            for (vector<int>::iterator fi = failures_rec_used[fe_index].begin();
                fi != failures_rec_used[fe_index].end(); ++fi){
                int type = -1;
                if (recomb_catalog[*fi].leaving == lv1){
                    type = 1;
                }
                else if (recomb_catalog[*fi].leaving == lv2){
                    type = 2;
                }
                else if (recomb_catalog[*fi].leaving == lv3){
                    type = 3;
                }
                // Get moved clade
                arg_node* movednode = NULL;
                for (arg_clademap::iterator n = sites_clade.equal_range(recomb_catalog[*fi].leaving).first;
                    n != sites_clade.equal_range(recomb_catalog[*fi].leaving).second; ++n){
                    if (site_ranges_overlap(recomb_catalog[*fi].left,
                        recomb_catalog[*fi].right,
                        n->second)){
                        movednode = n->second;
                        break;
                    }
                }
                if (type != -1 && movednode != NULL){
                    // Add solved edges.
                    got_new_edges = true;
                    add_solved_edges(f->first, f->second, movednode, type, true, 
                        newnodes, prop_dist, num_haplotypes);
                }
                else{
                    got_new_edges = true;
                    add_unsolved_edges(f->first, f->second, num_haplotypes);
                }
            }
            if (!got_new_edges){
                add_unsolved_edges(f->first, f->second, num_haplotypes);
            }
            del_connections_between(f->first, f->second, lv_grps, lv_grps_all, num_haplotypes);
            if (allnodes_l.find(f->first) != allnodes_l.end()){
                allnodes_l.erase(allnodes_l.find(f->first));
                for (vector<pair<long int, arg_node*> >::iterator ls = lnodes_sorted.begin();
                    ls != lnodes_sorted.end();){
                    if (ls->second == f->first){
                        lnodes_sorted.erase(ls);
                        break;
                    }
                    else{
                        ++ls;
                    }
                }
                // Expand range.
                long int newend = min(*f->first->sites.rbegin() + prop_dist,
                    left);
                //if (f->first->recomb_right != -1 && newend < f->first->recomb_right){
                //    newend = f->first->recomb_right-1;
                //}
                if (newend > f->first->end && (f->first->recomb_right == -1 || newend < f->first->recomb_right)){
                    
                    
                    expand_range_fixed(f->first, newend, 
                        sites_pos, sites_clade, lv_grps, prop_dist, root,
                        newnodes, false, num_haplotypes, true, nodes_external);
                    
                }
            }
            if (allnodes_r.find(f->second) != allnodes_r.end()){
                allnodes_r.erase(allnodes_r.find(f->second));
                for (vector<pair<long int, arg_node*> >::iterator rs = rnodes_sorted.begin();
                    rs != rnodes_sorted.end();){
                    if (rs->second == f->second){
                        rnodes_sorted.erase(rs);
                        break;
                    }
                    else{
                        ++rs;
                    }
                }
                // Expand range.
                long int newstart = max(*f->second->sites.begin() - prop_dist,
                   right);
                //if (f->second->recomb_left != -1 && newstart > f->second->recomb_left){
                //    newstart = f->second->recomb_left+1;
                //}
                if (newstart < f->second->start && (f->second->recomb_left == -1 || newstart > f->second->recomb_left)){
                    
                    expand_range_fixed(f->second, newstart, sites_pos, 
                        sites_clade, lv_grps, prop_dist, root, newnodes, false,
                        num_haplotypes, true, nodes_external);
                }
                
            }
            ++fe_index;
        }
        
        lv_grps_all.clear();
        for (set<arg_node*>::iterator l = allnodes_l.begin(); l != allnodes_l.end(); ++l){
            for (vector<rc_edge>::iterator from = (*l)->edges_from.begin(); from != (*l)->edges_from.end();
                ++from){
                if (from->final_dest != NULL && allnodes_r.find(from->final_dest) != allnodes_r.end()){
                    lv_grps_all.insert(from->node);
                }
            }
        }
        
        for (cladeset_set::iterator cl = 
            lv_grps_recomb.begin(); cl != lv_grps_recomb.end(); ++cl){
            for (set<arg_node*>::iterator lv = lv_grps_all.begin();
                lv != lv_grps_all.end();){
                if ((*lv)->clade == *cl){
                    lv_grps_all.erase(lv++);
                    break;
                }
                else{
                    ++lv;
                }
            }
        }
        
        // Remove any L & R nodes that no longer are connected by leaving 
        // groups in the remaining set.
        
        for (set<arg_node*>::iterator l = allnodes_l.begin(); l != allnodes_l.end();){
            bool connected = false;
            for (vector<rc_edge>::iterator from = (*l)->edges_from.begin();
                from != (*l)->edges_from.end(); ++from){
                if (lv_grps_all.find(from->node) != lv_grps_all.end() &&
                    from->final_dest != NULL && allnodes_r.find(from->final_dest) !=
                    allnodes_r.end()){
                    connected = true;
                    break;
                }
            }
            if (connected){
                ++l;
            }
            else{
                allnodes_l.erase(l++);
            }
        }
        for (set<arg_node*>::iterator r = allnodes_r.begin(); r != allnodes_r.end();){
            bool connected = false;
            for (vector<rc_edge>::iterator to = (*r)->edges_to.begin();
                to != (*r)->edges_to.end(); ++to){
                if (lv_grps_all.find(to->node) != lv_grps_all.end() &&
                    to->final_dest != NULL && allnodes_l.find(to->final_dest) != 
                    allnodes_l.end()){
                    connected = true;
                    break;
                }
            }
            if (connected){
                ++r;
            }
            else{
                allnodes_r.erase(r++);
            }
        }
        
    }
      
    for (vector<pair<long int, arg_node*> >::iterator ls = lnodes_sorted.begin();
        ls != lnodes_sorted.end();){
        if (allnodes_l.find(ls->second) != allnodes_l.end()){
            ++ls;
        }
        else{
            lnodes_sorted.erase(ls);
        }
    }
    for (vector<pair<long int, arg_node*> >::iterator rs = rnodes_sorted.begin();
        rs != rnodes_sorted.end();){
        if (allnodes_r.find(rs->second) != allnodes_r.end()){
            ++rs;
        }
        else{
            rnodes_sorted.erase(rs);
        }
    }
    
    if (allnodes_l.size() == 0 || allnodes_r.size() == 0 || lv_grps_all.size() == 0){
        main_left = NULL;
        main_right = NULL;
        
        return false;
    }
    else{
        main_left = lnodes_sorted.rbegin()->second;
        main_right = rnodes_sorted.begin()->second;
    }
    
    return true;
}

bool filter_recomb_nodes(set<arg_node*>& allnodes_l,
    set<arg_node*>& allnodes_r,
    vector<pair<long int, arg_node*> >& lnodes_sorted,
    vector<pair<long int, arg_node*> >& rnodes_sorted,
    set<arg_node*>& lv_grps_all,
    arg_node*& main_left,
    arg_node*& main_right,
    int num_haplotypes){
    //fprintf(stderr, "Rfrn\n");
    for (set<arg_node*>::iterator l = allnodes_l.begin(); l != allnodes_l.end();){
        bool connected_r = false;
        
        bool remove_recomb = false;
        if ((*l)->recomb_right != -1 && (*l)->recomb_right < *main_left->sites.rbegin()){
            remove_recomb = true;
        }
        else{
            for (vector<rc_edge>::iterator from = (*l)->edges_from.begin(); from != (*l)->edges_from.end();
                ++from){
                if (from->final_dest != NULL && allnodes_r.find(from->final_dest) != allnodes_r.end() &&
                    lv_grps_all.find(from->node) != lv_grps_all.end()){
                    connected_r = true;
                    break;
                }
            }
        }
        if (!connected_r || remove_recomb){
            for (vector<pair<long int, arg_node*> >::iterator ls = lnodes_sorted.begin();
                ls != lnodes_sorted.end();){
                if (ls->second == *l){
                    lnodes_sorted.erase(ls);
                    break;
                }
                else{
                    ++ls;
                }
            }
            allnodes_l.erase(l++);
        }
        else{
            ++l;
        }
    }

    if (allnodes_l.size() == 0 || lnodes_sorted.size() == 0){

        return false;
    }
    
    main_left = lnodes_sorted.rbegin()->second;

    for (set<arg_node*>::iterator r = allnodes_r.begin(); r != allnodes_r.end();){
        bool connected_l = false;
        
        bool remove_recomb = false;
        if ((*r)->recomb_left != -1 && (*r)->recomb_left > *main_right->sites.begin()){
            remove_recomb = true;
        }
        else{
            for (vector<rc_edge>::iterator to = (*r)->edges_to.begin(); to != (*r)->edges_to.end();
                ++to){
                if (to->final_dest != NULL && allnodes_l.find(to->final_dest) != allnodes_l.end() &&
                    lv_grps_all.find(to->node) != lv_grps_all.end()){
                    connected_l = true;
                    break;
                }
            }
        }
        if (!connected_l || remove_recomb){
            if (DEBUG_MODE){
                fprintf(stderr, "erasing R: %d %d\n", connected_l, remove_recomb);
                print_node_lite(*r, num_haplotypes);
            }
            for (vector<pair<long int, arg_node*> >::iterator rs = rnodes_sorted.begin();
                rs != rnodes_sorted.end();){
                if (rs->second == *r){
                    rnodes_sorted.erase(rs);
                    break;
                }
                else{
                    ++rs;
                }
            }
            allnodes_r.erase(r++);
        }
        else{
            ++r;
        }
    }

    if (allnodes_r.size() == 0 || rnodes_sorted.size() == 0){
        return false;
    }
    
    main_right = rnodes_sorted.begin()->second;

    for (set<arg_node*>::iterator lv = lv_grps_all.begin(); lv != lv_grps_all.end();){

        bool connected_l = false;

        for (vector<rc_edge>::iterator to = (*lv)->edges_to.begin(); to != (*lv)->edges_to.end();
            ++to){
            if (allnodes_l.find(to->node) != allnodes_l.end()){
                connected_l = true;
                break;
            }
        }
        bool connected_r = false;
        for (vector<rc_edge>::iterator from = (*lv)->edges_from.begin(); from != (*lv)->edges_from.end();
            ++from){
            if (allnodes_r.find(from->node) != allnodes_r.end()){
                connected_r = true;
                break;
            }
        }
        if (connected_l && connected_r){
            ++lv;
        }
        else{
            lv_grps_all.erase(lv++);
        }
    }
    if (lv_grps_all.size() == 0){
        return false;
    }
    return true;
}



void filter_lv_grps(set<arg_node*>& allnodes_l,
    set<arg_node*>& allnodes_r,    
    set<arg_node*>& lv_grps_all,
    vector<pair<int, arg_node*> >& lv_grps_degree,
    arg_sitemap& sites_pos,
    arg_clademap& sites_clade,
    arg_clademap& lv_grps,
    arg_node* root,
    arg_node* main_left,
    arg_node* main_right,
    long int prop_dist,
    int num_haplotypes){
    // Filter leaving groups.
    
    long int start = *main_left->sites.rbegin();
    long int end = *main_right->sites.begin();
    
    for (set<arg_node*>::iterator lv_grp = lv_grps_all.begin(); lv_grp != lv_grps_all.end();){
        
        //long int start = (*lv_grp)->start;
        //long int end = (*lv_grp)->end;
        
        if (DEBUG_MODE){
            fprintf(stderr, "lv grp:\n");
            print_node_lite(*lv_grp, num_haplotypes);
        }
        
        bool pass = true;
        bool deleted = false;
        
        if (has_node_invalidates(root, *main_left->sites.rbegin(), (*lv_grp)->clade)){
            pass = false;
        }
        else if (has_node_invalidates(root, *main_right->sites.begin(), (*lv_grp)->clade)){
            pass = false;
        }

        if (pass){
            
            pair<arg_clademap::iterator, arg_clademap::iterator> er = sites_clade.equal_range((*lv_grp)->clade);
            for (arg_clademap::iterator n = er.first; n != er.second; ++n){
                if (n->second->recomb_left != -1 && n->second->recomb_left >= start &&
                    n->second->recomb_left <= end){
                    pass = false;
                    break;
                }
                else if (n->second->recomb_right != -1 && n->second->recomb_right >= start &&
                    n->second->recomb_right <= end){
                    pass = false;
                    break;
                }
            }
            
            if (pass){
                for (arg_sitemap::iterator pos = sites_pos.begin();
                    pos != sites_pos.end(); ++pos){
                    if (pos->first >= start - prop_dist && pos->first <= end + prop_dist){
                        
                        if (site_ranges_overlap(start, end, pos->second)){
                            short comp = compare_bitsets((*lv_grp)->clade, pos->second->clade);
                            if (comp == -1){
                                pass = false;
                                break;
                            }
                        }
                    }
                    else if (pos->first > end + prop_dist){
                        break;
                    }
                    else if (!pass){
                        break;
                    }
                }
            }
            if (!pass){
                // It failed the test here; that means we need to delete it.
                if (DEBUG_MODE){
                    fprintf(stderr, "FILTERED OUT LV GRP\n");
                    print_node_lite(*lv_grp, num_haplotypes);
                }
                del_leaving_node(*lv_grp, num_haplotypes);
                pair<arg_clademap::iterator, arg_clademap::iterator> er = lv_grps.equal_range((*lv_grp)->clade);
                for (arg_clademap::iterator it = er.first; it != er.second;){
                    if (it->second == *lv_grp){
                        lv_grps.erase(it++);
                        //break;
                    }
                    else{
                        ++it;
                    }
                }
                deleted = true;
            }
        }
        if (pass){
            int degree_left = 0;
            int degree_right = 0;
            
            for (set<arg_node*>::iterator lnode = allnodes_l.begin(); lnode != allnodes_l.end();
                ++lnode){
                for (vector<rc_edge>::iterator from = (*lnode)->edges_from.begin();
                    from != (*lnode)->edges_from.end(); ++from){
                    
                    //if (from->node == *lv_grp){
                    if (from->node == *lv_grp && from->final_dest != NULL &&
                        allnodes_r.find(from->final_dest) != allnodes_r.end()){
                        degree_left++;
                    }
                }
            }
            for (set<arg_node*>::iterator rnode = allnodes_r.begin(); rnode != allnodes_r.end();
                ++rnode){
                for (vector<rc_edge>::iterator to = (*rnode)->edges_to.begin();
                    to != (*rnode)->edges_to.end(); ++to){
                    
                    //if (to->node == *lv_grp){
                    if (to->node == *lv_grp && to->final_dest != NULL && 
                        allnodes_l.find(to->final_dest) != allnodes_l.end()){
                        degree_right++;
                    }
                }
            }
            
            if (DEBUG_MODE && degree_left != degree_right){
                fprintf(stderr, "%d != %d\n", degree_left, degree_right);
                exit(1);
            }
            //lv_grps_degree.push_back(make_pair(-(degree_left+degree_right), *lv_grp));
            pair<int, arg_node*> p = make_pair(-max(degree_left, degree_right), *lv_grp);
            lv_grps_degree.push_back(p);
            //lv_grps_degree.push_back(make_pair(-max(degree_left, degree_right), *lv_grp));
        }
        
        if (deleted){
            delete *lv_grp;
            lv_grps_all.erase(lv_grp++);
        }
        else{
            ++lv_grp;
        }
    }
}

/**
 * true == managed to resolve a recomb 
 * false == did not
 */
bool resolve_recomb(arg_node* node,
    arg_clademap& lv_grps,
    arg_sitemap& sites_pos,
    arg_clademap& sites_clade,
    long int& prop_dist,
    arg_node* root,
    bool speculative,
    arg_node*& candidate_other,
    vector<recomb_event>& recomb_catalog,
    vector<arg_node*>& newnodes,
    int num_haplotypes,
    FILE* rcfile,
    long int hard_limit_l,
    long int hard_limit_r,
    set<arg_node*>& nodes_external){

    candidate_other = NULL;
    
    bool break_on_mistake = false;
    
    bool force_solve = false;
    
    short max_solve_attempts = 2;
    
    bool last_attempt = false;
    
    if (hard_limit_l != -1 && hard_limit_r != -1){
        force_solve = true;
    }
    else{
        if (node->solve_attempts == max_solve_attempts){
            return false;
            last_attempt = true;
            force_solve = true;
        }
        else{
            node->solve_attempts++;
            if (node->solve_attempts >= max_solve_attempts){
                force_solve = true;
            }
        }
    }
    
    if (DEBUG_MODE){
        fprintf(stderr, "**RESOLVE RECOMB**\n");
        print_node_lite(node, num_haplotypes);
        //print_recombs_right(node);
        fprintf(stderr, "\n");
    }
            
    // Keep track of all recombination events ever solved.
    //static vector<recomb_event> recomb_catalog;
    
    // Do preliminary checking to see if we really need to solve this.
    if (!recomb_check_prelim(node, lv_grps, num_haplotypes)){
        if (DEBUG_MODE){
            fprintf(stderr, "Prelim -> no need to solve\n");
        }
        return false;
    }
    
    vector<pair<int, arg_node*> > lv_grps_degree;
    
    arg_node* leaving_chosen = NULL;
    
    set<arg_node*> lv_grps_all;
    
    set<arg_node*> allnodes_l;
    set<arg_node*> allnodes_r;
    
    arg_node* main_left = NULL;
    arg_node* main_right = NULL;
    
    bool not_orig_left = false;
    
    gather_recomb_neighborhood(node, allnodes_l, allnodes_r, lv_grps_all, lv_grps, num_haplotypes);
    
    
    if (hard_limit_l != -1 && hard_limit_r != -1){
        for (set<arg_node*>::iterator ln = allnodes_l.begin(); ln != allnodes_l.end();){
            if ((*ln)->recomb_right < hard_limit_l){
                allnodes_l.erase(ln++);
            }
            else{
                ++ln;
            }
        }
        for (set<arg_node*>::iterator rn = allnodes_r.begin(); rn != allnodes_r.end();){
            if (*(*rn)->sites.begin() > hard_limit_r){
                allnodes_r.erase(rn++);
            }
            else{
                if (main_right == NULL || *(*rn)->sites.begin() < *main_right->sites.begin()){
                    main_right = *rn;
                }
                ++rn;
            }
        }
   
        for (set<arg_node*>::iterator rn = allnodes_r.begin(); rn != allnodes_r.end();){
            if ((*rn)->recomb_left >= *main_right->sites.begin() || main_right->recomb_right <= *(*rn)->sites.begin()){
                allnodes_r.erase(rn++);
            }
            else{
                ++rn;
            }
        }
        
        for (set<arg_node*>::iterator ln = allnodes_l.begin(); ln != allnodes_l.end();){
            bool connected = false;
            for (vector<rc_edge>::iterator from = (*ln)->edges_from.begin(); from != (*ln)->edges_from.end();
                ++from){
                if (from->final_dest != NULL && allnodes_r.find(from->final_dest) != allnodes_r.end()){
                    connected = true;
                    break;
                }
            }
            if (connected){
                ++ln;
            }
            else{
                allnodes_l.erase(ln++);
            }
        }
        for (set<arg_node*>::iterator rn = allnodes_r.begin(); rn != allnodes_r.end();){
            bool connected = false;
            for (vector<rc_edge>::iterator to = (*rn)->edges_to.begin(); to != (*rn)->edges_to.end();
                ++to){
                if (to->final_dest != NULL && allnodes_l.find(to->final_dest) != allnodes_l.end()){
                    connected = true;
                    break;
                }
            }
            if (connected){
                ++rn;
            }
            else{
                allnodes_r.erase(rn++);
            }
        }

    }
    
    
    // If no nodes fall on the left side of the recombination event, delete all
    // edges that caused this function to be called on the main node, and
    // bail out.
    if (allnodes_l.size() == 0){
        if (DEBUG_MODE){
            fprintf(stderr, "init allnodes_l size == 0\n");
        }
        set<arg_node*> rm;
        for (vector<rc_edge>::iterator from = node->edges_from.begin(); from != node->edges_from.end(); ++from){
            if (from->final_dest != NULL){
                rm.insert(from->final_dest);
            }
        }
        for (set<arg_node*>::iterator rmn = rm.begin(); rmn != rm.end(); ++rmn){
            //del_connections_between(node, *rmn, lv_grps, lv_grps_all);
        }

        return false;
    }
    
    
    
    // Sort everything.
    vector<pair<long int, arg_node*> > lnodes_sorted;
    for (set<arg_node*>::iterator l = allnodes_l.begin(); l != allnodes_l.end(); ++l){
        lnodes_sorted.push_back(make_pair(*(*l)->sites.rbegin(), *l));
    }
    sort(lnodes_sorted.begin(), lnodes_sorted.end());
    vector<pair<long int, arg_node*> > rnodes_sorted;
    for (set<arg_node*>::iterator r = allnodes_r.begin(); r != allnodes_r.end(); ++r){
        rnodes_sorted.push_back(make_pair(*(*r)->sites.begin(), *r));
    }
    sort(rnodes_sorted.begin(), rnodes_sorted.end());
    
    if (hard_limit_l == -1 || hard_limit_r == -1){
        // Remove nodes on the far left & far right that recombine with something before
        // the other side of this recombination event.
        if (!filter_far_edges(allnodes_l, allnodes_r, lv_grps_all, lnodes_sorted, 
            rnodes_sorted, main_left, main_right, node, force_solve, sites_pos, lv_grps, num_haplotypes)){
            return false;
        }
    
        // Determine what nodes can coexist to describe the same recombination event.
        if (!choose_nodeset(allnodes_l, allnodes_r, lv_grps_all, lnodes_sorted, 
            rnodes_sorted, node, main_left, main_right, force_solve, speculative, num_haplotypes)){
            //if (last_attempt){
            //    give_up(allnodes_l, allnodes_l, lv_grps_all, lv_grps,
            //        sites_pos, sites_clade, prop_dist);
            //}

            return false;
        }
    }

      
    main_left = lnodes_sorted.rbegin()->second;
    main_right = rnodes_sorted.begin()->second;
    
    set<arg_node*> removed_multi_l;
    set<arg_node*> removed_multi_r;
    
    if (!check_multi_recombs(allnodes_l, allnodes_r, lv_grps_all, lnodes_sorted, rnodes_sorted, 
        main_left, main_right, force_solve, candidate_other, num_haplotypes)){

        return false;
    }
    
    not_orig_left = main_left != node;
    
    // If the "real" left node is closer to the right side of the matrix
    // and can still be affected by new nodes, wait to solve anything.
    if (!force_solve && not_orig_left && *main_left->sites.rbegin() > root->end - 3 * prop_dist){

        return false;
    }
    
    if (!check_preexisting_recombs(allnodes_l, allnodes_r, lnodes_sorted,
        rnodes_sorted, lv_grps_all, main_left, main_right, lv_grps, recomb_catalog,
        sites_pos, sites_clade, prop_dist, root, newnodes, num_haplotypes, nodes_external)){

        return false;
    }
    
    if (allnodes_l.size() == 0 || allnodes_r.size() == 0 || 
        lv_grps_all.size() == 0 ||
        main_left == NULL || main_right == NULL ||
        allnodes_l.find(main_left) == allnodes_l.end() ||
        allnodes_r.find(main_right) == allnodes_r.end()){

        return false;
    }
    
    if (DEBUG_MODE){
        fprintf(stderr, "limited:\n");
        print_graph_sif_lim(lv_grps_all, allnodes_l, allnodes_r, num_haplotypes);
        fprintf(stderr, "done\n");
    }

    if (!filter_recomb_nodes(allnodes_l, allnodes_r, lnodes_sorted, rnodes_sorted, 
        lv_grps_all, main_left, main_right, num_haplotypes)){

        return false;
    }

    merge_lv_nodes_recomb(lv_grps_all, lv_grps);
    
    if (!filter_far_edges(allnodes_l, allnodes_r, lv_grps_all, lnodes_sorted, 
        rnodes_sorted, main_left, main_right, node, force_solve, sites_pos, lv_grps, num_haplotypes)){

        return false;
    }
    
    
    if (DEBUG_MODE){
        fprintf(stderr, "allnodes_l\n");
        for (set<arg_node*>::iterator l = allnodes_l.begin(); l != allnodes_l.end(); ++l){
            print_node_lite(*l, num_haplotypes);
            //print_recombs_right(*l);
        }
        fprintf(stderr, "\nallnodes_r\n");
        for (set<arg_node*>::iterator r = allnodes_r.begin(); r != allnodes_r.end(); ++r){
            print_node_lite(*r, num_haplotypes);
            //print_recombs_left(*r);
        }
        fprintf(stderr, "lv (unfiltered)\n");
        for (set<arg_node*>::iterator lv = lv_grps_all.begin(); lv != lv_grps_all.end(); ++lv){
            print_node_lite(*lv, num_haplotypes);
        }
        fprintf(stderr, "\n");
    }
    
    filter_lv_grps(allnodes_l, allnodes_r, lv_grps_all, lv_grps_degree, sites_pos, 
        sites_clade, lv_grps, root, main_left, main_right, prop_dist, num_haplotypes);
    
    if (lv_grps_degree.size() == 0){
        if (DEBUG_MODE){
            fprintf(stderr, "NO ELIGIBLE LEAVING GROUPS\n");
        }
        // Unsolvable. Delete all connections between main_left & main_right.

        return false;
    }

    // Sort by decreasing order of degree.
    sort(lv_grps_degree.begin(), lv_grps_degree.end());

    if (DEBUG_MODE){
        fprintf(stderr, "eligible leaving groups:\n");
        for (vector<pair<int, arg_node*> >::iterator lvd = lv_grps_degree.begin();
            lvd != lv_grps_degree.end(); ++lvd){
            fprintf(stderr, "%d:\n", -lvd->first);
            print_node_lite(lvd->second, num_haplotypes);
            fprintf(stderr, "\n");
        }
    }
    
    // If any of these leaving groups already exist as nodes in the ARG and are
    // tied to solved recombination events that would explain this one,
    // we can delete connections between the main nodes and bail out here.
    for (vector<pair<int, arg_node*> >::iterator lvd = lv_grps_degree.begin();
        lvd != lv_grps_degree.end();){
        bool rmlv = false;
        
        pair<arg_clademap::iterator, arg_clademap::iterator> er = sites_clade.equal_range(lvd->second->clade);
        for (arg_clademap::iterator n = er.first; n != er.second; ++n){
            if (site_ranges_overlap(*main_left->sites.rbegin(), *main_right->sites.begin(),
                n->second)){
                
                if (DEBUG_MODE){
                    fprintf(stderr, "Testing lv grp:\n");
                    print_node_lite(lvd->second, num_haplotypes);
                    print_node_lite(n->second, num_haplotypes);
                }
                
                bool has_inner_lnode = false;
                bool has_inner_rnode = false;
                
                for (vector<rc_edge>::iterator to = n->second->edges_to_solved.begin();
                    to != n->second->edges_to_solved.end(); ++to){
                    if (to->node->start > *main_left->sites.rbegin() && to->node->end < *main_right->sites.rbegin()){
                    /*
                    if (*to->node->sites.rbegin() >= *main_left->sites.rbegin() && 
                        *to->node->sites.rbegin() < *main_right->sites.begin()){
                    */    
                        if (!l_compat_r(main_left, to->node) || !l_compat_r(to->node, main_right)){
                            if (DEBUG_MODE){
                                fprintf(stderr, "hiln\n");
                                print_node_lite(to->node, num_haplotypes);
                            }
                            has_inner_lnode = true;
                            break;
                        }
                    }
                }
                for (vector<rc_edge>::iterator from = n->second->edges_from_solved.begin();
                    from != n->second->edges_from_solved.end(); ++from){
                    if (from->node->end < *main_right->sites.begin() && from->node->start > *main_left->sites.rbegin()){
                    /*
                    if (*from->node->sites.begin() <= *main_right->sites.begin() && 
                        *from->node->sites.begin() > *main_left->sites.rbegin()){
                    */    
                        if (!l_compat_r(main_left, from->node) || !l_compat_r(from->node, main_right)){
                            if (DEBUG_MODE){
                                fprintf(stderr, "hirn %d\n", from->type);
                                print_node_lite(from->node, num_haplotypes);
                                print_node_lite(main_left, num_haplotypes);
                                print_node_lite(main_right, num_haplotypes);
                            }
                            has_inner_rnode = true;
                            break;
                        }
                    }
                }
                if (has_inner_lnode || has_inner_rnode){
                    
                    // This recombination event can explain the one we're currently
                    // trying to solve. Bail.
                    if (DEBUG_MODE){
                        fprintf(stderr, "HAS INNER RECOMBS\n");
                    }
                    rmlv = true;
                }
            }
        }
        if (rmlv){
            if (DEBUG_MODE){
                fprintf(stderr, "removing\n");
            }
            lv_grps_degree.erase(lvd);
        }
        else{
            ++lvd;
        }
    }
    
    bool create_parent_both = false;

    // Choose a winning leaving group -- only resolve if there are no ties.
    // To do: when we need to choose one and break a tie, do that with precomputed 
    // genome-wide clade counts.
    //if (true){
    if (lv_grps_degree.size() == 1 || lv_grps_degree[0].first != lv_grps_degree[1].first){
        leaving_chosen = lv_grps_degree[0].second;
        
        if (DEBUG_MODE){
            fprintf(stderr, "CHOSEN:\n");
            print_node_lite(leaving_chosen, num_haplotypes);
        }
        
        // If the left node is still within propagation distance of new sites to
        // be read in, don't solve it yet. Just delete stuff that we know isn't
        // part of the recombination event.
        // Don't do this if the last site is the last site in the data set.
        
        
    }
    else if (lv_grps_degree.size() >= 2 && lv_grps_degree[0].first == -1){
        // If there's a lateral choice, choose it but also make a note to create the "parent of all"
        // clades.
        create_parent_both = true;
        for (vector<pair<int, arg_node*> >::iterator lvd = lv_grps_degree.begin();
           lvd != lv_grps_degree.end(); ++lvd){
            for (vector<rc_edge>::iterator to = lvd->second->edges_to.begin();
                to != lvd->second->edges_to.end(); ++to){
                if (allnodes_l.find(to->node) != allnodes_l.end() && to->type == 2){
                    leaving_chosen = lvd->second;
                    break;
                }
            }    
        }
        
        if (leaving_chosen == NULL){
            // Choose randomly.
            leaving_chosen = lv_grps_degree[0].second;
        }
        
    }
    else{
        if (force_solve){
            map<arg_node*, float> maxdists;
            clade_maxdists(maxdists, allnodes_l, allnodes_r, lnodes_sorted, rnodes_sorted,
                 lv_grps_degree, sites_pos, prop_dist, num_haplotypes);
            float maxdist = -1;
            for (map<arg_node*, float>::iterator md = maxdists.begin(); md != maxdists.end(); ++md){
                if (maxdist == -1 || md->second > maxdist){
                    maxdist = md->second;
                    leaving_chosen = md->first;
                }
            }
        }
        else{

            return false;
        }
    }
     
    if (leaving_chosen != NULL){
        
        main_left = NULL;
        main_right = NULL;
        // Find any recombs not explained by it and remove those nodes from the graph.
        for (set<arg_node*>::iterator l = allnodes_l.begin(); l != allnodes_l.end();){
            bool connected = false;
            for (vector<rc_edge>::iterator from = (*l)->edges_from.begin(); from != (*l)->edges_from.end(); ++from){
                if (from->node == leaving_chosen && from->final_dest != NULL && allnodes_r.find(from->final_dest) != allnodes_r.end()){
                    connected = true;
                    break;
                }
            }
            if (!connected){
                if (!force_solve){
                    return false;
                }
                else{
                    if (DEBUG_MODE){
                        fprintf(stderr, "ELIM L\n");
                        print_node_lite(*l, num_haplotypes);
                    }
                    allnodes_l.erase(l++);
                }
            }
            else{
                if (main_left == NULL || *(*l)->sites.rbegin() > *main_left->sites.rbegin()){
                    main_left = *l;
                }
                ++l;
            }
        }
        for (set<arg_node*>::iterator r = allnodes_r.begin(); r != allnodes_r.end();){
            bool connected = false;
            for (vector<rc_edge>::iterator to = (*r)->edges_to.begin(); to != (*r)->edges_to.end(); ++to){
                if (to->node == leaving_chosen && to->final_dest != NULL && allnodes_l.find(to->final_dest) != allnodes_l.end()){
                    connected = true;
                    break;
                }
            }
            if (!connected){
                if (!force_solve){
                    return false;
                }
                else{
                    if (DEBUG_MODE){
                        fprintf(stderr, "ELIM R\n");
                        print_node_lite(*r, num_haplotypes);
                    }
                    allnodes_r.erase(r++);
                }
            }
            else{
                if (main_right == NULL || *(*r)->sites.begin() < *main_right->sites.begin()){
                    main_right = *r;
                }
                ++r;
            }
        }
        
        for (vector<pair<long int, arg_node*> >::iterator ls = lnodes_sorted.begin(); ls != lnodes_sorted.end();){
            if (allnodes_l.find(ls->second) != allnodes_l.end()){
                ++ls;
            }
            else{
                lnodes_sorted.erase(ls);
            }
        }
        for (vector<pair<long int, arg_node*> >::iterator rs = rnodes_sorted.begin(); rs != rnodes_sorted.end();){
            if (allnodes_r.find(rs->second) != allnodes_r.end()){
                ++rs;
            }
            else{
                rnodes_sorted.erase(rs);
            }
        }
        
        if (!filter_far_edges(allnodes_l, allnodes_r, lv_grps_all, lnodes_sorted, 
            rnodes_sorted, main_left, main_right, node, force_solve, sites_pos, lv_grps, num_haplotypes)){

            return false;
        }
        
        bool recombs_before = false;
        // Do this again but for nodes other than the furthest L and R.
        for (set<arg_node*>::iterator ln = allnodes_l.begin(); ln != allnodes_l.end();){
            bool rm = false;
            if ((*ln)->recomb_right != -1 && (*ln)->recomb_right < 
                //rnodes_sorted.begin()->second->start){
                *rnodes_sorted.begin()->second->sites.begin()){
                //rm = true;
                recombs_before = true;
                break;
            }
            if (rm){
                for (vector<pair<long int, arg_node*> >::iterator ls = lnodes_sorted.begin();
                    ls != lnodes_sorted.end();){
                    if (ls->second == *ln){
                        lnodes_sorted.erase(ls);
                        break;
                    }
                    else{
                        ++ls;
                    }
                }
                allnodes_l.erase(ln++);
            }
            else{
                ++ln;
            }
        }
        if (allnodes_l.size() == 0){

            return false;
        }
        for (set<arg_node*>::iterator rn = allnodes_r.begin(); rn != allnodes_r.end();){
            bool rm = false;
            if ((*rn)->recomb_left != -1 && (*rn)->recomb_left > 
                //lnodes_sorted.rbegin()->second->end){
                *lnodes_sorted.rbegin()->second->sites.rbegin()){
                //rm = true;
                recombs_before = true;
                break;
            }
            if (rm){
                for (vector<pair<long int, arg_node*> >::iterator rs = rnodes_sorted.begin();
                    rs != rnodes_sorted.end();){
                    if (rs->second == *rn){
                        rnodes_sorted.erase(rs);
                        break;
                    }
                    else{
                        ++rs;
                    }
                }
                allnodes_r.erase(rn++);
            }
            else{
                ++rn;
            }
        }
        if (allnodes_r.size() == 0){

            return false;
        }
        if (recombs_before){

            return false;
        }
        
        set<arg_node*> alpha_left;
        set<arg_node*> beta_left;
        set<arg_node*> alpha_right;
        set<arg_node*> beta_right;
        
        // Track all nodes on the left & right side of this recombination event
        // These are used to expand the range of things after recombination 
        // event.
        set<arg_node*> nodes_left;
        set<arg_node*> nodes_right;
        
        // Build a recombination event
        recomb_event r;
        r.left = -1;
        r.right = -1;
        r.leaving = leaving_chosen->clade;
        r.leaving_known = true;
        
        long int widest_left = -1;
        long int widest_right = -1;
        
        long int shortest_alpha_left = -1;
        long int shortest_beta_right = -1;
        
        long int left_multi = -1;
        long int right_multi = -1;
        
        // Track the union of all clades involved in recomb in case it's up or down
        cladeset union_all;
        
        set<arg_node*> lv_grps_delete;
        
        set<arg_node*> lv_alt_del;
        
        bool mid_compat = l_compat_r(main_left, main_right);
        
        // First, determine recomb range/coordinates.
        
        // Left side:
        for (vector<rc_edge>::iterator edge = leaving_chosen->edges_to.begin();
            edge != leaving_chosen->edges_to.end(); ++edge){
            if (allnodes_l.find(edge->node) != allnodes_l.end()){
                
                if (r.left == -1 || *(edge->node->sites.rbegin()) > r.left){
                    r.left = *(edge->node->sites.rbegin());
                }
                if (widest_left == -1 || *(edge->node->sites.rbegin()) < widest_left){
                    widest_left = *(edge->node->sites.rbegin());
                }
                if (edge->type == 1 || edge->type == 2){
                    if (shortest_alpha_left == -1 || *(edge->node->sites.rbegin()) > shortest_alpha_left){
                        shortest_alpha_left = *(edge->node->sites.rbegin());
                    }
                }
                if (mid_compat && edge->node != main_left && !l_compat_r(edge->node, main_right)){
                    if (left_multi == -1 || *edge->node->sites.rbegin() > left_multi){
                        left_multi = *edge->node->sites.rbegin();
                    }
                }
            }
        }
        
        // Right side:
        for (vector<rc_edge>::iterator edge = leaving_chosen->edges_from.begin();
            edge != leaving_chosen->edges_from.end(); ++edge){
            if (allnodes_r.find(edge->node) != allnodes_r.end()){
               
                if (r.right == -1 || *(edge->node->sites.begin()) < r.right){
                    r.right = *(edge->node->sites.begin());
                }
                if (widest_right == -1 || *(edge->node->sites.begin()) > widest_right){
                    widest_right = *(edge->node->sites.begin());
                }
                if (edge->type == 2 || edge->type == 3){
                    if (shortest_beta_right == -1 || *(edge->node->sites.begin()) < shortest_beta_right){
                        shortest_beta_right = *(edge->node->sites.begin());
                    }
                }
                if (mid_compat && edge->node != main_right && !l_compat_r(main_left, edge->node)){
                    if (right_multi == -1 || *edge->node->sites.begin() < right_multi){
                        right_multi = *edge->node->sites.begin();
                    }
                }
            }
        }
        
        set<pair<arg_node*, arg_node*> > partners;
        
        // Now that we have recomb coordinates, we can exclude using clades
        // that recombine with the clades closest to the center of the recomb.
        
        // Next, determine clades.
        // Left side:
        for (vector<rc_edge>::iterator edge = leaving_chosen->edges_to.begin();
            edge != leaving_chosen->edges_to.end();){
            if (allnodes_l.find(edge->node) != allnodes_l.end()){

                if (edge->type == 1 || edge->type == 2){
                    
                    if (!r.alpha_known || issubset_bitset(edge->node->clade, r.alpha)){
                        r.alpha_known = true;
                        r.alpha = edge->node->clade;
                    }
                    alpha_left.insert(edge->node);
                    nodes_left.insert(edge->node);
                    if (union_all.size() == 0){
                        union_all = edge->node->clade;
                    }
                    else{
                        union_all |= edge->node->clade;
                    }
                }
                else if (edge->type == 3){
                    beta_left.insert(edge->node);
                    nodes_left.insert(edge->node);
                    if (union_all.size() == 0){
                        union_all = edge->node->clade;
                    }
                    else{
                        union_all |= edge->node->clade;
                    }
                }
            }
            
            
            bool erase = false;
            
            if (*edge->node->sites.rbegin() <= r.left){
                for (vector<rc_edge>::iterator from = edge->node->edges_from.begin(); 
                    from != edge->node->edges_from.end(); ++from){
                    if (from->type == edge->type && from->node == leaving_chosen &&
                        from->final_dest != NULL && *from->final_dest->sites.begin() >= r.right){
                        
                        partners.insert(make_pair(edge->node, from->final_dest));
                        
                    }
                }
            }
            

            if (erase){
                //del_connections_from_node(edge->node, leaving_chosen, false, lv_grps_delete);
            
                leaving_chosen->edges_to.erase(edge);
            }
            else{
                ++edge;
            }
        }
        
        // Right side:
        for (vector<rc_edge>::iterator edge = leaving_chosen->edges_from.begin();
            edge != leaving_chosen->edges_from.end();){
            
            if (allnodes_r.find(edge->node) != allnodes_r.end()){
                
                if (edge->type == 2 || edge->type == 3){
                    if (!r.beta_known || issubset_bitset(edge->node->clade, r.beta)){
                        r.beta_known = true;
                        r.beta = edge->node->clade;
                    }
                    beta_right.insert(edge->node);
                    nodes_right.insert(edge->node);
                    if (union_all.size() == 0){
                        union_all = edge->node->clade;
                    }
                    else{
                        union_all |= edge->node->clade;
                    }
                }
                else if (edge->type == 1){
                    alpha_right.insert(edge->node);
                    nodes_right.insert(edge->node);
                    if (union_all.size() == 0){
                        union_all = edge->node->clade;
                    }
                    else{
                        union_all |= edge->node->clade;
                    }
                }
            }
            
            bool erase = false;
            if (erase){
                //del_connections_to_node(edge->node, leaving_chosen, false, lv_grps_delete);
                leaving_chosen->edges_from.erase(edge);
            }
            else{
                ++edge;
            }
        }
        
        set<pair<arg_node*, arg_node*> > partners_edges;
        for (set<pair<arg_node*, arg_node*> >::iterator p = partners.begin(); p != partners.end(); ++p){
            if (allnodes_l.find(p->first) == allnodes_l.end() || 
                allnodes_r.find(p->second) == allnodes_r.end()){
                if (DEBUG_MODE){
                    fprintf(stderr, "DELETING PARTNERS:\n");
                    print_node_lite(p->first, num_haplotypes);
                    print_node_lite(p->second, num_haplotypes);
                    fprintf(stderr, "\n");
                }

                del_connections_between(p->first, p->second, lv_grps, lv_grps_all, num_haplotypes);
                
                partners_edges.insert(*p);
            }
        }
        
        if (!r.alpha_known){
            r.down = true;
            r.alpha = union_all;
            r.alpha_known = true;
        }
        else if (!r.beta_known){
            r.up = true;
            r.beta = union_all;
            r.beta_known = true;
        }
        

        // Set recomb range to a conservative value.
        vector<pair<long int, arg_node*> > lnodes_sorted;
        for (set<arg_node*>::iterator al = allnodes_l.begin(); al != allnodes_l.end(); ++al){
            lnodes_sorted.push_back(make_pair(*(*al)->sites.rbegin(), *al));
        }
        sort(lnodes_sorted.begin(), lnodes_sorted.end());
        vector<pair<long int, arg_node*> > rnodes_sorted;
        for (set<arg_node*>::iterator ar = allnodes_r.begin(); ar != allnodes_r.end(); ++ar){
            rnodes_sorted.push_back(make_pair(*(*ar)->sites.begin(), *ar));
        }
        sort(rnodes_sorted.begin(), rnodes_sorted.end());
        
        
        // == NEW == //
        vector<pair<arg_node*, arg_node*> > need_solved_edges;
        vector<int> need_solved_edges_types;
        
        lv_grps_all.clear();
        
        for (set<arg_node*>::iterator ln = allnodes_l.begin(); ln != allnodes_l.end(); ++ln){
            set<arg_node*> lpartners;
            for (vector<rc_edge>::iterator from = (*ln)->edges_from.begin(); from != (*ln)->edges_from.end();
                ++from){
                if (from->final_dest != NULL && *from->final_dest->sites.begin() > r.right &&
                    allnodes_r.find(from->final_dest) == allnodes_r.end() && 
                    !(l_compat_r(rnodes_sorted.rbegin()->second, from->final_dest))){
                    lpartners.insert(from->final_dest);
                }
            }
            for (set<arg_node*>::iterator partner = lpartners.begin(); partner != lpartners.end();
                ++partner){
                
                del_connections_between(*ln, *partner, lv_grps, lv_grps_all, num_haplotypes);
                
                int type = -1;
                if (r.leaving == set_diff_bitset((*ln)->clade, (*partner)->clade)){
                    type = 1;
                }
                else if (r.leaving == set_int_bitset((*ln)->clade, (*partner)->clade)){
                    type = 2;
                }
                else if (r.leaving == set_diff_bitset((*partner)->clade, (*ln)->clade)){
                    type = 3;
                }
               
                need_solved_edges.push_back(make_pair(*ln, *partner));
                need_solved_edges_types.push_back(type);
                
                if ((*partner)->edges_to.size() == 0){
                    long int left = -1;
                    if ((*partner)->recomb_left > r.right){
                        for (arg_sitemap::iterator it = sites_pos.equal_range((*partner)->recomb_left).second;
                            it != sites_pos.end(); ++it){
                            if (it->first > (*partner)->recomb_left){
                                left = it->first;
                                break;
                            }
                        }
                    }
                    else{
                        left = r.right;
                    }
                    if (left < (*partner)->start && left >= *(*partner)->sites.begin() - prop_dist &&
                        ((*partner)->recomb_left == -1 || left > (*partner)->recomb_left)){
                  
                        long int newsite = left;
                        
                        expand_range_fixed(*partner, newsite, sites_pos, sites_clade,
                            lv_grps, prop_dist, root, newnodes, false, num_haplotypes, true,
                            nodes_external); 
                        
                    }
                }
                
            }
            // Shrink range to conservative value so we can expand soon
            if ((*ln)->end > *(*ln)->sites.rbegin()){
                long int origend = (*ln)->end;
                (*ln)->end = *(*ln)->sites.rbegin();
                shorten_range(*ln, (*ln)->end - origend,
                    sites_pos, sites_clade, lv_grps, prop_dist, num_haplotypes, root, newnodes, -1);
            }
        }
        for (set<arg_node*>::iterator rn = allnodes_r.begin(); rn != allnodes_r.end(); ++rn){
            set<arg_node*> rpartners;
            for (vector<rc_edge>::iterator to = (*rn)->edges_to.begin(); to != (*rn)->edges_to.end();
                ++to){
                if (to->final_dest != NULL && *to->final_dest->sites.rbegin() < r.left &&
                    allnodes_l.find(to->final_dest) == allnodes_l.end() &&
                    !(l_compat_r(to->final_dest, lnodes_sorted.begin()->second))){
                    rpartners.insert(to->final_dest);
                }
            }
            for (set<arg_node*>::iterator partner = rpartners.begin(); partner != rpartners.end();
                ++partner){

                del_connections_between(*partner, *rn, lv_grps, lv_grps_all, num_haplotypes);
                
                int type = -1;
                if (r.leaving == set_diff_bitset((*partner)->clade, (*rn)->clade)){
                    type = 1;
                }
                else if (r.leaving == set_int_bitset((*partner)->clade, (*rn)->clade)){
                    type = 2;
                }
                else if (r.leaving == set_diff_bitset((*rn)->clade, (*partner)->clade)){
                    type = 3;
                }
                
                need_solved_edges.push_back(make_pair(*partner, *rn));
                need_solved_edges_types.push_back(type);
                
                if ((*partner)->edges_from.size() == 0){
                    long int right = -1;
                    if ((*partner)->recomb_right < r.left){
                        for (arg_sitemap::reverse_iterator it = 
                            arg_sitemap::reverse_iterator(sites_pos.equal_range((*partner)->recomb_right).first);
                            it != sites_pos.rend(); ++it){
                            if (it->first < (*partner)->recomb_right){
                                right = it->first;
                                break;
                            }
                        }
                    }
                    else{
                        right = r.left;
                    }
                    if (right > (*partner)->end && right <= *(*partner)->sites.rbegin() + prop_dist &&
                        ((*partner)->recomb_right == -1 || right < (*partner)->recomb_right)){
                        long int newsite = right;
                        
                        expand_range_fixed(*partner, newsite, sites_pos, sites_clade,
                            lv_grps, prop_dist, root, newnodes, false, num_haplotypes, true,
                            nodes_external); 
                        
                    }
                }
                
            }
            // Shrink to its first site.
            if ((*rn)->start < *(*rn)->sites.begin()){
                long int origstart = (*rn)->start;
                (*rn)->start = *(*rn)->sites.begin();
                shorten_range(*rn, (*rn)->start - origstart,
                    sites_pos, sites_clade, lv_grps, prop_dist, num_haplotypes, root, newnodes, -1);
            }
        }
        // ========= //
        

        long int left_fails = -1;
        long int right_fails = -1;
        for (vector<pair<long int, arg_node*> >::reverse_iterator ln = lnodes_sorted.rbegin();
            ln != lnodes_sorted.rend(); ++ln){
            short comp = compare_bitsets(ln->second->clade, main_right->clade);
            if (comp == -1){
                left_fails = ln->first;
                break;
            }
            else if (!speculative){
                
                if (allnodes_l.find(ln->second) != allnodes_l.end()){
                    allnodes_l.erase(allnodes_l.find(ln->second));
                }
                if (nodes_left.find(ln->second) != nodes_left.end()){
                    nodes_left.erase(nodes_left.find(ln->second));
                }
                if (beta_left.find(ln->second) != beta_left.end()){
                //    beta_left.erase(beta_left.find(ln->second));
                }
                
            }
        }

        for (vector<pair<long int, arg_node*> >::iterator rn = rnodes_sorted.begin(); rn != rnodes_sorted.end();
            ++rn){
            short comp = compare_bitsets(rn->second->clade, main_left->clade);
            if (comp == -1){
                right_fails = rn->first;
                break;
            }
            else if (!speculative){
                
                if (allnodes_r.find(rn->second) != allnodes_r.end()){
                    allnodes_r.erase(allnodes_r.find(rn->second));
                }
                if (nodes_right.find(rn->second) != nodes_right.end()){
                    nodes_right.erase(nodes_right.find(rn->second));
                }
                if (alpha_right.find(rn->second) != alpha_right.end()){
                //    alpha_right.erase(alpha_right.find(rn->second));
                }
                
            }
        }

        
        
        if (r.left == -1 || r.right == -1){
            
            return false;
        }
        
        
        long int newleft = -1;
        long int newright = -1;
        
        if (speculative){
            
            long int l_lim = -1;
            long int l_lim2 = -1;
            long int r_lim = -1;
            for (set<arg_node*>::iterator rn = allnodes_r.begin(); rn != allnodes_r.end(); ++rn){
                set<arg_node*> ltmp;
                set<arg_node*> rtmp;
                set<arg_node*> lvtmp;
                
                gather_recomb_neighborhood(*rn, ltmp, rtmp, lvtmp, lv_grps, num_haplotypes);
                for (set<arg_node*>::iterator lt = ltmp.begin(); lt != ltmp.end(); ++lt){
                    if (l_lim2 == -1 || *(*lt)->sites.rbegin() < l_lim2){
                        l_lim2 = *(*lt)->sites.rbegin();
                    }
                }
                
                if (l_lim == -1 || (*rn)->recomb_left > l_lim){
                    l_lim = (*rn)->recomb_left;
                }
            }
            for (set<arg_node*>::iterator ln = allnodes_l.begin(); ln != allnodes_l.end(); ++ln){
                if (r_lim == -1 || (*ln)->recomb_right < r_lim){
                    r_lim = (*ln)->recomb_right;
                }
            }
            
            if (l_lim != -1 && l_lim < 0){
                fprintf(stderr, "??? %ld\n", l_lim);
                exit(1);
            }
            
            long int lstart = r.left;
            if (l_lim != -1 && l_lim > lstart){
                lstart = l_lim;
            }
            long int rstart = r.right;
            if (r_lim != -1 && r_lim < rstart){
                rstart = r_lim;
            }
            if (l_lim2 > lstart && l_lim2 < rstart){
                lstart = l_lim2;
            }
            
            // Find tentative midpoint between L and R clades.
            //float midpt = ((float)r.right - (float)r.left)/2 + (float)r.left;
            
            //long int newleft = r.left;
            //long int newright = r.right;
            

            if (l_lim > r_lim){
            
            }
            vector<long int> positions;
            arg_sitemap::iterator sp = sites_pos.find(lstart);
            if (sp != sites_pos.end()){
                while(sp->first == lstart && sp != sites_pos.end()){
                    sp++;
                }
                for (; sp != sites_pos.end(); ++sp){
                    if (sp->first > lstart && sp->first < rstart){
                        bool found = false;
                        for (vector<long int>::iterator p = positions.begin(); p != positions.end(); ++p){
                            if (*p == sp->first){
                                found = true;
                                break;
                            }
                        }
                        if (!found){
                            positions.push_back(sp->first);
                        }
                    }
                    else{
                        break;
                    }
                }
                
                if (positions.size() > 1){
                    float midpt = lstart + (float)(rstart-lstart)/2.0;
                    vector<long int> lpos;
                    vector<long int> rpos;
                    for (vector<long int>::iterator pos = positions.begin(); pos != positions.end(); ++pos){
                        if (*pos < midpt){
                            lpos.push_back(*pos);
                        }
                        else if (*pos > midpt){
                            rpos.push_back(*pos);
                        }
                        else{
                            if (!RAND || (float)rand() / (float)RAND_MAX < 0.5){
                                lpos.push_back(*pos);  
                            }
                            else{
                                rpos.push_back(*pos);
                            }
                        }
                    }
                    if (lpos.size() == 0){
                        newleft = *rpos.begin();
                        newright = rpos[1];
                    }
                    else if (rpos.size() == 0){
                        newleft = lpos[lpos.size()-2];
                        newright = lpos[lpos.size()-1];
                    }
                    else{
                        newleft = *lpos.rbegin();
                        newright = *rpos.begin();
                    }
                    
                }
                else if (positions.size() == 1){
                    // Randomly assign L & R
                    if (!RAND || (float)rand() / (float)RAND_MAX < 0.5){
                        //newleft = r.left;
                        newleft = lstart;
                        newright = positions[0];
                    }
                    else{
                        newleft = positions[0];
                        //newright = r.right;
                        newright = rstart;
                    }
                }
            }
            
            if (newleft != -1 && newright != -1){
                r.left = newleft;
                r.right = newright;
            }
        }
        else{
            r.left = left_fails;
            r.right = right_fails;
        }

        if (DEBUG_MODE){
            
            for (set<arg_node*>::iterator al1 = allnodes_l.begin(); al1 != allnodes_l.end(); ++al1){
                for (set<arg_node*>::iterator al2 = allnodes_l.begin(); al2 != allnodes_l.end(); ++al2){
                    if (al1 != al2){
                        short comp = compare_bitsets((*al1)->clade, (*al2)->clade);
                        if (comp == -1){
                            fprintf(stderr, "lnodes fail:\n");
                            print_node_lite(*al1, num_haplotypes);
                            print_node_lite(*al2, num_haplotypes);
                            exit(1);
                        }
                    }
                }
            }
            for (set<arg_node*>::iterator ar1 = allnodes_r.begin(); ar1 != allnodes_r.end(); ++ar1){
                for (set<arg_node*>::iterator ar2 = allnodes_r.begin(); ar2 != allnodes_r.end(); ++ar2){
                    if (ar1 != ar2){
                        short comp = compare_bitsets((*ar1)->clade, (*ar2)->clade);
                        if (comp == -1){
                            fprintf(stderr, "rnodes fail:\n");
                            print_node_lite(*ar1, num_haplotypes);
                            print_node_lite(*ar2, num_haplotypes);
                            exit(1);
                        }
                    }
                }
            }
        
            fprintf(stderr, "==RECOMB== alpha_known %d beta_known %d\n", r.alpha_known, r.beta_known);
            fprintf(stderr, "(%ld %ld)\n", r.left, r.right);
            print_bitset_set(r.alpha, num_haplotypes);
            print_bitset_set(r.beta, num_haplotypes);
            print_bitset_set(r.leaving, num_haplotypes);
            fprintf(stderr, "\nalpha left:\n");
            for (set<arg_node*>::iterator al = alpha_left.begin(); al != alpha_left.end();
                ++al){
                print_node_lite(*al, num_haplotypes);
            }
            fprintf(stderr, "beta left:\n");
            for (set<arg_node*>::iterator bl = beta_left.begin(); bl != beta_left.end();
                ++bl){
                print_node_lite(*bl, num_haplotypes);
            }
            fprintf(stderr, "----------\nalpha right:\n");
            for (set<arg_node*>::iterator ar = alpha_right.begin(); ar != alpha_right.end();
                ++ar){
                print_node_lite(*ar, num_haplotypes);
                
            }
            fprintf(stderr, "beta right:\n");
            for (set<arg_node*>::iterator br = beta_right.begin(); br != beta_right.end();
                ++br){
                print_node_lite(*br, num_haplotypes);
            }
            fprintf(stderr, "\n");
            
            
            fprintf(stderr, "main left\n");
            print_node_lite(main_left, num_haplotypes);
            fprintf(stderr, "main right\n");
            print_node_lite(main_right, num_haplotypes);
            
        }
        
        recomb_catalog.push_back(r);
        
        long int l_boundary = r.left + (long int)floor((float)(r.right-r.left)/2);
        long int r_boundary = r.right - (long int)ceil((float)(r.right-r.left)/2);
        
        // print to recombination file in human-readable format.
        fprintf(rcfile, "seq\t%ld\t%ld\t", r.left, r.right);
        set<unsigned int> s1 = bitset2set(r.alpha, num_haplotypes);
        for (set<unsigned int>::iterator n = s1.begin(); n != s1.end(); ++n){
            fprintf(rcfile, "%d,", *n);
        }
        fprintf(rcfile, "\t");
        set<unsigned int> s2 = bitset2set(r.beta, num_haplotypes);
        for (set<unsigned int>::iterator n = s2.begin(); n != s2.end(); ++n){
            fprintf(rcfile, "%d,", *n);
        }
        fprintf(rcfile, "\t");
        set<unsigned int> s3 = bitset2set(r.leaving, num_haplotypes);
        for (set<unsigned int>::iterator n = s3.begin(); n != s3.end(); ++n){
            fprintf(rcfile, "%d,", *n);
        }
        fprintf(rcfile, "\n");
        
        
        if (DEBUG_MODE){
            if (r.right == 0){
                exit(1);
            }
            
            if (r.alpha == r.beta){
                fprintf(stderr, "recombination alpha and beta clades are identical\n");
                exit(1);
            }
            if (r.left == r.right){
                fprintf(stderr, "recombination left and right indices match\n");
                exit(1);
            }
        }

        // seg fault is somewhere after here
        
        for (arg_clademap::iterator l = lv_grps.begin();
            l != lv_grps.end(); ++l){
            
            if (DEBUG_MODE){
                if (l->second->edges_from.size() == 0 || l->second->edges_to.size() == 0){
                    fprintf(stderr, "leaving group has empty edges!\n");
                    print_node_lite(l->second, num_haplotypes);
                    exit(1);
                }
            }
            for (vector<rc_edge>::iterator t = l->second->edges_to.begin(); t != l->second->edges_to.end(); ++t){
                bool found = false;
                for (vector<rc_edge>::iterator f = t->node->edges_from.begin(); f != t->node->edges_from.end();
                    ++f){
                    if (f->node == l->second){
                        found = true;
                        break;
                    }
                }
                if (DEBUG_MODE && !found){
                    fprintf(stderr, "missing connection\n");
                    print_node_lite(l->second, num_haplotypes);
                    print_node_lite(t->node, num_haplotypes);
                    exit(1);
                }
            }
            for (vector<rc_edge>::iterator f = l->second->edges_from.begin(); f != l->second->edges_from.end(); ++f){
                bool found = false;
                for (vector<rc_edge>::iterator t = f->node->edges_to.begin(); t != f->node->edges_to.end(); ++t){
                    if (t->node == l->second){
                        found = true;
                        break;
                    }
                }
                if (DEBUG_MODE && !found){
                    fprintf(stderr, "missing connection2\n");
                    print_node_lite(l->second, num_haplotypes);
                    print_node_lite(f->node, num_haplotypes);
                    exit(1);
                }
            }
        }
        
        long int midpt = ((long int) round(((float)r.right - (float)r.left)/2)) + r.left;

        // Now adjust the rest of the graph. 
        
        // Create new nodes, test for four hap test failure with children of old nodes,
        // and inform parents of old nodes about these new nodes
        
        // Create leaving group node.
        arg_node* lv = new arg_node;
        lv->clade = r.leaving;
        
        //lv->start = max(min(r.left, midpt-prop_dist), (long int) 1);
        //lv->end = max(r.right, midpt + prop_dist);
        
        lv->start = max((long int)1, r.left-prop_dist);
        lv->end = r.right + prop_dist;
        
        //lv->start = max(r.left - prop_dist, (long int)1);
        //lv->end = r.right + prop_dist;
        lv->is_leaving_grp = true;
        lv->site_support = false;
        lv->recomb_left = -1;
        lv->recomb_right = -1;
        lv->sites.insert(r.left);
        lv->sites_lvgrp.insert(r.left);
        //lv->sites.insert(r.right);
        lv->fixed_left = false;
        lv->fixed_right = false;
        
        long int cml = -1;
        for (set<arg_node*>::iterator ln = allnodes_l.begin(); ln != allnodes_l.end(); ++ln){
            if (cml == -1 || (*ln)->closest_mut_l() > cml){
                cml = (*ln)->closest_mut_l();
            }
        }
        long int cmr = -1;
        for (set<arg_node*>::iterator rn = allnodes_r.begin(); rn != allnodes_r.end(); ++rn){
            if (cmr == -1 || (*rn)->closest_mut_r() < cmr){
                cmr = (*rn)->closest_mut_r();
            }
        }
        lv->closest_mut_left = cml;
        lv->closest_mut_right = cmr;
        
        if (r.left == -1 || r.right == -1){
            fprintf(stderr, "???1\n");
            print_node_lite(lv, num_haplotypes);
            exit(1);
        }

        arg_node* lvstored1 = create_node(lv, sites_pos, sites_clade, lv_grps, 
            prop_dist, root, newnodes, num_haplotypes, nodes_external);
        
        /*
        if (lvstored1 == NULL){
            fprintf(stderr, "null lv node\n");
            exit(1);
            
        }
        */

        if (lvstored1 != lv){
            /*
            if (lvstored == NULL){
                fprintf(stderr, "error: could not create leaving node???\n");
                exit(1);
            }
            */
            delete lv;
        }

        arg_node* lv2 = new arg_node;
        lv2->clade = r.leaving;
        
        //lv2->start = max(min(r.left, midpt-prop_dist), (long int) 1);
        //lv2->end = max(r.right, midpt + prop_dist);
        
        lv2->start = max((long int)1, r.left-prop_dist);
        lv2->end = r.right + prop_dist;
        
        //lv->start = max(r.left - prop_dist, (long int)1);
        //lv->end = r.right + prop_dist;
        lv2->is_leaving_grp = true;
        lv2->site_support = false;
        lv2->recomb_left = -1;
        lv2->recomb_right = -1;
        //lv->sites.insert(r.left);
        lv2->sites.insert(r.right);
        lv2->sites_lvgrp.insert(r.right);
        lv2->fixed_left = false;
        lv2->fixed_right = false;
        
        lv2->closest_mut_left = cml;
        lv2->closest_mut_right = cmr;

        arg_node* lvstored2 = create_node(lv2, sites_pos, sites_clade, lv_grps, 
            prop_dist, root, newnodes, num_haplotypes, nodes_external);

        /*
        if (lvstored2 == NULL){
            fprintf(stderr, "null lv node\n");
            exit(1);
        }
        */
        
        if (lvstored2 != lv2){
            delete lv2;
        }
        
        bool unsolvable = false;
        arg_node* lvstored = NULL;
        if (lvstored1 == lvstored2){
            lvstored = lvstored1;
        }
        else if (lvstored1 == NULL && lvstored2 != NULL){
            lvstored = lvstored2;
        }
        else if (lvstored2 == NULL && lvstored1 != NULL){
            lvstored = lvstored1;
        }
        else if (lvstored1 != NULL && lvstored2 != NULL){
            if (RAND==1){
                if((float) rand() / (float)(RAND_MAX) < 0.5){
                    lvstored = lvstored1;
                }
                else{
                    lvstored = lvstored2;
                }
            }
            else{
                lvstored = lvstored2;
            }
        }
        else{
            // Both are null? No idea how this could happen.
            unsolvable = true;
        }
        
        // Create finished edges that were already determined
        int nse_index = 0;
        for (vector<pair<arg_node*, arg_node*> >::iterator p = need_solved_edges.begin();
            p != need_solved_edges.end(); ++p){
            int type = need_solved_edges_types[nse_index];
            if (type == -1){
                add_unsolved_edges(p->first, p->second, num_haplotypes);
                //p->first->edges_from_unsolvable.push_back(rc_edge(-1, NULL, p->second));
                //p->second->edges_to_unsolvable.push_back(rc_edge(-1, NULL, p->first));
            }
            else{
                add_solved_edges(p->first, p->second, lvstored, type, true, newnodes, prop_dist, num_haplotypes);
            }
             
            ++nse_index;
        }
        
        /**
         * Delete obsolete edges
         */
        for (set<pair<arg_node*, arg_node*> >::iterator p = partners_edges.begin(); p != partners_edges.end(); ++p){
            cladeset lprime = set_diff_bitset(p->first->clade, p->second->clade);
            cladeset dd = p->first->clade & p->second->clade;
            cladeset rprime = set_diff_bitset(p->second->clade, p->first->clade);
            
            if (lvstored != NULL && (lvstored->clade == lprime || lvstored->clade == dd || lvstored->clade == rprime)){
                // Need to add solved edges to make up for this.
                int type;
                if (issuperset_bitset(p->first->clade, lvstored->clade)){
                    // L alpha
                    if (issuperset_bitset(p->second->clade, lvstored->clade)){
                        // R beta
                        type = 2;
                    }
                    else{
                        // R alpha
                        type = 1;
                    }
                }
                else{
                    // L beta; R must be beta
                    type = 3;
                }
               
                add_solved_edges(p->first, p->second, lvstored, type, true, newnodes, prop_dist, num_haplotypes);
                
                set<arg_node*> dummy;
                
                del_connections_between(p->first, p->second, lv_grps, dummy, num_haplotypes);
            }
            else{
                add_unsolved_edges(p->first, p->second, num_haplotypes);
                
                // is this bad?
                set<arg_node*> dummy;
                del_connections_between(p->first, p->second, lv_grps, dummy, num_haplotypes);
            }
        }
        if (unsolvable){
            for (set<arg_node*>::iterator ln = nodes_left.begin(); ln != nodes_left.end(); ++ln){
                for (set<arg_node*>::iterator rn = nodes_right.begin(); rn != nodes_right.end(); ++rn){
                    if (compare_bitsets((*ln)->clade, (*rn)->clade)){
                        add_unsolved_edges((*ln), (*rn), num_haplotypes);
                    }
                }
            }
            return false;
        }
        
        /**
         * UPDATE INDICES
         */
        for (set<arg_node*>::iterator node_left = nodes_left.begin(); node_left != nodes_left.end();
            ++node_left){
            
      
            if ((*node_left)->recomb_right == -1 || r.right < (*node_left)->recomb_right){
                //(*node_left)->recomb_right = r.right;

            }

        }
        
        for (set<arg_node*>::iterator node_right = nodes_right.begin(); node_right != nodes_right.end();
                ++node_right){   
            if ((*node_right)->recomb_left == -1 || r.left > (*node_right)->recomb_left){
                //(*node_right)->recomb_left = r.left;
            }
        }
        
        /**
         * ADJUST RANGES
         */
         
        long int ltarget = r.left;
        if (speculative && newleft != -1){
            ltarget = newleft;
        }
        //ltarget = l_boundary;
        
        bool created_failures = false;
        
        for (set<arg_node*>::iterator node_left = nodes_left.begin(); node_left != nodes_left.end();
            ++node_left){
            
            // Delete any potential recombs that are explained by this one OR
            // that conflict with this one.
            
            
            // Make propagation stop at rightmost point in recomb
            
            if (ltarget > (*node_left)->end){
                
                if (DEBUG_MODE){
                    fprintf(stderr, "expanding L node range:\n");
                    print_node_lite(*node_left, num_haplotypes);
                }
                
                long int orig_end = (*node_left)->end;
                
               
                if ((*node_left)->recomb_right != -1 && (*node_left)->recomb_right <= r.left){
                    //(*node_left)->end = (*node_left)->recomb_right-1;
                }
                
                set<arg_node*> failures;
                compile_failures(*node_left, root, failures,
                    orig_end, ltarget, prop_dist, sites_pos, num_haplotypes, false, false);
                
                expand_range_fixed(*node_left, ltarget, sites_pos, sites_clade,
                    lv_grps, prop_dist, root, newnodes, false, num_haplotypes, true, nodes_external);
                
                if ((*node_left)->recomb_right < (*node_left)->end){
                    fprintf(stderr, "recomb ind messup4\n");
                    print_node_lite(*node_left, num_haplotypes);
                    exit(1);
                }
                
                if ((*node_left)->end == ltarget){
                    elim_impossible_lv_grps(lv_grps, *node_left, num_haplotypes);
                    
                }
                else if (DEBUG_MODE){
                    fprintf(stderr, "unable to expand L:\n");
                    print_node_lite(*node_left, num_haplotypes);
                    print_recomb(r, num_haplotypes);
                    print_node_recursive_pos(root, r.left, 0, num_haplotypes);
                    //exit(1);
                }
            }
        }
         
       
        long int rtarget = r.right;
        if (speculative && newright != -1){
            rtarget = newright;
        }
        //rtarget = r_boundary;
        
        for (set<arg_node*>::iterator node_right = nodes_right.begin(); node_right != nodes_right.end();
            ++node_right){
            
            // Make propagation stop at leftmost point in recomb
            
            if (rtarget < (*node_right)->start){

                long int orig_start = (*node_right)->start;
                
                if ((*node_right)->recomb_left != -1 && (*node_right)->recomb_left >= r.right){
                //    (*node_right)->start = (*node_right)->recomb_left+1;
                }
                
                if (DEBUG_MODE){
                    fprintf(stderr, "expanding R node range; orig %ld:\n", orig_start);
                    print_node_lite(*node_right, num_haplotypes);
                    //print_node_upward_onelevel(*node_right, 0, num_haplotypes);
                }
                
                set<arg_node*> failures;
                compile_failures(*node_right, root, failures,
                    rtarget, orig_start, prop_dist, sites_pos, num_haplotypes, false, false);
                
                
                expand_range_fixed(*node_right, rtarget, sites_pos, sites_clade,
                    lv_grps, prop_dist, root, newnodes, false, num_haplotypes, true, nodes_external);
                
                if ((*node_right)->recomb_left > (*node_right)->start){
                    fprintf(stderr, "recomb ind messup5\n");
                    print_node_lite(*node_right, num_haplotypes);
                    exit(1);
                }
                
                if ((*node_right)->start == rtarget){
                    elim_impossible_lv_grps(lv_grps, *node_right, num_haplotypes);
                }
                else if (DEBUG_MODE){
                    fprintf(stderr, "unable to expand R:\n");
                    print_node_lite(*node_right, num_haplotypes);
                    print_recomb(r, num_haplotypes);
                    print_node_recursive_pos(root, r.right, 0, num_haplotypes);
                    //exit(1);
                }
            }
        }   
        
        
        /**
         * ADD "SOLVED" EDGES
         */
        for (set<arg_node*>::iterator al = alpha_left.begin(); al != alpha_left.end(); ++al){

            for (set<arg_node*>::iterator ar = alpha_right.begin(); ar != alpha_right.end(); ++ar){
                
                if (compare_bitsets((*al)->clade, (*ar)->clade) == -1 && 
                    set_diff_bitset((*al)->clade, (*ar)->clade) == lvstored->clade){
                
                    add_solved_edges(*al, *ar, lvstored, 1, false, newnodes, prop_dist, num_haplotypes);
                }
                
            }
            for (set<arg_node*>::iterator br = beta_right.begin(); br != beta_right.end(); ++br){
                if (set_int_bitset((*al)->clade, (*br)->clade) == lvstored->clade){
          
                    add_solved_edges(*al, *br, lvstored, 2, false, newnodes, prop_dist, num_haplotypes);
                }
            }
        }
        
        for (set<arg_node*>::iterator br = beta_right.begin(); br != beta_right.end(); ++br){
            for (set<arg_node*>::iterator bl = beta_left.begin(); bl != beta_left.end(); ++bl){
                if (compare_bitsets((*bl)->clade, (*br)->clade) == -1 && 
                    set_diff_bitset((*bl)->clade, (*br)->clade) == lvstored->clade){
                    add_solved_edges(*bl, *br, lvstored, 3, false, newnodes, prop_dist, num_haplotypes);
                }
            }
        }

        /**
         * CREATE NEW NODES
         */
        for (set<arg_node*>::iterator al = alpha_left.begin(); al != alpha_left.end(); ++al){
            if (issuperset_bitset((*al)->clade, r.alpha) &&
                *(*al)->sites.rbegin() + prop_dist > r.right){
                // Create a new node.
                
                cladeset newclade = set_diff_bitset((*al)->clade, r.leaving);
                
                if (newclade.count() > 1){
                    
                    //long int newstart = r_boundary;
                    long int newstart = r.right;
                    //long int newstart = r.left + 1;
                    //long int newend = *(*al)->sites.rbegin() + prop_dist;
                    long int newend = (*al)->closest_mut_l() + prop_dist;
                    
                    if (false){
                    //if (speculative){
                        newstart = midpt+1;
                    }
                    
                    if (newend >= newstart &&
                        safe_to_create_node(newclade, sites_pos, sites_clade, r.right, prop_dist, num_haplotypes)){
                        
                        arg_node* newnode = new arg_node;
                        newnode->recomb_left = -1;
                        newnode->recomb_right = -1;
                        newnode->is_leaving_grp = false;
                        newnode->skip = false;
                        newnode->site_support = false;
                        newnode->clade = newclade;
                        newnode->sites.insert(r.right);
                        newnode->start = newstart;
                        newnode->fixed_left = false;
                        newnode->fixed_right = false;
                        newnode->closest_mut_left = (*al)->closest_mut_l();

                        newnode->end = newend;
                        arg_node* stored = create_node(newnode, sites_pos, 
                            sites_clade, lv_grps, prop_dist, root, newnodes,
                            num_haplotypes, nodes_external);
                        
                        if (stored != newnode){
                            delete newnode;
                        }
                        
                        if (DEBUG_MODE){
                            fprintf(stderr, "CREATED new L -> R clade:\n");
                            print_node_lite(stored, num_haplotypes);
                            fprintf(stderr, "orig node:\n");
                            print_node_lite(*al, num_haplotypes);
                            print_bitset_set(newclade, num_haplotypes);
                            
                            if (stored->clade != newclade){
                                fprintf(stderr, "CLADES DON'T MATCH\n");
                                exit(1);
                            }
                            
                        }
                        
                    }
                }
            }
        }
        for (set<arg_node*>::iterator bl = beta_left.begin(); bl != beta_left.end(); ++bl){
            
            cladeset newclade = (*bl)->clade | r.leaving;
            
            if (newclade.count() > 1){
                long int newstart = r.right;
                //long int newstart = r_boundary;
                //long int newend = *(*bl)->sites.rbegin() + prop_dist;
                long int newend = (*bl)->closest_mut_l() + prop_dist;
                
                if (newend >= newstart &&
                    safe_to_create_node(newclade, sites_pos, sites_clade, r.right, prop_dist, num_haplotypes)){
                    
                    arg_node* newnode = new arg_node;
                    newnode->recomb_left = -1;
                    newnode->recomb_right = -1;
                    newnode->is_leaving_grp = false;
                    newnode->skip = false;
                    newnode->site_support = false;
                    newnode->clade = newclade;
                    newnode->sites.insert(r.right);
                    newnode->start = newstart;
                    newnode->end = newend;
                    newnode->fixed_left = false;
                    newnode->fixed_right = false;
                    newnode->closest_mut_left = (*bl)->closest_mut_l();
                    arg_node* stored = create_node(newnode, sites_pos, 
                        sites_clade, lv_grps, prop_dist, root, newnodes, num_haplotypes, nodes_external);
                    
                    if (stored != newnode){
                        delete newnode;
                    }
                    
                    if (DEBUG_MODE){
                        fprintf(stderr, "CREATED new L -> R clade from beta L:\n");
                        print_node_lite(stored, num_haplotypes);
                        fprintf(stderr, "orig node:\n");
                        print_node_lite(*bl, num_haplotypes);
                        print_bitset_set(newclade, num_haplotypes);
                        
                        if (stored->clade != newclade){
                            fprintf(stderr, "CLADES DON'T MATCH\n");
                            exit(1);
                        }
                    }
                    
                }
            }
        }

        for (set<arg_node*>::iterator ar = alpha_right.begin(); ar != alpha_right.end(); ++ar){
            cladeset newclade = (*ar)->clade | r.leaving;
            
            if (newclade.count() > 1){
                //long int newstart = *(*ar)->sites.begin() - prop_dist;
                long int newstart = (*ar)->closest_mut_r() - prop_dist;
                if (newstart < 1){
                    newstart = 1;
                }
                long int newend = r.left;
                //long int newend = l_boundary;
                
                if (false){
                //if (speculative){
                    newend = midpt;
                }
                if (newend >= newstart &&
                    safe_to_create_node(newclade, sites_pos, sites_clade, r.left, prop_dist, num_haplotypes)){
                    
                    arg_node* newnode = new arg_node;
                    newnode->recomb_left = -1;
                    newnode->recomb_right = -1;
                    newnode->is_leaving_grp = false;
                    newnode->skip = false;
                    newnode->site_support = false;
                    newnode->clade = newclade;
                    newnode->sites.insert(r.left);
                    newnode->start = newstart;
                    newnode->end = newend;
                    newnode->fixed_left = false;
                    newnode->fixed_right = false;
                    newnode->closest_mut_right = (*ar)->closest_mut_r();
                    arg_node* stored = create_node(newnode, sites_pos, 
                        sites_clade, lv_grps, prop_dist, root, newnodes,
                        num_haplotypes, nodes_external);  
                    
                    if (stored != newnode){
                        delete newnode;
                    }
                    
                    if (DEBUG_MODE){
                        fprintf(stderr, "CREATED new L <- R clade from alpha R:\n");
                        print_node_lite(stored, num_haplotypes);
                        print_bitset_set(newclade, num_haplotypes);
                        fprintf(stderr, "orig node:\n");
                        print_node_lite(*ar, num_haplotypes);
                        
                        if (stored->clade != newclade){
                            fprintf(stderr, "CLADES DON'T MATCH\n");
                            exit(1);
                        }
                            
                    }   
                    
                }
            }
        }

        for (set<arg_node*>::iterator br = beta_right.begin(); br != beta_right.end(); ++br){
            if (issuperset_bitset((*br)->clade, r.beta) &&
                *(*br)->sites.begin() - prop_dist < r.left){
                
                cladeset newclade = set_diff_bitset((*br)->clade, r.leaving);
                
                //long int newstart = *(*br)->sites.begin() - prop_dist;
                long int newstart = (*br)->closest_mut_r() - prop_dist;
                if (newstart < 1){
                    newstart = 1;
                }
                long int newend = r.left;
                //long int newend = l_boundary;
                //long int newend = r.right-1;
                if (false){
                //if (speculative){
                    newend = midpt;
                }
                if (newend >= newstart &&
                    safe_to_create_node(newclade, sites_pos, sites_clade, r.left, prop_dist, num_haplotypes)){
                    // Create a new node.
                    
                    arg_node* newnode = new arg_node;
                    newnode->recomb_left = -1;
                    newnode->recomb_right = -1;
                    newnode->is_leaving_grp = false;
                    newnode->skip = false;
                    newnode->site_support = false;
                    newnode->clade = newclade;
                    newnode->sites.insert(r.left);
                    newnode->start = newstart;
                    newnode->end = newend;
                    newnode->fixed_left = false;
                    newnode->fixed_right = false;
                    newnode->closest_mut_right = (*br)->closest_mut_r();
                    arg_node* stored = create_node(newnode, sites_pos, 
                        sites_clade, lv_grps, prop_dist, root, newnodes,
                        num_haplotypes, nodes_external);  
                    
                    if (stored != newnode){
                        delete newnode;
                    }
                    
                    if (DEBUG_MODE){
                        fprintf(stderr, "CREATED new L <- R clade:\n");
                        print_node_lite(stored, num_haplotypes);
                        print_bitset_set(newclade, num_haplotypes);
                        fprintf(stderr, "orig node:\n");
                        print_node_lite(*br, num_haplotypes);
                        print_bitset_set(stored->clade, num_haplotypes);
                        print_bitset_set(newclade, num_haplotypes);
                        
                        if (stored->clade != newclade){
                            fprintf(stderr, "CLADES DON'T MATCH\n");
                            exit(1);
                        }
                
                    }        
                }
            }
        }
        
        
        
        if (DEBUG_MODE){
            // Check if any two left or right nodes fail 4 hap test with each other.
            for (set<arg_node*>::iterator nl1 = nodes_left.begin(); nl1 != nodes_left.end(); ++nl1){
                for (set<arg_node*>::iterator nl2 = nodes_left.begin(); nl2 != nodes_left.end(); ++nl2){
                    if (nl1 != nl2){
                        short comp = compare_bitsets((*nl1)->clade, (*nl2)->clade);
                        if (comp == -1){
                            fprintf(stderr, "Two left nodes fail the four haplotype test\n");
                            print_node_lite(*nl1, num_haplotypes);
                            print_node_lite(*nl2, num_haplotypes);
                            exit(1);
                        }
                    }
                }
            }
            for (set<arg_node*>::iterator nr1 = nodes_right.begin(); nr1 != nodes_right.end(); ++nr1){
                for (set<arg_node*>::iterator nr2 = nodes_right.begin(); nr2 != nodes_right.end(); ++nr2){
                    if (nr1 != nr2){
                        short comp = compare_bitsets((*nr1)->clade, (*nr2)->clade);
                        if (comp == -1){
                            fprintf(stderr, "Two right nodes fail the four haplotype test\n");
                            print_node_lite(*nr1, num_haplotypes);
                            print_node_lite(*nr2, num_haplotypes);
                            exit(1);
                        }
                    }
                }
            }
        }

        if ((r.up || (!r.up && !r.down && create_parent_both)) && lvstored != NULL && 
            safe_to_create_node(r.beta, sites_pos, sites_clade, r.right, prop_dist, num_haplotypes)){
            // Make right clade
            
            
            arg_node* newnode = new arg_node;
            newnode->recomb_left = -1;
            newnode->recomb_right = -1;
            newnode->is_leaving_grp = false;
            newnode->skip = false;
            newnode->site_support = false;
            newnode->clade = r.beta;
            newnode->sites.insert(r.right);
            newnode->fixed_left = false;
            newnode->fixed_right = false;
            //newnode->closest_mut_right = cmr;
            newnode->closest_mut_right = r.right;
            newnode->closest_mut_left = r.right;
            
            //newnode->start = r.right;
            newnode->start = r.right - prop_dist;
            if (newnode->start < 1){
                newnode->start = 1;
            }
            newnode->end = r.right + prop_dist;

            arg_node* stored = create_node(newnode, sites_pos, sites_clade, 
                lv_grps, prop_dist, root, newnodes, num_haplotypes, nodes_external);
            
            if (stored != newnode){
                delete newnode;
            }
            
            // Create recombination edges.
            for (set<arg_node*>::iterator al = alpha_left.begin(); al != alpha_left.end();
                ++al){
                
                //add_solved_edges(*al, stored, lvstored, 1, false, newnodes, prop_dist);
                
               
            }
            
            if (DEBUG_MODE){
                fprintf(stderr, "created new upward R node\n");
                print_node_lite(stored, num_haplotypes);
                
                if (stored->clade != r.beta){
                    fprintf(stderr, "CLADES DON'T MATCH\n");
                    exit(1);
                }
                
            }
            
                        
        }
        else if ((r.down || (!r.up && !r.down && create_parent_both)) && lvstored != NULL && 
            safe_to_create_node(r.alpha, sites_pos, sites_clade, r.left, prop_dist, num_haplotypes)){
            // Make left clade
            arg_node* newnode = new arg_node;
            newnode->recomb_left = -1;
            newnode->recomb_right = -1;
            newnode->is_leaving_grp = false;
            newnode->skip = false;
            newnode->site_support = false;
            newnode->clade = r.alpha;
            newnode->sites.insert(r.left);
            newnode->start = r.left - prop_dist;
            newnode->fixed_left = false;
            newnode->fixed_right = false;
            //newnode->closest_mut_left = cml;
            newnode->closest_mut_left = r.left;
            newnode->closest_mut_right = r.left;
            
            if (newnode->start < 1){
                newnode->start = 1;
            }
            //newnode->end = r.left;
            newnode->end = r.left + prop_dist;
            arg_node* stored = create_node(newnode, sites_pos, sites_clade, 
                lv_grps, prop_dist, root, newnodes, num_haplotypes, nodes_external);
            
            if (stored != newnode){
                delete newnode;
            }
            
            if (DEBUG_MODE){
                
                fprintf(stderr, "created new downward L node\n");
                print_node_lite(stored, num_haplotypes);
                
                if (stored->clade != r.alpha){
                    fprintf(stderr, "CLADES DON'T MATCH\n");
                    exit(1);
                }
            }
            
            
            // Create recombination edges.
            for (set<arg_node*>::iterator br = beta_right.begin(); br != beta_right.end();
                ++br){
                
                //add_solved_edges(stored, *br, lvstored, 3, false, newnodes, prop_dist);
                
            }
            
        } 
    }
    
    if (DEBUG_MODE){
        fprintf(stderr, "RESOLVE RECOMB DONE\n");
    }
    
    return true;
}

