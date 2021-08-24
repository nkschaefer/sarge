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
#include "sets.h"
#include "treeNode.h"
#include "serialize.h"
#include "common.h"
#include "argnode.h"
#include "debug.h"

using std::cout;
using std::endl;
using namespace std;

/**
 * Allow sorting of vertical (parent/child) argnode edges
 */
bool operator<(const vert_edge& e1, const vert_edge& e2){
    if (e1.end < e2.end){
        return true;
    }
    else if (e2.end < e1.end){
        return false;
    }
    else{
        if (e1.start < e2.start){
            return true;
        }
        else if (e2.start < e1.start){
            return false;
        }
        else{
            return e1.to < e2.to;
        }
    }
}

bool operator<(const arg_node& n1, const arg_node& n2){
    if (n1.end < n2.end){
        return true;
    }
    else if (n2.end < n1.end){
        return false;
    }
    else{
        if (n1.start < n2.start){
            return true;
        }
        return false;
    }
}

bool sort_nodes_size(const arg_node* n1, const arg_node* n2){
    if (n1->clade.count() == n2->clade.count()){
        return n1 < n2;
    }
    else{
        return n1->clade.count() > n2->clade.count();
    }
}

struct rcsolved_info{
    arg_node* left;
    arg_node* right;
    arg_node* lv;
    int type;
    rcsolved_info(arg_node* l, arg_node* r, arg_node* mv, int t){
        this->left = l;
        this->right = r;
        this->lv = mv;
        this->type = t;
    };
};

bool operator<(const rcsolved_info& rc1, const rcsolved_info& rc2){
    if (rc1.left < rc2.left){
        return true;
    }
    else if (rc2.left < rc1.left){
        return false;
    }
    else{
        if (rc1.right < rc2.right){
            return true;
        }
        else if (rc2.right < rc1.right){
            return false;
        }
        else{
            if (rc1.lv < rc2.lv){
                return true;
            }
            else if (rc2.lv < rc1.lv){
                return false;
            }
            else{
                return rc1.type < rc2.type;
            }
        }
    }
};

void transfer_edges_from(arg_node* node, arg_node* newnode, int num_haplotypes){
    newnode->fixed_right = false;
    // Regular edges
    for (vector<rc_edge>::iterator edge_from = node->edges_from.begin(); 
        edge_from != node->edges_from.end(); ){
        if (edge_from->final_dest != NULL){
            edge_from->final_dest->fixed_left = false;
            for (vector<rc_edge>::iterator fd_to = edge_from->final_dest->edges_to.begin();
                fd_to != edge_from->final_dest->edges_to.end(); ++fd_to){
                if (fd_to->final_dest == node){
                    fd_to->final_dest = newnode;
                }
            }
        }
        for (vector<rc_edge>::iterator edge_to = edge_from->node->edges_to.begin();
            edge_to != edge_from->node->edges_to.end();){
            if (edge_to->node == node && edge_to->type == edge_from->type){
                edge_from->node->edges_to.erase(edge_to);
                //break;
            }
            else{
                ++edge_to;
            }
        }

        // Copy edge.
        rc_edge ef(*edge_from);

        newnode->edges_from.push_back(ef);
        rc_edge et(ef.type, newnode, NULL);
        edge_from->node->edges_to.push_back(et);
        node->edges_from.erase(edge_from);
    }
    
    // Solved edges
    
    // Store new edges to be created.
    set<rcsolved_info> solved_edges_new;
    
    for (vector<rc_edge>::iterator edge_from = node->edges_from_solved.begin(); 
        edge_from != node->edges_from_solved.end(); ++edge_from){
        if (edge_from->final_dest != NULL){
            solved_edges_new.insert(rcsolved_info(node, edge_from->final_dest,
                edge_from->node, edge_from->type));
        }
        else{
            for (vector<rc_edge>::iterator to = edge_from->node->edges_to_solved.begin();
                to != edge_from->node->edges_to_solved.end(); ++to){
                if (to->node == node && to->type == edge_from->type && to->final_dest != NULL){
                    solved_edges_new.insert(rcsolved_info(to->final_dest, edge_from->node,
                        node, edge_from->type));
                }
            }
        }

    }
    
    // Create new edges.
    vector<arg_node*> dummy;
    for (set<rcsolved_info>::iterator rcs = solved_edges_new.begin(); rcs != solved_edges_new.end();){
        del_solved_connections_between(rcs->left, rcs->right);
        if (rcs->left == node){

            add_solved_edges(newnode, rcs->right, rcs->lv, rcs->type, false, dummy, (long int)-1, num_haplotypes);
        }
        else if (rcs->lv == node){

            add_solved_edges(rcs->left, rcs->right, newnode, rcs->type, false, dummy, (long int)-1, num_haplotypes);
        }
        //add_solved_edges(rcs->left, rcs->right, rcs->lv, rcs->type, false, dummy, (long int)-1);
        solved_edges_new.erase(rcs++);
    }
    
    // Unsolvable edges
    set<arg_node*> dests;
    for (vector<rc_edge>::iterator edge_from = node->edges_from_unsolvable.begin();
        edge_from != node->edges_from_unsolvable.end();){
        edge_from->final_dest->fixed_left = false;
        dests.insert(edge_from->final_dest);
        for (vector<rc_edge>::iterator edge_to = edge_from->final_dest->edges_to_unsolvable.begin();
            edge_to != edge_from->final_dest->edges_to_unsolvable.end();){
            if (edge_to->final_dest == node){
                edge_from->final_dest->edges_to_unsolvable.erase(edge_to);
            }
            else{
                ++edge_to;
            }
        }
        node->edges_from_unsolvable.erase(edge_from);
    }
    for (set<arg_node*>::iterator dest = dests.begin(); dest != dests.end(); ++dest){
        bool found = false;
        for (vector<rc_edge>::iterator from = newnode->edges_from_unsolvable.begin();
            from != newnode->edges_from_unsolvable.end(); ++from){
            if (from->final_dest == *dest){
                found = true;
                break;
            }
        }
        if (!found){
            // Copy edge.
            add_unsolved_edges(newnode, *dest, num_haplotypes);
        }
    }

}

void transfer_edges_to(arg_node* node, arg_node* newnode, int num_haplotypes){
    newnode->fixed_left = false;
    // Regular edges
    for (vector<rc_edge>::iterator edge_to = node->edges_to.begin(); 
        edge_to != node->edges_to.end(); ){
        
        if (edge_to->final_dest != NULL){
            edge_to->final_dest->fixed_right = false;
            for (vector<rc_edge>::iterator fd_from = edge_to->final_dest->edges_from.begin();
                fd_from != edge_to->final_dest->edges_from.end(); ++fd_from){
                if (fd_from->final_dest == node){
                    fd_from->final_dest = newnode;
                }
            }
        }
        
        for (vector<rc_edge>::iterator edge_from = edge_to->node->edges_from.begin();
            edge_from != edge_to->node->edges_from.end();){
            if (edge_from->node == node && edge_from->type == edge_to->type){
                edge_to->node->edges_from.erase(edge_from);
                //break;
            }
            else{
                ++edge_from;
            }
        }
        
        // Copy edge.
        rc_edge et(*edge_to);
        newnode->edges_to.push_back(et);
        
        rc_edge ef(et.type, newnode);
        et.node->edges_from.push_back(ef);
        
        node->edges_to.erase(edge_to);
    }
    
    // Solved edges
    // Store new edges to be created.
    set<rcsolved_info> solved_edges_new;
    
    for (vector<rc_edge>::iterator edge_to = node->edges_to_solved.begin(); 
        edge_to != node->edges_to_solved.end(); ++edge_to){
        if (edge_to->final_dest != NULL){
            solved_edges_new.insert(rcsolved_info(edge_to->final_dest,
                node, edge_to->node, edge_to->type));
        }
        else{
            for (vector<rc_edge>::iterator from = edge_to->node->edges_from_solved.begin();
                from != edge_to->node->edges_from_solved.end(); ++from){
                if (from->node == node && from->type == edge_to->type &&
                    from->final_dest != NULL){
                    solved_edges_new.insert(rcsolved_info(edge_to->node, 
                        from->final_dest, node, edge_to->type));
                }
            }
        }
    }
    
    // Create new edges.
    vector<arg_node*> dummy;
    for (set<rcsolved_info>::iterator rcs = solved_edges_new.begin(); rcs != solved_edges_new.end();){
        del_solved_connections_between(rcs->left, rcs->right);
        if (rcs->right == node){

            add_solved_edges(rcs->left, newnode, rcs->lv, rcs->type, false, dummy, (long int)-1, num_haplotypes);
        }
        else if (rcs->lv == node){

            add_solved_edges(rcs->left, rcs->right, newnode, rcs->type, false, dummy, (long int)-1, num_haplotypes);
        }
        //add_solved_edges(rcs->left, rcs->right, rcs->lv, rcs->type, false, dummy, (long int)-1);
        solved_edges_new.erase(rcs++);
    }
    
    // Unsolvable edges
    set<arg_node*> dests;
    for (vector<rc_edge>::iterator edge_to = node->edges_to_unsolvable.begin();
        edge_to != node->edges_to_unsolvable.end();){
        edge_to->final_dest->fixed_left = false;
        dests.insert(edge_to->final_dest);
        for (vector<rc_edge>::iterator edge_from = edge_to->final_dest->edges_from_unsolvable.begin();
            edge_from != edge_to->final_dest->edges_from_unsolvable.end();){
            if (edge_from->final_dest == node){
                edge_to->final_dest->edges_from_unsolvable.erase(edge_from);
            }
            else{
                ++edge_from;
            }
        }
        node->edges_to_unsolvable.erase(edge_to);
    }
    for (set<arg_node*>::iterator dest = dests.begin(); dest != dests.end(); ++dest){
        bool found = false;
        for (vector<rc_edge>::iterator to = newnode->edges_to_unsolvable.begin();
            to != newnode->edges_to_unsolvable.end(); ++to){
            if (to->final_dest == *dest){
                found = true;
                break;
            }
        }
        if (!found){
            // Copy edge.
            add_unsolved_edges(*dest, newnode, num_haplotypes);
        }
    }
}

void add_solved_edges(arg_node* left, 
    arg_node* right, 
    arg_node* lv, 
    int type, 
    bool create_nodes, 
    vector<arg_node*>& newnodes,
    long int prop_dist,
    int num_haplotypes){
    if (*left->sites.begin() > *right->sites.rbegin()){
        fprintf(stderr, "flipped order?!\n");
        print_node_lite(left, num_haplotypes);
        print_node_lite(right, num_haplotypes);
        exit(1);
    }
    
    if (compare_bitsets(left->clade, right->clade) != -1){
        fprintf(stderr, "Adding solved recombination edges between two tree-compatible nodes\n");
        print_node_lite(left, num_haplotypes);
        print_node_lite(right, num_haplotypes);
        exit(1);
    }
    
    if (left->sites.size() == 0){
        fprintf(stderr, "ERROR: left node\n");
        print_node_lite(left, num_haplotypes);
        exit(1);
    }
    else if (right->sites.size() == 0){
        fprintf(stderr, "ERROR: right node\n");
        print_node_lite(right, num_haplotypes);
        exit(1);
    }
    else if (lv->sites.size() == 0){
        fprintf(stderr, "ERROR: lv node\n");
        print_node_lite(lv, num_haplotypes);
        exit(1);
    }
    
    if (lv->parents.size() == 0){
        fprintf(stderr, "no parents; lv\n");
        print_node_lite(lv, num_haplotypes);
        exit(1);
    }
    
    left->fixed_right = false;
    right->fixed_left = false;
    
    bool connected = false;
    if (left->edges_from_solved.size() <= right->edges_from_solved.size()){
        for (vector<rc_edge>::iterator from = left->edges_from_solved.begin();
            from != left->edges_from_solved.end(); ++from){
            if (from->final_dest != NULL && from->final_dest == right &&
                from->type == type){
                connected = true;
                break;
            }
        }
    }
    else{
        for (vector<rc_edge>::iterator to = right->edges_to_solved.begin();
            to != right->edges_to_solved.end(); ++to){
            if (to->final_dest != NULL && to->final_dest == left &&
                to->type == type){
                connected = true;
                break;   
            }
        }
    }
    
    if (type == 1){
        if (set_diff_bitset(left->clade, right->clade) != lv->clade){
            fprintf(stderr, "wrong clade (1)\n");
            exit(1);
        }
    }
    else if (type == 2){
        if (set_int_bitset(left->clade, right->clade) != lv->clade){
            fprintf(stderr, "wrong clade (2)\n");
            exit(1);
        }
    }
    else if (type == 3){
        if (set_diff_bitset(right->clade, left->clade) != lv->clade){
            fprintf(stderr, "wrong clade (3)\n");
            exit(1);
        }
    }
    
    if (!connected){
        left->edges_from_solved.push_back(rc_edge(type, lv, right));
        lv->edges_to_solved.push_back(rc_edge(type, left, NULL));
        lv->edges_from_solved.push_back(rc_edge(type, right, NULL));
        right->edges_to_solved.push_back(rc_edge(type, lv, left));
    }
    
    
    if (create_nodes){
        
        bool cr = left->recomb_right >= *right->sites.begin();
        
        bool cl = right->recomb_left <= *left->sites.rbegin();
        
        if (type == 1){
            // new right: left minus moved
            if (cr && left->closest_mut_l() != -1 && 
                left->closest_mut_l() + prop_dist >= *right->sites.begin()){
                arg_node* newnode = new arg_node;
                newnode->is_leaving_grp = false;
                newnode->site_support = false;
                newnode->recomb_left = -1;
                newnode->recomb_right = -1;
                newnode->skip = false;
                newnode->clade = set_diff_bitset(left->clade, lv->clade);
                newnode->sites.insert(*right->sites.begin());
                newnode->start = *right->sites.begin();
                newnode->end = *left->sites.rbegin() + prop_dist;
                newnode->fixed_left = false;
                newnode->fixed_right = false;
                newnode->closest_mut_left = left->closest_mut_l();
            
                newnodes.push_back(newnode);   
            }
            
            // new left: right plus moved
            if (cl && right->closest_mut_r() != -1 &&
                right->closest_mut_r() - prop_dist <= *left->sites.rbegin()){
                arg_node* newnode = new arg_node;
                newnode->is_leaving_grp = false;
                newnode->site_support = false;
                newnode->recomb_left = -1;
                newnode->recomb_right = -1;
                newnode->skip = false;
                newnode->clade = right->clade | lv->clade;
                newnode->sites.insert(*left->sites.rbegin());
                newnode->start = max((long int)1, *right->sites.begin()-prop_dist);
                newnode->end = *left->sites.rbegin();
                newnode->fixed_left = false;
                newnode->fixed_right = false;
                newnode->closest_mut_right = right->closest_mut_r();
                
                newnodes.push_back(newnode);
            }
        }
        else if (type == 2){
            // new right: left minus moved
            if (cr && left->closest_mut_l() != -1 && 
                left->closest_mut_l() + prop_dist >= *right->sites.begin()){
                arg_node* newnode = new arg_node;
                newnode->is_leaving_grp = false;
                newnode->site_support = false;
                newnode->recomb_left = -1;
                newnode->recomb_right = -1;
                newnode->skip = false;
                newnode->clade = set_diff_bitset(left->clade, lv->clade);
                newnode->sites.insert(*right->sites.begin());
                newnode->start = *right->sites.begin();
                newnode->end = *left->sites.rbegin() + prop_dist;
                newnode->fixed_left = false;
                newnode->fixed_right = false;
                newnode->closest_mut_left = left->closest_mut_l();
                
                newnodes.push_back(newnode);
            }
            // new left: right minus moved
            if (cl && right->closest_mut_r() != -1 && 
                right->closest_mut_r() - prop_dist <= *left->sites.rbegin()){
                arg_node* newnode = new arg_node;
                newnode->is_leaving_grp = false;
                newnode->site_support = false;
                newnode->recomb_left = -1;
                newnode->recomb_right = -1;
                newnode->skip = false;
                newnode->clade = set_diff_bitset(right->clade, lv->clade);
                newnode->sites.insert(*left->sites.rbegin());
                newnode->start = max((long int)1, *right->sites.begin()-prop_dist);
                newnode->end = *left->sites.rbegin();
                newnode->fixed_left = false;
                newnode->fixed_right = false;
                newnode->closest_mut_right = right->closest_mut_r();
                
                newnodes.push_back(newnode);
            }
        }
        else if (type == 3){
            // new right: left plus moved
            if (cr && left->closest_mut_l() != -1 &&
                 left->closest_mut_l() + prop_dist >= *right->sites.begin()){
                arg_node* newnode = new arg_node;
                newnode->is_leaving_grp = false;
                newnode->site_support = false;
                newnode->recomb_left = -1;
                newnode->recomb_right = -1;
                newnode->skip = false;
                newnode->clade = left->clade | lv->clade;
                newnode->sites.insert(*right->sites.begin());
                newnode->start = *right->sites.begin();
                newnode->end = *left->sites.rbegin() + prop_dist;
                newnode->fixed_left = false;
                newnode->fixed_right = false;
                newnode->closest_mut_left = left->closest_mut_l();
                
                newnodes.push_back(newnode);
            }
            // new left: right minus moved
            if (cl && right->closest_mut_r() != -1 && 
                right->closest_mut_r() - prop_dist <= *left->sites.rbegin()){
                arg_node* newnode = new arg_node;
                newnode->is_leaving_grp = false;
                newnode->site_support = false;
                newnode->recomb_left = -1;
                newnode->recomb_right = -1;
                newnode->skip = false;
                newnode->clade = set_diff_bitset(right->clade, lv->clade);
                newnode->sites.insert(*left->sites.rbegin());
                newnode->start = max((long int)1, *right->sites.begin()-prop_dist);
                newnode->end = *left->sites.rbegin();
                newnode->fixed_left = false;
                newnode->fixed_right = false;
                newnode->closest_mut_right = right->closest_mut_r();
                
                newnodes.push_back(newnode);
            }
        }
        
    }
    
}

void add_unsolved_edges(arg_node* left, arg_node* right, int num_haplotypes){
    if (*left->sites.begin() > *right->sites.rbegin()){
        fprintf(stderr, "unsolvable: flipped order?!\n");
        print_node_lite(left, num_haplotypes);
        print_node_lite(right, num_haplotypes);
        exit(1);
    }
    
    bool found = false;
    if (left->edges_from_unsolvable.size() <= right->edges_from_unsolvable.size()){
        for (vector<rc_edge>::iterator from = left->edges_from_unsolvable.begin();
            from != left->edges_from_unsolvable.end(); ++from){
            if (from->final_dest == right){
                found = true;
                break;
            }
        }
    }
    else{
        for (vector<rc_edge>::iterator to = right->edges_to_unsolvable.begin();
            to != right->edges_to_unsolvable.end(); ++to){
            if (to->final_dest == left){
                found = true;
                break;
            }
        }
    }
    if (!found){
        left->edges_from_unsolvable.push_back(rc_edge(-1, NULL, right));
        right->edges_to_unsolvable.push_back(rc_edge(-1, NULL, left));
    }
}

bool ranges_overlap(long int start1, long int end1, long int start2, long int end2){
    if ((start1 == -1 && end1 == -1) || (start2 == -1 && end2 == -1)){
        return false;
    }
    else if (start1 == end1 && start2 == end2){
        return start1 == start2;
    }
    else if (start1 == end1){
        return start1 >= start2 && end1 <= end2;
    }
    else if (start2 == end2){
        return start2 >= start1 && end2 <= end1;
    }
    else{
        return min(end1, end2) - max(start1, start2) >= 0;
    }
}

bool node_ranges_overlap(arg_node* n1, arg_node* n2){
    if ((n1->start == -1 && n1->end == -1) || (n2->start == -1 && n2->end == -1)){
        return false;
    }
    else if (n1->start == n1->end && n2->start == n2->end){
        return n1->start == n2->start;
    }
    else if (n1->start == n1->end){
        return n1->start >= n2->start && n1->end <= n2->end;
    }
    else if (n2->start == n2->end){
        return n2->start >= n1->start && n2->end <= n1->end;
    }
    else{
        return min(n1->end, n2->end) - max(n1->start, n2->start) >= 0;
    }
}

bool site_ranges_overlap(long int start, long int end, arg_node* n){
    if ((start == -1 && end == -1) || (n->start == -1 && n->end == -1)){
        return false;
    }
    else if (start == end && n->start == n->end){
        return start == n->start;
    }
    else if (start == end){
        return start >= n->start && end <= n->end;
    }
    else if (n->start == n->end){
        return n->start >= start && n->end <= end;
    }
    else{
        return min(n->end, end) - max(n->start, start) >= 0;
    }
}

bool has_node_invalidates(arg_node* parent, long int pos, cladeset& clade){
    for (set<vert_edge>::reverse_iterator child = parent->children.rbegin();
        child != parent->children.rend(); ++child){
        if (child->end < pos){
            break;
        }
        else if (child->start <= pos && child->end >= pos){
            short comp = compare_bitsets(clade, child->to->clade);
            if (comp == -1){
                return true;
            }
            else if (comp == 2){
                if (has_node_invalidates(child->to, pos, clade)){
                    return true;
                }
            }
            else if (comp == 3){
                return false;
            }
        }
    }
    return false;
}

bool has_node(arg_node* parent, long int pos, cladeset& clade){
    for (set<vert_edge>::reverse_iterator child = parent->children.rbegin(); 
        child != parent->children.rend(); ++child){
        if (child->end < pos){
            break;
        }
        else if (child->start <= pos && child->end >= pos){
            short comp = compare_bitsets(clade, child->to->clade);
            if (comp == 3){
                return true;
            }
            else if (comp == 2){
                if (has_node(child->to, pos, clade)){
                    return true;
                }
            }
        }
    }
    return false;
}

arg_node* get_smallest_containing(arg_node* parent, long int pos, cladeset& clade){
    for (set<vert_edge>::reverse_iterator child = parent->children.rbegin();
        child != parent->children.rend(); ++child){
        if (child->end < pos){
            break;
        }
        else if (child->start <= pos && child->end >= pos){
            if (issuperset_bitset(child->to->clade, clade)){
                return get_smallest_containing(child->to, pos, clade);
            }
        }
    }
    return parent;
}

void del_leaving_node(arg_node* lv, int num_haplotypes){
    set<pair<arg_node*, arg_node*> > recpairs;
    
    for (vector<rc_edge>::iterator e1 = lv->edges_to.begin(); e1 != lv->edges_to.end();){
        for (vector<rc_edge>::iterator e2 = e1->node->edges_from.begin(); e2 != e1->node->edges_from.end();){
            bool e2_rm = false;
            if (e2->node == lv){
                e2_rm = true;
                // See if there's an alternate route to final dest.
                bool store_pair = true;
                for (vector<rc_edge>::iterator e3 = e1->node->edges_from.begin();
                    e3 != e1->node->edges_from.end(); ++e3){
                    if (e3 != e2 && e3->final_dest != NULL && e3->node != lv && 
                        e3->final_dest == e2->final_dest){
                        store_pair = false;
                        break;
                    }
                }
                if (store_pair){
                    recpairs.insert(make_pair(e1->node, e2->final_dest));
                }
            }
            if (e2_rm){
                e1->node->edges_from.erase(e2);
            }
            else{
                ++e2;
            }
        }
        lv->edges_to.erase(e1);
    }
    for (vector<rc_edge>::iterator e1 = lv->edges_from.begin(); e1 != lv->edges_from.end();){
        for (vector<rc_edge>::iterator e2 = e1->node->edges_to.begin(); e2 != e1->node->edges_to.end();){
            bool e2_rm = false;
            if (e2->node == lv){
                e2_rm = true;
                // See if there's an alternate route to final dest.
                bool store_pair = true;
                for (vector<rc_edge>::iterator e3 = e1->node->edges_to.begin();
                    e3 != e1->node->edges_to.end(); ++e3){
                    if (e3 != e2 && e3->final_dest != NULL && e3->node != lv &&
                        e3->final_dest == e2->final_dest){
                        store_pair = false;
                        break;
                    }
                }
                if (store_pair){
                    recpairs.insert(make_pair(e2->final_dest, e1->node));
                }
            }
            if (e2_rm){
                e1->node->edges_to.erase(e2);
            }
            else{
                ++e2;
            }
        }
        lv->edges_from.erase(e1);
    }
    
    // Create "unsolvable" edges for nodes that can't be connected in any other
    // way.
    for (set<pair<arg_node*, arg_node*> >::iterator rp = recpairs.begin();
        rp != recpairs.end(); ++rp){
        add_unsolved_edges(rp->first, rp->second, num_haplotypes);
    }
}

/**
 * Returns negative for node2 left of node1, or positive for node2 right of node1.
 * Note: you MUST have already called handle_failure() on these nodes earlier to
 * split node2, if necessary. This assumes neither node is WITHIN the range 
 * of the other.
 */
long int node_dist(arg_node* node1, arg_node* node2){
    long int start1 = *node1->sites.begin();
    long int end1 = *node1->sites.rbegin();
    long int start2 = *node2->sites.begin();
    long int end2 = *node2->sites.rbegin();
    
    if (end1 < start2){
        return start2 - end1;
    }
    else if (start1 > end2){
        return end2 - start1;
    }
    else{
        fprintf(stderr, "can't compute node dist; was one node split?\n");
       
        exit(1);
    }
    return -1;
}

void compile_failures(arg_node* child, 
    arg_node* parent, 
    set<arg_node*>& failures,
    long int start, 
    long int end, 
    long int prop_dist,
    arg_sitemap& sites_pos,
    int num_haplotypes,
    bool exp_l,
    bool exp_r){
    bool expanding_only = false;
    
    if (start == -1){
        start = *child->sites.begin() - prop_dist;
    }
    else{
        expanding_only = true;
    }
    if (start < *child->sites.begin() - prop_dist){
        //start = *child->sites.begin() - prop_dist;
    }
    if (end == -1){
        end = *child->sites.rbegin() + prop_dist;
    }
    else{
        expanding_only = true;
    }
    if (end > *child->sites.rbegin() + prop_dist){
        //end = *child->sites.rbegin() + prop_dist;
    }
    
    if (exp_l && start > max(*child->sites.begin() - prop_dist, (long int)1)){
        start = max(*child->sites.begin() - prop_dist, (long int)1);
    }
    if (exp_r && end < *child->sites.rbegin() + prop_dist){
        end = *child->sites.rbegin() + prop_dist;
    }
    
    
    for (arg_sitemap::reverse_iterator site = sites_pos.rbegin(); 
        site != sites_pos.rend(); ++site){
        
        if (site->first + prop_dist < start){
            break;
        }
        else if (site->second != child && ranges_overlap(
            max(site->first - prop_dist, (long int)1), 
            site->first + prop_dist,
            start, 
            end)){
            
            if (DEBUG_MODE){
                fprintf(stderr, "checking site %ld\n", site->first);
                print_node_lite(site->second, num_haplotypes);
            }
            
            if (!expanding_only || site_ranges_overlap(start, end, site->second)){
                if (site->second->clade.count() > 1){
                    short comp = compare_bitsets(site->second->clade, child->clade);
                    if (comp == -1){
                        failures.insert(site->second);
                    }
                    
                }
            }
            else if (DEBUG_MODE){
                fprintf(stderr, "didn't check (%d); site %ld:\n", expanding_only, site->first);
                fprintf(stderr, "(%ld %ld) (%ld %ld) %ld\n", start, end, site->second->start, site->second->end,
                    min(end, site->second->end) - max(start, site->second->start));
                //print_node_lite(site->second, num_haplotypes);
            }
        }
    }
}

void gather_parents(arg_node* child, 
    arg_node* parent, 
    long int start, 
    long int end, 
    long int prop_dist,
    arg_sitemap& sites_pos,
    arg_clademap& sites_clade,
    arg_clademap& lv_grps,
    set<arg_node*, nodeptr_parent_sort>& nodes,
    int num_haplotypes){
    if (DEBUG_MODE){
        fprintf(stderr, "GATHER PARENTS\n");
        print_node_lite(child, num_haplotypes);
        print_node_lite(parent, num_haplotypes);
    }
    if (start == -1){
        start = child->start;
    }
    if (end == -1){
        end = child->end;
    }
    
    if (start < child->start){
        start = child->start;
    }
    if (end > child->end){
        end = child->end;
    }
    
    if (!site_ranges_overlap(start, end, child) || start < child->start || end > child->end){
        fprintf(stderr, "Invalid range detected in gather_parents(): %ld %ld\n", start, end);
        print_node_lite(child, num_haplotypes);
        exit(1);
    }

    for (set<vert_edge>::reverse_iterator edge = parent->children.rbegin();
        edge != parent->children.rend(); ++edge){
        if (edge->end < start - prop_dist){
            break;
        }
        if (site_ranges_overlap(start, end, edge->to)){
            short comp = compare_bitsets(child->clade, edge->to->clade);
            if (comp == 2){
                gather_parents(child, edge->to, start, end, prop_dist, sites_pos, 
                    sites_clade, lv_grps, nodes, num_haplotypes);
            }
        }
    }
    
    nodes.insert(parent);
    
    return;
}

// Strategy for dealing with failures should be this:
// 1. collect all failures in range
// 2. adjust indices of everything involved
// 3. add all parents/children of new clades/clades being expanded
// 4. set all 4haptest failures.

void handle_all_failures(arg_node* child,
    arg_node* parent,
    arg_sitemap& sites_pos,
    arg_clademap& sites_clade,
    arg_clademap& lv_grps,
    long int prop_dist,
    long int start,
    long int end,
    arg_node* root,
    vector<arg_node*>& newnodes,
    int num_haplotypes){
    bool has_failures = true;
    while(has_failures){
        set<arg_node*> failures;
        if (start != -1 && end != -1 && !site_ranges_overlap(start, end, child)){
            // child doesn't overlap range.
            return;
        }
        compile_failures(child, parent, failures, start, end, prop_dist, sites_pos, num_haplotypes, false, false);
        if (failures.size() == 0){
            has_failures = false;
        }
        else{
            set<arg_node*> newsplits;
            for (set<arg_node*>::iterator failure = failures.begin(); failure != failures.end();
                ++failure){
                arg_node* newsplit = handle_failure(child, *failure, sites_pos, sites_clade, lv_grps, prop_dist, false, root, newnodes, num_haplotypes);
                if (newsplit != NULL){
                    newsplits.insert(newsplit);
                }
            }
            for (set<arg_node*>::iterator newsplit = newsplits.begin(); newsplit != newsplits.end(); ++newsplit){
                failures.insert(*newsplit);
            }
            for (set<arg_node*>::iterator failure = failures.begin(); failure != failures.end(); ++failure){
                handle_failure(child, *failure, sites_pos, sites_clade, lv_grps, prop_dist, true, root, newnodes, num_haplotypes);
            }
            failures.clear();
        }
    }
}

void adjust_all_failure_indices(arg_node* child, 
    set<arg_node*>& failures,
    arg_sitemap& sites_pos,
    arg_clademap& sites_clade,
    arg_clademap& lv_grps,
    long int prop_dist,
    arg_node* root,
    vector<arg_node*>& newnodes,
    int num_haplotypes){
    set<arg_node*> newsplits;
    for (set<arg_node*>::iterator failure = failures.begin(); failure != failures.end();
        ++failure){
        arg_node* newsplit = handle_failure(child, *failure, sites_pos, sites_clade,
            lv_grps, prop_dist, false, root, newnodes, num_haplotypes);
        if (newsplit != NULL){
            newsplits.insert(newsplit);
        }
    }
    for (set<arg_node*>::iterator newsplit = newsplits.begin(); newsplit != newsplits.end(); ++newsplit){
        failures.insert(*newsplit);
    }
    
    // Indices should now be adjusted; safe to add parents/children.
}

void handle_failures(arg_node* child,
    set<arg_node*>& failures,
    arg_sitemap& sites_pos,
    arg_clademap& sites_clade,
    arg_clademap& lv_grps,
    long int prop_dist,
    arg_node* root,
    vector<arg_node*>& newnodes,
    int num_haplotypes){
    set<pair<long int, arg_node*> > dists_right;
    set<pair<long int, arg_node*> > dists_left;
    
    for (set<arg_node*>::iterator failure = failures.begin(); failure != failures.end();
        ++failure){
        
        long int nd = node_dist(child, *failure);
        if (nd < 0){
            dists_left.insert(make_pair(-nd, *failure));
        }
        else{
            dists_right.insert(make_pair(nd, *failure));
        }
    }
    
    // Now process things to the right, in order
    vector<arg_node*> right_processed;

    for (set<pair<long int, arg_node*> >::iterator dr = dists_right.begin(); dr != dists_right.end();
        ++dr){
        bool pass = true;
        for (vector<arg_node*>::iterator rp = right_processed.begin(); rp != right_processed.end();
            ++rp){
            
            short comp = compare_bitsets((*rp)->clade, dr->second->clade);
            if (comp == -1 && *(*rp)->sites.begin() < *dr->second->sites.begin()){
                pass = false;
                break;
            }
        }
        if (pass){
            right_processed.push_back(dr->second);
        }
    }
    
    // Process things to the left, in order
    vector<arg_node*> left_processed;
    for (set<pair<long int, arg_node*> >::iterator dl = dists_left.begin(); dl != dists_left.end();
        ++dl){
        bool pass = true;
        for (vector<arg_node*>::iterator lp = left_processed.begin(); lp != left_processed.end();
            ++lp){
            
            short comp = compare_bitsets((*lp)->clade, dl->second->clade);
            if (comp == -1 && *(*lp)->sites.rbegin() > *dl->second->sites.rbegin()){
                pass = false;
                break;
            }
        }
        if (pass){
            left_processed.push_back(dl->second);
        }
    }
    
    cladeset_set right_clades;
    for (vector<arg_node*>::iterator rp = right_processed.begin(); rp != right_processed.end();
        ++rp){
        if (right_clades.find((*rp)->clade) == right_clades.end()){
            handle_failure(child, *rp, sites_pos, sites_clade, lv_grps, prop_dist, true, root, newnodes, num_haplotypes);
            right_clades.insert((*rp)->clade);
        }
    }
    cladeset_set left_clades;
    for (vector<arg_node*>::iterator lp = left_processed.begin(); lp != left_processed.end();
        ++lp){
        if (left_clades.find((*lp)->clade) == left_clades.end()){
            handle_failure(child, *lp, sites_pos, sites_clade, lv_grps, prop_dist, true, root, newnodes, num_haplotypes);
            left_clades.insert((*lp)->clade);
        }
    }
    
    
    
    if (DEBUG_MODE){
        fprintf(stderr, "CHILD\n");
        print_node_lite(child, num_haplotypes);
        
        for (vector<rc_edge>::iterator to = child->edges_to.begin(); to != child->edges_to.end(); ++to){
            fprintf(stderr, "to:\n");
            print_node_lite(to->node, num_haplotypes);
            fprintf(stderr, "fd:\n");
            print_node_lite(to->final_dest, num_haplotypes);
            fprintf(stderr, "\n");
        }
    }
}

void insert_below_parent(arg_node* child, 
    arg_node* parent, 
    arg_sitemap& sites_pos,
    arg_clademap& sites_clade,
    arg_clademap& lv_grps,
    long int prop_dist,
    long int start, 
    long int end,
    arg_node* root,
    vector<arg_node*>& newnodes,
    int num_haplotypes){
    if (child == parent){
        return;
    }
    else if (child->clade == parent->clade){
        fprintf(stderr, "insert_below_parent(): child equals parent\n");
        print_node_lite(child, num_haplotypes);
        fprintf(stderr, "\n");
        print_node_lite(parent, num_haplotypes);
        exit(1);
    }
    
    if (DEBUG_MODE){
        fprintf(stderr, "ibp\n");
        print_node_lite(child, num_haplotypes);
        print_node_lite(parent, num_haplotypes);
    }
    
    set<arg_node*> failures;
    compile_failures(child, parent, failures, start, end, prop_dist, sites_pos, num_haplotypes, false, false);
    
    adjust_all_failure_indices(child, failures, sites_pos, sites_clade, lv_grps, 
        prop_dist, root, newnodes, num_haplotypes);
    
    set<arg_node*, nodeptr_parent_sort> parents;
    if (DEBUG_MODE){
        fprintf(stderr, "gp1 %ld %ld\n", start, end);
    }
    
    bool add_parents = true;
    bool leftward = false;
    
    if (start != -1 && start < child->start){
        // We were expanding range but hit a failure that actually contracted the range.
        start = child->start;
    }
    else if (end != -1 && end > child->end){
        // We were expanding range but hit a failure that actually contracted the range.
        end = child->end;
    }
    
    if (end < start){
        add_parents = false;
    }
    
    
    if (add_parents){
        gather_parents(child, parent, start, end, prop_dist, sites_pos, sites_clade, lv_grps, parents, num_haplotypes);
    
        if (DEBUG_MODE){
            fprintf(stderr, "parents:\n");
            for (set<arg_node*, nodeptr_parent_sort>::iterator parent = parents.begin(); parent != parents.end(); ++parent){
                print_node_downward_onelevel(*parent, 0, num_haplotypes);
            }
            fprintf(stderr, "\n");
        }
    
        // Now add parents as necessary. Smallest should be first.
        set<arg_node*, nodeptr_child_sort> children;
        for (set<arg_node*, nodeptr_parent_sort>::iterator node = parents.begin(); node != parents.end();
            ++node){
            set_parents(child, *node, false, prop_dist, sites_pos, sites_clade, lv_grps, children, num_haplotypes);
        }
    
        // Note: when we add a discarded child to the new node, there's no way we'll be knocking off anything
        // below the new node (as a new discarded child), since the node is new and we're adding
        // in the order of highest -> lowest level.
        set<arg_node*, nodeptr_child_sort> dummy;
        for (set<arg_node*, nodeptr_child_sort>::iterator discarded = children.begin(); discarded != children.end();
            ++discarded){
            set_parents(*discarded, child, false, prop_dist, sites_pos, sites_clade, lv_grps, dummy, num_haplotypes);
        }   
    } 
    
    // Finally, handle failures.
    handle_failures(child, failures, sites_pos, sites_clade, lv_grps, prop_dist, root, newnodes, num_haplotypes);
    
    return;
}

void elim_impossible_lv_grps(arg_clademap& lv_grps, arg_node* node, int num_haplotypes){
    return;
    
    for (arg_clademap::iterator lv = lv_grps.begin(); lv != lv_grps.end(); ){
        bool erase = false;
        short comp = compare_bitsets(lv->first, node->clade);
        if (comp == -1){
            if (node_ranges_overlap(lv->second, node)){
                if (DEBUG_MODE){
                    fprintf(stderr, "REMOVING:\n");
                    print_node_lite(lv->second, num_haplotypes);
                    
                    fprintf(stderr, "node:\n");
                    print_node_lite(node, num_haplotypes);
                }
                // Remove this.
                erase = true;
            }
        }
        else if (comp == 3 && (node->recomb_left != -1 || node->recomb_right != -1)){
            // Look for identical ones that intersect with sites that this site
            // says must have recombined from it.
            if (node->recomb_left != -1 && node->recomb_left >= lv->second->start && 
                node->recomb_left <= lv->second->end){
                erase = true;
            }
            else if (node->recomb_right != -1 && node->recomb_right >= lv->second->start && 
                node->recomb_right <= lv->second->end){
                erase = true;
            }
        }
        
        if (erase){
            del_leaving_node(lv->second, num_haplotypes);
            delete lv->second;
            lv_grps.erase(lv++);
        }
        else{
            ++lv;
        }
    }
}

/**
 * Returns false if the new clade would break any preexisting contiguous nodes.
 */
bool safe_to_create_node(cladeset& clade,
    arg_sitemap& sites_pos,
    arg_clademap& sites_clade,
    long int pos,
    long int prop_dist,
    int num_haplotypes){
    if (clade.count() == num_haplotypes){
        return false;
    }
    if (clade.count() == 1){
        // Singletons can't mess anything up.
        return true;
    }
    
    // Check to make sure it won't fail the four haplotype test with anything tied to
    // the same site OR be identical to something at the same site.
    if (sites_pos.count(pos) > 0){
        pair<arg_sitemap::iterator, arg_sitemap::iterator> er = sites_pos.equal_range(pos);
        for (arg_sitemap::iterator n = er.first; n != er.second; ++n){
            short comp = compare_bitsets(clade, n->second->clade);
            if (comp == -1){
                return false;
            }
            else if (comp == 3){
                return false;
            }
            // also check (solved) recombinations to make sure it's okay.
            
            for (vector<rc_edge>::iterator efs = n->second->edges_from_solved.begin();
                efs != n->second->edges_from_solved.end(); ++efs){
                if (efs->final_dest != NULL && efs->final_dest->clade == clade){
                    return false;
                }
            }
            for (vector<rc_edge>::iterator ets = n->second->edges_to_solved.begin();
                ets != n->second->edges_to_solved.end(); ++ets){
                if (ets->final_dest != NULL && ets->final_dest->clade == clade){
                    return false;
                }
            }
        }
    }
    
    // Ensure no identical nodes are marked as invalid by recombination at the site
    // at which this node will exist.
    pair<arg_clademap::iterator, arg_clademap::iterator> er = sites_clade.equal_range(clade);
    for (arg_clademap::iterator n = er.first; n != er.second; ++n){
        if (n->second->recomb_right == pos || n->second->recomb_left == pos){
            return false;
        }
    }
    
    //return true;
    
    arg_sitemap::iterator it = sites_pos.equal_range(pos).first;
    arg_sitemap::reverse_iterator leftward(it);
    ++leftward;
    for (; leftward != sites_pos.rend(); ++leftward){
        if (leftward->first < pos - prop_dist){
            break;   
        }
        
        if (leftward->second->clade.count() > 1 &&
            leftward->second->start <= pos && leftward->second->end >= pos){
            
            short comp = compare_bitsets(leftward->second->clade, clade);
            if (comp == -1 && pos >= *(leftward->second->sites.begin()) && pos <= *(leftward->second->sites.end())){
                return false;
            }
        }
    }
    arg_sitemap::iterator rightward = sites_pos.equal_range(pos).second;
    ++rightward;
    for (; rightward != sites_pos.end(); ++rightward){
        if (rightward->first > pos + prop_dist){
            break;
        }

        if (rightward->second->clade.count() > 1 && 
            rightward->second->start <= pos && rightward->second->end >= pos){
            
            short comp = compare_bitsets(rightward->second->clade, clade);
            //if (comp == -1){
            if (comp == -1 && pos >= *(rightward->second->sites.begin()) && pos <= *(rightward->second->sites.end())){
                return false;
            }
        }
    }
    
    return true;
}

arg_node* create_node(arg_node* node, 
    arg_sitemap& sites_pos,
    arg_clademap& sites_clade,
    arg_clademap& lv_grps,
    long int prop_dist,
    arg_node* root,
    vector<arg_node*>& newnodes,
    int num_haplotypes,
    set<arg_node*>& nodes_external){
    node->newnode = true;
    
    node->fixed_left = false;
    node->fixed_right = false;
    /*
    if (node->closest_mut_l() != -1 && node->start < node->closest_mut_l() - prop_dist){
        node->start = node->closest_mut_l() - prop_dist;
    }
    if (node->closest_mut_r() != -1 && node->end > node->closest_mut_r() + prop_dist){
        node->end = node->closest_mut_r() + prop_dist;
    }
    */
    
    if (node->end < node->start){
        fprintf(stderr, "create node: invalid range\n");
        print_node_lite(node, num_haplotypes);
        exit(1);
    }
    
    if (node->clade.count() == num_haplotypes){
        
        if (DEBUG_MODE){
            fprintf(stderr, "CREATING ROOT?\n");
        }
        for (set<long int>::iterator mut = node->mutations.begin(); mut != node->mutations.end(); ++mut){
            root->mutations.insert(*mut);
        }
    
        return root;
    }
    if (node->clade.count() == 1 && node->mutations.size() == 0 && !node->is_leaving_grp){
        // Don't bother creating -- it won't affect the range or branch lengths
        // of any other nodes (since no mutations support it)
        return NULL;
    }
    
    if (node->start > *node->sites.begin()){
        fprintf(stderr, "create node: illegal range 1\n");
        print_node_lite(node, num_haplotypes);
        exit(1);
    }
    if (node->end < *node->sites.rbegin()){
        fprintf(stderr, "create node: illegal range 2\n");
        print_node_lite(node, num_haplotypes);
        exit(1);
    }
    
    if (node->start < *node->sites.begin() - prop_dist){
        node->start = *node->sites.begin() - prop_dist;
    }
    if (node->end > *node->sites.rbegin() + prop_dist){
        node->end = *node->sites.rbegin() + prop_dist;
    }
    
    if (DEBUG_MODE){
        fprintf(stderr, "create node\n");
        print_node_lite(node, num_haplotypes);
    }
    
    if (node->sites.find(-1) != node->sites.end()){
        fprintf(stderr, "???0\n");
        print_node_lite(node, num_haplotypes);
        exit(1);
    }
    
    // If it fails the four haplotype test failure with another node at the same
    // site, we can't create it.
    for (set<long int>::iterator site = node->sites.begin(); site != node->sites.end();
        ++site){
        if (sites_pos.count(*site) > 0){
            pair<arg_sitemap::iterator, arg_sitemap::iterator> er = sites_pos.equal_range(*site);
            for (arg_sitemap::iterator prevnode = er.first; prevnode != er.second;
                ++prevnode){
                short comp = compare_bitsets(prevnode->second->clade, node->clade);
                if (comp == -1){
                    //delete node;
                    return NULL;
                }
            }
        }
    }
    
    bool range_set = false;
    
    bool merge_nodes = false;
    
    arg_node* nodetomerge = NULL;
    
    set<arg_node*> nodestomerge;
    // First check to see if we can merge with a preexisting node.
    if (sites_clade.count(node->clade) > 0){
        pair<arg_clademap::iterator, arg_clademap::iterator> er = sites_clade.equal_range(node->clade);
        for (arg_clademap::iterator prevnode = er.first; prevnode != er.second; ++prevnode){
            
            if (prevnode->second != node){
                
                bool pass = true;
                
                long int min1 = *(prevnode->second->sites.begin()) - prop_dist;
                long int max1 = *(prevnode->second->sites.rbegin()) + prop_dist;
                long int min2 = *(node->sites.begin()) - prop_dist;
                long int max2 = *(node->sites.rbegin()) + prop_dist;
                
                if (pass && min(max1, max2) - max(min1, min2) <= 0 && !node_ranges_overlap(prevnode->second, node)){
                    pass = false;
                    if (node_ranges_overlap(prevnode->second, node)){
                        pass = true;
                    }
                }
                
                set<arg_node*> parents_tmp;
                set<arg_node*> children_tmp;
                
                bool prevnode_left = false;
                bool prevnode_right = false;
                
                if (pass){

                    // Determine whether the new node is to the right or left.
                    if (*(node->sites.begin()) > *prevnode->second->sites.begin()){
                        if (prevnode->second->recomb_right != -1 && prevnode->second->recomb_right <= 
                            *(node->sites.rbegin())){
                            // Recombination between prev and current; can't merge
                            pass = false;
                            
                            if (node->sites.find(prevnode->second->recomb_right) != node->sites.end()){
                                fprintf(stderr, "create_node - node owns prevnode's rightward recomb site\n");
                                print_node_lite(prevnode->second, num_haplotypes);
                                print_node_lite(node, num_haplotypes);
                                exit(1);
                            }
                            
                            if (prevnode->second->end + 1 > node->start){
                                node->start = min(prevnode->second->end + 1, *node->sites.begin());
                                //node->start = prevnode->second->end + 1;
                            }
                            
                            if (*node->sites.begin() <= prevnode->second->recomb_right){
                                fprintf(stderr, "create_node - node's first site is <= prevnode's rightward recomb site\n");
                                print_node_lite(prevnode->second, num_haplotypes);
                                print_node_lite(node, num_haplotypes);
                                exit(1);
                            }
                        }
                        else{
                            // This ensures it's within propagation distance.
                            if (prevnode->second->recomb_left != -1 && node->start < prevnode->second->start){
                                node->start = min(prevnode->second->start, *node->sites.begin());
                                //node->start = prevnode->second->start;
                            }
                        
                            if (prevnode->second->end < *(node->sites.begin())){
                                // Check sites between old & new nodes.
                                arg_sitemap::iterator ps = sites_pos.find(*(node->sites.begin()));
                                arg_sitemap::reverse_iterator prevsite2(ps);
                                
                                while(prevsite2->first > prevnode->second->end){
                                    
                                    short comp = compare_bitsets(prevsite2->second->clade, node->clade);
                                    if (comp == -1){
                                        pass = false;
                                        break;
                                    }
                                        
                                    if (!pass){
                                        break;
                                    }
                                    ++prevsite2;
                                }
                                
                            }
                        
                        }
                    }
                    if (*(prevnode->second->sites.rbegin()) > *(node->sites.rbegin())){
                        if (prevnode->second->recomb_left != -1 && prevnode->second->recomb_left >=
                            *(node->sites.begin())){
                            pass = false;
                            if (node->sites.find(prevnode->second->recomb_left) != node->sites.end()){
                                fprintf(stderr, "create_node - node owns prevnode's leftward recomb site\n");
                                print_node_lite(prevnode->second, num_haplotypes);
                                print_node_lite(node, num_haplotypes);
                                fprintf(stderr, "other nodes at site:\n");
                                for (arg_sitemap::iterator n = sites_pos.equal_range(prevnode->second->recomb_left).first;
                                    n != sites_pos.equal_range(prevnode->second->recomb_left).second; ++n){ 
                                    print_node_lite(n->second, num_haplotypes);
                                }
                                exit(1);
                            }
                            
                            if (prevnode->second->start - 1 < node->end){
                                node->end = max(prevnode->second->start - 1, *node->sites.rbegin());
                                //node->end = prevnode->second->start-1;
                            }
                            
                            // If the node has sites that surround the other node, it's tough
                            // to know what to do.
                            if (*node->sites.rbegin() >= prevnode->second->recomb_left){
                                fprintf(stderr, "create_node - node's rightmost site is >= prevnode's leftward recomb site\n");
                                print_node_lite(prevnode->second, num_haplotypes);
                                print_node_lite(node, num_haplotypes);
                                exit(1);
                            }
                            
                        }
                        else{
                            // This ensures it's within propagation distance.
                            if (prevnode->second->recomb_right != -1 && node->end > prevnode->second->end){
                                node->end = max(prevnode->second->end, *node->sites.rbegin());
                                //node->end = prevnode->second->end;
                            }
                        
                            if (prevnode->second->start > *(node->sites.rbegin())){
                                arg_sitemap::iterator prevsite2 = sites_pos.find(*(node->sites.rbegin()));
                                
                                while(prevsite2->first < prevnode->second->start){
                                    
                                    short comp = compare_bitsets(prevsite2->second->clade, node->clade);
                                    
                                    if (comp == -1){
                                        pass = false;
                                        break;
                                    }
                                    
                                    if (!pass){
                                        break;
                                    }
                                    ++prevsite2;
                                }
                                
                            }
                        }
                    }
                }
                if (pass){
                    merge_nodes = true;
                    nodestomerge.insert(prevnode->second);
                    nodetomerge = prevnode->second;
                }
            }
        }
    }
    
    if (merge_nodes){
        // Merge these nodes.
        set<arg_node*> mergenodes_left;
        set<arg_node*> mergenodes_right;
        set<arg_node*> mergenodes_span;
        
        bool done = false;
        
        for (set<arg_node*>::iterator n = nodestomerge.begin(); n != nodestomerge.end(); ++n){
            
            // Look whether it's left, right, spanning, or within the new node's range.
            long int start1 = *(*n)->sites.begin();
            long int end1 = *(*n)->sites.rbegin();
            long int start2 = *node->sites.begin();
            long int end2 = *node->sites.rbegin();
            
            if (start1 > start2 && end1 < end2){
                // Node1 within Node2 range
                mergenodes_span.insert(*n);
            }
            else if (start2 >= start1 && end2 <= end1){
                // Node2 within Node1 range
                // Add sites to it and be done.
          
                for (set<long int>::iterator s = node->sites.begin(); s != node->sites.end(); ++s){
                    if (node->mutations.find(*s) != node->mutations.end()){
                        (*n)->mutations.insert(*s);
                    }
                    if (node->sites_lvgrp.find(*s) != node->sites_lvgrp.end()){
                        (*n)->sites_lvgrp.insert(*s);
                    }
                    (*n)->sites.insert(*s);
                    pair<arg_sitemap::iterator, arg_sitemap::iterator> er = sites_pos.equal_range(*s);
                    for (arg_sitemap::iterator it = er.first;
                        it != er.second; ){
                        if (it->second == node){
                            sites_pos.erase(it++);
                            break;
                        }
                        else{
                            ++it;
                        }
                    }
                    er = sites_pos.equal_range(*s);
                    bool found = false;
                    for (arg_sitemap::iterator it = er.first;
                        it != er.second; ++it){
                        if (it->second == *n){
                            found = true;
                            break;
                        }
                    }
                    if (!found){
                        sites_pos.insert(make_pair(*s, *n));
                    }
                }
                if ((*n)->start > *(*n)->sites.begin() || (*n)->end < *(*n)->sites.rbegin()){
                    fprintf(stderr, "merge_nodes - impossible node range\n");
                    print_node_lite(*n, num_haplotypes);
                    exit(1);
                }
                done = true;
                return *n;
                break;
            }
            else if ((start1 < start2 && end1 > start2 && end1 < end2) || end1 < start2){
                // Node1 hangs off to left of Node2
                mergenodes_left.insert(*n);
            }
            else if ((start2 < start1 && end2 > start1 && end2 < end1) || end2 < start1){
                // Node 1 hangs off to right of Node2
                mergenodes_right.insert(*n);
            }
        }
        
        /*
        if (mergenodes_left.size() > 1 || mergenodes_right.size() > 1 || 
            mergenodes_left.size() + mergenodes_right.size() > 2 ||
            mergenodes_span.size() >= 2){
      
            exit(1);
        }
        */
        if ((mergenodes_span.size() == 1 && 
            (mergenodes_left.size() == 1 || mergenodes_right.size() == 1)) ||
            (mergenodes_left.size() == 1 && mergenodes_right.size() == 1)){
            // Putting two existing nodes together.
   
            arg_node* left_node = NULL;
            arg_node* right_node = NULL;
            if (mergenodes_left.size() > 0){
                left_node = *mergenodes_left.begin();
            }
            else{
                left_node = *mergenodes_span.begin();
            }
            if (mergenodes_right.size() > 0){
                right_node = *mergenodes_right.begin();
            }
            else{
                right_node = *mergenodes_span.begin();
            }
        
            set<long int> sites_between;
            for (set<long int>::iterator site = node->sites.begin(); site != node->sites.end(); ++site){
                if (*site > *left_node->sites.rbegin() && *site < *right_node->sites.begin()){
                    sites_between.insert(*site);
                }
            }
            bool exp_r = false;
            
            if (*sites_between.rbegin() > right_node->start-1){
                // Expand right node leftward.
                exp_r = true;    
            }
            else{
                // Expand left node rightward.
                exp_r = false;
            }
            
            for (set<long int>::iterator site = node->sites.begin(); site != node->sites.end(); ++site){
                pair<arg_sitemap::iterator, arg_sitemap::iterator> er = sites_pos.equal_range(*site);
                for (arg_sitemap::iterator it = er.first;
                    it != er.second; ){
                    if (it->second == node){
                        sites_pos.erase(it++);
                        break;
                    }
                    else{
                        ++it;
                    }
                }
                if (exp_r){
                    right_node->sites.insert(*site);
                    bool found = false;
                    er = sites_pos.equal_range(*site);
                    for (arg_sitemap::iterator it = er.first;
                        it != er.second; ++it){
                        if (it->second == right_node){
                            found = true;
                            break;
                        }
                    }
                    if (!found){
                        sites_pos.insert(make_pair(*site, right_node));
                    }
                }
                else{
                    left_node->sites.insert(*site);
                    bool found = false;
                    er = sites_pos.equal_range(*site);
                    for (arg_sitemap::iterator it = er.first;
                        it != er.second; ++it){
                        if (it->second == left_node){
                            found = true;
                            break;
                        }
                    }
                    if (!found){
                        sites_pos.insert(make_pair(*site, left_node));
                    }
                }
            }
            for (set<long int>::iterator mut = node->mutations.begin(); mut != node->mutations.end(); ++mut){
                if (exp_r){
                    right_node->mutations.insert(*mut);
                }
                else{
                    left_node->mutations.insert(*mut);
                }
            }
            
            long int origstart1 = left_node->start;
            long int origend1 = left_node->end;
            long int origstart2 = right_node->start;
            long int origend2 = right_node->end;
            
            cladeset cl = left_node->clade;
            long int leftsitekey = *left_node->sites.begin();
            long int rightsitekey = *right_node->sites.rbegin();
            
            if (exp_r){
                if (left_node->end+1 < right_node->start){
                    if (left_node->end+1 < *right_node->sites.begin()-prop_dist){
                        long int origstart = right_node->start;
                        right_node->start = *right_node->sites.begin()-prop_dist;
                        expand_range(right_node, right_node->start-origstart, sites_pos, sites_clade,
                            lv_grps, prop_dist, root, newnodes, num_haplotypes, true);
                        if (right_node->start > *right_node->sites.begin()-prop_dist){
                            origstart = right_node->start;
                            right_node->start = *right_node->sites.begin()-prop_dist;
                            expand_range(right_node, right_node->start-origstart, sites_pos, 
                                sites_clade,
                                lv_grps, prop_dist, root, newnodes, num_haplotypes, true);
                            if (right_node->start > *right_node->sites.begin()-prop_dist){
                                fprintf(stderr, "didn't expand enough 1A\n");
                                print_node_lite(left_node, num_haplotypes);
                                print_node_lite(right_node, num_haplotypes);
                                exit(1);
                            }
                        }
                        // Bring up left node too.
                        long int origend = left_node->end;
                        left_node->end = right_node->start-1;
                        expand_range(left_node, left_node->end-origend, sites_pos, sites_clade,
                            lv_grps, prop_dist, root, newnodes, num_haplotypes, true);
                        if (left_node->end < right_node->start-1){
                            origend = left_node->end;
                            left_node->end = right_node->start-1;
                            expand_range(left_node, left_node->end-origend, sites_pos, sites_clade,
                                lv_grps, prop_dist, root, newnodes, num_haplotypes, true);
                            if (left_node->end < right_node->start-1){
                                fprintf(stderr, "didn't expand enough 1B\n");
                                print_node_lite(left_node, num_haplotypes);
                                print_node_lite(right_node, num_haplotypes);
                                exit(1);
                            }
                        }
                    }
                    else{
                        long int origstart = right_node->start;
                        right_node->start = left_node->end+1;
                        expand_range(right_node, right_node->start-origstart, sites_pos, sites_clade,
                            lv_grps, prop_dist, root, newnodes, num_haplotypes, true);
                        if (right_node->start > left_node->end+1){
                            origstart = right_node->start;
                            right_node->start = left_node->end+1;
                            expand_range(right_node, right_node->start-origstart, sites_pos, sites_clade,
                                lv_grps, prop_dist, root, newnodes, num_haplotypes, true);
                            if (right_node->start > left_node->end+1){
                                fprintf(stderr, "didn't expand enough 1C\n");
                                print_node_lite(left_node, num_haplotypes);
                                print_node_lite(right_node, num_haplotypes);
                                exit(1);
                            }
                        }
                    }
                }
                
            }
            else{
                if (right_node->start-1 > left_node->end){
                    if (right_node->start-1 > *left_node->sites.rbegin() + prop_dist){
                        long int origend = left_node->end;
                        left_node->end = *left_node->sites.rbegin() + prop_dist;
                        expand_range(left_node, left_node->end-origend, sites_pos, sites_clade,
                            lv_grps, prop_dist, root, newnodes, num_haplotypes, true);
                        if (left_node->end < *left_node->sites.rbegin() + prop_dist){
                            origend = left_node->end;
                            left_node->end = *left_node->sites.rbegin() + prop_dist;
                            expand_range(left_node, left_node->end-origend, sites_pos, sites_clade,
                                lv_grps, prop_dist, root, newnodes, num_haplotypes, true);
                            if (left_node->end < *left_node->sites.rbegin() + prop_dist){
                                fprintf(stderr, "didn't expand enough 2A\n");
                                print_node_lite(left_node, num_haplotypes);
                                print_node_lite(right_node, num_haplotypes);
                                exit(1);
                            }
                        }
                        // Need to bring down right node too.
                        long int origstart = right_node->start;
                        right_node->start = left_node->end+1;
                        expand_range(right_node, right_node->start-origstart, sites_pos, sites_clade,
                            lv_grps, prop_dist, root, newnodes, num_haplotypes, true);
                        if (right_node->start > left_node->end+1){
                            origstart = right_node->start;
                            right_node->start = left_node->end+1;
                            expand_range(right_node, right_node->start-origstart, sites_pos, sites_clade,
                                lv_grps, prop_dist, root, newnodes, num_haplotypes, true);
                            if (right_node->start > left_node->end+1){
                                fprintf(stderr, "didn't expand enough 2B\n");
                                print_node_lite(left_node, num_haplotypes);
                                print_node_lite(right_node, num_haplotypes);
                                exit(1);
                            }
                        }
                    }
                    else{
                        long int origend = left_node->end;
                        left_node->end = right_node->start-1;
                        expand_range(left_node, left_node->end-origend, sites_pos, sites_clade,
                            lv_grps, prop_dist, root, newnodes, num_haplotypes, true);
                        if (left_node->end < right_node->start-1){
                            origend = left_node->end;
                            left_node->end = right_node->start-1;
                            expand_range(left_node, left_node->end-origend, sites_pos, sites_clade,
                                lv_grps, prop_dist, root, newnodes, num_haplotypes, true);
                            if (left_node->end < right_node->start-1){
                                fprintf(stderr, "didn't expand enough 2C\n");
                                print_node_lite(left_node, num_haplotypes);
                                print_node_lite(right_node, num_haplotypes);
                                exit(1);
                            }
                        }
                    }
                }
            }
            
            // In most cases we can just join the two nodes here. It's possible, though, that
            // when calling expand_range() above, we hit a four haplotype test failure,
            // which was already solved using preexisting (solved) recombination edges.
            // This would result in expand_range_fixed() being called, which has the potential
            // to merge nodes. In that case, then, we need to re-look up both nodes and
            // see if there's still anything to merge.
            left_node = NULL;
            right_node = NULL;
            int leftnodecount = 0;
            int rightnodecount = 0;
            pair<arg_clademap::iterator, arg_clademap::iterator> er = sites_clade.equal_range(cl);
            for (arg_clademap::iterator it = er.first; it != er.second; ++it){
                if (it->second->sites.find(leftsitekey) != it->second->sites.end()){
                    leftnodecount++;
                    left_node = it->second;
                }
                if (it->second->sites.find(rightsitekey) != it->second->sites.end()){
                    rightnodecount++;
                    right_node = it->second;
                }
            }
            if (left_node != NULL && right_node != NULL && left_node != right_node &&
                left_node->end == right_node->start-1){
            
                join_nodes(left_node, right_node, sites_pos, sites_clade, lv_grps,
                    prop_dist, true, root, newnodes, num_haplotypes);
                for (arg_clademap::iterator it = sites_clade.equal_range(right_node->clade).first; 
                    it != sites_clade.equal_range(right_node->clade).second;){
                    if (it->second == right_node){
                        sites_clade.erase(it++);
                        break;
                    }
                    else{
                        ++it;
                    }
                }
                for (set<long int>::iterator site = right_node->sites.begin(); site != right_node->sites.end();
                    ++site){
                    pair<arg_sitemap::iterator, arg_sitemap::iterator> er = sites_pos.equal_range(*site);
                    for (arg_sitemap::iterator it = er.first;
                        it != er.second;){
                        if (it->second == right_node){
                            sites_pos.erase(it++);
                            break;
                        }
                        else{
                            ++it;
                        }
                    }
                }
                for (set<long int>::iterator site = left_node->sites.begin(); site != left_node->sites.end();
                    ++site){
                    pair<arg_sitemap::iterator, arg_sitemap::iterator> er = sites_pos.equal_range(*site);
                    for (arg_sitemap::iterator it = er.first;
                        it != er.second;){
                        if (it->second == right_node){
                            sites_pos.erase(it++);
                            break;
                        }
                        else{
                            ++it;
                        }
                    }
                }
                if (nodes_external.find(right_node) != nodes_external.end()){
                    nodes_external.erase(nodes_external.find(right_node));
                }
                delete right_node;
                // This might be unnecessary
                //update_recinds_left(left_node);
                //update_recinds_right(left_node);
                return left_node;
            }
            else{
                // Already merged?
                if (left_node != NULL){
                    return left_node;
                }
                else{
                    return right_node;
                }
            }
        }
        else{
          
            // Combining one existing node with newnode.
            if (mergenodes_left.size() > 0){
           
                arg_node* left_node = *mergenodes_left.begin();
                for (set<long int>::iterator site = node->sites.begin(); site != node->sites.end(); ++site){
                    pair<arg_sitemap::iterator, arg_sitemap::iterator> er = sites_pos.equal_range(*site);
                    for (arg_sitemap::iterator it = er.first;
                        it != er.second; ){
                        if (it->second == node){
                            
                            sites_pos.erase(it++);
                            break;
                        }
                        else{
                            ++it;
                        }
                    }
                    left_node->sites.insert(*site);
                    bool found = false;
                    er = sites_pos.equal_range(*site);
                    for (arg_sitemap::iterator it = er.first;
                        it != er.second; ++it){
                        if (it->second == left_node){
                            found = true;
                            break;
                        }
                    }
                    if (!found){
                        sites_pos.insert(make_pair(*site, left_node));
                    }
                }
                for (set<long int>::iterator mut = node->mutations.begin(); mut != node->mutations.end(); ++mut){
                    left_node->mutations.insert(*mut);
                }
                update_recinds_left(left_node);
                
                long int origend = left_node->end;
                if (node->end > left_node->end){
                    if (left_node->recomb_right == -1){
                        left_node->end = node->end;
                    }
                    else{
                        left_node->end = min(left_node->recomb_right-1, node->end);
                    }
          
                    expand_range(left_node, left_node->end-origend, sites_pos, sites_clade,
                        lv_grps, prop_dist, root, newnodes, num_haplotypes, true);
                    
                    // This block is probably unnecessary
                    if (left_node->end < *node->sites.begin()){
                        origend = left_node->end;
                        left_node->end = node->end;
                        expand_range(left_node, left_node->end-origend, sites_pos, sites_clade,
                            lv_grps, prop_dist, root, newnodes, num_haplotypes, true);
                        if (left_node->end < *node->sites.begin()){
                            fprintf(stderr, "nodes not merge-able after range expansion(1)\n");
                            print_node_lite(left_node, num_haplotypes);
                            print_node_lite(node, num_haplotypes);
                            exit(1);
                        }
                    }
                }
                
                
                //insert_below_parent(left_node, root, sites_pos, sites_clade, lv_grps,
                //    prop_dist, -1, -1, root, newnodes, num_haplotypes);
                if (left_node->end < *left_node->sites.rbegin() || left_node->start > *left_node->sites.begin()){
                    fprintf(stderr, "????1\n");
                    exit(1);
                }
                
                
                return left_node;
            }
            else if (mergenodes_right.size() > 0){
    
                arg_node* right_node = *mergenodes_right.begin();
                for (set<long int>::iterator site = node->sites.begin(); site != node->sites.end(); ++site){
                    pair<arg_sitemap::iterator, arg_sitemap::iterator> er = sites_pos.equal_range(*site);
                    for (arg_sitemap::iterator it = er.first;
                        it != er.second; ){
                        if (it->second == node){
                            sites_pos.erase(it++);
                            break;
                        }
                        else{
                            ++it;
                        }
                    }
                    right_node->sites.insert(*site);
                    bool found = false;
                    er = sites_pos.equal_range(*site);
                    for (arg_sitemap::iterator it = er.first;
                        it != er.second; ++it){
                        if (it->second == right_node){
                            found = true;
                            break;
                        }
                    }
                    if (!found){
                        sites_pos.insert(make_pair(*site, right_node));
                    }
                }
                for (set<long int>::iterator mut = node->mutations.begin(); mut != node->mutations.end(); ++mut){
                    right_node->mutations.insert(*mut);
                }
                update_recinds_right(right_node);
                long int origstart = right_node->start;
                if (node->start < right_node->start){
                    if (right_node->recomb_left == -1){
                        right_node->start = node->start;
                    }
                    else{
                        right_node->start = min(right_node->recomb_left+1, node->start);
                    }
                    expand_range(right_node, right_node->start-origstart, sites_pos, sites_clade,
                        lv_grps, prop_dist, root, newnodes, num_haplotypes, true);
                    
                    // This block is probably unnecessary
                    if (right_node->start > *node->sites.rbegin()){
                        origstart = right_node->start;
                        right_node->start = node->start;
                        expand_range(right_node, right_node->start-origstart, sites_pos, sites_clade,
                            lv_grps, prop_dist, root, newnodes, num_haplotypes, true);
                        if (right_node->start > *node->sites.rbegin()){
                            fprintf(stderr, "nodes not merge-able after range expansion(2)\n");
                            print_node_lite(node, num_haplotypes);
                            print_node_lite(right_node, num_haplotypes);
                            exit(1);
                        }
                    }
                }
                
                if (right_node->start > *right_node->sites.begin() || right_node->end < *right_node->sites.rbegin()){
                    fprintf(stderr, "?????2\n");
                    exit(1);
                }
                return right_node;
            }
            else{
                fprintf(stderr, "no merge-able nodes?\n");
                fprintf(stderr, "%ld %ld %ld\n", mergenodes_left.size(), mergenodes_right.size(), mergenodes_span.size());
                fprintf(stderr, "%ld\n", nodestomerge.size());
                print_node_lite(node, num_haplotypes);
                for (set<arg_node*>::iterator n = nodestomerge.begin(); n != nodestomerge.end(); ++n){
                    print_node_lite(*n, num_haplotypes);
                }
                exit(1);
            }
        }
    }
    else{                   
    
        // Store in maps for easy lookup.
        set<long int> sites_orig = node->sites;
        for (set<long int>::iterator site = node->sites.begin(); site != node->sites.end();
            ++site){
            bool found = false;
            pair<arg_sitemap::iterator, arg_sitemap::iterator> er = sites_pos.equal_range(*site);
            for (arg_sitemap::iterator it = er.first;
                it != er.second; ++it){
                if (it->second == node){
                    found = true;
                    break;
                }   
            }
            if (!found){
                sites_pos.insert(make_pair(*site, node));
            }
        }
        bool clfound = false;
        pair<arg_clademap::iterator, arg_clademap::iterator> cl_er = sites_clade.equal_range(node->clade);
        for (arg_clademap::iterator it = cl_er.first;
            it != cl_er.second; ++it){
            if (it->second == node){
                clfound = true;
                break;
            }
        }
        if (!clfound){
            sites_clade.insert(make_pair(node->clade, node));
        }
        
        if (!range_set){
            // Insert into the ARG.

            insert_below_parent(node, root, sites_pos, sites_clade, lv_grps, 
                prop_dist, -1, -1, root, newnodes, num_haplotypes);
        
            // If we survived this far, we should delete any leaving groups that overlap with & 
            // fail the four haplotype test failure with this node.
            elim_impossible_lv_grps(lv_grps, node, num_haplotypes);
        }
        
        node->newnode = false;
    
        // Allow any new stuff to get stored with the node.
        for (set<long int>::iterator site = node->sites.begin(); site != node->sites.end();
            ++site){
            
            if (sites_orig.find(*site) == sites_orig.end()){
                bool found = false;
                pair<arg_sitemap::iterator, arg_sitemap::iterator> er = sites_pos.equal_range(*site);
                for (arg_sitemap::iterator it = er.first;
                    it != er.second; ++it){
                    if (it->second == node){
                        found = true;
                        break;
                    }
                }
                if (!found){
                    sites_pos.insert(make_pair(*site, node));
                }
            }
        }
    
        if (node->end < node->start){
            fprintf(stderr, "create_node(): node end coord is less than start\n");
            print_node_lite(node, num_haplotypes);
            exit(1);
        }
    
        if (DEBUG_MODE){
            fprintf(stderr, "node created\n");
            print_node_upward_onelevel(node, 0, num_haplotypes);
        }
    
        return node;
    }
    

}

/**
 * This function finds where new edges connecting child and parent can fit in a currently existing
 * context.
 */
bool sp_aux(arg_node* child, 
    arg_node* parent, 
    vector<vert_edge>& edges_in,
    vector<vert_edge>& edges_new, 
    bool downward,
    bool make_changes,
    long int prop_dist,
    set<arg_node*>& nodes_blocking,
    arg_sitemap& sites_pos,
    arg_clademap& sites_clade,
    arg_clademap& lv_grps,
    set<arg_node*, nodeptr_child_sort>& discarded,
    int num_haplotypes){
    map<arg_node*, vector<vert_edge> > edges_discarded;
    
    vector<vert_edge> edges_old_adjusted;
    
    bool edges_already_found = false;
    if (edges_in.size() == 0){
        arg_node* ptr;
        if (downward){
            ptr = child;
        }
        else{
            ptr = parent;
        }
        vert_edge ei(max(child->start, parent->start), min(child->end, parent->end), ptr);
        edges_in.push_back(ei);
    }
    else{
        edges_already_found = true;
    }
    
    
    long int cur_context_start = max(child->start, parent->start);
    long int cur_context_end = min(child->end, parent->end);
    
    
    for (vector<vert_edge>::iterator ei = edges_in.begin(); ei != edges_in.end(); ++ei){
    
        long int new_start = ei->start;
        long int new_end = ei->end;
    
        set<vert_edge>::iterator edge;
        set<vert_edge>::iterator edge_end;
    
        if (downward){
            edge = parent->children.begin();
            edge_end = parent->children.end();
        }
        else{
            edge = child->parents.begin();
            edge_end = child->parents.end();
        }
    
        bool add_edge = true;
        
        bool edge_modified = false;
        
        for (; edge != edge_end;){
            bool remove_edge = false;
            
            if (ranges_overlap(edge->start, edge->end, new_start, new_end)){
                short comp;
                // Goal: if downward, replace edges with things child is a superset of
                // if upward, replace edges with things parent is a subset of
                // Strategy: if comp is 1, replace edge with node. If 2 or 3, leave edge in place.
                if (downward){
                    comp = compare_bitsets(child->clade, edge->to->clade);
                }
                else{
                    comp = compare_bitsets(edge->to->clade, parent->clade);
                }
                if (comp == -1){
                    fprintf(stderr, "sp_aux(): new child fails 4hap with old child\n");
                    print_node_upward_onelevel(child, 0, num_haplotypes);
                    print_node_downward_onelevel(parent, 0, num_haplotypes);
                    fprintf(stderr, "\n");
                    print_node_lite(edge->to, num_haplotypes);
                    exit(1);
                    return false;
                }
                else if (comp == 1 || comp == 2 || comp == 3){
                    // Determine type of overlap.
                    if (edge->start == new_start && edge->end == new_end){
                        // Identical ranges.
                        if (comp == 1){
                            remove_edge = true;
                            edges_discarded[edge->to].push_back(*edge);
                        }
                        else{
                            nodes_blocking.insert(edge->to);
                            add_edge = false;
                            break;
                        }
                    }
                    else if (edge->start >= new_start && edge->end <= new_end){
                        // edge within node range
                        if (comp == 1){
                            remove_edge = true;
                            edges_discarded[edge->to].push_back(*edge);
                        }
                        else{
                            nodes_blocking.insert(edge->to);
                            if (edge->start > new_start){
                                arg_node* ptr;
                                if (downward){
                                    ptr = child;
                                }
                                else{
                                    ptr = parent;
                                }
                                vert_edge edge_new(new_start, edge->start-1, ptr);
                                edges_new.push_back(edge_new);
                            }
                            new_start = edge->end + 1;
                            edge_modified = true;
                        }
                    }
                    else if (new_start >= edge->start && new_end <= edge->end){
                        // node within edge range
                        if (comp == 1){
                            // Replace for this part of the range
                            vert_edge edge_disc(*edge);
                            if (new_start > edge->start){
                                vert_edge edge_new(edge->start, new_start-1, edge->to);
                                edges_old_adjusted.push_back(edge_new);
                                edge_disc.start = new_start;
                            }
                            if (new_end < edge->end){
                                vert_edge edge_new(new_end+1, edge->end, edge->to);
                                edges_old_adjusted.push_back(edge_new);
                                edge_disc.end = new_end;
                            }
                            remove_edge = true;
                            
                            edges_discarded[edge_disc.to].push_back(edge_disc);
                        }
                        else{
                            nodes_blocking.insert(edge->to);
                            add_edge = false;
                            break;
                        }
                    }
                    else if (edge->start < new_start && edge->end >= new_start &&
                        edge->end < new_end){
                        // edge hangs off to left of node range  
                        if (comp == 1){
                            vert_edge edge_new(edge->start, new_start-1, edge->to);
                            edges_old_adjusted.push_back(edge_new);
                            remove_edge = true;
                            vert_edge edge_disc(*edge);
                            edge_disc.start = new_start;
                            
                            edges_discarded[edge_disc.to].push_back(edge_disc);
                        }
                        else{
                            nodes_blocking.insert(edge->to);
                            new_start = edge->end + 1;
                            edge_modified = true;
                        }
                    }
                    else if (new_start < edge->start && new_end >= edge->start &&
                        new_end < edge->end){
                        // edge hangs off to right of node range
                        if (comp == 1){
                            vert_edge edge_new(new_end+1, edge->end, edge->to);
                            edges_old_adjusted.push_back(edge_new);
                            remove_edge = true;
                            vert_edge edge_disc(*edge);
                            edge_disc.end = new_end;
                            
                            edges_discarded[edge_disc.to].push_back(edge_disc);
                        }
                        else{
                            nodes_blocking.insert(edge->to);
                            new_end = edge->start - 1;
                            edge_modified = true;
                        }
                    }
                }
            }

            if (remove_edge && make_changes){
                // Remove corresponding edge & delete edge
                if (downward){
                    
                    bool found = false;
                    for (set<vert_edge>::iterator e2 = edge->to->parents.begin(); e2 != edge->to->parents.end();){
                        if (e2->start == edge->start && e2->end == edge->end && e2->to == parent){
                            edge->to->parents.erase(e2++);
                            found = true;
                            break;
                        }
                        else{
                            ++e2;
                        }
                    }
                    if (!found){
                        fprintf(stderr, "!found1 (parent->child edge missing corresponding child->parent) <%ld %ld> \n", edge->start, edge->end);
                        fprintf(stderr, "this edge:\n");
                        print_node_lite(edge->to, num_haplotypes);
                        fprintf(stderr, "parent:\n");
                        print_node_lite(parent, num_haplotypes);
                        
                        exit(1);
                    }
                    
                    parent->children.erase(edge++);
                }
                else{
                    
                    bool found = false;
                    for (set<vert_edge>::iterator e2 = edge->to->children.begin(); e2 != edge->to->children.end();){
                        if (e2->start == edge->start && e2->end == edge->end && e2->to == child){
                            edge->to->children.erase(e2++);
                            found = true;
                            break;
                        }
                        else{
                            ++e2;
                        }
                    }
                    if (!found){
                        fprintf(stderr, "!found2 (child->parent edge missing corresponding parent->child) <%ld %ld>\n", edge->start, edge->end);
                        fprintf(stderr, "child:\n");
                        print_node_lite(child, num_haplotypes);
                        fprintf(stderr, "this edge:\n");
                        print_node_lite(edge->to, num_haplotypes);
                        exit(1);
                    }
                    
                    child->parents.erase(edge++);
                }
            }
            else{
                ++edge;
            }
        }
        
        if (edge_modified && edges_already_found && make_changes){
            fprintf(stderr, "sp_aux(): edge modified, but edges already set (%d):\n", downward);
            fprintf(stderr, "edges IN:\n");
            for (vector<vert_edge>::iterator e = edges_in.begin(); e != edges_in.end(); ++e){
                fprintf(stderr, "<%ld %ld>\n", e->start, e->end);
                print_node_lite(e->to, num_haplotypes);
            }
            fprintf(stderr, "edges OUT:\n");
            fprintf(stderr, "<%ld %ld>\n", new_start, new_end);
            print_node_lite(ei->to, num_haplotypes);
            fprintf(stderr, "edges new:\n");
            for (vector<vert_edge>::iterator e = edges_new.begin(); e != edges_new.end(); ++e){
                fprintf(stderr, "<%ld %ld>\n", e->start, e->end);
                print_node_lite(e->to, num_haplotypes);
            }
            fprintf(stderr, "child:\n");
            print_node_lite(child, num_haplotypes);
            //print_node_upward_onelevel(child, 0, num_haplotypes);
            fprintf(stderr, "parent:\n");
            print_node_lite(parent, num_haplotypes);
            fprintf(stderr, "nodes blocking:\n");
            for (set<arg_node*>::iterator nb = nodes_blocking.begin(); nb != nodes_blocking.end(); ++nb){
                print_node_lite(*nb, num_haplotypes);
                bool found = false;
                for (set<vert_edge>::iterator e = child->parents.begin(); e != child->parents.end(); ++e){
                    if (e->to == *nb){
                        fprintf(stderr, "found <%ld %ld>\n", e->start, e->end);
                        found = true;
                    }
                }
                fprintf(stderr, "found %d\n\n", found);
                //print_recombs_right(*nb, num_haplotypes);
                fprintf(stderr, "\n");
            }
            exit(1);
        }
    
        // Create "main" new edge.
        if (add_edge && new_end >= new_start){
            if (downward){
                vert_edge edge_new(new_start, new_end, child);
                edge_new.context_start = cur_context_start;
                edge_new.context_end = cur_context_end;
                edges_new.push_back(edge_new);
            }
            else{
                vert_edge edge_new(new_start, new_end, parent);
                //edge_new.context_start = context_start;
                //edge_new.context_end = context_end;
                edges_new.push_back(edge_new);
            }
        }
    }
    
    if (!make_changes){
        // Bail out early rather than editing the ARG.
        return true;
    }
        
    
    // Handle new pieces of old edges.
    for (vector<vert_edge>::iterator adj = edges_old_adjusted.begin(); adj != edges_old_adjusted.end(); ++adj){
        if (downward){
            long int c_start = max(adj->to->start, parent->start);
            long int c_end = min(adj->to->end, parent->end);
            adj->context_start = c_start;
            adj->context_end = c_end;
            parent->children.insert(*adj);
            vert_edge adj2(adj->start, adj->end, parent);
            adj2.context_start = c_start;
            adj2.context_end = c_end;
            adj->to->parents.insert(adj2);
        }
        else{
            long int c_start = max(adj->to->start, child->start);
            long int c_end = min(adj->to->end, child->end);
            child->parents.insert(*adj);
            vert_edge adj2(adj->start, adj->end, child);
            adj2.context_start = c_start;
            adj2.context_end = c_end;
            adj->to->children.insert(adj2);
        }
    }
    
    // Handle new edges.
    
    long int c_start_main = max(child->start, parent->start);
    long int c_end_main = min(child->end, parent->end);
    
    for (vector<vert_edge>::iterator edge_new = edges_new.begin(); edge_new != edges_new.end(); ++edge_new){
        // Check to see that this edge doesn't overlap with another edge to 
        // the same node.
        for (vector<vert_edge>::iterator en2 = edges_new.begin(); en2 != edges_new.end();
            ++en2){
            
            if (DEBUG_MODE && ranges_overlap(edge_new->start, edge_new->end, en2->start, en2->end) &&
                edge_new != en2 && edge_new->to == en2->to){
                fprintf(stderr, "new edges overlap?\n");
                fprintf(stderr, "<%ld %ld>\n", edge_new->start, edge_new->end);
                print_node_lite(edge_new->to, num_haplotypes);
                fprintf(stderr, "<%ld %ld>\n", en2->start, en2->end);
                print_node_lite(en2->to, num_haplotypes);
                exit(1);
            }
        }
        edge_new->context_start = c_start_main;
        edge_new->context_end = c_end_main;
        if (downward){
            if (DEBUG_MODE && (edge_new->start < parent->start || edge_new->end > parent->end)){
                fprintf(stderr, "new edge limits outside parent range\n");
                exit(1);
            }
            parent->children.insert(*edge_new);
        }
        else{
            if (DEBUG_MODE && (edge_new->start < child->start || edge_new->end > child->end)){
                fprintf(stderr, "new edge limits outside child range\n");
                exit(1);
            }
            child->parents.insert(*edge_new);
        }
    }
    vector<arg_node*> nodes_discarded;    
    for (map<arg_node*, vector<vert_edge> >::iterator d = edges_discarded.begin();
        d != edges_discarded.end(); ++d){
        nodes_discarded.push_back(d->first);
    }
    sort(nodes_discarded.begin(), nodes_discarded.end(), sort_nodes_size);
    
    if (DEBUG_MODE){
        fprintf(stderr, "discarded:\n");
        for (vector<arg_node*>::iterator n = nodes_discarded.begin(); n != nodes_discarded.end(); ++n){
            print_node_lite(*n, num_haplotypes);
        }
        fprintf(stderr, "\n");
    }
    
    for (vector<arg_node*>::iterator n = nodes_discarded.begin(); n != nodes_discarded.end(); ++n){
        
        if (downward){
            
            discarded.insert(*n);
        }
        else{
            
        }
    }
    
    return true;
}

/**
 * This function assumes a correct set of new child-parent edges (from above function)
 * and adjusts all other existing edges around them.
 */
void set_parents(arg_node* child, 
    arg_node* parent, 
    bool downward_first, 
    long int prop_dist,
    arg_sitemap& sites_pos,
    arg_clademap& sites_clade,
    arg_clademap& lv_grps,
    set<arg_node*, nodeptr_child_sort>& discarded,
    int num_haplotypes){
    if (child == parent){
        return;
    }
    
    if (child->start > child->end || child->start > *child->sites.begin() ||
        child->end < *child->sites.rbegin()){
        fprintf(stderr, "child range messed up\n");
        print_node_lite(child, num_haplotypes);
        fprintf(stderr, "parent\n");
        print_node_lite(parent, num_haplotypes);
        exit(1);
    }
    
    else if (DEBUG_MODE && child->clade == parent->clade){
        fprintf(stderr, "set_parents(): child clade matches parent clade?\n");
        print_node_lite(child, num_haplotypes);
        fprintf(stderr, "\n");
        print_node_lite(parent, num_haplotypes);
        exit(1);
    }
    
    if (DEBUG_MODE){
        fprintf(stderr, "set parents:\n");
        
        print_node_lite(child, num_haplotypes);
        fprintf(stderr, "\n");
        print_node_lite(parent, num_haplotypes);
        fprintf(stderr, "\n");
    }
    
    
    if (!DEBUG_MODE){
        vector<vert_edge> edges_in1;
        vector<vert_edge> edges_out1;
        vector<vert_edge> edges_out2;
        
        set<arg_node*> nodes_blocking1;
        set<arg_node*> nodes_blocking2;
        
        // Insert child -> parent direction first and store edges
        sp_aux(child, parent, edges_in1, edges_out1, false, true,
            prop_dist, nodes_blocking1, sites_pos, sites_clade, lv_grps, discarded, num_haplotypes);
        
        // Use edges just learned to create parent -> child edges
        sp_aux(child, parent, edges_out1, edges_out2, true, true, 
            prop_dist, nodes_blocking1, sites_pos, sites_clade, lv_grps, discarded, num_haplotypes);
    }
    else{
        // Strategy: look in both directions and use the more conservative set of edges.
        // To do: use this strategy ONLY when calling from shorten_range(). Otherwise, do 
        // one of the commented out strategies.
        vector<vert_edge> edges_in1;
        vector<vert_edge> edges_out1;
        vector<vert_edge> edges_in2;
        vector<vert_edge> edges_out2;
        
        set<arg_node*> nodes_blocking1;
        set<arg_node*> nodes_blocking2;
        
        bool needs_add = sp_aux(child, parent, edges_in1, edges_out1, false, false, 
            prop_dist, nodes_blocking1, sites_pos, sites_clade, lv_grps, discarded, num_haplotypes);
        if (!needs_add){
            return;
        }

        needs_add = sp_aux(child, parent, edges_in2, edges_out2, true, false, prop_dist, 
            nodes_blocking2, sites_pos, sites_clade, lv_grps, discarded, num_haplotypes);
        if (!needs_add){
            return;
        }
        
        // Reconcile edges.
        
        set<arg_node*> nb1_unique;
        set_difference(nodes_blocking1.begin(), nodes_blocking1.end(), nodes_blocking2.begin(),
            nodes_blocking2.end(), inserter(nb1_unique, nb1_unique.begin()));
        set<arg_node*> nb2_unique;
        set_difference(nodes_blocking2.begin(), nodes_blocking2.end(), nodes_blocking1.begin(),
            nodes_blocking1.end(), inserter(nb2_unique, nb2_unique.begin()));
        
        vector<vert_edge> edges_out_consensus;
        bool match = true;
        for (vector<vert_edge>::iterator eo1 = edges_out1.begin(); eo1 != edges_out1.end(); ++eo1){
            bool found = false;
            for (vector<vert_edge>::iterator eo2 = edges_out2.begin(); eo2 != edges_out2.end(); ++eo2){
                if (eo2->start == eo1->start && eo2->end == eo1->end){
                    found = true;
                }
                if (ranges_overlap(eo1->start, eo1->end, eo2->start, eo2->end)){
                    // Overlap.
                    // Take smallest possible interval.
                    vert_edge eo_cons(max(eo1->start, eo2->start), min(eo1->end, eo2->end), eo1->to);
                    edges_out_consensus.push_back(eo_cons);
                }
            }
            if (!found){
                match = false;
            }
        }
        
        if (!match){
            fprintf(stderr, "edges don't match\n");
            print_node_upward_onelevel(child, 0, num_haplotypes);
            print_node_downward_onelevel(parent, 0, num_haplotypes);
            fprintf(stderr, "\n");
            fprintf(stderr, "eo1\n");
            for (vector<vert_edge>::iterator eo1 = edges_out1.begin(); eo1 != edges_out1.end(); ++eo1){
                fprintf(stderr, "<%ld %ld>\n", eo1->start, eo1->end);
            }
            fprintf(stderr, "eo2\n");
            for (vector<vert_edge>::iterator eo2 = edges_out2.begin(); eo2 != edges_out2.end(); ++eo2){
                fprintf(stderr, "<%ld %ld>\n", eo2->start, eo2->end);
            }
            
            exit(1);
        }
        
        if (edges_out_consensus.size() == 0){
            return;
        }
        vector<vert_edge> edges_out3;
        vector<vert_edge> edges_out4;
        
        set<arg_node*> nodes_blocked3;
        set<arg_node*> nodes_blocked4;
        
        sp_aux(child, parent, edges_out_consensus, edges_out3, false, true, prop_dist, 
            nodes_blocked3, sites_pos, sites_clade, lv_grps, discarded, num_haplotypes);

        sp_aux(child, parent, edges_out_consensus, edges_out4, true, true, prop_dist, 
            nodes_blocked4, sites_pos, sites_clade, lv_grps, discarded, num_haplotypes);
    
        
        if (DEBUG_MODE){
            fprintf(stderr, "done\n");
            print_node_upward_onelevel(child, 0, num_haplotypes);
        }
    }
    return;
}

void expand_range_fixed(arg_node* node, long int site,
    arg_sitemap& sites_pos,
    arg_clademap& sites_clade,
    arg_clademap& lv_grps,
    long int prop_dist,
    arg_node* root,
    vector<arg_node*>& newnodes,
    bool fix_end,
    int num_haplotypes,
    bool no_failures,
    set<arg_node*>& nodes_external){
    
    // This function assumes that no failures will be created (these ranges have
    // been pre-tested and are "final." Therefore, if the expansion is going to cause
    // one or more new four haplotype test failures that would affect the ranges of
    // one or more nodes, we need to bail out.
    fix_end = false;
    
    if (fix_end){
        pair<arg_sitemap::iterator, arg_sitemap::iterator> er = sites_pos.equal_range(site);
        for (arg_sitemap::iterator sp = er.first;
            sp != er.second; ++sp){
            short comp = compare_bitsets(sp->second->clade, node->clade);
            if (comp == -1){
                fprintf(stderr, "EXPRF failure\n");
                print_node_lite(node, num_haplotypes);
                fprintf(stderr, "%ld\n", site);
                print_node_lite(sp->second, num_haplotypes);
                exit(1);
            }
        }
    }
    
    bool l = false;
    
    if (site > node->end){
        long int newend = site;
        if (site > *node->sites.rbegin() + prop_dist){
            newend = *node->sites.rbegin() + prop_dist;
            fix_end = false;
        }
        if (node->recomb_right != -1 && newend > node->recomb_right){
            for (arg_sitemap::reverse_iterator it = 
                arg_sitemap::reverse_iterator(sites_pos.equal_range(node->recomb_right).first);
                it != sites_pos.rend(); ++it){
                if (it->first < node->recomb_right){
                    newend = it->first;
                    fix_end = false;
                    //site = it->first;
                    break;
                }
            }
        }
        if (newend > node->end){
            // Check to see if this will create failures that will change one or more
            // node's coordinates; bail out if so.
            
            if (no_failures){
                /*
                set<arg_node*> failures;
                compile_failures(node, root, failures, node->end, newend, prop_dist, sites_pos,
                    num_haplotypes, false, false);
                for (set<arg_node*>::iterator f = failures.begin(); f != failures.end(); ++f){
                    bool found = false;
                    for (vector<rc_edge>::iterator from = node->edges_from.begin();
                        from != node->edges_from.end(); ++from){
                        if (from->final_dest != NULL && from->final_dest == *f){
                            found = true;
                            break;
                        }
                    }
                    if (!found){
                        for (vector<rc_edge>::iterator from = node->edges_from_solved.begin();
                            from != node->edges_from_solved.end(); ++from){
                            if (from->final_dest != NULL && from->final_dest == *f){
                                found = true;
                                break;
                            }
                        }
                    }
                    if (!found){
                        for (vector<rc_edge>::iterator from = node->edges_from_unsolvable.begin();
                            from != node->edges_from_unsolvable.end(); ++from){
                            if (from->final_dest != NULL && from->final_dest == *f){
                                found = true;
                                break;
                            }
                        }
                    }
                    if (!found){
                        return;
                    }
                }
                */
            }
            
            
            // Check to see if we can merge with another node.
            bool finished = false;
            arg_node* mergenode = NULL;
            for (arg_clademap::iterator prevnode = sites_clade.equal_range(node->clade).first;
                prevnode != sites_clade.equal_range(node->clade).second; ++prevnode){
                if (prevnode->second != node && site_ranges_overlap(node->start, newend, prevnode->second) &&
                    (node->recomb_right == -1 || node->recomb_right > *prevnode->second->sites.begin()) &&
                    (prevnode->second->recomb_left == -1 || prevnode->second->recomb_left < *node->sites.rbegin())){
                    mergenode = prevnode->second;
                    break;
                }
            }
            if (mergenode != NULL){
                bool expanded = true;
                // Expand range to the limit of the other node, 
                // then merge nodes and delete the old one.
                if (mergenode->start-1 > node->end){
                    long int origend = node->end;
                    node->end = mergenode->start-1;
                    expand_range(node, node->end-origend, sites_pos, sites_clade,
                        lv_grps, prop_dist, root, newnodes, num_haplotypes, false);
                    if (node->end < mergenode->start-1){
                        origend = node->end;
                        node->end = mergenode->start-1;
                        expand_range(node, node->end-origend, sites_pos, sites_clade,
                            lv_grps, prop_dist, root, newnodes, num_haplotypes, false);
                        if (node->end < mergenode->start-1){
                           
                            expanded = false;
                        }
                    }
                }
                if (expanded){
                    join_nodes(node, mergenode, sites_pos, sites_clade, lv_grps,
                        prop_dist, false, root, newnodes, num_haplotypes);
                    // Wipe out the other node
                    for (arg_clademap::iterator prevnode = sites_clade.equal_range(node->clade).first;
                        prevnode != sites_clade.equal_range(node->clade).second;){
                        if (prevnode->second == mergenode){
                            sites_clade.erase(prevnode++);
                        }
                        else{
                            ++prevnode;
                        }
                    }
                    if (nodes_external.find(mergenode) != nodes_external.end()){
                        nodes_external.erase(nodes_external.find(mergenode));
                    }
                    delete mergenode;
                    update_recinds_left(node);
                    finished = true;
                }
            }
            if (!finished){
                long int origend = node->end;
                node->end = newend;
                expand_range(node, node->end - origend, 
                    sites_pos, sites_clade, lv_grps, prop_dist, root,
                    newnodes, num_haplotypes, true);
            }
        }
    }
    else if (site < node->start){
        l = true;
        long int newstart = site;
        if (site < *node->sites.begin() - prop_dist){
            newstart = *node->sites.begin() - prop_dist;
            fix_end = false;
        }
        if (node->recomb_left != -1 && newstart < node->recomb_left){
            for (arg_sitemap::iterator it = sites_pos.equal_range(node->recomb_left).second;
                it != sites_pos.end(); ++it){
                if (it->first > node->recomb_left){
                    newstart = it->first;
                    fix_end = false;
                    //site = it->first;
                    break;
                }
            }
        }
        if (newstart < node->start){
        
            // Check to see if this will create failures that will change one or more
            // node's coordinates; bail out if so.
            
            /*
            if (no_failures){
                set<arg_node*> failures;
                compile_failures(node, root, failures, newstart, node->start, prop_dist, sites_pos,
                    num_haplotypes, false, false);
                for (set<arg_node*>::iterator f = failures.begin(); f != failures.end(); ++f){
                    bool found = false;
                    for (vector<rc_edge>::iterator to = node->edges_to.begin();
                        to != node->edges_to.end(); ++to){
                        if (to->final_dest != NULL && to->final_dest == *f){
                            found = true;
                            break;
                        }
                    }
                    if (!found){
                        for (vector<rc_edge>::iterator to = node->edges_to_solved.begin();
                            to != node->edges_to_solved.end(); ++to){
                            if (to->final_dest != NULL && to->final_dest == *f){
                                found = true;
                                break;
                            }
                        }
                    }
                    if (!found){
                        for (vector<rc_edge>::iterator to = node->edges_to_unsolvable.begin();
                            to != node->edges_to_unsolvable.end(); ++to){
                            if (to->final_dest != NULL && to->final_dest == *f){
                                found = true;
                                break;
                            }
                        }
                    }
                    if (!found){
                        return;
                    }
                }
            }
            */
            
            // Check to see if we can merge with another node.
            bool finished = false;
            arg_node* mergenode = NULL;
            for (arg_clademap::iterator prevnode = sites_clade.equal_range(node->clade).first;
                prevnode != sites_clade.equal_range(node->clade).second; ++prevnode){
                if (prevnode->second != node && site_ranges_overlap(newstart, node->end, prevnode->second) &&
                    (node->recomb_left == -1 || node->recomb_left < *prevnode->second->sites.rbegin()) &&
                    (prevnode->second->recomb_right == -1 || prevnode->second->recomb_right > *node->sites.begin())){
                    mergenode = prevnode->second;
                    break;
                }
            }
            if (mergenode != NULL){
                bool expanded = true;
                // Expand range to the limit of the other node,
                // then merge nodes and delete the old one.
                if (mergenode->end+1 < node->start){
                    long int origstart = node->start;
                    node->start = mergenode->end+1;
                    expand_range(node, node->start-origstart, sites_pos, sites_clade,
                        lv_grps, prop_dist, root, newnodes, num_haplotypes, false);
                    if (node->start < mergenode->end+1){
                        origstart = node->start;
                        node->start = mergenode->end+1;
                        expand_range(node, node->start-origstart, sites_pos, sites_clade,
                            lv_grps, prop_dist, root, newnodes, num_haplotypes, false);
                        if (node->start < mergenode->end +1){
                           
                            expanded = false;
                        }
                    }
                }
                if (expanded){
                    join_nodes(node, mergenode, sites_pos, sites_clade, lv_grps,
                        prop_dist, false, root, newnodes, num_haplotypes);
                    // Wipe out the other node.
                    for (arg_clademap::iterator prevnode = sites_clade.equal_range(node->clade).first;
                        prevnode != sites_clade.equal_range(node->clade).second;){
                        if (prevnode->second == mergenode){
                            sites_clade.erase(prevnode++);
                            break;
                        }
                        else{
                            ++prevnode;
                        }
                    }
                    if (nodes_external.find(mergenode) != nodes_external.end()){
                        nodes_external.erase(nodes_external.find(mergenode));
                    }
                    delete mergenode;
                    update_recinds_right(node);
                    finished = true;
                }
            }
            if (!finished){
                long int origstart = node->start;
                node->start = newstart;
                expand_range(node, node->start - origstart,
                    sites_pos, sites_clade, lv_grps, prop_dist, root,
                    newnodes, num_haplotypes, false);
            }
        }
    }
    
    // Can this be deleted?
    // ===================
    
    if (fix_end){
        if (!l && site == node->end){
            node->sites.insert(site);
            bool found = false;
            pair<arg_sitemap::iterator, arg_sitemap::iterator> er = sites_pos.equal_range(site);
            for (arg_sitemap::iterator it = er.first;
                it != er.second; ++it){
                if (it->second == node){
                    found = true;
                    break;
                }
            }
            if (!found){
                sites_pos.insert(make_pair(site, node));
            }
            
            for (vector<rc_edge>::iterator from = node->edges_from.begin();
                from != node->edges_from.end(); ++from){
                if (from->final_dest != NULL){
                    if (from->final_dest->recomb_left < site){
                        from->final_dest->recomb_left = site;
                    }
                }
            }
            for (vector<rc_edge>::iterator from = node->edges_from_solved.begin();
                from != node->edges_from_solved.end(); ++from){
                if (from->final_dest != NULL){
                    if (from->final_dest->recomb_left < site){
                        from->final_dest->recomb_left = site;
                        short comp = compare_bitsets(node->clade, from->final_dest->clade);
                        if (comp == -1){   
                            if (from->final_dest->recomb_left > from->final_dest->start){
                                fprintf(stderr, "messed up inds1\n");
                                print_node_lite(node, num_haplotypes);
                                print_recombs_right(node, num_haplotypes);
                                fprintf(stderr, "\n");
                                print_node_lite(from->final_dest, num_haplotypes);
                                print_recombs_left(from->final_dest, num_haplotypes);
                                exit(1);
                            }
                        }
                    }
                }
            }
            for (vector<rc_edge>::iterator from = node->edges_from_unsolvable.begin();
                from != node->edges_from_unsolvable.end(); ++from){
                if (from->final_dest != NULL && from->final_dest->recomb_left < site){
                    from->final_dest->recomb_left = site;
                }
            }
            
        }
        else if (l && site == node->start){
            node->sites.insert(site);
            bool found = false;
            pair<arg_sitemap::iterator, arg_sitemap::iterator> er = sites_pos.equal_range(site);
            for (arg_sitemap::iterator it = er.first;
                it != er.second; ++it){
                if (it->second == node){
                    found = true;
                    break;
                }
            }
            if (!found){
                sites_pos.insert(make_pair(site, node));
            }
            
            for (vector<rc_edge>::iterator to = node->edges_to.begin();
                to != node->edges_to.end(); ++to){
                if (to->final_dest != NULL){
                    if (to->final_dest->recomb_right > site){
                        to->final_dest->recomb_right = site;
                    }
                }
            }
            for (vector<rc_edge>::iterator to = node->edges_to_solved.begin();
                to != node->edges_to_solved.end(); ++to){
                if (to->final_dest != NULL){
                    if (to->final_dest->recomb_right > site){
                        short comp = compare_bitsets(node->clade, to->final_dest->clade);
                        if (comp == -1){
                            to->final_dest->recomb_right = site;
                            if (to->final_dest->recomb_right < to->final_dest->end){
                                fprintf(stderr, "messed up inds2\n");
                                print_node_lite(node, num_haplotypes);
                                print_recombs_left(node, num_haplotypes);
                                fprintf(stderr, "\n");
                                print_node_lite(to->final_dest, num_haplotypes);
                                print_recombs_right(to->final_dest, num_haplotypes);
                                exit(1);
                            }
                        }
                    }
                }   
            }
            for (vector<rc_edge>::iterator to = node->edges_to_unsolvable.begin();
                to != node->edges_to_unsolvable.end(); ++to){
                if (to->final_dest != NULL && to->final_dest->recomb_right > site){
                    to->final_dest->recomb_right = site;
                }
            }
            
        }
    }
    else{
        
    }

}

void update_recinds_left(arg_node* l){
    for (vector<rc_edge>::iterator from = l->edges_from.begin();
        from != l->edges_from.end(); ++from){
        if (from->final_dest != NULL){
            if (*l->sites.rbegin() > from->final_dest->recomb_left){
                from->final_dest->recomb_left = *l->sites.rbegin();;
            }
        }
    }
    for (vector<rc_edge>::iterator from = l->edges_from_solved.begin();
        from != l->edges_from_solved.end(); ++from){
        if (from->final_dest != NULL){
            if (*l->sites.rbegin() > from->final_dest->recomb_left &&
                compare_bitsets(l->clade, from->final_dest->clade) == -1){
                from->final_dest->recomb_left = *l->sites.rbegin();
            }
        }
    }
    for (vector<rc_edge>::iterator from = l->edges_from_unsolvable.begin();
        from != l->edges_from_unsolvable.end(); ++from){
        if (from->final_dest != NULL){
            if (*l->sites.rbegin() > from->final_dest->recomb_left){
                from->final_dest->recomb_left = *l->sites.rbegin();
            }
        }
    }
}

void update_recinds_right(arg_node* r){
    for (vector<rc_edge>::iterator to = r->edges_to.begin(); to != r->edges_to.end();
        ++to){
        if (to->final_dest != NULL && *r->sites.begin() < to->final_dest->recomb_right){
            to->final_dest->recomb_right = *r->sites.begin();
        }
    }
    for (vector<rc_edge>::iterator to = r->edges_to_solved.begin();
        to != r->edges_to_solved.end(); ++to){
        if (to->final_dest != NULL && *r->sites.begin() < to->final_dest->recomb_right &&
            compare_bitsets(r->clade, to->final_dest->clade) == -1){
            to->final_dest->recomb_right = *r->sites.begin();
        }
    }
    for (vector<rc_edge>::iterator to = r->edges_to_unsolvable.begin();
        to != r->edges_to_unsolvable.end(); ++to){
        if (to->final_dest != NULL && *r->sites.begin() < to->final_dest->recomb_right){
            to->final_dest->recomb_right = *r->sites.begin();
        }
    }
}

void unfix_edges_right(arg_node* node){
    node->fixed_right = false;
    for (vector<rc_edge>::iterator from = node->edges_from.begin(); from != node->edges_from.end(); ++from){
        if (from->final_dest != NULL){
            from->final_dest->fixed_left = false;
        }
    }
    for (vector<rc_edge>::iterator from = node->edges_from_solved.begin(); from != node->edges_from_solved.end(); ++from){
        if (from->final_dest != NULL){
            from->final_dest->fixed_left = false;
        }
    }
    for (vector<rc_edge>::iterator from = node->edges_from_unsolvable.begin(); from != node->edges_from_unsolvable.end(); ++from){
        if (from->final_dest != NULL){
            from->final_dest->fixed_left = false;
        }
    }
}

void unfix_edges_left(arg_node* node){
    node->fixed_left = false;
    for (vector<rc_edge>::iterator to = node->edges_to.begin(); to != node->edges_to.end(); ++to){
        if (to->final_dest != NULL){
            to->final_dest->fixed_right = false;
        }
    }
    for (vector<rc_edge>::iterator to = node->edges_to_solved.begin(); to != node->edges_to_solved.end(); ++to){
        if (to->final_dest != NULL){
            to->final_dest->fixed_right = false;
        }
    }
    for (vector<rc_edge>::iterator to = node->edges_to_unsolvable.begin(); to != node->edges_to_unsolvable.end(); ++to){
        if (to->final_dest != NULL){
            to->final_dest->fixed_right = false;
        }
    }
}


/**
 * amount == new index - old index
 * amount < 0: expanded to left
 * amount > 0: expanded to right
 */
void expand_range(arg_node* node, long int amount, 
    arg_sitemap& sites_pos,
    arg_clademap& sites_clade,
    arg_clademap& lv_grps,
    long int prop_dist,
    arg_node* root,
    vector<arg_node*>& newnodes,
    int num_haplotypes,
    bool newsites){
    if (node->end < node->start || node->start > *node->sites.begin() ||
        node->end < *node->sites.rbegin()){
        fprintf(stderr, "expand_range - invalid range\n");
        print_node_lite(node, num_haplotypes);
        bool found = false;
        for (arg_clademap::iterator sc = sites_clade.equal_range(node->clade).first;
            sc != sites_clade.equal_range(node->clade).second; ++sc){
            if (sc->second == node){
                found = true;
                break;
            }
        }
        fprintf(stderr, "found1 %d\n", found);
        bool found2 = false;
        for (arg_clademap::iterator lv = lv_grps.equal_range(node->clade).first;
            lv != lv_grps.equal_range(node->clade).second; ++lv){
            if (lv->second == node){
                found2 = true;
                break;
            }
        }
        fprintf(stderr, "found2 %d\n", found2);
        fprintf(stderr, "%ld parents %ld children\n", node->parents.size(), node->children.size());
        exit(1);   
    }
    
    if (amount < 0 && node->start < *node->sites.begin() - prop_dist){
        amount += (*node->sites.begin() - prop_dist - node->start);
        node->start = *node->sites.begin() - prop_dist;
    }
    if (amount > 0 && node->end > *node->sites.rbegin() + prop_dist){
        amount -= (node->end - (*node->sites.rbegin() + prop_dist));
        node->end = *node->sites.rbegin() + prop_dist;
    }
    
    if (DEBUG_MODE && amount > 0 && node->recomb_right != -1 && node->recomb_right < node->end){
        fprintf(stderr, "node's right edge extends past recombination\n");
        print_node_lite(node, num_haplotypes);
        exit(1);
    }
    else if (DEBUG_MODE && amount < 0 && node->recomb_left != -1 && node->recomb_left > node->start){
        fprintf(stderr, "node's left edge extends past recombination\n");
        print_node_lite(node, num_haplotypes);
        exit(1);
    }
    
    if (DEBUG_MODE){
        fprintf(stderr, "EXPAND RANGE; %ld\n", amount);
        print_node_lite(node, num_haplotypes);
    }

    long int start;
    long int end;
    
    bool exp_l = false;
    bool exp_r = false;
    
    if (amount > 0){
        start = node->end - amount;
        end = node->end;
        exp_r = true;
        unfix_edges_right(node);
    }
    else if (amount < 0){
        start = node->start;
        end = node->start - amount;
        exp_l = true;
        unfix_edges_left(node);
    }
    else{
        // zero? do nothing.
        return;
    }
    
    if (!newsites){
        exp_l = false;
        exp_r = false;
    }
    
    if (DEBUG_MODE){
        fprintf(stderr, "haf/gp2 %ld %ld\n", start, end);
    }
    
    long int origstart = node->start;
    long int origend = node->end;
    
    set<arg_node*> failures;
    compile_failures(node, root, failures, start, end, prop_dist, sites_pos, num_haplotypes, exp_l, exp_r);
    
    // Eliminate anything we already know about?
    
    adjust_all_failure_indices(node, failures, sites_pos, sites_clade, lv_grps, 
        prop_dist, root, newnodes, num_haplotypes);
    
    // Change amount, if we were blocked by a failure.
    if (amount < 0){
        long int before_func_start = origstart - amount;
        if (node->start > origstart){
            if (node->start > before_func_start){
                // We now are no longer actually expanding the range at all. 
                // No need to call shorten_range() since that will be taken care of
                // within handle_all_failures().
                return;
            }
            else{
                amount = node->start - before_func_start;
            }
        }
    }
    else if (amount > 0){
        long int before_func_end = origend - amount;
        if (node->end < origend){
            if (node->end < before_func_end){
                // We are now longer expanding the range at all.
                return;
            }
            else{
                amount = node->end - before_func_end;
            }
        }
    }
    
    // Handling failures might change coordinates and make gathering new parents unnecessary.
    bool need_new_parents = true;
    if (amount < 0){
        if (end < node->start){
            // Don't need new parents.
            need_new_parents = false;
        }
        else if (start < node->start){
            start = node->start;
        }
    }
    else if (amount > 0){
        if (start > node->end){
            need_new_parents = false;
        }
        else if (end > node->end){
            end = node->end;
        }
    }
    
    // Need to tell all parents & children about this node's expanded range.
    set<arg_node*, nodeptr_parent_sort> parents;
    for (set<vert_edge>::iterator edge = node->parents.begin(); edge != node->parents.end();){
        bool remove = false;
        if (amount < 0 && edge->start > edge->to->start && edge->to->start < end){
            remove = true;
        }
        else if (amount > 0 && edge->end < edge->to->end && edge->to->end > start){
            remove = true;
        }
        //remove = true;
        if (remove){
            parents.insert(edge->to);
            for (set<vert_edge>::iterator e2 = edge->to->children.begin(); e2 != edge->to->children.end();){
                if (e2->start == edge->start && e2->end == edge->end && e2->to == node){
                    edge->to->children.erase(e2++);
                    break;
                }
                else{
                    ++e2;
                }
            }
            node->parents.erase(edge++);
        }
        else{
            ++edge;
        }
    }
    
    set<arg_node*, nodeptr_child_sort> children;
    for (set<vert_edge>::iterator edge = node->children.begin(); edge != node->children.end();){
        bool remove = false;
        if (amount < 0 && edge->start > edge->to->start && edge->to->start < end){
            remove = true;
        }
        else if (amount > 0 && edge->end < edge->to->end && edge->to->end > start){
            remove = true;
        }
        //remove = true;
        if (remove){
            children.insert(edge->to);
            for (set<vert_edge>::iterator e2 = edge->to->parents.begin(); e2 != edge->to->parents.end();){
                if (e2->start == edge->start && e2->end == edge->end && e2->to == node){
                    edge->to->parents.erase(e2++);
                    break;
                }
                else{
                    ++e2;
                }
            }
            node->children.erase(edge++);
        }
        else{
            ++edge;
        }
    }
    
    if (DEBUG_MODE){
        fprintf(stderr, "gp2 %ld %ld\n", start, end);
    }
    
    if (need_new_parents){
        gather_parents(node, root, start, end, prop_dist, sites_pos, sites_clade, lv_grps, parents, num_haplotypes);
    }
    
    for (set<arg_node*, nodeptr_parent_sort>::iterator parent = parents.begin(); parent != parents.end();
        ++parent){
        if (node_ranges_overlap(node, *parent)){
            set_parents(node, *parent, true, prop_dist, sites_pos, sites_clade, lv_grps, children, num_haplotypes);
        }
    }
    set<arg_node*, nodeptr_child_sort> dummy;
    for (set<arg_node*, nodeptr_child_sort>::iterator child = children.begin(); child != children.end();
        ++child){
        if (node_ranges_overlap(node, *child)){
            set_parents(*child, node, false, prop_dist, sites_pos, sites_clade, lv_grps, dummy, num_haplotypes);
        }
    }
    
    if (DEBUG_MODE){
        fprintf(stderr, "expr handling failures\n");
    }
    
    handle_failures(node, failures, sites_pos, sites_clade, lv_grps, prop_dist, root, newnodes, num_haplotypes);
    
    if (DEBUG_MODE){
        if (dummy.size() > 0){
            fprintf(stderr, "?? dummy > 0\n");
            exit(1);
        }
        
        fprintf(stderr, "EXPAND RANGE DONE\n");
    }
    
}

void find_new_parents(arg_node* node, 
    arg_sitemap& sites_pos,
    arg_clademap& sites_clade,
    arg_clademap& lv_grps,
    long int prop_dist,
    long int start,
    long int end,
    arg_node* root,
    vector<arg_node*>& newnodes,
    int num_haplotypes){
    
    if (start < node->start){
        start = node->start;
    }
    if (end > node->end){
        end = node->end;
    }

    insert_below_parent(node, root, sites_pos, sites_clade, lv_grps, prop_dist, 
        -1, -1, root, newnodes, num_haplotypes);
}

void find_new_children(arg_node* node,
    arg_sitemap& sites_pos,
    arg_clademap& sites_clade,
    arg_clademap& lv_grps,
    long int prop_dist,
    long int start, long int end,
    arg_node* root,
    vector<arg_node*>& newnodes,
    int num_haplotypes){
    
    if (start < node->start){
        start = node->start;
    }
    if (end > node->end){
        end = node->end;
    }

    insert_below_parent(node, root, sites_pos, sites_clade, lv_grps, prop_dist, 
        -1, -1, root, newnodes, num_haplotypes);
}

/** 
 * amount == new index - old index
 * amount < 0: right side
 * amount > 0: left side
 */
void shorten_range(arg_node* node, long int amount, 
    arg_sitemap& sites_pos,
    arg_clademap& sites_clade,
    arg_clademap& lv_grps,
    long int prop_dist,
    int num_haplotypes,
    arg_node* root,
    vector<arg_node*>& newnodes,
    long int del_lim){
    //fprintf(stderr, "sr\n");
    if (!(node->start == -1 && node->end == -1)){
        if (node->end < *node->sites.rbegin()){
            fprintf(stderr, "shorten range - invalid range\n");
            print_node_lite(node, num_haplotypes);
            exit(1);
        }
        if (node->start > *node->sites.begin()){
            fprintf(stderr, "shorten range - invalid range\n");
            print_node_lite(node, num_haplotypes);
            exit(1);
        }
        if (node->end < node->start){
            fprintf(stderr, "shorten range - invalid range\n");
            print_node_lite(node, num_haplotypes);
            exit(1);
        }
    }
    
    if (amount == 0 && !(node->start == -1 && node->end == -1)){
        return;
    }
    
    if (DEBUG_MODE){
        fprintf(stderr, "SHORTEN RANGE (%ld)\n", amount);
        print_node_downward_onelevel(node, 0, num_haplotypes);
        print_node_upward_onelevel(node, 0, num_haplotypes);
    }
    
    // Store parents who are affected by the change in range
    set<arg_node*, nodeptr_parent_sort> parents;
    set<arg_node*, nodeptr_child_sort> children;
    
    for (set<vert_edge>::iterator parent = node->parents.begin(); parent != node->parents.end();){
        bool removed = false;
        
        if (node->start == -1 && node->end == -1){
            parents.insert(parent->to);
            removed = true;
        }
        else{
            if (amount > 0){
                if (parent->start < node->start){
                    // This edge is affected.
                    parents.insert(parent->to);
                    removed = true;
                }
            }
            else if (amount < 0){
                if (parent->end > node->end){
                    // This edge is affected.
                    parents.insert(parent->to);
                    removed = true;
                }
            }
        }
        if (removed){
            for (set<vert_edge>::iterator e2 = parent->to->children.begin(); e2 != parent->to->children.end();){
                if (e2->start == parent->start && e2->end == parent->end && e2->to == node){
                    parent->to->children.erase(e2++);
                    break;
                }
                else{
                    ++e2;
                }
            } 
            
            node->parents.erase(parent++);  
        }
        else{
            ++parent;
        }
    }
    
    
    
    for (set<vert_edge>::iterator child = node->children.begin(); child != node->children.end();){
        bool removed = false;
        if (node->start == -1 && node->end == -1){
            children.insert(child->to);
            removed = true;
        }
        else{
            if (amount > 0){
                if (child->start < node->start){
                    // This edge is affected.
                    children.insert(child->to);
                    removed = true;
                }
            }
            else if (amount < 0){
                if (child->end > node->end){
                    // This edge is affected.
                    children.insert(child->to);
                    removed = true;
                }
            }
        }
        if (removed){
            for (set<vert_edge>::iterator e2 = child->to->parents.begin(); e2 != child->to->parents.end();){
                if (e2->start == child->start && e2->end == child->end && e2->to == node){
                    child->to->parents.erase(e2++);
                    break;
                }
                else{
                    ++e2;
                }
            } 
            node->children.erase(child++);  
        }
        else{
            ++child;
        }
    }
    
    // Add back parents in proper order.
    for (set<arg_node*, nodeptr_parent_sort>::iterator parent = parents.begin(); parent != parents.end(); ++parent){
        if (amount > 0 && (*parent)->end < node->start){
            // Don't add back to node.
        }
        else if (amount < 0 && (*parent)->start > node->end){
            // Don't add back to node.
            if (DEBUG_MODE){
                fprintf(stderr, "skip parent\n");
                print_node_lite(*parent, num_haplotypes);
            }
        }
        else if (amount == 0 && node->start == -1 && node->end == -1){
            // Don't add back; deleting.
        }
        else{
            if (DEBUG_MODE){
                fprintf(stderr, "ADDING BACK PARENT:\n");
                print_node_lite(*parent, num_haplotypes);
            }
            //fprintf(stderr, "sp6\n");
            set_parents(node, *parent, false, prop_dist, sites_pos, sites_clade, lv_grps, children, num_haplotypes);
        }
        
    }
    set<arg_node*, nodeptr_child_sort> dummy;
    // Add back children in proper order.
    for (set<arg_node*, nodeptr_child_sort>::iterator child = children.begin(); child != children.end(); ++child){
        if (amount > 0 && (*child)->end < node->start){
            // Don't add back to node.
            if (DEBUG_MODE){
                fprintf(stderr, "skip child\n");
                print_node_lite(*child, num_haplotypes);
            }
        }
        else if (amount < 0 && (*child)->start > node->end){
            // Don't add back to node.
            if (DEBUG_MODE){
                fprintf(stderr, "skip child\n");
                print_node_lite(*child, num_haplotypes);
            }
        }
        else if (amount == 0 && node->start == -1 && node->end == -1){
            // Don't add back; deleting.
        }
        else{
            if (DEBUG_MODE){
                fprintf(stderr, "ADDING BACK CHILD:\n");
                print_node_lite(*child, num_haplotypes);
            }
            //fprintf(stderr, "sp7\n");
            set_parents(*child, node, false, prop_dist, sites_pos, sites_clade, lv_grps, dummy, num_haplotypes);
        }
    }
    if (dummy.size() > 0){
        fprintf(stderr, "??? dummy > 0 (2)\n");
        exit(1);
    }
    
    if (amount == 0 && node->start == -1 && node->end == -1){
        for (set<arg_node*, nodeptr_child_sort>::iterator child = children.begin();
            child != children.end(); ++child){
            if (del_lim == -1 || *(*child)->sites.rbegin() >= del_lim){
                insert_below_parent(*child, root, sites_pos, sites_clade, lv_grps,
                    prop_dist, -1, -1, root, newnodes, num_haplotypes);
            }
        }
    }
    
    for (set<arg_node*, nodeptr_parent_sort>::iterator parent = parents.begin(); parent != parents.end(); ++parent){
        // See if it needs to pick up one of the children on the end of its range.
        if (amount == 0 && node->start == -1 && node->end == -1){
            
            
        }
        else if (amount > 0 && (*parent)->start < node->start){
            for (set<arg_node*, nodeptr_child_sort>::iterator child = children.begin();
                child != children.end(); ++child){
                if ((*child)->start < node->start && node_ranges_overlap(*parent, *child)){
                    //fprintf(stderr, "sp8\n");
                    set_parents(*child, *parent, false, prop_dist, sites_pos, sites_clade, lv_grps, dummy, num_haplotypes);
                    //insert_below_parent(*child, *parent, sites_pos, sites_clade, lv_grps, 
                    //    prop_dist, -1, -1, root, newnodes, num_haplotypes);
                    if (dummy.size() > 0){
                        fprintf(stderr, "??? dummy > 0 (3)\n");
                        exit(1);
                    }
                }
            }
        }
        else if (amount < 0 && (*parent)->end > node->end){
            for (set<arg_node*, nodeptr_child_sort>::iterator child = children.begin();
                child != children.end(); ++child){
                if ((*child)->end > node->end && node_ranges_overlap(*parent, *child)){
                    //fprintf(stderr, "sp9\n");
                    set_parents(*child, *parent, false, prop_dist, sites_pos, sites_clade, lv_grps, dummy, num_haplotypes);
                    //insert_below_parent(*child, *parent, sites_pos, sites_clade, lv_grps, 
                    //    prop_dist, -1, -1, root, newnodes, num_haplotypes);
                    if (dummy.size() > 0){
                        fprintf(stderr, "??? dummy > 0 (4)\n");
                        exit(1);
                    }
                }   
            }
        }
    }
   
   
    if (DEBUG_MODE){
        fprintf(stderr, "SHORTEN RANGE DONE\n");
        print_node_downward_onelevel(node, 0, num_haplotypes);
        print_node_upward_onelevel(node, 0, num_haplotypes);
    }
    
    // Notify recombination partners.
    if (amount < 0){
        //unfix_edges_right(node);
        //node->fixed_right = false;
    }
    else if (amount > 0){
        //unfix_edges_left(node);
        //node->fixed_left = false;
    }
    
    return;
}

bool lv_grp_compatible(arg_node* left_node, 
    arg_node* right_node,
    int type,
    arg_node* lv_grp){
    // Collect all preexisting failures.
    vector<pair<arg_node*, arg_node*> > failures;
    vector<int> failuretypes;
    set<arg_node*> lnodes;
    for (vector<rc_edge>::iterator to = lv_grp->edges_to.begin(); to != lv_grp->edges_to.end(); ++to){
        lnodes.insert(to->node);
    }
    for (set<arg_node*>::iterator lnode = lnodes.begin(); lnode != lnodes.end(); ++lnode){
        for (vector<rc_edge>::iterator from = (*lnode)->edges_from.begin();
            from != (*lnode)->edges_from.end(); ++from){
            if (from->node == lv_grp && from->final_dest != NULL){
                failures.push_back(make_pair(*lnode, from->final_dest));
                failuretypes.push_back(from->type);
            }
        }
    }
    // Now check for this failure's compatibility with all preexisting failures.
    for (long unsigned int i = 0; i < failures.size(); ++i){
        int failuretype = failuretypes[i];
        
        if (type == 1){
            if (failuretype == 1){
                if (issubset_bitset(left_node->clade, failures[i].first->clade)){
                    short comp = compare_bitsets(left_node->clade, failures[i].second->clade);
                    if (comp != -1){
                        return false;
                    }
                }
                else if (issuperset_bitset(right_node->clade, failures[i].second->clade)){
                    short comp = compare_bitsets(failures[i].first->clade, right_node->clade);
                    if (comp != -1){
                        return false;
                    }
                }
            }
            else if (failuretype == 2){
                short comp = compare_bitsets(left_node->clade, failures[i].second->clade);
                if (comp != -1){
                    return false;
                }
            }
            else if (failuretype == 3){
                short comp = compare_bitsets(left_node->clade, failures[i].second->clade);
                if (comp != -1){
                    return false;
                }
            }
        }
        else if (type == 2){
            if (failuretype == 1){
                // check alpha/beta
                short comp = compare_bitsets(failures[i].first->clade, right_node->clade);
                if (comp != -1){
                    return false;
                }
                // check alpha/alpha
                if (issubset_bitset(left_node->clade, failures[i].first->clade)){
                    comp = compare_bitsets(left_node->clade, failures[i].second->clade);
                    if (comp != -1){
                        return false;
                    }
                }
            }
            else if (failuretype == 2){
                // check both alpha/betas
                short comp = compare_bitsets(left_node->clade, failures[i].second->clade);
                if (comp != -1){
                    return false;
                }
                comp = compare_bitsets(failures[i].first->clade, right_node->clade);
                if (comp != -1){
                    return false;
                }
            }
            else if (failuretype == 3){
                // check alpha/beta
                short comp = compare_bitsets(left_node->clade, failures[i].second->clade);
                if (comp != -1){
                    return false;
                }
                // check beta/beta
                if (issubset_bitset(right_node->clade, failures[i].second->clade)){
                    comp = compare_bitsets(failures[i].first->clade, right_node->clade);
                    if (comp != -1){
                        return false;
                    }
                }
            } 
        }
        else if (type == 3){
            if (failuretype == 1){
                short comp = compare_bitsets(failures[i].first->clade, right_node->clade);
                if (comp != -1){
                    return false;
                }
            }
            else if (failuretype == 2){
                short comp = compare_bitsets(failures[i].first->clade, right_node->clade);
                if (comp != -1){
                    return false;
                }
            }
            else if (failuretype == 3){
                if (issuperset_bitset(left_node->clade, failures[i].first->clade)){
                    short comp = compare_bitsets(left_node->clade, failures[i].second->clade);
                    if (comp != -1){
                        return false;
                    }
                }
                else if (issubset_bitset(right_node->clade, failures[i].second->clade)){
                    short comp = compare_bitsets(failures[i].first->clade, right_node->clade);
                    if (comp != -1){
                        return false;
                    }
                }
            }
        }
    }
    return true;
}

void merge_lv_nodes(arg_node* lv1, arg_node* lv2, arg_clademap& lv_grps){
    // Collapse ranges into one
    lv1->start = min(lv1->start, lv2->start);
    lv1->end = max(lv1->end, lv2->end);
    
    for (vector<rc_edge>::iterator to = lv2->edges_to.begin(); to != lv2->edges_to.end();){
        for (vector<rc_edge>::iterator from = to->node->edges_from.begin(); from != to->node->edges_from.end();
            ++from){
            // Edit exactly one edge to match this edge.
            if (from->node == lv2){
                from->node = lv1;
                break;
            }
        }
        lv1->edges_to.push_back(rc_edge(*to));
        lv2->edges_to.erase(to);
    }
    for (vector<rc_edge>::iterator from = lv2->edges_from.begin(); from != lv2->edges_from.end();){
        for (vector<rc_edge>::iterator to = from->node->edges_to.begin(); to != from->node->edges_to.end();
            ++to){
            // Edit exactly one edge to match this edge.
            if (to->node == lv2){
                to->node = lv1;
                break;
            }
        }
        lv1->edges_from.push_back(rc_edge(*from));
        lv2->edges_from.erase(from);
    }
    
    // Mark for deletion.
    pair<arg_clademap::iterator, arg_clademap::iterator> er = lv_grps.equal_range(lv2->clade);
    for (arg_clademap::iterator it = er.first; it != er.second;){
        if (it->second == lv2){
            lv_grps.erase(it++);
            //break;
        }
        else{
            ++it;
        }
    }
    
    delete lv2;
}

void merge_lv_nodes_recomb(set<arg_node*>& lv_grps_all,
    arg_clademap& lv_grps){
    // Merge any duplicate lv grps.
    bool lv_to_merge = true;
    while(lv_to_merge){
        arg_node* lvmerge1 = NULL;
        arg_node* lvmerge2 = NULL;
        for (set<arg_node*>::iterator lv1 = lv_grps_all.begin(); lv1 != lv_grps_all.end(); ++lv1){
            for (set<arg_node*>::iterator lv2 = lv_grps_all.begin(); lv2 != lv_grps_all.end(); ++lv2){
                if (lv2 != lv1 && (*lv1)->clade == (*lv2)->clade){
                    lvmerge1 = *lv1;
                    lvmerge2 = *lv2;
                    break;
                }
            }
            if (lvmerge1 != NULL && lvmerge2 != NULL){
                break;
            }
        }
        if (lvmerge1 != NULL && lvmerge2 != NULL){
            lv_grps_all.erase(lv_grps_all.find(lvmerge2));
            merge_lv_nodes(lvmerge1, lvmerge2, lv_grps);
        }
        else{
            lv_to_merge = false;
        }
    }
}

/**
 * Returns a pointer to the leaving group node.
 */
arg_node* add_lv_grp_graph(cladeset& clade,
    int type,
    arg_clademap& lv_grps,
    arg_node* left_node,
    arg_node* right_node,
    arg_sitemap& sites_pos,
    arg_clademap& sites_clade,
    long int prop_dist,
    int num_haplotypes){

    long int narrowest_left = *left_node->sites.rbegin();
    long int narrowest_right = *right_node->sites.begin();
    
    // Check to see if the leaving group node has already been built.
    arg_node* lv_grp = NULL;
    bool found = false;
    bool create = true;
    
    set<arg_node*> all_found;
    
    pair<arg_clademap::iterator, arg_clademap::iterator> er = lv_grps.equal_range(clade);
    for (arg_clademap::iterator node = er.first; node != er.second; ++node){
        
        // Don't allow overlap of 1 base
        if (min(narrowest_right, node->second->end) - max(narrowest_left, node->second->start) > 0){
            
            if (DEBUG_MODE){
                fprintf(stderr, "overlap? (%ld %ld) (%ld %ld)\n", narrowest_left,
                    narrowest_right, node->second->start, node->second->end);
            }
            
            found = true;
            
            bool compat = true;
            //bool compat = lv_grp_compatible(left_node, right_node, type, node->second);
            if (!compat){
                if (DEBUG_MODE){
                    fprintf(stderr, "NON-COMPAT\n");
                    print_bitset_set(clade, num_haplotypes);
                }
                found = false;
                create = false;
            }
            if (found){
                lv_grp = node->second;
                found = true;
                //break;
                all_found.insert(node->second);
            }
        }
    }

    if (!found){
        
        if (create){
            lv_grp = new arg_node;
            lv_grp->is_leaving_grp = true;
            lv_grp->start = narrowest_left;
            lv_grp->end = narrowest_right;
            lv_grp->site_support = false;
            lv_grp->clade = clade;
            lv_grps.insert(make_pair(clade, lv_grp));
        }
    }
    else{
        // Check to see if more than one match.
        
        if (all_found.size() > 1){
            merge_lv_nodes_recomb(all_found, lv_grps);
            lv_grp = *all_found.begin();
        }
        
        if (narrowest_left > lv_grp->start){
            lv_grp->start = narrowest_left;
        }
        if (narrowest_right < lv_grp->end){
            lv_grp->end = narrowest_right;
        }
        
    }
    
    if (lv_grp != NULL){
        // Add edges.

        rc_edge e1(type, lv_grp, right_node);
        left_node->edges_from.push_back(e1);
    
        rc_edge e2(type, left_node);
        lv_grp->edges_to.push_back(e2);
    
        rc_edge e3(type, right_node);
        lv_grp->edges_from.push_back(e3);
    
        rc_edge e4(type, lv_grp, left_node);
        right_node->edges_to.push_back(e4);

    }
    
    return lv_grp;
}

void del_connections_from_node(arg_node* node, 
    arg_node* lv, 
    bool erase_all,
    set<arg_node*>& lv_grps_delete){
    set<arg_node*> other_lv;

    set<arg_node*> dests_final;
    
    for (vector<rc_edge>::iterator edge_from = node->edges_from.begin();
        edge_from != node->edges_from.end();){
        if (edge_from->node == lv){
            dests_final.insert(edge_from->final_dest);
            if (erase_all){
                // Erase exactly one edge to match this edge.
                for (vector<rc_edge>::iterator edge_to = lv->edges_to.begin();
                    edge_to != lv->edges_to.end();){
                    if (edge_to->node == node && edge_from->type == edge_to->type){
                        lv->edges_to.erase(edge_to);
                        break;
                    }
                    else{
                        ++edge_to;
                    }
                }
            }
            node->edges_from.erase(edge_from);
        }
        else{
            ++edge_from;
        }
    }
    
    if (dests_final.size() > 0){
        for (vector<rc_edge>::iterator edge_from = node->edges_from.begin();
            edge_from != node->edges_from.end();){
            if (dests_final.find(edge_from->final_dest) != dests_final.end()){
                // Erase exactly one edge to match this edge.
                for (vector<rc_edge>::iterator edge_to = edge_from->node->edges_to.begin();
                    edge_to != edge_from->node->edges_to.end();){
                    if (edge_to->node == node && edge_to->type == edge_from->type){
                        edge_from->node->edges_to.erase(edge_to);
                        break;
                    }
                    else{
                        ++edge_to;
                    }
                }
                if (edge_from->node->edges_to.size() == 0 || edge_from->node->edges_from.size() == 0){
                    lv_grps_delete.insert(edge_from->node);
                }
                node->edges_from.erase(edge_from);
            }   
            else{
                ++edge_from;
            }
        }
    }    
    if (lv->edges_to.size() == 0 || lv->edges_from.size() == 0){
        lv_grps_delete.insert(lv);
    }
    
}

void del_connections_to_node(arg_node* node, 
    arg_node* lv,
    bool erase_all,
    set<arg_node*>& lv_grps_delete){
    set<arg_node*> other_lv;
    set<arg_node*> dests_final;
    
    for (vector<rc_edge>::iterator edge_to = node->edges_to.begin();
        edge_to != node->edges_to.end();){
        if (edge_to->node == lv){
            dests_final.insert(edge_to->final_dest);
            if (erase_all){
                // Erase exactly one edge to match this edge.
                for (vector<rc_edge>::iterator edge_from = lv->edges_from.begin();
                    edge_from != lv->edges_from.end();){
                    if (edge_from->node == node && edge_from->type == edge_to->type){
                        lv->edges_from.erase(edge_from);
                        break;
                    }
                    else{
                        ++edge_from;
                    }
                }
            }
            node->edges_to.erase(edge_to);
        }
        else{
            ++edge_to;
        }
    }
    
    if (dests_final.size() > 0){
        for (vector<rc_edge>::iterator edge_to = node->edges_to.begin();
            edge_to != node->edges_to.end();){
            if (dests_final.find(edge_to->final_dest) != dests_final.end()){
                // Erase exactly one edge to match this edge.
                for (vector<rc_edge>::iterator edge_from = edge_to->node->edges_from.begin();
                    edge_from != edge_to->node->edges_from.end();){
                    if (edge_from->node == node && edge_from->type == edge_to->type){
                        edge_to->node->edges_from.erase(edge_from);
                        break;
                    }
                    else{
                        ++edge_from;
                    }
                }
                if (edge_to->node->edges_from.size() == 0 || edge_to->node->edges_to.size() == 0){
                    lv_grps_delete.insert(edge_to->node);
                }
                node->edges_to.erase(edge_to);
            }   
            else{
                ++edge_to;
            }
        }
    }    
    if (lv->edges_to.size() == 0 || lv->edges_from.size() == 0){
        lv_grps_delete.insert(lv);
    }
}

void erase_lv_grps(set<arg_node*>& lv_grps_delete,
    arg_clademap& lv_grps,
    int num_haplotypes){
    set<arg_node*> already_deleted;
    
    while(lv_grps_delete.size() > 0){
        set<arg_node*> lv_grps_delete2;
        for (set<arg_node*>::iterator lv = lv_grps_delete.begin(); lv != lv_grps_delete.end();){
            
            bool found = false;
            pair<arg_clademap::iterator, arg_clademap::iterator> er = lv_grps.equal_range((*lv)->clade);
            for (arg_clademap::iterator it = er.first; it != er.second;){
                if (it->second == *lv){
                    lv_grps.erase(it++);
                    found = true;
                    //break;
                }
                else{
                    ++it;
                }
            }
            
            if (!found){
                fprintf(stderr, "OOPS; leaving group not found when deleting!\n");
                exit(1);
            }
            
            del_leaving_node(*lv, num_haplotypes);
            
            delete *lv;
            lv_grps_delete.erase(lv++);
        }
        lv_grps_delete = lv_grps_delete2;
    }
    
    
}

void del_connections_between(arg_node* main_left, arg_node* main_right,
    arg_clademap& lv_grps, set<arg_node*>& lv_grps_set, int num_haplotypes){
    main_left->fixed_right = false;
    main_right->fixed_left = false;
    
    set<arg_node*> lv_del;
    
    for (vector<rc_edge>::iterator from = main_left->edges_from.begin();
        from != main_left->edges_from.end();){
        bool remove = false;
        if (from->final_dest != NULL && from->final_dest == main_right){
            bool found = false;
            for (vector<rc_edge>::iterator from2 = from->node->edges_from.begin();
                from2 != from->node->edges_from.end();){
                if (from2->node == main_right && from2->type == from->type){
                    from->node->edges_from.erase(from2);
                    found = true;
                    break;
                }
                else{
                    ++from2;
                }
            }
            if (!found){
                fprintf(stderr, "!found1lv\n");
                exit(1);
            }
            remove = true;
        }
        if (remove){
            main_left->edges_from.erase(from);
        }
        else{
            ++from;
        }
    }
    for (vector<rc_edge>::iterator to = main_right->edges_to.begin();
        to != main_right->edges_to.end();){
        bool remove = false;
        if (to->final_dest != NULL && to->final_dest == main_left){
            bool found = false;
            for (vector<rc_edge>::iterator to2 = to->node->edges_to.begin();
                to2 != to->node->edges_to.end();){
                if (to2->node == main_left && to2-> type == to->type){
                    to->node->edges_to.erase(to2);
                    found = true;
                    break;
                }
                else{
                    ++to2;
                }
            }
            if (to->node->edges_from.size() == 0 || to->node->edges_to.size() == 0){
                lv_del.insert(to->node);
            }
            if (!found){
                fprintf(stderr, "!found2lv\n");
                exit(1);
            }
            remove = true;
        }
        if (remove){
            main_right->edges_to.erase(to);
        }
        else{
            ++to;
        }
    }
    
    for (set<arg_node*>::iterator lvd = lv_del.begin(); lvd != lv_del.end(); ++lvd){
        if (lv_grps_set.find(*lvd) != lv_grps_set.end()){
            lv_grps_set.erase(lv_grps_set.find(*lvd));
        }
    }
    
    erase_lv_grps(lv_del, lv_grps, num_haplotypes);
}

void del_solved_connections_between(arg_node* main_left, arg_node* main_right){
    main_left->fixed_right = false;
    main_right->fixed_left = false;
    
    for (vector<rc_edge>::iterator from = main_left->edges_from_solved.begin();
        from != main_left->edges_from_solved.end();){
        bool remove = false;
        if (from->final_dest != NULL && from->final_dest == main_right){
            bool found = false;
            for (vector<rc_edge>::iterator from2 = from->node->edges_from_solved.begin();
                from2 != from->node->edges_from_solved.end();){
                if (from2->node == main_right && from2->type == from->type){
                    from->node->edges_from_solved.erase(from2);
                    found = true;
                    break;
                }
                else{
                    ++from2;
                }
            }
            if (!found){
                fprintf(stderr, "!found1lvS\n");
                exit(1);
            }
            remove = true;
        }
        if (remove){
            main_left->edges_from_solved.erase(from);
        }
        else{
            ++from;
        }
    }
    for (vector<rc_edge>::iterator to = main_right->edges_to_solved.begin();
        to != main_right->edges_to_solved.end();){
        bool remove = false;
        if (to->final_dest != NULL && to->final_dest == main_left){
            bool found = false;
            for (vector<rc_edge>::iterator to2 = to->node->edges_to_solved.begin();
                to2 != to->node->edges_to_solved.end();){
                if (to2->node == main_left && to2-> type == to->type){
                    to->node->edges_to_solved.erase(to2);
                    found = true;
                    break;
                }
                else{
                    ++to2;
                }
            }
            if (!found){
                fprintf(stderr, "!found2lv\n");
                exit(1);
            }
            remove = true;
        }
        if (remove){
            main_right->edges_to_solved.erase(to);
        }
        else{
            ++to;
        }
    }
}

bool path_exists_between(arg_node* left, arg_node* right){
    for (vector<rc_edge>::iterator from = left->edges_from_solved.begin(); from != left->edges_from_solved.end();
        ++from){
        if (from->final_dest != NULL && from->final_dest == right){
            return true;
        }
    }
    // If not directly connected, search recursively.
    for (vector<rc_edge>::iterator from = left->edges_from_solved.begin(); from != left->edges_from_solved.end();
        ++from){
        if (from->final_dest != NULL && path_exists_between(from->final_dest, right)){
            return true;
        }
    }
    return false;
}

bool l_compat_r(arg_node* left, arg_node* right){
    if (min(left->end, right->end) - max(left->start, right->start) > 0){
        return true;
    }
    if (left->recomb_right == -1 || left->recomb_right > *right->sites.begin()){
        if (right->recomb_left == -1 || right->recomb_left < *left->sites.rbegin()){
            short comp = compare_bitsets(left->clade, right->clade);
            if (comp == -1){
                return false;
            }
            else{
                return true;
            }
        }
        else{
            return false;
        }
    }
    else{
        return false;
    }
}

bool l_compat_r_conservative(arg_node* left, arg_node* right, bool recomb_left){
    if (min(left->end, right->end) - max(left->start, right->start) > 0){
        return true;
    }
    if (left->recomb_right == -1 || left->recomb_right > *right->sites.begin()){
        if (right->recomb_left == -1 || right->recomb_left < *left->sites.rbegin()){
            short comp = compare_bitsets(left->clade, right->clade);
            if (comp == -1){
                return false;
            }
            else{
                if (recomb_left && left->recomb_right != -1 && left->recomb_right <= *right->sites.rbegin()){
                    return false;
                }
                else if (!recomb_left && right->recomb_left != -1 && right->recomb_left >= *left->sites.begin()){
                    return false;
                }
                else{
                    return true;
                }
            }
        }
        else{
            return false;
        }
    }
    else{
        return false;
    }
}

void set_4haptest_failure(arg_node* left_node, 
    arg_node* right_node,
    long int left_site,
    long int right_site,
    arg_sitemap& sites_pos,
    arg_clademap& sites_clade,
    arg_clademap& lv_grps,
    long int prop_dist,
    arg_node* root,
    vector<arg_node*>& newnodes,
    int num_haplotypes){

    if (DEBUG_MODE){
        fprintf(stderr, "s4ht %ld %ld\n", left_site, right_site);
        print_node_lite(left_node, num_haplotypes);
        print_node_lite(right_node, num_haplotypes);
    }
    
    if (right_site - left_site > 2*prop_dist){
        if (DEBUG_MODE){
            fprintf(stderr, "too far apart!\n");
        }
        //exit(1);
        return;
    }
    if (DEBUG_MODE && right_site <= left_site){
        fprintf(stderr, "impossible right/left site combo %ld %ld\n", left_site, right_site);
        exit(1);
    }
    
    if (left_node->recomb_right == -1 || *(right_node->sites.begin()) < left_node->recomb_right){
        left_node->recomb_right = *(right_node->sites.begin());
        //unfix_edges_right(left_node);
    }
    if (right_node->recomb_left == -1 || *(left_node->sites.rbegin()) > right_node->recomb_left){
        right_node->recomb_left = *(left_node->sites.rbegin());
        //unfix_edges_left(right_node);
    }
    
    if (left_node->recomb_right < left_node->end || right_node->recomb_left > right_node->start){
        fprintf(stderr, "recomb indices messup1\n");
        print_node_lite(left_node, num_haplotypes);
        print_node_lite(right_node, num_haplotypes);
        exit(1);
    }
    
    if (DEBUG_MODE && left_node->clade == right_node->clade){
        fprintf(stderr, "Four haplotype test failure being set between identical clades!\n");
        print_node_lite(left_node, num_haplotypes);
        print_node_lite(right_node, num_haplotypes);
        exit(1);
    }
    
    // If there are already recombination edges connecting these two nodes, then do nothing.
    if (left_node->edges_from.size() <= right_node->edges_to.size()){
        for (vector<rc_edge>::iterator edge_from = left_node->edges_from.begin(); edge_from != left_node->edges_from.end();
            ++edge_from){
            if (edge_from->final_dest == right_node){
                return;
            }
        }
    }
    else{
        for (vector<rc_edge>::iterator to = right_node->edges_to.begin(); to != right_node->edges_to.end(); ++to){
            if (to->final_dest == left_node){
                return;
            }
        }
    }
    if (left_node->edges_from_solved.size() <= right_node->edges_to_solved.size()){
        for (vector<rc_edge>::iterator edge_from_solved = left_node->edges_from_solved.begin();
            edge_from_solved != left_node->edges_from_solved.end(); ++edge_from_solved){
            if (edge_from_solved->final_dest == right_node){
                return;
            }
        }
    }
    else{
        for (vector<rc_edge>::iterator to = right_node->edges_to_solved.begin();
            to != right_node->edges_to_solved.end(); ++to){
            if (to->final_dest == left_node){
                return;
            }
        }
    }
    if (left_node->edges_from_unsolvable.size() <= right_node->edges_to_unsolvable.size()){
        for (vector<rc_edge>::iterator from = left_node->edges_from_unsolvable.begin(); 
            from != left_node->edges_from_unsolvable.end(); ++from){
            if (from->final_dest == right_node){
                return;
            }
        }
    }
    else{
        for (vector<rc_edge>::iterator to = right_node->edges_to_unsolvable.begin();
            to != right_node->edges_to_unsolvable.end(); ++to){
            if (to->final_dest == left_node){
                return;
            }
        }
    }
    
    // Determine which 3 clades could represent possible
    // leaving groups. For each one that does not fail the
    // four haplotype test within its range, create a leaving
    // node joining the left & right clade.
    
    bool already_solved = false;

    arg_node* moved = NULL;
    bool up = false;
    bool down = false;
    bool left_alpha = false;
    bool right_alpha = false;
    bool left_beta = false;
    bool right_beta = false;
    bool already_connected = false;
    bool l_fd_before = false;
    bool r_fd_before = false;
    
    bool left_new = false;
    
    set<arg_node*> lv_between;
    
    for (vector<rc_edge>::iterator lr = left_node->edges_from_solved.begin(); lr != left_node->edges_from_solved.end();
        ++lr){
        if (lr->final_dest != NULL){
            if (lr->final_dest == right_node){
                already_connected = true;
                break;
            }
            if (*lr->final_dest->sites.begin() <= *right_node->sites.begin()){
                
                if (compare_bitsets(lr->final_dest->clade, right_node->clade) == -1){
                    l_fd_before = true;
                }
                
                lv_between.insert(lr->node);
            
                if (lr->type == 1){
                    // Determine if right node is alpha (if beta, then this is a completely
                    // different recombination event).
                    if (lr->node->clade == set_diff_bitset(left_node->clade, right_node->clade)){
                        // right node is alpha
                        already_solved = true;
                        moved = lr->node;
                        right_alpha = true;
                        up = true;
                        left_alpha = true;
                        //break;
                    }
                
                }
                else if (lr->type == 2){
                    // Determine if right node is alpha or beta.
                    short comp = compare_bitsets(right_node->clade, lr->final_dest->clade);
                    if (comp == 1 || comp == 3){
                        // right node is beta
                        already_solved = true;
                        moved = lr->node;
                        right_beta = true;
                        left_alpha = true;
                        //break;
                    }
                    else{
                        cladeset testclade = right_node->clade | lr->node->clade;
                        if (issuperset_bitset(testclade, left_node->clade)){
                            // right node is alpha
                            already_solved = true;
                            right_alpha = true;
                            moved = lr->node;
                            right_alpha = true;
                            left_alpha = true;
                            //break;
                        }
                    }
                
                }
                else if (lr->type == 3){
                    // Here, the left node is the superset of everything (leaving clade
                    // moved down). If it fails the four haplotype test with right node,
                    // it's not a beta/beta.
                }
            }
        }
    }

    if (!already_connected){
    //if (!already_solved){
        for (vector<rc_edge>::iterator rl = right_node->edges_to_solved.begin(); rl !=
            right_node->edges_to_solved.end(); ++rl){
            if (rl->final_dest != NULL){
                if (rl->final_dest == left_node){
                    already_connected = true;
                    break;
                }
                if (*rl->final_dest->sites.rbegin() >= *left_node->sites.rbegin()){
                    
                    if (compare_bitsets(rl->final_dest->clade, left_node->clade) == -1){
                        r_fd_before = true;
                    }
                    
                    lv_between.insert(rl->node);
                    
                    if (rl->type == 1){
                        // Here, the right node is the superset of everything (leaving clade
                        // moved up). If it fails the four haplotype test with left node,
                        // it's not alpha/alpha.
                    }
                    else if (rl->type == 2){
                        // Determine if left node is alpha or beta.
                        short comp = compare_bitsets(left_node->clade, rl->final_dest->clade);
                        if (comp == 1 || comp == 3){
                            already_solved = true;
                            moved = rl->node;
                            left_alpha = true;
                            //break;
                            right_beta = true;
                            left_new = true;
                        }
                        else{
                            cladeset testclade = left_node->clade | rl->node->clade;
                            if (issuperset_bitset(testclade, right_node->clade)){
                                // left node is beta
                                already_solved = true;
                                moved = rl->node;
                                left_beta = true;
                                //break;
                                right_beta = true;
                                left_new = true;
                            }
                        }
                    }
                    else if (rl->type == 3){
                        // Determine if left node is beta (if alpha, then this is a 
                        // completely different recombination event).
                        if (rl->node->clade == set_diff_bitset(right_node->clade, left_node->clade)){
                            already_solved = true;
                            moved = rl->node;
                            left_beta = true;
                            down = true;
                            //break;
                            right_beta = true;
                            left_new = true;
                        }
                    }
                }
            }
        }
    }

    if (already_connected){
        return;
    }
    if (l_fd_before && r_fd_before){
        
        for (set<arg_node*>::iterator lv = lv_between.begin(); lv != lv_between.end();
            ++lv){
            cladeset lprime = 
                set_diff_bitset(left_node->clade, right_node->clade);
            cladeset rprime = 
                set_diff_bitset(right_node->clade, left_node->clade);
            cladeset dd = 
                left_node->clade & right_node->clade;
            
            int type = -1;
            if ((*lv)->clade == lprime){
                type = 1;
            }
            else if ((*lv)->clade == rprime){
                type = 3;
            }
            else if ((*lv)->clade == dd){
                type = 2;
            }
            if (type != -1){
                add_solved_edges(left_node, right_node, *lv, type, false, newnodes, prop_dist, num_haplotypes);
            }
            else{
                // Add unsolvable edges.
                add_unsolved_edges(left_node, right_node, num_haplotypes);
            }
        }
    
        return;
    }
    
    if (already_solved){
        // Adjust start/end limits as necessary.
        
        // Add recombination edges as necessary. Don't need to add edges if recomb is
        // neither up nor down and left_beta or right_alpha.
        if (true){
        //if (up || down || (!up && !down && !left_beta && !right_alpha)){
            cladeset lprime = 
                set_diff_bitset(left_node->clade, right_node->clade);
            cladeset rprime = 
                set_diff_bitset(right_node->clade, left_node->clade);
            cladeset dd = 
                left_node->clade & right_node->clade;
                
            int type = -1;
            if (up && moved->clade == lprime){
                type = 1;
            }
            else if (down && moved->clade == rprime){
                type = 3;
            }
            else if (left_alpha && right_alpha && moved->clade == lprime){
                type = 1;
            }
            else if (left_alpha && right_beta && moved->clade == dd){
                type = 2;
            }
            else if (left_beta && right_beta && moved->clade == rprime){
                type = 3;
            }
            if (type != -1){
 
                add_solved_edges(left_node, right_node, moved, type, false, newnodes, prop_dist, num_haplotypes);
            }
            else{
                add_unsolved_edges(left_node, right_node, num_haplotypes);
            }
            
        }
        
        
        // Expand range to limits.
        if ((left_new) && !l_fd_before && (left_node->recomb_right == -1 || left_node->recomb_right >= *right_node->sites.begin())){
            // Need to expand left node rightward.
            long int newend = right_node->recomb_left;
            set<arg_node*> dummy;
            expand_range_fixed(left_node, newend, sites_pos, sites_clade, lv_grps,
                prop_dist, root, newnodes, false, num_haplotypes, true, dummy);
            
        }
        else if ((!left_new) && !r_fd_before && (right_node->recomb_left == -1 || right_node->recomb_left <= *left_node->sites.rbegin())){
            // Need to expand right node leftward.
            long int newstart = left_node->recomb_right;
            set<arg_node*> dummy;
            expand_range_fixed(right_node, newstart, sites_pos, sites_clade, lv_grps,
                prop_dist, root, newnodes, false, num_haplotypes, true, dummy);
            
        }
        
        // Check to see if we should go ahead and create new left and/or right nodes.
        
        bool cr = left_node->recomb_right >= *right_node->sites.begin();
        bool cl = right_node->recomb_left <= *left_node->sites.begin();
        
        // Create new nodes as necessary.
        //if (false){
        if (left_new && cr && left_alpha && *left_node->sites.rbegin() + prop_dist >= *right_node->sites.begin() &&
            (left_node->recomb_right == -1 || left_node->recomb_right >= *right_node->sites.begin())){
            
            cladeset clade = set_diff_bitset(left_node->clade, moved->clade);
            
            if (safe_to_create_node(clade, sites_pos, sites_clade, 
                *right_node->sites.begin(), prop_dist, num_haplotypes) &&
                left_node->closest_mut_l() != -1 &&
                left_node->closest_mut_l() + prop_dist >= *right_node->sites.begin()){

                arg_node* newnode = new arg_node;
                newnode->is_leaving_grp = false;
                newnode->site_support = false;
                newnode->recomb_left = -1;
                newnode->recomb_right = -1;
                newnode->skip = false;
                newnode->clade = clade;
                newnode->sites.insert(*right_node->sites.begin());
                newnode->start = *right_node->sites.begin();
                newnode->end = *left_node->sites.rbegin() + prop_dist;
                newnode->fixed_left = false;
                newnode->fixed_right = false;
                newnode->closest_mut_left = left_node->closest_mut_l();
                newnodes.push_back(newnode);  
            }
        }
        else if (left_new && cr && left_beta && *left_node->sites.rbegin() + prop_dist >= *right_node->sites.begin() &&
            (left_node->recomb_right == -1 || left_node->recomb_right >= *right_node->sites.begin())){
            cladeset clade = left_node->clade | moved->clade;
            if (safe_to_create_node(clade, sites_pos, sites_clade, *right_node->sites.begin(),
                prop_dist, num_haplotypes) &&
                left_node->closest_mut_l() != -1 && 
                left_node->closest_mut_l() + prop_dist >= *right_node->sites.begin()){

                arg_node* newnode = new arg_node;
                newnode->is_leaving_grp = false;
                newnode->site_support = false;
                newnode->recomb_left = -1;
                newnode->recomb_right = -1;
                newnode->skip = false;
                newnode->clade = clade;
                newnode->sites.insert(*right_node->sites.begin());
                newnode->start = *right_node->sites.begin();
                newnode->end = *left_node->sites.rbegin() + prop_dist;
                newnode->fixed_left = false;
                newnode->fixed_right = false;
                newnode->closest_mut_left = left_node->closest_mut_l();
                newnodes.push_back(newnode);
            }
        }
        else if (!left_new && cl && right_alpha && *right_node->sites.begin() - prop_dist <= *left_node->sites.rbegin() &&
            (right_node->recomb_left == -1 || right_node->recomb_left <= *left_node->sites.rbegin())){
            cladeset clade = right_node->clade | moved->clade;
            if (safe_to_create_node(clade, sites_pos, sites_clade, *left_node->sites.rbegin(),
                prop_dist, num_haplotypes) &&
                right_node->closest_mut_r() != -1 && 
                right_node->closest_mut_r() - prop_dist <= *left_node->sites.rbegin()){

                arg_node* newnode = new arg_node;
                newnode->is_leaving_grp = false;
                newnode->site_support = false;
                newnode->recomb_left = -1;
                newnode->recomb_right = -1;
                newnode->skip = false;
                newnode->clade = clade;
                newnode->sites.insert(*left_node->sites.rbegin());
                newnode->start = max((long int)1, *right_node->sites.begin()-prop_dist);
                newnode->end = *left_node->sites.rbegin();
                newnode->fixed_left = false;
                newnode->fixed_right = false;
                newnode->closest_mut_right = right_node->closest_mut_r();
                newnodes.push_back(newnode);
            }
        }
        else if (!left_new && cl && right_beta && *right_node->sites.begin()-prop_dist <= *left_node->sites.rbegin() &&
            (right_node->recomb_left == -1 || right_node->recomb_left <= *left_node->sites.rbegin())){
            
            cladeset clade = set_diff_bitset(right_node->clade, moved->clade);
            
            if (safe_to_create_node(clade, sites_pos, sites_clade,
                *left_node->sites.rbegin(), prop_dist, num_haplotypes) &&
                right_node->closest_mut_r() != -1 &&
                right_node->closest_mut_r() - prop_dist <= *left_node->sites.rbegin()){

                arg_node* newnode = new arg_node;
                newnode->is_leaving_grp = false;
                newnode->site_support = false;
                newnode->recomb_left = -1;
                newnode->recomb_right = -1;
                newnode->skip = false;
                newnode->clade = clade;
                newnode->sites.insert(*left_node->sites.rbegin());
                newnode->start = max((long int)1, *right_node->sites.begin()-prop_dist);
                newnode->end = *left_node->sites.rbegin();
                newnode->fixed_left = false;
                newnode->fixed_right = false;
                newnode->closest_mut_right = right_node->closest_mut_r();
                
                newnodes.push_back(newnode);
            }
        }
        return;
    }
    
    cladeset leftprime = set_diff_bitset(left_node->clade, right_node->clade);
    cladeset rightprime = set_diff_bitset(right_node->clade, left_node->clade);
    cladeset dd = set_int_bitset(left_node->clade, right_node->clade);
    
    bool leftprime_valid = true;
    bool rightprime_valid = true;
    bool dd_valid = true;
    
    
    // Check to see if a previous recombination event has been added to the graph that
    // shares a leaving group with this recombination event, but whose left or
    // Build graph nodes.

    
    int num_edges = 0;
    
    if (leftprime_valid){
        add_lv_grp_graph(leftprime, 1, lv_grps,
            left_node, right_node, sites_pos, sites_clade, prop_dist, num_haplotypes);
        num_edges++;
    }
    if (rightprime_valid){
        
        add_lv_grp_graph(rightprime, 3, lv_grps,
            left_node, right_node, sites_pos, sites_clade, prop_dist, num_haplotypes);
        num_edges++;  
    }
    if (dd_valid){
        add_lv_grp_graph(dd, 2, lv_grps, 
            left_node, right_node, sites_pos, sites_clade, prop_dist, num_haplotypes);
        num_edges++;
    }
    
    if (DEBUG_MODE){
        for (arg_clademap::iterator l = lv_grps.begin();
            l != lv_grps.end(); ++l){
            if (l->second->edges_from.size() == 0 || l->second->edges_to.size() == 0){
                fprintf(stderr, "s4ht end: leaving group has empty edges!\n");
                print_node_lite(l->second, num_haplotypes);
                exit(1);
            }
        }
        
        if (left_node->end < left_node->start){
            fprintf(stderr, "left node end < start\n");
            print_node_lite(left_node, num_haplotypes);
            exit(1);
        }
        if (right_node->end < right_node->start){
            fprintf(stderr, "right node end < start\n");
            print_node_lite(right_node, num_haplotypes);
            exit(1);
        }
    }
}

/**
 * Eliminates second node.
 */
void join_nodes(arg_node* node1,
    arg_node* node2,
    arg_sitemap& sites_pos,
    arg_clademap& sites_clade,
    arg_clademap& lv_grps,
    long int prop_dist,
    bool exp_range,
    arg_node* root,
    vector<arg_node*>& newnodes,
    int num_haplotypes){

    if (DEBUG_MODE){
        if (node1->clade != node2->clade){
            fprintf(stderr, "joining nodes with different clades\n");
            print_node_lite(node1, num_haplotypes);
            print_node_lite(node2, num_haplotypes);
            exit(1);
        }
        
        if (node1 == root || node2 == root){
            fprintf(stderr, "joining to root\n");
            print_node_lite(node1, num_haplotypes);
            print_node_lite(node2, num_haplotypes);
            exit(1);
        }
        
        fprintf(stderr, "JOIN NODES\n");
        
        print_node_upward_onelevel(node1, 0, num_haplotypes);
        print_node_upward_onelevel(node2, 0, num_haplotypes);
    }
    
    long int node2start = node2->start;
    long int node2end = node2->end;
    
    long int node2sitestart = *(node2->sites.begin());
    long int node2siteend = *(node2->sites.rbegin());
     
    // Handle recombination stuff.
    transfer_edges_from(node2, node1, num_haplotypes);
    transfer_edges_to(node2, node1, num_haplotypes);

    if (node2->recomb_left != -1 && (node1->recomb_left == -1 || node2->recomb_left > node1->recomb_left)){
        node1->recomb_left = node2->recomb_left;
    }
    if (node2->recomb_right != -1 && (node1->recomb_right == -1 || node2->recomb_right < node1->recomb_right)){
        node1->recomb_right = node2->recomb_right;
    }
    
    if ((node1->recomb_left != -1 && node1->recomb_left > node1->start) || 
        (node1->recomb_right != -1 && node1->recomb_right < node1->end)){
        fprintf(stderr, "recomb indices messup2\n");
        print_node_lite(node1, num_haplotypes);
        print_node_lite(node2, num_haplotypes);
        exit(1);
    }

    for (set<long int>::iterator site = node2->sites.begin(); site != node2->sites.end();){
        pair<arg_sitemap::iterator, arg_sitemap::iterator> er = sites_pos.equal_range(*site);
        
        for (arg_sitemap::iterator it = er.first; it != er.second; ){
            if (it->second == node2){
                sites_pos.erase(it++);
                break;
            }
            else{
                ++it;
            }
        }
        
        if (node1->sites.find(*site) == node1->sites.end()){
            node1->sites.insert(*site);
            bool found = false;
            er = sites_pos.equal_range(*site);
            for (arg_sitemap::iterator it = er.first;
                it != er.second; ++it){
                if (it->second == node1){
                    found = true;
                    break;
                }
            }
            if (!found){
                sites_pos.insert(make_pair(*site, node1));
            }
        }
        
        node2->sites.erase(site++);
    }

    for (set<long int>::iterator mut = node2->mutations.begin(); mut != node2->mutations.end();){
        node1->mutations.insert(*mut);
        node2->mutations.erase(mut++);
    }
    for (set<long int>::iterator site = node2->sites_lvgrp.begin(); site != node2->sites_lvgrp.end();){
        node1->sites_lvgrp.insert(*site);
        node2->sites_lvgrp.erase(site++);
    }

    // Wipe out node2 connections to parents & children
    set<arg_node*, nodeptr_parent_sort> parents;
    set<arg_node*, nodeptr_child_sort> children;
    
    for (set<vert_edge>::iterator parent = node2->parents.begin(); parent != node2->parents.end();){
        for (set<vert_edge>::iterator child = parent->to->children.begin(); child != parent->to->children.end();){
            if (parent->start == child->start && parent->end == child->end && child->to == node2){
                parent->to->children.erase(child++);
                break;
            }
            else{
                ++child;
            }
        }
        parents.insert(parent->to);
        node2->parents.erase(parent++);
    }
    for (set<vert_edge>::iterator child = node2->children.begin(); child != node2->children.end();){
        for (set<vert_edge>::iterator parent = child->to->parents.begin(); parent != child->to->parents.end();){
            if (parent->start == child->start && parent->end == child->end && parent->to == node2){
                child->to->parents.erase(parent++);
                break;
            }
            else{
                ++parent;
            }
        }
        children.insert(child->to);
        node2->children.erase(child++);
    }
    
    node1->start = min(node1->start, node2->start);
    node1->end = max(node1->end, node2->end);
    
    for (set<arg_node*, nodeptr_parent_sort>::iterator p = parents.begin(); p != parents.end(); ++p){
        set_parents(node1, *p, true, prop_dist, sites_pos, sites_clade, lv_grps, children, num_haplotypes);
    }
    for (set<arg_node*, nodeptr_child_sort>::iterator c = children.begin(); c != children.end(); ++c){
        set_parents(*c, node1, true, prop_dist, sites_pos, sites_clade, lv_grps, children, num_haplotypes);

    }
    return;
    
    
    long int prevstart = -1;
    long int prevend = -1;
    
    if (node2start < node1->start && (node1->recomb_left == -1 || node1->recomb_left < node2start) &&
        !(node1->start == *(node1->sites.begin()) && node2start < node2sitestart)){
        prevstart = node1->start;
        node1->start = node2start;
    }
    else if (node2sitestart < node1->start && (node1->recomb_left == -1 || node1->recomb_left < node2sitestart)){
        prevstart = node1->start;
        node1->start = node2sitestart;
    }

    if (node2end > node1->end && (node1->recomb_right == -1 || node1->recomb_right > node2end) &&
        !(node1->end == *(node1->sites.rbegin()) && node2end > node2siteend)){
        prevend = node1->end;
        node1->end = node2end;
    }
    else if (node2siteend > node1->end && (node1->recomb_right == -1 || node1->recomb_right > node2siteend)){
        prevend = node1->end;
        node1->end = node2siteend;
    }
    
    // Figure out what new sites need to be checked
    set<arg_node*> failures;
    
    if (prevstart != -1 && node1->start < prevstart){

        if (DEBUG_MODE){
            fprintf(stderr, "haf/gp3 %ld %ld\n", node1->start, prevstart-1);
        }
        compile_failures(node1, root, failures, node1->start, prevstart-1, prop_dist, sites_pos, num_haplotypes, false, false);

        adjust_all_failure_indices(node1, failures, sites_pos, sites_clade, 
            lv_grps, prop_dist, root, newnodes, num_haplotypes);
        
        // Gathering failures might change the start coordinate and make adding new parents
        // unnecessary.
        if (node1->start < prevstart){
            if (DEBUG_MODE){
                fprintf(stderr, "gp3 %ld %ld\n", node1->start, prevstart-1);
            }
            gather_parents(node1, root, node1->start, prevstart-1, prop_dist, 
                sites_pos, sites_clade, lv_grps, parents, num_haplotypes);
        }
    }

    if (prevend != -1 && node1->end > prevend){

        if (DEBUG_MODE){
            fprintf(stderr, "haf/gp4 %ld %ld\n", prevend+1, node1->end);
        }
        compile_failures(node1, root, failures, prevend+1, node1->end, prop_dist, sites_pos, num_haplotypes, false, false);

        adjust_all_failure_indices(node1, failures, sites_pos, sites_clade, 
            lv_grps, prop_dist, root, newnodes, num_haplotypes);
        
        // Gathering failures might change the end coordinate and make adding new parents 
        // unnecessary.
        if (node1->end > prevend){
            if (DEBUG_MODE){
                fprintf(stderr, "gp4 %ld %ld\n", prevend+1, node1->end);
            }
            gather_parents(node1, root, prevend+1, node1->end, prop_dist, 
                sites_pos, sites_clade, lv_grps, parents, num_haplotypes);
        }
    }
    
    for (set<vert_edge>::iterator edge = node2->parents.begin();
        edge != node2->parents.end();){
        for (set<vert_edge>::iterator e2 = edge->to->children.begin();
            e2 != edge->to->children.end();){
            if (e2->to == node2){
                edge->to->children.erase(e2++);
            }
            else{
                ++e2;
            }
        }
        if (DEBUG_MODE){
            fprintf(stderr, "erased parent <%ld %ld>\n", edge->start, edge->end);
            print_node_lite(edge->to, num_haplotypes);
        }
        parents.insert(edge->to);
        node2->parents.erase(edge++);
    }

    for (set<vert_edge>::iterator edge = node1->parents.begin();
        edge != node1->parents.end();){
        for (set<vert_edge>::iterator e2 = edge->to->children.begin();
            e2 != edge->to->children.end();){
            if (e2->to == node1){
                edge->to->children.erase(e2++);
            }
            else{
                ++e2;
            }
        }
        if (DEBUG_MODE){
            fprintf(stderr, "erased parent <%ld %ld>\n", edge->start, edge->end);
            print_node_lite(edge->to, num_haplotypes);
        }
        parents.insert(edge->to);
        node1->parents.erase(edge++);
    }
    
    for (set<vert_edge>::iterator edge = node2->children.begin();
        edge != node2->children.end();){
        for (set<vert_edge>::iterator e2 = edge->to->parents.begin();
            e2 != edge->to->parents.end();){
            if (e2->to == node2){
                edge->to->parents.erase(e2++);
            }
            else{
                ++e2;
            }
        }   
        children.insert(edge->to);
        node2->children.erase(edge++);
    }
 
    for (set<vert_edge>::iterator edge = node1->children.begin();
        edge != node1->children.end();){
        for (set<vert_edge>::iterator e2 = edge->to->parents.begin();
            e2 != edge->to->parents.end();){
            if (e2->to == node1){
                edge->to->parents.erase(e2++);
            }
            else{
                ++e2;
            }
        }
        children.insert(edge->to);
        node1->children.erase(edge++);
    }

    for (set<arg_node*, nodeptr_parent_sort>::iterator parent = parents.begin(); parent != parents.end(); ++parent){
        if (*parent != root){
            set_parents(node1, *parent, true, prop_dist, sites_pos, sites_clade, lv_grps, children, num_haplotypes);
        }
    }
    
    set<arg_node*, nodeptr_child_sort> dummy;
    
    for (set<arg_node*, nodeptr_child_sort>::iterator child = children.begin(); child != children.end(); ++child){
        set_parents(*child, node1, false, prop_dist, sites_pos, sites_clade, lv_grps, dummy, num_haplotypes);
        if (DEBUG_MODE && dummy.size() > 0){
            fprintf(stderr, "??? dummy > 0 (5)\n");
            exit(1);
        }
    }
    handle_failures(node1, failures, sites_pos, sites_clade, lv_grps, prop_dist, 
        root, newnodes, num_haplotypes);
    

    if (DEBUG_MODE && (node1->start > *(node1->sites.begin()) || 
        node1->end < *(node1->sites.rbegin()))){
        fprintf(stderr, "node1 range is impossible\n");
        exit(1);
    }
    
    if (DEBUG_MODE){
        fprintf(stderr, "JOIN NODES DONE\n");
        print_node_downward_onelevel(node1, 0, num_haplotypes);
        
        for (set<vert_edge>::iterator c = node1->children.begin(); c != node1->children.end(); ++c){
            check_for_holes(c->to, root, num_haplotypes);
        }
    }
    
    node1->fixed_left = false;
    node1->fixed_right = false;
}

/**
 * Returns pointer to new node created (to the right)
 */
arg_node* split_node(arg_node* node, 
    long int split_pos, 
    arg_sitemap& sites_pos,
    arg_clademap& sites_clade,
    arg_clademap& lv_grps,
    long int prop_dist,
    arg_node* root,
    vector<arg_node*>& newnodes,
    int num_haplotypes){
    if (DEBUG_MODE){
        fprintf(stderr, "SPLIT NODE\n");
    }
    
    long int oldnodeend = node->end;
    
    set<long int> sites_right;
    set<long int> muts_right;
    for (set<long int>::iterator site = node->sites.begin(); site != node->sites.end();){
        if (*site < split_pos){
            ++site;
        }
        else if (*site > split_pos){
            pair<arg_sitemap::iterator, arg_sitemap::iterator> er = sites_pos.equal_range(*site);
            for (arg_sitemap::iterator it = er.first; it != er.second; ){
                if (it->second == node){
                    sites_pos.erase(it++);
                }
                else{
                    ++it;
                }
            }
            sites_right.insert(*site);
            if (node->mutations.find(*site) != node->mutations.end()){
                node->mutations.erase(node->mutations.find(*site));
                muts_right.insert(*site);
            }
            node->sites.erase(site++);  
        }
    }
    
    long int prev_end = node->end;
    node->end = *(node->sites.rbegin());

    // Shorten the range of the existing node.
    shorten_range(node, node->end - prev_end, sites_pos, sites_clade, lv_grps, prop_dist, 
        num_haplotypes, root, newnodes, -1);
    
    node->recomb_right = split_pos;
    
    node->fixed_right = false;
    
    // Stick new node into ARG and maps.
    arg_node* newnode = new arg_node;
    newnode->is_leaving_grp = node->is_leaving_grp;
    newnode->skip = node->skip;
    newnode->site_support = node->site_support;
    newnode->clade = node->clade;
    newnode->sites = sites_right;
    newnode->mutations = muts_right;
    newnode->start = *(sites_right.begin());
    newnode->end = oldnodeend;
    newnode->recomb_left = split_pos;
    newnode->fixed_left = false;
    newnode->fixed_right = false;
    
    set<arg_node*> dummy;
    
    arg_node* stored = create_node(newnode, sites_pos, sites_clade, lv_grps, 
        prop_dist, root, newnodes, num_haplotypes, dummy);
    
    
    if (stored == NULL){
        fprintf(stderr, "split_node() could not create??\n");
        print_node_lite(node, num_haplotypes);
        exit(1);
    }
    else if (stored != newnode){
        delete newnode;
    }
    
    set<rcsolved_info> edges_add;
    for (vector<rc_edge>::iterator from = node->edges_from_solved.begin(); from != node->edges_from_solved.end(); ++from){
        if (from->final_dest != NULL){
            edges_add.insert(rcsolved_info(node, from->final_dest, from->node, from->type));
        }
        else{
            for (vector<rc_edge>::iterator to = from->node->edges_to_solved.begin();
                to != from->node->edges_to_solved.end(); ++to){
                if (to->final_dest != NULL && to->node == node){
                    edges_add.insert(rcsolved_info(to->final_dest, from->node, node, to->type));
                }
            }
        }
    }
    
    vector<arg_node*> dummy2;
    for (set<rcsolved_info>::iterator ea = edges_add.begin(); ea != edges_add.end(); ++ea){
        del_solved_connections_between(ea->left, ea->right);
        if (ea->left == node){

            add_solved_edges(newnode, ea->right, ea->lv, ea->type, false, dummy2, (long int)-1, num_haplotypes);
        }
        else if (ea->lv == node){

            add_solved_edges(ea->left, ea->right, newnode, ea->type, false, dummy2, (long int)-1, num_haplotypes);
        }
    }
    
    node->fixed_right = false;
    for (vector<rc_edge>::iterator from = node->edges_from.begin(); from != node->edges_from.end(); ++from){
        if (from->final_dest != NULL){
            from->final_dest->fixed_left = false;
        }
    }
    
    for (vector<rc_edge>::iterator from = node->edges_from_unsolvable.begin();
        from != node->edges_from_unsolvable.end();){
        from->final_dest->fixed_left = false;
        for (vector<rc_edge>::iterator to = from->final_dest->edges_to_unsolvable.begin();
            to != from->final_dest->edges_to_unsolvable.end();){
            if (to->final_dest == node){
                from->final_dest->edges_to_unsolvable.erase(to);
                //break;
            }
            else{
                ++to;
            }
        }
        node->edges_from_unsolvable.erase(from);
    }
    
    if (DEBUG_MODE){
        fprintf(stderr, "SPLIT NODE DONE\n");
    }
    
    return stored;
}

/**
 * Returns NULL (normal/success) OR a pointer to a new node created by splitting node2.
 *
 * parameter finalize = when true, sets 4hap test failures (which can trigger building
 * new recombination events). when false, does not do this - just sets ranges
 * appropriately.
 */
arg_node* handle_failure(arg_node* node1, arg_node* node2, 
    arg_sitemap& sites_pos,
    arg_clademap& sites_clade,
    arg_clademap& lv_grps,
    long int prop_dist,
    bool finalize,
    arg_node* root,
    vector<arg_node*>& newnodes,
    int num_haplotypes){
    
    if (DEBUG_MODE){
        fprintf(stderr, "handle failure; finalize = %d\n", finalize);
        print_node_lite(node1, num_haplotypes);
        print_node_lite(node2, num_haplotypes);

        for (set<long int>::iterator site1 = node1->sites.begin(); site1 != node1->sites.end(); ++site1){
            if (node2->sites.find(*site1) != node2->sites.end()){
                fprintf(stderr, "two clades that fail 4hap test exist at the same site\n");
                print_node_lite(node1, num_haplotypes);
                print_node_lite(node2, num_haplotypes);
                exit(1);
            }
        }
        
        if (node1->clade == node2->clade){
            fprintf(stderr, "failure between identical nodes?\n");
            exit(1);
        }
        
        // If there's any overlap in sites, we need to merge the two together.
        for (set<long int>::iterator site1 = node1->sites.begin(); site1 != node1->sites.end(); ++site1){
            if (node2->sites.find(*site1) != node2->sites.end()){
                //return true;
                fprintf(stderr, "failure between nodes with overlapping sites?!\n");
                exit(1);
            }
        }
    }
    
    long int start1 = *(node1->sites.begin());
    long int end1 = *(node1->sites.rbegin());
    long int start2 = *(node2->sites.begin());
    long int end2 = *(node2->sites.rbegin());

    if (start1 > start2 && end1 < end2){
        // Node1 within Node2 range
        
        if (DEBUG_MODE){
            // There should not be any sites belonging to node2 between node1 start & end
            for (set<long int>::iterator site = node2->sites.begin(); site != node2->sites.end(); ++site){
                if (*site > start1 && *site < end1){
                    fprintf(stderr, "node2 has site between node1 start & end\n");
                    print_node_lite(node1, num_haplotypes);
                    print_node_lite(node2, num_haplotypes);
                    exit(1);
                }
            }
        }
        
        // Shorten range on both ends
        long int prevstart = node1->start;
        if (prevstart != start1){
            node1->start = start1;
            shorten_range(node1, node1->start - start1, sites_pos, sites_clade, lv_grps, 
                prop_dist, num_haplotypes, root, newnodes, -1);
            node1->fixed_left = false;
            unfix_edges_left(node1);
        }
        long int prevend = node1->end;
        if (prevend != end1){
            node1->end = end1;
            shorten_range(node1, node1->end - prevend, sites_pos, sites_clade, lv_grps, 
                prop_dist, num_haplotypes, root, newnodes, -1);
            node1->fixed_right = false;
            unfix_edges_right(node1);
        }
        
        
        arg_node* newnode = split_node(node2, start1, sites_pos, sites_clade, 
            lv_grps, prop_dist, root, newnodes, num_haplotypes);
        
        node2->fixed_right = false;
        unfix_edges_right(node2);
        
        if (finalize){
            set_4haptest_failure(node2, node1, node2->end, node1->start, sites_pos, sites_clade,
                lv_grps, prop_dist, root, newnodes, num_haplotypes);
            set_4haptest_failure(node1, newnode, end1, newnode->start, sites_pos, sites_clade,
                lv_grps, prop_dist, root, newnodes, num_haplotypes);
        }
        return newnode;
    }
    else if (start2 > start1 && end2 < end1){
        // Node2 within Node1 range
        
        if (DEBUG_MODE){
            // There should not be any sites belonging to node1 between node2 start & end
            for (set<long int>::iterator site = node1->sites.begin(); site != node1->sites.end(); ++site){
                if (*site > start2 && *site < end2){
                    fprintf(stderr, "node1 has site between node2 start & end\n");
                    print_node_lite(node1, num_haplotypes);
                    print_node_lite(node2, num_haplotypes);
                    exit(1);
                }
            }
        }
        
        // Shorten range on both ends
        long int prevstart = node2->start;
        if (prevstart != start2){
            node2->start = start2;
            shorten_range(node2, node2->start - start2, sites_pos, sites_clade, lv_grps, 
                prop_dist, num_haplotypes, root, newnodes, -1);
            node2->fixed_left = false;
            unfix_edges_left(node2);
        }
        long int prevend = node2->end;
        if (prevend != end2){
            node2->end = end2;
            shorten_range(node2, node2->end - prevend, sites_pos, sites_clade, lv_grps, 
                prop_dist, num_haplotypes, root, newnodes, -1);
            node2->fixed_right = false;
            unfix_edges_right(node2);
        }
        
        // NOTE: node1 (the new node) should NEVER be split. This is because new nodes either
        // come from single sites (can't be split), or they are leaving groups (2 sites).
        // But for a leaving group to be valid/usable, it must not fail the 4hap test with
        // any other clade within its range.
        fprintf(stderr, "about to split node1! this should never happen\n");
        print_node_lite(node1, num_haplotypes);
        print_node_lite(node2, num_haplotypes);
        exit(1);
        /*
        arg_node* newnode = split_node(node1, start2, sites_pos, sites_clade, 
            lv_grps, prop_dist, root, newnodes, num_haplotypes);
        node1->fixed_right = false;
        unfix_edges_right(node1);
        
        if (finalize){
            set_4haptest_failure(node1, node2, node1->end, node2->start, sites_pos, sites_clade,
               lv_grps, prop_dist, root, newnodes, num_haplotypes);
            set_4haptest_failure(node2, newnode, end2, newnode->start, sites_pos, sites_clade,
                lv_grps, prop_dist, root, newnodes, num_haplotypes);
        }
        return newnode;
        */
        return NULL;
    }
    else if ((start1 < start2 && end1 > start2 && end1 < end2) || end1 < start2){
        
        node1->fixed_right = false;
        unfix_edges_right(node1);
        node2->fixed_left = false;
        unfix_edges_left(node2);
            
        // Node1 hangs off to left of Node2
        if (node1->recomb_right != -1 && node1->recomb_right <= start2 && node1->end < node2->start){
            // Don't alter node1 range.
        }
        else{
            long int prevend = node1->end;
            node1->end = end1;
            if (end1 != prevend){
                shorten_range(node1, node1->end - prevend, sites_pos, sites_clade, lv_grps,
                     prop_dist, num_haplotypes, root, newnodes, -1);
            }
            
        }
        if (node2->recomb_left != -1 && node2->recomb_left >= end1 && node1->end < node2->start){
            // Don't alter node2 range.
        }
        else{
            long int prevstart = node2->start;
            node2->start = start2;
            
            if (start2 != prevstart){
                shorten_range(node2, node2->start - prevstart, sites_pos, sites_clade, 
                    lv_grps, prop_dist, num_haplotypes, root, newnodes, -1);
            }
        }
        
        if (DEBUG_MODE && (node1->start > *node1->sites.begin() || node1->end < *node1->sites.rbegin() ||
            node2->start > *node2->sites.begin() || node2->end < *node2->sites.rbegin())){
            fprintf(stderr, "impossible range(1)\n");
            print_node_lite(node1, num_haplotypes);
            print_node_lite(node2, num_haplotypes);
            exit(1);
        }

        if (finalize){
            set_4haptest_failure(node1, node2, end1, start2, sites_pos, sites_clade,
                lv_grps, prop_dist, root, newnodes, num_haplotypes);
        }
       
        if (DEBUG_MODE && (node1->recomb_left > node1->start || (node1->recomb_right != -1 && node1->recomb_right < node1->end) ||
            node2->recomb_left > node2->start || (node2->recomb_right != -1 && node2->recomb_right < node2-> end))){
            fprintf(stderr, "invalid recomb ind(s)1\n");
            print_node_lite(node1, num_haplotypes);
            print_node_lite(node2, num_haplotypes);
            exit(1);
        }
        
        
        return NULL;
        
    }
    else if ((start2 < start1 && end2 > start1 && end2 < end1) || end2 < start1){
        
        node1->fixed_left = false;
        unfix_edges_left(node1);
        node2->fixed_right = false;
        unfix_edges_right(node2);
            
        // Node 1 hangs off to right of Node2
        if (node1->recomb_left != -1 && node1->recomb_left >= end2 && node2->end < node1->start){
            // Don't alter node1 range.
        }
        else{
            long int prevstart = node1->start;
            node1->start = start1;
            if (start1 != prevstart){
                shorten_range(node1, node1->start - prevstart, sites_pos, sites_clade, 
                    lv_grps, prop_dist, num_haplotypes, root, newnodes, -1);
            }
        }
        
        if (node2->recomb_right != -1 && node2->recomb_right <= start1 && node2->end < node1->start){
            // Don't alter node2 range.
        }
        else{
            long int prevend = node2->end;
            node2->end = end2;
            if (end2 != prevend){
                shorten_range(node2, node2->end - prevend, sites_pos, sites_clade, 
                    lv_grps, prop_dist, num_haplotypes, root, newnodes, -1);
            }
            
        }
        
        
        if (DEBUG_MODE && (node1->start > *node1->sites.begin() || node1->end < *node1->sites.rbegin() ||
            node2->start > *node2->sites.begin() || node2->end < *node2->sites.rbegin())){
            fprintf(stderr, "impossible range(2)\n");
            print_node_lite(node1, num_haplotypes);
            print_node_lite(node2, num_haplotypes);
            exit(1);
        }
        if (finalize){
            set_4haptest_failure(node2, node1, node2->end, node1->start, sites_pos, sites_clade,
                lv_grps, prop_dist, root, newnodes, num_haplotypes);
        }
        
        
        if (DEBUG_MODE && (node1->recomb_left > node1->start || (node1->recomb_right != -1 && node1->recomb_right < node1->end) ||
            node2->recomb_left > node2->start || (node2->recomb_right != -1 && node2->recomb_right < node2-> end))){
            fprintf(stderr, "invalid recomb ind(s)2\n");
            print_node_lite(node1, num_haplotypes);
            print_node_lite(node2, num_haplotypes);
            exit(1);
        }
        
        return NULL;
    }
    else{
        fprintf(stderr, "didn't set 4haptest failure %ld %ld | %ld %ld\n", start1, end1, start2, end2);
        print_node_lite(node1, num_haplotypes);
        print_node_lite(node2, num_haplotypes);
        
        for (set<long int>::iterator site = node1->sites.begin(); site != node1->sites.end(); ++site){
            pair<arg_sitemap::iterator, arg_sitemap::iterator> er = sites_pos.equal_range(*site);
            bool n1found = false;
            for (arg_sitemap::iterator n = er.first; n != er.second; ++n){
                if (n->second == node1){
                    n1found = true;
                    break;
                }
            }
            if (!n1found){
                fprintf(stderr, "!n1found; site %ld\n", *site);
            }
        }
        for (set<long int>::iterator site = node2->sites.begin(); site != node2->sites.end(); ++site){
            pair<arg_sitemap::iterator, arg_sitemap::iterator> er = sites_pos.equal_range(*site);
            bool n2found = false;
            for (arg_sitemap::iterator n = er.first; n != er.second; ++n){
                if (n->second == node2){
                    n2found = true;
                    break;
                }
            }
            if (!n2found){
                fprintf(stderr, "!n2found; site %ld\n", *site);
            }
        }
        exit(1);
    }
    
    if (DEBUG_MODE){
        fprintf(stderr, "hf DONE\n");
    }
    
    if (node1->recomb_left > node1->start || (node1->recomb_right != -1 && node1->recomb_right < node1->end) ||
        node2->recomb_left > node2->start || (node2->recomb_right != -1 && node2->recomb_right < node2-> end)){
        fprintf(stderr, "invalid recomb ind(s)3\n");
        print_node_lite(node1, num_haplotypes);
        print_node_lite(node2, num_haplotypes);
        exit(1);
    }
        
    return NULL;
}

/**
 * amount = newpos - oldpos
 * amount > 0 == to the right
 * amount < 0 == to the left
 */
void expand_range_newsites(arg_node* node, 
    long int amount, 
    arg_sitemap& sites_pos,
    arg_clademap& sites_clade,
    arg_clademap& lv_grps,
    long int prop_dist,
    bool newnode,
    arg_node* root,
    vector<arg_node*>& newnodes,
    int num_haplotypes){

    if (amount == 0){
        return;
    }
    else if (amount > 0){
        node->fixed_right = false;
    }
    else if (amount < 0){
        node->fixed_left = false;
    }

    if (amount > 0 && node->recomb_right != -1 && node->recomb_right < node->end){
        node->end = *(node->sites.rbegin());
        return;
    }
    else if (amount < 0 && node->recomb_left != -1 && node->recomb_left > node->start){
        node->start = *(node->sites.begin());
        return;
    }
    
    if (DEBUG_MODE){
        fprintf(stderr, "EXPRNS\n");
        print_node_lite(node, num_haplotypes);
        
        fprintf(stderr, "prevs:\n");
        pair<arg_clademap::iterator, arg_clademap::iterator> er = sites_clade.equal_range(node->clade);
        for (arg_clademap::iterator prev = er.first; prev != er.second; ++prev){
            print_node_lite(prev->second, num_haplotypes);
        }
    }
    
    // See if another node already exists that can be joined with this one.
    pair<arg_clademap::iterator, arg_clademap::iterator> er = sites_clade.equal_range(node->clade);
    for (arg_clademap::iterator prev = er.first; prev != er.second; ){
    
        if (DEBUG_MODE && prev->second->clade != node->clade){
            fprintf(stderr, "prevnode clade does not match node clade?\n");
            print_node_lite(prev->second, num_haplotypes);
            print_node_lite(node, num_haplotypes);
            exit(1);
        }
        bool removed = false;
        if (prev->second != node && node_ranges_overlap(prev->second, node)){
            
            if (DEBUG_MODE){
                fprintf(stderr, "checking prev:\n");
                print_node_lite(prev->second, num_haplotypes);
            }
            
            bool can_join = true;
            // See if there's a failure between them.
            if (amount < 0){

                if (prev->second->recomb_right != -1 && prev->second->recomb_right <= *node->sites.begin()){
                    if (node->sites.find(prev->second->recomb_right) != node->sites.end()){
                        fprintf(stderr, "Already-known recombination between otherwise-merge-able nodes\n");
                        print_node_lite(prev->second, num_haplotypes);
                        print_node_lite(node, num_haplotypes);
                        pair<arg_sitemap::iterator, arg_sitemap::iterator> er = sites_pos.equal_range(prev->second->recomb_right);
                        for (arg_sitemap::iterator n = er.first; n != er.second; ++n){
                            print_node_lite(n->second, num_haplotypes);
                        }
                        exit(1);
                    }
                    can_join = false;
                }
                else{

                    // Test sites in between.
                    arg_sitemap::iterator between = sites_pos.equal_range(*prev->second->sites.rbegin()).second;
                    while (between->first < *node->sites.begin() && between != sites_pos.end()){
                        
                        short comp = compare_bitsets(between->second->clade, node->clade);
                        if (comp == -1){
                            can_join = false;
                            break;
                        }
                        ++between;
                        if (between == sites_pos.end()){
                            fprintf(stderr, "went off the end of sites_pos(1)\n");
                            exit(1);
                        }
                    } 
                }
                
            }
            else{
                if (prev->second->recomb_left != -1 && prev->second->recomb_left > *node->sites.rbegin()){
                    if (node->sites.find(prev->second->recomb_left) != node->sites.end()){
                        fprintf(stderr, "Already-known recombination between otherwise merge-able nodes (2)\n");
                        print_node_lite(prev->second, num_haplotypes);
                        print_node_lite(node, num_haplotypes);
                        pair<arg_sitemap::iterator, arg_sitemap::iterator> er2 = sites_pos.equal_range(prev->second->recomb_left);
                        for (arg_sitemap::iterator n = er2.first; n != er2.second; ++n){
                            print_node_lite(n->second, num_haplotypes);
                        }
                        exit(1);
                    }
                    can_join = false;
                }
                else{

                    // Test sites in between.
                    arg_sitemap::iterator between = sites_pos.equal_range(*node->sites.rbegin()).second;
                    between++;
                    //arg_sitemap::iterator between = sites_pos.find(*node->sites.rbegin());
                    while (between != sites_pos.end() && between->first < *prev->second->sites.begin()){
                        if (DEBUG_MODE){
                            fprintf(stderr, "checking site %ld\n", between->first);

                        }
                        
                        short comp = compare_bitsets(between->second->clade, node->clade);
                        if (comp == -1){
                            can_join = false;
                            break;
                        }
          
                        ++between;
                        if (between == sites_pos.end()){
                            fprintf(stderr, "went off the end of sites_pos (2)\n");
                            exit(1);
                        }
                    }
                    if (DEBUG_MODE){
                        fprintf(stderr, "done checking\n");
                    }
                }
            }
            
            if (can_join){
                
                arg_node* othernode = prev->second;
                
                
                if (DEBUG_MODE){
                    fprintf(stderr, "CAN JOIN:\n");
                    print_node_lite(node, num_haplotypes);
                    print_node_lite(othernode, num_haplotypes);
                    fprintf(stderr, "%ld\n", amount);
                }
                
                
                // First let the expanding node pick up sites leading up to the other node
                // that will be joined.
                if (amount > 0){
                    long int before_func_end = node->end - amount;
                    
                
                    if (DEBUG_MODE){
                        fprintf(stderr, "before_func_end %ld\n", before_func_end);
                    }
                    if (othernode->start > before_func_end && othernode->start-1 > before_func_end){
     
                        long int prevend = node->end;
                        node->end = othernode->start-1;
                        // Take off sites too far for this range; add back afterward.
                        
                        set<long int> sites_stash;
                        for (set<long int>::iterator site = node->sites.begin();
                            site != node->sites.end();){
                            if (*site > node->end){
                                sites_stash.insert(*site);
                                node->sites.erase(site++);
                            }
                            else{
                                ++site;
                            }
                        }
                        
                        if (DEBUG_MODE){
                            fprintf(stderr, "expR1before\n");
                        }

                        expand_range(node, node->end-before_func_end, sites_pos, 
                            sites_clade, lv_grps, prop_dist, root, newnodes, num_haplotypes, false);
                        if (DEBUG_MODE){
                            fprintf(stderr, "expR1here\n");
                        }
                        
                        // Add back stashed sites.
                        for (set<long int>::iterator site = sites_stash.begin();
                            site != sites_stash.end();){
                            node->sites.insert(*site);
                            sites_stash.erase(site++);
                        }
                        
                        for (vector<rc_edge>::iterator from = node->edges_from.begin();
                            from != node->edges_from.end(); ++from){
                            if (from->final_dest != NULL && from->final_dest->recomb_left != -1 &&
                                from->final_dest->recomb_left < *node->sites.rbegin()){
                                from->final_dest->recomb_left = *node->sites.rbegin();   
                            }
                        }
                        for (vector<rc_edge>::iterator from = node->edges_from_solved.begin();
                            from != node->edges_from_solved.end(); ++from){
                            if (from->final_dest != NULL && from->final_dest->recomb_left != -1 &&
                                from->final_dest->recomb_left < *node->sites.rbegin()){
                                from->final_dest->recomb_left = *node->sites.rbegin();   
                            }
                        }
                        for (vector<rc_edge>::iterator from = node->edges_from_unsolvable.begin();
                            from != node->edges_from_unsolvable.end(); ++from){
                            if (from->final_dest != NULL && from->final_dest->recomb_left != -1 &&
                                from->final_dest->recomb_left < *node->sites.rbegin()){
                                from->final_dest->recomb_left = *node->sites.rbegin();   
                            }
                        }
                    
                        if (prevend > othernode->end){
                            //amount = prevend - othernode->end;
                            long int rightside_amount = prevend - othernode->end;
                            othernode->end = prevend;
                            expand_range(othernode, rightside_amount, sites_pos, 
                                sites_clade, lv_grps, prop_dist, root, newnodes, num_haplotypes, false);
                            
                        }
                        amount = 0;
                    }
                    else if (DEBUG_MODE){
                        fprintf(stderr, "AMOUNT %ld node %ld %ld  before_end %ld othernode %ld %ld\n", 
                            amount, node->start, node->end, before_func_end, othernode->start, othernode->end);
                    }
                }
                else if (amount < 0){
                    long int before_func_start = node->start - amount;
                    if (DEBUG_MODE){
                        fprintf(stderr, "before_func_start %ld\n", before_func_start);
                    }
                    if (othernode->end < before_func_start && othernode->end+1 < before_func_start){
                        long int prevstart = node->start;
                        node->start = othernode->end+1;
                        if (DEBUG_MODE){
                            fprintf(stderr, "expR2before\n");
                        }
                        
                        set<long int> sites_stash;
                        for (set<long int>::iterator site = node->sites.begin();
                            site != node->sites.end();){
                            if (*site < node->start){
                                sites_stash.insert(*site);
                                node->sites.erase(site++);
                            }
                            else{
                                ++site;
                            }
                        }

                        expand_range(node, node->start-before_func_start, 
                            sites_pos, sites_clade, lv_grps, prop_dist, 
                            root, newnodes, num_haplotypes, false);
                        if (DEBUG_MODE){
                            fprintf(stderr, "expR2after\n");
                        }
                        
                            
                        // Add back stashed sites.
                        for (set<long int>::iterator site = sites_stash.begin();
                            site != sites_stash.end();){
                            node->sites.insert(*site);
                            sites_stash.erase(site++);
                        }
                        
                        for (vector<rc_edge>::iterator to = node->edges_to.begin();
                            to != node->edges_to.end(); ++to){
                            if (to->final_dest != NULL && to->final_dest->recomb_right != -1 &&
                                to->final_dest->recomb_right > *node->sites.begin()){
                                to->final_dest->recomb_right = *node->sites.begin();   
                            }
                        }
                        for (vector<rc_edge>::iterator to = node->edges_to_solved.begin();
                            to != node->edges_to_solved.end(); ++to){
                            if (to->final_dest != NULL && to->final_dest->recomb_right != -1 &&
                                to->final_dest->recomb_right > *node->sites.begin()){
                                to->final_dest->recomb_right = *node->sites.begin();   
                            }
                        }
                        for (vector<rc_edge>::iterator to = node->edges_to_unsolvable.begin();
                            to != node->edges_to_unsolvable.end(); ++to){
                            if (to->final_dest != NULL && to->final_dest->recomb_right != -1 &&
                                to->final_dest->recomb_right > *node->sites.begin()){
                                to->final_dest->recomb_right = *node->sites.begin();   
                            }
                        }
                            
                        if (prevstart < othernode->start){
                            
                            long int leftside_amount = prevstart - othernode->start;
                            othernode->start = prevstart;

                            expand_range(othernode, leftside_amount, sites_pos, 
                                sites_clade, lv_grps, prop_dist, root, newnodes,
                                num_haplotypes, false);
                            
                            
                        }
                        
                        amount = 0;
                    }
                    else if (DEBUG_MODE){
                        fprintf(stderr, "AMOUNT %ld node %ld %ld before_start %ld othernode %ld %ld\n", 
                            amount, node->start, node->end, before_func_start, othernode->start, othernode->end);
                    }
                }
                else if (DEBUG_MODE){
                    fprintf(stderr, "AMOUNT %ld node %ld %ld othernode %ld %ld\n", 
                        amount, node->start, node->end, othernode->start, othernode->end);
                }
                
                join_nodes(node, othernode, sites_pos, sites_clade, lv_grps, 
                    prop_dist, false, root, newnodes, num_haplotypes);
                
                update_recinds_left(node);
                update_recinds_right(node);
                
                // Need to get back the iterator (we've just messed with sites_pos and sites_clade).
                removed = true;
            }
            if (DEBUG_MODE){
                fprintf(stderr, "checking done\n");
            }
        }
        if (removed){
            arg_node* nodeptr = prev->second;
            sites_clade.erase(prev++);
            
            delete nodeptr;
        }
        else{
            ++prev;
        }
    }
    
    if (DEBUG_MODE){
        fprintf(stderr, "did not merge\n");
    }
    
    if (amount != 0){
        // Adjust current children/parents to inform them of the change.
        expand_range(node, amount, sites_pos, sites_clade, lv_grps, prop_dist, 
            root, newnodes, num_haplotypes, false);
    }
    
    return;
    
}


void print_node_recursive(arg_node* n, int indent_level, int num_haplotypes){
    string indent_str = "";
    for (int i = 0; i < indent_level; ++i){
        indent_str += "\t";
    }
    fprintf(stderr, "%s(%ld, %ld)", indent_str.c_str(), n->start, n->end);
    if (n->is_new_left){
        fprintf(stderr, " L ");
    }
    if (n->is_leaving_grp){
        fprintf(stderr, " * ");
    }
    if (n->is_new_right){
        fprintf(stderr, " R ");
    }
    fprintf(stderr, " {%ld} {%ld}", n->recomb_left, n->recomb_right);
    fprintf(stderr, "\n");
    fprintf(stderr, "%s[ ", indent_str.c_str());
    for (set<long int>::iterator site = n->sites.begin(); site != n->sites.end(); ++site){
        fprintf(stderr, "%ld ", *site);
    }
    fprintf(stderr, "]\n");
    fprintf(stderr, "%s{ ", indent_str.c_str());
    for (set<long int>::iterator site = n->mutations.begin(); site != n->mutations.end(); ++site){
        fprintf(stderr, "%ld ", *site);
    }
    fprintf(stderr, "}\n");
    string members = indent_str + "(";
    for (unsigned int i = 0; i < num_haplotypes; ++i){
        if (n->clade.test(num_haplotypes-i-1)){
            char buf[4];
            sprintf(buf, "%d", i);
            string bufstr(buf);
            members += bufstr + ",";
        }
    }
    members = members.substr(0, members.length()-1);
    members += ")";
    fprintf(stderr, "%s\n", members.c_str());
    for (set<vert_edge>::iterator child = n->children.begin(); child != n->children.end();
        ++child){
        fprintf(stderr, "%s\t<%ld %ld>\n", indent_str.c_str(), child->start, child->end);
        print_node_recursive(child->to, indent_level+1, num_haplotypes);
    }
}

void print_node_recursive_upward(arg_node* n, int indent_level, int num_haplotypes){
    string indent_str = "";
    for (int i = 0; i < indent_level; ++i){
        indent_str += "\t";
    }
    fprintf(stderr, "%s(%ld, %ld)", indent_str.c_str(), n->start, n->end);
    if (n->is_new_left){
        fprintf(stderr, " L ");
    }
    if (n->is_leaving_grp){
        fprintf(stderr, " * ");
    }
    if (n->is_new_right){
        fprintf(stderr, " R ");
    }
    fprintf(stderr, " {%ld} {%ld}", n->recomb_left, n->recomb_right);
    fprintf(stderr, "\n");
    fprintf(stderr, "%s[ ", indent_str.c_str());
    for (set<long int>::iterator site = n->sites.begin(); site != n->sites.end(); ++site){
        fprintf(stderr, "%ld ", *site);
    }
    fprintf(stderr, "]\n");
    fprintf(stderr, "%s{ ", indent_str.c_str());
    for (set<long int>::iterator site = n->mutations.begin(); site != n->mutations.end(); ++site){
        fprintf(stderr, "%ld ", *site);
    }
    fprintf(stderr, "}\n");
    string members = indent_str + "(";
    for (unsigned int i = 0; i < num_haplotypes; ++i){
        if (n->clade.test(num_haplotypes-i-1)){
            char buf[4];
            sprintf(buf, "%d", i);
            string bufstr(buf);
            members += bufstr + ",";
        }
    }
    members = members.substr(0, members.length()-1);
    members += ")";
    fprintf(stderr, "%s\n", members.c_str());
    for (set<vert_edge>::iterator parent = n->parents.begin(); parent != n->parents.end();
        ++parent){
        fprintf(stderr, "%s\t<%ld %ld>\n", indent_str.c_str(), parent->start, parent->end);
        print_node_recursive_upward(parent->to, indent_level+1, num_haplotypes);
    }
}

void print_node_recursive_pos(arg_node* n, long int pos, int indent_level, int num_haplotypes){
    string indent_str = "";
    for (int i = 0; i < indent_level; ++i){
        indent_str += "\t";
    }
    fprintf(stderr, "%s(%ld, %ld)", indent_str.c_str(), n->start, n->end);
    if (!n->site_support){
        fprintf(stderr, " N ");
    }
    if (n->is_new_left){
        fprintf(stderr, " L ");
    }
    if (n->is_leaving_grp){
        fprintf(stderr, " * ");
    }
    if (n->is_new_right){
        fprintf(stderr, " R ");
    }
    fprintf(stderr, " {%ld} {%ld}", n->recomb_left, n->recomb_right);
    fprintf(stderr, "\n");
    fprintf(stderr, "%s[ ", indent_str.c_str());
    for (set<long int>::iterator site = n->sites.begin(); site != n->sites.end(); ++site){
        fprintf(stderr, "%ld ", *site);
    }
    fprintf(stderr, "]\n");
    fprintf(stderr, "%s{ ", indent_str.c_str());
    for (set<long int>::iterator site = n->mutations.begin(); site != n->mutations.end(); ++site){
        fprintf(stderr, "%ld ", *site);
    }
    fprintf(stderr, "}\n");
    string members = indent_str + "(";
    for (unsigned int i = 0; i < num_haplotypes; ++i){
        if (n->clade.test(num_haplotypes-i-1)){
            char buf[4];
            sprintf(buf, "%d", i);
            string bufstr(buf);
            members += bufstr + ",";
        }
    }
    members = members.substr(0, members.length()-1);
    members += ")";
    fprintf(stderr, "%s\n", members.c_str());
    for (set<vert_edge>::iterator child = n->children.begin(); child != n->children.end();
        ++child){
        if (child->start <= pos && child->end >= pos){
            fprintf(stderr, "%s\t<%ld %ld>\n", indent_str.c_str(), child->start, child->end);
            print_node_recursive_pos(child->to, pos, indent_level+1, num_haplotypes);
        }
    }
}

void print_node_recursive_upward_pos(arg_node* n, long int pos, int indent_level, int num_haplotypes){
    string indent_str = "";
    for (int i = 0; i < indent_level; ++i){
        indent_str += "\t";
    }
    fprintf(stderr, "%s(%ld, %ld)", indent_str.c_str(), n->start, n->end);
    if (n->is_new_left){
        fprintf(stderr, " L ");
    }
    if (n->is_leaving_grp){
        fprintf(stderr, " * ");
    }
    if (n->is_new_right){
        fprintf(stderr, " R ");
    }
    fprintf(stderr, " {%ld} {%ld}", n->recomb_left, n->recomb_right);
    fprintf(stderr, "\n");
    fprintf(stderr, "%s[ ", indent_str.c_str());
    for (set<long int>::iterator site = n->sites.begin(); site != n->sites.end(); ++site){
        fprintf(stderr, "%ld ", *site);
    }
    fprintf(stderr, "]\n");
    fprintf(stderr, "%s{ ", indent_str.c_str());
    for (set<long int>::iterator site = n->mutations.begin(); site != n->mutations.end(); ++site){
        fprintf(stderr, "%ld ", *site);
    }
    fprintf(stderr, "}\n");
    string members = indent_str + "(";
    for (unsigned int i = 0; i < num_haplotypes; ++i){
        if (n->clade.test(num_haplotypes-i-1)){
            char buf[4];
            sprintf(buf, "%d", i);
            string bufstr(buf);
            members += bufstr + ",";
        }
    }
    members = members.substr(0, members.length()-1);
    members += ")";
    fprintf(stderr, "%s\n", members.c_str());
    for (set<vert_edge>::iterator parent = n->parents.begin(); parent != n->parents.end();
        ++parent){
        if (parent->start <= pos && parent->end >= pos){
            fprintf(stderr, "%s\t<%ld %ld>\n", indent_str.c_str(), parent->start, parent->end);
            print_node_recursive_upward_pos(parent->to, pos, indent_level+1, num_haplotypes);
        }
    }
}

void print_node_upward_onelevel(arg_node* n, int indent_level, int num_haplotypes){
    string indent_str = "";
    for (int i = 0; i < indent_level; ++i){
        indent_str += "\t";
    }
    fprintf(stderr, "%s(%ld, %ld)", indent_str.c_str(), n->start, n->end);
    if (n->is_new_left){
        fprintf(stderr, " L ");
    }
    if (n->is_leaving_grp){
        fprintf(stderr, " * ");
    }
    if (n->is_new_right){
        fprintf(stderr, " R ");
    }
    fprintf(stderr, " {%ld} {%ld}", n->recomb_left, n->recomb_right);
    fprintf(stderr, "\n");
    fprintf(stderr, "%s[ ", indent_str.c_str());
    for (set<long int>::iterator site = n->sites.begin(); site != n->sites.end(); ++site){
        fprintf(stderr, "%ld ", *site);
    }
    fprintf(stderr, "]\n");
    fprintf(stderr, "%s{ ", indent_str.c_str());
    for (set<long int>::iterator site = n->mutations.begin(); site != n->mutations.end(); ++site){
        fprintf(stderr, "%ld ", *site);
    }
    fprintf(stderr, "}\n");
    string members = indent_str + "(";
    for (unsigned int i = 0; i < num_haplotypes; ++i){
        if (n->clade.test(num_haplotypes-i-1)){
            char buf[4];
            sprintf(buf, "%d", i);
            string bufstr(buf);
            members += bufstr + ",";
        }
    }
    members = members.substr(0, members.length()-1);
    members += ")";
    fprintf(stderr, "%s\n", members.c_str());
    
    if (indent_level == 0){
        for (set<vert_edge>::iterator parent = n->parents.begin(); parent != n->parents.end();
            ++parent){
            fprintf(stderr, "%s\t<%ld %ld>\n", indent_str.c_str(), parent->start, parent->end);
            print_node_upward_onelevel(parent->to, indent_level+1, num_haplotypes);
            fprintf(stderr, "\n");
        }
    }
}

void print_node_downward_onelevel(arg_node* n, int indent_level, int num_haplotypes){
    string indent_str = "";
    for (int i = 0; i < indent_level; ++i){
        indent_str += "\t";
    }
    fprintf(stderr, "%s(%ld, %ld)", indent_str.c_str(), n->start, n->end);
    if (n->is_new_left){
        fprintf(stderr, " L ");
    }
    if (n->is_leaving_grp){
        fprintf(stderr, " * ");
    }
    if (n->is_new_right){
        fprintf(stderr, " R ");
    }
    fprintf(stderr, " {%ld} {%ld}", n->recomb_left, n->recomb_right);
    fprintf(stderr, "\n");
    fprintf(stderr, "%s[ ", indent_str.c_str());
    for (set<long int>::iterator site = n->sites.begin(); site != n->sites.end(); ++site){
        fprintf(stderr, "%ld ", *site);
    }
    fprintf(stderr, "]\n");
    fprintf(stderr, "%s{ ", indent_str.c_str());
    for (set<long int>::iterator site = n->mutations.begin(); site != n->mutations.end(); ++site){
        fprintf(stderr, "%ld ", *site);
    }
    fprintf(stderr, "}\n");
    string members = indent_str + "(";
    for (unsigned int i = 0; i < num_haplotypes; ++i){
        if (n->clade.test(num_haplotypes-i-1)){
            char buf[4];
            sprintf(buf, "%d", i);
            string bufstr(buf);
            members += bufstr + ",";
        }
    }
    members = members.substr(0, members.length()-1);
    members += ")";
    fprintf(stderr, "%s\n", members.c_str());
    
    if (indent_level == 0){
        for (set<vert_edge>::iterator child = n->children.begin(); child != n->children.end();
            ++child){
            fprintf(stderr, "%s\t<%ld %ld>\n", indent_str.c_str(), child->start, child->end);
            print_node_downward_onelevel(child->to, indent_level+1, num_haplotypes);
            fprintf(stderr, "\n");
        }
    }
}


void print_node_lite(arg_node* n, int num_haplotypes){
    fprintf(stderr, "(%ld, %ld)", n->start, n->end);
    if (!n->site_support){
        fprintf(stderr, " N");
    }
    if (n->is_leaving_grp){
        fprintf(stderr, " * (%ld) (%ld)", n->edges_to.size(), n->edges_from.size());
    }
    fprintf(stderr, " {%ld} {%ld}", n->recomb_left, n->recomb_right);
    fprintf(stderr, "\n");
    fprintf(stderr, "[ ");
    for (set<long int>::iterator site = n->sites.begin(); site != n->sites.end(); ++site){
        fprintf(stderr, "%ld ", *site);
    }
    fprintf(stderr, "]\n");
    // Don't print all mutations for the root node
    if (n->clade.count() != num_haplotypes){
        fprintf(stderr, "{ ");
        for (set<long int>::iterator mut = n->mutations.begin(); mut != n->mutations.end(); ++mut){
            fprintf(stderr, "%ld ", *mut);
        }
        fprintf(stderr, "}\n");
    }
    string members = "(";
    for (unsigned int i = 0; i < num_haplotypes; ++i){
        if (n->clade.test(num_haplotypes-i-1)){
            char buf[4];
            sprintf(buf, "%d", i);
            string bufstr(buf);
            members += bufstr + ",";
        }
    }
    members = members.substr(0, members.length()-1);
    members += ")";
    fprintf(stderr, "%s\n", members.c_str());
    
}

void print_graph_sif(set<arg_node*>& lv_grps_all, int num_haplotypes){
    fprintf(stderr, "\n");
    set<arg_node*> lnodes;
    set<arg_node*> rnodes;
    
    for (set<arg_node*>::iterator lv = lv_grps_all.begin(); lv != lv_grps_all.end(); ++lv){
        string lvname = "* " + node2name(*lv, num_haplotypes);
        for (vector<rc_edge>::iterator to = (*lv)->edges_to.begin(); to != (*lv)->edges_to.end();
            ++to){
            if (true){
            //if (*to->node->sites.rbegin() <= *main_left->sites.rbegin()){
                for (vector<rc_edge>::iterator from = to->node->edges_from.begin();
                    from != to->node->edges_from.end(); ++from){
                    if (from->node == *lv){
                        if (from->final_dest != NULL){
                            string lname = "L " + node2name(to->node, num_haplotypes);
                            string fname = "R " + node2name(from->final_dest, num_haplotypes);
                            fprintf(stderr, "%s\trecomb_%d\t%s\n", lname.c_str(), to->type,
                                lvname.c_str());
                            fprintf(stderr, "%s\tfinal_%d\t%s\n", lname.c_str(), to->type,
                                fname.c_str());
                        }
                        else{
                            fprintf(stderr, "OOPS; final dest is null\n");
                            print_node_lite(to->node, num_haplotypes);
                            print_node_lite(*lv, num_haplotypes);
                            fprintf(stderr, "%d\n", to->type);
                            exit(1);
                        }
                    }
                }
            }
            
        }
        for (vector<rc_edge>::iterator from = (*lv)->edges_from.begin(); from != (*lv)->edges_from.end();
            ++from){
            if (true){
            //if (*from->node->sites.begin() >= *main_right->sites.begin()){
                for (vector<rc_edge>::iterator to = from->node->edges_to.begin();
                    to != from->node->edges_to.end(); ++to){
                    if (to->node == *lv){
                        if (to->final_dest != NULL){
                            string rname = "R " + node2name(from->node, num_haplotypes);
                            string fname = "L " + node2name(to->final_dest, num_haplotypes);
                            fprintf(stderr, "%s\trecomb_%d\t%s\n", lvname.c_str(), from->type,
                                rname.c_str());
                            fprintf(stderr, "%s\tfinal_%d\t%s\n", rname.c_str(), from->type,
                                fname.c_str());
                        }
                        else{
                            fprintf(stderr, "OOPS; final dest is null2\n");
                            print_node_lite(from->node, num_haplotypes);
                            print_node_lite(*lv, num_haplotypes);
                            fprintf(stderr, "%d\n", from->type);
                            exit(1);
                        }
                    }
                }
            }
        }
    }
    fprintf(stderr, "\n");
}

void print_graph_sif_lim(set<arg_node*>& lv_grps_all, set<arg_node*>& lnodes, 
    set<arg_node*>& rnodes, int num_haplotypes){
    fprintf(stderr, "\n");
    
    for (set<arg_node*>::iterator lv = lv_grps_all.begin(); lv != lv_grps_all.end(); ++lv){
        string lvname = "* " + node2name(*lv, num_haplotypes);
        for (vector<rc_edge>::iterator to = (*lv)->edges_to.begin(); to != (*lv)->edges_to.end();
            ++to){
            if (lnodes.find(to->node) != lnodes.end()){
                for (vector<rc_edge>::iterator from = to->node->edges_from.begin();
                    from != to->node->edges_from.end(); ++from){
                    if (from->node == *lv && from->final_dest != NULL && rnodes.find(from->final_dest) != rnodes.end()){
                        string lname = "L " + node2name(to->node, num_haplotypes);
                        string fname = "R " + node2name(from->final_dest, num_haplotypes);
                        fprintf(stderr, "%s\trecomb_%d\t%s\n", lname.c_str(), to->type,
                            lvname.c_str());
                        fprintf(stderr, "%s\tfinal_%d\t%s\n", lname.c_str(), to->type,
                            fname.c_str());
                    }
                }
            }
        }
        for (vector<rc_edge>::iterator from = (*lv)->edges_from.begin(); from != (*lv)->edges_from.end();
            ++from){
            if (rnodes.find(from->node) != rnodes.end()){
                for (vector<rc_edge>::iterator to = from->node->edges_to.begin();
                    to != from->node->edges_to.end(); ++to){
                    if (to->node == *lv && to->final_dest != NULL && lnodes.find(to->final_dest) != lnodes.end()){
                        string rname = "R " + node2name(from->node, num_haplotypes);
                        string fname = "L " + node2name(to->final_dest, num_haplotypes);
                        fprintf(stderr, "%s\trecomb_%d\t%s\n", lvname.c_str(), from->type,
                            rname.c_str());
                        fprintf(stderr, "%s\tfinal_%d\t%s\n", rname.c_str(), from->type,
                            fname.c_str());
                    }
                }
            }
        }
    }
    fprintf(stderr, "done\n");
    fprintf(stderr, "\n");
}

void print_recombs_right(arg_node* n, int num_haplotypes){
    for (vector<rc_edge>::iterator from = n->edges_from.begin(); from != n->edges_from.end();
        ++from){
        fprintf(stderr, "FROM (%d):\n", from->type);
        print_node_lite(from->node, num_haplotypes);
        if (from->final_dest != NULL){
            fprintf(stderr, "fd:\n");
            print_node_lite(from->final_dest, num_haplotypes);
        }
    }
    for (vector<rc_edge>::iterator from = n->edges_from_solved.begin(); from != n->edges_from_solved.end();
        ++from){
        fprintf(stderr, "FROM solved (%d):\n", from->type);
        print_node_lite(from->node, num_haplotypes);
        if (from->final_dest != NULL){
            fprintf(stderr, "fd:\n");
            print_node_lite(from->final_dest, num_haplotypes);
        }
    }
    for (vector<rc_edge>::iterator from = n->edges_from_unsolvable.begin(); from != n->edges_from_unsolvable.end();
        ++from){
        fprintf(stderr, "FROM unsolvable:\n");
        print_node_lite(from->final_dest, num_haplotypes);
    }
    fprintf(stderr, "\n");
}

void print_recombs_right_solved(arg_node* n, int num_haplotypes){
    for (vector<rc_edge>::iterator from = n->edges_from_solved.begin(); from != n->edges_from_solved.end();
        ++from){
        fprintf(stderr, "FROM solved (%d):\n", from->type);
        print_node_lite(from->node, num_haplotypes);
        if (from->final_dest != NULL){
            fprintf(stderr, "fd:\n");
            print_node_lite(from->final_dest, num_haplotypes);
        }
    }
    fprintf(stderr, "\n");
}

void print_recombs_left(arg_node* n, int num_haplotypes){
    for (vector<rc_edge>::iterator to = n->edges_to.begin(); to != n->edges_to.end(); ++to){
        fprintf(stderr, "TO (%d):\n", to->type);
        print_node_lite(to->node, num_haplotypes);
        if (to->final_dest != NULL){
            fprintf(stderr, "fd:\n");
            print_node_lite(to->final_dest, num_haplotypes);
        }
    }
    for (vector<rc_edge>::iterator to = n->edges_to_solved.begin(); to != n->edges_to_solved.end(); ++to){
        fprintf(stderr, "TO solved (%d):\n", to->type);
        print_node_lite(to->node, num_haplotypes);
        if (to->final_dest != NULL){
            fprintf(stderr, "fd:\n");
            print_node_lite(to->final_dest, num_haplotypes);
        }
    }
    for (vector<rc_edge>::iterator to = n->edges_to_unsolvable.begin(); to != n->edges_to_unsolvable.end(); ++to){
        fprintf(stderr, "TO unsolvable:\n");
        print_node_lite(to->final_dest, num_haplotypes);
    }
    fprintf(stderr, "\n");
}

void print_recombs_left_solved(arg_node* n, int num_haplotypes){
    for (vector<rc_edge>::iterator to = n->edges_to_solved.begin(); to != n->edges_to_solved.end(); ++to){
        fprintf(stderr, "TO solved (%d):\n", to->type);
        print_node_lite(to->node, num_haplotypes);
        if (to->final_dest != NULL){
            fprintf(stderr, "fd:\n");
            print_node_lite(to->final_dest, num_haplotypes);
        }
    }
    fprintf(stderr, "\n");
}



string node2name(arg_node* node, int num_haplotypes){
    char name[100];
    sprintf(name, "(%ld %ld) [", node->start, node->end);
    string namestr(name);
    for (int i = num_haplotypes-1; i >= 0; --i){
        if (node->clade.test(i)){
            char indv[5];
            sprintf(indv, "%d,", num_haplotypes-i-1);
            string indvstr(indv);
            namestr += indvstr;
        }
    }
    namestr += "]";
    return namestr;
}

