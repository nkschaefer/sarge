#include <getopt.h>
#include <argp.h>
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
#include "argnode.h"
#include "recomb.h"
#include "debug.h"
#include "rc_trees.h"

using std::cout;
using std::endl;
using namespace std;

//arg_node* root;
//long int last_deleted;

std::vector<std::string> indvs;

/**
 * Function to read in a chunk from the input file.
 * Returns the number of new sites added; modifies hapList, infile, and eof by reference.
 */
long int expand_matrix(deque<cladeset >& site_bitsets,
    long int bitsets_popped,
    gzFile &infile, 
    long int bufsize,
    const unsigned int num_haplotypes,
    bool &eof,
    long int prop_dist){
    
    static long int site_index = 0;
    
    long unsigned int new_sites = 0;
    char* chunk = (char*) malloc(bufsize+1);
    //char chunk[bufsize+1];
    
    int chunk_size;
    chunk_size = gzread(infile, chunk, bufsize);
    chunk[chunk_size] = '\0';
    
    eof = gzeof(infile);
    
    if (!eof){
        // Find the last newline in the chunk and rewind to just before that position.
        int steps_back = 0;
        // Should not have to do this if bufsize is properly set.
        if (chunk[chunk_size-1] != '\n'){
            for (int i=chunk_size-1; i >= 0; i--){
                if (chunk[i] == '\n'){
                    // Truncate the string here.
                    chunk[i] = '\0';
                    gzseek(infile, -steps_back, SEEK_CUR);
                    break;
                }
                steps_back++;
            }
        }
        else{
            chunk[chunk_size-1] = '\0';
        }
    }
    
    // Split lines   
    istringstream linesplitter(chunk);
    string line;
    int token_index = 0;

    while(std::getline(linesplitter, line, '\n')){
        if (line.size() == num_haplotypes){
            new_sites++;
                
            // Convert to bit string and store
            cladeset bitstr(line);
            
            site_bitsets.push_back(bitstr);
            
            // Check to see if this is marked "deleted"
            if (del_sites.size() > 0 && *del_sites.begin() == site_index){
                // Set current site to all ancestral so it will be ignored.
                site_bitsets[site_index-bitsets_popped].reset();
                del_sites.erase(del_sites.begin());
            }
        }
        ++site_index;
    }
    
    // Avoid letting garbage into this buffer next time
    memset(chunk, 0, bufsize);
    
    free(chunk);
    
    return new_sites;
}

void rcnodes2tree(arg_node* node, treeNode* tree, long int pos, long int prop_dist,
    int num_haplotypes){
    
    for (set<vert_edge>::iterator child = node->children.begin(); child != node->children.end();
        ++child){
        if (child->start <= pos && child->end >= pos){
            treeNode* childNode = new treeNode(num_haplotypes, child->to->clade);
            for (set<long int>::iterator site = child->to->sites.begin(); site != child->to->sites.end();
                ++site){
                if (*site >= pos - prop_dist && *site <= pos + prop_dist){
                    childNode->dist++;
                }
            }
            //childNode->dist += child->to->sites.size();
            tree->leaves = set_diff_bitset(tree->leaves, childNode->leaves);
            tree->children.push_back(childNode);
            childNode->parent = tree;
            tree->clear_cache();
            rcnodes2tree(child->to, childNode, pos, prop_dist, num_haplotypes);
        }
    }
}

void print_site_tree(arg_node* root, long int pos, long int prop_dist, int num_haplotypes){
    treeNode* t = new treeNode();
    t->set_haps(num_haplotypes);
    for (int i = 0; i < num_haplotypes; ++i){
        t->leaves.set(num_haplotypes-i-1);
    }
    rcnodes2tree(root, t, pos, prop_dist, num_haplotypes);
    vector<string> dummy;
    fprintf(stderr, "%s\n", t->newick(false, false, dummy).c_str());
    delete t;
}

bool clade_exists(arg_node* node, long int pos, const cladeset& clade){
    if (clade.count() <= 1){
        return true;
    }
    if (node->clade == clade){
        return true;
    }
    else{
        short comp = compare_bitsets(node->clade, clade);
        if (comp == -1){
            return true;
        }
        else if (comp == 2){
            return false;
        }
    }
    for (set<vert_edge>::iterator child = node->children.begin();
        child != node->children.end(); ++child){
        if (child->start <= pos && child->end >= pos){
            if (clade_exists(child->to, pos, clade)){
                return true;
            }
        }
    }
    return false;
}


void get_tree_flat(arg_node* node, 
    long int pos,
    unordered_set<cladeset>& tree_flat){
    tree_flat.insert(node->clade);
    for (set<vert_edge>::iterator child = node->children.begin(); child != node->children.end(); ++child){
        if (child->start <= pos && child->end >= pos){
            get_tree_flat(child->to, pos, tree_flat);
        }
    }
}

bool has_clade_invalidates(unordered_set<cladeset>& flat, const cladeset& cl){
    for (unordered_set<cladeset>::iterator f = flat.begin(); f != flat.end(); ++f){
        if (compare_bitsets(*f, cl) == -1){
            return true;
        }
    }
    return false;
}

long int new_closest_l(arg_node* l_rm, arg_node* r){
    //fprintf(stderr, "Mncl\n");
    long int closest_l = -1;
    for (vector<rc_edge>::iterator to = r->edges_to.begin(); to != r->edges_to.end(); ++to){
        if (to->final_dest != NULL && to->final_dest != l_rm && 
            (closest_l == -1 || *to->final_dest->sites.rbegin() > closest_l)){
            closest_l = *to->final_dest->sites.rbegin();
        }
    }
    for (vector<rc_edge>::iterator to = r->edges_to_solved.begin(); to != r->edges_to_solved.end(); ++to){
        if (to->final_dest != NULL && to->final_dest != l_rm &&
            (closest_l == -1 || *to->final_dest->sites.rbegin() > closest_l)){
            closest_l = *to->final_dest->sites.rbegin();
        }
    }
    for (vector<rc_edge>::iterator to = r->edges_to_unsolvable.begin(); to != r->edges_to_unsolvable.end(); ++to){
        if (to->final_dest != NULL && to->final_dest != l_rm &&
            (closest_l == -1 || *to->final_dest->sites.rbegin() > closest_l)){
            closest_l = *to->final_dest->sites.rbegin();
        }
    }
    return closest_l;
}

long int new_closest_r(arg_node* l, arg_node* r_rm){
   // fprintf(stderr, "Mncr\n");
    long int closest_r = -1;
    for (vector<rc_edge>::iterator from = l->edges_from.begin(); from != l->edges_from.end(); ++from){
        if (from->final_dest != NULL && from->final_dest != r_rm &&
            (closest_r == -1 || *from->final_dest->sites.begin() < closest_r)){
            closest_r = *from->final_dest->sites.begin();
        }
    }
    for (vector<rc_edge>::iterator from = l->edges_from_solved.begin(); from != l->edges_from_solved.end(); ++from){
        if (from->final_dest != NULL && from->final_dest != r_rm &&
            (closest_r == -1 || *from->final_dest->sites.begin() < closest_r)){
            closest_r = *from->final_dest->sites.begin();
        }
    }
    for (vector<rc_edge>::iterator from = l->edges_from_unsolvable.begin(); from != l->edges_from_unsolvable.end(); ++from){
        if (from->final_dest != NULL && from->final_dest != r_rm &&
            (closest_r == -1 || *from->final_dest->sites.begin() < closest_r)){
            closest_r = *from->final_dest->sites.begin();
        }
    }
    return closest_r;
}

bool is_fixed_right(arg_node* node, long int pos){
    //fprintf(stderr, "Mifr\n");
    if (node->fixed_right){
        return true;
    }
    else{
        for (vector<rc_edge>::iterator from = node->edges_from.begin();
            from != node->edges_from.end(); ++from){
            if (from->final_dest != NULL && from->final_dest->start <= pos &&
                from->final_dest->end >= pos){
                node->fixed_right = true;
                return true;
            }
        }
        for (vector<rc_edge>::iterator from = node->edges_from_solved.begin();
            from != node->edges_from_solved.end(); ++from){
            if (from->final_dest != NULL && from->final_dest->start <= pos &&
                from->final_dest->end >= pos &&
                compare_bitsets(node->clade, from->final_dest->clade) == -1){
                node->fixed_right = true;
                return true;
            }
        }
        for (vector<rc_edge>::iterator from = node->edges_from_unsolvable.begin();
            from != node->edges_from_unsolvable.end(); ++from){
            if (from->final_dest != NULL && from->final_dest-> start <= pos &&
                from->final_dest->end >= pos){
                node->fixed_right = true;
                return true;
            }
        }
    }
    return false;
}

bool is_fixed_left(arg_node* node, arg_sitemap& sites_pos){
    //fprintf(stderr, "Mifl\n");
    if (node->fixed_left){
        return true;
    }
    else{
        // Get pos before node
        long int pos = -1;
        for (arg_sitemap::reverse_iterator sp = 
            arg_sitemap::reverse_iterator(sites_pos.equal_range(*node->sites.begin()).first);
            sp != sites_pos.rend(); ++sp){
            if (sp->first < node->start){
                pos = sp->first;
                break;
            }
        }
        
        if (pos != -1){
            for (vector<rc_edge>::iterator to = node->edges_to.begin(); to != node->edges_to.end();
                ++to){
                if (to->final_dest != NULL && to->final_dest->start <= pos &&
                    to->final_dest->end >= pos){
                    node->fixed_left = true;
                    return true;
                }
            }
            for (vector<rc_edge>::iterator to = node->edges_to_solved.begin();
                to != node->edges_to_solved.end(); ++to){
                if (to->final_dest != NULL && to->final_dest->start <= pos &&
                    to->final_dest->end >= pos &&
                    compare_bitsets(node->clade, to->final_dest->clade) == -1){
                    node->fixed_left = true;
                    return true;   
                }
            }
            for (vector<rc_edge>::iterator to = node->edges_to_unsolvable.begin();
                to != node->edges_to_unsolvable.end(); ++to){
                if (to->final_dest != NULL && to->final_dest->start <= pos &&
                    to->final_dest->end >= pos){
                    node->fixed_left = true;
                    return true;
                }
            }
        }
    }
    return false;
}

void serialize_arg_tree(arg_node* node, 
    long int pos, 
    long int nextpos, 
    long int prop_dist, 
    gzFile& out, 
    arg_sitemap& sites_pos, 
    arg_clademap& sites_clade,
    arg_clademap& lv_grps, 
    arg_node* root, 
    short expdir,
    int num_haplotypes){
    //fprintf(stderr, "Msat\n");
    
    // Write out a placeholder for the node name.
    string name = "";
    serialize_str(out, name);
    
    // Write out a placeholder for the clade's persistence.
    long int persistence = -1;
    serialize_longint(out, persistence);
    
    // CHANGE THIS BACK TO 0 EVENTUALLY
    float brlen = 0.0;
    set<long int> muts_in_range;
    long int start = -1;
    long int end = -1;
    if (node == root){
        start = max(pos - prop_dist, (long int)1);
        end = pos + prop_dist;
    }
    else{
        // Figure out parent at current site.
        for (set<vert_edge>::iterator parent = node->parents.begin();
            parent != node->parents.end(); ++parent){
            if (parent->start <= pos && parent->end >= pos){
                start = parent->start;
                end = parent->end;
            }
        }
    }
    // Limit to propagation distance.
    if (start < pos - prop_dist){
        start = pos - prop_dist;
    }
    if (end > pos + prop_dist){
        end = pos + prop_dist;
    }
    
    // DELETE LATER
    //start = pos - prop_dist;
    //end = pos + prop_dist;
    
    set<long int> muts_in_range2;
    
    for (set<long int>::iterator site = node->mutations.begin(); site != node->mutations.end(); ++site){
    //for (set<long int>::iterator site = node->sites.begin(); site != node->sites.end(); ++site){
        if (*site >= start && *site <= end){
            brlen++;
            // Save space by not storing sites for root
            if (node != root){
                muts_in_range.insert(*site);
            }
        }
        if (*site >= pos - prop_dist && *site <= pos + prop_dist){
            muts_in_range2.insert(*site);
        }
        else if (*site > pos + prop_dist){
            break;
        }
    }
    
    /*
    for (set<long int>::iterator site = node->sites.begin(); site != node->sites.end(); ++site){
        if (*site >= start && *site <= end && muts_in_range2.find(*site) == muts_in_range2.end()){
            sites_in_range.insert(*site);
        }
        else if (*site > pos + prop_dist){
            break;
        }
    }
    */
    
    float brlen_lv = 0;
    for (set<long int>::iterator site = node->sites_lvgrp.begin(); site != node->sites_lvgrp.end();
        ++site){
        if (*site >= pos - prop_dist && *site <= pos + prop_dist){
            brlen_lv += 0.5;
        }
    } 
    
    // Normalize by genomic span of the parent edge at this site
    brlen = brlen/(float)(end-start+1);
    serialize_float(out, brlen);
    
    // Store another version normalized by the full range of the node
    float brlen2 = (float)muts_in_range2.size()/(float)(min(node->end, pos + prop_dist) - max(max(node->start, pos - prop_dist), (long int)1)+1);
    serialize_float(out, brlen2);
    
    // Store a "branch length" of leaving group sites normalized by the full range of the node
    float brlen_lv_norm = brlen_lv/(float)(min(node->end, pos + prop_dist) - max(max(node->start, pos - prop_dist), (long int)1)+1);
    serialize_float(out, brlen_lv_norm);
    
    serialize_longint(out, muts_in_range2.size());
    for (set<long int>::iterator mut = muts_in_range2.begin(); mut != muts_in_range2.end(); ++mut){
        serialize_longint(out, *mut);
    }
    
    /*
    serialize_longint(out, sites_in_range.size());
    for (set<long int>::iterator site = sites_in_range.begin(); site != sites_in_range.end(); ++site){
        serialize_longint(out, *site);
    }
    */

    // Get leaves attached to this node (its bitset minus all child bitsets)
    cladeset this_leaves = node->clade;
    
    
    map<arg_node*, long int> new_edge_indices;
    map<arg_node*, long int> new_edge_anchors;
    
    short expdir_this = 0;
    
    long int cr = -1;
    
    // Count up children overlapping position
    set<arg_node*> node_children;
    int numchildren = 0;
    for (set<vert_edge>::iterator child = node->children.begin(); child != node->children.end();
        ++child){
        if (child->start <= pos && child->end >= pos){
            this_leaves = set_diff_bitset(this_leaves, child->to->clade);
            node_children.insert(child->to);
            numchildren++;
        }
    }
    
    serialize_bitset(out, this_leaves, num_haplotypes);
    
    // Number of children
    serialize_int(out, numchildren);
    
    for (set<arg_node*>::iterator child = node_children.begin(); child != node_children.end();
        ++child){
        serialize_arg_tree(*child, pos, nextpos, prop_dist, out, sites_pos, sites_clade, lv_grps, root, expdir_this, num_haplotypes);
    }
}

void print_node(gzFile& main_out, 
    string& chrom, 
    arg_sitemap& sites_pos, 
    arg_clademap& sites_clade, 
    arg_clademap& lv_grps, 
    long int pos, 
    long int prop_dist, 
    arg_node* root, 
    long int nextpos,
    int num_haplotypes){
    //fprintf(stderr, "Mpn\n");
    // Don't create entries for "SNPs" that were just fixed derived alleles.
    if (root->mutations.find(pos) != root->mutations.end()){
        return;
    }
    
    static string prevchrom = "";
    if (chrom.compare(prevchrom) != 0){
        // Need to store string data.
        serialize_bool(main_out, true);
        serialize_str(main_out, chrom);
        prevchrom = chrom;
    }
    else{
        // No need to store string data.
        serialize_bool(main_out, false);
    }
    
    // Write position
    serialize_longint(main_out, pos);
    
    // Write ARG at position as tree.
    serialize_arg_tree(root, pos, nextpos, prop_dist, main_out, sites_pos, sites_clade,
        lv_grps, root, 0, num_haplotypes);

}

pair<long int, long int> inds_between(long int left, long int right, arg_sitemap& sites_pos){
    //fprintf(stderr, "Mib\n");
    // Find tentative midpoint between L and R clades.
    float midpt = ((float)right - (float)left)/2 + (float)left;
    long int newleft = -1;
    long int newright = -1;
    vector<long int> positions;
    arg_sitemap::iterator sp = sites_pos.find(left);
    if (sp != sites_pos.end()){
        while(sp->first == left && sp != sites_pos.end()){
            sp++;
        }
        for (; sp != sites_pos.end(); ++sp){
            if (sp->first > left && sp->first < right){
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
            int l_ind;
            int r_ind;
            if (positions.size() % 2 == 0){
                r_ind = (int)((float)positions.size()/2.0);
                l_ind = r_ind - 1;
            }
            else{
                // 2 ways to divide - choose randomly
                if (!RAND || (float)rand() / (float)RAND_MAX < 0.5){
                    l_ind = (int) ((float)positions.size()-1)/2.0;
                    r_ind = l_ind + 1;
                }
                else{
                    r_ind = (int) ((float)positions.size()-1)/2.0;
                    l_ind = r_ind-1;
                }
            }
            
            newleft = positions[l_ind];
            newright = positions[r_ind];
        }
        else if (positions.size() == 1){
            // Randomly assign L & R
            if (!RAND || (float)rand() / (float)RAND_MAX < 0.5){
                newleft = left;
                newright = positions[0];
            }
            else{
                newleft = positions[0];
                newright = right;
            }
        }
    }
    return make_pair(newleft, newright);
}


void get_nodes_rec(arg_node* node,
    long int pos,
    long int nextpos, 
    long int prop_dist,
    set<arg_node*>& nodes_l,
    set<arg_node*>& nodes_l_all,
    set<arg_node*>& nodes_r,
    multimap<arg_node*, arg_node*>& rnodes_lv){
    //fprintf(stderr, "Mgnr\n");
    if (node->end < nextpos && node->recomb_right != -1){
        // Check to see if one of this node's recombination partners covers the site.
        
        bool covered_from_r = is_fixed_right(node, nextpos);

        if (!covered_from_r && node->recomb_right > nextpos && *node->sites.rbegin() + prop_dist >= nextpos){   
            nodes_l.insert(node);
        }
        nodes_l_all.insert(node);

        for (vector<rc_edge>::iterator from = node->edges_from.begin();
            from != node->edges_from.end(); ++from){
            if (from->final_dest != NULL && from->final_dest->start > nextpos){
                //(from->final_dest->recomb_left <= *node->sites.rbegin() ||
                //node->recomb_right == *from->final_dest->sites.begin())){
                //&& (from->final_dest->recomb_left == -1 || from->final_dest->recomb_left < nextpos) &&
                //nextpos >= *from->final_dest->sites.begin() - prop_dist){
                nodes_r.insert(from->final_dest);
                rnodes_lv.insert(make_pair(from->final_dest, from->node));
            }
        }
        for (vector<rc_edge>::iterator from = node->edges_from_solved.begin();
            from != node->edges_from_solved.end(); ++from){
            if (from->final_dest != NULL && from->final_dest->start > nextpos){
                //(from->final_dest->recomb_left <= *node->sites.rbegin() || 
                //node->recomb_right == *from->final_dest->sites.begin()){
                short comp = compare_bitsets(from->final_dest->clade, node->clade);
                if (comp == -1){
                    //&& (from->final_dest->recomb_left == -1 || from->final_dest->recomb_left < nextpos) &&
                    //nextpos >= *from->final_dest->sites.begin() - prop_dist){
                    nodes_r.insert(from->final_dest);
                    rnodes_lv.insert(make_pair(from->final_dest, from->node));
                }
            }
        }
        for (vector<rc_edge>::iterator from = node->edges_from_unsolvable.begin();
            from != node->edges_from_unsolvable.end(); ++from){
            if (from->final_dest != NULL && from->final_dest->start > nextpos){
                nodes_r.insert(from->final_dest);
            }
        }
    }
    
    for (set<vert_edge>::iterator child = node->children.begin(); child != node->children.end();
        ++child){
        if (child->start <= pos && child->end >= pos){
            get_nodes_rec(child->to, pos, nextpos, prop_dist, nodes_l, nodes_l_all, nodes_r, rnodes_lv);
        }
    }
}


void get_nodes_rec_nextpos(arg_node* node,
    long int pos,
    long int nextpos, 
    long int prop_dist,
    set<arg_node*>& nodes_l,
    set<arg_node*>& nodes_r,
    set<arg_node*>& lv_grps){
    //fprintf(stderr, "Mgnrnp\n");
    if (node->end < nextpos && node->recomb_right != -1){
        // Check to see if one of this node's recombination partners covers the site.
        
        bool covered_from_r = is_fixed_right(node, nextpos);

        if (covered_from_r){
            nodes_l.insert(node);
            
            for (vector<rc_edge>::iterator from = node->edges_from.begin();
                from != node->edges_from.end(); ++from){
                if (from->final_dest != NULL && from->final_dest->start <= nextpos && from->final_dest->end >= nextpos){
                    lv_grps.insert(from->node);
                    nodes_r.insert(from->final_dest);
                }
            }
            for (vector<rc_edge>::iterator from = node->edges_from_solved.begin();
                from != node->edges_from_solved.end(); ++from){
                if (from->final_dest != NULL && from->final_dest->start <= nextpos && from->final_dest->end >= nextpos){
                    short comp = compare_bitsets(from->final_dest->clade, node->clade);
                    if (comp == -1){
                        lv_grps.insert(from->node);
                        nodes_r.insert(from->final_dest);
                    }
                }
            }
            for (vector<rc_edge>::iterator from = node->edges_from_unsolvable.begin();
                from != node->edges_from_unsolvable.end(); ++from){
                if (from->final_dest != NULL && from->final_dest->start <= nextpos && from->final_dest->end >= nextpos){
                    nodes_r.insert(from->final_dest);
                }
            }
        }
    }
    
    for (set<vert_edge>::iterator child = node->children.begin(); child != node->children.end();
        ++child){
        if (child->start <= pos && child->end >= pos){
            get_nodes_rec_nextpos(child->to, pos, nextpos, prop_dist, nodes_l, nodes_r, lv_grps);
        }
    }
}

void fill_in_recomb(arg_node* root,
    long int pos,
    long int nextpos, 
    arg_sitemap& sites_pos,
    arg_clademap& sites_clade, 
    arg_clademap& lv_grps,
    long int prop_dist,
    long int right_site,
    vector<arg_node*>& newnodes,
    int num_haplotypes){
    //fprintf(stderr, "Mfir\n");
    set<arg_node*> nodes_l;
    set<arg_node*> nodes_l_all;
    set<arg_node*> nodes_r;
    multimap<arg_node*, arg_node*> rnodes_lv;
    
    get_nodes_rec(root, pos, nextpos, prop_dist, nodes_l, nodes_l_all, nodes_r, rnodes_lv);
    
    for (set<arg_node*>::iterator r = nodes_r.begin(); r != nodes_r.end(); ++r){
        for (vector<rc_edge>::iterator to = (*r)->edges_to.begin(); to != (*r)->edges_to.end();
            ++to){
            if (to->final_dest != NULL && to->final_dest->end >= pos){
                nodes_l.insert(to->final_dest);
            }
        }
        for (vector<rc_edge>::iterator to = (*r)->edges_to_solved.begin(); to != (*r)->edges_to_solved.end();
            ++to){
            if (to->final_dest != NULL && to->final_dest->end >= pos){
                short comp = compare_bitsets(to->final_dest->clade, (*r)->clade);
                if (comp == -1){
                    nodes_l.insert(to->final_dest);
                }
            }
        }
        for (vector<rc_edge>::iterator to = (*r)->edges_to_unsolvable.begin(); to != (*r)->edges_to_unsolvable.end();
            ++to){
            if (to->final_dest != NULL && to->final_dest->end >= pos){
                nodes_l.insert(to->final_dest);
            }
        }
    }
    
    set<long int> sites_between;
    
    if (nodes_l.size() > 0 || nodes_r.size() > 0){
        
        if (DEBUG_MODE){
            fprintf(stderr, "L\n");
            for (set<arg_node*>::iterator l = nodes_l.begin(); l != nodes_l.end(); ++l){
                print_node_lite(*l, num_haplotypes);
            }
            fprintf(stderr, "R\n");
            for (set<arg_node*>::iterator r = nodes_r.begin(); r != nodes_r.end(); ++r){
                print_node_lite(*r, num_haplotypes);
                
            }
            fprintf(stderr, "\n");
        }
        
        if (nodes_r.size() > 0){
            // Determine rightmost boundary in R nodes
            long int rmost = -1;
            for (set<arg_node*>::iterator r = nodes_r.begin(); r != nodes_r.end(); ++r){
                if (*(*r)->sites.begin() > rmost){
                    rmost = *(*r)->sites.begin();
                }
            }
            
            // Determine all sites between L and R nodes
            for (arg_sitemap::iterator sp = sites_pos.equal_range(nextpos).second;
                sp != sites_pos.end(); ++sp){
                if (sp->first >= rmost){
                    break;
                }
                else{
                    sites_between.insert(sp->first);
                }
            }
            
        }
        
        sites_between.insert(nextpos);
        
        // Map destination sites to left nodes
        map<arg_node*, long int> lnode_dests;
        // Map destination sites to right nodes
        map<arg_node*, long int> rnode_dests;
        
        // Now try to fix things up at each site.
        for (set<long int>::iterator site = sites_between.begin(); site != sites_between.end(); ++site){
            
            if (DEBUG_MODE){
                fprintf(stderr, "\n=====SITE %ld=====\n", *site);
            }
            
            set<arg_node*> nodes_l_site;
            for (set<arg_node*>::iterator l = nodes_l.begin(); l != nodes_l.end(); ){
                bool erase = false;
                if ((*l)->start <= *site && (*l)->end < *site && (*l)->recomb_right > *site && 
                    *(*l)->sites.rbegin() + prop_dist >= *site){
                    if (!is_fixed_right(*l, *site)){
                        nodes_l_site.insert(*l);
                    }
                    else{
                        // It's out of the running for the rest of the sites.
                        erase = true;
                    }
                }
                if (erase){
                    nodes_l.erase(l++);
                }
                else{
                    ++l;
                }
            }
            
            set<arg_node*> nodes_r_site;
            for (set<arg_node*>::iterator r = nodes_r.begin(); r != nodes_r.end(); ++r){
                if ((*r)->start > *site && (*r)->recomb_left < *site && 
                    *(*r)->sites.begin() - prop_dist <= *site){
                    if (!is_fixed_left(*r, sites_pos)){
                        // See if node has recomb final dest in between this site and 
                        // its own end position
                        bool has_rec_at_site = false;
                        for (vector<rc_edge>::iterator to = (*r)->edges_to.begin();
                            to != (*r)->edges_to.end(); ++to){
                            if (to->final_dest != NULL && to->final_dest->end >= *site){
                                has_rec_at_site = true;
                                break;
                            }
                        }
                        if (!has_rec_at_site){
                            for (vector<rc_edge>::iterator to = (*r)->edges_to_solved.begin();
                                to != (*r)->edges_to_solved.end(); ++to){
                                if (to->final_dest != NULL &&
                                    to->final_dest->end >= *site &&
                                    compare_bitsets(to->final_dest->clade, (*r)->clade) == -1){
                                    has_rec_at_site = true;
                                    break;
                                }
                            }
                            if (!has_rec_at_site){
                                for (vector<rc_edge>::iterator to = (*r)->edges_to_unsolvable.begin();
                                    to != (*r)->edges_to_unsolvable.end(); ++to){
                                    if (to->final_dest != NULL &&
                                        to->final_dest->end >= *site){
                                        has_rec_at_site = true;
                                        break;
                                    }
                                }
                            }
                        }
                        if (!has_rec_at_site){
                            nodes_r_site.insert(*r);
                        }
                    }
                }
            }
            
            
            // Determine if we should expand left or right nodes.
            bool exp_l;
            
            if (nodes_r_site.size() == 0){
                // Let L nodes take over.
                exp_l = true;
            }
            else if (nodes_l_site.size() == 0){

                // Let R nodes take over.
                exp_l = false;
            }
            else{

                // Determine which is closest.
                long int l_maxr = -1;
                long int r_minl = -1;
                for (set<arg_node*>::iterator l = nodes_l_site.begin(); l != nodes_l_site.end(); ++l){
                    if (l_maxr == -1 || *(*l)->sites.rbegin() > l_maxr){
                        l_maxr = *(*l)->sites.rbegin();
                    }
                }
                for (set<arg_node*>::iterator r = nodes_r_site.begin(); r != nodes_r_site.end(); ++r){
                    if (r_minl == -1 || *(*r)->sites.begin() < r_minl){
                        r_minl = *(*r)->sites.begin();
                    }
                }
                
                if (*site - l_maxr < r_minl - *site){
                    exp_l = true;
                }
                else if (*site - l_maxr > r_minl - *site){
                    exp_l = false;
                }
                else{
                    // Flip a coin.
                    if (!RAND || (float)rand() / (float)RAND_MAX < 0.5){
                        exp_l = true;
                    }
                    else{
                        exp_l = false;
                    }
                }
            }
            if (exp_l){
                for (set<arg_node*>::iterator l = nodes_l_site.begin(); l != nodes_l_site.end(); ++l){

                    if (DEBUG_MODE){
                        fprintf(stderr, "expand L -> R\n");
                        print_node_lite(*l, num_haplotypes);
                        print_recombs_right(*l, num_haplotypes);
                        fprintf(stderr, "= %ld =\n", (*l)->edges_from_solved.size());
                    }
                    
                    lnode_dests[*l] = *site;
                }
            }
            else{
                for (set<arg_node*>::iterator r = nodes_r_site.begin(); r != nodes_r_site.end(); ++r){

                    if (DEBUG_MODE){
                        fprintf(stderr, "expand L <- R\n");
                        print_node_lite(*r, num_haplotypes);
                    }
                    
                    rnode_dests[*r] = *site;
                    
                    // If we expanded a right node left, we don't have to keep checking
                    // it, since it's now at its max range.
                    (*r)->fixed_left = true;
                    nodes_r.erase(nodes_r.find(*r));
                    // We can also erase all L nodes that are tied to it.
                    for (vector<rc_edge>::iterator to = (*r)->edges_to.begin(); to != (*r)->edges_to.end(); ++to){
                        if (to->final_dest != NULL && nodes_l.find(to->final_dest) != nodes_l.end()){
                            nodes_l.erase(nodes_l.find(to->final_dest));
                        }
                    }
                    for (vector<rc_edge>::iterator to = (*r)->edges_to_solved.begin(); to != (*r)->edges_to_solved.end(); ++to){
                        if (to->final_dest != NULL && nodes_l.find(to->final_dest) != nodes_l.end()){
                            nodes_l.erase(nodes_l.find(to->final_dest));
                        }
                    }
                    for (vector<rc_edge>::iterator to = (*r)->edges_to_unsolvable.begin();
                        to != (*r)->edges_to_unsolvable.end(); ++to){
                        if (to->final_dest != NULL && nodes_l.find(to->final_dest) != nodes_l.end()){
                            nodes_l.erase(nodes_l.find(to->final_dest));
                        }
                    }
                    if (nodes_l.size() == 0){
                        break;
                    }
                }
            }
        }
        for (map<arg_node*, long int>::iterator lnd = lnode_dests.begin(); lnd != lnode_dests.end();
            ++lnd){
            expand_range_fixed(lnd->first, lnd->second, sites_pos, sites_clade, lv_grps, prop_dist,
                root, newnodes, false, num_haplotypes, false, nodes_l);
            /*
            long int origend = lnd->first->end;
            lnd->first->end = lnd->second;                   
            expand_range_newsites(lnd->first, lnd->first->end-origend, sites_pos, sites_clade,
                lv_grps, prop_dist, false, root, newnodes, num_haplotypes);
            */
        }
        for (map<arg_node*, long int>::iterator rnd = rnode_dests.begin(); rnd != rnode_dests.end();
            ++rnd){
            expand_range_fixed(rnd->first, rnd->second, sites_pos, sites_clade, lv_grps, prop_dist,
                root, newnodes, false, num_haplotypes, false, nodes_r);
            /*
            long int origstart = rnd->first->start;
            rnd->first->start = rnd->second;
            expand_range_newsites(rnd->first, rnd->first->start-origstart, sites_pos, sites_clade,
                lv_grps, prop_dist, false, root, newnodes, num_haplotypes);
            */
        }
    }
    
}

/**
 * Builds and compares nodes, creating relationships between them.
 */
void compare_all_sites(
    deque<cladeset >& site_bitsets,
    arg_sitemap& sites_pos,
    arg_clademap& sites_clade,
    arg_clademap& lv_grps,
    long int &rightIndex, 
    gzFile &infile, 
    long int bufsize,
    const unsigned int num_haplotypes, 
    long unsigned int &num_sites, 
    bool &eof,
    long int prop_dist,
    FILE* recomb_out,
    arg_node* root,
    gzFile& main_out,
    vector<loc>& locs){
    
    long unsigned int newsites;
    
    vector<recomb_event> recomb_catalog;
    
    long int progress = 100;
    long int prev = 0;
    
    if (num_sites < 2){
        newsites = expand_matrix(site_bitsets, 0,
            infile, bufsize, num_haplotypes, eof, prop_dist);
        num_sites += newsites;
    }
    
    bool matrix_end_reached = false;
    
    long int last_site_printed = 0;
    
    long int last_deleted = 0;
    
    vector<arg_node*> newnodes;
    
    long int bitsets_popped = 0;
    
    string prevchrom;
    
    while(!matrix_end_reached){
    
        // Grab the chromosome position of the right index.
        long int right_site = locs[rightIndex].pos;
        string curchrom = locs[rightIndex].chrom;
        
        // Track propagation limits.
        long int right_start = max(right_site - prop_dist, (long int) 1);
        long int right_end = min(right_site + prop_dist, locs[locs.size()-1].pos);
        
        cladeset right_clade = site_bitsets[0];
        
        set<arg_node*> dummy;
        
        // Skip all-ancestral sites
        if (right_clade.count() > 0){
            
            if (right_clade.count() == num_haplotypes){
                root->mutations.insert(right_site);
            }
            else{
                
                arg_node* newnode = new arg_node();
                newnode->sites.insert(right_site);
                newnode->mutations.insert(right_site);
                newnode->start = max(right_site - prop_dist, (long int) 1);
                newnode->end = right_site + prop_dist;
                newnode->recomb_left = -1;
                newnode->recomb_right = -1;
                newnode->site_support = true;
                newnode->is_leaving_grp = false;
                newnode->skip = false;
                newnode->clade = right_clade;
                newnode->fixed_left = false;
                newnode->fixed_right = false;
                
                if (DEBUG_MODE){
                    fprintf(stderr, "CREATING %ld\n", right_site);
                }
                arg_node* stored = create_node(newnode, sites_pos, sites_clade, 
                    lv_grps, prop_dist, root, newnodes, num_haplotypes, dummy);
                
                if (stored != newnode){
                    delete newnode;
                }
                
            }

            if (DEBUG_MODE){
                check_solved_connections(sites_pos, num_haplotypes);
                check_everything(root, root, sites_pos, sites_clade, num_haplotypes);
                check_cp_edges(sites_pos, num_haplotypes);
                check_solved_connections(sites_pos, num_haplotypes);
            }
        }
        
        // Create any new nodes.
        bool all_checked = false;
        
        
        //if (true){
        while(newnodes.size() > 0 && !all_checked){
            vector<arg_node*> nncpy;
            for (vector<arg_node*>::iterator nn = newnodes.begin(); nn != newnodes.end();){
                
                bool deleted = true;
                
                // Eliminate any sites too far from the current site to matter.
                for (set<long int>::iterator site = (*nn)->sites.begin(); site != (*nn)->sites.end();){
                    if (*site <= last_site_printed){
                        if ((*nn)->mutations.find(*site) != (*nn)->mutations.end()){
                            (*nn)->mutations.erase((*nn)->mutations.find(*site));
                        }
                        (*nn)->sites.erase(site++);
                    }
                    else{
                        ++site;
                    }
                }
                
                if ((*nn)->sites.size() > 0 && 
                    *(*nn)->sites.rbegin() < right_site - 2*prop_dist &&
                    (*nn)->end > right_site - 4*prop_dist &&
                    safe_to_create_node((*nn)->clade, sites_pos, sites_clade, *(*nn)->sites.begin(),
                    prop_dist, num_haplotypes)){
                    
                    if ((*nn)->end < (*nn)->start){
                        fprintf(stderr, "nn end less than start\n");
                        print_node_lite(*nn, num_haplotypes);
                        exit(1);
                    }
                
                    if ((*nn)->start < last_site_printed){
                        (*nn)->start = last_site_printed; 
                    }
                    
                    arg_node* stored = create_node(*nn, sites_pos, sites_clade, 
                        lv_grps, prop_dist, root, nncpy, num_haplotypes, dummy);
                    if (stored != *nn){
                        delete *nn;
                    }
                    
                }
                else if ((*nn)->end >= right_site - 2*prop_dist){
                    // Don't create yet, but don't delete.
                    deleted = false;
                }
                else{
                    delete *nn;
                    deleted = true;
                }
                if (deleted){
                    newnodes.erase(nn);
                }
                else{
                    ++nn;
                }
            }
            all_checked = true;
            for (vector<arg_node*>::iterator nn_new = nncpy.begin(); nn_new != nncpy.end();){
                newnodes.push_back(*nn_new);
                nncpy.erase(nn_new);
                all_checked = false;
            }
        }
        
        // Solve all recombination events > 2x prop dist away from current site.
        // Sort in reverse order to make solving narrowest recomb events more
        // efficient.
        set<arg_node*> to_solve;
        
        for (arg_sitemap::iterator site_dat = sites_pos.begin(); site_dat != sites_pos.end(); ++site_dat){
            if (site_dat->first < right_site - (2*prop_dist)){
                if (site_dat->first > right_site - (4*prop_dist)){
                    if (site_dat->first == *site_dat->second->sites.rbegin()){
                        if (site_dat->second->edges_from.size() > 0){
                            to_solve.insert(site_dat->second);
                        }
                    }
                }
            }
            else{
                break;
            }
        }
        /*
        fprintf(stderr, "prelim to solve:\n");
        for (int x = 0; x < to_solve.size(); ++x){
            fprintf(stderr, "%d:\n", x);
            print_node_lite(to_solve[x], num_haplotypes);
        }
        */
        
        while (to_solve.size() > 0){
            
            arg_node* n = *to_solve.begin();
            
            bool recomb_solvable = n->edges_from.size() > 0;
            while(recomb_solvable){
                //fprintf(stderr, "attempting\n");
                if (DEBUG_MODE){
                    fprintf(stderr, "about to resolve recomb:\n");
                    print_node_lite(n, num_haplotypes);
                    print_recombs_right(n, num_haplotypes);
                }

                arg_node* try_other = NULL;
                
                recomb_solvable = resolve_recomb(n, lv_grps, 
                    sites_pos, sites_clade, prop_dist, root, 
                    true, try_other, recomb_catalog, 
                    newnodes, num_haplotypes, recomb_out, -1, -1, to_solve,
                    curchrom);

                while (!recomb_solvable && try_other != NULL && 
                    *try_other->sites.rbegin() < right_site - 2*prop_dist && 
                    *try_other->sites.rbegin() > right_site - 4*prop_dist){

                    recomb_solvable = resolve_recomb(try_other, lv_grps,
                        sites_pos, sites_clade, prop_dist, root,
                        true, try_other, recomb_catalog, 
                        newnodes, num_haplotypes, recomb_out, -1, -1, to_solve,
                        curchrom);
                    
                    
                }
                if (n->edges_from.size() == 0){
                    recomb_solvable = false;
                }
            }
            if (to_solve.find(n) != to_solve.end()){
                to_solve.erase(to_solve.find(n));
            }
        } 
        
        // Fill in recombination gaps > 4x prop_dist away from current site and
        // determine what sites need to be printed.
        set<long int> to_print;
        
        for (arg_sitemap::iterator site_dat = sites_pos.begin(); site_dat != sites_pos.end(); ++site_dat){
            if (site_dat->first <= last_site_printed || site_dat->first <= last_deleted){
                // Do nothing.
            }
            else if (site_dat->first < right_site - 4*prop_dist && site_dat->first >= last_site_printed){
                to_print.insert(site_dat->first);
            }
            else if (site_dat->first >= right_site - 4*prop_dist){
                break;
            }
        }
        
        
        long int printlen = to_print.size();
        
        for (set<long int>::iterator printsite = to_print.begin(); printsite != to_print.end();){
            long int nextpos = -1;
            for (arg_sitemap::iterator it = sites_pos.equal_range(*printsite).second; 
                it != sites_pos.end(); ++it){
                if (it->first > *printsite){
                    nextpos = it->first;
                    break;
                }
            }
            long int prevpos = -1;
            for (arg_sitemap::reverse_iterator it = arg_sitemap::reverse_iterator(sites_pos.equal_range(*printsite).first);
                it != sites_pos.rend(); ++it){
                if (it->first < *printsite){
                    prevpos = it->first;
                    break;
                }
            }
            
            fill_in_recomb(root, *printsite, nextpos, sites_pos, sites_clade,
                lv_grps, prop_dist, right_site, newnodes, num_haplotypes);
            
            
            print_node(main_out, locs[rightIndex].chrom, sites_pos, sites_clade, 
                lv_grps, *printsite, prop_dist, root, nextpos, num_haplotypes);
            last_site_printed = *printsite;
            to_print.erase(printsite++);
            
        }
        
        set<arg_node*> to_delete;
        set<arg_node*> to_adjust;
        map<arg_node*, long int> adj_inds;
        
        for (arg_sitemap::iterator site_dat = sites_pos.begin(); site_dat != sites_pos.end(); ){
            if (site_dat->first < right_site - (6*prop_dist)){
                
                if (DEBUG_MODE){
                    fprintf(stderr, "DELETING %ld\n", site_dat->first);
                    print_node_lite(site_dat->second, num_haplotypes);
                }
                
                last_deleted = site_dat->first;
                
                if (root->sites.find(site_dat->first) != root->sites.end()){
                    root->sites.erase(root->sites.find(site_dat->first));
                    if (root->mutations.find(site_dat->first) != root->mutations.end()){
                        root->mutations.erase(root->mutations.find(site_dat->first));
                    }
                }
                
                // Erase all nodes.

                if (site_dat->second->sites.size() > 0 &&
                    *(site_dat->second->sites.rbegin()) > site_dat->first){
                    if (DEBUG_MODE){
                        fprintf(stderr, "preserving...\n");
                        print_node_lite(site_dat->second, num_haplotypes);
                    }
                    // If this node is still needed by other sites, adjust its range and avoid deleting it.
                    
                    to_adjust.insert(site_dat->second);
                    
                    if (site_dat->second->sites.find(site_dat->first) != site_dat->second->sites.end()){
                        site_dat->second->sites.erase(site_dat->second->sites.find(site_dat->first));
                    }
                    if (site_dat->second->mutations.find(site_dat->first) != site_dat->second->mutations.end()){
                        site_dat->second->mutations.erase(site_dat->second->mutations.find(site_dat->first));
                    }
                    if (site_dat->second->sites_lvgrp.find(site_dat->first) != site_dat->second->sites_lvgrp.end()){
                        site_dat->second->sites_lvgrp.erase(site_dat->second->sites_lvgrp.find(site_dat->first));
                    }
                    long int prevstart = site_dat->second->start;
                    adj_inds[site_dat->second] = prevstart;
                    
                     if (site_dat->second->recomb_left != -1 && site_dat->second->recomb_left < site_dat->first){
                        site_dat->second->recomb_left = -1;
                    }
                    
                }
                else if (site_dat->second->sites.size() == 0 ||
                    site_dat->first == *site_dat->second->sites.rbegin()){
                    // Erase it.
                    if (to_adjust.find(site_dat->second) != to_adjust.end()){
                        to_adjust.erase(to_adjust.find(site_dat->second));
                    }
                    to_delete.insert(site_dat->second);
                }

                if (DEBUG_MODE){
                    fprintf(stderr, "erasing from data structure\n");
                }
                
                sites_pos.erase(site_dat++);
                if (DEBUG_MODE){
                    fprintf(stderr, "done\n");
                }
            }
            else{
                break;
            }
        }
        
        for (set<arg_node*>::iterator adj = to_adjust.begin(); adj != to_adjust.end();){

            (*adj)->start = max(*((*adj)->sites.begin()), (long int) 1);
            
            shorten_range(*adj, (*adj)->start - adj_inds[*adj], sites_pos, 
                sites_clade, lv_grps, prop_dist, num_haplotypes, root, newnodes, -1);
            
            to_adjust.erase(adj++);
        }
        
        set<arg_node*, nodeptr_child_sort> children_save;
        
        for (set<arg_node*>::iterator del = to_delete.begin(); del != to_delete.end();){
            if (DEBUG_MODE){
                fprintf(stderr, "erasing...\n");
            }
            
            // Wipe out connections to parents & children
            (*del)->start = -1;
            (*del)->end = -1;
            
            for (set<vert_edge>::iterator e = (*del)->parents.begin(); e != (*del)->parents.end();){
                
                for (set<vert_edge>::iterator e2 = e->to->children.begin(); e2 != e->to->children.end();){
                    if (e2->start == e->start && e2->end == e->end && e2->to == *del){
                        e->to->children.erase(e2++);
                        break;
                    }
                    else{
                        ++e2;
                    }
                }
                
                (*del)->parents.erase(e++);
            }
            for (set<vert_edge>::iterator e = (*del)->children.begin(); e != (*del)->children.end();){
                if (to_delete.find(e->to) == to_delete.end()){
                    children_save.insert(e->to);
                }
                for (set<vert_edge>::iterator e2 = e->to->parents.begin(); e2 != e->to->parents.end();){
                    if (e2->start == e->start && e2->end == e->end && e2->to == *del){
                        e->to->parents.erase(e2++);
                        break;
                    }
                    else{
                        ++e2;
                    }
                }
                (*del)->children.erase(e++);
            }
            
            // Make sure it's gone from sites_pos.
            for (set<long int>::iterator site = (*del)->sites.begin();
                site != (*del)->sites.end(); ++site){
                pair<arg_sitemap::iterator, arg_sitemap::iterator> er1 = 
                    sites_pos.equal_range(*site);
                for (arg_sitemap::iterator it = er1.first; it != er1.second;){
                    if (it->second == *del){
                        sites_pos.erase(it++);
                        break;
                    }
                    else{
                        ++it;
                    }
                }
            }
            
            // Erase from other data structure too.
            pair<arg_clademap::iterator, arg_clademap::iterator> er = sites_clade.equal_range((*del)->clade);
            for (arg_clademap::iterator it = er.first;
                it != er.second; ){
                if (it->second == *del){
                    sites_clade.erase(it++);
                    break;
                }
                else{
                    ++it;
                }
            }
            
            vector<pair<arg_node*, arg_node*> > nodepairs;
            for (vector<rc_edge>::iterator from = (*del)->edges_from.begin();
                from != (*del)->edges_from.end(); ++from){
                if (from->final_dest != NULL){
                    nodepairs.push_back(make_pair(*del, from->final_dest));
                }
            }
            for (vector<rc_edge>::iterator to = (*del)->edges_to.begin();
                to != (*del)->edges_to.end(); ++to){
                if (to->final_dest != NULL){
                    nodepairs.push_back(make_pair(to->final_dest, *del));
                }
            }
            set<arg_node*> dummy;
            for (vector<pair<arg_node*, arg_node*> >::iterator np = 
                nodepairs.begin(); np != nodepairs.end(); ++np){

                del_connections_between(np->first, np->second, lv_grps, dummy, num_haplotypes);
            }
            nodepairs.clear();
            if (to_adjust.find(*del) != to_adjust.end()){
                to_adjust.erase(to_adjust.find(*del));
            }
            if (adj_inds.count(*del) > 0){
                adj_inds.erase(*del);
            }
            for (vector<arg_node*>::iterator nn = newnodes.begin(); nn != newnodes.end();){
                if (*nn == *del){
                    newnodes.erase(nn);
                }
                else{
                    ++nn;
                }
            }
            delete *del;
            to_delete.erase(del++);
        }
        
        
        for (set<arg_node*, nodeptr_child_sort>::iterator child = children_save.begin(); child !=
            children_save.end(); ++child){            
            
            insert_below_parent(*child, root, sites_pos, sites_clade, lv_grps, prop_dist,
                (*child)->start, (*child)->end, root, newnodes, num_haplotypes);    
       
        }
        
        site_bitsets.pop_front();
        bitsets_popped++;
        
        ++rightIndex;
        
        // If the right index is larger than the size of the matrix, attempt to 
        // read in more data from the input file.
        if (site_bitsets.size() == 0){
            if (!eof){
                newsites = expand_matrix(site_bitsets, bitsets_popped,
                    infile, bufsize, num_haplotypes, eof, prop_dist);
                num_sites += newsites;
                
                if (newsites == 0){
                    // We didn't actually recover anything new from the file;
                    // bail out here.
                    --rightIndex;
                    matrix_end_reached = true;
                    eof = true;
                }
            }
            else{
                // We've reached the end of the file and the end of the matrix.
                --rightIndex;
                matrix_end_reached = true;
            }
        }
        
        if (!DEBUG_MODE && rightIndex % progress == 0){
            fprintf(stderr, "Processed %s %ld\r", locs[rightIndex-1].chrom.c_str(), locs[rightIndex-1].pos);
        }
        
        prevchrom = curchrom;
        
    }
    
    // === Clean up === //
    
    // Resolve all last recomb events.
    // Do this in two passes so we can try to "force solve" everything.
    for (int attempt = 0; attempt < 2; ++attempt){
        set<arg_node*> solvenodes;
        for (arg_sitemap::iterator site_dat = sites_pos.begin(); site_dat != sites_pos.end(); ++site_dat){
            if (site_dat->first > last_site_printed){
                bool recomb_solvable = site_dat->second->edges_from.size() > 0;
                if (recomb_solvable){
                    solvenodes.insert(site_dat->second);
                }
            }
        }
        while(solvenodes.size() > 0){
            arg_node* n = *solvenodes.begin();
            //for (set<arg_node*>::iterator n = solvenodes.begin(); n != solvenodes.end(); ++n){
            bool recomb_solvable = n->edges_from.size() > 0;
            while(recomb_solvable){
                
                if (DEBUG_MODE){
                    fprintf(stderr, "about to resolve recomb:\n");
                    print_node_lite(n, num_haplotypes);
                    print_recombs_right(n, num_haplotypes);
                }
                arg_node* try_other = NULL;

                recomb_solvable = resolve_recomb(n, lv_grps, 
                    sites_pos, sites_clade, prop_dist, root, 
                    true, try_other, recomb_catalog, 
                    newnodes, num_haplotypes, recomb_out, -1, -1, solvenodes,
                    prevchrom);
                while (!recomb_solvable && try_other != NULL){

                    recomb_solvable = resolve_recomb(try_other, lv_grps,
                        sites_pos, sites_clade, prop_dist, root,
                        true, try_other, recomb_catalog, 
                        newnodes, num_haplotypes, recomb_out, -1, -1, solvenodes,
                        prevchrom);
                }
                if (n->edges_from.size() == 0){
                    recomb_solvable = false;
                }
            }
            if (solvenodes.find(n) != solvenodes.end()){
                solvenodes.erase(solvenodes.find(n));
            }
        }
    }
    
    // Create any new nodes
    for (vector<arg_node*>::iterator nn = newnodes.begin(); nn != newnodes.end();){
        delete *nn;
        newnodes.erase(nn);
    }
    
    // Print out last sites.
    set<long int> printsites;
    for (arg_sitemap::iterator site_dat = sites_pos.begin(); site_dat != sites_pos.end(); ++site_dat){
        if (site_dat->first > last_site_printed){
            printsites.insert(site_dat->first);
        }
    }
    
    for (set<long int>::iterator printsite = printsites.begin(); printsite != printsites.end();){
        long int nextpos = -1;
        for (arg_sitemap::iterator it = sites_pos.equal_range(*printsite).second; 
            it != sites_pos.end(); ++it){
            if (it->first > *printsite){
                nextpos = it->first;
                break;
            }
        }
        
        fill_in_recomb(root, *printsite, nextpos, sites_pos, sites_clade,
            lv_grps, prop_dist, -1, newnodes, num_haplotypes);
                
        print_node(main_out, locs[rightIndex].chrom, sites_pos, sites_clade, lv_grps,
            *printsite, prop_dist, root, nextpos, num_haplotypes);
            
        last_site_printed = *printsite;
        printsites.erase(printsite++);
    }
    
    // Delete everything.
    set<arg_node*> to_delete;
    
    for (arg_sitemap::iterator site_dat = sites_pos.begin(); site_dat != sites_pos.end();){
        
        to_delete.insert(site_dat->second);
        
        sites_pos.erase(site_dat++);
    }
    
    for (set<arg_node*>::iterator del = to_delete.begin(); del != to_delete.end();){
        // Erase from other data structure too.
        pair<arg_clademap::iterator, arg_clademap::iterator> er = sites_clade.equal_range((*del)->clade);
        for (arg_clademap::iterator it = er.first;
            it != er.second; ){
            if (it->second == *del){
                sites_clade.erase(it++);
                break;
            }
            else{
                ++it;
            }
        }
        
        // Delete parent/child edges
        for (set<vert_edge>::iterator e = (*del)->parents.begin(); e != (*del)->parents.end();){
            for (set<vert_edge>::iterator e2 = e->to->children.begin(); e2 != e->to->children.end();){
                if (e2->start == e->start && e2->end == e->end && e2->to == *del){
                    e->to->children.erase(e2++);
                    break;
                }
                else{
                    ++e2;
                }
            }
            
            (*del)->parents.erase(e++);
        }
        for (set<vert_edge>::iterator e = (*del)->children.begin(); e != (*del)->children.end();){
            for (set<vert_edge>::iterator e2 = e->to->parents.begin(); e2 != e->to->parents.end();){
                if (e2->start == e->start && e2->end == e->end && e2->to == *del){
                    e->to->parents.erase(e2++);
                    break;
                }
                else{
                    ++e2;
                }
            }
            (*del)->children.erase(e++);
        }
        
        // Delete recombination edges
        vector<pair<arg_node*, arg_node*> > nodepairs;
        for (vector<rc_edge>::iterator from = (*del)->edges_from.begin();
            from != (*del)->edges_from.end(); ++from){
            if (from->final_dest != NULL){
                nodepairs.push_back(make_pair(*del, from->final_dest));
            }
        }
        for (vector<rc_edge>::iterator to = (*del)->edges_to.begin();
            to != (*del)->edges_to.end(); ++to){
            if (to->final_dest != NULL){
                nodepairs.push_back(make_pair(to->final_dest, *del));
            }
        }
        set<arg_node*> dummy;
        for (vector<pair<arg_node*, arg_node*> >::iterator np = 
            nodepairs.begin(); np != nodepairs.end(); ++np){
            del_connections_between(np->first, np->second, lv_grps, dummy, 
                num_haplotypes);
        }
        nodepairs.clear();
        
        delete *del;
        to_delete.erase(del++);
    }
    fprintf(stderr, "Processed %s %ld\r", locs[rightIndex-1].chrom.c_str(), locs[rightIndex-1].pos);
    
    fprintf(stderr, "\n");
}


/**
 * Print a help message to the terminal and exit.
 */
void help(int code){
   fprintf(stderr, "sarge [OPTIONS]\n");
   fprintf(stderr, "Main program. Takes input files and runs ARG inference. Outputs gzipped, \
binary data that can be viewed and analyzed by other programs.\n");
   fprintf(stderr, "[OPTIONS]:\n");
   fprintf(stderr, "    --sites -s a file listing all sites in the main input \
file, where each line is a chromosome and position (in bp), tab separated \n");
   fprintf(stderr, "    --bufsize -b The number of characters to read from the \
input file at a time\n");
   fprintf(stderr, "    --input -i The main input file; a matrix of genotype data \
where rows are sites and columns are individuals. It must contain 1s and 0s only, \
with no spaces\n");
    fprintf(stderr, "   --out -o The name of the output file to create. Output will be \
in gzipped binary format.\n");
    fprintf(stderr, "   --recombfile -r Name of file to use for writing recombination \
events (default = INPUTFILENAME.recomb).\n");
    fprintf(stderr, "   --bufsize -b The number of characters to read from the \
input file at a time\n");
    fprintf(stderr, "   --exclude_bed -e Provide a BED file containing genomic regions \
to exclude from ARG construction (i.e. CpG sites)\n");
    fprintf(stderr, "   --prop_dist -p A number of bases for which each clade should be \
allowed to extend on either side of any given site that tags it\n");
    fprintf(stderr, "   --noheader -H If you are going to concatenate multiple output files \
(i.e. this file is only part of a chromosome and not the first part), set this option to \
keep from printing a header. The output file can then be combined with other output files \
using zcat, provided that the first output file has a header.\n");
exit(code);
}

int main(int argc, char *argv[]) {    
    // Define arguments 
    // http://www.gnu.org/software/libc/manual/html_node/Getopt-Long-Options.html#Getopt-Long-Options
    // http://www.gnu.org/software/libc/manual/html_node/Getopt-Long-Option-Example.html#Getopt-Long-Option-Example
    // Fields for each argument: name, has_arg (values: no_argument, required_argument,
    //     optional_argument)
    // flag = int value to store flag for the option, or NULL if option is string
    // val = short name for string option, or NULL
    
    static struct option long_options[] = {
       {"sites", required_argument, 0, 's'},
       {"bufsize", optional_argument, 0, 'b'},
       {"input", required_argument, 0, 'i'},
       {"recombfile", optional_argument, 0, 'r'},
       {"exclude_bed", required_argument, 0, 'e'},
       {"prop_dist", required_argument, 0, 'p'},
       {"out", required_argument, 0, 'o'},
       {"noheader", no_argument, 0, 'H'},
       {0, 0, 0, 0} 
    };
    
    // Set default values
    //std::vector<loc> locs;
    //std::vector<std::string> indvs;
    //int bufsize = 1024;
    // Read 1 MB at a time by default
    long int bufsize = 1048576;
    //long int bufsize = 8388608;
    std::string infilename;
    
    int option_index = 0;
    int ch;
    
    string recombfilename;
    bool recomb_filename_given = false;
    
    string locsfile;
    bool locs_given = false;
    
    bool exclude_bed_given = false;
    string exclude_bed_name;
    
    long int prop_dist = 10000;
    
    string mainoutfn;
    
    vector<loc> locs;
    
    bool no_header = false;
    
    if (argc == 1){
        help(0);
    }
    while((ch = getopt_long(argc, argv, "s:b:i:r:e:p:o:H", long_options, &option_index )) != -1){
        switch(ch){
            case 0:
                // This option set a flag. No need to do anything here.
                break;
            case 's':
                locsfile = optarg;
                locs_given = true;
                //locs = parse_locs(optarg);
                break;
            case 'b':
                bufsize = atol(optarg);
                //bufsize = strtol(optarg, NULL, 0);
                break;
            case 'i':
                infilename = optarg;
                break;
            case 'r':
                recombfilename = optarg;
                recomb_filename_given = true;
                break;
            case 'e':
                exclude_bed_given = true;
                exclude_bed_name = optarg;
                break;
            case 'p':
                prop_dist = atol(optarg);
                break;
            case 'o':
                mainoutfn = optarg;
                break;
            case 'H':
                no_header = true;
                break;
            case '?':
                //help(0);
                break;
            default:
                help(0);
        }    
    }
    
    // Error check.
    if (bufsize <= 0){
        fprintf(stderr, "ERROR: invalid buffer size of %ld provided. Exiting.\n", bufsize);
        exit(1);
    }
    if (prop_dist <= 0){
        fprintf(stderr, "ERROR: you must provide a propagation distance greater than zero\n");
        exit(1);
    }
    if (mainoutfn.length() == 0){
        fprintf(stderr, "ERROR: you must provide an output file name.\n");
        exit(1);
    }
    
    map<string, set<pair<long int, long int> > > exclude_bed;
    if (exclude_bed_given){
        parse_exclude_bed(exclude_bed_name, exclude_bed);
    }
    
    if (locs_given){
        locs = parse_locs(locsfile, exclude_bed);
    }
    else{
        fprintf(stderr, "ERROR: you must provide a file listing sites in the input file.\n");
        exit(1);
    }
    
    // Get "base name" for input file. Other file names will be based on this.
    string infilebase;
    int suffixstart = infilename.find_last_of('.');
    if (suffixstart > 0){
        infilebase = infilename.substr(0, suffixstart);
        // Check to see if this was just ".gz" -- if so, need to strip off the
        // extension too.
        if (infilename.substr(suffixstart, infilename.length()).compare(".gz") == 0){
            suffixstart = infilebase.find_last_of('.');
            if (suffixstart > 0){
                infilebase = infilebase.substr(0, suffixstart);
            }
        }
    }
    else{
        infilebase = infilename;
    }
    
    // Set up recombination output file
    if (!recomb_filename_given){
        recombfilename = infilebase + ".recomb";
    }
    
    // Open main input file for reading (binary format since it only contains 1s and 0s)
    gzFile infile;
    infile = gzopen(infilename.c_str(), "r");
    if (!infile){
        fprintf(stderr, "ERROR: the input file %s does not exist. Exiting.\n", infilename.c_str());
        exit(1);
    }
    
    // Find out how many haplotypes there are by reading the first line of the input
    // file
    
    // Ensure that the buffer is big enough to read in at least one whole row of 
    // haplotype data from the input file at a time
    const unsigned int num_haplotypes = get_num_haplotypes(infile, bufsize);
    
    // Check that input genotype data file and input sites file contain the
    // same number of sites.
    int numsites = get_num_sites(infile, bufsize, num_haplotypes);
    
    fprintf(stderr, "Input genotype file contains %d sites\n", numsites);
    
    if (numsites != locs.size()){
        fprintf(stderr, "ERROR: read %d sites from genotype file and %ld from sites file. \
Files should contain the same number of sites.\n", numsites, locs.size());
        exit(1);
    }
    
    // Create root node, parent of all others
    arg_node* root = new arg_node();
    for (int i = 0; i < num_haplotypes; ++i){
        root->clade.set(num_haplotypes-i-1);
    }
        
    // Make it cover all positions
    root->start = 1;
    root->end = locs[locs.size()-1].pos + prop_dist;
    
    bool eof = false;
    
    FILE* recomb_out;
    recomb_out = fopen(recombfilename.c_str(), "w");
    if (!recomb_out){
        fprintf(stderr, "ERROR: unable to open file %s for writing. Exiting.\n", recombfilename.c_str());
        exit(1);
    }
    
    gzFile main_out;
    main_out = gzopen(mainoutfn.c_str(), "wb");
    if (!main_out){
        fprintf(stderr, "ERROR opening file %s for writing.\n", mainoutfn.c_str());
        exit(1);
    }
    
    if (!no_header){
        // In the long run, this should be removed. For now, keep it in the
        // output file for backward-compatibility
        cladeset mask;
        mask.reset();
        for (int i = 0; i < num_haplotypes; ++i){
            mask.set(num_haplotypes-i-1);
        }
        serialize_header(main_out, num_haplotypes, mask);
    }
    
    // We'll now make sure that bufsize is a multiple of num_haplotypes+1 so 
    // we never have to gzseek().
    if (bufsize < num_haplotypes + 1){
        bufsize = num_haplotypes + 1;
    }
    else{
        bufsize = (long int) round( (float) (num_haplotypes+1) * round((float) bufsize / (float) (num_haplotypes+1)));
    }
    
    // Init random number seed
    srand(time(NULL));
    //srand(1);
    
    arg_sitemap sites_pos;
    arg_clademap sites_clade;
    arg_clademap lv_grps;
    
    deque<cladeset > site_bitsets;
    
    long int rightIndex = 0;
    
    // Read in the first chunk of the haplotype matrix.
    long unsigned int num_sites = expand_matrix(site_bitsets, 0, infile, bufsize,
        num_haplotypes, eof, prop_dist);
    
    // Build graph, reading more data & writing out data as it goes.
    compare_all_sites(site_bitsets, sites_pos, sites_clade, lv_grps, rightIndex,
        infile, bufsize, num_haplotypes, num_sites, eof,
        prop_dist, recomb_out, root, main_out, locs);
    
    fclose(recomb_out);
    gzclose(infile);
    gzclose(main_out);
    
    return 0;
}
