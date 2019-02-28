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
#include <fstream>
#include <bitset>
#include "sets.h"
#include "treeNode.h"
#include "serialize.h"
#include "common.h"
#include "argnode.h"
#include "debug.h"

using std::cout;
using std::endl;
using namespace std;

void check_for_holes(arg_node* node, arg_node* root, int num_haplotypes){
    if (node == root){
        return;
    }
    
    // See if there's any piece of the node's range that's not covered by parents.
    long int prog = node->start;

    for (set<vert_edge>::iterator p = node->parents.begin(); p != node->parents.end(); ++p){
        if (p->start > prog ){
            fprintf(stderr, "\n");
            fprintf(stderr, "hole in parents %ld %ld\n", prog, p->start);
            //print_node_upward_onelevel(node, 0, num_haplotypes);
            print_node_lite(node, num_haplotypes);  
            exit(1);
        }
        prog = p->end + 1;
    }
    if (prog < node->end){
        fprintf(stderr, "hole in parents? %ld %ld\n", prog, node->end);
        //print_node_upward_onelevel(node, 0, num_haplotypes);
        print_node_lite(node, num_haplotypes);
        exit(1);
    }
}

void check_ov_parents(arg_node* node, int num_haplotypes){
    for (set<vert_edge>::iterator p = node->parents.begin(); p != node->parents.end(); ++p){
        for (set<vert_edge>::iterator p2 = node->parents.begin(); p2 != node->parents.end(); ++p2){
            if (p != p2){
                if (ranges_overlap(p->start, p->end, p2->start, p2->end) && p->to == p2->to){
                    fprintf(stderr, "overlapping range in parents\n");
                    fprintf(stderr, "<%ld %ld> <%ld %ld>\n", p->start, p->end, p2->start, p2->end);
                    print_node_lite(p->to, num_haplotypes);
                    print_node_upward_onelevel(node, 0, num_haplotypes);
                    exit(1);
                }
            }
        }
    }
}

void check_ov_children(arg_node* node, int num_haplotypes){
    for (set<vert_edge>::iterator c = node->children.begin(); c != node->children.end(); ++c){
        for (set<vert_edge>::iterator c2 = node->children.begin(); c2 != node->children.end(); ++c2){
            if (c != c2){
                if (ranges_overlap(c->start, c->end, c2->start, c2->end) && c->to == c2->to){
                    fprintf(stderr, "overlapping range in children\n");
                    fprintf(stderr, "<%ld %ld> <%ld %ld>\n", c->start, c->end, c2->start, c2->end);
                    print_node_lite(c->to, num_haplotypes);
                    print_node_downward_onelevel(node, 0, num_haplotypes);
                    exit(1);
                }
            }
        }
    }
}

void check_solved_connected(arg_sitemap& sites_pos, long int left, long int right, int num_haplotypes){
    arg_node* lnode = NULL;
    arg_node* rnode = NULL;
    for (arg_sitemap::iterator it = sites_pos.equal_range(left).first;
        it != sites_pos.equal_range(left).second; ++it){
        if (it->second->mutations.find(left) != it->second->mutations.end()){
            lnode = it->second;
            break;
        }
    }
    for (arg_sitemap::iterator it = sites_pos.equal_range(right).first;
        it != sites_pos.equal_range(right).second; ++it){
        if (it->second->mutations.find(right) != it->second->mutations.end()){
            rnode = it->second;
            break;
        }
    }
    if (lnode != NULL && rnode != NULL){
        bool found_l = false;
        bool found_r = false;
        for (vector<rc_edge>::iterator from = lnode->edges_from_solved.begin();
            from != lnode->edges_from_solved.end(); ++from){
            if (from->final_dest != NULL && from->final_dest == rnode){
                found_l = true;
                break;
            }
        }
        for (vector<rc_edge>::iterator to = rnode->edges_to_solved.begin();
            to != rnode->edges_to_solved.end(); ++to){
            if (to->final_dest != NULL && to->final_dest == lnode){
                found_r = true;
                break;
            }
        }
        if (!found_l || !found_r){
            fprintf(stderr, "FOUND L %d R %d\n", found_l, found_r);
            print_node_lite(lnode, num_haplotypes);
            print_recombs_right_solved(lnode, num_haplotypes);
            fprintf(stderr, "\n");
            print_node_lite(rnode, num_haplotypes);
            print_recombs_left_solved(rnode, num_haplotypes);
            exit(1);
        }
    }
}

void check_solved_connections(arg_sitemap& sites_pos, int num_haplotypes){
    set<arg_node*> allnodes;
    for (arg_sitemap::iterator it = sites_pos.begin(); it != sites_pos.end(); ++it){
        allnodes.insert(it->second);
    }
    for (set<arg_node*>::iterator it = allnodes.begin(); it != allnodes.end(); ++it){
        for (vector<rc_edge>::iterator from = (*it)->edges_from_solved.begin();
            from != (*it)->edges_from_solved.end(); ++from){
            bool conn1 = false;
            for (vector<rc_edge>::iterator to = from->node->edges_to_solved.begin();
                to != from->node->edges_to_solved.end(); ++to){
                //if (to->node == *it){
                if (to->node == (*it) && to->type == from->type){
                    conn1 = true;
                    break;
                }
            }
            if (!conn1){
                fprintf(stderr, "L -> lv solved; missing L <- lv solved\n");
                print_node_lite(*it, num_haplotypes);
                print_node_lite(from->node, num_haplotypes);
                exit(1);
            }
            if (from->final_dest != NULL){
                bool connected = false;
                for (vector<rc_edge>::iterator to = from->final_dest->edges_to_solved.begin();
                    to != from->final_dest->edges_to_solved.end(); ++to){
                    if (to->final_dest != NULL && to->final_dest == *it){
                        connected = true;
                        break;
                    }
                }
                if (!connected){
                    fprintf(stderr, "L -> R solved; missing L <- R solved\n");
                    print_node_lite(*it, num_haplotypes);
                    print_node_lite(from->final_dest, num_haplotypes);
                    exit(1);
                }
                bool conn3 = false;
                for (vector<rc_edge>::iterator from2 = from->node->edges_from_solved.begin();
                    from2 != from->node->edges_from_solved.end(); ++from2){
                    if (from2->type == from->type && from2->node == from->final_dest){
                        conn3 = true;
                        break;
                    }   
                }
                if (!conn3){
                    fprintf(stderr, "L -> R solved; missing lv -> R solved\n");
                    print_node_lite(*it, num_haplotypes);
                    print_node_lite(from->final_dest, num_haplotypes);
                    print_node_lite(from->node, num_haplotypes);
                    fprintf(stderr, "=====\n");
                    print_node_lite(*it, num_haplotypes);
                    print_recombs_right_solved(*it, num_haplotypes);
                    fprintf(stderr, "\n");
                    print_node_lite(from->node, num_haplotypes);
                    print_recombs_right_solved(from->node, num_haplotypes);
                    fprintf(stderr, "\n");
                    print_node_lite(from->final_dest, num_haplotypes);
                    print_recombs_left_solved(from->final_dest, num_haplotypes);
                    exit(1);
                }
            }
        }
        for (vector<rc_edge>::iterator to = (*it)->edges_to_solved.begin();
            to != (*it)->edges_to_solved.end(); ++to){
            bool conn1 = false;
            for (vector<rc_edge>::iterator from = to->node->edges_from_solved.begin();
                from != to->node->edges_from_solved.end(); ++from){
                //if (from->node == *it){
                if (from->node == *it && from->type == to->type){
                    conn1 = true;
                    break;
                }
            }
            if (!conn1){
                fprintf(stderr, "lv <- R solved; missing lv -> R solved\n");
                
                if (to->final_dest != NULL){
                    fprintf(stderr, "fd:\n");
                    print_node_lite(to->final_dest, num_haplotypes);
                    print_recombs_right_solved(to->final_dest, num_haplotypes);
                }
                fprintf(stderr, "\nnode:\n");
                print_node_lite(*it, num_haplotypes);
                print_recombs_left_solved(*it, num_haplotypes);
                fprintf(stderr, "\n");
                fprintf(stderr, "lv:\n");
                print_node_lite(to->node, num_haplotypes);
                print_recombs_left_solved(to->node, num_haplotypes);
                print_recombs_right_solved(to->node, num_haplotypes);
                exit(1);
            }
            if (to->final_dest != NULL){
                bool connected = false;
                for (vector<rc_edge>::iterator from = to->final_dest->edges_from_solved.begin();
                    from != to->final_dest->edges_from_solved.end(); ++from){
                    if (from->final_dest != NULL && from->final_dest == *it){
                        connected = true;
                        break;
                    }
                }
                if (!connected){
                    fprintf(stderr, "L <- R solved; missing L -> R solved\n");
                    print_node_lite(*it, num_haplotypes);
                    print_node_lite(to->final_dest, num_haplotypes);
                    exit(1);
                }
                bool conn3 = false;
                for (vector<rc_edge>::iterator to2 = to->node->edges_to_solved.begin();
                    to2 != to->node->edges_to_solved.end(); ++to2){
                    if (to2->type == to->type && to2->node == to->final_dest){
                        conn3 = true;
                        break;
                    }
                }
                if (!conn3){
                    fprintf(stderr, "L <- R solved; missing L <- lv solved\n");
                    print_node_lite(*it, num_haplotypes);
                    print_node_lite(to->node, num_haplotypes);
                    exit(1);
                }
            }
        }
        for (vector<rc_edge>::iterator from = (*it)->edges_from_unsolvable.begin();
            from != (*it)->edges_from_unsolvable.end(); ++from){
            bool found = false;
            for (vector<rc_edge>::iterator to = from->final_dest->edges_to_unsolvable.begin();
                to != from->final_dest->edges_to_unsolvable.end(); ++to){
                if (to->final_dest == *it){
                    found = true;
                    break;
                }
            }
            if (!found){
                fprintf(stderr, "unsolved L -> R -- L <- R not found\n");
                print_node_lite(*it, num_haplotypes);
                print_node_lite(from->final_dest, num_haplotypes);
                exit(1);
            }
        }
        for (vector<rc_edge>::iterator to = (*it)->edges_to_unsolvable.begin();
            to != (*it)->edges_to_unsolvable.end(); ++to){
            bool found = false;
            for (vector<rc_edge>::iterator from = to->final_dest->edges_from_unsolvable.begin();
                from != to->final_dest->edges_to_unsolvable.end(); ++from){
                if (from->final_dest == *it){
                    found = true;
                    break;
                }
            }
            if (!found){
                fprintf(stderr, "unsolved L <- R -- L -> R not found\n");
                print_node_lite(to->final_dest, num_haplotypes);
                print_node_lite(*it, num_haplotypes);
                exit(1);
            }
        }
    }
}

void check_cp_edges(arg_sitemap& sites_pos, int num_haplotypes){
    set<arg_node*> allnodes;
    for (arg_sitemap::iterator it = sites_pos.begin(); it != sites_pos.end(); ++it){
        allnodes.insert(it->second);
    }
    for (set<arg_node*>::iterator it = allnodes.begin(); it != allnodes.end(); ++it){
        for (set<vert_edge>::iterator p = (*it)->parents.begin(); p != (*it)->parents.end(); ++p){
            if (p->to == *it){
                fprintf(stderr, "ERROR: node has parent edge to itself!\n");
                print_node_lite(*it, num_haplotypes);
                exit(1);
            }
        }
        for (set<vert_edge>::iterator c = (*it)->children.begin(); c != (*it)->children.end(); ++c){
            if (c->to == *it){
                fprintf(stderr, "ERROR: node has child edge to itself!\n");
                print_node_lite(*it, num_haplotypes);
                exit(1);
            }
        }
    }
}

void check_everything(arg_node* n, arg_node* root, arg_sitemap& sites_pos,
    arg_clademap& sites_clade, int num_haplotypes){
    check_for_holes(n, root, num_haplotypes);
    check_ov_parents(n, num_haplotypes);
    if (n != root && n->parents.size() == 0){
        fprintf(stderr, "NO PARENTS\n");
        print_node_recursive(n, 0, num_haplotypes);
        exit(1);
    }
    if (n != root){
        for (set<long int>::iterator site = n->sites.begin(); site != n->sites.end(); ++site){
            pair<arg_sitemap::iterator, arg_sitemap::iterator> er = sites_pos.equal_range(*site);
            bool found = false;
            for (arg_sitemap::iterator it = er.first; it != er.second; ++it){
                if (it->second == n){
                    found = true;
                    break;
                }
            }
            if (!found){
                fprintf(stderr, "Node missing key in site map for %ld\n", *site);
                print_node_lite(n, num_haplotypes);
                exit(1);
            }
        }
        pair<arg_clademap::iterator, arg_clademap::iterator> er = sites_clade.equal_range(n->clade);
        bool found = false;
        for (arg_clademap::iterator it = er.first; it != er.second; ++it){
            if (it->second == n){
                found = true;
                break;
            }
        }
        if (!found){
            fprintf(stderr, "Node missing key in clade map:\n");
            print_node_lite(n, num_haplotypes);
            exit(1);
        }
    }
    for (set<vert_edge>::iterator p = n->parents.begin(); p != n->parents.end(); ++p){
        bool found = false;
        for (set<vert_edge>::iterator c = p->to->children.begin(); c != p->to->children.end(); ++c){
            if (c->start == p->start && c->end == p->end && c->to == n){
                found = true;
                break;
            }
        }
        if (!found){
            fprintf(stderr, "missing parent -> node <%ld %ld>\n", p->start, p->end);
            print_node_upward_onelevel(n, 0, num_haplotypes);
            fprintf(stderr, "-----\n");
            print_node_downward_onelevel(p->to, 0, num_haplotypes);
            exit(1);
        }
    }
    for (set<vert_edge>::iterator c = n->children.begin(); c != n->children.end(); ++c){
        bool found = false;
        for (set<vert_edge>::iterator p = c->to->parents.begin(); p != c->to->parents.end(); ++p){
            if (p->start == c->start && p->end == c->end && p->to == n){
                found = true;
                break;
            }
        }
        if (!found){
            fprintf(stderr, "missing child -> node <%ld %ld>\n", c->start, c->end);
            print_node_downward_onelevel(n, 0, num_haplotypes);
            fprintf(stderr, "-----\n");
            print_node_upward_onelevel(c->to, 0, num_haplotypes);
            exit(1);
        }
        else{
            check_everything(c->to, root, sites_pos, sites_clade, num_haplotypes);
        }
    }
}


void check_node_recursive(arg_node* n, long int pos, int num_haplotypes){
    // If any layer of this node has two children that fail 4 hap test with each other
    // overlapping the same position, it's a problem.
    for (set<vert_edge>::iterator child = n->children.begin(); child != n->children.end(); ++child){
        if (child->start <= pos && child->end >= pos){
            for (set<vert_edge>::iterator child2 = n->children.begin(); child2 != n->children.end(); ++child2){
                if (child2->to != child->to && child2->start <= pos && child2->end >= pos){
                    short comp = compare_bitsets(child->to->clade, child2->to->clade);
                    if (comp == -1){
                        fprintf(stderr, "FAILED\n");
                        print_node_recursive_pos(n, pos, 0, num_haplotypes);
                        exit(1);
                    }
                }
            }
        }
    }
    for (set<vert_edge>::iterator child = n->children.begin(); child != n->children.end(); ++child){
        if (child->start <= pos && child->end >= pos){
            check_node_recursive(child->to, pos, num_haplotypes);
        }
    }
}

void check_node_parents_recursive(arg_node* n, int num_haplotypes){
    for (set<vert_edge>::iterator edge = n->parents.begin(); edge != n->parents.end(); ++edge){
        for (long int i = edge->to->start; i <= edge->to->end; ++i){
            check_node_recursive(edge->to, i, num_haplotypes);
        }
    }
}
