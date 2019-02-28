/**
 * Contains code used to solve recombination events between two fully known trees.
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
#include "sets.h"
#include "treeNode.h"
#include "serialize.h"
#include "common.h"
#include "recomb.h"
#include "rc_trees.h"

using namespace std;
using std::cout;
using std::endl;

void get_failures(vector<cladeset >& t1flat,
    vector<cladeset >& t2flat,
    vector<failure>& failures){
    
    for (vector<cladeset >::iterator clade1 = t1flat.begin();
        clade1 != t1flat.end(); ++clade1){
        for (vector<cladeset >::iterator clade2 = t2flat.begin(); 
            clade2 != t2flat.end(); ++clade2){
            short comp = compare_bitsets(*clade1, *clade2);
            if (comp == -1){
                failure f;
                f.l = *clade1;
                f.r = *clade2;
                f.lprime = set_diff_bitset(*clade1, *clade2);
                f.rprime = set_diff_bitset(*clade2, *clade1);
                f.dd = set_int_bitset(*clade1, *clade2);
                failures.push_back(f);
            }
        }
    }   
}

recomb_event bestclade2recomb(vector<failure>& failures,
    cladeset& bestclade){
    
    cladeset alpha;
    cladeset beta;
    cladeset un;
    bool alpha_set = false;
    bool beta_set = false;
    bool un_set = false;
    
    for (vector<failure>::iterator f = failures.begin(); f != failures.end(); ++f){
        if (f->lprime == bestclade){
            if (!alpha_set || issubset_bitset(f->l, alpha)){
                alpha_set = true;
                alpha = f->l;
            }
            if (!un_set){
                un = f->l | f->r;
            }
            else{
                un |= (f->l | f->r);
            }
            
        }
        else if (f->rprime == bestclade){
            if (!beta_set || issubset_bitset(f->r, beta)){
                beta_set = true;
                beta = f->r;
            }
            if (!un_set){
                un = f->l | f->r;
            }
            else{
                un |= (f->l | f->r);
            }
        }
        else if (f->dd == bestclade){
            if (!alpha_set || issubset_bitset(f->l, alpha)){
                alpha_set = true;
                alpha = f->l;
            }
            if (!beta_set || issubset_bitset(f->r, beta)){
                beta_set = true;
                beta = f->r;
            }
            if (!un_set){
                un = f->l | f->r;
            }
            else{
                un |= (f->l | f->r);
            }
        }
    }
    recomb_event r;
    r.leaving = bestclade;
    r.leaving_known = true;
    if (alpha_set){
        r.alpha = alpha;
        r.alpha_known = true;
    }
    else{

        r.down = true;
        r.alpha = un;
        r.alpha_known = true;
    }
    if (beta_set){
        r.beta = beta;
        r.beta_known = true;
    }
    else{

        r.up = true;
        r.beta = un;
        r.beta_known = true;
    }
    
    return r;
}

void redo_failures(vector<failure>& failures,
    cladeset& bestclade){
    // Remove anything explained by this recombination event; eliminate
    // the leaving group from all those that remain.
    for (vector<failure>::iterator f = failures.begin(); f != failures.end();){
        if (f->lprime == bestclade || f->rprime == bestclade || f->dd == bestclade){
            failures.erase(f);
        }
        else{
            f->lprime = set_diff_bitset(f->lprime, bestclade);
            f->rprime = set_diff_bitset(f->rprime, bestclade);
            f->dd = set_diff_bitset(f->dd, bestclade);
            ++f;
        }
    }
}

/** 
 * Given two trees, finds all recombination events between them.
 */
std::vector<std::vector<recomb_event> > solve_recombs_trees(treeNode& left, treeNode& right){
    // Get all failures
    vector<cladeset > t1flat;
    vector<cladeset > t2flat;
    left.flatten(t1flat);
    right.flatten(t2flat);
    vector<failure> failures;
    get_failures(t1flat, t2flat, failures);
    
    vector<vector<recomb_event> > results;
    
    if (failures.size() == 0){
        return results;
    }
    
    while(failures.size() > 0){
        vector<string> d;
        
        map<cladeset, int, cladesetcomp> lv_counts;
        map<cladeset, int, cladesetcomp> lv_types;
        for (vector<failure>::iterator f = failures.begin(); f != failures.end();
            ++f){
            if (f->lprime.count() > 0){
                if (lv_counts.count(f->lprime) == 0){
                    lv_counts.insert(make_pair(f->lprime, 0));
                }
                lv_counts[f->lprime]++;
                if (lv_types.count(f->lprime) == 0){
                    lv_types.insert(make_pair(f->lprime, 1));
                }
            }
            if (f->dd.count() > 0){
                if (lv_counts.count(f->dd) == 0){
                    lv_counts.insert(make_pair(f->dd, 0));
                }
                lv_counts[f->dd]++;
                if (lv_types.count(f->dd) == 0){
                    lv_types.insert(make_pair(f->dd, 2));
                }
                else{
                    lv_types[f->dd] = 2;
                }
            }
            if (f->rprime.count() > 0){
                if (lv_counts.count(f->rprime) == 0){
                    lv_counts.insert(make_pair(f->rprime, 0));
                }
                lv_counts[f->rprime]++;
                if (lv_types.count(f->rprime) == 0){
                    lv_types.insert(make_pair(f->dd, 3));
                }
            }
        }
        // Determine & choose best clade.
        // If there's a tie, choose the one that is lateral (not up or down).
        int maxcount = -1;
        vector<cladeset > maxclades;
        for (map<cladeset, int>::iterator lc = 
            lv_counts.begin(); lc != lv_counts.end(); ++lc){
            if (lc->second > maxcount){
                maxcount = lc->second;
                maxclades.clear();
                maxclades.push_back(lc->first);
            }
            else if (lc->second == maxcount){
                maxclades.push_back(lc->first);
            }
        }
        
        if (maxclades.size() > 0){
            vector<recomb_event> s;
            for (vector<cladeset >::iterator cl = 
                maxclades.begin(); cl != maxclades.end(); ++cl){
                s.push_back(bestclade2recomb(failures, *cl));
            }
            redo_failures(failures, maxclades[0]);
            results.push_back(s);
        }
    }
    return results;
}

