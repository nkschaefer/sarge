#ifndef SARGE_DEBUG_H
#define SARGE_DEBUG_H

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


/**
 * Debugging related functions
 */
void check_for_holes(arg_node* node, arg_node* root, int num_haplotypes);
void check_ov_parents(arg_node* node, int num_haplotypes);
void check_ov_children(arg_node* node, int num_haplotypes);
void check_solved_connections(arg_sitemap& sites_pos, int num_haplotypes);
void check_cp_edges(arg_sitemap& sites_pos, int num_haplotypes);
void check_everything(arg_node* n, arg_node* root, arg_sitemap& sites_pos, arg_clademap& sites_clade, int num_haplotypes);
void check_node_recursive(arg_node* n, long int pos, int num_haplotypes);
void check_node_parents_recursive(arg_node* n, int num_haplotypes);
void check_solved_connected(arg_sitemap&, long int, long int, int num_haplotypes);

#endif
