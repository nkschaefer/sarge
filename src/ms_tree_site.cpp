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
#include <fstream>
#include <zlib.h>
#include "sets.h"
#include "treeNode.h"
#include "serialize.h"
#include "common.h"

using std::cout;
using std::endl;
using namespace std;

int main(int argc, char *argv[]) {    
    int bufsize = 1048576;
    
    if (argc < 3){
        fprintf(stderr, "USAGE: ms_tree_site msfile site\n");
        fprintf(stderr, "Where msfile is an output file from MS and site is an integer index\n");
        exit(1);
    }
    
    string msfile = string(argv[1]);
    
    long int pos = atol(argv[2]);

    // Parse trees from MS run
    map<pair<long int, long int>, treeNode> ms_trees;
    fprintf(stderr, "Parsing trees from MS run...\n");
    int nhaps = parse_ms_file(ms_trees, msfile);
    fprintf(stderr, "Finished\n");
    
    for (map<pair<long int, long int>, treeNode>::iterator ms_it = ms_trees.begin();
        ms_it != ms_trees.end(); ++ms_it){
        
        if (ms_it->first.first > pos){
            break;
        }
        if (ms_it->first.first <= pos && ms_it->first.second >= pos){
            
            vector<string> dummy;
            fprintf(stdout, "%s\n", ms_it->second.newick(false, false, dummy).c_str());
        }
    }
    
}
