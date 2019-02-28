#ifndef SPARGERCTREES_H
#define SPARGERCTREES_H

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

struct failure{
    cladeset l;
    cladeset r;
    cladeset lprime;
    cladeset rprime;
    cladeset dd;
    int type;
};

void redo_failures(std::vector<failure>& failures, cladeset& bestclade);
std::vector<std::vector<recomb_event> > solve_recombs_trees(treeNode&, treeNode&);

#endif
