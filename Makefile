# http://gcc.gnu.org/onlinedocs/cpp/Environment-Variables.html
# These variables must be set, if stuff is in a weird location:
# CPLUS_INCLUDE_PATH = anything that should be included with -I
# LIBRARY_PATH = anything that should be included with -L

OBJDEPS = serialize.o sets.o treeNode.o common.o argnode.o recomb.o debug.o rc_trees.o
CPPDEPS = src/serialize.cpp src/sets.cpp src/treeNode.cpp src/common.cpp src/argnode.cpp src/recomb.cpp src/debug.cpp src/rc_trees.cpp
HDEPS = src/serialize.h src/sets.h src/treeNode.h src/common.h src/argnode.h src/recomb.h src/debug.h src/rc_trees.h
DEBUG_MODE ?= 0
MAXHAPS ?= 5000
RAND ?= 1
SHELL=/bin/bash
OPTS = -D DEBUG_MODE=$(DEBUG_MODE) -D MAXHAPS=$(MAXHAPS) -D RAND=$(RAND) -D NODESETSIZE=$$(( $(MAXHAPS) * 5 ))
ifeq ($(HTSLIB), 1)
HTSOPTS = -D HTSLIB -lhts
endif
COMP = g++
ifeq ($(CLANG), 1)
COMP = clang++
OPTS += -stdlib=libc++
endif
CCOMP = gcc

all: bin/sarge bin/trees2newick bin/index bin/index_geno bin/ms_tree_site bin/admix_scan bin/tmrcas_scan bin/dstat bin/upgma bin/vcf2sarge bin/vcf_depths bin/cpg_bed bin/ms2sarge bin/haplens bin/dump_haplens bin/bait_scan bin/rth bin/rel_ils

bin/sarge: src/main.cpp $(OBJDEPS)
	$(COMP) -std=c++11 $(OPTS) -g src/main.cpp -o bin/sarge $(OBJDEPS) -lz 

bin/trees2newick: src/trees2newick.cpp $(OBJDEPS)
	$(COMP) -std=c++11 $(OPTS) src/trees2newick.cpp -o bin/trees2newick $(OBJDEPS) -lz

bin/index: src/index.cpp $(OBJDEPS) src/serialize.h
	$(COMP) -std=c++11 $(OPTS) src/index.cpp -o bin/index $(OBJDEPS) -lz

bin/index_geno: src/index_geno.cpp $(OBJDEPS)
	$(COMP) -std=c++11 $(OPTS) src/index_geno.cpp -o bin/index_geno $(OBJDEPS) -lz

bin/ms_tree_site: src/ms_tree_site.cpp src/common.cpp $(OBJDEPS) src/common.h
	$(COMP) -std=c++11 $(OPTS) src/ms_tree_site.cpp -o bin/ms_tree_site $(OBJDEPS) -lz

bin/admix_scan: src/admix_scan.cpp $(OBJDEPS)
	$(COMP) -std=c++11 $(OPTS) src/admix_scan.cpp -o bin/admix_scan $(OBJDEPS) -lz

bin/bait_scan: src/bait_scan.cpp $(OBJDEPS)
	$(COMP) -std=c++11 $(OPTS) src/bait_scan.cpp -o bin/bait_scan $(OBJDEPS) -lz

bin/rth: src/rth.cpp $(OBJDEPS)
	$(COMP) -std=c++11 $(OPTS) src/rth.cpp -o bin/rth $(OBJDEPS) -lz

bin/rel_ils: src/rel_ils.cpp $(OBJDEPS)
	$(COMP) -std=c++11 $(OPTS) src/rel_ils.cpp -o bin/rel_ils $(OBJDEPS) -lz

bin/tmrcas_scan: src/tmrcas_scan.cpp $(OBJDEPS)
	$(COMP) -std=c++11 $(OPTS) src/tmrcas_scan.cpp -o bin/tmrcas_scan $(OBJDEPS) -lz
	
bin/dstat: src/dstat.cpp $(OBJDEPS)
	$(COMP) -std=c++11 $(OPTS) src/dstat.cpp -o bin/dstat $(OBJDEPS) -lz

bin/upgma: src/upgma.cpp $(OBJDEPS)
	$(COMP) -std=c++11 $(OPTS) -pthread src/upgma.cpp -o bin/upgma $(OBJDEPS) -lz

bin/vcf2sarge: src/vcf2sarge.cpp
	$(COMP) -std=c++11 src/vcf2sarge.cpp -o bin/vcf2sarge $(HTSOPTS) -lz

bin/vcf_depths: src/vcf_depths.cpp
	$(COMP) -std=c++11 src/vcf_depths.cpp -o bin/vcf_depths $(HTSOPTS) -lz

bin/cpg_bed: src/cpg_bed.c
	$(CCOMP) src/cpg_bed.c -o bin/cpg_bed $(HTSOPTS) -lz

bin/ms2sarge: src/ms2sarge.cpp src/common.cpp src/treeNode.cpp $(OBJDEPS)
	$(COMP) -std=c++11 $(OPTS) src/ms2sarge.cpp -o bin/ms2sarge $(OBJDEPS) -lz

bin/haplens: src/haplens.cpp src/common.cpp src/treeNode.cpp src/serialize.cpp $(OBJDEPS)
	$(COMP) -std=c++11 $(OPTS) src/haplens.cpp -o bin/haplens $(OBJDEPS) -lz

bin/dump_haplens: src/dump_haplens.cpp src/common.cpp src/treeNode.cpp  $(OBJDEPS)
	$(COMP) -std=c++11 $(OPTS) src/dump_haplens.cpp -o bin/dump_haplens $(OBJDEPS) -lz

serialize.o: src/serialize.cpp src/serialize.h
	$(COMP) -std=c++11 $(OPTS) -g -c src/serialize.cpp

sets.o: src/sets.cpp
	$(COMP) -std=c++11 $(OPTS) -g -c src/sets.cpp

treeNode.o: src/treeNode.cpp src/serialize.h src/sets.h
	$(COMP) -std=c++11 $(OPTS) -c src/treeNode.cpp

common.o: src/common.cpp $(HDEPS)
	$(COMP) -std=c++11 $(OPTS) -g -c src/common.cpp

argnode.o: src/argnode.cpp $(HDEPS) common.o sets.o serialize.o debug.o
	$(COMP) -std=c++11 -Wall $(OPTS) -g -c src/argnode.cpp

recomb.o: src/recomb.cpp $(HDEPS) common.o argnode.o
	$(COMP) -std=c++11 -Wall $(OPTS) -g -c src/recomb.cpp

debug.o: src/debug.cpp $(HDEPS)
	$(COMP) -std=c++11 $(OPTS) -g -c src/debug.cpp

rc_trees.o: src/rc_trees.cpp src/rc_trees.h
	$(COMP) -std=c++11 $(OPTS) -g -c src/rc_trees.cpp

clean:
	rm *.o
