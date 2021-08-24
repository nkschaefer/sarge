# sarge

Fast, heuristic parsimony-based ancestral recombination graph inference

## Warning
This repository is under heavy development. Output file format might change slightly in the future, and the organization of the code is subject to change as well. The code that is on here right now might be somewhat less than beautiful.

## Requirements:

### Essential stuff
* A C++ compiler capable of supporting the C++11 standard; for example:
  * [g++](https://gcc.gnu.org/) version >= 4.9.0
  * [clang++](https://clang.llvm.org/cxx_status.html/) version >= 3.3
* [GNU Make](https://www.gnu.org/software/make/)
* [zlib](https://zlib.net/)

### Optional stuff
* For creating BED files of CpG sites (to exclude from analysis), a C compiler:
  * [gcc](https://gcc.gnu.org/)
* For efficiently creating input files:
  * [HTSLib](http://www.htslib.org/download/)
  * [BCFTools](http://www.htslib.org/download/)
* For splitting input files for parallelization (utilities/split_input_gaps.py), [Python 2.7](https://www.python.org/download/releases/2.7/) installed and in your ``$PATH``
* For running jobs automatically in parallel (utilities/sarge_parallel.sh), [GNU Parallel](https://www.gnu.org/software/parallel/) installed and in your ``$PATH``
* For visualizing trees (utilities/view_trees.py), the Phylo package from [BioPython](https://biopython.org/) must be available

## Compilation instructions
After downloading, programs can be compiled by running `make` in the root directory. This will create executable programs in the `bin` directory.

### Number of haplotypes

**Important**: due to sarge's use of bitsets, and the fact that bitset sizes must be defined at compile time, you *must* define the maximum number of haplotypes you plan to use at the time you compile sarge. This means that if your dataset contains 100 diploid genomes, you must tell sarge at the time of compilation that you plan to use 200 haplotypes. You can set this number higher than the actual size of your dataset, but this can cause the program to run slower than necessary. To set the number of haplotypes:

``make MAXHAPS=<maximum number of haplotypes>``

Where ``<maximum number of haplotypes>`` should be replaced with an integer number. If you compile for a certain number of haplotypes and later want to change it, simply go back to the root directory and run

``make clean``

``make MAXHAPS=<new number>``

There is another place in the program where bitsets are used to store choices of clades to include (or exclude) from a recombination event. By default, it allows a number of clades up to 5 times the maximum number of haplotypes. If you experience a problem with this, you should receive an error message telling you what went wrong and how to recompile. The upshot is that you will need to ``make clean`` and then run ``make`` again, increasing either ``NODESETSIZE`` to a number greater than 5 times whatever you set ``MAXHAPS`` to before, or to increase ``MAXHAPS`` further. 

### Input data conversion

[HTSLib](http://www.htslib.org/download/) is a useful set of C-language libraries used for processing standard high-throughput sequence data types. It is used by [samtools](http://www.htslib.org/doc/samtools.html) and [bcftools](https://samtools.github.io/bcftools/), among other programs. By default, SARGE does not depend on HTSLib. The program ``bin/vcf2sarge``, however, can use HTSLib to make the process of converting input data files (VCF or BCF format) into SARGE format much faster. If you have genotype data in VCF format, the program ``bin/vcf2sarge`` can be used to process the VCF file directly, but this will take a long time. To make processing much faster, install HTSLib and compile with HTSLib support:

``make clean``

``make HTSLIB=1``

And this will allow ``bin/vcf2sarge`` to read BCF, rather than VCF, input. If you still want to work with VCF files, simply convert them to BCF using ``bcftools view`` and pipe to ``bin/vcf2sarge`` (yes, this is still faster than parsing VCF directly):

``bcftools view -O u <file.vcf> | bin/vcf2sarge ...``

Where ``<file.vcf>`` is the VCF file and ``...`` represents other options to ``bin/vcf2sarge.``

### Random choices

There are times when SARGE has no principled way of selecting between two equivalent choices, except to randomly choose one. This should not be a problem, except in the case of debugging, when replicable output is desired (this can also be achieved by setting the random number seed in ``src/sarge_main.cpp`` to a predetermined number rather than the current time). If for some reason random choices need to be disabled, recompile with RAND=0:

``make clean``

``make RAND=0``

### Parsing large strings

This will not be relevant for most users, but SARGE has a function for parsing [ms](http://home.uchicago.edu/rhudson1/source.html) output for the sake of testing. This works fine, unless many haplotypes (thousands) are included in the simulation, leading trees in the simulated ARG to consist of very long [Newick-format](https://en.wikipedia.org/wiki/Newick_format) strings. The default C++ library g++ uses (libstdc++) has a bug that prevents its regular expression library from being able to process very large strings without crashing. For this reason, a workaround in this very specific case is to install the [clang](https://clang.llvm.org/cxx_status.html) compiler, along with [libc++](https://libcxx.llvm.org/) and tell ``make`` to use these when compiling. For the sake of simplicity, I've added a variable for that:

``make clean``

``make CLANG=1``

### Combining multiple options

All of these above options can be combined during compilation. In other words, if you want to support up to 500 haplotypes and build-in support for HTSLib (supposing it's installed and available in your ``$LIBRARY_PATH`` environment variable:

``make MAXHAPS=500 HTSLIB=1``

If you make a mistake, just ``make clean`` and start over.

## Preparing data

SARGE expects 3 main input files, in a very similar format to other ARG programs. The main one, a \*.geno.gz file, is a gzip-compressed text file where each line represents a variable site and each position in the line is 1 or 0, representing if an individual (phased) haplotype has an ancestral (0) or derived (1) allele. All lines must be the same length, and there can be no missing data. The outgroup you use should not be included as one of the sequences -- instead, distance from the outgroup is measured in the number of fixed derived mutations in your input data (i.e. if you are studying humans with the chimpanzee as outgroup, then every human/chimpanzee difference should be represented as a row of all "1" entries in this file).

Next, you should include a \*.sites file. Lines in this file correspond to variable sites in the \*.geno.gz file. Each entry in the sites file is a tab-separated chromosome and position of the variable site.

Finally, you should have a \*.haps file. Each line of this file should correspond to a character (column) in the \*.geno.gz file. This file lists which phased haplotype each column in the \*.geno.gz refers to. In other words, if you have 10 diploid individuals (20 phased haplotypes), each line in \*.geno.gz should have 20 characters, and the \*.haps file should have 20 lines, like "individual1-1," "individual1-2," "individual2-1," "individual2-2," and so on.

Some programs within SARGE can use a \*.alleles file as well. This file corresponds to the \*.sites file; each line in the \*.alleles file represents a variable site. There should be two tab-separated characters per line, the ancestral allele (the one matching "0" in the input file) and the derived allele (matching "1" in the input file).

There is a program (bin/vcf2sarge) that can help with this conversion. It requires that you have [HTSLib](http://www.htslib.org/download/) installed, compile SARGE with the HTSLIB=1 option, and pipe output from [BCFTools](http://www.htslib.org/download/) to create input data.

One parameter that bin/vcf2sarge can take is a list of (minimum and maximum) coverage depth cutoffs per sample. This is extremely important to set when working with ancient data; a histogram of coverage depths can be produced from BCF/VCF using bin/vcf_depths. These can then be inspected to determine the correct coverage cutoffs before running bin/vcf2sarge.

## Running SARGE

The main program in SARGE is bin/sarge. It requires a \*.geno.gz file, \*.sites file, and a propagation distance. The propagation distance parameter describes how many bases (upstream and downstream) a given site is allowed to communicate its existence. Increasing this number will cause slower execution and more memory usage but will make somewhat better trees. Typically, there is a point of diminishing returns beyond which it is not useful to increase this parameter (it will only slow things down). SARGE creates an output file (gzipped binary format), along with a "recombination" file, which is a tab separated file of chromosome, start position, end position, upstream clade, downstream clade, and moving clade. 

SARGE output can be indexed (bin/index) and retrieved using region syntax similar to the UCSC Genome Browser and SAMTools. This output can then be piped to the other programs (or just zcat the file to run on everything) for downstream analyses. If you want to export Newick-format trees that can be used in another program, run the program bin/trees2newick, which will output tab separated chromosome, position, and Newick-format tree. For example, to dump chromosome 12, site 100,000 to 200,000 in Newick format: 

``bin/index [OUTPUTFILE] 12:100000-200000 | bin/trees2newick -v [HAPSFILE]``

Currently, SARGE is designed to expect only one chromosome per set of input files, but this may change in the future.

All programs should display a usage message if run with no arguments.

## Utility Programs

Some utility programs are provided for help with specific tasks. 

*<code>utilities/sarge_parallel.sh</code> will split up a set of input files into chunks that can be run in parallel (dividing where there are large gaps in a sequence, such as centromeres). It will then run as much as possible in parallel and join output files together when finished.

*<code>utilities/view_trees.py</code> can be used to look at individual trees (must select all trees (<code>zcat [output file]</code>) or a specific region (<code>bin/sarge_index [output file] chr3:5000000-10000000</code>) and pipe to <code>bin/trees2newick</code> and then to <code>utilities/view_trees.py</code>

*<code>utilities/run_relate.sh</code> takes SARGE-format input files, converts to Relate format, runs Relate, and converts output to SARGE-format output files (requires Relate binaries, last tested with Relate v1.0.17)

*<code>utilities/run_tsinfer.sh</code> takes SARGE-format input files, runs tsinfer and tsdate, and converts output to SARGE-format output files (requires tsinfer and tsdate to be installed in python in the current $PATH)

*<code>utilities/run_rentplus.sh</code> takes SARGE-format intput files, converts to Rent+ format, runs Rent+, and converts output to SARGE-format output files (requires RentPlus.jar file)

## Testing programs

There are several programs for testing output by comparing to trees from an MS simulation (with the -T parameter enabled for outputting the simulation ARG). Each program takes uncompressed SARGE output format as stdin: i.e. <code> zcat testdata.sarge.gz | [program]</code>

These programs are:

*<code>test/score_trees_ms</code> to gauge accuracy of trees

*<code>test/ms_branch_compare</code> to compare branch lengths on correct clades to those from the simulation

*<code>test/ms_corr_dist</code> to determine, when a clade is incorrect, the distance to a position at which the clade would have been correctly inferred

*<code>test/ms_missing_clades</code> to check the sizes (in number of leaf haplotypes) of true clades (from a simulation) that were not captured in inferred/output data
 

