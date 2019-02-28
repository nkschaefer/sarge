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
* For splitting input files for parallelization, [Python 2.7](https://www.python.org/download/releases/2.7/) installed and in your ``$PATH``
* For running jobs automatically in parallel, [GNU Parallel](https://www.gnu.org/software/parallel/) installed and in your ``$PATH``

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

 

