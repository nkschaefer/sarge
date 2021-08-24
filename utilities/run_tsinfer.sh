#! /bin/bash

if [ $# -lt 3 ]; then
    >&2 echo "USAGE: run_tsinfer.sh <basename> <N> <mu>"
    >&2 echo "Given a <basename> for SARGE input files, runs tsinfer and converts"
    >&2 echo "the output to SARGE format. Output file will be <basename>_tsinfer.sarge.gz"
    >&2 echo "N is the effective population size (for branch length inference by tsdate module)"
    >&2 echo "mu is the mutation rate (for branch length inference by tsdate module)"
    exit 1
fi

basename=$1
N=$2
mu=$3

UTDIR=$(realpath "$0")
UTDIR=$(dirname "${UTDIR}")
SARGEDIR=$(dirname "${UTDIR}")

"${UTDIR}/sarge_tsinfer.py" --N $N --m $mu --basename $basename | "${SARGEDIR}/test/treedat2sarge" -i "${basename}.geno.gz" -o "${basename}_tsinfer.sarge.gz" 
