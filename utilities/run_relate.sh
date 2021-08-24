#! /bin/bash

if [ $# -lt 3 ]; then
    >&2 echo "USAGE: run_relate.sh <basename> <N> <mu> <Relate_path>"
    >&2 echo "Converts SARGE input files to Relate format, runs Relate, and "
    >&2 echo "converts Relate output to SARGE output format. Output file name "
    >&2 echo "will be <basename>_relate.sarge.gz"
    >&2 echo "<basename> is the base file name of all SARGE input files"
    >&2 echo "<N> is the effective population size for Relate (2x normal effective population size)"
    >&2 echo "<mu> is the mutation rate"
    >&2 echo "<Relate_path> is the root folder where Relate was downloaded"
    exit 1
fi

basename=$1
N=$2
mu=$3
relate_path=$4

UTDIR=$(realpath "$0")
UTDIR=$(dirname "${UTDIR}")
SARGEDIR=$(dirname "${UTDIR}")

# Convert input files
"${UTDIR}/sarge_relate_input.py" -i "${basename}" -o "${basename}_relate"

# Run Relate -- needs to be in same directory as files
DATDIR=$(dirname $(realpath $basename))
DATNAME=$(basename $basename)
cd $DATDIR
if [ -d "${DATNAME}_relate" ]; then
    rm -r "${DATNAME}_relate"
fi



"${relate_path}/bin/Relate" \
    --mode All \
    -m $mu \
    -N $N \
    --map "${DATNAME}_relate.map" \
    --haps "${DATNAME}_relate.haps" \
    --sample "${DATNAME}_relate.samples" \
    --output "${DATNAME}_relate"

# Convert output
"${UTDIR}/sarge_relate_output.py" -b "${basename}_relate" | "${SARGEDIR}/test/treedat2sarge" -i "${basename}.geno.gz" -o "${basename}_relate.sarge.gz"


