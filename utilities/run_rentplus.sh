#! /bin/bash

if [ $# -lt 2 ]; then
    >&2 echo "USAGE: run_rentplus.sh <basename> <Rentplus_path>"
    >&2 echo "Given a <basename> for SARGE input files, runs Rent+ and converts"
    >&2 echo "the output to SARGE format."
    >&2 echo "Output file will be <basename>_rentplus.sarge.gz"
    >&2 echo "<Rentplus_path> is the directory containing the RentPlus.jar file."
    >&2 echo "A file called <basename>_rentplus.sarge.gz will be created containing results."
    exit 1
fi

basename=$1
rentplus_path=$2

UTDIR=$(realpath "$0")
UTDIR=$(dirname "${UTDIR}")
SARGEDIR=$(dirname "${UTDIR}")

"${UTDIR}/sarge_rentplus_input.py" "${basename}" "${basename}_rentplus"
java -jar "${rentplus_path}/RentPlus.jar" -t "${basename}_rentplus" 
cat "${basename}_rentplus.trees" | sed 's/^/seq\t/' | "${SARGEDIR}/test/newick2sarge" -i "${basename}.geno.gz" -d -o "${basename}_rentplus.sarge.gz"
