#! /bin/bash

# Splits up input files where there are gaps between sites large enough to do so.
# Runs everything in parallel (with GNU parallel), then puts pieces back
# together.

# NOTE: requires both python2.7 and GNU parallel to be installed and in 
# current $PATH.

# Arguments:
# Directory containing input files
# Propagation distance

if [ $# -lt 2 ]; then
    >&2 echo "Usage:"
    >&2 echo "sarge_parallel.sh <directory> <propdist>"
    >&2 echo "<directory> is a directory containing all input files. Each chromosome should \
have a *.geno.gz and *.sites file. These will be split into multiple files where \
possible; output will be re-joined at the end"
    exit 1
fi

# Detect root directory of SARGE installation
sarge_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null && pwd )"
sarge_dir=$(dirname $sarge_dir)

dir=$1
propdist=$2

for fn in $(find $dir -mindepth 1 -maxdepth 1 -type f -name "*.geno.gz" -not -name "*split*" | sort -V); do
    fn="${fn##*/}"
    fnbase="${fn%.geno.gz}"
    if [ ! -e "${dir}/${fnbase}.sites" ]; then
        >&2 echo "WARNING: file ${dir}/${fnbase}.sites not found"
    elif [ -e "${dir}/${fnbase}_split_${propdist}.1.geno.gz" ]; then
        >&2 echo "Files already split for ${fnbase}"
    else
        >&2 echo "Splitting input file ${fn}..."
        "${sarge_dir}/utilities/split_input_gaps.py" -g "${dir}/${fn}" -s "${dir}/${fnbase}.sites" -p $propdist -o "${dir}/${fnbase}_split_${propdist}"
    fi
done

if [ -e "${dir}/sarge_parallel.jobs" ]; then
    rm "${dir}/sarge_parallel.jobs"
fi

# Create jobs file
>&2 echo "Creating jobs file..."
for fn in $(find $dir -mindepth 1 -maxdepth 1 -type f -name "*_split_${propdist}*.geno.gz" | sort -V); do
    fn="${fn##*/}"
    fnbase="${fn%.geno.gz}"
    noheader=""
    index="${fnbase%.geno.gz}"
    index="${index##*_split*.}"
    if [ "$index" -eq 1 ]; then
        noheader=""
    else
        noheader="-H "
    fi
    exclbed=""
    if [ -e "${dir}/${fnbase%_split*}_exclude.bed" ]; then
        exclbed="-e ${dir}/${fnbase%_split*}_exclude.bed "
    fi
    echo "nice ${sarge_dir}/bin/sarge ${noheader}${exclbed}-i ${dir}/${fnbase}.geno.gz -s ${dir}/${fnbase}.sites -p ${propdist} -o ${dir}/${fnbase}.sarge.gz" >> "${dir}/sarge_parallel.jobs"
done

# Run jobs using GNU parallel
time parallel < "${dir}/sarge_parallel.jobs"

if [ -e "${dir}/sarge_parallel.jobs" ]; then
    rm "${dir}/sarge_parallel.jobs"
fi

# Put results back together and index results

if [ -e "${dir}/sarge_index.jobs" ]; then
    rm "${dir}/sarge_index.jobs"
fi

for fnbase in $(find $dir -type f -name "*_split_${propdist}*.sarge.gz" -printf '%f\n' | sed "s/_split_${propdist}./\&/" | \
    cut -d'&' -f1 | sort -V | uniq); do
    >&2 echo "Joining files for ${fnbase}"
    filescat=$(find "${dir}" -mindepth 1 -maxdepth 1 -type f -name "${fnbase}_split_${propdist}*.sarge.gz" | sort -V)
    >&2 echo $filescat
    zcat $filescat | gzip -c - > "${dir}/${fnbase}_joined_${propdist}.sarge.gz"
    echo "${sarge_dir}/bin/index ${dir}/${fnbase}_joined_${propdist}.sarge.gz" >> "${dir}/sarge_index.jobs"
done

parallel < "${dir}/sarge_index.jobs"


