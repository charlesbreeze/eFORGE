#!/bin/bash

key=$1
bedFn=$2
parentDir=$3
pattern=$4

dsArgs=()

mkdir -p ${key}

module add bedops/2.4.29-typical

#
# build lexicographically-ordered map list
#
find ${parentDir} -name ${pattern} -print0 | sort -z | xargs -r0 echo | sed 's/ /\n/g' > ${key}/fns.txt

#
# map against files in map list
#
while IFS='' read -r fn || [[ -n "$fn" ]]; do
    dsKey=`echo "$fn" | awk 'match($0, /DS[0-9]+/) { print substr( $0, RSTART, RLENGTH );}'`
    if [[ -z "${dsKey// }" ]]; then
        exit -1
    fi
    echo ${fn}
    echo ${dsKey}
    dsArgs+=("${key}/$dsKey.indicator")
    convert2bed --input=wig < ${fn} \
        | bedmap --indicator ${bedFn} - \
        > ${key}/${dsKey}.indicator
done < ${key}/fns.txt

#
# build matrix
#
paste <(cut -f4 ${bedFn}) "${dsArgs[@]}" \
    | gzip -c - \
    > ${key}/bits.mtx.gz

#
# collapse bits to bitstring and apply summation
#
gunzip -c ${key}/bits.mtx.gz \
    | awk -v OFS="\t" '{ s=0; b=""; for(i=2;i<=NF;i++){ s+=$i; b=b$i; } print $1,s,b; }' - \
    | gzip -c - \
    > ${key}/bitstrings.txt.gz

#
# cleanup
#
rm ${key}/*.indicator
