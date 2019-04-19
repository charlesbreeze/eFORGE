#!/bin/bash

key=$1
bedFn=$2

dsArgs=()

mkdir -p ${key}

module add bedops/2.4.30-typical
ulimit -n 2048

#
# build lexicographically-ordered map list
#
cp ${key}-files_fns.txt ${key}/fns.txt

#
# map against files in map list
#
while IFS='' read -r fn || [[ -n "$fn" ]]; do
    fn=${fn}
    sampleKey=`echo "$fn" | awk 'match($0, /\/[A-Za-z0-9_\-]+./) { print substr( $0, RSTART+1, RLENGTH-2 );}'`
    if [[ -z "${sampleKey// }" ]]; then
        exit -1
    fi
    echo ${fn}
    echo ${sampleKey}
    sampleArgs+=("${key}/$sampleKey.indicator")
    tempBedFile=$(mktemp /tmp/bitstrings.XXXXXX)
    convertCmd="convert2bed --input=wig < ${fn} > ${tempBedFile}"
    eval ${convertCmd}
    convertReturnCode=$?
    #
    # treat input as BED, but strip headers
    #
    if [ $convertReturnCode != 0 ]; then
        echo "treating input as BED and stripping headers"
        sort-bed ${fn} > ${tempBedFile}
        head ${tempBedFile}
    fi
    bedmap --indicator ${bedFn} ${tempBedFile} > ${key}/${sampleKey}.indicator
    echo "${key}/${sampleKey}.indicator"
    rm ${tempBedFile}
done < ${key}/fns.txt

#
# build matrix
#
paste <(cut -f4 ${bedFn}) "${sampleArgs[@]}" \
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
