#!/bin/bash

input_path=$1
output_path=$2
num_jobs=$3


# Check if the correct number of arguments is provided
if [ $# -ne 3 ]; then
echo "Usage: $0 <source_path> <output_path><num_jobs>"
exit 1
fi

# Function to process a single gzip file
process_assembly() {
local F="$1"

R=${F%_*}_2.fastq.gz
BASE=${F##*/}
SAMPLE=${BASE%_*}
echo $SAMPLE

if [ -e $output_path/${SAMPLE}/${SAMPLE}.fa ]; then
    echo "$SAMPLE SKIP"
else
    megahit -1 $F -2 $R -m 100 -t 64 -o $output_path/${SAMPLE} --out-prefix ${SAMPLE} --continue
    echo "$SAMPLE assembly DONE"
fi

}

# Use parallel to process gzip files in parallel with specified number of jobs
for F in $input_path/*_1.fastq.gz; do
process_assembly "$F" &
((++processed_files))
[ $((processed_files % num_jobs)) -eq 0 ] && wait
done
