#!/bin/bash

#bash ~/meta_script/Meta_qc_multi.sh rawdata qc_data xxx
contig_path=$1
seq_path=$2
depth_path=$3
bin_path=$4
num_jobs=$5

# Check if the correct number of arguments is provided
if [ $# -ne 5 ]; then
echo "Usage: $0 <contig_path> <seq_path> <depth_path> <bin_path> <num_jobs>"
exit 1
fi

# Function to process a single file
process_bwa_index() {
local seq="$1"
echo $seq
BASE=${seq##*/}
SAMPLE=${BASE%%.*}
#echo $SAMPLE

if [ -e $contig_path/$BASE.sa ]; then
    echo "$SAMPLE SKIP"
else
    bwa index $seq
    echo "$SAMPLE index DONE"
fi

}

# Function to process a single file
process_bwa_map() {
local contig="$1"

SAMPLE=$(basename $contig .fa)
F=$seq_path/${SAMPLE}_1.fastq.gz
R=$seq_path/${SAMPLE}_2.fastq.gz

echo $F
#echo $SAMPLE

if [ -e $depth_path/${SAMPLE}_depth.txt ]; then
    echo "$SAMPLE SKIP"
else
    bwa mem -t 128 $contig $F $R > $depth_path/${SAMPLE}_depth.sam
    samtools sort -@ 128 $depth_path/${SAMPLE}_depth.sam > $depth_path/${SAMPLE}_depth.bam
    samtools index  -@ 128 $depth_path/${SAMPLE}_depth.bam
    jgi_summarize_bam_contig_depths --outputDepth $depth_path/${SAMPLE}_depth.txt $depth_path/${SAMPLE}_depth.bam
    echo "$SAMPLE map DONE"
    rm $depth_path/${SAMPLE}_depth.sam
    rm $depth_path/${SAMPLE}_depth.bam
    rm $depth_path/${SAMPLE}_depth.bam.bai
fi

}


# Function to process a single file
process_bin() {
local seq_bin="$1"

BASE=${seq_bin##*/}
SAMPLE=${BASE%%.*}
echo $SAMPLE

if [ -e $bin_path/${SAMPLE}_done ]; then
    echo "$SAMPLE SKIP"
else
    metabat2 -m 1500 -t 64 -i $seq_bin -a $depth_path/${SAMPLE}_depth.txt -o $bin_path/${SAMPLE}
    echo "$SAMPLE bin DONE" > $bin_path/${SAMPLE}_done
    echo "$SAMPLE bin DONE"
fi

}

# Use parallel to index
for seq in $contig_path/*.fa; do
process_bwa_index "$seq" &
((++processed_files))
[ $((processed_files % num_jobs)) -eq 0 ] && wait
done

wait


# Use parallel to map
for contig in $contig_path/*.fa; do
process_bwa_map "$contig" &
((++processed_files))
[ $((processed_files % num_jobs)) -eq 0 ] && wait
done

wait


# Use parallel to bin
for seq_bin in $contig_path/*.fa; do
process_bin "$seq_bin" &
((++processed_files))
[ $((processed_files % num_jobs)) -eq 0 ] && wait
done

wait







