#!/bin/bash

########################
# MAIN ANALYSIS SCRIPT #
########################
# args:
# 1. full path to working directory
# 2. full path to data directory
# 3. sample names
## these will probably never change, so i'm just going to copy paste them here.
# WORKING DIRECTORY: /home/kchan_mb18/p1
# DATA DIRECTORY: /home/kchan_mb18/project_1 <-- note this is a symlink!
# SAMPLE NAMES: /home/kchan_mb18/p1/sample_names.txt
# ./master.sh /home/kchan_mb18/p1 /home/kchan_mb18/project_1 /home/kchan_mb18/p1/sample_names.txt

WORK_DIR=$1
DATA_DIR=$2
FASTQC=$(which fastqc)
BWA=$(which bwa)
THREADS=5

usage(){
echo "
MICB405: Project 1
Written by Kevin Chan

Description:  main script to run all analyses for project 1.

Usage:        master.sh working_dir data_dir sample_names

positional arguments:
 working_dir       full path to the working directory.
 data_dir          full path to the data directory.
 sample_names	   text file with each line containing a single sample name.
"   
}

if [ -z "$1" ] || [[ $1 == -h ]] || [[ $1 == --help ]]; then
	usage
	exit
fi

echo
echo "#################################"
echo "# STARTING MAIN ANALYSIS SCRIPT #"
echo "#################################"
echo 

if [[ -n "$WORK_DIR" && -n "$DATA_DIR" ]]; then
	echo "working directory is $WORK_DIR"
	echo "data directory is $DATA_DIR"
else
	echo "ERROR: please pass in the full path to your working and data directories!"
	exit 1
fi

cd $WORK_DIR

##########
# FASTQC #
##########
echo
echo "################################"
echo "# (POTENTIALLY) RUNNING FASTQC #"
echo "################################"
echo 

mkdir -p $WORK_DIR/fastqc

if [ -z "$(ls $WORK_DIR/fastqc)" ]; then
	echo "running fastqc for all samples..."
	$FASTQC -t $THREADS -o $WORK_DIR/fastqc
else
	echo "detected fastqc output, skipping..."
fi

#############
# ALIGNMENT #
#############
echo
echo "######################"
echo "# ALIGNMENT WITH BWA #"
echo "######################"
echo 

mkdir -p $WORK_DIR/align/index $WORK_DIR/align/sam $WORK_DIR/align/logs

# should only need to do this once
# cp $DATA_DIR/ref_genome.fasta $WORK_DIR/align/index
# echo "generating index from reference fasta..."
# $BWA index -p $WORK_DIR/align/index/ref_genome $WORK_DIR/align/index/ref_genome.fasta

while read name; do
	echo "aligning for sample $name"
	# might not make sense to align non-human to human ref genome...
	$BWA mem -t $THREADS $WORK_DIR/align/index/ref_genome $DATA_DIR/"$name"_1.fastq.gz $DATA_DIR/"$name"_1.fastq.gz \
		1> $WORK_DIR/align/sam/$name.sam 2> $WORK_DIR/align/logs/$name_log.txt
done < $WORK_DIR/sample_names.txt

echo
echo "########"
echo "# DONE #"
echo "########"
echo
