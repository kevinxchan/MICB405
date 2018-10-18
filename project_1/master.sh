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
SAMPLE_NAMES=$3

FASTQC=$(which fastqc)
TRIMMOMATIC=$(which trimmomatic)
BWA=$(which bwa)
SAMTOOLS=$(which samtools)

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

if [[ -n "$WORK_DIR" && -n "$DATA_DIR" && -n "$SAMPLE_NAMES" ]]; then
	echo "working directory is $WORK_DIR"
	echo "data directory is $DATA_DIR"
	echo "sample names file: $SAMPLE_NAMES"
else
	echo "ERROR: please pass in the full path to your working and data directories, and sample_names.txt!"
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

################
# QUALITY TRIM #
################
# echo
# echo "#################################"
# echo "# QUALITY TRIM WITH TRIMMOMATIC #"
# echo "#################################"
# echo 

# mkdir -p $WORK_DIR/qual_trimmed/paired $WORK_DIR/qual_trimmed/unpaired
# # some obscure location. Nextera adapters (found by looking at fastQC output)
# adapter_path=/home/linuxbrew/.linuxbrew/Cellar/trinity/2.8.3/libexec/trinity-plugins/Trimmomatic-0.36/adapters/NexteraPE-PE.fa

# while read name; do
# 	echo "trimming sample $name"
# 	$TRIMMOMATIC PE -threads $THREADS -phred33 $DATA_DIR/"$name"_1.fastq.gz $DATA_DIR/"$name"_2.fastq.gz \
# 	$WORK_DIR/qual_trimmed/paired/"$name"_1_paired.fastq.gz $WORK_DIR/qual_trimmed/unpaired/"$name"_1_unpaired.fastq.gz \
# 	$WORK_DIR/qual_trimmed/paired/"$name"_2_paired.fastq.gz $WORK_DIR/qual_trimmed/unpaired/"$name"_2_unpaired.fastq.gz \
# 	ILLUMINACLIP:"$adapter_path":2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
# done < $SAMPLE_NAMES

#############
# ALIGNMENT #
#############
echo
echo "######################"
echo "# ALIGNMENT WITH BWA #"
echo "######################"
echo 

mkdir -p $WORK_DIR/align/index $WORK_DIR/align/sam $WORK_DIR/align/logs $WORK_DIR/align/bam

# should only need to do this once
# cp $DATA_DIR/ref_genome.fasta $WORK_DIR/align/index
# echo "generating index from reference fasta..."
# $BWA index -p $WORK_DIR/align/index/ref_genome $WORK_DIR/align/index/ref_genome.fasta

while read name; do
	echo "aligning for sample $name"
	$BWA mem -t $THREADS $WORK_DIR/align/index/ref_genome $WORK_DIR/qual_trimmed/paired/"$name"_1_paired.fastq.gz $WORK_DIR/qual_trimmed/paired/"$name"_2_paired.fastq.gz \
		1> $WORK_DIR/align/sam/"$name".sam 2> $WORK_DIR/align/logs/"$name"_log.txt
	$SAMTOOLS view -b $WORK_DIR/align/sam/"$name".sam | \
		$SAMTOOLS sort --threads $THREADS -o $WORK_DIR/align/bam/"$name"_sorted.bam && $SAMTOOLS index $WORK_DIR/align/bam/"$name"_sorted.bam
done < $SAMPLE_NAMES



echo
echo "########"
echo "# DONE #"
echo "########"
echo
