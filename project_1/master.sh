#!/bin/bash

########################
# MAIN ANALYSIS SCRIPT #
########################
# args:
# 1. full path to working directory
# 2. full path to data directory
## these will probably never change, so i'm just going to copy paste them here.
# WORKING DIRECTORY: /home/kchan_mb18/p1
# DATA DIRECTORY: /home/kchan_mb18/project_1 <-- note this is a symlink!

WORK_DIR=$1
DATA_DIR=$2
FASTQC=$(which fastqc)
BWA=$(which bwa)

usage(){
echo "
MICB405: Project 1
Written by Kevin Chan

Description:  main script to run all analyses for project 1.

Usage:        master.sh working_dir data_dir

positional arguments:
 working_dir       full path to the working directory.
 data_dir          full path to the data directory.
"   
}

if [ -z "$1" ] || [[ $1 == -h ]] || [[ $1 == --help ]]; then
	usage
	exit
fi

echo -e "\n > starting main analysis script...\n"

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
echo -e "\n > 1. fastqc..."
mkdir -p $WORK_DIR/fastqc
if [ -z "$(ls $WORK_DIR/fastqc)" ]; then
	echo -e "\n > running fastqc for all samples\n"
	$FASTQC -t 10 -o $WORK_DIR/fastqc
else
	echo -e "\n > detected fastqc output, skipping...\n"
fi

#############
# ALIGNMENT #
#############

echo -e "\n > DONE\n"
