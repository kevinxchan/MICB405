#!/bin/bash

########################
# MAIN ANALYSIS SCRIPT #
########################
# args:
# 1. full path to working directory
# 2. full path to data directory
# 3. sample names
## these will probably never change, so i'm just going to copy paste them here.
# WORKING DIRECTORY: /home/kchan_mb18/p2
# DATA DIRECTORY: /home/kchan_mb18/SaanichInlet_200m <-- note this is a symlink!
# /home/kchan_mb18/scripts/project_2/master.sh /home/kchan_mb18/p2 /home/kchan_mb18/SaanichInlet_200m

WORK_DIR=$1
DATA_DIR=$2

# dependencies. assumed to be installed
PROKKA=$(which prokka)

THREADS=5

usage(){
echo "
MICB405: Project 2
Written by Kevin Chan

Description:  main script to run all analyses for project 2.

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

echo
echo "##################"
echo "# PRE-PROCESSING #"
echo "##################"
echo 

cd $WORK_DIR/MedQPlus_MAGs
for f in *.fa; do
	echo "$f"
done

cd $DATA_DIR/MetaBAT2_SaanichInlet_200m/gtdbtk_output/ # assumed to be here
awk '{print $1}' gtdbtk.ar122.classification_pplacer.tsv

##########
# PROKKA #
##########
echo
echo "##########"
echo "# PROKKA #"
echo "##########"
echo 

# mkdir -p $WORK_DIR/fastqc

# if [ -z "$(ls $WORK_DIR/fastqc)" ]; then
# 	echo "running fastqc for all samples..."
# 	$FASTQC -t $THREADS -o $WORK_DIR/fastqc
# else
# 	echo "detected fastqc output, skipping..."
# fi





echo
echo "########"
echo "# DONE #"
echo "########"
echo
