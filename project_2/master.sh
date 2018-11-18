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
PYTHON3=$(which python3)
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

if [[ ! -f $WORK_DIR/match_mag_taxa.py ]]; then
	wget https://raw.githubusercontent.com/kevinxchan/MICB405/master/project_2/match_mag_taxa.py -O match_mag_taxa.py
fi

$PYTHON3 match_mag_taxa.py -a $DATA_DIR/MetaBAT2_SaanichInlet_200m/gtdbtk_output/gtdbtk.ar122.classification_pplacer.tsv -b $DATA_DIR/MetaBAT2_SaanichInlet_200m/gtdbtk_output/gtdbtk.bac120.classification_pplacer.tsv -m $DATA_DIR/MetaBAT2_SaanichInlet_200m/MedQPlus_MAGs -o $WORK_DIR

echo
echo "##########"
echo "# PROKKA #"
echo "##########"
echo 

mkdir -p $WORK_DIR/prokka_output

while read line; do
	name=$(echo $line | awk '{print $1}')
	name=${name//.fa}
	taxa=$(echo $line | awk '{print $2}')
	mkdir -p $WORK_DIR/prokka_output/"$name"

	echo
	echo "processing MAG $name"
	echo

	$PROKKA --force --outdir $WORK_DIR/prokka_output/"$name" --prefix $name --kingdom $taxa --cpus $THREADS $DATA_DIR/MetaBAT2_SaanichInlet_200m/MedQPlus_MAGs/$line
done < $WORK_DIR/id_taxa_map.txt




echo
echo "########"
echo "# DONE #"
echo "########"
echo
