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
META_T_DIR=/projects/micb405/resources/project\_2/2018/Metatranscriptomes/
RESOURCE_DIR=/projects/micb405/resources/project_2/2018/
ASSIGNED_DEPTH=200m

# dependencies. assumed to be installed
PYTHON3=$(which python3)
PROKKA=$(which prokka)
BWA=$(which bwa)

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

printf "checking if custom python script (match_mag_taxa.py) exists..."
if [[ ! -f $WORK_DIR/match_mag_taxa.py ]]; then
	echo "no. downloading"
	wget https://raw.githubusercontent.com/kevinxchan/MICB405/master/project_2/match_mag_taxa.py -O match_mag_taxa.py
else
	echo "yes"
fi

$PYTHON3 match_mag_taxa.py -a $DATA_DIR/MetaBAT2_SaanichInlet_200m/gtdbtk_output/gtdbtk.ar122.classification_pplacer.tsv -b $DATA_DIR/MetaBAT2_SaanichInlet_200m/gtdbtk_output/gtdbtk.bac120.classification_pplacer.tsv -m $DATA_DIR/MetaBAT2_SaanichInlet_200m/MedQPlus_MAGs -o $WORK_DIR

echo
echo "##########"
echo "# PROKKA #"
echo "##########"
echo 

# mkdir -p $WORK_DIR/prokka_output

# while read line; do
# 	name=$(echo $line | awk '{print $1}')
# 	name=${name//.fa}
# 	taxa=$(echo $line | awk '{print $2}')
# 	mkdir -p $WORK_DIR/prokka_output/"$name"

# 	echo
# 	echo "processing MAG $name"
# 	echo

# 	$PROKKA --force --outdir $WORK_DIR/prokka_output/"$name" --prefix $name --kingdom $taxa --cpus $THREADS $DATA_DIR/MetaBAT2_SaanichInlet_200m/MedQPlus_MAGs/$line
# done < $WORK_DIR/id_taxa_map.txt

echo
printf "concatenating all ORFs (.faa files) and .ffn files separately from each MAG..."

while read line; do
	name=$(echo $line | awk '{print $1}')
	name=${name//.fa}
	cat $WORK_DIR/prokka_output/$name/"$name.faa" >> "$WORK_DIR/SaanichInlet_200m_all_MAGs_ORFs.faa"
	cat $WORK_DIR/prokka_output/$name/"$name.ffn" >> "$WORK_DIR/SaanichInlet_200m_all_ref.ffn"
done < $WORK_DIR/id_taxa_map.txt
echo "done"

# at this points, user needs to upload the .faa file to KAAS and get annotations. see note below

echo
echo "####################################"
echo "# ALIGNING METATRANSCRIPTOME READS #"
echo "####################################"
echo

mkdir -p $WORK_DIR/align/index $WORK_DIR/align/sam $WORK_DIR/align/logs $WORK_DIR/align/bam
mv $WORK_DIR/SaanichInlet_200m_all_ref.ffn $WORK_DIR/align/index

echo "generating index from reference fasta..."
$BWA index -p $WORK_DIR/align/index/SaanichInlet_200m_all_ref $WORK_DIR/align/index/SaanichInlet_200m_all_ref.ffn

for f in $META_T_DIR/*"$ASSIGNED_DEPTH"*; do
	name=$(basename $f)
	name=${name//.gz}
	name=${name//.fastq}
	echo "aligning for file $name"
	$BWA mem -t $THREADS -p $WORK_DIR/align/index/SaanichInlet_200m_all_ref $f \
		1> $WORK_DIR/align/sam/"$name".sam 2> $WORK_DIR/align/logs/"$name"_log.txt
done


echo
echo "########"
echo "# DONE #"
echo "########"
echo

echo "LAST STEPS:"
echo "upload the concatenated .faa file (found at $WORK_DIR/SaanichInlet_200m_all_MAGs_ORFs.faa)"
echo "to KAAS for annotation at THIS LINK: https://www.genome.jp/kaas-bin/kaas_main?prog=GHOSTX&way=s"
echo "use all defaults. then upload this back to the server at $WORK_DIR, and save the file as"
echo "SaanichInlet_200m_all_MAGs_ORFs_ko.txt. finally, run this command:"
echo
echo "grep '\sK' $WORK_DIR/SaanichInlet_200m_all_MAGs_ORFs_ko.txt > $WORK_DIR/SaanichInlet_200m_all_MAGs_ORFs_ko.cleaned.txt"
echo
echo "and then you're golden."
echo
