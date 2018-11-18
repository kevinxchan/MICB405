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

# dependencies. assumed to be installed
FASTQC=$(which fastqc)
TRIMMOMATIC=$(which trimmomatic)
BWA=$(which bwa)
SAMTOOLS=$(which samtools)
BCFTOOLS=$(which bcftools)
VCFTOOLS=$(which vcftools)
PYTHON=$(which python3)
MAFFT=$(which mafft)
MUSCLE=$(which muscle)
TRIMAL=$(which trimal)
RAXML_NG=$(which raxml-ng)

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
echo
echo "#################################"
echo "# QUALITY TRIM WITH TRIMMOMATIC #"
echo "#################################"
echo 

mkdir -p $WORK_DIR/qual_trimmed/paired $WORK_DIR/qual_trimmed/unpaired
# some obscure location. Nextera adapters (found by looking at fastQC output)
adapter_path=/home/linuxbrew/.linuxbrew/Cellar/trinity/2.8.3/libexec/trinity-plugins/Trimmomatic-0.36/adapters/NexteraPE-PE.fa

while read name; do
	echo "trimming sample $name"
	$TRIMMOMATIC PE -threads $THREADS -phred33 $DATA_DIR/"$name"_1.fastq.gz $DATA_DIR/"$name"_2.fastq.gz \
	$WORK_DIR/qual_trimmed/paired/"$name"_1_paired.fastq.gz $WORK_DIR/qual_trimmed/unpaired/"$name"_1_unpaired.fastq.gz \
	$WORK_DIR/qual_trimmed/paired/"$name"_2_paired.fastq.gz $WORK_DIR/qual_trimmed/unpaired/"$name"_2_unpaired.fastq.gz \
	ILLUMINACLIP:"$adapter_path":2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
done < $SAMPLE_NAMES

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
cp $DATA_DIR/ref_genome.fasta $WORK_DIR/align/index
echo "generating index from reference fasta..."
$BWA index -p $WORK_DIR/align/index/ref_genome $WORK_DIR/align/index/ref_genome.fasta

while read name; do
	echo "aligning for sample $name"
	$BWA mem -t $THREADS $WORK_DIR/align/index/ref_genome $WORK_DIR/qual_trimmed/paired/"$name"_1_paired.fastq.gz $WORK_DIR/qual_trimmed/paired/"$name"_2_paired.fastq.gz \
		1> $WORK_DIR/align/sam/"$name".sam 2> $WORK_DIR/align/logs/"$name"_log.txt
	$SAMTOOLS view -b $WORK_DIR/align/sam/"$name".sam | \
		$SAMTOOLS sort --threads $THREADS -o $WORK_DIR/align/bam/"$name"_sorted.bam && $SAMTOOLS index $WORK_DIR/align/bam/"$name"_sorted.bam
	# clean up
	rm $WORK_DIR/align/sam/"$name".sam
done < $SAMPLE_NAMES

###################
# VARIANT CALLING #
###################
echo
echo "###################"
echo "# VARIANT CALLING #"
echo "###################"
echo 

echo "NOTE:"
echo "For this and the next portion of the pipeline, I got the parameters from the TB demo."
echo "In the words of Dr. Jen Gardy: 'Developing a proper filtering protocol is MSc/PhD-level serious bioinformatician stuff.'"
echo "I'm neither, and I trust her, so I'll just use these params."
echo

mkdir -p $WORK_DIR/variants/bcf $WORK_DIR/variants/vcf_var_only
cp $DATA_DIR/ref_genome.fasta $WORK_DIR/variants/bcf

while read name; do
	echo "piling up $name"
	$SAMTOOLS mpileup -q 30 -u -f $WORK_DIR/variants/bcf/ref_genome.fasta $WORK_DIR/align/bam/"$name"_sorted.bam > $WORK_DIR/variants/bcf/"$name".bcf -I
	echo "outputting variant sites for $name"
	$BCFTOOLS call -O v -mv $WORK_DIR/variants/bcf/"$name".bcf > $WORK_DIR/variants/vcf_var_only/"$name".vcf
done < $SAMPLE_NAMES

#####################
# VARIANT FILTERING #
#####################
echo
echo "#####################"
echo "# VARIANT FILTERING #"
echo "#####################"
echo 

mkdir -p $WORK_DIR/variants/vcf_filtered

while read name; do
	echo "filtering low quality variants for $name"
	# filter out DP < 5, "weird DP4", qual score >= 150
	$VCFTOOLS --vcf $WORK_DIR/variants/vcf_var_only/"$name".vcf --minDP 5.0 --stdout --recode | \
		$BCFTOOLS filter --exclude "QUAL < 150" > $WORK_DIR/variants/vcf_filtered/"$name"_DP5_QUAL150.vcf
done < $SAMPLE_NAMES

printf "checking if vcf_to_fasta_het.py exists..."
if [ -e $WORK_DIR/vcf_to_fasta_het.py ]; then
	echo "yes"
else
	echo "no"
	echo "retrieving script from github"
	echo
	wget https://raw.githubusercontent.com/jlgardy/tb_demo/master/vcf_to_fasta_het.py  -O $WORK_DIR/vcf_to_fasta_het.py
fi

echo "running custom python script for variant FASTA and tabular output..."
$PYTHON $WORK_DIR/vcf_to_fasta_het.py -x $WORK_DIR/variants/vcf_filtered/ dp5_qual150 && mv $WORK_DIR/variants/vcf_filtered/dp5_qual150* $WORK_DIR/variants/

#########################################
# MULTIPLE SEQUENCE ALIGNMENT AND TREES #
#########################################
echo
echo "###############################"
echo "# MULTIPLE SEQUENCE ALIGNMENT #"
echo "###############################"
echo 

mkdir -p $WORK_DIR/tree

echo "performing MSA with mafft and muscle..."
$MAFFT --auto --thread $THREADS $WORK_DIR/variants/dp5_qual150.fasta > $WORK_DIR/tree/dp5_qual150_mafft.mfa
$MUSCLE -in $WORK_DIR/variants/dp5_qual150.fasta -out $WORK_DIR/tree/dp5_qual150_muscle.mfa

echo "cleaning up with trimal..."
$TRIMAL -in $WORK_DIR/tree/dp5_qual150_mafft.mfa -out $WORK_DIR/tree/dp5_qual150_mafft_trimal_auto.mfa -automated1
$TRIMAL -in $WORK_DIR/tree/dp5_qual150_muscle.mfa -out $WORK_DIR/tree/dp5_qual150_muscle_trimal_auto.mfa -automated1

echo
echo "#################"
echo "# TREE BUILDING #"
echo "#################"
echo 

echo "builing tree with default params..."
$RAXML_NG --all --msa $WORK_DIR/tree/dp5_qual150_mafft_trimal_auto.mfa --model GTR+G4 --threads 1 --prefix mafft_trimal_GTR_G4_DEFAULT && mv mafft_trimal_GTR_G4_DEFAULT* $WORK_DIR/tree

echo
echo "########"
echo "# DONE #"
echo "########"
echo
