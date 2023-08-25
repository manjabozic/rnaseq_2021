#!/usr/bin/env bash

##This is a script for checking the quality of the raw reads using fastqc.

##Create a directory where the results will be stored or do that beforehand

	#mkdir ~/rnaseq_2021/01_quality_control/results

##Then, change directories to where the fastq files are stored

cd ~/rnaseq_2021/test_files

##Load the conda environment where fastqc is installed

source ~/anaconda3/etc/profile.d/conda.sh
conda activate fastqc

#Do the analysis

declare -a runqc=(*.fastq.gz)

for file in ${runqc[@]}

do
	echo Checking the quality of ${runqc}
	fastqc -o ~/rnaseq_2021/01_quality_control/results --noextract ${file}
done



