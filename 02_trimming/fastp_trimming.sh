#!/usr/bin/env bash

##This is a script for trimming thr raw sequences based on quality, as well as adapter content for paired end (PE) reads using trimmomatic.

##Create a directory where the results will be stored, or do that beforehand

	#mkdir ~/rnaseq_2021/02_trimming/results

##Then, change directories to where the fastq files are stored

cd ~/rnaseq_2021/00_raw_data

##Load the conda environment where fastqc is installed

source ~/anaconda3/etc/profile.d/conda.sh
conda activate fastp_env

## Enter the analysis parameters:

# Tell fastp what adapters you want to use

echo "Write down the path to the FASTA file containing the adapters you want to cut out. The easiest way is to have the file in the same folder as the raw data."
read -r adapters

# Set the sliding quality

echo "Quality threshold is:"
read -r q

#Set the minimum length of reads

echo "Length minimum is:"
read -r l


##Do the analysis:

##First we tell the loop what reads to use, starting with the forward reads from the directory (the ones whose names finish with "_1.fq.gz")

declare -a runtrim=(*_1.fq.gz)

#Then for each file that fits with the previous selection, do the following:

for file1 in ${runtrim[@]}

do
	file2=$(echo $file1|sed 's/_1/_2/') #This will extract the reverse reads, or all the reads whose names end with "_2.fq.gz" and add them to the variable file2:
	output=$(echo $file1|sed 's/_1.fq.gz//') #This will create the appropriate output file name;
	echo Started trimming ${file1} and ${file2} 
	fastp -i ${file1} -I ${file2} -o ~/rnaseq_2021/02_trimming/results/${output}_1P.fq.gz -O ~/rnaseq_2021/02_trimming/results/${output}_2P.fq.gz -q ${q} --detect_adapter_for_pe --adapter_fasta ${adapters} -l ${l} -f 8 -F 8 -c -p -R "${output}" -j ~/rnaseq_2021/02_trimming/results/${output}.json -h ~/rnaseq_2021/02_trimming/results/${output}.html -w 8
	echo Finished trimming ${file1} and ${file2}
done







