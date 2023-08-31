#!/usr/bin/env bash

##This is a script for trimming thr raw sequences based on quality, as well as adapter content for paired end (PE) reads using trimmomatic.

##Create a directory where the results will be stored, or do that beforehand

	#mkdir ~/rnaseq_2021/01_trimming/results

##Then, change directories to where the fastq files are stored

cd ~/rnaseq_2021/00_raw_data

##Load the conda environment where fastqc is installed

source ~/anaconda3/etc/profile.d/conda.sh
conda activate trimmomatic_env

# Java path to trimmomatic is:

#~/anaconda3/envs/trimmomatic_env/share/trimmomatic-0.39-2/trimmomatic.jar


## Enter the analysis parameters:

#First give the CPU number the process will use

echo "The number of threads you would like the tool to use:"
read -r t

# Tell trimmomatic what adapters you want to use

echo "Write down the path to the FASTA file containing the adapters you want to cut out. The easiest way is to have the file in the same folder as the raw data."
read -r adapters

# Set the number of seed mismatches (usually it is 2)

echo "The number of seed mismatches is:"
read -r seed_mismatches

# Set the palindrom clip threshold (default is 30)
echo "Palindrom cliip threshold is:"
read -r pct

# Set the simple clip threshold (default is 10)
echo "Simple cliip threshold is:"
read -r sct

# Set the sliding window info (size and quality)
echo "Sliding window size is:"
read -r sws
echo "Quality threshold is:"
read -r q


##Do the analysis:

##First we tell the loop what reads to use, starting with the forward reads from the directory (the ones whose names finish with "_1.fq.gz")

declare -a runtrim=(*_1.fq.gz)

#Then for each file that fits with the previous selection, do the following:

for file1 in ${runtrim[@]}

do
	file2=$(echo $file1|sed 's/_1/_2/') #This will extract the reverse reads, or all the reads whose names end with "_2.fq.gz" and add them to the variable file2:
	output=$(echo $file1|sed 's/_1.fq.gz//') #This will create the appropriate output file name;
	echo Started trimming ${file1} and ${file2} 
	java -jar ~/anaconda3/envs/trimmomatic_env/share/trimmomatic-0.39-2/trimmomatic.jar PE -threads ${t} ${file1} ${file2} -baseout ~/rnaseq_2021/02_trimming/results/${output}.fq.gz ILLUMINACLIP:${adapters}:${seed_mismatches}:${pct}:${sct} SLIDINGWINDOW:${sws}:${q} MINLEN:30
	echo Finished trimming ${file1} and ${file2}
done





