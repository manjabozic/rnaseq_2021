#!/bin/bash
#SBATCH --job-name=TRIM
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --mem=20gb
#SBATCH --time=05:00:00               # Time limit hrs:min:sec
#SBATCH --output=TRIM.%j.out
#SBATCH --error=TRIM.%j.err

#SBATCH --mail-type=NONE
#SBATCH --mail-user=djl35334@uga.edu

#module load Trimmomatic/0.39-Java-1.8.0_144

#java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE -threads 15 SRR10676751_1.fastq SRR10676751_2.fastq -baseout SRR10676751.TRIM ILLUMINACLIP:ADAPTERS.fa:2:30:10 SLIDINGWINDOW:5:20
#java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE -threads 15 SRR10676752_1.fastq SRR10676752_2.fastq -baseout SRR10676752.TRIM ILLUMINACLIP:ADAPTERS.fa:2:30:10 SLIDINGWINDOW:5:20

#module purge
module load FastQC/0.11.8-Java-1.8.0_144
fastqc SRR10676751.TRIM_1.fastq
fastqc SRR10676751.TRIM_2.fastq
fastqc SRR10676752.TRIM_1.fastq
fastqc SRR10676752.TRIM_2.fastq
