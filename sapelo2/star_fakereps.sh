#!/bin/bash

#SBATCH --job-name=STAR
#SBATCH --partition=batch
#SBATCH --ntasks=16
#SBATCH --mem=128G
#SBATCH --time=72:00:00               # Time limit hrs:min:sec
#SBATCH --output=STAR.%j.out    # Standard output log
#SBATCH --error=STAR.%j.err     # Standard error log
#SBATCH --mail-type=ALL         #Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=mb62095@uga.edu  # Where to send mail


module load STAR/2.7.10b-GCC-11.3.0

## This is a script for aligning reads to a reference genome using STAR aligner. Make sure you know exactly where your input files are and what format they are. Also, a folder with the
## indexes has to be in the working directory when you start star.


STAR --runThreadN 16  --genomeDir /scratch/mb62095/star/indexes/ --readFilesCommand zcat --quantMode GeneCounts --outFilterMismatchNmax 15 --outFilterMultimapNmax 50 --outFileNamePrefix /scratch/mb62095/fakerep_mapping/C_24_III_1_ --readFilesIn /scratch/mb62095/data/fake_replicates/C_24_III_1_1P.fq.gz /scratch/mb62095/data/fake_replicates/C_24_III_1_2P.fq.gz
STAR --runThreadN 16  --genomeDir /scratch/mb62095/star/indexes/ --readFilesCommand zcat --quantMode GeneCounts --outFilterMismatchNmax 15 --outFilterMultimapNmax 50 --outFileNamePrefix /scratch/mb62095/fakerep_mapping/C_24_III_2_ --readFilesIn /scratch/mb62095/data/fake_replicates/C_24_III_2_1P.fq.gz /scratch/mb62095/data/fake_replicates/C_24_III_2_2P.fq.gz
STAR --runThreadN 16  --genomeDir /scratch/mb62095/star/indexes/ --readFilesCommand zcat --quantMode GeneCounts --outFilterMismatchNmax 15 --outFilterMultimapNmax 50 --outFileNamePrefix /scratch/mb62095/fakerep_mapping/C_6_III_1_ --readFilesIn /scratch/mb62095/data/fake_replicates/C_6_III_1_1P.fq.gz /scratch/mb62095/data/fake_replicates/C_6_III_1_2P.fq.gz
STAR --runThreadN 16  --genomeDir /scratch/mb62095/star/indexes/ --readFilesCommand zcat --quantMode GeneCounts --outFilterMismatchNmax 15 --outFilterMultimapNmax 50 --outFileNamePrefix /scratch/mb62095/fakerep_mapping/C_6_III_2_ --readFilesIn /scratch/mb62095/data/fake_replicates/C_6_III_2_1P.fq.gz /scratch/mb62095/data/fake_replicates/C_6_III_2_2P.fq.gz

STAR --runThreadN 16  --genomeDir /scratch/mb62095/star/indexes/ --readFilesCommand zcat --quantMode GeneCounts --outFilterMismatchNmax 15 --outFilterMultimapNmax 50 --outFileNamePrefix /scratch/mb62095/fakerep_mapping/C_24_II_1_ --readFilesIn /scratch/mb62095/data/fake_replicates/C_24_II_1_1P.fq.gz /scratch/mb62095/data/fake_replicates/C_24_II_1_2P.fq.gz
STAR --runThreadN 16  --genomeDir /scratch/mb62095/star/indexes/ --readFilesCommand zcat --quantMode GeneCounts --outFilterMismatchNmax 15 --outFilterMultimapNmax 50 --outFileNamePrefix /scratch/mb62095/fakerep_mapping/C_24_II_2_ --readFilesIn /scratch/mb62095/data/fake_replicates/C_24_II_2_1P.fq.gz /scratch/mb62095/data/fake_replicates/C_24_II_2_2P.fq.gz
STAR --runThreadN 16  --genomeDir /scratch/mb62095/star/indexes/ --readFilesCommand zcat --quantMode GeneCounts --outFilterMismatchNmax 15 --outFilterMultimapNmax 50 --outFileNamePrefix /scratch/mb62095/fakerep_mapping/C_6_II_1_ --readFilesIn /scratch/mb62095/data/fake_replicates/C_6_II_1_1P.fq.gz /scratch/mb62095/data/fake_replicates/C_6_II_1_2P.fq.gz
STAR --runThreadN 16  --genomeDir /scratch/mb62095/star/indexes/ --readFilesCommand zcat --quantMode GeneCounts --outFilterMismatchNmax 15 --outFilterMultimapNmax 50 --outFileNamePrefix /scratch/mb62095/fakerep_mapping/C_6_II_2_ --readFilesIn /scratch/mb62095/data/fake_replicates/C_6_II_2_1P.fq.gz /scratch/mb62095/data/fake_replicates/C_6_II_2_2P.fq.gz


STAR --runThreadN 16  --genomeDir /scratch/mb62095/star/indexes/ --readFilesCommand zcat --quantMode GeneCounts --outFilterMismatchNmax 15 --outFilterMultimapNmax 50 --outFileNamePrefix /scratch/mb62095/fakerep_mapping/T_24_III_1_ --readFilesIn /scratch/mb62095/data/fake_replicates/T_24_III_1_1P.fq.gz /scratch/mb62095/data/fake_replicates/T_24_III_1_2P.fq.gz
STAR --runThreadN 16  --genomeDir /scratch/mb62095/star/indexes/ --readFilesCommand zcat --quantMode GeneCounts --outFilterMismatchNmax 15 --outFilterMultimapNmax 50 --outFileNamePrefix /scratch/mb62095/fakerep_mapping/T_24_III_2_ --readFilesIn /scratch/mb62095/data/fake_replicates/T_24_III_2_1P.fq.gz /scratch/mb62095/data/fake_replicates/T_24_III_2_2P.fq.gz
STAR --runThreadN 16  --genomeDir /scratch/mb62095/star/indexes/ --readFilesCommand zcat --quantMode GeneCounts --outFilterMismatchNmax 15 --outFilterMultimapNmax 50 --outFileNamePrefix /scratch/mb62095/fakerep_mapping/T_6_III_1_ --readFilesIn /scratch/mb62095/data/fake_replicates/T_6_III_1_1P.fq.gz /scratch/mb62095/data/fake_replicates/T_6_III_1_2P.fq.gz
STAR --runThreadN 16  --genomeDir /scratch/mb62095/star/indexes/ --readFilesCommand zcat --quantMode GeneCounts --outFilterMismatchNmax 15 --outFilterMultimapNmax 50 --outFileNamePrefix /scratch/mb62095/fakerep_mapping/T_6_III_2_ --readFilesIn /scratch/mb62095/data/fake_replicates/T_6_III_2_1P.fq.gz /scratch/mb62095/data/fake_replicates/T_6_III_2_2P.fq.gz

STAR --runThreadN 16  --genomeDir /scratch/mb62095/star/indexes/ --readFilesCommand zcat --quantMode GeneCounts --outFilterMismatchNmax 15 --outFilterMultimapNmax 50 --outFileNamePrefix /scratch/mb62095/fakerep_mapping/T_24_II_1_ --readFilesIn /scratch/mb62095/data/fake_replicates/T_24_II_1_1P.fq.gz /scratch/mb62095/data/fake_replicates/T_24_II_1_2P.fq.gz
STAR --runThreadN 16  --genomeDir /scratch/mb62095/star/indexes/ --readFilesCommand zcat --quantMode GeneCounts --outFilterMismatchNmax 15 --outFilterMultimapNmax 50 --outFileNamePrefix /scratch/mb62095/fakerep_mapping/T_24_II_2_ --readFilesIn /scratch/mb62095/data/fake_replicates/T_24_II_2_1P.fq.gz /scratch/mb62095/data/fake_replicates/T_24_II_2_2P.fq.gz
STAR --runThreadN 16  --genomeDir /scratch/mb62095/star/indexes/ --readFilesCommand zcat --quantMode GeneCounts --outFilterMismatchNmax 15 --outFilterMultimapNmax 50 --outFileNamePrefix /scratch/mb62095/fakerep_mapping/T_6_II_1_ --readFilesIn /scratch/mb62095/data/fake_replicates/T_6_II_1_1P.fq.gz /scratch/mb62095/data/fake_replicates/T_6_II_1_2P.fq.gz
STAR --runThreadN 16  --genomeDir /scratch/mb62095/star/indexes/ --readFilesCommand zcat --quantMode GeneCounts --outFilterMismatchNmax 15 --outFilterMultimapNmax 50 --outFileNamePrefix /scratch/mb62095/fakerep_mapping/T_6_II_2_ --readFilesIn /scratch/mb62095/data/fake_replicates/T_6_II_2_1P.fq.gz /scratch/mb62095/data/fake_replicates/T_6_II_2_2P.fq.gz


