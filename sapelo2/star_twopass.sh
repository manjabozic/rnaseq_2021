#!/bin/bash

#SBATCH --job-name=STAR_control_2p
#SBATCH --partition=batch
#SBATCH --ntasks=16
#SBATCH --mem=128G
#SBATCH --time=72:00:00               # Time limit hrs:min:sec
#SBATCH --output=STAR_control_2p.%j.out    # Standard output log
#SBATCH --error=STAR_control_2p.%j.err     # Standard error log
#SBATCH --mail-type=ALL         #Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=mb62095@uga.edu  # Where to send mail


module load STAR/2.7.10b-GCC-11.3.0

## This is a script for aligning reads to a reference genome using STAR aligner. Make sure you know exactly where your input files are and what format they are. Also, a folder with the
## indexes has to be in the working directory when you start star.


STAR --runThreadN 16  --genomeDir /scratch/mb62095/star/indexes/ --readFilesCommand zcat --quantMode GeneCounts --outFilterMismatchNmax 15 --outFilterMultimapNmax 50 --outFileNamePrefix /scratch/mb62095/star/two_pass/C_6_III_ --readFilesIn /scratch/mb62095/data/C_6_III_1P.fq.gz /scratch/mb62095/data/C_6_III_2P.fq.gz
STAR --runThreadN 16  --genomeDir /scratch/mb62095/star/indexes/ --readFilesCommand zcat --quantMode GeneCounts --outFilterMismatchNmax 15 --outFilterMultimapNmax 50 --sjdbFileChrStartEnd /scratch/mb62095/star/two_pass/C_6_III_SJ.out.tab --outFileNamePrefix /scratch/mb62095/star/two_pass/C_6_III_ --readFilesIn /scratch/mb62095/data/C_6_III_1P.fq.gz /scratch/mb62095/data/C_6_III_2P.fq.gz

STAR --runThreadN 16  --genomeDir /scratch/mb62095/star/indexes/ --readFilesCommand zcat --quantMode GeneCounts --outFilterMismatchNmax 15 --outFilterMultimapNmax 50 --outFileNamePrefix /scratch/mb62095/star/two_pass/C_24_III_ --readFilesIn /scratch/mb62095/data/C_24_III_1P.fq.gz /scratch/mb62095/data/C_24_III_2P.fq.gz
STAR --runThreadN 16  --genomeDir /scratch/mb62095/star/indexes/ --readFilesCommand zcat --quantMode GeneCounts --outFilterMismatchNmax 15 --outFilterMultimapNmax 50 --sjdbFileChrStartEnd /scratch/mb62095/star/two_pass/C_24_III_SJ.out.tab --outFileNamePrefix /scratch/mb62095/star/two_pass/C_24_III_ --readFilesIn /scratch/mb62095/data/C_24_III_1P.fq.gz /scratch/mb62095/data/C_24_III_2P.fq.gz


STAR --runThreadN 16  --genomeDir /scratch/mb62095/star/indexes/ --readFilesCommand zcat --quantMode GeneCounts --outFilterMismatchNmax 15 --outFilterMultimapNmax 50 --outFileNamePrefix /scratch/mb62095/star/two_pass/C_6_II_ --readFilesIn /scratch/mb62095/data/C_6_II_1P.fq.gz /scratch/mb62095/data/C_6_III_2P.fq.gz
STAR --runThreadN 16  --genomeDir /scratch/mb62095/star/indexes/ --readFilesCommand zcat --quantMode GeneCounts --outFilterMismatchNmax 15 --outFilterMultimapNmax 50 --sjdbFileChrStartEnd /scratch/mb62095/star/two_pass/C_6_II_SJ.out.tab --outFileNamePrefix /scratch/mb62095/star/two_pass/C_6_II_ --readFilesIn /scratch/mb62095/data/C_6_II_1P.fq.gz /scratch/mb62095/data/C_6_II_2P.fq.gz

STAR --runThreadN 16  --genomeDir /scratch/mb62095/star/indexes/ --readFilesCommand zcat --quantMode GeneCounts --outFilterMismatchNmax 15 --outFilterMultimapNmax 50 --outFileNamePrefix /scratch/mb62095/star/two_pass/C_24_II_ --readFilesIn /scratch/mb62095/data/C_24_II_1P.fq.gz /scratch/mb62095/data/C_24_III_2P.fq.gz
STAR --runThreadN 16  --genomeDir /scratch/mb62095/star/indexes/ --readFilesCommand zcat --quantMode GeneCounts --outFilterMismatchNmax 15 --outFilterMultimapNmax 50 --sjdbFileChrStartEnd /scratch/mb62095/star/two_pass/C_24_II_SJ.out.tab --outFileNamePrefix /scratch/mb62095/star/two_pass/C_24_II_ --readFilesIn /scratch/mb62095/data/C_24_II_1P.fq.gz /scratch/mb62095/data/C_24_II_2P.fq.gz







