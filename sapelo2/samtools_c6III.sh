#!/bin/bash

#SBATCH --job-name=samtools_c6III
#SBATCH --partition=batch
#SBATCH --ntasks=16
#SBATCH --mem=128G
#SBATCH --time=48:00:00               # Time limit hrs:min:sec
#SBATCH --output=samtools_c6III.%j.out    # Standard output log
#SBATCH --error=samtools_c6III.%j.err     # Standard error log
#SBATCH --mail-type=ALL         #Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=mb62095@uga.edu  # Where to send mail


module load SAMtools/1.16.1-GCC-11.3.0

#FIRST, WE TRANSFORM THE OUTPUT FILES TO A BAM FORMAT AND SORT THEM ACCORDING TO CHROMOSOMAL POSITIONS

samtools view -b -@ 16 /scratch/mb62095/star/results/C_6_III/C_6_IIIAligned.out.sam | samtools sort -@ 16 > /scratch/mb62095/star/results/bamfiles/C_6_III.bam


#THEN WE CAN INDEX THE BAM FILES

samtools index -b -@ 16 /scratch/mb62095/star/results/bamfiles/C_6_III.bam > /scratch/mb62095/star/results/bamfiles/C_6_III.bam.bai

