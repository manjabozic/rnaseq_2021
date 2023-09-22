#!/bin/bash

#SBATCH --job-name=samtools_t6II
#SBATCH --partition=batch
#SBATCH --ntasks=16
#SBATCH --mem=128G
#SBATCH --time=48:00:00               # Time limit hrs:min:sec
#SBATCH --output=samtools_t6II.%j.out    # Standard output log
#SBATCH --error=samtools_t6II.%j.err     # Standard error log
#SBATCH --mail-type=ALL         #Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=mb62095@uga.edu  # Where to send mail


module load SAMtools/1.16.1-GCC-11.3.0

#FIRST, WE TRANSFORM THE OUTPUT FILES TO A BAM FORMAT AND SORT THEM ACCORDING TO CHROMOSOMAL POSITIONS

samtools view -b -@ 16 /scratch/mb62095/star/results/T_6_II/T_6_IIAligned.out.sam | samtools sort -@ 16 > /scratch/mb62095/star/results/bamfiles/T_6_II.bam


#THEN WE CAN INDEX THE BAM FILES

samtools index -b -@ 16 /scratch/mb62095/star/results/bamfiles/T_6_II.bam > /scratch/mb62095/star/results/bamfiles/T_6_II.bam.bai

