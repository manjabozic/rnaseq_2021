#!/bin/bash

#SBATCH --job-name=bambam
#SBATCH --partition=batch
#SBATCH --ntasks=16
#SBATCH --mem=64G
#SBATCH --time=48:00:00               # Time limit hrs:min:sec
#SBATCH --output=bambam.%j.out    # Standard output log
#SBATCH --error=bambam.%j.err     # Standard error log
#SBATCH --mail-type=ALL         #Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=mb62095@uga.edu  # Where to send mail


module load SAMtools/1.16.1-GCC-11.3.0

#FIRST, WE TRANSFORM THE OUTPUT FILES TO A BAM FORMAT AND SORT THEM ACCORDING TO CHROMOSOMAL POSITIONS

samtools view -Sb /scratch/mb62095/star/results/C_6_II/C_6_IIAligned.out.sam | samtools sort > /scratch/mb62095/star/results/bamfiles/C_6_II.bam
samtools view -Sb /scratch/mb62095/star/results/C_24_II/C_24_IIAligned.out.sam | samtools sort >/scratch/mb62095/star/results/bamfiles/ C_24_II.bam
samtools view -Sb /scratch/mb62095/star/results/C_6_III/C_6_IIIAligned.out.sam | samtools sort > /scratch/mb62095/star/results/bamfiles/C_6_III.bam
samtools view -Sb /scratch/mb62095/star/results/C_24_III/C_24_IIIAligned.out.sam | samtools sort > /scratch/mb62095/star/results/bamfiles/C_24_III.bam


samtools view -Sb /scratch/mb62095/star/results/T_6_II/T_6_IIAligned.out.sam | samtools sort > /scratch/mb62095/star/results/bamfiles/T_6_II.bam
samtools view -Sb /scratch/mb62095/star/results/T_24_II/T_24_IIAligned.out.sam | samtools sort > /scratch/mb62095/star/results/bamfiles/T_24_II.bam
samtools view -Sb /scratch/mb62095/star/results/T_6_III/T_6_IIIAligned.out.sam | samtools sort > /scratch/mb62095/star/results/bamfiles/T_6_III.bam
samtools view -Sb /scratch/mb62095/star/results/T_24_III/T_24_IIIAligned.out.sam | samtools sort > /scratch/mb62095/star/results/bamfiles/T_24_III.bam

#THEN WE CAN INDEX THE BAM FILES

samtools index /scratch/mb62095/star/results/bamfiles/C_6_II.bam
samtools index /scratch/mb62095/star/results/bamfiles/C_24_II.bam
samtools index /scratch/mb62095/star/results/bamfiles/C_6_III.bam
samtools index /scratch/mb62095/star/results/bamfiles/C_24_III.bam

samtools index /scratch/mb62095/star/results/bamfiles/T_6_II.bam
samtools index /scratch/mb62095/star/results/bamfiles/T_24_II.bam
samtools index /scratch/mb62095/star/results/bamfiles/T_6_III.bam
samtools index /scratch/mb62095/star/results/bamfiles/T_24_III.bam
