#!/bin/bash

#SBATCH --job-name=stringtie
#SBATCH --partition=batch
#SBATCH --ntasks=16
#SBATCH --mem=128G
#SBATCH --time=72:00:00               # Time limit hrs:min:sec
#SBATCH --output=stringtie.%j.out    # Standard output log
#SBATCH --error=stringtie.%j.err     # Standard error log
#SBATCH --mail-type=ALL         #Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=mb62095@uga.edu  # Where to send mail


module load StringTie/2.2.1-GCC-11.2.0

stringtie /scratch/mb62095/star/results/bamfiles/C_6_II.bam -G /scratch/mb62095/reference/Zea_mays_B73_RefGen_v5.exon.gtf -p 24 -o /scratch/mb62095/star/results/bamfiles/stringtie/temp_C_6_II.gtf
stringtie /scratch/mb62095/star/results/bamfiles/C_24_II.bam -G /scratch/mb62095/reference/Zea_mays_B73_RefGen_v5.exon.gtf -p 24 -o /scratch/mb62095/star/results/bamfiles/stringtie/temp_C_24_II.gtf

stringtie /scratch/mb62095/star/results/bamfiles/C_6_III.bam -G /scratch/mb62095/reference/Zea_mays_B73_RefGen_v5.exon.gtf -p 24 -o /scratch/mb62095/star/results/bamfiles/stringtie/temp_C_6_III.gtf
stringtie /scratch/mb62095/star/results/bamfiles/C_24_III.bam -G /scratch/mb62095/reference/Zea_mays_B73_RefGen_v5.exon.gtf -p 24 -o /scratch/mb62095/star/results/bamfiles/stringtie/temp_C_24_III.gtf

stringtie /scratch/mb62095/star/results/bamfiles/T_6_II.bam -G /scratch/mb62095/reference/Zea_mays_B73_RefGen_v5.exon.gtf -p 24 -o /scratch/mb62095/star/results/bamfiles/stringtie/temp_T_6_II.gtf
stringtie /scratch/mb62095/star/results/bamfiles/T_24_II.bam -G /scratch/mb62095/reference/Zea_mays_B73_RefGen_v5.exon.gtf -p 24 -o /scratch/mb62095/star/results/bamfiles/stringtie/temp_T_24_II.gtf

stringtie /scratch/mb62095/star/results/bamfiles/T_6_III.bam -G /scratch/mb62095/reference/Zea_mays_B73_RefGen_v5.exon.gtf -p 24 -o /scratch/mb62095/star/results/bamfiles/stringtie/temp_T_6_III.gtf
stringtie /scratch/mb62095/star/results/bamfiles/T_24_III.bam -G /scratch/mb62095/reference/Zea_mays_B73_RefGen_v5.exon.gtf -p 24 -o /scratch/mb62095/star/results/bamfiles/stringtie/temp_T_24_III.gtf

stringtie --merge /scratch/mb62095/star/results/bamfiles/stringtie/*.gtf -G ~/rnaseq_2021/03_mapping/reference_genome/Zea_mays_B73_RefGen_v5.exon.gtf -p 16 -o /scratch/mb62095/star/results/bamfiles/stringtie/merged.gtf




