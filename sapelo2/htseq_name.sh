#!/bin/bash

#SBATCH --job-name=htseq_name
#SBATCH --partition=batch
#SBATCH --ntasks=16
#SBATCH --mem=128G
#SBATCH --time=48:00:00               # Time limit hrs:min:sec
#SBATCH --output=htseq_name.%j.out    # Standard output log
#SBATCH --error=htseq_name.%j.err     # Standard error log
#SBATCH --mail-type=ALL         #Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=mb62095@uga.edu  # Where to send mail


module load HTSeq/2.0.2-foss-2022a


htseq-count --stranded=no -s no -r pos -m union -o /scratch/mb62095/fakerep_mapping/htseq/C_6_II_1.bam -p bam -c /scratch/mb62095/fakerep_mapping/htseq/C_6_II_1.tsv -n 16 /scratch/mb62095/fakerep_mapping/bamfiles_fakereps/C_6_II_1.bam /scratch/mb62095/reference/Zea_mays_B73_RefGen_v5.exon.gtf
htseq-count --stranded=no -s no -r pos -m union -o /scratch/mb62095/fakerep_mapping/htseq/C_6_II_2.bam -p bam -c /scratch/mb62095/fakerep_mapping/htseq/C_6_II_2.tsv -n 16 /scratch/mb62095/fakerep_mapping/bamfiles_fakereps/C_6_II_2.bam /scratch/mb62095/reference/Zea_mays_B73_RefGen_v5.exon.gtf

htseq-count --stranded=no -s no -r pos -m union -o /scratch/mb62095/fakerep_mapping/htseq/C_24_II_1.bam -p bam -c /scratch/mb62095/fakerep_mapping/htseq/C_24_II_1.tsv -n 16 /scratch/mb62095/fakerep_mapping/bamfiles_fakereps/C_24_II_1.bam /scratch/mb62095/reference/Zea_mays_B73_RefGen_v5.exon.gtf
htseq-count --stranded=no -s no -r pos -m union -o /scratch/mb62095/fakerep_mapping/htseq/C_24_II_2.bam -p bam -c /scratch/mb62095/fakerep_mapping/htseq/C_24_II_2.tsv -n 16 /scratch/mb62095/fakerep_mapping/bamfiles_fakereps/C_24_II_2.bam /scratch/mb62095/reference/Zea_mays_B73_RefGen_v5.exon.gtf


htseq-count --stranded=no -s no -r pos -m union -o /scratch/mb62095/fakerep_mapping/htseq/C_6_III_1.bam -p bam -c /scratch/mb62095/fakerep_mapping/htseq/C_6_III_1.tsv -n 16 /scratch/mb62095/fakerep_mapping/bamfiles_fakereps/C_6_III_1.bam /scratch/mb62095/reference/Zea_mays_B73_RefGen_v5.exon.gtf
htseq-count --stranded=no -s no -r pos -m union -o /scratch/mb62095/fakerep_mapping/htseq/C_6_III_2.bam -p bam -c /scratch/mb62095/fakerep_mapping/htseq/C_6_III_2.tsv -n 16 /scratch/mb62095/fakerep_mapping/bamfiles_fakereps/C_6_III_2.bam /scratch/mb62095/reference/Zea_mays_B73_RefGen_v5.exon.gtf

htseq-count --stranded=no -s no -r pos -m union -o /scratch/mb62095/fakerep_mapping/htseq/C_24_III_1.bam -p bam -c /scratch/mb62095/fakerep_mapping/htseq/C_24_III_1.tsv -n 16 /scratch/mb62095/fakerep_mapping/bamfiles_fakereps/C_24_III_1.bam /scratch/mb62095/reference/Zea_mays_B73_RefGen_v5.exon.gtf
htseq-count --stranded=no -s no -r pos -m union -o /scratch/mb62095/fakerep_mapping/htseq/C_24_III_2.bam -p bam -c /scratch/mb62095/fakerep_mapping/htseq/C_24_III_2.tsv -n 16 /scratch/mb62095/fakerep_mapping/bamfiles_fakereps/C_24_III_2.bam /scratch/mb62095/reference/Zea_mays_B73_RefGen_v5.exon.gtf





htseq-count --stranded=no -s no -r pos -m union -o /scratch/mb62095/fakerep_mapping/htseq/T_6_II_1.bam -p bam -c /scratch/mb62095/fakerep_mapping/htseq/T_6_II_1.tsv -n 16 /scratch/mb62095/fakerep_mapping/bamfiles_fakereps/T_6_II_1.bam /scratch/mb62095/reference/Zea_mays_B73_RefGen_v5.exon.gtf
htseq-count --stranded=no -s no -r pos -m union -o /scratch/mb62095/fakerep_mapping/htseq/T_6_II_2.bam -p bam -c /scratch/mb62095/fakerep_mapping/htseq/T_6_II_2.tsv -n 16 /scratch/mb62095/fakerep_mapping/bamfiles_fakereps/T_6_II_2.bam /scratch/mb62095/reference/Zea_mays_B73_RefGen_v5.exon.gtf

htseq-count --stranded=no -s no -r pos -m union -o /scratch/mb62095/fakerep_mapping/htseq/T_24_II_1.bam -p bam -c /scratch/mb62095/fakerep_mapping/htseq/T_24_II_1.tsv -n 16 /scratch/mb62095/fakerep_mapping/bamfiles_fakereps/T_24_II_1.bam /scratch/mb62095/reference/Zea_mays_B73_RefGen_v5.exon.gtf
htseq-count --stranded=no -s no -r pos -m union -o /scratch/mb62095/fakerep_mapping/htseq/T_24_II_2.bam -p bam -c /scratch/mb62095/fakerep_mapping/htseq/T_24_II_2.tsv -n 16 /scratch/mb62095/fakerep_mapping/bamfiles_fakereps/T_24_II_2.bam /scratch/mb62095/reference/Zea_mays_B73_RefGen_v5.exon.gtf


htseq-count --stranded=no -s no -r pos -m union -o /scratch/mb62095/fakerep_mapping/htseq/T_6_III_1.bam -p bam -c /scratch/mb62095/fakerep_mapping/htseq/T_6_III_1.tsv -n 16 /scratch/mb62095/fakerep_mapping/bamfiles_fakereps/T_6_III_1.bam /scratch/mb62095/reference/Zea_mays_B73_RefGen_v5.exon.gtf
htseq-count --stranded=no -s no -r pos -m union -o /scratch/mb62095/fakerep_mapping/htseq/T_6_III_2.bam -p bam -c /scratch/mb62095/fakerep_mapping/htseq/T_6_III_2.tsv -n 16 /scratch/mb62095/fakerep_mapping/bamfiles_fakereps/T_6_III_2.bam /scratch/mb62095/reference/Zea_mays_B73_RefGen_v5.exon.gtf

htseq-count --stranded=no -s no -r pos -m union -o /scratch/mb62095/fakerep_mapping/htseq/T_24_III_1.bam -p bam -c /scratch/mb62095/fakerep_mapping/htseq/T_24_III_1.tsv -n 16 /scratch/mb62095/fakerep_mapping/bamfiles_fakereps/T_24_III_1.bam /scratch/mb62095/reference/Zea_mays_B73_RefGen_v5.exon.gtf
htseq-count --stranded=no -s no -r pos -m union -o /scratch/mb62095/fakerep_mapping/htseq/T_24_III_2.bam -p bam -c /scratch/mb62095/fakerep_mapping/htseq/T_24_III_2.tsv -n 16 /scratch/mb62095/fakerep_mapping/bamfiles_fakereps/T_24_III_2.bam /scratch/mb62095/reference/Zea_mays_B73_RefGen_v5.exon.gtf


