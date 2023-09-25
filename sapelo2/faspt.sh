#!/bin/bash

#SBATCH --job-name=fastp
#SBATCH --partition=batch
#SBATCH --ntasks=16
#SBATCH --mem=128G
#SBATCH --time=48:00:00               # Time limit hrs:min:sec
#SBATCH --output=fastp.%j.out    # Standard output log
#SBATCH --error=fastp.%j.err     # Standard error log
#SBATCH --mail-type=ALL         #Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=mb62095@uga.edu  # Where to send mail


module load fastp/0.23.2-GCC-11.2.0

fastp -i /scratch/mb62095/raw_data/C_6_II_1.fq.gz -I /scratch/mb62095/raw_data/C_6_II_2.fq.gz -o /scratch/mb62095/trimming/C_6_II_1P.fq.gz -O /scratch/mb62095/trimming/C_6_II_2P.fq.gz -q 30 --detect_adapter_for_pe --adapter_fasta /scratch/mb62095/trimming/adapters.fa -l 30 -f 8 -F 8 -c -p -R "C_6_II" -j /scratch/mb62095/trimming/C_6_II.json -h /scratch/mb62095/trimming/C_6_II.html -w 16
fastp -i /scratch/mb62095/raw_data/C_24_II_1.fq.gz -I /scratch/mb62095/raw_data/C_24_II_2.fq.gz -o /scratch/mb62095/trimming/C_24_II_1P.fq.gz -O /scratch/mb62095/trimming/C_24_II_2P.fq.gz -q 30 --detect_adapter_for_pe --adapter_fasta /scratch/mb62095/trimming/adapters.fa -l 30 -f 8 -F 8 -c -p -R "C_24_II" -j /scratch/mb62095/trimming/C_24_II.json -h /scratch/mb62095/trimming/C_24_II.html -w 16
fastp -i /scratch/mb62095/raw_data/C_6_III_1.fq.gz -I /scratch/mb62095/raw_data/C_6_III_2.fq.gz -o /scratch/mb62095/trimming/C_6_III_1P.fq.gz -O /scratch/mb62095/trimming/C_6_III_2P.fq.gz -q 30 --detect_adapter_for_pe --adapter_fasta /scratch/mb62095/trimming/adapters.fa -l 30 -f 8 -F 8 -c -p -R "C_6_III" -j /scratch/mb62095/trimming/C_6_III.json -h /scratch/mb62095/trimming/C_6_III.html -w 16
fastp -i /scratch/mb62095/raw_data/C_24_III_1.fq.gz -I /scratch/mb62095/raw_data/C_24_III_2.fq.gz -o /scratch/mb62095/trimming/C_24_III_1P.fq.gz -O /scratch/mb62095/trimming/C_24_III_2P.fq.gz -q 30 --detect_adapter_for_pe --adapter_fasta /scratch/mb62095/trimming/adapters.fa -l 30 -f 8 -F 8 -c -p -R "C_24_III" -j /scratch/mb62095/trimming/C_24_III.json -h /scratch/mb62095/trimming/C_24_III.html -w 16

fastp -i /scratch/mb62095/raw_data/T_6_II_1.fq.gz -I /scratch/mb62095/raw_data/T_6_II_2.fq.gz -o /scratch/mb62095/trimming/T_6_II_1P.fq.gz -O /scratch/mb62095/trimming/T_6_II_2P.fq.gz -q 30 --detect_adapter_for_pe --adapter_fasta /scratch/mb62095/trimming/adapters.fa -l 30 -f 8 -F 8 -c -p -R "T_6_II" -j /scratch/mb62095/trimming/T_6_II.json -h /scratch/mb62095/trimming/T_6_II.html -w 16
fastp -i /scratch/mb62095/raw_data/T_24_II_1.fq.gz -I /scratch/mb62095/raw_data/T_24_II_2.fq.gz -o /scratch/mb62095/trimming/T_24_II_1P.fq.gz -O /scratch/mb62095/trimming/T_24_II_2P.fq.gz -q 30 --detect_adapter_for_pe --adapter_fasta /scratch/mb62095/trimming/adapters.fa -l 30 -f 8 -F 8 -c -p -R "T_24_II" -j /scratch/mb62095/trimming/T_24_II.json -h /scratch/mb62095/trimming/T_24_II.html -w 16
fastp -i /scratch/mb62095/raw_data/T_6_III_1.fq.gz -I /scratch/mb62095/raw_data/T_6_III_2.fq.gz -o /scratch/mb62095/trimming/T_6_III_1P.fq.gz -O /scratch/mb62095/trimming/T_6_III_2P.fq.gz -q 30 --detect_adapter_for_pe --adapter_fasta /scratch/mb62095/trimming/adapters.fa -l 30 -f 8 -F 8 -c -p -R "T_6_III" -j /scratch/mb62095/trimming/T_6_III.json -h /scratch/mb62095/trimming/T_6_III.html -w 16
fastp -i /scratch/mb62095/raw_data/T_24_III_1.fq.gz -I /scratch/mb62095/raw_data/T_24_III_2.fq.gz -o /scratch/mb62095/trimming/T_24_III_1P.fq.gz -O /scratch/mb62095/trimming/T_24_III_2P.fq.gz -q 30 --detect_adapter_for_pe --adapter_fasta /scratch/mb62095/trimming/adapters.fa -l 30 -f 8 -F 8 -c -p -R "T_24_III" -j /scratch/mb62095/trimming/T_24_III.json -h /scratch/mb62095/trimming/T_24_III.html -w 16