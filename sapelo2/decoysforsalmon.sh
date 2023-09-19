#!/bin/bash

#SBATCH --job-name=decoysforSalmon
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --mem=128G
#SBATCH --time=10:00:00               # Time limit hrs:min:sec
#SBATCH --output=decoysforSalmon.%j.out    # Standard output log
#SBATCH --error=decoysforSalmon.%j.err     # Standard error log
#SBATCH --mail-type=ALL         #Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=mb62095@uga.edu  # Where to send mail


module load BEDTools/2.30.0-GCC-12.2.0
module load 
