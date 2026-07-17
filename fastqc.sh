#!/bin/bash
#SBATCH -p normal
#SBATCH --mem=16g
#SBATCH -N 1
#SBATCH -n 4

echo 'fastqc.sh' 
date

module load FastQC
d='./original_samples/'
samples=($(ls ${d}*.gz))
fastqc -t 4 -o $d "${samples[@]}"


echo 'done' 
date