#!/bin/bash
#SBATCH --job-name=gff2SAF
#SBATCH --output=logs/gff2SAF_%A_%a.out
#SBATCH --error=logs/gff2SAF_%A_%a.err
#SBATCH -p normal
#SBATCH --mem=16g
#SBATCH -N 1
#SBATCH -n 2

echo 'generating genome SAF files'
date

genome='SPAdes_contigSynGen_mergedSamples'
echo $genome
# convert gff file to saf
python gff2SAF.py $genome

echo 'done'
date
