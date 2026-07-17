#!/bin/bash
#SBATCH --job-name=synGenomeIndex
#SBATCH --output=logs/synGenomeIndex_%A_%a.out
#SBATCH --error=logs/synGenomeIndex_%A_%a.err
#SBATCH -N 1
#SBATCH -n 2
#SBATCH --mem=32G
#SBATCH -p normal

echo 'indexing synthetic genome'
date

module swap Python/3.10.1-IGB-gcc-8.2.0 Python/3.7.2-IGB-gcc-8.2.0
module load Bowtie2

genome='SPAdes_contigSynGen_newSamples'

bowtie2-build $genome.fa $genome --offrate 1 --large-index -p 2

echo 'done'
date
