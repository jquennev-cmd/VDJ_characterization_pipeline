#!/bin/bash
#SBATCH -p normal
#SBATCH --mem=32g
#SBATCH -N 1
#SBATCH -n 2
#SBATCH --job-name=syn_genome_constructor
#SBATCH --output=logs/syn_genome_constructor_%A_%a.out
#SBATCH --error=logs/syn_genome_constructor_%A_%a.err

echo 'syn_genome_constructor'
date

module load Biopython
module load BLAST+

# build samples array
mapfile -t samples < samples.txt


# NOTE: blast db needs to be updated manually in script
python contig_synGenome_generator_v2.py "${samples[@]}" newSamples


echo 'done'
date