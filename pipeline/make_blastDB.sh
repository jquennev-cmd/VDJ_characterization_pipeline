#!/bin/bash
#SBATCH --job-name=make_blastdb
#SBATCH --output=logs/make_blastdb_%A.out
#SBATCH --error=logs/make_blastdb_%A.err
#SBATCH -p normal
#SBATCH --mem=16g
#SBATCH -N 1
#SBATCH -n 2


module load Biopython
module load ncbi-db
module load BLAST+

# requires manual updating for db file location and species/keyword searches
python make_blastdb.py

# make sure to update this command with the correct fasta file
makeblastdb -in nr_Sscrofa_Ig.fasta -parse_seqids -blastdb_version 5 -title "Cookbook demo" -dbtype prot
