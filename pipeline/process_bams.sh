#!/bin/bash
#SBATCH -p normal
#SBATCH --mem=32g
#SBATCH -N 1
#SBATCH -n 2

echo "$0"
date

module load Biopython

samples=(INSERT SAMPLE FILE PATHS HERE)
genome=[INSERT GENOME ANNOTATION FILE NAME HERE]

echo "${samples[@]}"

py_cmd="process_bam_coverage.py $genome ${samples[@]}"
py_cmd2="process_bam_expression.py readCounts" 
echo $py_cmd
python $py_cmd
echo $py_cmd2
python $py_cmd2

echo 'done'
date
