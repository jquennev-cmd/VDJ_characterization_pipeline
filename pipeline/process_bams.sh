#!/bin/bash
#SBATCH -p normal
#SBATCH --mem=32g
#SBATCH -N 1
#SBATCH -n 2

echo "$0"
date

module load Biopython

samples=('./coverage/67_0dpi_S4_newSamples+originalSamples.coverage' './coverage/68_0dpi_S8_newSamples+originalSamples.coverage' './coverage/68_14dpi_S7_newSamples+originalSamples.coverage' './coverage/68_21dpi_S27_newSamples+originalSamples.coverage' './coverage/68_7dpi_S28_newSamples+originalSamples.coverage' './coverage/77_0dpc_S17_newSamples+originalSamples.coverage' './coverage/77_7dpc_S35_newSamples+originalSamples.coverage')
genome='SPAdes_contigSynGen_mergedSamples.gff'

echo "${samples[@]}"

py_cmd="process_bam_coverage.py $genome ${samples[@]}"
py_cmd2="process_bam_expression.py readCounts" 
echo $py_cmd
python $py_cmd
echo $py_cmd2
python $py_cmd2

echo 'done'
date
