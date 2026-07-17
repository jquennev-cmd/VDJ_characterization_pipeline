#!/bin/bash
#SBATCH --job-name=spades_array
#SBATCH --output=logs/spades_%A_%a.out
#SBATCH --error=logs/spades_%A_%a.err
#SBATCH --array=0-6           # Adjust based on number of samples. 0 is first sample. 5# is max concorant jobs
#SBATCH --cpus-per-task=16
#SBATCH --mem=40G
#SBATCH -p normal

echo "Hello from array task ID: $SLURM_ARRAY_TASK_ID"
date


module swap Python/3.10.1-IGB-gcc-8.2.0 Python/3.7.2-IGB-gcc-8.2.0
module load SPAdes


sample_dir='./trimmomatic_output'
# Read sample name from list using task ID
sample=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" samples.txt)

spades.py --rna -t 16 --pe1-1 ${sample_dir}/${sample}_R1_trimmed.fastq.gz --pe1-2 ${sample_dir}/${sample}_R2_trimmed.fastq.gz -k 21,31,41,51 --only-assembler -o spades_${sample}



echo "done"
date
