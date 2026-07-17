#!/bin/bash
#SBATCH --job-name=trim_array
#SBATCH --output=logs/trim_%A_%a.out
#SBATCH --error=logs/trim_%A_%a.err
#SBATCH --array=0-6%3           # Adjust based on number of samples. 0 is first sample. 5# is max concorant jobs
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH -p normal

echo "Hello from array task ID: $SLURM_ARRAY_TASK_ID"
date

module swap Java/15.0.1 Java/1.8.0_201
module load Trimmomatic

# Read sample name from list using task ID
sample=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" samples.txt)
echo $sample
# Define input/output paths
d="./original_samples"
out_dir="./trimmomatic_output"
mkdir -p "$out_dir"

r1="$d/${sample}_R1_001.fastq.gz"
r2="$d/${sample}_R2_001.fastq.gz"

out_r1="$out_dir/${sample}_R1_trimmed.fastq.gz"
out_r2="$out_dir/${sample}_R2_trimmed.fastq.gz"
unp_r1="$out_dir/${sample}_R1_unpaired.fastq.gz"
unp_r2="$out_dir/${sample}_R2_unpaired.fastq.gz"

# Run Trimmomatic
trimmomatic PE -threads 4 -phred33 \
  "$r1" "$r2" \
  "$out_r1" "$unp_r1" \
  "$out_r2" "$unp_r2" \
  ILLUMINACLIP:./original_samples/nextera_adaptor.txt:2:30:10:1:true MINLEN:40

# trim_cmd="trimmomatic PE -threads 4 \
#   "$r1" "$r2" \
#   "$out_r1" "$unp_r1" \
#   "$out_r2" "$unp_r2" \
#   ILLUMINACLIP:adaptor.txt:2:30:10:1:true MINLEN:40"
# echo $trim_cmd
# $trim_cmd


echo 'done'
date