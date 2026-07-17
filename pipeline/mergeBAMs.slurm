#!/bin/bash
#SBATCH --job-name=bamMerge_array
#SBATCH --output=logs/bamMerge_array_%A_%a.out
#SBATCH --error=logs/bamMerge_array_%A_%a.err
#SBATCH --mail-user=jordanq@illinois.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --array=0-3           # Adjust based on number of samples. 0 is first sample. add %[#] is max concorant jobs
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16 # based on https://hpc.nih.gov/apps/bowtie2
#SBATCH --mem=60G # samtools merge is RAM heavy, assuming 30G per sorted bam file
#SBATCH -p normal

echo 'merging bams'
echo "Hello from array task ID: $SLURM_ARRAY_TASK_ID"
date

module load SAMtools

samples=($(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" bams2merge2.txt))
sample_dir='bams'
echo ${samples[@]}

##### Use the following code if samples have the same chromosome ID. 
##### Merge will not work if both bams have the same chromosome identifier
##### big ups to Pierre https://www.biostars.org/p/13462/

# new_ids=('42' '43') # chromosome name in first and second file
# i=0
# for file in ${samples[@]}
#     do
#     file=${file%.*}
#     samtools view -H "$file" | sed -e "s/SN:\([^ \t]*\)/SN:${new_ids[$i]}/" | samtools reheader - $file > ${file}_updated_header.bam
#     ((i++))
#     done

##### note this code has not been tested.
##### End of chromosome ID updater section

##### start of bam merge section 
base=${samples[0]%_SP*}
t=${samples[0]#*Gen_}
genome1=${t%%.sorted.bam}
t=${samples[1]#*Gen_}
genome2=${t%%.sorted.bam}
out_name=${base}_${genome1}+${genome2}'.sorted.bam'
echo $out_name



samtools merge -f -o "$out_name" ./${sample_dir}/${samples[0]} ./${sample_dir}/${samples[1]} -@ $SLURM_CPUS_PER_TASK
samtools index "$out_name" -@ $SLURM_CPUS_PER_TASK

##### end of bam merge section

echo 'done'
date
