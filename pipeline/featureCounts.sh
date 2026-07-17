#!/bin/bash
#SBATCH --job-name=FC_array
#SBATCH --output=logs/FC_array_%A_%a.out
#SBATCH --error=logs/FC_array_%A_%a.err
#SBATCH --mail-user=jordanq@illinois.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --array=0-6           # Adjust based on number of samples. 0 is first sample. add %[#] is max concorant jobs
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16 # based on https://hpc.nih.gov/apps/bowtie2
#SBATCH --mem=16G
#SBATCH -p normal

echo 'featureCounts.sh'
date

module load Subread

# make local scratch dir, accounting for job array
SCRATCH_DIR=/scratch/${USER}_${SLURM_JOB_ID}_$SLURM_ARRAY_TASK_ID
mkdir -p $SCRATCH_DIR
echo "$SCRATCH_DIR"

starting_dir=$(pwd)
saf='SPAdes_contigSynGen_mergedSamples.saf'
sample_dir='bams'
sample_l='FC_sampleList.txt'
sample=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" $sample_l)
echo $sample
out_name=${sample:0:-4}_${saf:0:-4}.'readCounts'

cp $sample_dir/$sample $SCRATCH_DIR
cp $saf $SCRATCH_DIR

cd $SCRATCH_DIR
featureCounts -p --countReadPairs -F SAF -T $SLURM_CPUS_PER_TASK -M --fraction -a $saf -o $out_name $sample

mkdir -p $starting_dir/readCounts
mv $out_name $starting_dir/readCounts

echo 'done'
date
