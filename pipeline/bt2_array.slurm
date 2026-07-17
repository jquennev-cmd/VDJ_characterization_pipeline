#!/bin/bash
#SBATCH --job-name=N2O_bt2_array
#SBATCH --output=logs/bt2_array_%A_%a.out
#SBATCH --error=logs/bt2_array_%A_%a.err
#SBATCH --mail-user=jordanq@illinois.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --array=0-6           # Adjust based on number of samples. 0 is first sample. add %[#] is max concorant jobs
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16 # based on https://hpc.nih.gov/apps/bowtie2
#SBATCH --mem=16G
#SBATCH -p normal

echo "Hello from array task ID: $SLURM_ARRAY_TASK_ID"
date

module swap Python/3.10.1-IGB-gcc-8.2.0 Python/3.7.2-IGB-gcc-8.2.0 
module load Bowtie2
module load SAMtools

sample=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" samples.txt)
echo $sample

# make local scratch dir, accounting for job array
SCRATCH_DIR=/scratch/${USER}_${SLURM_JOB_ID}_$SLURM_ARRAY_TASK_ID
mkdir -p $SCRATCH_DIR
echo "$SCRATCH_DIR"

# need to collect and copy over fastq, genome, genome index files to local scratch
sample_dir='trimmomatic_output'
fastq_r1=${sample}'_R1_trimmed.fastq.gz'
fastq_r2=${sample}'_R2_trimmed.fastq.gz'
genome='SPAdes_contigSynGen_newSamples'
genome_dir='/home/a-m/jordanq/ngs_dataAnalysis/pig_RNAseq/priority_sample_reanalysis/'
genome_fis=()
for ext in fa gff bt2l; do
    genome_fis+=(${genome}*.${ext})
done
# file collection complete
starting_dir=$(pwd)
echo 'copying files to scratch'
cp $starting_dir/$sample_dir/$fastq_r1 $SCRATCH_DIR
cp $starting_dir/$sample_dir/$fastq_r2 $SCRATCH_DIR

for f in "${genome_fis[@]}"
    do
    cp $starting_dir/$f $SCRATCH_DIR
    done

echo 'copying complete'

# move to right dir
cd $SCRATCH_DIR


out="${sample}_${genome}"
bt_cmd="bowtie2 --time \
-q --phred33 --very-sensitive -p 16 -I 35 -X 400 --fr -L 32 \
--mp 100,90 --rdg 100,90 --rfg 100,90 -a --no-unal \
-x $genome \
-1 $fastq_r1 \
-2 $fastq_r2 \
-S ${out}.sam " 

echo "$bt_cmd"
eval "$bt_cmd"

# now make samtools calls
samtools view -bS $out.sam -@ $SLURM_CPUS_PER_TASK > $out.bam
samtools sort $out.bam -o $out.sorted.bam -@ $SLURM_CPUS_PER_TASK
samtools index $out.sorted.bam -@ $SLURM_CPUS_PER_TASK
samtools depth $out.sorted.bam -o $out.coverage -@ $SLURM_CPUS_PER_TASK

#need to cp all outputs back to home

cp $out.sam $starting_dir
cp $out.bam $starting_dir
cp $out.sorted.bam $starting_dir
cp $out.sorted.bam.bai $starting_dir
cp $out.coverage $starting_dir

# clean up starting dir
mkdir -p $starting_dir/sam_files
mv *.sam $starting_dir/sam_files

mkdir -p $starting_dir/bams
mv *.bam* $starting_dir/bams

mkdir -p $starting_dir/coverage
mv *.coverage $starting_dir/coverage

echo 'done'
date
