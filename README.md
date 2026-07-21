# VDJ construction and characterization pipeline
Pipeline for generating VDJ sequences from NGS data in any species. Best run on a dedicated compute cluster, but can be run on other achitecture with adjustment. Note: pipeline is written to run on a SLURM-based cluster.\
For questions, feel free to reach out at: jquennev@gmail.com

## Required softwares:
  ### Stand-alone tools:
   - fastQC (v 0.12.1): used for verifying read qualities.
   - Trimmomatic (v 0.40): used for removing low-qaulity reads and adaptor sequences
   - SPADes (v 3.15.5 requires Python-3.7.2): used for generating kmer contigs from trimmed reads
   - ncbi-db: pipeline uses ncbi databases for constructing reference files used in synthetic genome construction steps
   - BLAST+ (v 2.9.0): used in synthetic genome construction steps
   - bowtie2 (v 2.5.4): used for read alignment to synthetic genome
   - SAMtools (v 1.9): used for bam file manipulation
   - Subread (v 2.0.4): used in VDJ contig expression quantification
   - Python (v3.10): used for synthetic genome construction, contig expression quantification, other custom scripts 
     - Biopython (v1.76)
     - Pandas (v2.2.0)

## Step 0: Prepare working directory
- Install all required softwares.
- Make 'logs' folder in working_directory.
- Copy read files to working_directory/original_samples.
  - Note: all scripts use "original_samples" as directory name containing original sample read files.
- Make samples.txt file
  - In working_directory, make a text file of sample root names, one sample per line.
  - Used by steps 3 and onward for defining samples and subsequent files.
- Copy all script files in github pipeline folder to working_directory  

## Step 1: Read quality control
- run fastqc.sh script
  - Generates read qualtiy reports for all samples in "original_samples" folder. Use reports to remove samples from analysis and guage expectations for read trimming.

## Step 2: Read trimming
- move "nextera_adaptor.txt" file to "original_samples" folder.
  - contains adaptor sequences used for samples library prep. Update with your own sample library's adaptor sequences.
- run trim_array.slurm
  - runs trimmomatic in parallel to remove adaptor sequences and low quality read pairs.
  - Currently set up to run 7 jobs with 4 CPUs and 8G of RAM per job. Adjust to your preferences.

## Step 3: Generate kmer contigs
- run spades.slurm
  - runs SPADes as a 7 job array, 16 CPUs and 40G of RAM per job. These resources are more than sufficient for all processed samples, with only high-diversity samples requiring >10G RAM. Adjust as needed.

## Step 4: Generate synthetic genome
### Preparation step: Generate custom BLAST database
To reduce computation time when running following BLAST analyses, generating a custom BLAST database containing only your species & keywords (example: "antibody", "IGH", "immunoglobin", etc...) will significantly reduce computation times. Only needs to be run once per species and keyword pairing.
- run make_blastDB.sh
  - runs on one 1 CPU, 16G RAM.
  - provides access to ncbi database files using ncbi-db
  - Calls python script make_blastDB.py
    - NOTE: loads BLAST nr database file staticly to "db_fi" variable on line 10. Adjust as needed to point to your BLAST NR file location.
    - NOTE: species and keywords are staticly defined on lines 11 and 12. Redefine as needed for your use-case.
    - NOTE: outputs static database name, currently set to "nr_Sscrofa_Ig.fasta" on line 22. This is used statically by the next command in make_blastdb.sh. Adjust as preffered.
      - [OTHER SCRIPTS THAT USE THE DB FILE] contig_synGenome_generator_v2.py , bt2_array.slurm, merge_genomes.py, gff2SAF.py, 
### Main step: construct synthetic genome
- run syn_genome_generator.sh
- requires sample.txt file
- calls contig_synGenome_generator_v2.py script
  - NOTE: contig_synGenome_generator_v2.py requires manual update for BLAST db file. This is on line 14, "blast_db_fi" variable. Update if BLAST db file was changed.
  - NOTE: contig_synGenome_generator_v2.py contains static variables for defining intended contig size and size threshold variables, defined on lines 16 and 17. Update as needed for your use-case.
  - NOTE: synthetic genome root name is defined in syn_genome_generator.sh as the last variable in the contig_synGenome_generator_v2.py call (currently "newSamples"). Update as preferred.
- Outputs a synthetic genome fasta file and a genome annotation (gff) file.

## Step 5: Align reads to synthetic genome
- run index_syn_genome.sh
  - NOTE: requires update with synthetic genome file name on line 16.
  - NOTE: settings are optimized for following alignment step. Maximum density is recommended.
- run bt2_array.slurm
  - Launches 7 job array using 16 CPU and 16G RAM per job. Given longer compute time, also sends emails for job completion, requiring update.
  - To reduce computation times, read and genome files are copied over to the local scratch. Update this this as needed to run on your system.
  - genome root name and location variables on lines 32 & 33 need to be updated to user genome_name & genome_location

## Optional step: Merging alignments across sample sets
In the case of appending new samples to a sampleset, or if a comparison between datasets wants to be done, run the following steps to update and merge alignments. Perform this step after running previous steps on new samples and alignment of old samples to synthetic genome derived from new samples. 
- Merge genome files
  - run merge_genomes.py
    - expected call: python merge_genomes.py [GENOMES_ROOTNAMES] output_name
    - outputs a merged fasta and annotation file. Each genome becomes a separate fasta entry, or "synthetic chromosome"
- Merge BAM files 
  - make bams2merge.txt sample file
    - each line, pair of files to merge
    - run mergeBAMs.slurm
      - assumes bam files have different chromosome numbers (essential for merging)
        - if not, run commented out code in .slurm file (lines 27-34). Updates chromosome ID info
      - assumes bam files are sorted

## Step 6: Quantification of of contig expression and final output tables
### Prepatory step: generate SAF formatted genome annotation file
- run gff2SAF.sh
  - requires update to relevant genome file rootname
  - runs gff2SAF.py
### main step: Quantifying expression 
- prepare FC_samplelist.txt file
  - text file containing .bam files to be quantified, 1 file per line.
- run featureCounts.sh
  - Calculates read counts for all samples
  - launches 7 job array using 16 CPUs and 16G of RAM per job.
  - requires copying files to local scratch, update on line 19
  - requires updating line 24 to correct SAF file, generated in prepatory step
- run process_bams.sh
  - update script with samples to be processed and genome annotation file for quantification, on lines 12 and 13
  - runs process_bam_coverage.py and process_bam_expression.py
    - in process_bam_coverage.py, there are 2 variables for defining coverage regions when calculating coverage_adjusted CPM.
      - max_allowed_coverage_segments = 1
      - percent_cov_threshold = 0.8
      - change these at your discretion. these were optimized for VDJ contig coverage settings.

# Final output
Final output of the pipeline is a table of VDJ-like contigs, their nucleotide sequence, best amino acid sequence, their expression calculated in CPM, and their expression in coverage-adjusted CPM\
Coverage-adjusted CPM is CPM; but if the contig is not covered in a single area within the contig (within a predefined threshold), the expression is converted to 0. This is because we only want extremely well covered contigs for further analyses. Those with multiple blocks of coverage (or gaps in coverage) should be discarded. 

# Downstream analyses
Contigs can be further characterized using IgBLAST and a custom pig VDJ database derived from IMGT.\
Setting up custom IgBLAST databases can be found here: https://ncbi.github.io/igblast/cook/How-to-set-up.html\
The IgBLAST commands and the custom pig VDJ files are found in the IgBLAST folder.
