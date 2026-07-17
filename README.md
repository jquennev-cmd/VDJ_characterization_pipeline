# VDJ construction and characterization pipeline
Pipeline for generating VDJ sequences from NGS data in any species. Best run on a dedicated compute cluster. 

Required softwares:
  Stand-alone tools:
    fastQC (v 0.12.1): used for verifying read qualities.
    Trimmomatic (v 0.40): used for removing low-qaulity reads and adaptor sequences
    SPADes (v 3.15.5 requires Python-3.7.2): used for generating kmer contigs from trimmed reads
    ncbi-db: pipeline uses ncbi databases for constructing reference files used in synthetic genome construction steps
    BLAST+ (v 2.9.0): used in synthetic genome construction steps
    bowtie2 (v 2.5.4): used for read alignment to synthetic genome
    SAMtools (v 1.9): used for bam file manipulation
    Subread (v 2.0.4): used in VDJ contig expression quantification
  Python Packages:
    Biopython (v)
    Pandas (v)
    
