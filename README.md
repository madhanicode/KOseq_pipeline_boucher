# KOseq_pipeline_boucher

Contains python script, KOseq_pipeline.py, to process KO-seq fastq files and output .csv files with read counts for each junction (upstream or downstream) for each mutant.

This script was implemented in Python version 3.7.1 using the conda environment described in the KOseq_pipeline.yml file.

KOseq_pipeline.py performs two tasks.

First, "build_index" creates a custom bowtie2 index that can be used to map reads to the upstream or downstream junctions of each gene in the C. neoformans deletion library, using the expected KO junctions denoted in KO_primers_coords.csv.

Second, "align" takes KO-seq data in fastq format and uses that custom index to count reads at each junction (upstream or downstream) for each mutant.

A sample "build_index" call is:
run ./KOseq_pipeline.py build_index --input ./KO_primers_coords.csv --format csv --genome_fasta ./CNA3-gobs.fa --flank_size 300

A sample "align" call is:
run ./KOseq_pipeline.py align --index ./index/KO_primers_coords_bowtie2 --directory [path to directory containing fastq files] --threads 12
