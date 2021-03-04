# Scripts

Here, all scripts are described.

## trimm_map.sh

The script loops through the fastq-files, trims off adapters using BBDuk, maps against the reference oligo pool (reference fasta and gff files can be found in [../data/reference_sequences/](../data/reference_sequences/)) and counts the reads mapped to coding sequences and sRNAs using featureCounts.

Trimming, mapping and counting statistics are stored in the log file [../analysis/logfile.log](../analysis/logfile.log) . The directory [../data/rna_align](../data/rna_align) includes all bam-alignment files as well as the count tables for our oligos [../data/rna_align/counttable.txt](../data/rna_align/counttable.txt) as well as count tables for rRNA and tRNA reads only. These will be imported into R for further downstream analysis.

## get_mismatches.sh

Here I use the seqmap tool to identify all possible zero- one- and two- mismatches which might map to other oligos in the pool.  

## getreadlengthsofalignedfiles.sh

Does exactly what its name suggests and stores the result in [../analysis/readlengths](../analysis/readlengths) . Uses samtools. 

## ribopool.Rmd

Does all downstream analysis. Further description can be found within the script.

