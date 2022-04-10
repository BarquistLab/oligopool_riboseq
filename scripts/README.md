# Scripts

Here, all scripts are described.

## trimm_map.sh

The script loops through the fastq-files, trims off adapters using BBDuk, maps against the reference oligo pool (reference fasta and gff files can be found in [../data/reference_sequences/](../data/reference_sequences/)) and counts the reads mapped using featureCounts.

Trimming, mapping and counting statistics are stored in the log file [../analysis/logfile.log](../analysis/logfile.log) . The directory [../data/rna_align](../data/rna_align) includes all bam-alignment files as well as the count tables for our oligos [../data/rna_align/counttable.txt](../data/rna_align/counttable.txt) as well as count tables for rRNA and tRNA reads only. These will be imported into R for further downstream analysis. the [ref](ref) directort contains processed reference sequences for the BBtools commands. 

## get_mismatches.sh

Here I use the seqmap tool to identify all possible zero- one- and two- mismatches which might map to other oligos in the pool.  Results were later used for visualizing off-targets in the script [ribopool.Rmd](ribopool.Rmd). 

## getreadlengthsofalignedfiles.sh

Does exactly what its name suggests and stores the result in [../analysis/readlengths](../analysis/readlengths) . Uses samtools. 

## coverage.sh

bash script using bamCoverage command to map 3' ends of reads to the reference genomes to get exact position of ribosomes.

## bigwigtowig.sh

Used to convert bigwig files to wig format using bigWigToWig command.

## ribopool.Rmd

Does downstream analysis on count table. Further description can be found within the R markdown script. 

## identify_initiationsites.py

functions which help to find initiation sites using the reference oligopool and wigglefiles.

## make_oligos.py

script to generate all oligomers  from *E. coli* K12

## check_annot_start_sites.py

script used to generate predictions of start sites using a peak-calling algorithm described in the manuscript.  

## rewrite_wiggles.py

script that rewrites wiggle files to compare them with E. coli published data of meydan et al. 

## fasta_add_genenames.py

Here, I just add genenames to fasta for convenience.

## create_fasta_gff_oligopool.py

here I generate fasta and gff files for the oligopools from K12

## metagene_plots.R

Script to create metagene plots showing where TIS are. 
