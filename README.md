# INRI-seq analysis to Hör et al. 2022

- Project name: INRI-seq
- Experiments: Jens Hör, Svetlana D-M
- Supervision: Lars Barquist, Jörg Vogel
- Data analysis: Jakob Jung
- Start: January 2021

## Introduction

This project directory contains the analysis of INRI-Seq manuscript.

 

## Directory structure

The project is divided in 3 main directories:

- data: contains all raw, intermediate and final data of the project.  
- analysis: contains analysis files, such as figures, plots, tables, etc. 
- scripts: contains scripts to process and analyze data from the data directory.

Some directories have their own README.md file with information on the respective files. 



## Workflow

Here I describe the workflow, which can be followed to fully reproduce the INRI-Seq results.



### 1. Prerequisites

For running the whole IRNI-seq analysis, one needs following packages/tools/software:

- BBMap (v38.84) & BBDuk
- R (v.4.1.1) along with packages from Bioconductor/CRAN 
- Linux shell (we used Ubuntu 20.04) for commands & bash scripts
- featureCounts (v2.0.1) from Subread package
- bedtools (v2.26.0) 
- samtools (v1.12)
- bamCoverage (v3.4.3) for generation of wiggle files and 3' end mapping
- Python (v3.7) along with packages from BioPython and Anaconda/Bioconda. 



### 2. Oligo design

To design an oligo for each gene in the e. coli genome (shown in Figure 1), we downloaded the *E. coli* K12 substr. MG1655 genome and gff files from https://www.ncbi.nlm.nih.gov/nuccore/556503834 to  [./data/reference_sequences/](./data/reference_sequences/). Then we ran our custom python script to design the oligos:

```bash
python create_fasta_gff_oligopool.py
```

This script designed the complete oligopool (fasta and gff files) and saves it in the [./data/reference_sequences/](./data/reference_sequences/) directory. 

### 3. Mapping

The raw FastQ files of the INRI-Seq experiment are accessible via the Gene expresiion Omnibus under accession number: GSE190954. To re-run the analysis, the FastQ files need to be downloaded to the directory [./data/fastq](data/fastq).

All raw FastQ files should then be located in the folder [./data/fastq](data/fastq) . Details on samples and setup of the experiment can be found in the methods section of the manuscript. Navigate to [analysis/fastqc](./analysis/fastqc) to find fastQC quality statistics of the raw reads. To run the mapping, run the bash script [./scripts/trimm_map_BB.sh](./scripts/trimm_map_BB.sh):

```bash
cd scripts
sh timm_map_BB.sh
```

The script loops through the fastq-files, trims off adapters using BBDuk, maps against the reference oligo sequences (reference fasta and gff files can be found in [./data/reference_sequences/](./data/reference_sequences/)) and counts of mapped reads using featureCounts.

Trimming, mapping and counting statistics are stored in the log file [./analysis/logfile.log](./analysis/logfile.log) . The directory [./data/rna_align](./data/rna_align) includes all bam-alignment files as well as the count table [./data/rna_align/counttable.txt](./data/rna_align/counttable.txt) . This count table will be imported into R for Differential expression analysis. After mapping, bamCoverage are used to do the 3' end mapping of  the reads. bigwig files are stored in [./analysis/bamcoverage/](./analysis/bamcoverage/). The trimm_map script was run on a an external server, but it can also be run on a local computer, with an expected run-times of around 3 hours. 

We also looked at the distribution of length of reads sequenced, using the custom script:

```bash
sh getreadlengthsofalignedfiles.sh
```

It uses samtools view, perl, and the bash tools "sort", "uniq" and "sed" to loop through the aligned bam-files and store lengths of reads in   [./analysis/readlengths/](./analysis/readlengths/). Plots for read length distribution are shown in the R markdown script [ribopool.Rmd](./scripts/ribopool.Rmd). 



### 4. Coverage analysis

To run coverage analysis, we first converted bigwig-files to wig files, as wig files are easier to parse manually. In order to do this we ran the script:

```bash
sh bigwig_to_wig.sh
```

We saved resulting coverage files under [./data/wigglefiles/](./data/wigglefiles/). Next, we created metagene plots using the R script [metagene_plots.R](./scripts/metagene_plots.R)  (Figure 2A). 

For most coverage plots in the manuscript's figures (e.g. Figure 2C, Figure 3C-F, Figure S), we used the Integrative Genomics Viewer (IGV). For this, we loaded the fasta and gff files from the reference genomes ([./data/reference_sequences/](./data/reference_sequences/) ) and then loaded the respective coverage files (from [./data/wigglefiles/](./data/wigglefiles/) ) for peak visualization. 



### 5. TIS analysis

For the metagene plots, we compared our data to that of Meydan et al. (2019). We first downloaded their data from GEO (accession number: GSE122129) to [data/wigglefiles/wiggles-papercomparison_2](./data/wigglefiles/wiggles-papercomparison_2) and modified our INRI-Seq wiggle files so that they align to the *E. coli* genome using:

```bash
python rewrite_wiggles.py
```

To identify novel TIS, we wrote a custom python script  [check_annot_start_sites.py](./scripts/check_annot_start_sites.py). it can be run with:

```bash
python check_annot_start_sites.py
```

The script first checks whether annotated TIS can be found in INRI-Seq data and the RIBO-RET data from Meydan et al. Then it also screens for possibly novel TIS. Output data are stored in [./analysis/annotated_sites.csv](./analysis/annotated_sites.csv) and [./analysis/alt_tis_w_threshold.csv](./analysis/alt_tis_w_threshold.csv).  The  algorithm is described in the methods section and Figure S4 of the manuscript. We used the data to generate Figure 2D. The code to generate the figure is included in the script [ribopool.Rmd](./scripts/ribopool.Rmd). 



### 6. Differential translation analysis

To see the concentration dependent effects of PNA on the translation levels, we did differential translation analysis (Figures 4, 5). Prior to analysis, we checked for off-targets iin the oligo pool using:

```bash
sh get_mismatches.sh
```

This script uses the PNA sequence and screens for matching regions other than the target one, using seqmap. the output files are saved in  [./data/mismatches](./data/mismatches). 

To run the differential translation analysis, we used the R markdown script [./scripts/ribopool.Rmd](./scripts/ribopool.Rmd) . This outputs data for the manuscripts Figures 4,5, S2. Figures are stored in the [./analysis](./analysis) directory. The figures on the manuscript are not all in the github Repository, but they were  created using the ouput data in the repo. It can take up to 10 minutes to run this script on a low-memory laptop. 



