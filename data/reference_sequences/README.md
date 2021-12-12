# reference_sequences

here, all reference sequences are stored. 

## Workflow to generate oligo sequences

1. Download the reference genome of *Escherichia coli* str. K-12 substr. MG1655 from [RefSeq](https://www.ncbi.nlm.nih.gov/nuccore/NC_000913.2) (download date: 03-11-2020). saved as: NC_000913.3.fasta , NC_000913.3.gff
2. run script [scripts/make_oligos.py](scripts/make_oligos.py) to generate 



Here I created new fasta & gff records to only include the oligos we created, as well as trnas and rrnas.

To merge the gff with oligos with the one for trnas and rrnas, execute:

```bash
cat trnas_rrnas_mod.gff >> oligos_cds.gff 
```

