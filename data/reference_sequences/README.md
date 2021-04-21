# oligos_ref_seqs

Here I created new fasta & gff records to only include the oligos we created, as well as trnas and rrnas.



To merge the gff with oligos with the one for trnas and rrnas, execute:

```bash
cat trnas_rrnas_mod.gff >> oligos_cds.gff 
```

