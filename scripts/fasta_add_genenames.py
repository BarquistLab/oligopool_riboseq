import pandas as pd

link = pd.read_csv("../data/link_lt_gn.tab", sep="\t", names=["genename", "locustag"])

old_fasta_file = open("../data/reference_sequences/oligos_cds.fasta", "r")
new_fasta_file = open("../data/reference_sequences/oligos_cds_with_name.fasta", "w")

for line in old_fasta_file:
    if line[0] != ">":
        new_fasta_file.write(line)
    else:
        splitted_line = line.split("_")
        lt = splitted_line[0][1:]
        gname = link.loc[link["locustag"] == lt]
        if len(gname) > 0:
            line = line[:-1] + "_" + gname.iloc[0, 0] + "\n"
        new_fasta_file.write(line)


old_fasta_file.close()
new_fasta_file.close()