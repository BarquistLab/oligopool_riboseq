import subprocess
import re
from ushuffle import shuffle
import random
import pandas as pd


f_out = open("../data/MFE_3primeends.csv", "w")
f_out.write("TIS_type,MFE\n")

pandas_gff = pd.read_table("../data/reference_sequences/startpeaks_cds.gff", header=None)
list_entries = pandas_gff[0].to_list()

f_gff = open("../data/reference_sequences/startpeaks_cds.gff", "a")

with open("../analysis/alt_tis_w_threshold.csv", "r") as f:
    for line in f:
        column = line.split(",")
        if column[1] == "locus_tag":
            continue
        pos_from_s = int(column[5])

        # generate gff:
        lt = column[1]
        entry_gff = pandas_gff.loc[pandas_gff[0].str.contains(lt)]
        entry_gff = entry_gff.values[0].tolist()
        entry_gff = [str(i) for i in entry_gff]

        entry_gff[3] = str(pos_from_s + 58 - 20)
        entry_gff[4] = str(pos_from_s + 58 + 20)
        entry_gff[2] = "alt_tis"

        print(entry_gff)
        newline = "\t".join(entry_gff)
        print(newline)
        f_gff.write(newline + "\n")


        # do mfe stuff:
        upstream_seq = column[16]
        hoer = column[10]

        if hoer == "True" and 1 < pos_from_s < 50:
            print(type(pos_from_s), pos_from_s)
            three_prime = upstream_seq[-25:]
            print(three_prime, len(three_prime))
            cmd = " ".join(["echo", "'" + three_prime + "'", "|", "RNAfold", "|", "grep", "-Ev", "'A|U|G|C'", "|",
                            "sed", "-E", "'s/.* \\(([^\\)]+).*$/\\1/'"])
            ps = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            MFE = float(ps.communicate()[0].decode()[:-2].strip())
            print(MFE)
            f_out.write("novel" + "," + str(MFE)+"\n")

f_gff.close()

with open("../analysis/annotated_sites.csv", "r") as f:
    for line in f:
        column = line.split(",")
        if column[1] == "locus_tag":
            continue
        upstream_seq = column[17]
        hoer = column[11]
        if hoer == "True":
            print(type(pos_from_s), pos_from_s)
            three_prime = upstream_seq[-25:]
            print(three_prime, len(three_prime))
            cmd = " ".join(["echo", "'" + three_prime + "'", "|", "RNAfold", "|", "grep", "-Ev", "'A|U|G|C'", "|",
                            "sed", "-E", "'s/.* \\(([^\\)]+).*$/\\1/'"])
            ps = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            MFE = float(ps.communicate()[0].decode()[:-2].strip())
            print(MFE)
            f_out.write("annotated" + "," + str(MFE)+"\n")

f_out.close()