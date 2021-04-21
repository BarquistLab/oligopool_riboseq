import os
import pandas as pd

ecoli_gff = pd.read_csv("../data/reference_sequences/NC_000913.3.gff", sep="\t", comment="#", header=None)


directory = '../data/wigglefiles'

for filename in os.listdir(directory):
    path = os.path.join(directory, filename)
    if path.endswith("_mod.wig"):
        print(path)
        f = open(path)
        new_forward = open(path[:-8] + "_new_F.wig", "w")
        new_reverse = open(path[:-8] + "_new_R.wig", "w")
        for line in f:
            line_splitted = line.split("\t")
            lt = line_splitted[0].split("_")[0]
            if int(line_splitted[1]) < 28 or int(line_splitted[1]) > 220:
                continue
            wanted_row = ecoli_gff[ecoli_gff[8].str.contains(lt)].iloc[0, :]
            if wanted_row[6] == "+":
                old_pos_start_end = [int(line_splitted[1]), int(line_splitted[2])]
                start = wanted_row[3]
                new_pos_start_end = a = [x - 58 - 1 + start for x in old_pos_start_end]
                new_forward.write("\t".join(["NC_000913.3", str(new_pos_start_end[0]), str(new_pos_start_end[1]),
                                             line_splitted[3]]))
            else:
                old_pos_start_end = [int(line_splitted[2]), int(line_splitted[1])]
                start = wanted_row[4]
                new_pos_start_end = a = [(x - 58)*(-1) + start for x in old_pos_start_end]
                new_reverse.write("\t".join(["NC_000913.3", str(new_pos_start_end[0]), str(new_pos_start_end[1]),
                                             str(float(line_splitted[3][:-1]) * -1) + "\n"]))
        new_forward.close()
        new_reverse.close()

    else:
        continue

# df_neg = pd.read_csv("../data/wigglefiles/GSE122129_RAW/GSM3455900_RET_BWK_U00096_3_R.wig", sep="\t", header=None, skiprows=2)
# d f_neg[1] = df_neg[1] * (-1)
# df_neg.to_csv("../data/wigglefiles/GSE122129_RAW/GSM3455900_RET_BWK_U00096_3_mod_R.wig", sep="\t", header=None, index=False)


# create wiggles:
wiggle_paper_F = pd.read_csv("../data/wigglefiles/wiggles-papercomparison_2/GSM3455900_RET_BWK_U00096_3_F_new.wig",
                             sep="\t", skiprows=2, header=None)

startsites_F = pd.read_csv("../data/reference_sequences/startsites__ecoli_F.tab",
                             sep="\t", header=None)

startsite = wiggle_paper_F.iloc[0, 0]



