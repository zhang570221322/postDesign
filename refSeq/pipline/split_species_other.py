k_number = 31
file_dir = "/home/wlzhang/classfication/refseq/viral_database/kmer_count_k" + \
    str(k_number)+"_union"
f_species = open(file_dir+"_species", "w")
f_others = open(file_dir+"_others", "w")
for line in open(file_dir, "r"):
    if len(line.split("\t")) > 2:
        f_others.writelines(line)
    else:
        f_species.writelines(line)
