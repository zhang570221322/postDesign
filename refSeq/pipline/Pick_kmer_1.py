import pickle
import math
import random
import numpy as np
# Pick_kmer.py k_number
k_number = 31
basename = "/home/wlzhang/classfication/refseq/viral_database/"
file_dir = basename+"kmer_count_k" + \
    str(k_number)
kmer_count_k28_union_bin_dir = file_dir+"_union_species_bin"
_kmer_taxid_count = {}
with open(kmer_count_k28_union_bin_dir, 'rb') as f:
    _kmer_taxid_count = pickle.load(f)
print("[DONE] kmer_count_k28_union_bin")
_taxid_bps = {}
with open("taxid_bps_bin", 'rb') as f:
    _taxid_bps = pickle.load(f)


def handle(pick_divisor):
    _OUT_DIC = {}
    first = True
    taxid = None
    temp_kmer = set()
    for line in open(file_dir, "r"):
        if line.startswith("taxid"):
            current_taxid = line.split(":|:")[1]
            if first:
                taxid = current_taxid
                first = False
                continue
            pick_number = math.ceil(len(temp_kmer)*pick_divisor)
            _OUT_DIC[int(taxid)] = random.sample(
                temp_kmer, pick_number)
            taxid = current_taxid
            temp_kmer = set()
        else:
            kmer, count = line.replace("\n", "").split("\t")
            if _kmer_taxid_count.get(kmer):
                temp_kmer.add(kmer+","+count)
    _OUT_DIC[int(taxid)] = random.sample(
        temp_kmer, math.ceil(len(temp_kmer)*pick_divisor))
    print(np.array([len(v) for k, v in _OUT_DIC.items()]).mean())
    # 555.8358715685323
    out_dir = basename+"taxid_kmer"+str(k_number)+"_pick_1"
    with open(out_dir+"_taxid_lines_"+str(pick_divisor), 'wb') as f:
        pickle.dump(_OUT_DIC, f, pickle.HIGHEST_PROTOCOL)


handle(0.1866)
handle(0.207)
