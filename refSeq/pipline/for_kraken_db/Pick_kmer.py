import pickle
import math
import random
import numpy as np
# Pick_kmer.py k_number
k_number = 31
basename = "/home/wlzhang/classfication/refseq/viral_database_kraken2/"
# file_dir =
_SPECIES_kmer_taxid_count_bin_file = basename+"_TOTAL_kmer_taxid_count_bin"
_SPECIES_kmer_taxid_count = {}
with open(_SPECIES_kmer_taxid_count_bin_file, 'rb') as f:
    _SPECIES_kmer_taxid_count = pickle.load(f)
print("[DONE] _SPECIES_kmer_taxid_count")


def randomboolean(pick_divisor):


def handle(pick_divisor):
    _OUT_DIC = {}
    for k, v in _SPECIES_kmer_taxid_count.items():
        pick_number = math.ceil(len(temp_kmer)*pick_divisor)
        _OUT_DIC[int(taxid)] = random.sample(
            temp_kmer, pick_number)
        taxid = current_taxid
        temp_kmer = set()
    _OUT_DIC[int(taxid)] = random.sample(
        temp_kmer, math.ceil(len(temp_kmer)*pick_divisor))
    print(np.array([len(v) for k, v in _OUT_DIC.items()]).mean())
    # 555.8358715685323
    out_dir = basename+"taxid_kmer"+str(k_number)+"_pick_1"
    with open(out_dir+"_taxid_lines_"+str(pick_divisor), 'wb') as f:
        pickle.dump(_OUT_DIC, f, pickle.HIGHEST_PROTOCOL)


handle(0.1866)
handle(0.207)
