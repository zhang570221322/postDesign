import pickle
from ete3 import NCBITaxa
import numpy as np


# 一个k-mer的共同祖先表示是那个祖先的k-mer,然后对其挑选
k_number = 28
file_union = "/home/wlzhang/classfication/refseq/viral_database/kmer_count_k" + \
    str(k_number)+"_union"

# tax目录
ncbi = NCBITaxa()
tree = ncbi.get_topology(ncbi.get_descendant_taxa(
    10239), intermediate_nodes=True)
print("load taxonomy database done")
all_taxid = [str(taxid).replace("\n", "") for taxid in open(
    # windows TEST FILE
    "/home/wlzhang/classfication/refseq/viral_database/all_taxid.db", "r").readlines()]
merge_dic = ncbi._translate_merged(all_taxid)[1]
#  输出dic
_Total_0 = {}
for i, cur_line in enumerate(open(file_union, "r")):
    tax_count_array = np.array([tax_count.split(":") for tax_count in cur_line.replace(
        "\n", "").split("\t")[1:]], dtype=np.int64)
    for k, v in merge_dic.items():
        tax_count_array[tax_count_array[:, 0] == k, 0] = v
    tax_count_array_1 = tax_count_array[:, 0]
    size = len(tax_count_array_1)
    if size > 1:
        common_ancestor = tree.get_common_ancestor(
            list(map(str, tax_count_array_1))).name
        if common_ancestor:
            temp_set = _Total_0.setdefault(int(common_ancestor), set())
            temp_set.add(cur_line.split("\t")[
                         0]+","+str(tax_count_array[:, 1].mean()))
    if i % 1000000 == 0:
        print(i)
with open("/home/wlzhang/classfication/refseq/viral_database/com_ancestor_kmer28_1", 'wb') as f:
    pickle.dump(_Total_0, f, pickle.HIGHEST_PROTOCOL)
