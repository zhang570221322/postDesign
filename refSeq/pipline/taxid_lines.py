import sys
import pickle
import math
import os

# taxid_lines.py k_number
_taxid_lines = {}
k_number = sys.argv[1]
file_dir = "/home/wlzhang/classfication/refseq/viral_database/kmer_count_k" + \
    str(k_number)


def taxid_count_yeild(file_dir):
    for i, current_line in enumerate(open(file_dir, "r")):
        if current_line.startswith("taxid"):
            taxid = current_line.split(":|:")[1].replace('\n', '')
            yield taxid, i


def handle_taxid_lines_dic(file_dir):
    _loc_dic = {}
    first = True
    taxid = ""
    index_start = -1
    for current_taxid, index in taxid_count_yeild(file_dir):
        if first:
            index_start = index+1
            taxid = current_taxid
            first = False
            continue
        _loc_dic[taxid] = str(index_start)+","+str(index)
        index_start = index+1
        taxid = current_taxid
    _loc_dic[taxid] = str(index_start)+","
    return _loc_dic


# taixd:"line_start,line_end"字典
_taxid_lines = handle_taxid_lines_dic(file_dir)
with open(file_dir+"_taxid_lines", 'wb') as f:
    pickle.dump(_taxid_lines, f, pickle.HIGHEST_PROTOCOL)

#  内存化
_taxid_lines = {}
with open(file_dir+"_taxid_lines", 'rb') as f:
    _taxid_lines = pickle.load(f)
