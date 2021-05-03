import sys
import pickle

#   union2bindb.py
_DIC1 = {}
k_number = 31
file_union = "/home/wlzhang/classfication/refseq/viral_database/kmer_count_k" + \
    str(k_number)+"_union_species"
for cur_line in open(file_union, "r"):
    try:
        temp = cur_line.replace("\n", "").split("\t")
        k = temp[0]
        v = temp[1:][0]
        _DIC1[k] = v
    except Exception as e:
        print(e)
        print(cur_line)
with open(file_union+"_bin", 'wb') as f:
    pickle.dump(_DIC1, f, pickle.HIGHEST_PROTOCOL)
print("DONE species")
_DIC1 = {}
file_union = "/home/wlzhang/classfication/refseq/viral_database/kmer_count_k" + \
    str(k_number)+"_union_others"
for cur_line in open(file_union, "r"):
    try:
        temp = cur_line.replace("\n", "").split("\t")
        k = temp[0]
        v = temp[1:]
        _DIC1[k] = v
    except Exception as e:
        print(e)
        print(cur_line)
with open(file_union+"_bin", 'wb') as f:
    pickle.dump(_DIC1, f, pickle.HIGHEST_PROTOCOL)
print("DONE others")


_DIC1 = {}
for cur_line in open("/home/wlzhang/classfication/refseq/viral_database/kmer_count_k" + str(k_number)+"_union", "r"):
    try:
        temp = cur_line.replace("\n", "").split("\t")
        k = temp[0]
        v = temp[1:]
        _DIC1[k] = v
    except Exception as e:
        print(e)
        print(cur_line)
with open("/home/wlzhang/classfication/refseq/viral_database/kmer_count_k" + str(k_number)+"_union"+"_bin", 'wb') as f:
    pickle.dump(_DIC1, f, pickle.HIGHEST_PROTOCOL)
