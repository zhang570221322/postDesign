import numpy as np
from ete3 import NCBITaxa
import multiprocessing as mp
import sys
import math


def my_setdefault(mydict, k, default):
    return mydict.get(k) if mydict.get(k) else mydict.setdefault(k, default)


def handle_martix_distance(distance_matrix1, myMatrix1, myMatrix2, tree, _stroe_Short_Path, i, j, err_count):
    # f(x_i,y_j) implements
    distance = 0
    if myMatrix1[i, 0] != myMatrix2[j, 0]:
        try:
            distance = my_setdefault(_stroe_Short_Path, str(
                myMatrix1[i, 0])+"_"+str(myMatrix2[j, 0]),  tree.get_distance(
                str(myMatrix1[i, 0]), str(myMatrix2[j, 0])))
        except Exception as identifier:
            err_count[0] += 1
            print("WARNING, Ignore the species. ", str(identifier))
    else:
        distance = 1
    distance_matrix1[i, j, 0] = distance
    distance_matrix1[i, j, 1] = (max(
        myMatrix1[i, 1], myMatrix2[j, 1])+abs(myMatrix2[j, 1]-myMatrix1[i, 1])/2) * distance


def getAver_distance(myMatrix1, myMatrix2, _stroe_Short_Path, tree):
    try:
        return_value_0, return_value_1 = -1, -1
        err_count = [0]
        # 对于自身,减少循环，加快速度， 其实用第二种也可以
        if myMatrix1 is myMatrix2:
            n = myMatrix1.shape[0]
            # with[n,m,1] or without[n,m,0] the number of k-mer appear in species
            distance_matrix1 = np.full((n, n, 2), 0, dtype=np.float32)
            for i in range(n):
                for j in range(i, n):
                    handle_martix_distance(
                        distance_matrix1, myMatrix1, myMatrix2, tree, _stroe_Short_Path, i, j, err_count)
            # windows TEST PRINT
            # print(distance_matrix1)
            a_0 = distance_matrix1[:, :, 0].sum(
            )*2 - np.trace(distance_matrix1[:, :, 0])
            a_1 = distance_matrix1[:, :, 1].sum(
            )*2 - np.trace(distance_matrix1[:, :, 1])
            b = n**2-err_count[0]
            return_value_0, return_value_1 = a_0/b, a_1/b
        else:
            # 对于任意情况
            n = myMatrix1.shape[0]
            m = myMatrix2.shape[0]
            # with[n,m,1] or without[n,m,0] the number of k-mer appear in species
            distance_matrix1 = np.full((n, m, 2), 1, dtype=np.float32)
            for i in range(n):
                for j in range(m):
                    if myMatrix1[i, 0] != myMatrix2[j, 0]:
                        handle_martix_distance(
                            distance_matrix1, myMatrix1, myMatrix2, tree, _stroe_Short_Path, i, j, err_count)
            a_0 = distance_matrix1[:, :, 0].sum()
            a_1 = distance_matrix1[:, :, 1].sum()
            b = n*m-err_count[0]
            return_value_0, return_value_1 = a_0/b, a_1/b
        return round(return_value_0, 4), round(return_value_1, 4)
    except Exception as e:
        raise e


def my_process(temp_total, _kmer_union, tree, _stroe_Short_Path, merge_dic,  _LCA_matrix):
    for LCA, temp_LCA_kmer_dic in temp_total:
        try:
            n = len(temp_LCA_kmer_dic)
            distance_matrix = np.full((n, n, 1), 0, dtype=np.float32)
            kmers = sorted(list(temp_LCA_kmer_dic.keys()))
            for i in range(n):
                for j in range(i+1, n):
                    kmer1 = kmers[i]
                    kmer2 = kmers[j]
                    union1 = np.array([tax_count.split(",")
                                       for tax_count in _kmer_union[kmer1]], dtype=np.int64)
                    union2 = np.array([tax_count.split(",")
                                       for tax_count in _kmer_union[kmer2]], dtype=np.int64)
                    for k, v in merge_dic.items():
                        union2[union2[:, 0] == k, 0] = v
                        union1[union1[:, 0] == k, 0] = v
                    # i:kmer1
                    # j:[kmer2_1,kmer2_2]
                    dis = getAver_distance(
                        union1, union2, _stroe_Short_Path, tree)[1]
                    distance_matrix[i, j] = dis
                    distance_matrix[j, i] = dis
            _LCA_matrix[LCA] = distance_matrix
        except Exception as e:
            print(e)
            print(LCA)
            break


# def test(kmers):

# n = len(kmers)
# distance_matrix = np.full((n, n, 1), 'ttt', dtype='<U10')
# for i in range(len(kmers)):
#     for j in range(i+1, len(kmers)):
#         kmer1 = kmers[i]
#         kmer2 = kmers[j]
#         distance_matrix[i][j] = kmer1+","+kmer2
#         distance_matrix[j][i] = kmer1+","+kmer2
# print(distance_matrix.T)
# 'A', 'T', 'C', 'G'
# AT AC AG
# TC TG CG
# [[['ttt' 'A,T' 'A,C' 'A,G']
#   ['A,T' 'ttt' 'T,C' 'T,G']
#   ['A,C' 'T,C' 'ttt' 'C,G']
#   ['A,G' 'T,G' 'C,G' 'ttt']]]
if __name__ == '__main__':
k_number = 28
file_union = "/home/wlzhang/classfication/refseq/viral_database_kraken2/library.fna" + \
    str(k_number)+"_union"

ncbi = NCBITaxa()
tree = ncbi.get_topology(ncbi.get_descendant_taxa(
    10239), intermediate_nodes=True)
print("load taxonomy database done")
all_taxid = set([taxid.split("|")[1] for taxid in open(
    "/home/wlzhang/classfication/refseq/viral_database_kraken2/prelim_map.txt", "r").readlines()])
merge_dic = ncbi._translate_merged(all_taxid)[1]
print("load merge  database done")

_kmer_union = {}
_LCA_kmer = {}
for line in open(file_union+"_LCA", "r"):
    temp = line.replace("\n", "").split("\t")
    LCA, LCA_count = temp[0].split(",")  # 123,1
    kmer = temp[1]  # ATCGACTG
    Union_Taxid = temp[2:]  # ['123,1','123,1','123,1']
    _kmer_union[kmer] = Union_Taxid
    temp_LCA_kmer_dic = _LCA_kmer.setdefault(LCA, {})
    temp_LCA_kmer_dic[kmer] = LCA_count
print("load DATA  done")
total = list(_LCA_kmer.items())
length = len(total)
num_cores = int(sys.argv[1])
_stroe_Short_Path = mp.Manager().dict()
_LCA_matrix = mp.Manager().dict()
print("multiprocessing " + str(num_cores))
pool = mp.Pool(int(num_cores))
for i in range(int(num_cores)):
    temp_total = total[math.floor(
        i / num_cores * length):math.floor((i + 1) / num_cores * length)]
    pool.apply_async(my_process, args=(temp_total, _kmer_union,
                                       tree, _stroe_Short_Path, merge_dic, _LCA_matrix,))
