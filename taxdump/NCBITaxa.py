import multiprocessing as mp
from ete3 import NCBITaxa
import math
import os
import numpy as np
import sys


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
    # windows TEST PRINT
    # print("t1:", myMatrix1[i, 0], " coe:", myMatrix1[i, 1], " t2:", myMatrix1[j, 0], " coe:",
    #       myMatrix1[j, 1], " distance:",  distance_matrix1[i, j, 0], " 带系数:", distance_matrix1[i, j, 1])


def getAver_distance(myMatrix1, myMatrix2, _stroe_Short_Path, tree):
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
        distance_matrix1 = np.full((n, m), 1, dtype=np.float32)
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
    return round(return_value_0, 3), round(return_value_1, 3)


def my_process(lines, tree, _stroe_Short_Path, i, file_union, merge_dic):
    try:
        k_mer_self_distance_0 = {}
        k_mer_self_distance_1 = {}
        for cur_line in lines:
            tax_count_array = np.array([tax_count.split(":") for tax_count in cur_line.replace(
                "\n", "").split("\t")[1:]], dtype=np.int64)
            for k, v in merge_dic.items():
                tax_count_array[tax_count_array[:, 0] == k, 0] = v
            size = len(tax_count_array)
            if size > 1:
                temp_array_0 = k_mer_self_distance_0.setdefault(size, [])
                temp_array_1 = k_mer_self_distance_1.setdefault(size, [])
                temp = getAver_distance(
                    tax_count_array, tax_count_array, _stroe_Short_Path, tree)
                temp_array_0.append(temp[0])
                temp_array_1.append(temp[1])
                # window TEST CODE
                # break
        print("Subprocess "+str(i), "DONE ,Start to write")
        with open(file_union+"1_kmer_self_distance"+str(i), "w") as f:
            f.writelines(str(k_mer_self_distance_0)+"\n")
            f.writelines(str(k_mer_self_distance_1)+"\n")
        print("Subprocess "+str(i), "DONE ,End to write")
    except Exception as identifier:
        print("ERROR, subProcess:"+str(i))
        print(identifier)


if __name__ == '__main__':
    # 多线程
    num_cores = int(sys.argv[1])
    _stroe_Short_Path = mp.Manager().dict()
    ncbi = NCBITaxa()
    tree = ncbi.get_topology(ncbi.get_descendant_taxa(
        10239), intermediate_nodes=True)
    print("load taxonomy database done")
    all_taxid = [str(taxid).replace("\n", "") for taxid in open(
        # windows TEST FILE
        # "/home/wlzhang/classfication/refseq/viral_database/all_taxid.db", "r").readlines()]
        "E:\\Github\\classification\\refSeq\\all_taxid.db", "r").readlines()]
    merge_dic = ncbi._translate_merged(all_taxid)[1]
    print("load merge  database done")
    for k_number in [26]:
        # windows TEST FILE
        file_union = "E:\\Github\\classification\\refSeq\\test\\database_fna_k26_union"
        # file_union = "/home/wlzhang/classfication/refseq/viral_database/kmer_count_k" + \
        # str(k_number)+"_union"
        lines = open(file_union, "r").readlines()
        length = len(lines)
        print("Read data end :" + file_union)
        print("multiprocessing " + str(num_cores))
        pool = mp.Pool(int(num_cores))
        for i in range(int(num_cores)):
            pool.apply_async(my_process, args=(lines[math.floor(
                i / num_cores * length):math.floor((i + 1) / num_cores * length)], tree, _stroe_Short_Path, i, file_union, merge_dic,))
        print('Waiting for all subprocesses done...')
        pool.close()
        pool.join()
        print('All subprocesses done.')
        _Total_0 = {}
        _Total_1 = {}
        for i in range(num_cores):
            with open(file_union+"1_kmer_self_distance"+str(i), "r") as f:
                lines_data = f.readlines()
                _sub_data_0 = eval(lines_data[0].replace("\n", ""))
                _sub_data_1 = eval(lines_data[1].replace("\n", ""))
                for k, v in _sub_data_0.items():
                    temp_array = _Total_0.setdefault(k, [])
                    for score in v:
                        temp_array.append(score)
                for k, v in _sub_data_1.items():
                    temp_array = _Total_1.setdefault(k, [])
                    for score in v:
                        temp_array.append(score)
                # 删除子进程文件
            os.system("rm " + file_union+"1_kmer_self_distance"+str(i))
        with open(file_union+"1_kmer_self_distance", "w") as f:
            f.writelines(str(_Total_0)+"\n")
            f.writelines(str(_Total_1)+"\n")
        print("Union result done.")
