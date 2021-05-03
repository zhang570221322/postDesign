from ete3 import NCBITaxa
import numpy as np
import multiprocessing as mp
import sys
import math
from functools import reduce
# 一个k-mer的共同祖先表示是那个祖先的k-mer,然后对其挑选


def my_process(line_temp1, file_union1, thread_number1, merge_dic1, tree1):
    try:
        f_lca = open(file_union1+"_LCA"+str(thread_number1), "w")
        f_species = open(file_union1+"_SPECIES"+str(thread_number1), "w")
        for cur_line in line_temp1:
            tax_count_array = np.array([tax_count.split(",") for tax_count in cur_line.replace(
                "\n", "").split("\t")[1:]], dtype=np.int64)
            if len(tax_count_array) > 1:
                for k, v in merge_dic1.items():
                    tax_count_array[tax_count_array[:, 0] == k, 0] = v
                common_ancestor = tree1.get_common_ancestor(
                    list(map(str, tax_count_array[:, 0]))).name
                if common_ancestor:
                    f_lca.writelines("%s,%s\t%s" % (str(common_ancestor), str(
                        tax_count_array[:, 1].sum()), cur_line))
            else:
                f_species.writelines(cur_line)
    except Exception as e:
        print(e)


if __name__ == '__main__':
    ncbi = NCBITaxa()
    tree = ncbi.get_topology(ncbi.get_descendant_taxa(
        10239), intermediate_nodes=True)
    print("load taxonomy database done")
    all_taxid = set([taxid.split("|")[1] for taxid in open(
        "/home/wlzhang/classfication/refseq/viral_database_kraken2/prelim_map.txt", "r").readlines()])
    merge_dic = ncbi._translate_merged(all_taxid)[1]
    print("load merge  database done")
    num_cores = int(sys.argv[1])
    for k_number in [28]:
        file_union = "/home/wlzhang/classfication/refseq/viral_database_kraken2/library.removeNX.fna" + \
            str(k_number)+"_union"
        lines = open(file_union, "r").readlines()
        length = len(lines)
        print("Read data end :" + file_union)
        print("multiprocessing " + str(num_cores))
        pool = mp.Pool(int(num_cores))
        for i in range(int(num_cores)):
            line_temp = lines[math.floor(
                i / num_cores * length):math.floor((i + 1) / num_cores * length)]
            pool.apply_async(my_process, args=(
                line_temp, file_union, i, merge_dic, tree,))
        print('Waiting for all subprocesses done...')
        pool.close()
        pool.join()
        print('All subprocesses done.')
        print("cat %s > %s" % (reduce(lambda x, y: x+" "+y, [file_union+"_LCA"+str(thread_number)
                                                             for thread_number in range(num_cores)]), file_union+"_LCA"))
        print("cat %s > %s" % (reduce(lambda x, y: x+" "+y, [file_union+"_SPECIES"+str(thread_number)
                                                             for thread_number in range(num_cores)]), file_union+"_LCA"))
# nohup  python Com_ancestor.py 24  > Com_ancestor.py.log 2>Com_ancestor.py.log &
