from Bio.SeqIO.FastaIO import SimpleFastaParser
from collections import Counter
from functools import reduce
print("library load DONE!!")
k_number = 28


def split_dna(dna, kmer_size):
    for start in range(0, len(dna)-(kmer_size-1), 1):
        kmer = dna[start:start+kmer_size]
        yield kmer


file = "/home/wlzhang/classfication/refseq/viral_database_kraken2/library.removeNX.fna"
_kmer_count = {}
# 把数据库当成一个库, kmer作为key, 出现的taxid作为value, 格式为taxid\ttaxid\ttaxid
with open(file) as in_handle:
    for title, seq in SimpleFastaParser(in_handle):
        taxid = title.split("|")[1]
        for k_mer in split_dna(seq, k_number):
            _kmer_count[k_mer] = _kmer_count.get(k_mer, "")+"\t"+taxid
print("Union DONE!!")
with open(file + str(k_number)+"_union", "w") as f:
    for k, v in _kmer_count.items():
        f.writelines("%s\t%s\n" % (k, reduce(lambda x, y: x+"\t"+y,
                                             [k2+","+str(v2) for k2, v2 in Counter(v.split("\t")[1:]).items()])))
print("Counter and save DONE!!")
# 最终文件为这样的格式：
# CAAAGCAGCATCGGACGCTATCGCCGGT	2182394,1	2656569,1	2776873,1	2656605,1	2027893,1	1541823,1
# nohup  python Union_Counter_save.py  > Union_Counter_save.py.log 2>Union_Counter_save.py.log &
