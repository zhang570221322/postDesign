import pickle

LCA = "/home/wlzhang/classfication/refseq/viral_database_kraken2/library.removeNX.fna28_union_LCA"
SPECIES = "/home/wlzhang/classfication/refseq/viral_database_kraken2/library.removeNX.fna28_union_SPECIES"
library_fna28_union_file = "/home/wlzhang/classfication/refseq/viral_database_kraken2/library.removeNX.fna28_union"
_taxid_kmers_count = {}  # {taxid:['ATCG,count','ATCG,count','ATCG,count']}
_kmer_taxid_count = {}  # {kmer:taxid,count}
# {LCA_kmer:['2182360,1', '2250319,1']}  第一个是LCA:count
_LCAkmer_taxids_counts = {}
_SPECIESkmer_taxid_count = {}  # {SPECIES_kmer:taxid,count}
for cur_line in open(LCA, "r"):
    try:
        temp = cur_line.replace("\n", "").split("\t")
        LCA_str, LCA_count = temp[0].split(",")  # 'LCA_taxid,count'
        kmer = temp[1]  # 'ATTTCG'
        taxids = temp[2:]  # ['taxid,count','taxid:count']
        _taxid_kmers_count.setdefault(LCA_str, []).append(kmer)
        _kmer_taxid_count[kmer] = temp[0]
        _LCAkmer_taxids_counts[kmer] = [temp[0]] + taxids
    except Exception as e:
        print(e)
        print(cur_line)
for cur_line in open(SPECIES, "r"):
    try:
        temp = cur_line.replace("\n", "").split("\t")
        kmer = temp[0]  # 'ATTTCG'
        taxid, count = temp[1].split(",")  # 'taxid,count'
        _taxid_kmers_count.setdefault(taxid, []).append(kmer)
        _kmer_taxid_count[kmer] = temp[1]
        _SPECIESkmer_taxid_count[kmer] = temp[1]
    except Exception as e:
        print(e)
        print(cur_line)
_library_fna28_union = {}   # {'ATTTCG':['taxid,count','taxid,count']}
for cur_line in open(library_fna28_union_file, "r"):
    try:
        temp = cur_line.replace("\n", "").split("\t")
        _library_fna28_union[temp[0]] = temp[1:]
    except Exception as e:
        print(e)
        print(cur_line)


def save(file_name, dic):
    with open(file_name+"_bin", 'wb') as f:
        pickle.dump(dic, f, pickle.HIGHEST_PROTOCOL)


save("_library_fna28_union", _library_fna28_union)
save("_taxid_kmers_count", _taxid_kmers_count)
save("_TOTAL_kmer_taxid_count", _kmer_taxid_count)
save("_LCAk_mer_taxids_counts", _LCAkmer_taxids_counts)
save("_SPECIES_kmer_taxid_count", _SPECIESkmer_taxid_count)
