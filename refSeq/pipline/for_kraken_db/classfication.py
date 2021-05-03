import numpy as np
# from ete3 import NCBITaxa
# from collections import Counter
from Bio.SeqIO.FastaIO import SimpleFastaParser
import pickle
import datetime
# import numpy as np


def get_max_dic(d):
    max_value = max(d.values())
    return [key for key, value in d.items() if value == max_value]


def split_dna(dna, kmer_size):
    for start in range(0, len(dna)-(kmer_size-1), 1):
        yield dna[start:start+kmer_size]


_kmer_taxid_1 = {}
with open("/home/wlzhang/classfication/refseq/viral_database_kraken2/_library_fna28_union_bin", 'rb') as f:
    _kmer_taxid_1 = pickle.load(f)


def getComplement(rev):
    rev = rev.replace("A", "X")
    rev = rev.replace("T", "A")
    rev = rev.replace("X", "T")
    rev = rev.replace("C", "X")
    rev = rev.replace("G", "C")
    rev = rev.replace("X", "G")
    return rev


k_number = 28


def test_once(_kmer_taxid, k_number):
    # performance["out"]
    hit_information = []
    with open("test.fa") as in_handle:
        for title, seq in SimpleFastaParser(in_handle):
            read_id = title.split(" ")[0]
            real_taxid = title.split("|")[0].split("reference")[1][1:]
            prediction = [read_id]
            is_Empty = True
            for kmer in split_dna(seq, k_number):
                temp = _kmer_taxid.get(kmer)
                if temp:
                    # real_taxid,count
                    prediction.append(str(temp))
                    is_Empty = False
                else:
                    temp = _kmer_taxid.get(getComplement(kmer[::-1]))
                    if temp:
                        prediction.append(str(temp))
                        is_Empty = False
                    else:
                        prediction.append("None")
            prediction.append(real_taxid)
            if is_Empty:
                print(read_id, ",", real_taxid)
                continue
            hit_information.append(prediction)
    return hit_information


hit_information1 = test_once(_kmer_taxid_1, k_number)


with open("t.csv", "w") as f:
    for prediction in hit_information1:
        f.writelines("\t".join(prediction[1:])+"\n")


# test


# ncbi = NCBITaxa()
# tree = ncbi.get_topology(ncbi.get_descendant_taxa(
#     10239), intermediate_nodes=True)
# print("load taxonomy database done")


def test():
    starttime = datetime.datetime.now()
    h = []
    label = []
    with open("t1.csv", "r") as f:
        for line in f:
            hit_lists = [eval(hit_str)
                         for hit_str in line.replace("\n", "").split("\t")]
            _taxid_count = {}
            for hit in hit_lists[:-1]:
                if hit:
                    for temp in hit:
                        taxid, count = temp.split(",")
                        _taxid_count[taxid] = _taxid_count.get(
                            taxid, 0)+int(count)
            merge_dic = ncbi._translate_merged(list(_taxid_count.keys()))[1]
            if merge_dic:
                for k, v in merge_dic.items():
                    _taxid_count[str(v)] = _taxid_count.pop(str(k))
            h.append(_taxid_count)
            label.append(hit_lists[-1])
    endtime = datetime.datetime.now()
    print((endtime - starttime).total_seconds())
    return h, label


h, label = test()


def get_max_dic(d):
    max_value = max(d.values())
    return [key for key, value in d.items() if value == max_value]


def result_rule_read(prediction_list):
    #  return the prediction Under certain rules
    _stastic = {}
    for temp in prediction_list:
        taxid, weight = temp.split(",")
        _stastic[taxid] = _stastic.get(taxid, 0)+int(weight)
    return get_max_dic(_stastic)


def is_right(prediction, real_taxid, performance, res, error, read_id):
    if len(prediction) == 0:
        # print("[empty] real taxid is %s, the res is empty %s" %
        #       (real_taxid, str(prediction)))
        performance['empty'] += 1
        res['empty'].append(real_taxid)
        error['empty'].append(read_id)
        return False
    pp = result_rule_read(prediction)
    # 如果多个预测。
    if len(pp) != 1:
        # print("[multiple] real taxid is %s, multiple predicting results %s" %
        #       (real_taxid, str(res)))
        performance['multiple'] += 1
        res['multiple'].append(pp)
        error['multiple'].append(read_id)
        return "multiple"
    if pp[0] == real_taxid:
        performance['right'] += 1
        res['right'].append(real_taxid)
        error['right'].append(read_id)
        return True
    else:
        # print("[error] real taxid is %s, error results is  %s" %
        #       (real_taxid, str(res[0])))
        performance['error'] += 1
        res['error'].append([pp[0], real_taxid])
        error['error'].append(read_id)
        return False


def handel_releative():
    performance = {"right": 0, "total": 0,
                   "empty": 0, "multiple": 0, "error": 0, "hit": 0}
    res = {"right": [],
           "empty": [], "multiple": [], "error": []}
    _read_id = {"right": [],
                "empty": [], "multiple": [], "error": [], "his_total_minus_dis_hit": []}
    with open("t1.csv", "r") as f:
        for i, line in enumerate(f):
            prediction = []
            total_hits = 0
            test = set()
            hit_lists = [eval(hit_str)
                         for hit_str in line.replace("\n", "").split("\t")]
            for hit in hit_lists[:-1]:
                if hit:
                    prediction += hit
                    performance["hit"] += len(prediction)
                    total_hits += 1
                    test.add(str(hit))
            if is_right(prediction, str(
                    hit_lists[-1]), performance, res, _read_id, i) == "multiple":
                _read_id['his_total_minus_dis_hit'].append(
                    [total_hits, len(test)])
    return performance, res, _read_id


p, r, e = handel_releative()
all_data = []
label = []
with open("t1.csv", "r") as f:
    for line in f:
        hits_list = []
        hit_lists = [eval(hit_str)
                     for hit_str in line.replace("\n", "").split("\t")]
        for hits in hit_lists[:-1]:
            if hits:
                hits_list.append([temp.split(",") for temp in hits])
            else:
                hits_list.append([['0', -1]])
        all_data.append(hits_list)
        label.append(hit_lists[-1])
with open("test.csv_bin", 'wb') as f:
    pickle.dump(all_data, f, pickle.HIGHEST_PROTOCOL)
with open("train.csv_bin", 'wb') as f:
    pickle.dump(all_data, f, pickle.HIGHEST_PROTOCOL)
with open("train.csv_bin", 'rb') as f:
    all_data = pickle.load(f)
# [一共9904行数据,每行数据123匹配,每个匹配的节点,节点的属性(两个维度)]
with open("t1.csv_bin", 'rb') as f:
    all_data = pickle.load(f)
len(label)
