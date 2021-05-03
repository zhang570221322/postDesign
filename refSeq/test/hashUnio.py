from efficient_apriori import apriori
import sys
_Total_dic = {}
file_name = "/home/wlzhang/classfication/refseq/viral_database/kmer_count_k26"

# file_name = "E:/Github/classification/refSeq/test/test"


# def handle_muti_in_name(l: list) -> list:
#     # 因为我用的是:分割，所以遇到如taxid:429564,species_taxid:77811,id:GCF_000869345.1,name:Euphorbia mosaic virus - A [Mexico:Yucatan:2004]
#     # 的情况无法正常转为dict
#     temp = l[0:3][:]
#     temp.append(list((l[3:][0][0], ":".join(l[3:][0][1:]))))
#     return temp


process_count = 0
# file_name = "refSeq/test/test"
with open(file_name, "r") as f:
    current_line = f.readline()
    tax_id = ""
    while current_line:
        if current_line.startswith("taxid"):
            tax_id = current_line.split(":|:")[1].replace('\n', '')
            process_count += 1
            if process_count % 1000 == 0:
                print(process_count)
        else:
            k_mer, k_mer_count = current_line.split("\t")
            # {k_mer: {taxid: count, taxid: count, taxid: count}, }的形式  mode_1
            # 如果_Total_dic不存在k_mer, 设置为{}, 并返回 get(k_mer), 如果存在 get(k_mer)
            _taxId_K_mer_count = _Total_dic.setdefault(k_mer, {})
            # 取得_taxId_K_mer_count[tax_id]与当前k_mer_count比较, 更新为最大的. 解决1个taxid对应多个GCF文件
            _taxId_K_mer_count[tax_id] = max(
                int(k_mer_count),  _taxId_K_mer_count.get(tax_id, -1))
            # {k_mer:[[taxid_1,taxid_2],[count_1,count_2]]}的形式
            # K_mer_count_list = _Total_dic.setdefault(k_mer, [[], []])
            # K_mer_count_list[0].append(tax_id)
            # K_mer_count_list[1].append(k_mer_count)
        current_line = f.readline()
with open(file_name+"_union", "w") as f:
    for k, v in _Total_dic.items():
        # for mode_1
        f.write(k+"\t"+"\t".join([k2+":"+str(v2)
                                  for k2, v2 in v.items()])+"\n")
        # f.write(k+"\t"+"\t".join([v[0][i]+":"+v[1][i]
        #                           for i in range(len(v[0]))]))
# unio数据统计
k_number = sys.argv[1]
file_union = "/home/wlzhang/classfication/refseq/viral_database/kmer_count_k" + \
    str(k_number)+"_union"
_count_taxid = {}
_count_kmer = {}
with open("E:\\Github\\classification\\refSeq\\test\\database_fna_k26_union", "r") as f:
    # with open(file_union, "r") as f:
    current_line = f.readline()
    while current_line:
        kmer_taxid_list = current_line.replace("\n", "").split("\t")[1:]
        _count_taxid[len(kmer_taxid_list)] = _count_taxid.get(
            len(kmer_taxid_list), 0)+1
        if len(kmer_taxid_list) > 6:
            print(current_line)
            break
        current_line = f.readline()
        taxid_count = sum([int(temp1.split(":")[1])
                           for temp1 in kmer_taxid_list])
        _count_kmer[taxid_count] = _count_kmer.get(
            taxid_count, 0)+1
        current_line = f.readline()
with open(file_union+"_statistics", "w") as f:
    f.writelines(str(_count_taxid)+"\n")
    f.writelines(str(_count_kmer)+"\n")

_data = {}
for i in [10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 31]:
    i = str(i)
    with open("/home/wlzhang/classfication/refseq/viral_database/kmer_count_k"+i+"_union_statistics") as f:
        lines = f.readlines()
        _data["k"+i+"_1"] = eval(lines[0].replace("\n", ""))
        _data["k"+i+"_2"] = eval(lines[1].replace("\n", ""))

l1 = []
l2 = []
l3 = []
for i in [10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 31]:
    k = str(i)
    b = sum([int(v) for k2, v in _data['k'+k+'_2'].items()])
    print("--------------k="+k+"-------------")
    print("dis_kmer:", b)
    a = sum([k2*v for k2, v in _data['k'+k+'_2'].items()])/b
    print("average:", a)
    c = sum([int(v) for k2, v in _data['k'+k+'_1'].items() if int(k2) > 50])
    print("k>1000 count:", c)
    l1.append(b)
    l2.append(a)
    l3.append(c)
print(l1)
print(l2)
print(l3)

data_set = []
k_number = str(sys.argv[1])
with open("/home/wlzhang/classfication/refseq/viral_database/kmer_count_k"+k_number+"_union", "r") as f:
    current_line = f.readline()
    while current_line:
        data_set.append(tuple(set(
            [temp.split(":")[0] for temp in current_line.replace("\n", "").split("\t")[1:]])))
        current_line = f.readline()
itemsets, rules = apriori(data_set, min_support=float(
    sys.argv[2]), min_confidence=0.95)
open("/home/wlzhang/classfication/refseq/viral_database/kmer_count_k" +
     k_number+"_union_apriori_min_support"+str(sys.argv[2]), "w").writelines(rules)
