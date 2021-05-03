from read_kmer.Taxonomy import Tax
# all_taxid = "/home/wlzhang/classfication/refseq/viral_database_kraken2/test_data/all_taxId.db"
all_taxid = "/home/cszhang/workSpace/classification/refSeq/pipline2/all_taxId.db"
tax = Tax(10239, all_taxid)


def my_print(a, b):
    print(a, round(b*100, 6))


# result = "/home/wlzhang/classfication/refseq/viral_database_kraken2/test_data/grinder-reads.fa.result.handle"
result = "/home/cszhang/workSpace/classification/refSeq/pipline2/grinder-reads.fa.result.handle_2"
total = 0
error_pricision = 0
pricision = 0
kmer_found_pricision = 0
mutiple_pricison = 0
empty_pricision = 0
LCA_pricision = {}
# 存储ID信息
error_pricision_id = []
pricision_id = []
kmer_found_pricision_id = []
mutiple_pricison_id = []
empty_pricision_id = []
LCA_pricision_id = {}
with open(result, "r") as f:
    line = f.readline()
    while line:
        # 进度
        if total % 10000 == 0:
            print("process:", total)
        # 去除>和\n
        title = line[1:-1]
        out_list = f.readline()[:-1]
        non_set = f.readline()[:-1]
        line = f.readline()
        total += 1
        title_split = title.split('\t')
        # 如果预测只有一个.
        if len(title_split[-1].split(',')) == 1:
            # 如果一个也没命中
            if int(title_split[-1]) == 0:
                empty_pricision += 1
                empty_pricision_id.append(title_split[0])
                continue
                # 如果预测命中
            if int(title_split[1]) == int(title_split[-1]):
                pricision += 1
                pricision_id.append(title_split[0])
            # 如果真实taxid存在kmer的命中列表
            # TODO需要降低干扰命中的权重使其选择到.
            elif title_split[1] in out_list:
                kmer_found_pricision += 1
                kmer_found_pricision_id.append(title_split[0])
            # 如果不在里面, 看是否是上升到LCA了.
            else:
                # 与数据库保持最新
                real = str(tax.merge_dic.get(
                    int(title_split[1]), title_split[1]))
                predi = str(tax.merge_dic.get(
                    int(title_split[-1]), title_split[-1]))
                # 查询共同祖先
                common_node = (tax.tree & predi).get_common_ancestor(
                    [real, predi])
                # 去LCA之后为节点本身,说明上升到LCA了.
                if common_node.name == predi:
                    LCA_pricision.setdefault(common_node.rank, 0)
                    LCA_pricision[common_node.rank] += 1
                    temp = LCA_pricision_id.setdefault(common_node.rank, [])
                    temp.append(title_split[0])
                else:
                    error_pricision += 1
                    error_pricision_id.append(title_split[0])
        elif title_split[1] in title_split[-1]:
            mutiple_pricison += 1
            mutiple_pricison_id.append(title_split[0])
my_print("直接命中精度", pricision/total)
my_print("一个也未命中的精度", empty_pricision/total)
my_print("kmer命中精度", kmer_found_pricision/total)
my_print("多个预测值", mutiple_pricison/total)
my_print("错误:", error_pricision/total)
for k, v in LCA_pricision.items():
    my_print(k+"水平精度", v/total)


# 2.
# 直接命中精度 69.516
# 一个也未命中的精度 0.326
# kmer命中精度 6.452
# 多个预测值 0.0
# 错误: 0.0
# order水平精度 0.52
# no rank水平精度 3.392
# genus水平精度 16.538
# family水平精度 1.686
# subfamily水平精度 0.553
# species水平精度 0.987
# clade水平精度 0.027
# superkingdom水平精度 0.003
