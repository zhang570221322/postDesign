
"""统计kraken2的结果
"""
from read_kmer.Taxonomy import Tax
all_taxid = "/home/cszhang/workSpace/classification/refSeq/pipline2/all_taxId.db"
tax = Tax(10239, all_taxid)


def my_print(a, b):
    print(a, round(b*100, 6))


#  原始数据
my_result = "/home/cszhang/workSpace/classification/refSeq/pipline2/grinder-reads.fa"
# 正确的序列id和对应的真是taxid.
# 用来和kraken2检测
data_result = {}
for line in open(my_result, 'r'):
    if line.startswith(">"):
        title = line
        split_title = title.split(" ")
        ID = split_title[0][1:]
        taxid = split_title[1].split("|")[1]
        data_result[ID] = taxid
print("start")
result = "/home/cszhang/workSpace/classification/refSeq/pipline2/grinder-reads.fa.output"
total = 0
pricision = 0
kmer_found_pricision = 0
mutiple_pricison = 0
LCA_pricision = {}
empty_pricision = 0
error_pricision = 0
error_pricision_id = []
with open(result, "r") as f:
    line = f.readline()
    while line:
        # 进度
        if total % 1000 == 0:
            print("process:", total)
        line_split = line[:-1].split('\t')
        is_classification = line_split[0] == 'C'
        seq_id = line_split[1]
        predi = line_split[2]
        out_list = line_split[4:]
        line = f.readline()
        total += 1
        # 如果一个也没命中
        if not is_classification:
            empty_pricision += 1
            continue
        # 如果预测命中
        if int(predi) == int(data_result[seq_id]):
            pricision += 1
        # 如果真实taxid存在kmer的命中列表
        # 需要降低干扰命中的权重使其选择到.
        elif data_result[seq_id] in out_list:
            kmer_found_pricision += 1
        # 如果不在里面, 看是否是上升到LCA了.
        else:
            # 与数据库保持最新
            real = str(tax.merge_dic.get(
                int(data_result[seq_id]), data_result[seq_id]))
            predi = str(tax.merge_dic.get(
                int(predi), predi))
            # 查询共同祖先
            common_node = ''
            # 去LCA之后为节点本身,说明上升到LCA了.
            try:
                common_node = (tax.tree & predi).get_common_ancestor(
                    [real, predi])
            except Exception as e:
                common_node = tax.tree.get_common_ancestor([real, predi])
            if common_node:
                if common_node.name == predi:
                    temp = LCA_pricision.setdefault(common_node.rank, 0)
                    LCA_pricision[common_node.rank] += 1
                else:
                    error_pricision += 1
                    error_pricision_id.append(seq_id)
            else:
                print([seq_id, real, predi])

my_print("直接命中精度", pricision/total)
my_print("一个也未命中的精度", empty_pricision/total)
my_print("kmer命中精度", kmer_found_pricision/total)
my_print("多个预测值", mutiple_pricison/total)
my_print("错误:", error_pricision/total)
for k, v in LCA_pricision.items():
    my_print(k+"水平精度", round(v/total, 5))
# 直接命中精度 70.854
# 一个也未命中的精度 0.0
# kmer命中精度 0.0
# 多个预测值 0.0
# 错误: 0.034
# no rank水平精度 6.778
# genus水平精度 19.878
# subfamily水平精度 1.608
# family水平精度 0.79
# order水平精度 0.058

# 输出错误

error_set = set()
for i in [data_result[i] for i in error_pricision_id]:
    error_set.add(i)
    # 找到存在于grinder-reads.fa.output的相关seq_id的taxid
for line in open('grinder-reads.fa.output', 'r'):
    line_split = line.split('\t')
    seq_id = line_split[1]
    predi = line_split[2]
    if seq_id in error_pricision_id:
        OUTs = set([i.split(':')[0] for i in line_split[4].strip().split(' ')])
        a = set(OUTs)
        try:
            a.remove('0')
        except Exception:
            pass
        if data_result[seq_id] in a:
            print(seq_id)
        print("seq_id:{0},真实taxid:{1},预测tax_id{2}".format(
            seq_id, data_result[seq_id], predi))
        print('|'.join([data_result[seq_id]]+list(a)))
        print([i for i in line_split[4].strip().split(' ')])


with open('error_id', 'w') as test:
    for i in error_set:
        test.write(i+'\n')
