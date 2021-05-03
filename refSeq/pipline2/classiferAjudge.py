import pickle
from read_kmer import read
from rule import weight
from read_kmer.Taxonomy import Tax
# 分类树信息
all_taxid = "/home/wlzhang/classfication/refseq/viral_database_kraken2/Script/all_taxId.db"
# 10239为病毒根节点.
tax = Tax(10239, all_taxid)

_kmer_taxid = {}
with open("/home/wlzhang/classfication/refseq/viral_database_kraken2/_TOTAL_kmer_taxid_count_bin", 'rb') as f:
    _kmer_taxid = pickle.load(f)
# 测试数据路径
test_data_dir = "/home/wlzhang/classfication/refseq/viral_database_kraken2/test_data/"
# 测试数据文件名
test_data_name = 'grinder-reads.fa'
result_csv = open(test_data_dir+test_data_name+".result.handle", 'w')
kmer_size = 28
for title, seq in read.fasta_generator(test_data_dir+test_data_name):
    # 格式化title
    # 最终信息为二维表
    # 表头为 [ID, taxid, NC_NO, position, description ]
    format_title = read.format_title(title)
    # 输出表头
    result_csv.write(">"+"\t".join(format_title))
    # 更改152835..152984为152835的格式
    format_title[3] = int(format_title[3].split(".")[0].split("(")[-1])
    kmer_hit_informations = []
    # 处理一个read
    # ATCGGGAGGC 3'-> 5'
    # TAGCCCTCCG 另一条链
    # CGGAGGGCTA 5'-> 3'
    # GCCTCCCGAT 另一条链
    for kmer in read.kmer_generator(seq, kmer_size):
        raw_kmer = kmer
        # 对 5'->3' 做检测
        OUT = _kmer_taxid.get(kmer)
        if not OUT:
            # 获取另一条链
            complementary_kmer = read.getComplement(kmer)
            OUT = _kmer_taxid.get(complementary_kmer)
        # 对 3'->5'做检测
        if not OUT:
            kmer = kmer[::-1]
            OUT = _kmer_taxid.get(kmer)
            if not OUT:
                complementary_kmer = read.getComplement(kmer)
                OUT = _kmer_taxid.get(complementary_kmer)
        if not OUT:
            OUT = "0"
        # 存储命中信息
        kmer_hit_informations.append([format_title[3], raw_kmer, OUT])
        # 更新碱基坐标信息.
        format_title[3] = format_title[3]+kmer_size
    # predi: 预测taxid
    # weigh_hit: 命中位置的加权每个taxid的权重
    # non_hit: 没有命中的坐标位置
    predi, weigh_hit, non_hit = weight.judge1(kmer_hit_informations)
    # 如果为多个判断
    # 应该降低LCA的权重
    if len(predi) > 1:
        # 得到非叶子节点的taxid
        Leaf_taxid = []
        LCA_taxid = []
        for taxid in predi:
            if (tax.tree & taxid).is_leaf():
                Leaf_taxid.append(taxid)
            else:
                LCA_taxid.append(taxid)
        # 得到对应的kmers
        LCA_kmers = [kmer for pos, kmer,

                     out in kmer_hit_informations if out != '0' and out.split(",")[0] in LCA_taxid]
        # 降低权重,罚分为0.1
        for kmer in LCA_kmers:
            taxid, kmer_weight = _kmer_taxid[kmer].split(',')
            _kmer_taxid[kmer] = taxid+','+str(int(kmer_weight)+1)
        # 如果预测的叶子节点只有一个
        if len(Leaf_taxid) == 1:
            predi = Leaf_taxid
        else:
            # 如果预测的叶子节点有多个, 即A物种一半,B物种一半.
            # 可以采用kraken2的方式,root-to-leaf的方式来进行甄别.
            print("kmer命中有多个叶子节点,且权重一样{seq_id}".format(seq_id=format_title[0]))
    result_csv.write("\t{0}\n".format(
        ",".join(predi)))
    result_csv.write("\t".join(weigh_hit)+"\n")
    result_csv.write("\t".join(map(str, non_hit))+"\n")
result_csv.flush()
result_csv.close()
