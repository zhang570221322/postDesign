import pickle
from read_kmer import read
# library是kraken2下载的.
# 把数据库加载进内存.
# _library_fna28_union_bin数据形式是
# ATACTAAAGGGGGCCC  10699:1 20978:3    10777:4  物种内重复和物种间重复.
# _kmer_taxid = {}
# with open("/home/wlzhang/classfication/refseq/viral_database_kraken2/_library_fna28_union_bin", 'rb') as f:
#     _kmer_taxid = pickle.load(f)
# _TOTAL_kmer_taxid_count_bin的数据格式是
# 'AGGATTAATGCACTATCTCAAGCGTTTG':'10662,6' (后面是出现的LCA,如果是叶子节点则是物种本身,然后是出现次数) #TODO 但是这里的次数是物种间次数和物种内次数的总和.
_kmer_taxid = {}
with open("/home/wlzhang/classfication/refseq/viral_database_kraken2/_TOTAL_kmer_taxid_count_bin", 'rb') as f:
    _kmer_taxid = pickle.load(f)
# 测试数据路径
test_data_dir = "/home/wlzhang/classfication/refseq/viral_database_kraken2/test_data/"
# 测试数据文件名
test_data_name = 'grinder-reads.fa'
result_csv = open(test_data_dir+test_data_name+".result", 'w')

kmer_size = 28
for title, seq in read.fasta_generator(test_data_dir+test_data_name):
    # 格式化title
    # 最终信息为二维表
    # 表头为 [ID, taxid, NC_NO, position, description ]
    format_title = read.format_title(title)
    # 输出表头
    result_csv.write(">"+"\t".join(format_title)+"\n")
    # 更改152835..152984为152835的格式
    format_title[3] = int(format_title[3].split(".")[0])
    for kmer in read.kmer_generator(seq, kmer_size):
        # 对正链做检测
        OUT = _kmer_taxid.get(kmer)
        if not OUT:
            # 获取负链
            inverse_kmer = read.getComplement(kmer)
            OUT = _kmer_taxid.get(inverse_kmer)
        # 检测命中
        if not OUT:
            OUT = "0"
        # 写入信息
        result_csv.write("{start_index}\t{out}\n".format(
            start_index=format_title[3],
            out=OUT
        ))
        # 更新碱基坐标信息.
        format_title[3] = format_title[3]+kmer_size
result_csv.flush()
result_csv.close()
