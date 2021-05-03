## Download raw fna and assembly_summary.txt
```bash
# assembly_summary.txt
wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt
wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/{category}/assembly_summary.txt
# rasyc download raw seq
# rsync_download.py
```
##  Merge into one file $1(dir including *.fna.gz) $2(assembly_summary_refseq.txt dir) $3 (out file's name,*.fna)
```bash
python MergeIntoFna.py E:\\Github\\classification\\refSeq\\pipline\\raw_data  E:\\Github\\classification\\refSeq\\stimu\\assembly_summary.txt 100_seq.fna
```
## 
    1. DNA2k-mer with taixd
    2. k-mer count  k-mer\ttaxid:count\ttaxid:count
    3. pick k-mer
    4.,利用k-mer count来聚类, 利用Taxonomy Group来递归
        4.1  利用DNA2k-mer with taixd组成taixd:"line_start:line_end"字典 , 可以知道哪些k-mer是哪个taxid DONE. 从这些k-mer pick
        4.2  一个k-mer的共同祖先表示是那个祖先的k-mer,然后对其挑选
        4.3 制定挑选规则
        4.4 机器学习进行二分类。