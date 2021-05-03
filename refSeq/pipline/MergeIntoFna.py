import os
import os.path
import re
import gzip
import sys

basename = sys.argv[1]
assembly_summary = sys.argv[2]
out_file_name = sys.argv[3]


def read_gz_file(path):
    if os.path.exists(path):
        with gzip.open(path, 'r') as pf:
            for line in pf:
                yield line
    else:
        print('the path [{}] is not exist!'.format(path))


def searchFile(pathname, filename):
    matchedFile = []
    for root, dirs, files in os.walk(pathname):
        for file in files:
            if re.match(filename, file):
                fname = os.path.abspath(os.path.join(root, file))
                # print(fname)
                matchedFile.append(fname)
    return matchedFile


_GCF_TAXID = {}
for line in open(assembly_summary, "rb").readlines()[2:]:
    line_list = str(line, "utf-8").split("\t")
    _GCF_TAXID[line_list[0]] = line_list[5]
# _GCF_TAXID1 = {}
# for line in open("E:\\Github\\classification\\refSeq\\pipline\\database\\seq2taxid2.db", "r").readlines():
#     line_list = line.split(",")
#     _GCF_TAXID1[line_list[2]] = line_list[0]

wr_f = open(out_file_name, "w")
list_s = searchFile(basename, r'.+\.gz')
for fasta in list_s:
    taxid = _GCF_TAXID.get(os.path.dirname(fasta).split("/")[-1])
    if taxid:
        for line in read_gz_file(fasta):
            str_line = str(line, "utf-8")
            if str_line.startswith(">"):
                # each seq metaInfo
                str_line = str_line[0]+str(taxid)+"|"+str_line[1:]
            wr_f.writelines(str_line)
    else:
        print("WARNING :", fasta)
