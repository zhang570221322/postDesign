from typing import Iterable, List, Union
from Bio.SeqIO.FastaIO import SimpleFastaParser


def getComplement(rev: str):
    """输入一个kmer,返回它的逆序
    """
    rev = rev.replace("A", "X")
    rev = rev.replace("T", "A")
    rev = rev.replace("X", "T")
    rev = rev.replace("C", "X")
    rev = rev.replace("G", "C")
    rev = rev.replace("X", "G")
    return rev


def kmer_generator(dna, kmer_size) -> str:
    """根据核苷酸序列分为kemr
    kmer_size为k的大小
    """
    for start in range(0, len(dna)-(kmer_size-1), 1):
        temp = dna[start:start+kmer_size]
        yield temp


def fasta_generator(file_name: str) -> Iterable[Union[str, str]]:
    """根据fasta文件获取到信息
    返回(tile,seq)类型的元组的生成器
    """
    with open(file_name) as in_handle:
        for title, seq in SimpleFastaParser(in_handle):
            yield (title, seq)


def format_title(title: str) -> List[str]:
    description = title[title.index('description')+13:-1]
    split_title = title.split(" ")
    ID = split_title[0]
    taxid, NC_NO = split_title[1].split("|")[1:]
    position = split_title[2].split("=")[-1]
    return [ID, taxid, NC_NO, position, description]
