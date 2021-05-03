from ete3 import NCBITaxa
import os


class Tax:
    ncbi = NCBITaxa()
    tree = None
    merge_dic = None

    def __init__(self, root_number, all_taxid):

        if not os.path.exists(all_taxid):
            raise FileNotFoundError("{0} is not found!".format(all_taxid))
        Tax.tree = Tax.ncbi.get_topology(Tax.ncbi.get_descendant_taxa(
            root_number), intermediate_nodes=True)
        print("load taxonomy database done")
        all_taxid = [str(taxid).replace("\n", "")
                     for taxid in open(all_taxid, "r").readlines()]
        # NCBI树更新了, 但是本地没有更新,因此需要找到新的与之对应
        Tax.merge_dic = Tax.ncbi._translate_merged(all_taxid)[1]
        print("load merge  database done")
