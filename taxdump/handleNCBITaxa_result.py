import numpy as np
_Total_0 = {}
_Total_1 = {}
file_union = "/home/wlzhang/classfication/refseq/viral_database/kmer_count_k28_union1_kmer_self_distance"
with open(file_union, "r") as f:
    temp = f.readlines()
    _Total_0 = eval(temp[0].replace("\n", ""))
    _Total_1 = eval(temp[1].replace("\n", ""))
a = np.array([[k, np.array(v).mean()] for k, v in _Total_0.items()])

b = np.array([[k, np.array(v).mean()] for k, v in _Total_1.items()])
a = np.around(np.sort(a, 0), 3)
b = np.around(np.sort(b, 0), 3)
for i in range(len(a)):
    print("{product: '%s', 'L(x)': %s ,'Cost(x,x)': %s }," % (
        str(a[i, 0]), str(a[i, 1]), str(b[i, 1])))
