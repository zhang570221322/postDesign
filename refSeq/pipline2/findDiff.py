# %%
f1 = open("grinder-reads.fa.result.handle", 'r')
f2 = open("grinder-reads.fa.result.handle_2", 'r')
while True:
    line1_1, line1_2, line1_3 = f1.readline(), f1.readline(), f1.readline()
    line2_1, line2_2, line2_3 = f2.readline(), f2.readline(), f2.readline()

    if not (line1_1 and line1_1):
        break
    seq_id = None

    seq_id = line1_1.split('\t')[0]
    taxid_weight1 = line1_2[:-1].split("\t")
    taxid_weight2 = line2_2[:-1].split("\t")
    flag = False
    try:
        for i in range(len(taxid_weight1)):
            taxid1, weight1 = taxid_weight1[i].split(':')
            taxid2, weight2 = taxid_weight2[i].split(':')
            if taxid1 != taxid2:
                flag = True
                break
            if round(float(weight1), 3) != float(weight2):
                flag = True
                break
        if flag:
            print(seq_id)
        flag = False
    except Exception:
        pass
