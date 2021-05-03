# 加权.
from typing import List, Union


def final_out(ls: List):
    """策略组合
    """
    _dic = strategy_weight(ls)
    return _dic, get_max_dic(_dic)


def strategy_weight(ls: List):
    """根据
    [[66797,1],[66797,1],[66797,3]]
    这样的命中信息得到count,且对其做加权,如66797,3,则权值为0.333
    最终结果为: {'66797':'2.333'}
    """
    _dic = {}
    for out, fre in ls:
        # 获取到了值,如果值不存在设置为fre
        temp = _dic.get(out)
        fre = float(fre)
        if temp:
            _dic[out] += round(1/fre, 3)
        else:
            _dic[out] = round(1/fre, 3)
    return _dic


def get_max_dic(d: dict):
    """获取字典值最大的key
    如果不存在,则返回0
    """
    if not d:
        return ['0']
    max_value = max(d.values())
    return [key for key, value in d.items() if value == max_value]


def judge1(kmer_hit_information: List) -> Union[List, List, List]:
    """根据 kmer_hit_information[['碱基位置','taxid,次数'],...]来判断read的最终标签.
    Output:
        1. out: 预测taxid,多个以,号分割
        2. 命中位置的加权每个taxid的权重
        3. 没有命中的坐标位置
    """
    non_set = []
    out_list = []

    for kmer_position, kmer, temp in kmer_hit_information:
        if temp.strip() != '0':
            out_list.append(temp.split(','))
        else:
            non_set.append(kmer_position)
    _dic, out = final_out(out_list)
    return out, [k+':'+str(round(v, 3)) for k, v in _dic.items()], non_set
