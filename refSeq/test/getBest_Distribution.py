import numpy as np
import scipy.stats as st
# 分布拟合
_data1 = {1: 237414592, 2: 24751998, 3: 7122753, 4: 3623556, 5: 1975215, 6: 1231667, 7: 719756, 8: 413919,
          9: 309696, 10: 255862, 11: 231553, 12: 179577, 13: 148377, 14: 134557, 15: 101203, 16: 109067, 17: 89980, 18: 84974, 19: 85988, 20: 121412, 21: 78282, 22: 44154, 23: 46891, 24: 26221, 25: 36870, 26: 13364, 27: 13063, 28: 13558, 29: 13753, 30: 10277, 31: 5118, 32: 5219, 33: 4678, 34: 4131, 35:
          3817, 36: 3801, 37: 3732, 38: 2879, 39: 2467, 40: 2339, 41: 1806, 42: 1706, 43: 1732, 44: 1620, 45: 1561, 46: 1476, 47: 1609, 48: 1386, 49: 1490, 50: 1663, 51: 1840, 52: 1088, 53: 1390, 54: 1158,
          55: 295, 56: 350, 57: 299, 58: 263, 59: 252, 60: 245, 61: 220, 62: 167, 63: 153, 64: 200, 65: 146, 66: 162, 67: 191, 68: 118, 69: 84, 70: 83, 71: 95, 72: 80, 73: 60, 74: 71, 75: 59, 76: 99, 77: 49, 78: 71, 79: 68, 80: 47, 81: 58, 82: 31, 83: 51, 84: 25, 85: 34, 86: 25, 87: 40, 88: 42, 89: 52,
          90: 29, 91: 9, 92: 22, 93: 35, 94: 23, 95: 12, 96: 8, 98: 2, 99: 7, 100: 1, 101: 4, 102: 1, 104: 1, 106: 2, 113: 2, 120: 1, 122: 2, 125: 1, 131: 1, 136: 2, 137: 1, 138: 1, 149: 1, 154: 1, 155: 1,
          158: 1, 159: 2, 162: 1}
_data1_density = {}
# 1: 237414592 只出现在1个物种里的k_mer有237414592个(占distinct k-mer 84%)
# Best fit reached using laplace, MLE value: 193810574.6741538
_data2 = {1: 236051083, 2: 25748257, 3: 7268837, 4: 3724741, 5: 2005247, 6: 1261319, 7: 729550, 8: 425599,
          9: 315367, 10: 260498, 11: 234369, 12: 183635, 13: 151081, 14: 137448, 15: 101770, 16: 110606, 17: 90347, 18: 85215, 19: 86795, 20: 122268, 21: 78543, 22: 44926, 23: 48550, 24: 27664, 25: 38369, 26: 14491, 27: 13832, 28: 14297, 29: 14309, 30: 10666, 31: 5341, 32: 5346, 33: 4759, 34: 4199, 35:
          3939, 36: 3855, 37: 3893, 38: 3037, 39: 2630, 40: 2376, 41: 1882, 42: 1794, 43: 1728, 44: 1647, 45: 1563, 46: 1534, 47: 1628, 48: 1418, 49: 1493, 50: 1679, 51: 1860, 52: 1102, 53: 1415, 54: 1222,
          55: 345, 56: 371, 57: 354, 58: 280, 59: 285, 60: 272, 61: 220, 62: 180, 63: 197, 64: 219, 65: 175, 66: 180, 67: 194, 68: 127, 69: 92, 70: 86, 71: 104, 72: 82, 73: 65, 74: 68, 75: 84, 76: 108, 77:
          61, 78: 96, 79: 83, 80: 56, 81: 75, 82: 39, 83: 64, 84: 35, 85: 41, 86: 30, 87: 41, 88: 47, 89: 49, 90: 37, 91: 13, 92: 25, 93: 28, 94: 38, 95: 30, 96: 10, 97: 4, 98: 4, 99: 9, 100: 4, 101: 9, 102: 1, 104: 6, 105: 4, 106: 2, 107: 1, 108: 2, 109: 1, 110: 2, 111: 6, 113: 2, 115: 3, 117: 1, 118:
          1, 120: 3, 121: 4, 122: 13, 123: 2, 124: 15, 125: 4, 127: 3, 129: 6, 130: 9, 131: 3, 133: 1, 134:
          1, 136: 11, 137: 6, 138: 4, 139: 7, 141: 2, 142: 1, 145: 1, 146: 6, 149: 1, 150: 1, 151: 1, 154: 1, 155: 1, 158: 3, 159: 2, 162: 7, 163: 6, 164: 5, 166: 1, 168: 7, 170: 8, 172: 1, 174: 10, 175: 1, 178: 3, 179: 1, 182: 4, 183: 2, 186: 2, 205: 1, 206: 1, 207: 1, 209: 1, 221: 1, 228: 1, 231: 1, 237: 1, 245: 1, 268: 1, 294: 1, 326: 1, 338: 1, 339: 1, 341: 1, 345: 1, 347: 2, 348: 1, 349: 1, 354: 1, 355: 1, 372: 3, 374: 1, 375: 2, 377: 2, 378: 5, 379: 2, 381: 3, 408: 1, 425: 1, 435: 1, 446:
          2, 447: 1, 450: 1, 451: 1, 459: 1, 489: 2, 519: 1, 734: 1, 741: 1, 872: 1, 882: 1, 889: 1, 903: 1, 909: 1, 913: 1, 932: 1, 934: 1, 2812: 1}  # 在整个病毒参考基因组中,只重复1次的k_mer有236051083个(占all k-mer 61.3%)
# Best fit reached using norm,so MLE value: 535859008.3514898
count = 0
for k, v in _data1:
    _data1[k]


def getBest_Distribution(_data):
    data = np.array(_data)
    distributions = [st.laplace, st.norm]
    mles = []
    for distribution in distributions:
        pars = distribution.fit(data)
        mle = distribution.nnlf(pars, data)
        mles.append(mle)
    best_fit = sorted(zip(distributions, mles), key=lambda d: d[1])[0]
    print('Best fit reached using {}, MLE value: {}'.format(
        best_fit[0].name, best_fit[1]))


# _data1_list = sum([[k]*v for k, v in _data1.items()], [])
_data2_list = sum([[k]*v for k, v in _data2.items()], [])
getBest_Distribution(_data2_list)
