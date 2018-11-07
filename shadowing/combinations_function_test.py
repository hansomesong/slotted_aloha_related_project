__author__ = 'qsong'

from itertools import combinations
import numpy as np
from operator import mul
from scipy.special import lambertw

def cal_prob_n_occur(prob_vector, n):
    tmp_p = list(list(e) for e in combinations(prob_vector, n))
    rest_p = [list(set(prob_vector) - set(element)) for element in tmp_p]
    new_p = [[1-x for x in e[0]]+e[1] for e in zip(tmp_p, rest_p)]
    print new_p
    tmp_p = [[1-e for e in element] for element in tmp_p]
    # rest_p = []
    return sum([reduce(mul, element) for element in new_p])
    # print





if __name__ == "__main__":
    alpha = 0.9
    test_p = [1.0, 0.8416, 0.7083, 0.5961, 0.5017, 0.42226017342572847]
    # Each element is the prob for I=0
    test_p = [np.exp(-alpha*element) for element in test_p]
    print "input", test_p

    total = 0.0

    for i in range(len(test_p)+1):
        tmp  = cal_prob_n_occur(test_p, i)
        total += tmp
        print tmp

    print total

