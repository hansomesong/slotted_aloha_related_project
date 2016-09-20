# -*- coding: utf-8 -*-
__author__ = 'qsong'

import numpy as np
import csv
import scipy.special
from itertools import combinations
from operator import mul

def cal_prob_n_occur(prob_vector, n):
    tmp_p = list(list(e) for e in combinations(prob_vector, n))
    rest_p = [list(set(prob_vector) - set(element)) for element in tmp_p]
    new_p = [[1-x for x in e[0]]+e[1] for e in zip(tmp_p, rest_p)]
    # print new_p
    return sum([reduce(mul, element) for element in new_p])

def trans_failure_p(thrld, P, alpha, v, k, rho, sigma):
    '''
        本方法用于计算随机变量 总干扰功率强度值为N 的累计概率值，基本思路是通过对 特征函数的数值计算，计算反傅里叶变换系数
        本方法用于计算在考虑 Narrow-band fading and shadowing 的情况下，第 k 次重传的失败概率
        @thrld:     SINR 的门限值, a scalar number，SINR 大于等于该值的传输才算是成功的
        @P:         概率矢量, a matrix of shape (1, X)
        @alpha:     任意slot的包达到率平均值
        @:          power increment factor
        @k:         第 k th 重传
    '''
    # 概率矢量的维度 等于 最大传输次数 （最大重传次数 + 1）
    K = P.size
    BETA = np.log(10)/10
    p_zero_inf = np.exp(-1.0*alpha*sum(P))

    # mean_Y_plus = 0.0
    # mean_Y_2_plus = 0.0
    # for n in range(1, K+1, 1):
    #     tmp_proba = cal_prob_n_occur(P, n)
    #     mean_Y_plus_n = np.exp(np.power(BETA*sigma, 2))*sum([np.power(v, m-k)*alpha*P[m] for m in range(n)])
    #     mean_Y_plus += mean_Y_plus_n*tmp_proba
    #     var_Y_plus_n = np.exp(np.power(2*BETA*sigma, 2))*sum([np.power(v, 2*(m-k))*alpha*P[m] for m in range(n)])
    #     mean_Y_2_plus_n = var_Y_plus_n + mean_Y_plus_n**2
    #     mean_Y_2_plus += mean_Y_2_plus_n*tmp_proba


    # mean_Y = np.exp(np.power(BETA*sigma, 2))*sum([np.power(v, m-k)*alpha*P[m] for m in range(K)])
    mean_Y = (1+BETA**2*sigma**2)*sum([np.power(v, m-k)*alpha*P[m] for m in range(K)])
    # mean_Y_2 = (7*BETA**4*sigma**4+1+3*BETA**2*sigma**2)*sum([np.power(v, 2*(m-k))*alpha*P[m] for m in range(K)]) + mean_Y**2
    var_Y = np.exp(np.power(2*BETA*sigma, 2))*sum([np.power(v, 2*(m-k))*alpha*P[m] for m in range(K)])
    mean_Y_2 = np.power(mean_Y, 2) + var_Y
    mean_Y_plus = mean_Y/(1-p_zero_inf)
    mean_Y_2_plus = mean_Y_2/(1-p_zero_inf)


    mu_eq = 2*np.log(mean_Y_plus) - 0.5*np.log(mean_Y_2_plus)
    sigma_eq = np.log(mean_Y_2_plus) - 2*np.log(mean_Y_plus)


    # 注意：这里的 thrld 一定是 以 分贝为单位的，千万不要弄错了。。。
    ERROR = (-BETA*thrld - mu_eq)/(np.sqrt(2*sigma_eq**2))

    print "k", k, "mu_eq", mu_eq, "sigma_eq", sigma_eq, "Failure prob", 0.5*(1-scipy.special.erf(ERROR))*(1-p_zero_inf), "Input P", P

    return 0.5*(1-scipy.special.erf(ERROR))*(1-p_zero_inf)


def solve_fxp(P, delta, alpha, v, sigma, rho, thrld):
    '''
        @P          : start probability vector
        @threld     : the SINR threshold, in unit of dB
    '''
    # i is the iteration number
    i = 0
    MAX_I = 1000
    # The difference between two consecutive iterations. The defaut value is a large number.
    # When ecart is less than a value, the iteration will stop
    ecart = 100
    K = P.size
    tmp = np.copy(P)

    # print "Current power increment factor is:", v
    # print "For intensity:", alpha
    while i < MAX_I and ecart > delta:
        # print "In iteration ", i, "The proba vector P is:\t\t", ["{:10.5f}".format(x) for x in list(P)]
        # P[:-1] is a slicing operation in numpy, meaning taking elements except the last one.
        failure = [trans_failure_p(thrld, P, alpha, v, k, rho, sigma) for k in range(K-1)]
        # print "In iteration ", i, "New failure proba:\t\t\t", ["{:10.5f}".format(x) for x in failure]
        # print "In iteration ", i, "New success proba:\t\t\t", ["{:10.5f}".format(1-x) for x in failure]

        tmp[1:] = [p_f*P[x] for x, p_f in enumerate(failure)]
        ecart = sum(abs(tmp-P))
        P = np.copy(tmp)
        i += 1
    # Note that np.floor(l**(K-1)/threld-M[-1])
    result = []
    p_loss = P[-1]*trans_failure_p(thrld, P, alpha, v, K-1, rho, sigma)
    result.extend([float("{0:.4g}".format(element)) for element in list(P)])
    result.append(p_loss)
    result.append(alpha)

    return result

def do_analytic(P, delta, start, end, l, m, sigma, rho, thrld, step):
    result = []
    while start <= end:
        vector_p = solve_fxp(P, delta, start, 1.0*l/m, sigma, rho, thrld)
        start += step
        result.append(vector_p)
    return result

if __name__ == "__main__":

    MAX_TRANS = 5
    DELTA = 0.0000001
    SIGMA = 1
    P = np.zeros(MAX_TRANS)
    P[0] = 1
    # 注意： 这里 门限值一定是 分贝单位
    THRLD = -3
    alpha_start = 1.12
    l, m = 1, 1
    alpha_end = 1.12
    step = 0.005
    rho = 1

    result = do_analytic(P, DELTA, alpha_start, alpha_end, l, m, SIGMA, rho, THRLD, step)

    print result

    # result_f = "shadowing_analytical_result_threshold={0}dB_l={1}_m={2}_sigma={3}.csv".format(THRLD, l, m, SIGMA)
    #
    # with open(result_f, 'w') as f_handler:
    #     spamwriter = csv.writer(f_handler, delimiter=',')
    #     for n, vector_p in enumerate(result, 1):
    #         print n, vector_p
    #         spamwriter.writerow(vector_p)


