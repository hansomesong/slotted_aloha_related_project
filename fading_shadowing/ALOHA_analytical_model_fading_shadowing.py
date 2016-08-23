# -*- coding: utf-8 -*-
__author__ = 'qsong'

import numpy as np
import csv

def trans_failure_p(thrld, P, alpha, v, k, sigma):
    '''
        本方法用于计算随机变量 总干扰功率强度值为N 的累计概率值，基本思路是通过对 特征函数的数值计算，计算反傅里叶变换系数
        本方法用于计算在考虑 Narrow-band fading and shadowing 的情况下，第 k 次重传的失败概率
        @thrld:     SINR 的门限值, a scalar number，SINR 大于等于该值的传输才算是成功的
        @P:         概率矢量, a matrix of shape (1, X)
        @alpha:     任意slot的包达到率平均值
        @:          power increment factor
    '''

    # 返回值初始化为0
    result = 0.0
    # 概率矢量的维度 等于 最大传输次数 （最大重传次数 + 1）
    K = P.size
    # Be careful, in python, array type is not matrix. THe operators such as *  / are still element-wise
    MULTIPLIER_1 = np.array(
        [
            np.power(
                1.0+np.power(
                    v**(m-k)*thrld,
                    np.power(
                        1.0+0.125*np.pi*sigma**2,
                        -0.5
                    )
                ),
                -1.0
            ) for m in range(K)
        ]
    ).reshape(K, 1) #相当于转置

    return 1-np.exp(alpha*(P.dot(MULTIPLIER_1)-sum(P)))

def solve_fxp(P, delta, alpha, v, sigma, thrld):
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
    # print "For intensity:", alpha
    while i < MAX_I and ecart > delta:
        # print "In iteration ", i, "The proba vector P is:\t\t", ["{:10.5f}".format(x) for x in list(P)]
        # P[:-1] is a slicing operation in numpy, meaning taking elements except the last one.
        failure = [trans_failure_p(thrld, P, alpha, v, k, sigma) for k in range(K-1)]
        # print "In iteration ", i, "New failure proba:\t\t\t", ["{:10.5f}".format(x) for x in failure]
        # print "In iteration ", i, "New success proba:\t\t\t", ["{:10.5f}".format(1-x) for x in failure]

        tmp[1:] = [p_f*P[x] for x, p_f in enumerate(failure)]
        ecart = sum(abs(tmp-P))


        P = np.copy(tmp)
        i += 1
    # Note that np.floor(l**(K-1)/threld-M[-1])
    result = []
    p_loss = P[-1]*trans_failure_p(thrld, P, alpha, v, K-1, sigma)
    result.extend([float("{0:.4g}".format(element)) for element in list(P)])
    result.extend(list(p_loss))
    result.append(alpha)

    return result

def do_analytic(P, delta, start, end, l, m, sigma, thrld, step):
    result = []
    while start <= end:
        vector_p = solve_fxp(P, delta, start, 1.0*l/m, sigma, thrld)
        start += step
        result.append(vector_p)
    return result

if __name__ == "__main__":

    MAX_TRANS = 5
    DELTA = 0.0000001
    SIGMA = np.sqrt(2)*np.log(10)/10.0
    P = np.zeros(MAX_TRANS)
    P[0] = 1
    THRLD = 0.0  #unit dB
    THRLD = 10**(THRLD/10.0) # decimal
    alpha_start = 0.2
    l, m = 1, 1
    alpha_end = 2.02
    step = 0.005

    result = do_analytic(P, DELTA, alpha_start, alpha_end, l, m, SIGMA, THRLD, step)

    result_f = "fading_shadowing_analytical_result_threshold={0}dB_l={1}_m={2}_sigma=1.csv".format(10*np.log10(THRLD), l, m)
    with open(result_f, 'w') as f_handler:
        spamwriter = csv.writer(f_handler, delimiter=',')
        for n, vector_p in enumerate(result, 1):
            print n, vector_p
            spamwriter.writerow(vector_p)


