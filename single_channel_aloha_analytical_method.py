# -*- coding: utf-8 -*-
__author__ = 'qsong'

import numpy as np
import csv

def cul_power_chf2cdf(N, P, alpha, l, m):
    '''
        本方法用于计算随机变量 总干扰功率强度值为N 的累计概率值，基本思路是通过对 特征函数的数值计算，计算反傅里叶变换系数
        @N:         总干扰功率值为N, a scalar number
        @P:         概率矢量, a matrix of shape (1, X)
        @alpha:     任意slot的包达到率平均值
        @l@M:       归一化后 接收功率的表示值
    '''

    # 返回值初始化为0
    result = 0.0
    #构建乘子
    # 概率矢量的维度 等于 最大传输次数 （最大重传次数 + 1）
    K = P.size
    # Be careful, in python, array type is not matrix. THe operators such as *  / are still element-wise
    LM = np.array([l**k*m**(K-1-k) for k in range(K)]).reshape(K, 1)
    # t 是特征函数的自变量，取值范围从0到PI
    t = np.linspace(0.0001, np.pi, 1000)
    E1 = np.exp(1j*LM*t)
    chf = np.exp(alpha*(P.dot(E1)-sum(P)))
    # 乘子 E2
    E2 = np.exp(-1j*t*N*0.5)
    Y = np.sin((N+1)*t*0.5) * np.real(chf*E2) / np.sin(0.5*t)
    result = np.trapz(Y, x=t)/np.pi
    return result


def solve_fxp(P, delta, alpha, l, m, threld):
    '''
        @P          : start probability vector
    '''
    # i is the iteration number
    i = 0
    MAX_I = 1000
    # The difference between two consecutive iterations. The defaut value is a large number.
    # When ecart is less than a value, the iteration will stop
    ecart = 100

    K = P.size
    M = [np.floor(l**k*m**(K-1-k)/threld) -1 for k in range(K-1)]

    tmp = np.copy(P)
    # print "For intensity:", alpha
    while i < MAX_I and ecart > delta:
        # print "In iteration ", i, "The proba vector P is:\t\t", ["{:10.5f}".format(x) for x in list(P)]
        # P[:-1] is a slicing operation in numpy, meaning taking elements except the last one.
        failure = [1-cul_power_chf2cdf(x, P, alpha, l, m) for x in M]
        # print "In iteration ", i, "New failure proba:\t\t\t", ["{:10.5f}".format(x) for x in failure]
        # print "In iteration ", i, "New success proba:\t\t\t", ["{:10.5f}".format(1-x) for x in failure]

        tmp[1:] = [p_f*P[x] for x, p_f in enumerate(failure)]
        ecart = sum(abs(tmp-P))


        P = np.copy(tmp)
        i += 1
    # Note that np.floor(l**(K-1)/threld-M[-1])
    p_loss = P[-1]*(1-cul_power_chf2cdf(np.floor(l**(K-1)/threld-l**(K-1)), P, alpha, l, m))
    # p_loss = (1-np.exp(-sum(P)*alpha) - sum(P)*alpha*np.exp(-sum(P)*alpha))**MAX_TRANS
    # p_loss = (1-np.exp(-sum(P)*alpha))**MAX_TRANS
    # print "alpha", alpha, "P[-1]", P[-1], "\t\tFailed Proba", 1-cul_power_chf2cdf(np.floor(l**(K-1)/threld), P, alpha, l, m), "proba_loss", p_loss
    # print alpha, [float("{0:.4g}".format(element)) for element in list(P)], p_loss
    return alpha, [float("{0:.4g}".format(element)) for element in list(P)], p_loss

def do_analytic(P, DELTA, start, end, l, m, thresld, step):
    result = []
    while start <= end:
        alpha, vector_p, p_loss = solve_fxp(P, DELTA, start, l, m, thresld)
        start += step
        vector_p.append(p_loss)
        vector_p.append(alpha)
        result.append(vector_p)
    return result

if __name__ == "__main__":
    """
        result 1.5 ['0.5524', '1', '0.9958', '0.9841', '0.9542', '0.8599', '1', '0.792', '0.5847', '0.3795', '0.1801']
    """

    np.set_printoptions(precision=4)
    MAX_TRANS = 5
    alpha_start = 0.1
    # l, m = 2, 1
    l, m = 1, 1
    P = np.zeros(MAX_TRANS)
    P[0] = 1
    DELTA = 0.0000001
    THRESLD = 0.5
    alpha_end =2.05
    step = 0.05
    print [cul_power_chf2cdf(n, P, alpha_start, l, m) for n in range(50)]
    # print solve_fxp(P, DELTA, alpha, l, m, THRESLD)

    # solve_fxp(P, DELTA, 1.5, l, m, 0.33)

    result = do_analytic(P, DELTA, alpha_start, alpha_end, l, m, THRESLD, step)


    result_f = "/Users/qsong/Documents/MATLAB/analytical_result_0p5.csv"
    with open(result_f, 'w') as f_handler:
        spamwriter = csv.writer(f_handler, delimiter=',')
        for n, vector_p in enumerate(result, 1):
            print n, vector_p
            spamwriter.writerow(vector_p)