# -*- coding: utf-8 -*-
__author__ = 'qsong'

import numpy as np
import csv
import scipy.special


def trans_failure_p(thrld, P, alpha, v, k, rho, sigma):
    '''
        本方法用于计算随机变量 总干扰功率强度值为N 的累计概率值，基本思路是通过对 特征函数的数值计算，计算反傅里叶变换系数
        本方法用于计算在考虑 Narrow-band fading and shadowing 的情况下，第 k 次重传的失败概率
        @thrld:     SINR 的门限值, a scalar number，SINR 大于等于该值的传输才算是成功的
        @P:         概率矢量, a matrix of shape (1, X)
        @alpha:     任意slot的包达到率平均值
        @:          power increment factor
    '''
    # 概率矢量的维度 等于 最大传输次数 （最大重传次数 + 1）
    K = P.size
    BETA = np.log(10)/10
    p_zero_inf = np.exp(-1.0*alpha*sum(P))

    mean_Y = np.exp(0.5*(BETA*sigma)**2)*sum([(v**k)*rho*alpha*P[k] for k in range(K)])
    var_Y = np.exp(2*(BETA*sigma)**2)*sum([v**(2*k)*(rho**2)*alpha*P[k] for k in range(K)])
    # mean_Y2 = (rho**2)*np.exp((BETA*sigma)**2)*(
    #     sum([v**(2*k)*alpha*P[k] for k in range(K)]) + np.power(sum([(v**k)*alpha*P[k] for k in range(K)]), 2)
    # )
    mean_Y_plus = mean_Y/(1-p_zero_inf)
    mean_Y_2_plus = (mean_Y**2 + var_Y) / (1-p_zero_inf)

    mu_eq = 2*np.log(mean_Y_plus) - 0.5*np.log(mean_Y_2_plus)
    sigma_eq = np.log(mean_Y_2_plus) - 2*np.log(mean_Y_plus)

    # 注意：这里的 thrld 一定是 以 分贝为单位的，千万不要弄错了。。。
    ERROR = (np.log(v**k*rho) - BETA*thrld - mu_eq)/(np.sqrt(2*((BETA*sigma)**2+sigma_eq**2)))
    return 0.5*(1-scipy.special.erf(ERROR))*(1-p_zero_inf)
    # return 0.5*(1-scipy.special.erf(ERROR))


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
    SIGMA = 4.0
    P = np.zeros(MAX_TRANS)
    P[0] = 1
    # 注意： 这里 门限值一定是 分贝单位
    THRLD = 3
    alpha_start = 0.2
    l, m = 1, 1
    alpha_end = 1.02
    step = 0.005
    rho = 1

    result = do_analytic(P, DELTA, alpha_start, alpha_end, l, m, SIGMA, rho, THRLD, step)

    result_f = "shadowing_analytical_result_threshold=3dB_l=1_m=1_sigma=1.csv"
    with open(result_f, 'w') as f_handler:
        spamwriter = csv.writer(f_handler, delimiter=',')
        for n, vector_p in enumerate(result, 1):
            print n, vector_p
            spamwriter.writerow(vector_p)


