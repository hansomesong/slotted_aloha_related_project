# -*- coding: utf-8 -*-
__author__ = 'qsong'

import numpy as np
import csv

def chf_compound_sln(alpha, P, omega, l, m, k):
    '''
        Y_k = \sum_{m=0}^{K} v^{m-k} \sum_{j=1}^{N_m} e^{\beta\epsilion}
        计算 Y_k 的 CDF
        @alpha      : the fresh packet arrival intensity, unit per second
        @P          : the probability vector, initial value is e.g. [1,0,0,0,0]
        @sigma      : power control error variance, usually between 1dB and 4dB
        @v          : power increment factor, it may take 1, 2, 0,5
        @k          : the k th retransmission, varying from 0 to MAX_K
        @omega      : the values in frequency space, type np.array
    '''

    DIMENSION = P.size
    W_MULTIPLIER = np.array([np.exp(1j*omega*l**x*m**(DIMENSION-1-x)) for x in range(DIMENSION)])
    return np.exp(alpha*(P.dot(W_MULTIPLIER)-sum(P)))

def trans_failure_p(thrld, P, alpha, l, m, k):
    thrld = 10**(thrld/10.0)
    omega_end = 50
    N = 1000
    omega = np.linspace(0, omega_end, N)
    eta = 0.8
    EPS = 0.0001
    shift_omega = 1j*eta + omega + EPS
    Y = ((1+EPS*1j)/(eta-1j*omega))*np.exp(-1j*omega/thrld)*chf_compound_sln(alpha, P, shift_omega, l, m, k)

    result = np.exp(eta*1.0/thrld)*np.real(np.trapz(y=Y, x=omega))/np.pi
    # print "proba:", 1-result

    return 1-result

def solve_fxp(P, delta, alpha, l, m, thrld):
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
        failure = [trans_failure_p(thrld, P, alpha, l, m, k) for k in range(K-1)]
        # print "In iteration ", i, "New failure proba:\t\t\t", ["{:10.5f}".format(x) for x in failure]
        # print "In iteration ", i, "New success proba:\t\t\t", ["{:10.5f}".format(1-x) for x in failure]

        tmp[1:] = [p_f*P[x] for x, p_f in enumerate(failure)]
        ecart = sum(abs(tmp-P))
        P = np.copy(tmp)
        i += 1
    # Note that np.floor(l**(K-1)/threld-M[-1])
    result = []
    p_loss = P[-1]*trans_failure_p(thrld, P, alpha, l, m, K-1)
    result.extend([float("{0:.4g}".format(element)) for element in list(P)])
    result.append(p_loss)
    result.append(alpha)

    return result

def do_analytic(P, delta, start, end, l, m, thrld, step):
    result = []
    while start <= end:
        vector_p = solve_fxp(P, delta, start, l, m, thrld)
        start += step
        result.append(vector_p)
        print vector_p
    return result

if __name__ == "__main__":
    """
        result 1.5 ['0.5524', '1', '0.9958', '0.9841', '0.9542', '0.8599', '1', '0.792', '0.5847', '0.3795', '0.1801']
    """
    np.set_printoptions(precision=4)
    MAX_TRANS = 5
    DELTA = 0.000001
    P = np.zeros(MAX_TRANS)
    P[0] = 1
    # 注意： 这里 门限值一定是 分贝单位
    # THRLD = -3.0102999566398121
    THRLD = 6.0

    alpha_start = 0.30
    l, m = 1, 2
    alpha_end = 0.70
    step = 0.005
    rho = 1

    result = do_analytic(P, DELTA, alpha_start, alpha_end, l, m, THRLD, step)

    print result
    #
    result_f = "analytical_result_threshold={0}dB_l={1}_m={2}_sigma={3}.csv".format(THRLD, l, m, 0)

    with open(result_f, 'w') as f_handler:
        spamwriter = csv.writer(f_handler, delimiter=',')
        for n, vector_p in enumerate(result, 1):
            print n, vector_p
            spamwriter.writerow(vector_p)