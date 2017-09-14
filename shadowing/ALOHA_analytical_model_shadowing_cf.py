# -*- coding: utf-8 -*-
__author__ = 'qsong'

import numpy as np
import matplotlib.pyplot as plt
import csv
import scipy.special
from scipy.special import lambertw

def chf_compound_sln(alpha, P, sigma, omega, v, k):
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
    BETA = np.log(10)/10.0
    DIMENSION = P.size
    W_MULTIPLIER = np.array([
        np.power(1+lambertw(-2j*omega*np.power(v, m-k)*np.power(BETA*sigma, 2), 0), -0.5)*
        np.exp(
            -1.0*(
                2*lambertw(-2j*omega*np.power(v, m-k)*np.power(BETA*sigma, 2), 0)
                + np.power(lambertw(-2j*omega*np.power(v, m-k)*np.power(BETA*sigma, 2), 0), 2)
            )/np.power(2*BETA*sigma, 2)
        )

        for m in range(DIMENSION)
    ])

    # W_MULTIPLIER = np.array([
    #     np.power(1-lambertw(2j*omega*np.power(v, m-k)*np.power(BETA*sigma, 2), 0), -0.5)*
    #     np.exp(
    #         -1.0*(
    #             -2*lambertw(2j*omega*np.power(v, m-k)*np.power(BETA*sigma, 2), 0)
    #             + np.power(lambertw(2j*omega*np.power(v, m-k)*np.power(BETA*sigma, 2), 0), 2)
    #         )/np.power(2*BETA*sigma, 2)
    #     )
    #
    #     for m in range(DIMENSION)
    # ])

    # print W_MULTIPLIER
    # print P.shape

    # W_MULTIPLIER = np.transpose(W_MULTIPLIER)

    return np.exp(alpha*(P.dot(W_MULTIPLIER)-sum(P)))

def cdf_from_char_func_FFT(alpha, P, sigma, v, k, omega_end, eta, N):
    '''
    % cdf_from_char_func_FFT produces CDF via FFT of alpha-stable characteristic function
    % Decay factor eta necessary to avoid divergence at zero
        @eta            : the damping function parameter, e^{-eta}
        @N              : the frequency samples number
        @omega_end      : the maximum value of omega (frequency)
    '''
    omega = np.linspace(0, omega_end, N)
    delta_omega = omega[1]-omega[0]
    EPS = np.finfo(float).eps

    dx = 2*np.pi/(delta_omega*N)
    xpoint = (dx*(N-1)/2.0)
    x = np.arange(-xpoint, xpoint+dx, dx)

    shift_omega = eta +1j*omega+ EPS*1j

    cdf_phi = ((1+EPS*1j)/(eta-1j*omega))*chf_compound_sln(alpha, P, sigma, omega, v, k)

def chf_poisson(w, alpha):
    return np.exp(alpha*(np.exp(1j*w)-1))

def chf_log_normal(w, k, v, sigma):
    BETA = np.log(10)/10.0
    return np.power(1+lambertw(-2j*w*np.power(v, m-k)*np.power(BETA*sigma, 2), 0), -0.5)*np.exp(-1.0*(2*lambertw(-2j*w*np.power(v, m-k)*np.power(BETA*sigma, 2), 0)+ np.power(lambertw(-2j*w*np.power(v, m-k)*np.power(BETA*sigma, 2), 0), 2))/np.power(2*BETA*sigma, 2))

def chf_cumu_interference(alpha, w, P, k, v, sigma):
    '''
        本方法用于计算某一给定设备在第 k th 重新传输时候，其遭遇的总干涉强度的特征函数
    '''
    BETA = np.log(10)/10.0
    DIMENSION = P.size
    W_MULTIPLIER = np.array([
        np.power(1+lambertw(-2j*w*np.power(v, m-k)*np.power(BETA*sigma, 2), 0), -0.5)*
        np.exp(
            -1.0*(
                2*lambertw(-2j*w*np.power(v, m-k)*np.power(BETA*sigma, 2), 0)
                + np.power(lambertw(-2j*w*np.power(v, m-k)*np.power(BETA*sigma, 2), 0), 2)
            )/np.power(2*BETA*sigma, 2)
        )

        for m in range(DIMENSION)
    ])

    # print W_MULTIPLIER
    # print P.shape

    # W_MULTIPLIER = np.transpose(W_MULTIPLIER)

    return np.exp(alpha*(P.dot(W_MULTIPLIER)-sum(P)))

def trans_failure_p_1(thrld, P, alpha, v, k, rho, sigma):
    thrld = 10**(thrld/10.0)
    w_array = np.linspace(0.0001, 60, 2000)


    # Y = [
    #     (np.exp(1j*w/thrld)*chf_cumu_interference(alpha, -w, P, k, v, sigma) -
    #     np.exp(-1j*w/thrld)*chf_cumu_interference(alpha, w, P, k, v, sigma))/(1j*w)
    #     for w in w_array
    # ]
    #
    # Y = np.real(Y)
    #
    # result = 0.5+np.trapz(Y, x=w_array)/(2*np.pi)

    Y = [
        (np.imag(chf_cumu_interference(alpha, w, P, k, v, sigma))*np.cos(w/thrld)
         - np.real(chf_cumu_interference(alpha, w, P, k, v, sigma))*np.sin(w/thrld))/w
        for w in w_array
    ]


    result = 0.5 - np.trapz(Y, x=w_array)/np.pi
    # print "proba:", 1-result

    return 1-result

def trans_failure_p_2(thrld, P, alpha, v, k, sigma):
    thrld = np.power(10, thrld/10.0)
    omega_end = 50
    N = 2000
    omega = np.linspace(0, omega_end, N)
    eta = 1.0
    EPS = np.finfo(float).eps
    shift_omega = 1j*eta + omega + EPS
    Y = ((1+EPS*1j)/(eta-1j*omega))*np.exp(-1j*omega/thrld)*chf_compound_sln(alpha, P, sigma, shift_omega, v, k)

    result = np.exp(eta*1.0/thrld)*np.real(np.trapz(y=Y, x=omega))/np.pi
    # print "proba:", 1-result

    return 1-result

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
    failure_prob = 0.0001
    # print "Current power increment factor is:", v
    # print "For intensity:", alpha
    while i < MAX_I and ecart > delta:
        # print "In iteration ", i, "The proba vector P is:\t\t", ["{:10.5f}".format(x) for x in list(P)]
        # P[:-1] is a slicing operation in numpy, meaning taking elements except the last one.
        failure = [trans_failure_p_2(thrld, P, alpha, v, k, sigma) for k in range(K-1)]
        # print "In iteration ", i, "New failure proba:\t\t\t", ["{:10.5f}".format(x) for x in failure]
        # print "In iteration ", i, "New success proba:\t\t\t", ["{:10.5f}".format(1-x) for x in failure]

        tmp[1:] = [p_f*P[x] for x, p_f in enumerate(failure)]
        # ecart = sum(abs(tmp-P))
        P = np.copy(tmp)
        i += 1
        tmp_failure_p = P[-1]*trans_failure_p_2(thrld, P, alpha, v, K-1, sigma)
        ecart = abs(tmp_failure_p - failure_prob)/failure_prob
        failure_prob = tmp_failure_p + np.finfo(float).eps

        # print "vector P", P, P[-1]*trans_failure_p_2(thrld, P, alpha, v, K-1, sigma)

    # Note that np.floor(l**(K-1)/threld-M[-1])
    result = []
    p_loss = P[-1]*trans_failure_p_2(thrld, P, alpha, v, K-1, sigma)
    result.extend([float("{0:.4g}".format(element)) for element in list(P)])
    result.append(p_loss)
    result.append(alpha)

    return result

def do_analytic(P, delta, start, end, l, m, sigma, thrld, step):
    result = []
    while start <= end:
        vector_p = solve_fxp(P, delta, start, 1.0*l/m, sigma, thrld)
        start += step
        result.append(vector_p)
        print vector_p
    return result

if __name__ == "__main__":

    MAX_TRANS = 5
    DELTA = 0.00001
    SIGMA = 2.0 # 3dB
    P = np.zeros(MAX_TRANS)
    P[0] = 1
    # 注意： 这里 门限值一定是 分贝单位
    THRLD = -3.0
    alpha_start = 0.60
    l, m = 2, 1
    alpha_end = 2.05
    step = 0.005

    result = do_analytic(P, DELTA, alpha_start, alpha_end, l, m, SIGMA, THRLD, step)


    # plt.figure()
    # P = np.array(result[0][0:5])
    # x_axis = np.arange(0+0.000000000000000000000000001, 6.0, 0.02)
    # x_axis_2 = np.array([-10.0*np.log10(x) for x in x_axis])
    # y = [1-trans_failure_p_2(x, P, alpha_start, 1, 0, SIGMA) for x in x_axis_2]
    # print trans_failure_p_2(-3.0, P, alpha_start, 1, 0, SIGMA)
    # plt.plot(x_axis, y)
    #
    # plt.show()
    # print result
    #
    result_f = "analytical_result_K={4}_threshold={0}dB_l={1}_m={2}_sigma={3}.csv".format(THRLD, l, m, SIGMA, MAX_TRANS)

    with open(result_f, 'w') as f_handler:
        spamwriter = csv.writer(f_handler, delimiter=',')
        for n, vector_p in enumerate(result, 1):
            print n, vector_p
            spamwriter.writerow(vector_p)

    # k = 0
    # w = 0
    # v = 1.0
    # sigma = 1.0
    # alpha = 0.4
    # # # print chf_cumu_interference(0.4, 1000, P, k, v, sigma)
    # # # print trans_failure_p_1(3.0, P, 0.4, v, k, rho, sigma)
    # #
    # #
    # plt.figure()
    # plt.subplot(2, 1, 1)
    # t= np.linspace(0.0001, 100, 10000)
    # chf_1 = [chf_cumu_interference(alpha, w, P, k, v, sigma) for w in t]
    # real_part_1 = np.real(chf_1)
    # imag_part_1 = np.imag(chf_1)
    #
    # # real_part2 = np.real([chf_log_normal(w, k, v, sigma)for w in t])
    # # imag_part2 = np.imag([chf_log_normal(w, k, v, sigma)for w in t])
    # #
    # # real_part3 = np.real([chf_poisson(w, 0.4) for w in t])
    # # imag_part3 = np.imag([chf_poisson(w, 0.4) for w in t])
    #
    # plt.plot(t, real_part_1)
    # plt.plot(t, imag_part_1)
    # plt.title("1")
    # plt.subplot(2, 1, 2)
    #
    # chf = chf_compound_sln(alpha, P, sigma, t, v, k)
    # real_part_2 = np.real(chf)
    # imag_part_2 = np.imag(chf)
    # print real_part_2
    # print imag_part_2
    # plt.title("2")
    # plt.plot(t, real_part_2)
    # plt.plot(t, imag_part_2)
    # plt.show()


