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

    # print W_MULTIPLIER
    return np.exp(alpha*(P.dot(W_MULTIPLIER)-sum(P)))

def cdf_from_char_func_FFT(thld, alpha, P, sigma, v, k, omega_end, eta, N):
    '''
    % cdf_from_char_func_FFT produces CDF via FFT of alpha-stable characteristic function
    % Decay factor eta necessary to avoid divergence at zero
        @eta            : the damping function parameter, e^{-eta}
        @N              : the frequency samples number
        @omega_end      : the maximum value of omega (frequency)
        @thld           : SINR threshold, in unit dB
    '''
    thld = 10**(thld/10.0)
    omega = np.linspace(0, omega_end, N)
    delta_omega = omega[1]-omega[0]
    EPS = np.finfo(float).eps

    dx = 2*np.pi/(delta_omega*N)
    xpoint = dx*(N-1)/2.0
    x = np.arange(-xpoint, xpoint+dx, dx)

    shift_omega = 1j*eta + omega + EPS
    print "omega", omega
    print "shift_omega:", shift_omega

    cdf_phi = ((1+EPS*1j)/(eta-1j*omega))*chf_compound_sln(alpha, P, sigma, shift_omega, v, k)

    fft_input = delta_omega*np.exp(-1j*x[0]*omega)*cdf_phi

    # Trapezoidal Rule: 1st and last divided by 2
    fft_input[0] = fft_input[0]/2.0
    fft_input[-1] = fft_input[-1]/2.0
    print "cdf special:", (np.exp(eta*0.51)/np.pi)*np.real(np.fft.fft(fft_input))
    cdf = (np.exp(eta*x)/np.pi)*np.real(np.fft.fft(fft_input))
    print cdf[0:10]
    print "cdf max", max(cdf), "cdf min", min(cdf)
    # Normalize
    cdf = (cdf-min(cdf))/max(cdf)
    print cdf[0:10]
    x_close_2_thrld = min(x, key=lambda element: abs(element-1.0/thld))
    print x
    print "cdf", cdf
    print x_close_2_thrld
    index_x = np.where(x == x_close_2_thrld)[0][0]
    result = 1.0 - cdf[index_x]
    print "Failure Proba", result
    return result


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
    omega_end, eta, N = 314, 1.0, 2048
    # print "Current power increment factor is:", v
    # print "For intensity:", alpha
    while i < MAX_I and ecart > delta:
        # print "In iteration ", i, "The proba vector P is:\t\t", ["{:10.5f}".format(x) for x in list(P)]
        # P[:-1] is a slicing operation in numpy, meaning taking elements except the last one.
        failure = [cdf_from_char_func_FFT(thrld, alpha, P, sigma, v, k, omega_end, eta, N) for k in range(K-1)]
        # print "In iteration ", i, "New failure proba:\t\t\t", ["{:10.5f}".format(x) for x in failure]
        # print "In iteration ", i, "New success proba:\t\t\t", ["{:10.5f}".format(1-x) for x in failure]

        tmp[1:] = [p_f*P[x] for x, p_f in enumerate(failure)]
        ecart = sum(abs(tmp-P))
        P = np.copy(tmp)
        i += 1
    # Note that np.floor(l**(K-1)/threld-M[-1])
    result = []
    p_loss = P[-1]*cdf_from_char_func_FFT(thrld, alpha, P, sigma, v, K-1, omega_end, eta, N)
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
        print vector_p
    return result

if __name__ == "__main__":

    MAX_TRANS = 5
    DELTA = 0.000001
    SIGMA = 1
    P = np.zeros(MAX_TRANS)
    P[0] = 1
    # 注意： 这里 门限值一定是 分贝单位
    THRLD = 3.0
    alpha_start = 0.8
    l, m = 2, 1
    alpha_end = 0.804
    step = 0.005
    rho = 1

    # result = do_analytic(P, DELTA, alpha_start, alpha_end, l, m, SIGMA, rho, THRLD, step)
    #
    # print result
    #
    # result_f = "fft_shadowing_analytical_result_threshold={0}dB_l={1}_m={2}_sigma={3}.csv".format(THRLD, l, m, SIGMA)
    #
    # with open(result_f, 'w') as f_handler:
    #     spamwriter = csv.writer(f_handler, delimiter=',')
    #     for n, vector_p in enumerate(result, 1):
    #         print n, vector_p
    #         spamwriter.writerow(vector_p)



    k = 0
    v = 1.0
    sigma = 1.0
    alpha = 1.5
    omega_end, eta, N = 50, 0.8, 256
    thld = 3.0
    # # # print chf_cumu_interference(0.4, 1000, P, k, v, sigma)
    # # # print trans_failure_p_1(3.0, P, 0.4, v, k, rho, sigma)
    # #
    # #
    # plt.figure()
    # plt.subplot(2, 1, 1)
    # t= np.linspace(0.0001, 60, 10000)
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

    # chf = chf_compound_sln(alpha, P, sigma, t, v, k)
    # real_part_2 = np.real(chf)
    # imag_part_2 = np.imag(chf)
    # print real_part_2
    # print imag_part_2
    # plt.title("2")
    # plt.plot(t, real_part_2)
    # plt.plot(t, imag_part_2)
    # plt.grid()
    # plt.show()


    cdf_from_char_func_FFT(thld, alpha, P, sigma, v, k, omega_end, eta, N)


