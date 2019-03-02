# -*- coding: utf-8 -*-
# This script collects all involved method for maximum ratio combining.

__author__ = 'qsong'
import numpy as np
from scipy.special import gamma as gamma_f
from scipy.special import erf as erf # Import error function
from scipy.special import binom as binom
from scipy.optimize import leastsq  # 引入最小二乘法算法

import scipy.stats as st
import scipy.integrate as integrate
import pandas as pd


SIM_FILE_PATH = "/Users/qsong/Documents/slotted_aloha_related_project/logs/SimuAvecFittageDifTetas"


def cumu_sir_lt(s, lambda_m, lambda_b, gamma, p, sigma_dB, pure=False, itf_mean=True):
    """
    This method is used to calculate the value of Laplace transform, given a complex value s
    :param s:
    :param lambda_m:
    :param lambda_b:
    :param gamma:
    :param p:
    :param sigma_dB:
    :param pure:
    :param itf_mean:
    :return:
    """
    # Calculate the fractional moment of cumulative interference
    fm_cumu_itf = itf_frac_moment_calculator(lambda_m, p, sigma_dB, gamma, pure, itf_mean)

    BETA = np.log(10.0)/10.0
    sigma = BETA*sigma_dB
    sigma_X = 2.0*sigma/gamma

    # Fractional moment of shadowing effect
    fm_shadowing = np.exp(0.5*sigma_X**2)

    C = np.pi*gamma_f(1-2.0/gamma)*gamma_f(1+2.0/gamma)

    return np.exp(
        -lambda_b * fm_shadowing * C * fm_cumu_itf * np.power(s, 2.0/gamma)
    )


def itf_frac_moment_calculator(lambda_m, p, sigma_dB, gamma, pure, itf_mean):
    # 2017-10-12. I have found the general expression for fractional moment calculation problem.
    BETA = np.log(10.0)/10.0
    sigma = BETA*sigma_dB
    sigma_X = 2.0*sigma/gamma

    # Fractional moment of shadowing effect
    fm_shadowing = np.exp(0.5*sigma_X**2)
    A = gamma_f(1-2.0/gamma)*gamma_f(1+2.0/gamma)

    if pure:
        if itf_mean:
            A = 2*gamma/(2+gamma) * A
        else:
            A = 2.0 * A

    return gamma*np.power(2*gamma_f(2.0/gamma)*p*lambda_m*np.pi*fm_shadowing*A, -1.0)


def lt2plr(thetha_dB, lambda_m, lambda_b, gamma, p, sigma_dB, pure, itf_mean):
    '''
    This method is used to numerically calculate packet loss rate from laplace transform.
    The underlying formula originiates from Eq.(35) of reference: "Unified Stochastic Geometry Modeling and Analysis of
    Cellular Networks in LOS/NLOS and Shadowed Fading".

    :param thetha_dB:       scalar, in unit of dB, capture ratio threshold
    :param lambda_m:        scalar,
    :param lambda_b:
    :param gamma:
    :param p:
    :param sigma_dB:
    :param pure:
    :param itf_mean:
    :return:
    '''
    T = np.power(10.0, thetha_dB/10.0)
    #
    # To avoid any name collision with our namespace. For parameters used in above mentioned is added with append _num
    A_NUM = 18.4
    M_NUM = 11
    N_NUM = 15

    def real_part_lt_of_cdf_theta(x, lambda_m, lambda_b, gamma, p, sigma_dB, pure, itf_mean):
        return np.real(
            cumu_sir_lt(x, lambda_m, lambda_b, gamma, p, sigma_dB, pure, itf_mean)/x
        )

    binom_coeff = np.array([binom(M_NUM, k) for k in range(M_NUM + 1)])
    s_n_k_t_list = []
    for k in range(M_NUM+1):
        a_coeff = np.exp(0.5*A_NUM) * np.array(
            [
                np.power(2*T, -1)
                    if l == 0
                    else
                np.power(T, -1)
                    for l in range(N_NUM + k + 1)
            ]
        )

        b_coeff = np.array(
            [
                np.power(-1.0, l)*
                real_part_lt_of_cdf_theta(
                    (A_NUM+2j*l*np.pi)/2/T,
                    lambda_m, lambda_b, gamma, p, sigma_dB, pure, itf_mean
                )
                for l in range(N_NUM + k + 1)
            ]
        )

        s_n_k_t_list.append(np.sum(a_coeff * b_coeff))

    s_n_k_t = np.array(s_n_k_t_list)
    result = np.sum(binom_coeff * s_n_k_t) * np.power(2.0, -1.0*M_NUM)

    return result




def file_parse(csv_file):
    df = pd.read_csv(csv_file, sep=";", decimal=',', skiprows=9, nrows=15).dropna(axis=1, how='all')

    print [float(x.replace(',', '.')) for x in df.columns[1:].values.tolist()]
    print df


def analytical_mrc_writter(thetha_dB, gamma, sigma_dB, pure, itf_mean):
    p = 0.01
    lambda_b = 150.0
    charges = np.array([0.02, 0.04, 0.06, 0.08, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5])
    lambda_m = lambda_b * charges/p
    plr = []
    for element in lambda_m:
        tmp = lt2plr(thetha_dB, element, lambda_b, gamma, p, sigma_dB, pure, itf_mean)
        plr.append(tmp)

    return np.array(plr)

# Methods for curve fitting for Maximum Ratio Combining techniques
## The target function to fit.
## Given that
## 1) when path-loss exponent \gamma is 4.0, the packet loss rate when MRC is applied is a function of erfc function
## 2) Laplace transform of output SIR in case of MRC is in terms of exponential function and depends
##      on A*L*theta^{2/\gamma}
#
def func(k, x):
    return k*x

## 偏差函数：x,y都是列表:这里的x,y更上面的Xi,Yi中是一一对应的
def error(k, x, y):
    return func(k, x) - y





if __name__ == "__main__":
    gamma = 3.8
    sigma_dB = 5.0
    thetha_dB = 4.0
    lambda_m = 3000.0
    lambda_b = 150.0
    p = 0.01
    # L = 0.2
    pure = True
    itf_mean = True
    print lt2plr(thetha_dB, lambda_m, lambda_b, gamma, p, sigma_dB, pure, itf_mean)



    first_empirical_results = analytical_mrc_writter(thetha_dB, gamma, sigma_dB, pure, itf_mean)
    print "first empirical results:", first_empirical_results[5:]

    print analytical_mrc_writter(thetha_dB, gamma, sigma_dB, pure, itf_mean)

    csv_file = "/Users/qsong/Documents/slotted_aloha_related_project/logs/SimuAvecFittageDifTetas/AlohaMultg33s8t4.csv"

    file_parse(csv_file)

    tmp = np.array([0.01444, 0.06472, 0.22408,	0.3752, 	0.48454])




    reference = np.array([6.01044850e-03,   4.82264534e-02,   2.11996231e-01, 3.65075870e-01,   4.78565375e-01])


    print np.sqrt(np.sum(np.power((reference - tmp)/reference, 2.0)))
    print np.sqrt(np.sum(np.power((reference - first_empirical_results[5:])/reference, 2.0)))


    Para=leastsq(error, p0, args=(Xi,Yi))



