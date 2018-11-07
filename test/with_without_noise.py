__author__ = 'qsong'

import numpy as np
from scipy.special import gamma as gamma_f
from scipy.special import erf as erf

from analytical_model import sgam
import matplotlib.pyplot as plt


FIGSIZE = (8, 6)


def psr(x, chi):
    # x, the distance
    # chi, np.power()
    return np.exp(-p*intensity_m * A * K * x**2 * (np.exp(-chi))**(2.0/gamma))

def pdf(x):
    """
    pdf is actually the probability density function of a RV of normal distribution. Its standard deviation error is
    noted as \sigma_X and mean as \log(B), where B = A * K * p * \lambda_m / (pi * \lambda_b ),
    A = gamma( 1 + 2.0/\gamma) * np.exp(0.5 \sigma_X ^ 2) * \theta ^{ 2.0/\gamma},
    K depends on slotted or pure ALOHA. In slotted ALOHA, k = pi * gamma(1-2.0/\gamma)
    :param x:           a vector
    :return:
    """
    tmp1 = np.exp(-0.5 * np.power((x - np.log(B)) / sigma_X, 2))
    tmp2 = np.power(np.sqrt(2 * np.pi) * sigma_X, -1.0)
    return tmp2*tmp1


def logistical(x):
    return np.power(1+np.exp(x),  -1)

def logistical2(x):
    return np.power(1+np.exp(-x),  -1)

def q(x):
    return 0.5 + 0.5*erf(x/1.595769/np.sqrt(2))

def logist(x):
    return np.power(1+np.exp(-x),  -1)


def repara_logist(x):
    gamma_3 = np.sqrt(1+np.pi*sigma_X**2/8)
    return np.power(1 + np.exp(-x/gamma_3), -1)

def repara_logist2(x):
    gamma_3 = -np.power(1+np.pi*sigma_X**2/8, -0.5)
    return np.power(1 + x**gamma_3, -1)






if __name__ == "__main__":

    intensity_m = 0.01
    gamma = 4.0
    thetha_dB = 3.0
    threshold = np.power(10, thetha_dB/10)
    intensity_b = 0.08
    p = 0.008
    N = 100
    intensity_m = np.linspace(0.21, 3, N)
    sigma_dB = 8.0
    BETA = np.log(10)/10
    sigma = sigma_dB * BETA
    sigma_X = 2.0/gamma*sigma
    length = 50000
    A = gamma_f(1+2.0/gamma) * np.exp(0.5*np.power(sigma_X, 2)) * np.power(threshold, 2.0/gamma)
    K = np.pi * gamma_f(1-2.0/gamma)
    B = np.tile(A*K*p*intensity_m/(np.pi*intensity_b), (length, 1))
    print B.shape


    p_f_bs_nst_att_8 = sgam.bs_nearest_atch_op(intensity_m, intensity_b, gamma, p, thetha_dB, 8)

    chi = np.random.lognormal(0, sigma, size=(length, N))
    f_chi = np.power(B*np.power(chi, -2.0/gamma) + 1, -1.0)
    p_f_mont_carlo = 1-np.mean(f_chi, axis=0)

    diff = np.abs(p_f_bs_nst_att_8-p_f_mont_carlo)

    # fig, axes = plt.subplots(1, 1, figsize=FIGSIZE, sharey=False)

    print "p_f_bs_nst_att_8"
    print p_f_bs_nst_att_8
    print "p_f_mont_carlo"
    print p_f_mont_carlo
    print "=========="
    print np.max(diff)
    print "The ratio..."
    print "maximum ratio:", np.max(100*np.array(diff/p_f_mont_carlo))
    print 100*np.array(diff/p_f_mont_carlo)
    max_idx = np.argmax(diff)
    print max_idx
    print p*intensity_m[int(max_idx)]/intensity_b
    plt.yscale('log')
    # plt.title('Comparison of Analytical Packet Loss Rate and Monte-Carlo Simulation Results')
    plt.grid(True, which="both", axis="both")
    plt.plot(p*intensity_m/intensity_b, p_f_bs_nst_att_8, 'b', label="Analytical Results")
    plt.plot(p*intensity_m/intensity_b, p_f_mont_carlo, 'r', label="Monte-Carlo Simulation Results")
    plt.ylabel("Packet loss rate")
    plt.xlabel("Normalized load")
    plt.xticks(np.arange(0.02, 0.3, 0.03))
    # plt.yticks(0.07, 0.08, 0.09, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7)
    plt.axis([0.02, 0.3, 0.07, 0.7])
    plt.legend(loc='best')
    plt.savefig('comparison_monte_carlo_approximation.eps', format='eps', dpi=300)
    plt.show()









