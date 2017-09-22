__author__ = 'qsong'

import numpy as np
import os
from scipy.special import gamma as gamma_f
from scipy.special import erf as erf

from analytical_model import sgam
import matplotlib.pyplot as plt

from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
from mpl_toolkits.axes_grid1.inset_locator import inset_axes


import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
# from analytical_model.sgam import bs_nearest_atch_op, bs_rx_div_op

params = {
    'legend.fontsize': 20,
    'xtick.labelsize': 20,
    'ytick.labelsize': 20,
    'axes.labelsize': 20,
}
plt.rcParams.update(params)

# gamma = 4.0
#
# r = 1.0
#
# thetha_dB = 3.0
# threshold = np.power(10, thetha_dB/10)
# intensity_b = 0.004
# p = 0.008
# intensity_m = 0.01
# sigma_dB = 8.0
# BETA = np.log(10)/10
# sigma = sigma_dB * BETA
# sigma_X = 2.0/gamma*sigma
# # length = 30000
# # chi = np.random.lognormal(0, sigma, length)
# A = gamma_f(1+2.0/gamma) * np.exp(0.5*np.power(sigma_X, 2)) * np.power(threshold, 2.0/gamma)
# K = np.pi * gamma_f(1-2.0/gamma)


FIGSIZE = (10, 8)


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



# x = np.linspace(-100, 100, 100000)
# y = logistical(x)*pdf(x)
#
# # plt.plot(x, logistical(x), 'b')
# # plt.plot(x, logistical2(x), 'r')
#
# gamma_2 = np.sqrt(1+np.pi*sigma_X**2/8)

# print "K", K, "B", B
# print "p_{f,n}=", 1-np.mean(f_chi)
# print "p_{f,n}, trap, ",  1-np.trapz(y, x)
# print "p_{f,n}, trap, ",  np.trapz(logistical2(x)*pdf(x), x)
# print "p_{f,n}, trap, logis-normal-app ",  np.trapz(q(x)*pdf(x), x)
# print "p_{f,n}, trap, logis-normal-app-further ",  0.5 + 0.5*erf(np.log(B)*np.sqrt(np.pi)/(4*gamma_2))
# print "p_{f,n}, repara",  repara_logist(np.log(B))
# print "p_{f,n}, repara2",  repara_logist2(B)
# print "p_{f,n}, lower bound",  1-np.power(1+np.exp(np.log(B) + 0.5*sigma_X), -1)
#
# print "p_{f,n}, upper bound",  1-np.power(1+np.exp(np.log(B) - 0.5*sigma_X), -1)
#
# print "p_{s}(r)",  np.mean(psr(1, chi))
# print "p_{f}(best)",  1 - np.power(B+1, -1)
#
#
# print "p_{f}(best)",  1 - np.power(B+1, -1)


if __name__ == "__main__":

    FIG_DEST = "/Users/qsong/Documents/Final_Thesis/phd-thesis-template-2.2.2/Chapter5/Figures"
    FIG_NAME = 'comparison_monte_carlo_approximation.eps'
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

    fig, ax1 = plt.subplots(1, 1, figsize=FIGSIZE, sharey=False)

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
    ax1.set_yscale('log')
    # plt.title('Comparison of Analytical Packet Loss Rate and Monte-Carlo Simulation Results')
    ax1.grid(True, which="both", axis="both")
    ax1.plot(p*intensity_m/intensity_b, p_f_bs_nst_att_8, 'b', linestyle='-',label="Analytical Results")
    ax1.plot(p*intensity_m/intensity_b, p_f_mont_carlo, 'r', linestyle=':', label="Monte-Carlo Simulation Results")
    ax1.set_ylabel("Packet loss rate")
    ax1.set_xlabel("Normalized load")
    ax1.set_xticks(np.arange(0.02, 0.3, 0.03))
    # plt.yticks(0.07, 0.08, 0.09, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7)
    ax1.set_xlim([0.02, 0.29])
    ax1.set_ylim([0.07, 0.6])




    axins1 = inset_axes(ax1, 3, 4, bbox_to_anchor=[650, 250], loc=7)
    axins1.plot(p*intensity_m/intensity_b, p_f_bs_nst_att_8, 'b', linestyle='-', label="Analytical Results")
    axins1.plot(p*intensity_m/intensity_b, p_f_mont_carlo, 'r', linestyle=':', label="Monte-Carlo Simulation Results")

    x1, x2, y1, y2 = 0.02, 0.04, 0.09, 0.14 # specify the limits
    axins1.set_xlim(x1, x2) # apply the x-limits
    axins1.grid()
    axins1.set_ylim(y1, y2) # apply the y-limits
    mark_inset(ax1, axins1, loc1=2, loc2=4, fc="none", ec="0.5")

    ax1.legend(loc='best')

    plt.savefig(os.path.join(FIG_DEST, FIG_NAME), format='eps', transparent=False, dpi=300)

    plt.show()









