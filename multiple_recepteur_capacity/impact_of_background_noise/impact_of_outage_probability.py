__author__ = 'qsong'


# This script is the latest one (15/09/2017).
import numpy as np
from numpy import sqrt, exp, pi
from scipy import optimize
from sympy.solvers import solve
from sympy import Symbol
import matplotlib.pyplot as plt
import matplotlib
import scipy.integrate as integrate

from scipy.special import erfc as erfc
from analytical_model import sgam
from scipy.special import gamma as gamma_f

import os


# Some common constant declaration
EPS = np.finfo(np.float64).eps
BETA = np.log(10.0)/10.0
FIGSIZE = (12, 8)

def sc_diversity_packet_loss_rate(r, thin_lambda_m, lambda_b, normalized_n, sigma_dB, thetha_dB):

    THETA = np.power(10, thetha_dB/10.0)
    sigma = BETA*sigma_dB
    sigma_X = 2.0*sigma/4.0
    fm_shadowing = np.exp(0.5*sigma_X**2)


    eta_T = normalized_n * THETA
    # B = p lamda_m pi A theta^{2/gamma} exp(sigma**2/8)
    # Attention! here A is for pure ALOHA, maximum interference
    A = 2 * gamma_f(0.5) * gamma_f(1.5)
    epsilon_T =thin_lambda_m * pi * np.power(THETA, .5) * fm_shadowing * A
    # Term 1 is the packet loss rate with only the nearest one.
    term1 = 1 - np.exp(-eta_T*r**4 - epsilon_T*r**2)

    term2 = np.exp(-np.pi * lambda_b * fm_shadowing * integrate.quad(lambda x: x*exp(-eta_T*x**4 - epsilon_T * x**2), r, 100)[0])

    # print "r:", r, "term1:", term1, "term2:", term2

    return term1 * term2

def best_bs_attach_packet_loss_rate(r, thin_lambda_m, lambda_b, normalized_n, sigma_dB, thetha_dB):

    THETA = np.power(10, thetha_dB/10.0)
    sigma = BETA*sigma_dB
    sigma_X = 2.0*sigma/4.0

    fm_shadowing = np.exp(0.5*sigma_X**2)
    eta_T = normalized_n * THETA
    # epsilon_T = p lamda_m pi A theta^{2/gamma} exp(sigma**2/8)
    # Take the value of A for best case,
    A = gamma_f(0.5) * gamma_f(1.5)
    epsilon_T =thin_lambda_m * pi * np.power(THETA, 0.5) * fm_shadowing * A
    # Term 1 is the packet loss rate with only the nearest one.
    term1 = 1 - np.exp(-eta_T*np.power(r, 4.0))*np.exp(-epsilon_T*np.power(r, 2.0))

    return term1


def calculate_critical_distance(lambda_b, thin_lambda_m, normalized_n, sigma_dB, thetha_dB, p_f_target, plr_func):
    """

        r:          intial value
        plr_func:       the packet loss rate function
    """
    MAX_ITR_NB = 500
    RELATIVE_ERROR = 1.0
    n = 0
    r_max = 10.0
    r_min = 0.0
    r_mean = 0.0

    while n < MAX_ITR_NB and RELATIVE_ERROR > 0.001:

        r_mean = (r_max + r_min) / 2.0

        tmp = plr_func(r_mean, thin_lambda_m, lambda_b, normalized_n, sigma_dB, thetha_dB)

        # print "r_mean:", r_mean, "tmp:", tmp

        if tmp > p_f_target:
            r_max = r_mean
        elif tmp < p_f_target:
            r_min = r_mean

        # update relative error
        RELATIVE_ERROR = np.abs(tmp - p_f_target) / p_f_target
        n += 1

    return r_mean




if __name__ == "__main__":

    # Maximum supported normalized load with respect to outage probabilities
    FIG_DST = "/Users/qsong/Documents/Final_Thesis/phd-thesis-template-2.2.2/Chapter5/Figures"
    FIG_NAME = 'min_bs_intensity_new.eps'
    MARKER_EVERY = 1

    # The normalized noise power level
    normalized_n_dB = -5.4
    normalized_n = np.power(10, normalized_n_dB/10.0)

    p_f_target, sigma_dB, thetha_dB = 0.1, 8.0, 3.0

    sigma = BETA*sigma_dB
    THETA = np.power(10, thetha_dB/10.0)
    sigma_X = 2.0*sigma/4.0
    fm_shadowing = np.exp(0.5*sigma_X**2)


    # define the charges level
    thin_lambda_m = 0.005

    min_bs_intensity_array = []

    START = 0.001
    END = 0.1
    p_outages = np.linspace(START, END, 50)

    # Now we change another strategy to get the relationship between outage probability and minimum BS intensity
    # We first guess that an interval of minimum BS intensity, then for each given BS intensity, calculate the
    # calculate the critical
    # distance from the packet loss rate formula, then calculate the outage probability.
    bs_intensity = np.linspace(0.4, 80, 1000)

    proba_l_sc = []
    proba_l_best = []

    proba_l_sc_no_noise = []
    proba_l_best_no_noise = []


    thin_lambda_m_1 = 0.5
    thin_lambda_m_2 = 0.5

    for element in bs_intensity:
        # c_d_src refers to the critical distance when using SC macro reception diversity
        c_d_src = calculate_critical_distance(element, thin_lambda_m_1, normalized_n, sigma_dB, thetha_dB, p_f_target, sc_diversity_packet_loss_rate)
        c_d_best = calculate_critical_distance(element, thin_lambda_m_1, normalized_n, sigma_dB, thetha_dB,p_f_target, best_bs_attach_packet_loss_rate)
        # normalized_n = 0.0 => no noise
        normalized_n = 0.0
        c_d_src_no_noise = calculate_critical_distance(element, thin_lambda_m_2, normalized_n, sigma_dB, thetha_dB, p_f_target, sc_diversity_packet_loss_rate)
        c_d_best_no_noise = calculate_critical_distance(element, thin_lambda_m_2, normalized_n, sigma_dB, thetha_dB, p_f_target, best_bs_attach_packet_loss_rate)
        proba_o_sc = np.exp(-np.pi * element * fm_shadowing * c_d_src**2)
        proba_o_best = np.exp(-np.pi * element * fm_shadowing * c_d_best**2)

        proba_o_sc_no_noise = np.exp(-np.pi * element * fm_shadowing * c_d_src_no_noise**2)
        proba_o_best_no_noise = np.exp(-np.pi * element * fm_shadowing * c_d_best_no_noise**2)
        print "sc macro:", "element:", element, "critical distance:", c_d_src, "outage_probability:", proba_o_sc
        print "best BS attach:", "element:", element, "critical distance:", c_d_best, "outage_probability:", proba_o_best

        proba_l_sc.append(proba_o_sc)
        proba_l_best.append(proba_o_best)

        proba_l_sc_no_noise.append(proba_o_sc_no_noise)
        proba_l_best_no_noise.append(proba_o_best_no_noise)


    print "proba_l", proba_l_sc

    slotted_A = gamma_f(0.5) * gamma_f(1.5)

    # Section: Draw the performance curve.
    fig, axes = plt.subplots(1, 1, figsize=FIGSIZE, sharey=False)
    axes.plot(
        proba_l_sc,
        bs_intensity,
        color='r',  marker='d', linestyle='-.', linewidth=2, markevery=MARKER_EVERY+3, label="SC Macro diversity, pure ALOHA, \n maximum interference, $\lambda_m = {0}$".format(thin_lambda_m_1)
    )
    axes.plot(
        proba_l_best,
        bs_intensity,
        color='b',  marker='o', linestyle='-.', linewidth=2, markevery=MARKER_EVERY+3, label="Best BS attach, slotted ALOHA.$\lambda_m = {0}$".format(thin_lambda_m_1)
    )

    axes.plot(
        proba_l_sc_no_noise,
        bs_intensity,
        color='r',  marker='+', linestyle='--', linewidth=2, markevery=MARKER_EVERY+5, label="SC Macro diversity, pure ALOHA, \n maximum interference, $\lambda_m = {0}$".format(thin_lambda_m_2)
    )
    axes.plot(
        proba_l_best_no_noise,
        bs_intensity,
        color='b',  marker='s', linestyle='--', linewidth=2, markevery=MARKER_EVERY+5, label="Best BS attach, slotted ALOHA, $\lambda_m = {0}$".format(thin_lambda_m_2)
    )

    # Study the impact of background noise
    # axes.plot(
    #     proba_l_sc,
    #     thin_lambda_m_2/bs_intensity,
    #     color='r',  marker='', linestyle='--', linewidth=2, markevery=MARKER_EVERY+3, label="SC Macro diversity, pure ALOHA, \n maximum interference, $\lambda_m = {0}$".format(thin_lambda_m_1)
    # )
    # axes.plot(
    #     proba_l_best,
    #     thin_lambda_m_2/bs_intensity,
    #     color='b',  marker='o', linestyle='--', linewidth=2, markevery=MARKER_EVERY+3, label="Best BS attach, slotted ALOHA.$\lambda_m = {0}$".format(thin_lambda_m_1)
    # )
    #
    # axes.plot(
    #     proba_l_sc_no_noise,
    #     thin_lambda_m_2/bs_intensity,
    #     color='r',  marker='+', linestyle='--', linewidth=2, markevery=MARKER_EVERY+5, label="SC Macro diversity, pure ALOHA, \n maximum interference, $\lambda_m = {0}$".format(thin_lambda_m_2)
    # )
    # axes.plot(
    #     proba_l_best_no_noise,
    #     thin_lambda_m_2/bs_intensity,
    #     color='b',  marker='s', linestyle='--', linewidth=2, markevery=MARKER_EVERY+5, label="Best BS attach, slotted ALOHA, $\lambda_m = {0}$".format(thin_lambda_m_2)
    # )

    # axes.plot(
    #     p_outages,
    #     np.power(slotted_A*np.power(THETA, 2.0/4.0), -1) * np.log(1-p_f_target)/np.log(p_outages),
    #     color='b',  marker='s', linestyle='-', linewidth=2, markevery=MARKER_EVERY, label="Best BS attach, slotted ALOHA, no noise"
    # )




    axes.set_yscale("linear")
    axes.set_xscale('log')
    axes.set_xlim([START, END])
    # axes.set_ylim([0.0, 2.5])
    axes.set_xlabel(r"Outage Probability")
    axes.set_ylabel(r"Minimum BS intensity")
    axes.grid(which='major')
    axes.grid(which='minor')
    # plt.minorticks_on()
    plt.legend(loc='best')
    fig.savefig(os.path.join(FIG_DST, FIG_NAME), format='eps', dpi=300)


    plt.show()
