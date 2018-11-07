__author__ = 'qsong'
"""
    This script is used to draw the maximum supported load with respect to the outage probability constraint.
    The generated plot is inserted into the final thesis.
"""

from scipy.special import gamma as gamma_f
import numpy as np
from scipy import optimize
import os
import scipy.integrate as integrate

import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
params = {
    'legend.fontsize': 15,
    "lines.markersize" : 8,
    'lines.linewidth' : 2,
    'xtick.labelsize': 20,
    'ytick.labelsize': 20,
    'axes.labelsize': 20,
    'legend.numpoints': 1
}
plt.rcParams.update(params)

# Some common constant declaration
EPS = np.finfo(np.float64).eps
BETA = np.log(10.0)/10.0

FIGSIZE = (12, 8)

def outage_distance(p_outage, lambda_b, sigma_dB, thetha_dB, gamma):

    THETA = np.power(10, thetha_dB/10)
    sigma = BETA*sigma_dB
    sigma_X = 2.0*sigma/gamma
    fm_shadowing = np.exp(0.5*sigma_X**2)

    return np.sqrt(np.log(1.0/(p_outage)) / (lambda_b*np.pi*fm_shadowing))

def func(lambda_m, p, thetha_dB, gamma, r, lambda_b, p_outage, pure=False, itf_mean=True):
    """

        Fixed function to be solved.
    """
    THETA = np.power(10, thetha_dB/10)
    sigma = BETA*sigma_dB
    sigma_X = 2.0*sigma/gamma
    fm_shadowing = np.exp(0.5*sigma_X**2)
    L = p*lambda_m/lambda_b
    if not pure:
        A = gamma_f(1-2.0/gamma)*gamma_f(1+2.0/gamma)
    else:
        if itf_mean:
            factor = 2*gamma/(2+gamma)
        else:
            factor = 2.0
        A = factor*gamma_f(1-2.0/gamma)*gamma_f(1+2.0/gamma)

    term_1 = np.power(p*np.pi*A*np.power(THETA, 2.0/gamma)*fm_shadowing*np.power(r, 2), -1)
    term_2 = np.log(1-np.exp(-p*lambda_m*np.pi*A*np.power(THETA, 2.0/gamma)*fm_shadowing*np.power(r, 2)*np.power(r, 2))) - np.log(p_outage)
    term_3 = np.log(A*np.power(THETA, 2.0/gamma)*L*term_2)
    result = -1.0*term_1*term_3
    return result


def func_plr(plr, p, thetha_dB, gamma, lambda_b, lambda_m, N, pure=False, itf_mean=True):

    THETA = np.power(10, thetha_dB/10)
    sigma = BETA*sigma_dB
    sigma_X = 2.0*sigma/gamma
    fm_shadowing = np.exp(0.5*sigma_X**2)
    thetha_term = np.power(THETA, 2.0/gamma)
    L = p*lambda_m/lambda_b
    G = (1-plr)*L/(1-np.power(plr, 1.0/N))
    if not pure:
        A = gamma_f(1-2.0/gamma)*gamma_f(1+2.0/gamma)
    else:
        if itf_mean:
            factor = 2*gamma/(2+gamma)
        else:
            factor = 2.0
        A = factor*gamma_f(1-2.0/gamma)*gamma_f(1+2.0/gamma)

    return np.exp(np.sum([np.power(-1, i)*np.power(i*A*thetha_term*G, -1) for i in range(1, N+1)]))


def single_attach(p, thetha_dB, gamma, r, p_outage, pure=False, itf_mean=True):
    THETA = np.power(10, thetha_dB/10)
    sigma = BETA*sigma_dB
    sigma_X = 2.0*sigma/gamma
    fm_shadowing = np.exp(0.5*sigma_X**2)
    if not pure:
        A = gamma_f(1-2.0/gamma)*gamma_f(1+2.0/gamma)
    else:
        if itf_mean:
            factor = 2*gamma/(2+gamma)
        else:
            factor = 2.0
        A = factor*gamma_f(1-2.0/gamma)*gamma_f(1+2.0/gamma)

    term_2 = -1.0*np.log(1-p_outage)
    term_1 = np.power(p*np.pi*A*np.power(THETA, 2.0/gamma)*fm_shadowing*np.power(r, 2), -1)
    return term_1 * term_2


def sc_diversity_packet_loss_rate(r, lambda_b, normalized_n, thetha_dB, gamma, thin_lambda_m):

    THETA = np.power(10, thetha_dB/10.0)
    sigma = BETA*sigma_dB
    sigma_X = 2.0*sigma/gamma

    fm_shadowing = np.exp(0.5*sigma_X**2)
    eta_T = normalized_n * THETA
    # B = p lamda_m pi A theta^{2/gamma} exp(sigma**2/8)
    # Attention! here A is for pure ALOHA, maximum interference
    A = 2 * gamma_f(0.5) * gamma_f(1.5)
    epsilon_T =thin_lambda_m * np.pi * np.power(THETA, .5) * fm_shadowing * A
    # Term 1 is the packet loss rate with only the nearest one.
    term1 = 1 - np.exp(-eta_T*r**4 - epsilon_T*r**2)

    term2 = np.exp(-np.pi * lambda_b * fm_shadowing * integrate.quad(lambda x: np.exp(-eta_T*x**4 - epsilon_T * x**2), r, 10)[0])

    # print "r:", r, "term1:", term1, "term2:", term2

    return term1 * term2


def calculate_critical_distance(lambda_b, thin_lambda_m, normalized_n, thetha_dB, gamma, p_f, plr_func):
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

        tmp = plr_func(r_mean, lambda_b, normalized_n, thetha_dB, gamma, thin_lambda_m)

        # print "r_mean:", r_mean, "tmp:", tmp

        if tmp > p_f:
            r_max = r_mean
        elif tmp < p_f:
            r_min = r_mean

        # update relative error
        RELATIVE_ERROR = np.abs(tmp - p_f) / p_f
        n += 1

    return r_mean


if __name__ == "__main__":

    FIG_DST = "/Users/qsong/Documents/Final_Thesis/phd-thesis-template-2.2.2/Chapter5/Figures"
    FIG_NAME = 'max_normalized_load_with_outage_probability.eps'

    MARKER_EVERY = 4

    # vert strange: when fixed piont analysis does not work when thetha_dB =3.0 instead of 3!
    p_f_target, lambda_b, sigma_dB, thetha_dB, gamma = 0.1, 0.08, 8.0, 3.0, 4.0


    THETA = np.power(10, thetha_dB/10.0)


    sigma = BETA*sigma_dB
    sigma_X = 2.0*sigma/gamma
    fm_shadowing = np.exp(0.5*sigma_X**2)
    p = 0.008
    r = outage_distance(p_f_target, lambda_b, sigma_dB, thetha_dB, gamma)
    print r

    # print optimize.fixed_point(func, [1.2, 1.3], args=(p, thetha_dB, gamma, r, lambda_b, plr, True, False))
    print func(0.79123752, p, thetha_dB, gamma, r, lambda_b, p_f_target, True, False)

    fig, axes = plt.subplots(1, 1, figsize=FIGSIZE, sharey=False)


    X_END = 0.1
    p_outages = np.linspace(0.001, X_END, 50)
    r_array = outage_distance(p_outages, lambda_b, sigma_dB, thetha_dB, gamma)

    print "r_array", r_array

    lambda_m_array = []

    # for r in r_array:
    #     tmp = optimize.fixed_point(func, [1.0], args=(p, thetha_dB, gamma, r, lambda_b, p_f_target, True, False))
    #     lambda_m_array.append(tmp[0])
    #
    #
    # lambda_m_array = np.array(lambda_m_array)



    bs_intensity = np.linspace(0.4, 100, 100)
    thin_lambda_m = 0.5 # thus, the normalized L varies from 1.25 to 0.005
    normalized_n = 0.0 # ignorance of noise
    proba_l_sc = []
    proba_l_best = []

    for element in bs_intensity:
        # c_d_src refers to the critical distance when using SC macro reception diversity
        c_d_src = calculate_critical_distance(element, thin_lambda_m, normalized_n, thetha_dB, gamma, p_f_target, sc_diversity_packet_loss_rate)

        proba_o_sc = np.exp(-np.pi * element * fm_shadowing * c_d_src**2)
        print "sc macro:", "element:", element, "critical distance:", c_d_src, "outage_probability:", proba_o_sc

        proba_l_sc.append(proba_o_sc)


    print "proba_l", proba_l_sc

    slotted_A = gamma_f(1-2.0/gamma)*gamma_f(1+2.0/gamma)

    print lambda_m_array

    # axes.plot(
    #     p_outages,
    #     p*lambda_m_array/lambda_b,
    #     color='r',  marker='', linestyle='--', linewidth=2, label="Diversity,ANA, MAX_ITF"
    # )
    #

    # axes.plot(
    #     p_outages,
    #     p*single_attach(p, thetha_dB, gamma, r, p_outages, pure=False, itf_mean=True)/lambda_b,
    #     color='k',  marker='', linestyle='--', linewidth=2, label="Diversity,ANA, MAX_ITF"
    # )


    # axes.plot(
    #     np.linspace(0, 1, 50),
    #     1 - (1-np.linspace(0, 1, 50))**(1.0/(slotted_A*np.power(THETA, 2.0/gamma)*0.08)),
    #     color='b',  marker='', linestyle='--', linewidth=2, label="Cumulative Distribution Function of $P_{f, b}$"
    # )

    axes.plot(
        p_outages,
        np.power(slotted_A*np.power(THETA, 2.0/gamma), -1) * np.log(1-p_f_target)/np.log(p_outages),
        color='b',  marker='o', linestyle='-', linewidth=2, markevery=MARKER_EVERY, label="Best BS attach, slotted ALOHA"
    )

    # axes.plot(
    #     p_outages,
    #     0.05*np.ones(p_outages.size),
    #     color='b', linestyle='-', linewidth=2, markevery=MARKER_EVERY, label="Best BS attach, slotted ALOHA, no outage constraint"
    # )

    axes.plot(
        proba_l_sc,
        thin_lambda_m/bs_intensity ,
        color='r',  marker='*', linestyle='--', linewidth=2, markevery=MARKER_EVERY, label="SC Macro diversity, pure ALOHA, \n maximum interference"
    )

    # axes.plot(
    #     p_outages,
    #     0.1**np.ones(p_outages.size),
    #     color='r',  linestyle='--', linewidth=2, markevery=MARKER_EVERY, label="SC Macro diversity, pure ALOHA, \n maximum interference, no outage constraint"
    # )


    lambda_m_list  = np.linspace(0.00001, 3, 50)

    print slotted_A*np.power(THETA, 2.0/gamma)
    print "xx", 1.0/slotted_A*np.power(THETA, 2.0/gamma)
    print "L:", np.power(slotted_A*np.power(THETA, 2.0/gamma), -1) * np.log(1-p_f_target)/np.log(p_outages)

    axes.set_yscale("linear")
    axes.set_xscale('log')
    axes.set_xlim([0.001, X_END])
    axes.set_ylim([0.0, 0.10])
    axes.set_xlabel(r"Outage Probability")
    axes.set_ylabel(r"Maximum Supported Normalized Load")
    axes.grid(which='major')
    axes.grid(which='minor')
    plt.legend(loc='best')
    plt.savefig(os.path.join(FIG_DST, FIG_NAME), format='eps', dpi=300)
    plt.show()






