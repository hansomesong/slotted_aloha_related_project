__author__ = 'qsong'

from scipy.special import gamma as gamma_f
import numpy as np
from scipy import optimize
from sympy.solvers import solve
from sympy import Symbol
import matplotlib.pyplot as plt
import matplotlib

from scipy.special import erfc as erfc
from analytical_model import sgam
from scipy.special import gamma as gamma_f


import os

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

def func_min_bs_intensity(lambda_b_min, normalized_n, p_outage, thetha_dB, gamma, p_f_target):
    """
    This is a fixed point function about minimum BS intensity, which is to be solved.
    :param lambda_b_min:            the minimum BS intensity to be solved
    :param normalized_n:            the normalized background noise
    :param thetha_dB:
    :param gamma:
    :param r:
    :param p_f_target:
    :return:lambda_b_min            the fixed point equation
    """
    THETA = np.power(10, thetha_dB/10.0)
    sigma = BETA*sigma_dB

    sigma_X = 2.0*sigma/gamma
    fm_shadowing = np.exp(0.5*sigma_X**2)

    A = normalized_n * THETA

    r_0_square = np.log(1.0/p_outage) * np.power(lambda_b_min*np.pi*fm_shadowing, -1.0)

    term_1 = np.pi * fm_shadowing * erfc(np.sqrt(A) * r_0_square) * np.sqrt(np.pi) / 2.0
    term_2 = np.log(
        (1 - np.exp(-A*np.power(r_0_square, 2)))
        / p_f_target
    )

    return np.sqrt(A) * term_2 / term_1


def func_critical_distance(r, lambda_b_min, normalized_n, thetha_dB, gamma, p_f_target):
    THETA = np.power(10, thetha_dB/10.0)
    sigma = BETA*sigma_dB
    sigma_X = 2.0*sigma/gamma

    fm_shadowing = np.exp(0.5*sigma_X**2)
    c_lambda_b = lambda_b_min * fm_shadowing
    A = normalized_n * THETA

    term1 = np.power(A, -1.0)
    term2 = p_f_target*np.exp(np.pi * c_lambda_b * np.sqrt(term1) * erfc(np.sqrt(A) * r))
    term3 = -np.log(1 - term2)

    # Attention: the return value is actually is the square of critical distance.
    return np.sqrt(-1.0 * term1 * term3)

def packet_loss_rate(r, lambda_b_min, normalized_n, thetha_dB, gamma, sigma_dB):
    THETA = np.power(10, thetha_dB/10.0)
    sigma = BETA*sigma_dB
    sigma_X = 2.0*sigma/gamma

    fm_shadowing = np.exp(0.5*sigma_X**2)
    c_lambda_b = lambda_b_min * fm_shadowing
    A = normalized_n * THETA

    term1 = 1 - np.exp(-A*r**4)
    term2 = np.exp(-np.pi * lambda_b_min * fm_shadowing * 0.5 * np.sqrt(np.pi) * np.power(A, -0.5) * erfc(np.sqrt(A) * r**2))

    return term1 * term2

def outage_prob_with_bs_intensity():
    pass


def calculate_critical_distance(lambda_b_min, normalized_n, thetha_dB, gamma, p_f, sigma_dB):
    """
        r:          intial value
    """
    MAX_ITR_NB = 500
    RELATIVE_ERROR = 1.0
    n = 0
    r_max = 10.0
    r_min = 0.0
    r_mean = 0.0

    while n < MAX_ITR_NB and RELATIVE_ERROR > 0.001:

        r_mean = (r_max + r_min) / 2.0

        tmp = packet_loss_rate(r_mean, lambda_b_min, normalized_n, thetha_dB, gamma, sigma_dB)

        # print "r_mean:", r_mean, "tmp:", tmp

        if tmp > p_f:
            r_max = r_mean
        else:
            r_min = r_mean

        # update relative error
        RELATIVE_ERROR = np.abs(tmp - p_f) / p_f
        n += 1

    return r_mean




if __name__ == "__main__":

    FIG_DST = "/Users/qsong/Documents/Final_Thesis/phd-thesis-template-2.2.2/Chapter5/Figures"
    FIG_NAME = 'bs_intensity_vs_outage_proba.eps'
    MARKER_EVERY = 4
    normalized_n_dB = -5.4
    normalized_n = np.power(10, normalized_n_dB/10.0)

    p_f_target, sigma_dB, thetha_dB, gamma = 0.1, 8.0, 3.0, 4.0
    sigma = BETA*sigma_dB
    THETA = np.power(10, thetha_dB/10.0)
    sigma_X = 2.0*sigma/gamma
    fm_shadowing = np.exp(0.5*sigma_X**2)

    min_bs_intensity_array = []

    START = 0.001
    END = 0.1
    p_outages = np.linspace(START, END, 50)



    # Section: Draw the performance curve.
    fig, axes = plt.subplots(1, 1, figsize=FIGSIZE, sharey=False)


    # Now we change another strategy to get the relationship between outage probability and minimum BS intensity
    # We first guess that an interval of minimum BS intensity, then for each given BS intensity, calculate the
    # calculate the critical
    # distance from the packet loss rate formula, then calculate the outage probability.
    bs_intensity = np.linspace(0.4, 1.5, 50)

    proba_l = []
    proba_l_0 = []

    for element in bs_intensity:

        c_d_8 = calculate_critical_distance(element, normalized_n, thetha_dB, gamma, p_f_target, sigma_dB)
        c_d_0 = calculate_critical_distance(element, normalized_n, thetha_dB, gamma, p_f_target, 0)
        proba_8 = np.exp(-np.pi * element * fm_shadowing * c_d_8**2)
        proba_0 = np.exp(-np.pi * element * fm_shadowing * c_d_0**2)

        print "element:", element, "critical distance 8dB:", c_d_8, "outage_probability 8dB:", proba_8
        print "element:", element, "critical distance 0dB:", c_d_0, "outage_probability 0dB:", proba_0

        proba_l.append(proba_8)
        proba_l_0.append(proba_0)


    x = np.linspace(0.4, 3.5, 100)

    y_8 = sgam.min_bs_intensity_best_case(p_outages, p_f_target, normalized_n, gamma, thetha_dB, sigma_dB)
    y_0 = sgam.min_bs_intensity_best_case(p_outages, p_f_target, normalized_n, gamma, thetha_dB, 0)
    axes.plot(
        y_8,
        p_outages,
        color='b',  marker='o', linestyle='-', linewidth=2, markevery=MARKER_EVERY, label="Best BS attach, slotted ALOHA, shadowing: 8 dB"
    )

    axes.plot(
        y_0,
        p_outages,
        color='b',  marker='', linestyle='-', linewidth=2, markevery=MARKER_EVERY, label="Best BS attach, slotted ALOHA, no shadowing"
    )


    axes.plot(
        bs_intensity,

        proba_l,

        color='r',  marker='*', linestyle='-.', linewidth=2, markevery=MARKER_EVERY, label="SC Macro diversity, pure ALOHA, \n maximum interference"
    )

    axes.plot(
        bs_intensity,

        proba_l_0,

        color='r',  marker='', linestyle='-.', linewidth=2, markevery=MARKER_EVERY, label="SC Macro diversity, pure ALOHA, \n maximum interference, no shadowing"
    )

    axes.set_xscale("linear")
    axes.set_yscale('log')
    axes.set_ylim([START, END])
    axes.set_xlim([0.0, 6])
    axes.set_ylabel(r"Outage Probability")
    axes.set_xlabel(r"Minimum BS intensity")
    axes.grid(which='major')
    axes.grid(which='minor')
    # plt.minorticks_on()
    plt.legend(loc='best')
    fig.savefig(os.path.join(FIG_DST, FIG_NAME), format='eps', dpi=300)

    plt.show()





