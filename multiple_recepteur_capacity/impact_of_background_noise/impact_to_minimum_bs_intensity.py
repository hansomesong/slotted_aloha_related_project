__author__ = 'qsong'

from scipy.special import gamma as gamma_f
import numpy as np
from scipy import optimize

import matplotlib.pyplot as plt
import matplotlib
from scipy.special import erf as erf
from scipy.special import erf as erfc
from analytical_model import sgam

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

    r_0_square = np.log(1.0/p_outage) *np.power(lambda_b_min*np.pi*fm_shadowing, -1.0)

    term_1 = np.pi * fm_shadowing * erfc(np.sqrt(A) * np.power(r_0_square, 2))
    term_2 = np.log(
        (1 - np.exp(-A*np.power(r_0_square, 2)))
        / p_f_target
    )

    return np.sqrt(A) * term_2 / term_1

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


if __name__ == "__main__":

    FIG_DST = "/Users/qsong/Documents/Final_Thesis/phd-thesis-template-2.2.2/Chapter5/Figures"
    FIG_NAME = 'min_bs_intensity.eps'
    MARKER_EVERY = 4
    normalized_n_dB = -5.4
    normalized_n = np.power(10, normalized_n_dB/10.0)

    p_f_target, sigma_dB, thetha_dB, gamma = 0.1, 8.0, 3.0, 4.0

    THETA = np.power(10, thetha_dB/10.0)

    min_bs_intensity_array = []

    p_outages = np.linspace(0.001, 0.1, 50)

    for p_outage in p_outages:
        tmp = optimize.fixed_point(func_min_bs_intensity, [1.0], args=(normalized_n, p_outage, thetha_dB, gamma, p_f_target))

        min_bs_intensity_array.append(tmp[0])


    min_bs_intensity_array = np.array(min_bs_intensity_array)

    print "min_bs_intensity_array:", min_bs_intensity_array

    # Section: Draw the performance curve.
    fig, axes = plt.subplots(1, 1, figsize=FIGSIZE, sharey=False)

    axes.plot(
        p_outages,
        sgam.min_bs_intensity_best_case(p_outages, p_f_target, normalized_n, gamma, thetha_dB, sigma_dB),
        color='b',  marker='o', linestyle='--', linewidth=2, markevery=MARKER_EVERY, label="Best BS attach, slotted ALOHA"
    )

    axes.plot(
        p_outages,
        min_bs_intensity_array,
        color='r',  marker='*', linestyle='--', linewidth=2, markevery=MARKER_EVERY, label="SC Macro diversity, pure ALOHA, \n maximum interference"
    )

    axes.set_yscale("linear")
    axes.set_xscale('log')
    axes.set_xlim([0.001, 0.1])
    axes.set_ylim([0.0, 3.5])
    axes.set_xlabel(r"Outage Probability")
    axes.set_ylabel(r"Minimum BS intensity")
    axes.grid(which='major')
    axes.grid(which='minor')
    # plt.minorticks_on()
    plt.legend(loc='best')
    plt.savefig(os.path.join(FIG_DST, FIG_NAME), format='eps', dpi=300)
    plt.show()






