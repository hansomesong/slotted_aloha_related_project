__author__ = 'qsong'

from scipy.special import gamma as gamma_f
import numpy as np
from scipy import optimize
from sympy.solvers import solve
from sympy import Symbol
import matplotlib.pyplot as plt
import matplotlib
from scipy.special import erf as erf
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


if __name__ == "__main__":

    FIG_DST = "/Users/qsong/Documents/Final_Thesis/phd-thesis-template-2.2.2/Chapter5/Figures"
    FIG_NAME = 'min_bs_intensity.eps'
    MARKER_EVERY = 4
    normalized_n_dB = -5.4
    normalized_n = np.power(10, normalized_n_dB/10.0)

    p_f_target, sigma_dB, thetha_dB, gamma = 0.1, 8.0, 3.0, 4.0
    sigma = BETA*sigma_dB
    THETA = np.power(10, thetha_dB/10.0)

    min_bs_intensity_array = []

    START = 0.001
    END = 0.1
    p_outages = np.linspace(START, END, 50)

    for p_outage in p_outages:
        tmp = optimize.fixed_point(func_min_bs_intensity, [5.0], args=(normalized_n, p_outage, thetha_dB, gamma, p_f_target))

        min_bs_intensity_array.append(tmp[0])


    min_bs_intensity_array = np.array(min_bs_intensity_array)

    print "min_bs_intensity_array:", min_bs_intensity_array
    print sgam.min_bs_intensity_best_case(p_outages, p_f_target, normalized_n, gamma, thetha_dB, sigma_dB)

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


    l = 5.0
    p_outage = 0.001
    sigma_X = 2.0*sigma/gamma
    fm_shadowing = np.exp(0.5*sigma_X**2)
    r = np.log(1.0/p_outage) * np.power(l*np.pi*fm_shadowing, -1.0)

    print "l:", l, "fixed point:", func_min_bs_intensity(l, normalized_n, p_outage, thetha_dB, gamma, p_f_target)
    # print np.exp(-np.pi*0.4105*np.exp(sigma**2/8)*np.sqrt(np.pi)/2/np.sqrt(normalized_n))
    # print (1-np.exp(-normalized_n*THETA*np.power(r, 2))) * np.exp(-np.pi * 0.4105 * 0.5* np.sqrt(np.pi/normalized_n/THETA) * erfc( np.sqrt(normalized_n*THETA))*r)
    #
    # r = 5.0
    # print (1-np.exp(-normalized_n*THETA*np.power(r, 2))) * np.exp(-np.pi * 0.4105 * 0.5* np.sqrt(np.pi/normalized_n/THETA) * erfc( np.sqrt(normalized_n*THETA))*r)
    #
    # print (1-np.exp(-normalized_n*THETA*np.power(r, 2)))
    axes.set_yscale("linear")
    axes.set_xscale('log')
    axes.set_xlim([START, END])
    axes.set_ylim([0.0, 3.5])
    axes.set_xlabel(r"Outage Probability")
    axes.set_ylabel(r"Minimum BS intensity")
    axes.grid(which='major')
    axes.grid(which='minor')
    # plt.minorticks_on()
    plt.legend(loc='best')
    fig.savefig(os.path.join(FIG_DST, FIG_NAME), format='eps', dpi=300)



    fig, axes = plt.subplots(1, 1, figsize=FIGSIZE, sharey=False)

    l = 2.5

    r = np.linspace(0, 5, 100)

    p = 0.008
    l_m = 6.0
    l_b = 0.08
    L = p*l_m/l_b
    plr1 = (1-np.exp(-normalized_n*THETA*np.power(r, 4))) * np.exp(
        -np.pi * l * fm_shadowing* 0.5 * np.sqrt(np.pi)* erfc(
            np.sqrt(normalized_n*THETA)*r**2
        )
        /np.sqrt(normalized_n*THETA)
    )
    plr2 = (1-np.exp(-normalized_n*THETA*np.power(r, 4)))
    A = 2 * gamma_f(0.5) * gamma_f(1.5)
    plr3 = (1-np.exp(-p*l_m*np.pi*A*np.power(THETA, 0.5) * fm_shadowing * np.power(r, 2))) * np.exp(-np.exp(-p*l_m*np.pi*A*np.power(THETA, 0.5)*fm_shadowing*np.power(r, 2))/A/np.power(THETA, 0.5)/L)

    print "term in plr1:", np.exp(
        -np.pi *
        l *
        fm_shadowing*
        0.5 *
        np.sqrt(np.pi)*
        erfc(
            np.sqrt(normalized_n*THETA)*0**2
        ) / np.sqrt(normalized_n*THETA)
    )

    print np.exp(
        -np.pi *
        l *
        fm_shadowing*
        0.5 *
        np.sqrt(np.pi) *
        erfc(0.0)
        / np.sqrt(normalized_n*THETA)
    )



    # x = Symbol('x')
    # y = solve((1-np.e**(-normalized_n*THETA*np.power(np.log(1.0/p_outage) * np.power(x*np.pi*fm_shadowing, -1.0), 2)))*np.e**(-np.pi * l * fm_shadowing* 0.5 * np.sqrt(np.pi)* erfc(np.sqrt(normalized_n*THETA)* np.log(1.0/p_outage) * np.power(x*np.pi*fm_shadowing, -1.0))/np.sqrt(normalized_n*THETA))- 0.1, x)
    # print "y:", y
    # print "term in plr3:", np.exp(-np.exp(-p*l_m*np.pi*A*np.power(THETA, 0.5)*fm_shadowing*np.power(r, 2))/A/np.power(THETA, 0.5)/L)
    axes.plot(r, plr1)
    axes.plot(r, plr2)
    # axes.plot(r, plr3)

    # plt.show()



