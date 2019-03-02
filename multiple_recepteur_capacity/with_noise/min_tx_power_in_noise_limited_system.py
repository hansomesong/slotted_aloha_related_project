__author__ = 'qsong'

from scipy.special import gamma as gamma_f
import numpy as np
from scipy import optimize

import matplotlib.pyplot as plt
import matplotlib
from scipy.special import erf as erf
from scipy.special import erf as erfc
from analytical_model import sgam

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


def min_tx_power_fixed_func(p_t, lambda_b, thetha_dB, gamma, r, p_f_target):
    """
        This is a fixed point function about minimum transmit power, which is to be solved.
        Note that the background noise is regarded as 1.0
    """
    THETA = np.power(10, thetha_dB/10)
    sigma = BETA*sigma_dB
    sigma_X = 2.0*sigma/gamma
    fm_shadowing = np.exp(0.5*sigma_X**2)

    A = THETA * np.power(r, 2) / p_t

    term_1 = THETA * np.power(r, 2)
    term_2 = np.pi * lambda_b * fm_shadowing * erfc(np.sqrt(A) * np.power(r, 2))
    term_3 = np.log(( 1 - np.exp(-A*np.power(r, 2))) / p_f_target )

    return term_1 * np.power(term_3/term_2, 2)


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

    MARKER_EVERY = 4

    plr, lambda_b, sigma_dB, thetha_dB, gamma = 0.1, 0.08, 8, 3, 4

    THETA = np.power(10, thetha_dB/10)

    fig, axes = plt.subplots(1, 1, figsize=FIGSIZE, sharey=False)


    p_outages = np.linspace(0.01, 0.1, 50)

    r_array = outage_distance(p_outages, lambda_b, sigma_dB, thetha_dB, gamma)

    print "the distance:", r_array

    min_tx_power_array = []

    for r in r_array:
        tmp = optimize.fixed_point(min_tx_power_fixed_func, [1.0], args=(lambda_b, thetha_dB, gamma, r, plr))
        min_tx_power_array.append(tmp[0])


    min_tx_power_array = np.array(min_tx_power_array)
    slotted_A = gamma_f(1-2.0/gamma)*gamma_f(1+2.0/gamma)

    print "min_tx_power_array ", min_tx_power_array

    print sgam.min_tx_power_best_case(p_outages, plr, lambda_b, gamma, thetha_dB, sigma_dB)

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
        sgam.min_tx_power_best_case(p_outages, plr, lambda_b, gamma, thetha_dB, sigma_dB),
        color='b',  marker='o', linestyle='--', linewidth=2, markevery=MARKER_EVERY, label="Best BS attach"
    )

    axes.plot(
        p_outages,
        min_tx_power_array,
        color='r',  marker='*', linestyle='--', linewidth=2, markevery=MARKER_EVERY, label="SC Macro diversity"
    )



    print slotted_A*np.power(THETA, 2.0/gamma)
    print "xx", 1.0/slotted_A*np.power(THETA, 2.0/gamma)

    axes.set_yscale("linear")
    axes.set_xscale('log')

    axes.set_xlabel(r"Outage Probability")
    axes.set_ylabel(r"Minimum Tx Power")
    axes.grid()
    plt.legend(loc='best')
    plt.savefig('min_tx_power.eps', format='eps', dpi=300)
    plt.show()






