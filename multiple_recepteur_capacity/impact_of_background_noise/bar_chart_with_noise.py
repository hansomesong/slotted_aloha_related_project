__author__ = 'qsong'

from scipy.special import gamma as gamma_f
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

from matplotlib import rc
rc('text', usetex=True)

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
params = {
    'legend.fontsize': 15,
    "lines.markersize" : 8,
    'lines.linewidth' : 2,
    'xtick.labelsize': 15,
    'ytick.labelsize': 15,
    'axes.labelsize': 15,
    'legend.numpoints': 1
}
plt.rcParams.update(params)

FIGSIZE = (12, 8)


def cdistance2bs_intensity(r_c, shadowing_f, p_o):
    """
        Calculate the minimum BS intensity from a critical distance and outage probability threshold.
    :return:
    """

    return np.log(1.0/p_o) * np.power(np.pi * shadowing_f * r_c**2, -1.0)

def fixed_p_fuc_best_bs(r_c, p_f, eta_t, epsilon_t):
    """
    The fixed point equation to solve the critical distance for best BS attach method. Consider the background noise
    The formula:
        r_c = sqrt(1.0/epsilon_T * { -log(1-p_f) - eta_T * r_c^4}

    Attention: the initial value of r_c is very important! For example, let r_c = 1.0, finally we won't get the solution.
    if we let r_c = 0.001, we get solution...
    :return:
    """
    return np.sqrt((-np.log(1-p_f) - eta_t * r_c**4) / epsilon_t)

def sc_diversity_packet_loss_rate(r, p_o, eta, epsilon, shadowing_f):

    lambda_b = np.log(1.0/p_o) * np.power(np.pi * shadowing_f * r**2, -1.0)

    term1 = 1 - np.exp(-eta*r**4.0 - epsilon*r**2.0)

    term2 = exp(-2 * pi * lambda_b * shadowing_f * integrate.quad(lambda x: x*exp(-eta*x**4 - epsilon * x**2), r, 20)[0])

    # print "r:", r, "term1:", term1, "term2:", term2

    return term1 * term2


def calculate_critical_distance4sc_macro(eta, epsilon, p_f, shadowing_f):
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

        tmp = sc_diversity_packet_loss_rate(r_mean, p_o, eta, epsilon, shadowing_f)

        # print "r_mean:", r_mean, "tmp:", tmp

        if tmp > p_f:
            r_max = r_mean
        elif tmp < p_f:
            r_min = r_mean

        # update relative error
        RELATIVE_ERROR = np.abs(tmp - p_f) / p_f
        n += 1

    return r_mean

def best_bs_attach_packet_loss_rate(r, eta, epsilon):

    return 1 - np.exp(-eta*np.power(r, 4.0) - epsilon*np.power(r, 2.0))

if __name__ == "__main__":
    FIG_DST = "/Users/qsong/Documents/Final_Thesis/phd-thesis-template-2.2.2/Chapter5/Figures"
    p_f = 0.1   #target network level packet loss rate
    p_o = 0.1   #outage probability threshold
    FIG_NAME = 'bar_chart_min_bs_intensity_load_{0}'.format(p_o)

    BETA = np.log(10.0)/10.0
    sigma_dB = 8.0 # unit of dB
    sigma = sigma_dB * BETA

    thetha_dB = 3.0 # unit of dB
    theta = np.power(10, thetha_dB/10.0)
    # To analyze the case without noise, just let norm_noise_power = 0.0
    norm_noise_power_dB = -5.4 # unit of dB
    norm_noise_power = np.power(10.0, norm_noise_power_dB/10.0)
    # The same no matter which association policy.
    eta_t = norm_noise_power * theta
    # eta_t =0.0
    shadowing_f = np.exp(sigma**2/8.0)


    # for macro diversity, constant A is 2*gamma(0.5)*gamma(1.5) = pi
    loads = np.array([0.05, 0.1, 0.15, 0.2, 0.25])  #load = p * lambda_m
    epsilon_t_sc = loads * np.pi * np.pi * shadowing_f * np.power(theta, 0.5)
    # for best BS attach, constant A is gamma(0.5)*gamma(1.5) = pi/2.0
    epsilon_t_best = loads * np.pi * (np.pi/2.0) * shadowing_f * np.power(theta, 0.5)



    critical_dists =  optimize.fixed_point(
        fixed_p_fuc_best_bs,
        [0.001]*loads.size,
        args=(
            p_f*np.ones(loads.size),
            eta_t*np.ones(loads.size),
            epsilon_t_best
        )
    )

    critical_dists_no_noise = optimize.fixed_point(
        fixed_p_fuc_best_bs,
        [0.001]*loads.size,
        args=(
            p_f*np.ones(loads.size),
            0.0*np.ones(loads.size),
            epsilon_t_best
        )
    )
    print "Critical Distance vector:", critical_dists, "Packet loss rate vector:", \
        best_bs_attach_packet_loss_rate(critical_dists, eta_t, epsilon_t_best)

    mini_bs_intensity_best = cdistance2bs_intensity(critical_dists, shadowing_f, p_o)
    mini_bs_intensity_best_no_noise = cdistance2bs_intensity(critical_dists_no_noise, shadowing_f, p_o)

    print mini_bs_intensity_best

    mini_bs_intensity_sc = []
    mini_bs_intensity_sc_no_noise = []

    for element in epsilon_t_sc:
        sc_c_dist = calculate_critical_distance4sc_macro(eta_t, element, p_f, shadowing_f)
        sc_c_dist_no_noise = calculate_critical_distance4sc_macro(0.0, element, p_f, shadowing_f)

        tmp = cdistance2bs_intensity(sc_c_dist, shadowing_f, p_o)
        tmp_no_noise = cdistance2bs_intensity(sc_c_dist_no_noise, shadowing_f, p_o)

        print "SC critical distance:", sc_c_dist, "target loss rate:", \
            sc_diversity_packet_loss_rate(sc_c_dist, p_o, eta_t, element, shadowing_f),\
            "bs intensity:", tmp
        mini_bs_intensity_sc.append(tmp)
        mini_bs_intensity_sc_no_noise.append(tmp_no_noise)

    mini_bs_intensity_sc = np.array(mini_bs_intensity_sc)
    mini_bs_intensity_sc_no_noise = np.array(mini_bs_intensity_sc_no_noise)

    # Linear fitting
    z1 = np.polyfit(loads, mini_bs_intensity_best, 1)
    p1 = np.poly1d(z1)

    print p1


    width = 0.23       # the width of the bars


    # Plotting part. The tutorial: https://matplotlib.org/examples/api/barchart_demo.html
    fig, ax = plt.subplots(figsize=FIGSIZE)
    ind = np.arange(loads.size)

    rects1 = ax.bar(ind, mini_bs_intensity_best, width, edgecolor = "k", color='w')
    rects1_no_noise = ax.bar(ind+width, mini_bs_intensity_best_no_noise, width, edgecolor = "k", color='b')
    rects2 = ax.bar(ind+2*width, mini_bs_intensity_sc, width, edgecolor = "k", color='r')
    rects2_no_noise = ax.bar(ind+3*width, mini_bs_intensity_sc_no_noise, width, edgecolor = "k", color='m')

    for element in rects1:
        element.set_hatch('*')
    for element in rects2:
        element.set_hatch('o')
    #

    # add some text for labels, title and axes ticks
    ax.set_ylabel(r'Mininum BS spatial density $\lambda_b$ (per $km^2$)')
    ax.set_xlabel(r'Network load $p\lambda_m$ (per slot)')
    ax.set_xticks(ind + 1.5*width)
    ax.set_xticklabels(['%.2f' % load for load in loads])
    ax.set_ylim([0, int(np.max(mini_bs_intensity_best)+1)])

    ax.legend(
        (
            rects1[0],
            rects1_no_noise[0],
            rects2[0],
            rects2_no_noise[0]
        ),
        (
            'Best BS attach, \nslotted ALOHA',
            'Best BS attach, \nslotted ALOHA, no noise',
            'SC macro diversity,\nPure ALOHA \nwith maximum interference',
            'SC macro diversity,\nPure ALOHA \nwith maximum interference, no noise'
        )
    )

    def autolabel(rects):
        """
        Attach a text label above each bar displaying its height
        """
        for rect in rects:
            height = rect.get_height()
            ax.text(rect.get_x() + rect.get_width()/2., 1.00*height,
                    '%.2f' % height, fontsize=15,
                    ha='center', va='bottom')

    autolabel(rects1)
    autolabel(rects2)
    autolabel(rects1_no_noise)
    autolabel(rects2_no_noise)

    plt.tight_layout()
    fig.savefig(os.path.join(FIG_DST, FIG_NAME), format='eps', dpi=300)
    plt.show()