# -*- coding: utf-8 -*-
__author__ = 'qsong'

import csv
import matplotlib.pyplot as plt
import os
import numpy as np
import matplotlib.ticker as ticker
import numpy as np
import glob
import scipy.stats as st
import pandas as pd
import glob
import re
from scipy.special import gamma as gamma_f
import scipy.special as ss
from analytical_model import mrc_curve_fitting
from analytical_model import sgam

from scipy.special import erf as erf
# from analytical_model.sgam import bs_nearest_atch_op, bs_rx_div_op

import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

params = {
    'legend.fontsize': 15,
    'xtick.labelsize': 20,
    'ytick.labelsize': 20,
    'axes.labelsize': 20,
    'legend.numpoints': 1
}
plt.rcParams.update(params)

FIGSIZE = (8, 6)

if __name__ == '__main__':

    # 生成 lambda_m 的 ndarray
    X_END = 5.0

    lambda_m = np.linspace(0, X_END, 190)
    # 原则上我们让 基站的密度是个肯定值，毕竟这个东西投资大，没必要变化 lambda_b
    lambda_b = 0.08
    # gamma, path loss component. Usually take 4 as value.
    p = 0.008

    L = p*lambda_m/lambda_b

    gamma = 3.7
    thetha_dB = 3.0 # unit dB
    theta = np.power(10.0, thetha_dB/10)

    mu_shadow = 0.0
    X_END = p*max(lambda_m)/lambda_b

    X_START = 0.0
    X_STEP = 0.002
    Y_START = 1e-3
    Y_END = 1.0
    Y_STEP = 0.1

    MAX_TRANS = 1
    LINEWIDTH = 3

    SCALE = ["log", "linear"]
    FIG_DST = "/Users/qsong/Documents/Communication_Letter_MRC_curve_fitting"
    FIG_NAME = 'pure_slot_packet_loss_rate_mpr_gamma_{0}_theta_{1}.eps'.format(3, thetha_dB)
    FIG_NAME = 'slot_ana_fit_and_sim_comp.eps'


    SIM_LOG_DIR = '/Users/qsong/Documents/slotted_aloha_related_project/logs/SimuNov23PlusdePoints'
    PURE, ITF_MEAN = True, True
    if PURE == False:
        FIG_TITLE = 'Slotted ALOHA'
    else:
        if ITF_MEAN == True:
            FIG_TITLE = 'Pure ALOHA, mean'
        else:
            FIG_TITLE = 'Pure ALOHA, max'

    SC_SHOW = False # Whether the curves about SC should be present in the figure.
    SIM_SC_SHOW = False
    MRC_SHOW = True # Whether the curves about MRC should be present in the figure.
    SIM_MRC_SHOW = True
    FIT_MRC_SHOW = False # Whether the fitted curves about MRC should be present in the figure.
    SIM_FIT_MRC_SHOW = False # Whether the fitted curves about MRC should be present in the figure.



    sim_intensity, sim_plr_sc_divers, sim_plr_sc_divers_semi_ci, sim_plr_mrc_divers, sim_plr_mrc_divers_semi_ci = \
        mrc_curve_fitting.sim_parser(SIM_LOG_DIR, PURE, ITF_MEAN)


    gammas = [3.0, 3.3, 3.5, 3.7, 4.0, 4.2, 4.5, 4.7, 5.0, 6.0]
    p_f_rx_div, p_f_mrc_div, fit_p_f_mrc_div, empirical_p_f_mrc_div \
        = mrc_curve_fitting.sc_mrc_anayltical_parser(lambda_m, lambda_b, p, thetha_dB, gammas, PURE, ITF_MEAN)

    sim_fit_p_f_mrc_div ={}
    sim_fit_log = '/Users/qsong/Documents/slotted_aloha_related_project/analytical_model/sim_fit_result_theta_3.csv'

    for gamma in gammas:
        label = str(int(gamma*10))
        sim_fit_p_f_mrc_div[label] = mrc_curve_fitting.sim_fitted_function(sim_fit_log, thetha_dB, L, gamma, PURE, ITF_MEAN)


    fig, axes = plt.subplots(1, 1, figsize=FIGSIZE, sharey=False)

    axes.set_yscale("log")

    MARCO_DIV_TYPE = ['MRC']

    # plot curves for sc
    ana_slt_gammas = [3.0, 3.3, 4.5, 5.0, 6.0]
    lines_styles = ['-', '-.', '--', ':', '-']
    for i, gamma in enumerate(ana_slt_gammas):
        label = str(int(gamma*10))
        if SC_SHOW:
            axes.plot(
            L,
            p_f_rx_div[label],
            color='r',  marker='', linestyle=lines_styles[i], linewidth=LINEWIDTH, label=r"SC,ANA,$\gamma={0}$".format(gamma)
            )

    # plot curves for MRC
        if MRC_SHOW:
            axes.plot(
                L,
                p_f_mrc_div[label],
                color='m',  marker='', linestyle=lines_styles[i], linewidth=LINEWIDTH, label="ANA,$\gamma={0}$".format(gamma)
            )

        if FIT_MRC_SHOW:
            axes.plot(
                L,
                fit_p_f_mrc_div[label],
                color='b',  marker='', linestyle=lines_styles[i], linewidth=LINEWIDTH, label="FIT_ANA,$\gamma={0}$".format(gamma)
            )

        if SIM_FIT_MRC_SHOW:
            axes.plot(
                L,
                sim_fit_p_f_mrc_div[label],
                color='b',  marker='', linestyle=lines_styles[i], linewidth=LINEWIDTH, label="FIT_SIM,$\gamma={0}$".format(gamma)
            )


    sim_slt_gammas = [3.0, 3.3, 4.5, 5.0, 6.0]
    sc_marker_styles = ['s', 'o', '*', '+', 'd', 'D']
    mrc_marker_styles = ['v', '^', '>', '<', 'h', 'p']

    for i, gamma in enumerate(sim_slt_gammas):
        label = str(int(gamma*10))
        if SIM_SC_SHOW:
            axes.errorbar(
                sim_intensity[label],
                sim_plr_sc_divers[label],
                yerr=[sim_plr_sc_divers_semi_ci[label],  sim_plr_sc_divers_semi_ci[label]],
                fmt=sc_marker_styles[i],
                mfc='none',
                ecolor='r',
                capthick=2,
                label="SC,SIM,$\gamma={0}$".format(gamma)
            )
    for i, gamma in enumerate(sim_slt_gammas):
        label = str(int(gamma*10))
        if SIM_MRC_SHOW:
            axes.errorbar(
                sim_intensity[label],
                sim_plr_mrc_divers[label],
                yerr=[sim_plr_mrc_divers_semi_ci[label],  sim_plr_mrc_divers_semi_ci[label]],
                fmt=mrc_marker_styles[i],
                mfc='none',
                ecolor='b',
                capthick=2,
                label=r"SIM,$\gamma={0}$".format(gamma)
            )

    axes.grid()
    axes.axis([X_START, X_END, Y_START, Y_END])
    axes.set_xlabel(r"Normalized Load")
    axes.set_ylabel("Packet Loss Rate")
    axes.set_title(FIG_TITLE)

    # lines = axes.get_lines()
    # legend1 = plt.legend([lines[i] for i in [0,1,2]], ["algo1", "algo2", "algo3"], loc=1)
    # legend2 = plt.legend([lines[i] for i in [0,3,6]], parameters, loc=4)
    # axes.add_artist(legend1)
    # axes.add_artist(legend2)

    # print zip(p*lambda_m/lambda_b, p_f_rx_div_0)
    print mrc_curve_fitting.empirical_plr_mrc(thetha_dB, 2.0, lambda_b, 4.0, p, PURE, ITF_MEAN)
    print mrc_curve_fitting.empirical_plr_mrc(thetha_dB, 2.0, lambda_b, 3.3, p, PURE, ITF_MEAN)
    # plt.legend(bbox_to_anchor=(0.119, 0.01, 0.79, 1), loc=1, ncol=3, mode="expand", bbox_transform=plt.gcf().transFigure)
    plt.legend(loc='best')
    # plt.legend(bbox_to_anchor=(0.119, 0.02, 0.79, 1), loc=1, ncol=3, mode="expand", bbox_transform=plt.gcf().transFigure)
    plt.savefig(os.path.join(FIG_DST, FIG_NAME), format='eps', dpi=300)

    plt.show()

