# -*- coding: utf-8 -*-
__author__ = 'qsong'

import csv
import matplotlib.pyplot as plt
import os
import numpy as np
import matplotlib.ticker as ticker
import numpy as np
import scipy.stats as st
import pandas as pd
import glob
from scipy.special import erf as erf

from scipy.special import gamma as gamma_f
import scipy.special as ss
from analytical_model import sgam, mrc_curve_fitting
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

LOG_DIR = 'logs'
SUB_DIR = 'multiple_reception'
CASE_DIR = 'fading'
SUB_CASE_DIR = "bs_0.01_p_0.002_R_40|100"
DATA_FOLDED = '.'
FIGSIZE = (8, 6)

SCALE = ["log", "linear"]

if __name__ == '__main__':

    X_START = 0.0
    X_END = 5.0
    X_STEP = 0.002
    Y_START = 1e-3
    Y_END = 1.0
    Y_STEP = 0.1

    MAX_TRANS = 1
    LINEWIDTH = 3
    # 生成 lambda_m 的 ndarray
    lambda_m = np.linspace(0, X_END, 190)
    # 原则上我们让 基站的密度是个肯定值，毕竟这个东西投资大，没必要变化 lambda_b
    lambda_b = 0.08
    # gamma, path loss component. Usually take 4 as value.

    p = 0.008
    L = p*lambda_m / lambda_b

    mu_shadow = 0.0


    X_END = p*max(lambda_m)/lambda_b

    gamma = 3.7
    thetha_dB = 3.0 # unit dB
    theta = np.power(10, thetha_dB/10)

    FIG_DST = "/Users/qsong/Documents/Communication_Letter_MRC"
    FIG_NAME = 'pure_max_packet_loss_rate_mpr_gamma_{0}_theta_{1}.eps'.format(3, thetha_dB)

    # Analytical results for SC macro diversity
    pure_max_p_f_rx_div ={}
    pure_max_p_f_rx_div['30'] = sgam.bs_rx_div_op(lambda_m, lambda_b, 3.0, p, thetha_dB, 8, True, False)
    pure_max_p_f_rx_div['33'] = sgam.bs_rx_div_op(lambda_m, lambda_b, 3.3, p, thetha_dB, 8, True, False)
    pure_max_p_f_rx_div['35'] = sgam.bs_rx_div_op(lambda_m, lambda_b, 3.5, p, thetha_dB, 8, True, False)
    pure_max_p_f_rx_div['37'] = sgam.bs_rx_div_op(lambda_m, lambda_b, 3.7, p, thetha_dB, 8, True, False)
    pure_max_p_f_rx_div['40'] = sgam.bs_rx_div_op(lambda_m, lambda_b, 4.0, p, thetha_dB, 8, True, False)
    pure_max_p_f_rx_div['60'] = sgam.bs_rx_div_op(lambda_m, lambda_b, 6.0, p, thetha_dB, 8, True, False)

    # Analytical results for MRC macro diversity
    pure_max_p_f_mrc_div = {}
    pure_max_p_f_mrc_div['30'] = []
    pure_max_p_f_mrc_div['33'] = []
    pure_max_p_f_mrc_div['35'] = []
    pure_max_p_f_mrc_div['37'] = []
    pure_max_p_f_mrc_div['40'] = []
    pure_max_p_f_mrc_div['60'] = []
    for l in L:
        pure_max_p_f_mrc_div['30'].append(mrc_curve_fitting.lt2plr(thetha_dB, l, gamma=3.0, pure=True, itf_mean=False))
        pure_max_p_f_mrc_div['33'].append(mrc_curve_fitting.lt2plr(thetha_dB, l, gamma=3.3, pure=True, itf_mean=False))
        pure_max_p_f_mrc_div['35'].append(mrc_curve_fitting.lt2plr(thetha_dB, l, gamma=3.5, pure=True, itf_mean=False))
        pure_max_p_f_mrc_div['37'].append(mrc_curve_fitting.lt2plr(thetha_dB, l, gamma=3.7, pure=True, itf_mean=False))
        pure_max_p_f_mrc_div['40'].append(mrc_curve_fitting.lt2plr(thetha_dB, l, gamma=4.0, pure=True, itf_mean=False))
        pure_max_p_f_mrc_div['60'].append(mrc_curve_fitting.lt2plr(thetha_dB, l, gamma=6.0, pure=True, itf_mean=False))

    pure_max_p_f_mrc_div['30'] = np.array(pure_max_p_f_mrc_div['30'])
    pure_max_p_f_mrc_div['33'] = np.array(pure_max_p_f_mrc_div['33'])
    pure_max_p_f_mrc_div['35'] = np.array(pure_max_p_f_mrc_div['35'])
    pure_max_p_f_mrc_div['37'] = np.array(pure_max_p_f_mrc_div['37'])
    pure_max_p_f_mrc_div['40'] = np.array(pure_max_p_f_mrc_div['40'])
    pure_max_p_f_mrc_div['60'] = np.array(pure_max_p_f_mrc_div['60'])


    # Simulation results from Xavier
    sim_intensity = np.array([0.02,	0.04,	0.06,	0.08,	0.1,	0.15,	0.2,	0.3,	0.4,	0.5])

    sim_plr_sc_divers_pure_max ={}
    sim_plr_sc_divers_pure_max['33'] = np.array([0.00022,	0.00568,	0.0289,	    0.0695,	    0.12134,	0.25206,	0.3595,	    0.51378,	0.60664,	0.67144])
    sim_plr_sc_divers_pure_max['35'] = np.array([0.00012,	0.0054,	    0.02654,	0.06574,	0.11288,	0.23792,	0.33928,	0.49264,	0.58266,	0.65498])
    sim_plr_sc_divers_pure_max['37'] = np.array([0.00012,	0.00504,	0.0249,	    0.05996,	0.10294,	0.21926,	0.31906,	0.4689,	    0.5607,	    0.632])
    sim_plr_sc_divers_pure_max['40'] = np.array([0.00012,	0.00436,	0.02366,	0.05458,	0.09734,	0.2042,	0.30282,	    0.44526,	0.53832,	0.6156])

    sim_plr_sc_divers_pure_max_semi_ci = {}
    sim_plr_sc_divers_pure_max_semi_ci['33'] = np.array([0.00013904,	0.00066626,	0.0014544,	0.0025873,	0.0026012,	0.0047677,	0.0036176,	0.0042529,	0.0047698,	0.0040795])
    sim_plr_sc_divers_pure_max_semi_ci['35'] = np.array([0.00010577,	0.0006022,	0.0011296,	0.002036,	0.0022707,	0.0030058,	0.0047599,	0.0043963,	0.0041246,	0.0039331])
    sim_plr_sc_divers_pure_max_semi_ci['37'] = np.array([9.01e-05,	0.00053737,	0.0015001,	0.0016057,	0.0023207,	0.0032936,	0.0037813,	0.0046561,	0.0044045,	0.0046428])
    sim_plr_sc_divers_pure_max_semi_ci['40'] = np.array([0.00010577,	0.00045634,	0.0012661,	0.0020286,	0.002698,	0.003538,	0.0037602,	0.003696,	0.0044676,	0.004157])





    sim_plr_mrc_divers_pure_max ={}

    sim_plr_mrc_divers_pure_max['33'] = np.array([0,	0,	0,	0.00028,	0.0019,	0.02968,	0.10302,	0.30052,	0.4556,	0.56456])
    sim_plr_mrc_divers_pure_max['35'] = np.array([0,	0,	2.00e-05,	0.00036,	0.00304,	0.03504,	0.11296,	0.30564,	0.45162,	0.55778])
    sim_plr_mrc_divers_pure_max['37'] = np.array([0,	0,	0,	0.0007,	    0.0042,	    0.0426,	    0.12232,	0.31026,	0.44458,	0.54676])
    sim_plr_mrc_divers_pure_max['40'] = np.array([0,	0,	0.00012,	0.00098,	0.00566,	0.04852,	0.13198,	0.30608,	0.43796,	0.53782])

    sim_plr_mrc_divers_pure_max_semi_ci ={}

    sim_plr_mrc_divers_pure_max_semi_ci['33'] = np.array([0,	0,	0,	0.00014709,	0.00035605,	0.0015879,	0.0024801,	0.004156,	0.0046076,	0.004349])
    sim_plr_mrc_divers_pure_max_semi_ci['35'] = np.array([0,	0,	3.88e-05,	0.00017319,	0.00045364,	0.0014403,	0.0030202,	0.0045112,	0.0041415,	0.0042623])
    sim_plr_mrc_divers_pure_max_semi_ci['37'] = np.array([0,	0,	0,	0.00022347,	0.00060981,	0.0018545,	0.002809,	0.0036797,	0.0038787,	0.0044185])
    sim_plr_mrc_divers_pure_max_semi_ci['40'] = np.array([0,	0,	9.01e-05,	0.00023184,	0.00059864,	0.0018908,	0.0033682,	0.0034784,	0.0041117,	0.0043899])

    fig, axes = plt.subplots(1, 1, figsize=FIGSIZE, sharey=False)
    axes.set_yscale("log")

    axes.plot(
        L,
        pure_max_p_f_rx_div['35'],
        color='r',  marker='', linestyle='-', linewidth=LINEWIDTH, label=r"SC Diversity,ANA, $\gamma=3.5$"
    )
    # axes.plot(
    #     L,
    #     pure_max_p_f_rx_div['37'],
    #     color='r',  marker='', linestyle='-.', linewidth=LINEWIDTH, label=r"SC Diversity,ANA, $\gamma=3.7$"
    # )
    axes.plot(
        L,
        pure_max_p_f_rx_div['40'],
        color='r',  marker='', linestyle='--', linewidth=LINEWIDTH, label=r"SC Diversity,ANA, $\gamma=4.0$"
    )

    axes.plot(
        L,
        pure_max_p_f_mrc_div['35'],
        color='m',  marker='', linestyle='-', linewidth=LINEWIDTH, label="MRC Diversity,ANA, $\gamma=3.5$"
    )

    # axes.plot(
    #     L,
    #     pure_max_p_f_mrc_div['37'],
    #     color='m',  marker='', linestyle='-.', linewidth=LINEWIDTH, label="MRC Diversity,ANA, $\gamma=3.7$"
    # )
    axes.plot(
        L,
        pure_max_p_f_mrc_div['40'],
        color='m',  marker='', linestyle='--', linewidth=LINEWIDTH, label="MRC Diversity,ANA, $\gamma=4.0$"
    )

    axes.errorbar(
        sim_intensity,
        sim_plr_sc_divers_pure_max['40'],
        yerr=[sim_plr_sc_divers_pure_max_semi_ci['40'],  sim_plr_sc_divers_pure_max_semi_ci['40']],
        fmt='s',
        mfc='none',
        ecolor='r',
        capthick=2,
        label="SC Diversity,SIM,$\gamma=4.0$"
    )

    axes.errorbar(
        sim_intensity,
        sim_plr_sc_divers_pure_max['35'],
        yerr=[sim_plr_sc_divers_pure_max_semi_ci['35'],  sim_plr_sc_divers_pure_max_semi_ci['35']],
        fmt='o',
        mfc='none',
        ecolor='r',
        capthick=2,
        label="SC Diversity,SIM,$\gamma=3.5$"
    )
    # axes.errorbar(
    #     sim_intensity,
    #     sim_plr_sc_divers_pure_max['33'],
    #     yerr=[sim_plr_sc_divers_pure_max_semi_ci['33'],  sim_plr_sc_divers_pure_max_semi_ci['33']],
    #     fmt='*',
    #     mfc='none',
    #     ecolor='r',
    #     capthick=2,
    #     label="SC Diversity,SIM,8dB shadowing, $\gamma=3.3$"
    # )


    axes.errorbar(
        sim_intensity,
        sim_plr_mrc_divers_pure_max['40'],
        yerr=[sim_plr_mrc_divers_pure_max_semi_ci['40'],  sim_plr_mrc_divers_pure_max_semi_ci['40']],
        fmt='v',
        mfc='none',
        ecolor='b',
        capthick=2,
        label=r"MRC Diversity,SIM,$\gamma=4.0$"
    )

    axes.errorbar(
        sim_intensity,
        sim_plr_mrc_divers_pure_max['35'],
        yerr=[sim_plr_mrc_divers_pure_max_semi_ci['35'],  sim_plr_mrc_divers_pure_max_semi_ci['35']],
        fmt='^',
        mfc='none',
        ecolor='b',
        capthick=2,
        label=r"MRC Diversity,SIM,$\gamma=3.5$"
    )




    #===================================================MAYBE DELETE THIS PART =========================================
    # lambda_m = np.linspace(1.0, 3.0, 100)
    # axes.plot(
    #     p*lambda_m/lambda_b,
    #     sgam.mrc_bs_rx_div_op(lambda_m, lambda_b, gamma, p, thetha_dB, 8.0, pure=True, itf_mean=False),
    #     color='k',  marker='', linestyle=':', linewidth=LINEWIDTH, label="MRC Diversity,ANA numerical"
    # )
    #===================================================MAYBE DELETE THIS PART =========================================


    #=====================================Gamma = 3.5 simulation ======================================================

    axes.grid()
    axes.axis([X_START, X_END, Y_START, Y_END])
    axes.set_xlabel(r"Normalized Load")
    axes.set_ylabel("Packet Loss Rate")
    plt.legend(loc='best')
    # plt.legend(bbox_to_anchor=(0.119, 0.02, 0.79, 1), loc=1, ncol=3, mode="expand", bbox_transform=plt.gcf().transFigure)
    plt.savefig(os.path.join(FIG_DST, FIG_NAME), format='eps', dpi=300)

    plt.show()

