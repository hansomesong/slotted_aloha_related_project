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

LOG_DIR = 'logs'
SUB_DIR = 'multiple_reception'
CASE_DIR = 'fading'
SUB_CASE_DIR = "bs_0.01_p_0.002_R_40|100"
DATA_FOLDED = '.'
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
    thetha_dB = 7.0 # unit dB
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
    FIG_DST = "/Users/qsong/Documents/Communication_Letter_MRC"
    FIG_NAME = 'pure_slot_packet_loss_rate_mpr_gamma_{0}_theta_{1}.eps'.format(3, thetha_dB)

    # Analytical results for SC macro diversity
    pure_p_f_rx_div_8 = sgam.bs_rx_div_op(lambda_m, lambda_b, gamma, p, thetha_dB, 8, True)

        # Analytical results for SC macro diversity
    pure_mean_p_f_rx_div ={}
    pure_mean_p_f_rx_div['30'] = sgam.bs_rx_div_op(lambda_m, lambda_b, 3.0, p, thetha_dB, 8, True, True)
    pure_mean_p_f_rx_div['33'] = sgam.bs_rx_div_op(lambda_m, lambda_b, 3.3, p, thetha_dB, 8, True, True)
    pure_mean_p_f_rx_div['35'] = sgam.bs_rx_div_op(lambda_m, lambda_b, 3.5, p, thetha_dB, 8, True, True)
    pure_mean_p_f_rx_div['37'] = sgam.bs_rx_div_op(lambda_m, lambda_b, 3.7, p, thetha_dB, 8, True, True)
    pure_mean_p_f_rx_div['40'] = sgam.bs_rx_div_op(lambda_m, lambda_b, 4.0, p, thetha_dB, 8, True, True)
    pure_mean_p_f_rx_div['60'] = sgam.bs_rx_div_op(lambda_m, lambda_b, 6.0, p, thetha_dB, 8, True, True)

    # Analytical results for MRC macro diversity
    pure_mean_p_f_mrc_div = {}
    pure_mean_p_f_mrc_div['30'] = []
    pure_mean_p_f_mrc_div['33'] = []
    pure_mean_p_f_mrc_div['35'] = []
    pure_mean_p_f_mrc_div['37'] = []
    pure_mean_p_f_mrc_div['40'] = []
    pure_mean_p_f_mrc_div['60'] = []
    for l in L:
        pure_mean_p_f_mrc_div['30'].append(mrc_curve_fitting.lt2plr(thetha_dB, l, gamma=3.0, pure=True, itf_mean=True))
        pure_mean_p_f_mrc_div['33'].append(mrc_curve_fitting.lt2plr(thetha_dB, l, gamma=3.3, pure=True, itf_mean=True))
        pure_mean_p_f_mrc_div['35'].append(mrc_curve_fitting.lt2plr(thetha_dB, l, gamma=3.5, pure=True, itf_mean=True))
        pure_mean_p_f_mrc_div['37'].append(mrc_curve_fitting.lt2plr(thetha_dB, l, gamma=3.7, pure=True, itf_mean=True))
        pure_mean_p_f_mrc_div['40'].append(mrc_curve_fitting.lt2plr(thetha_dB, l, gamma=4.0, pure=True, itf_mean=True))
        pure_mean_p_f_mrc_div['60'].append(mrc_curve_fitting.lt2plr(thetha_dB, l, gamma=6.0, pure=True, itf_mean=True))

    pure_mean_p_f_mrc_div['30'] = np.array(pure_mean_p_f_mrc_div['30'])
    pure_mean_p_f_mrc_div['33'] = np.array(pure_mean_p_f_mrc_div['33'])
    pure_mean_p_f_mrc_div['35'] = np.array(pure_mean_p_f_mrc_div['35'])
    pure_mean_p_f_mrc_div['37'] = np.array(pure_mean_p_f_mrc_div['37'])
    pure_mean_p_f_mrc_div['40'] = np.array(pure_mean_p_f_mrc_div['40'])
    pure_mean_p_f_mrc_div['60'] = np.array(pure_mean_p_f_mrc_div['60'])


    # Simulation results from Xavier
    sim_intensity = np.array([0.02,	0.04,	0.06,	0.08,	0.1,	0.15,	0.2,	0.3,	0.4,	0.5])


    sim_plr_sc_divers_pure_mean ={}

    sim_plr_sc_divers_pure_mean['33'] = np.array([0.0001,	0.00098,	0.00756,	0.02356,	0.05076,	0.14212,    0.2389,	    0.39282,	0.5013,	    0.57616])
    sim_plr_sc_divers_pure_mean['35'] = np.array([0.00044,	0.01102,	0.04534,	0.10048,	0.164,	    0.29942,	0.41212,	0.55842,	0.64626,	0.70212])
    sim_plr_sc_divers_pure_mean['37'] = np.array([0,	    0.00092,	0.00654,	0.02166,	0.04422,	0.1251,	    0.21082,    0.35776,	0.45316,	0.53576])
    sim_plr_sc_divers_pure_mean['40'] = np.array( [0.00036,0.0095,0.03786,0.07948,0.13594,0.25344,0.35672,0.49996,0.59108,0.65616])
    sim_plr_sc_divers_pure_mean['60'] = np.array([0,	0.00014,	0.00264,	0.01198,	0.02954,	0.1027,	 0.17998,	0.3263,	0.43388,	0.51678])




    sim_plr_sc_divers_pure_mean_semi_ci ={}
    sim_plr_sc_divers_pure_mean_semi_ci['33'] = np.array([9.99e-05,	0.00033488,	0.00084286,	0.0013086,	0.0018366,	0.0038045,	0.0036188,	0.0039937,	0.0046482,	0.0040863])
    sim_plr_sc_divers_pure_mean_semi_ci['35'] = np.array([0.0001767,	0.0009416,	0.0016301,	0.0030147,	0.0040222,	0.0037803,	0.0042595,	0.0038252,	0.0035824,	0.0041142])
    sim_plr_sc_divers_pure_mean_semi_ci['37'] = np.array([0,	    0.00025308,	0.0006817,	0.001209,	0.0017366,	0.0026157,	0.0033865,	0.00393,	0.0038966,	0.0044171])
    sim_plr_sc_divers_pure_mean_semi_ci['40'] = np.array([0,	    0.00021959,	0.00074213,	0.0010926,	0.0017278,	0.0026366,	0.0032863,	0.0041037,	0.0040214,	0.0043918])
    sim_plr_sc_divers_pure_mean_semi_ci['60'] = np.array([0,	    9.62e-05,	0.00054249,	0.00074478,	0.001625,	0.0026157,	0.0036363,	0.0036231,	0.005317,	0.0038629])




    sim_plr_mrc_divers_pure_mean ={}
    sim_plr_mrc_divers_pure_mean['33'] = np.array([0,	0,	0,	2.00e-05,	6.00e-05,	0.00282,	0.02198,	0.13558,	0.28186,	0.40376])
    sim_plr_mrc_divers_pure_mean['35'] = np.array([0,	0,	0.00014,	0.0022,	0.01024,	0.07658,	0.19246,	0.40558,	0.54364,	0.63082])
    sim_plr_mrc_divers_pure_mean['37'] = np.array([0,	0,	0,	0,	0.00018,	0.00684,	0.03562,	0.16374,	0.29302,	0.40534])
    sim_plr_mrc_divers_pure_mean['40'] = np.array([0,	0,	0.0005,	0.0037,	0.01562,	0.08542,	0.19438,	0.38124,	0.50978,	0.59434])
    sim_plr_mrc_divers_pure_mean['60'] = np.array([0,	0,	0.00014,	0.00158,	0.00538,	0.03572,	0.08954,	0.2081,	0.30984,	0.39702])



    sim_plr_mrc_divers_pure_mean_semi_ci ={}
    sim_plr_mrc_divers_pure_mean_semi_ci['33'] = np.array([0,	0,	0,	        3.88e-05,	6.58e-05,	0.00046939,	0.0011714,	0.002595,	0.0040964,	0.0043537])
    sim_plr_mrc_divers_pure_mean_semi_ci['35'] = np.array([0,	0,	9.62e-05,	0.00032797,	0.0010305,	0.0023245,	0.0037094,	0.0046027,	0.0042513,	0.0036675])
    sim_plr_mrc_divers_pure_mean_semi_ci['37'] = np.array([0,	0,	0,	        0,	        0.00012006,	0.00062563,	0.0015008,	0.0029574,	0.0040819,	0.0038141])
    sim_plr_mrc_divers_pure_mean_semi_ci['40'] = np.array([0.00018185,	0.00082646,	0.0016924,	0.0023801,	0.0028351,	0.003261,	0.004622,	0.0046263,	0.0039582,	0.0052032])
    sim_plr_mrc_divers_pure_mean_semi_ci['60'] = np.array([0,	0,	9.62e-05,	0.00037636,	0.00059162,	0.0015885,	0.0024185,	0.0034905,	0.0039816,	0.003827])



    fig, axes = plt.subplots(1, 1, figsize=FIGSIZE, sharey=False)
    axes.set_yscale("log")


    axes.plot(
        L,
        pure_mean_p_f_rx_div['35'],
        color='r',  marker='', linestyle='-', linewidth=LINEWIDTH, label=r"SC Diversity,ANA, $\gamma=3.5$"
    )
    axes.plot(
        L,
        pure_mean_p_f_rx_div['60'],
        color='r',  marker='', linestyle='-.', linewidth=LINEWIDTH, label=r"SC Diversity,ANA, $\gamma=6.0$"
    )
    axes.plot(
        L,
        pure_mean_p_f_rx_div['40'],
        color='r',  marker='', linestyle='--', linewidth=LINEWIDTH, label=r"SC Diversity,ANA, $\gamma=4.0$"
    )

    axes.plot(
        L,
        pure_mean_p_f_mrc_div['35'],
        color='m',  marker='', linestyle='-', linewidth=LINEWIDTH, label="MRC Diversity,ANA, $\gamma=3.5$"
    )

    axes.plot(
        L,
        pure_mean_p_f_mrc_div['60'],
        color='m',  marker='', linestyle='-.', linewidth=LINEWIDTH, label="MRC Diversity,ANA, $\gamma=6.0$"
    )

    axes.plot(
        L,
        pure_mean_p_f_mrc_div['40'],
        color='m',  marker='', linestyle='--', linewidth=LINEWIDTH, label="MRC Diversity,ANA, $\gamma=4.0$"
    )

    axes.errorbar(
        sim_intensity,
        sim_plr_sc_divers_pure_mean['60'],
        yerr=[sim_plr_sc_divers_pure_mean_semi_ci['60'],  sim_plr_sc_divers_pure_mean_semi_ci['60']],
        fmt='*',
        mfc='none',
        ecolor='r',
        capthick=2,
        label="SC Diversity,SIM,$\gamma=6.0$"
    )

    axes.errorbar(
        sim_intensity,
        sim_plr_sc_divers_pure_mean['40'],
        yerr=[sim_plr_sc_divers_pure_mean_semi_ci['40'],  sim_plr_sc_divers_pure_mean_semi_ci['40']],
        fmt='s',
        mfc='none',
        ecolor='r',
        capthick=2,
        label="SC Diversity,SIM,$\gamma=4.0$"
    )

    axes.errorbar(
        sim_intensity,
        sim_plr_sc_divers_pure_mean['35'],
        yerr=[sim_plr_sc_divers_pure_mean_semi_ci['35'],  sim_plr_sc_divers_pure_mean_semi_ci['35']],
        fmt='o',
        mfc='none',
        ecolor='r',
        capthick=2,
        label="SC Diversity,SIM,$\gamma=3.5$"
    )
    # axes.errorbar(
    #     sim_intensity,
    #     sim_plr_sc_divers_pure_mean['33'],
    #     yerr=[sim_plr_sc_divers_pure_mean_semi_ci['33'],  sim_plr_sc_divers_pure_mean_semi_ci['33']],
    #     fmt='*',
    #     mfc='none',
    #     ecolor='r',
    #     capthick=2,
    #     label="SC Diversity,SIM,8dB shadowing, $\gamma=3.3$"
    # )


    axes.errorbar(
        sim_intensity,
        sim_plr_mrc_divers_pure_mean['40'],
        yerr=[sim_plr_mrc_divers_pure_mean_semi_ci['40'],  sim_plr_mrc_divers_pure_mean_semi_ci['40']],
        fmt='v',
        mfc='none',
        ecolor='b',
        capthick=2,
        label=r"MRC Diversity,SIM,$\gamma=4.0$"
    )

    axes.errorbar(
        sim_intensity,
        sim_plr_mrc_divers_pure_mean['35'],
        yerr=[sim_plr_mrc_divers_pure_mean_semi_ci['35'],  sim_plr_mrc_divers_pure_mean_semi_ci['35']],
        fmt='^',
        mfc='none',
        ecolor='b',
        capthick=2,
        label=r"MRC Diversity,SIM,$\gamma=3.5$"
    )

    axes.errorbar(
        sim_intensity,
        sim_plr_mrc_divers_pure_mean['60'],
        yerr=[sim_plr_mrc_divers_pure_mean_semi_ci['60'],  sim_plr_mrc_divers_pure_mean_semi_ci['60']],
        fmt='>',
        mfc='none',
        ecolor='b',
        capthick=2,
        label=r"MRC Diversity,SIM,$\gamma=6.0$"
    )


    axes.grid()
    axes.axis([X_START, X_END, Y_START, Y_END])
    axes.set_xlabel(r"Normalized Load")
    axes.set_ylabel("Packet Loss Rate")
    #
    # # #===================================================MAYBE DELETE THIS PART =========================================
    # lambda_m = np.linspace(1.0, 5.0, 100)
    # axes.plot(
    #     p*lambda_m/lambda_b,
    #     # sgam.mrc_bs_rx_div_op(lambda_m, lambda_b, 3.0, p, thetha_dB, 8.0, pure=True, itf_mean=True),
    #     [sgam.lt2plr(thetha_dB, x, lambda_b, 3.0, p, 8.0, pure=True, itf_mean=True) for x in lambda_m],
    #     color='k',  marker='', linestyle=':', linewidth=LINEWIDTH, label="MRC Diversity,ANA numerical, gamma=3"
    # )
    #
    # plr_gamma_3 = [0, 0.000925, 0.01235, 0.12442, 0.28232, 0.41652]
    # L = np.array([25,	37.5, 50, 	75, 100,	125])/250.0
    #
    # plt.scatter(L, plr_gamma_3, label="simulation, gamma=3.0")
    # axes.plot(
    #     p*lambda_m/lambda_b,
    #     # sgam.mrc_bs_rx_div_op(lambda_m, lambda_b, 4.0, p, thetha_dB, 8.0, pure=True, itf_mean=True),
    #     [sgam.lt2plr(thetha_dB, x, lambda_b, 3.5, p, 8.0, pure=True, itf_mean=True) for x in lambda_m],
    #     color='y',  marker='', linestyle=':', linewidth=LINEWIDTH, label="MRC Diversity,ANA numerical, gamma=3.5"
    # )
    #
    # plr_gamma_35 = np.array([5e-005,	0.0038, 	0.0295, 	0.15515, 	0.29825, 	0.41218])
    # L_gamma35 = np.array([25,	37.5,	50,	 75,	100,	125])/250.0
    #
    # plt.scatter(L_gamma35, plr_gamma_35, label="simulation, gamma=3.5")
    #
    #
    #
    # axes.plot(
    #     p*lambda_m/lambda_b,
    #     # sgam.mrc_bs_rx_div_op(lambda_m, lambda_b, 5.0, p, thetha_dB, 8.0, pure=True, itf_mean=True),
    #     [sgam.lt2plr(thetha_dB, x, lambda_b, 6.0, p, 8.0, pure=True, itf_mean=True) for x in lambda_m],
    #     color='b',  marker='', linestyle=':', linewidth=LINEWIDTH, label="MRC Diversity,ANA numerical, gamma=6"
    # )
    #
    # plr_gamma_6 = np.array([0.00705, 0.088075,	0.20797,	0.31498])
    #
    # L_gamma6 = np.array([25,	50,	 75, 100])/250.0
    #
    # plt.scatter(L_gamma6, plr_gamma_6, label="simulation, gamma=6")
    # # #===================================================MAYBE DELETE THIS PART =========================================

    # print zip(p*lambda_m/lambda_b, p_f_rx_div_0)

    # plt.legend(bbox_to_anchor=(0.119, 0.01, 0.79, 1), loc=1, ncol=3, mode="expand", bbox_transform=plt.gcf().transFigure)
    plt.legend(loc='best')
    # plt.legend(bbox_to_anchor=(0.119, 0.02, 0.79, 1), loc=1, ncol=3, mode="expand", bbox_transform=plt.gcf().transFigure)
    plt.savefig(os.path.join(FIG_DST, FIG_NAME), format='eps', dpi=300)

    plt.show()

