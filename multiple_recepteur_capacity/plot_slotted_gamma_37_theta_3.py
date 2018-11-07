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
from scipy.special import erf as erf
from analytical_model import sgam, mrc_curve_fitting


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

A_P_DECREMENT ="analytical, power decrement"
A_P_IDENTIC = "analytical, identical power"
A_P_INCREMENT = "analytical, power increment"
S_P_IDENTIC = "simulation, identical power"
S_P_INCREMENT = "simulation, power increment"
S_P_DECREMENT = "simulation, power decrement"



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
    lambda_m = np.linspace(.5, X_END, 100)
    # 原则上我们让 基站的密度是个肯定值，毕竟这个东西投资大，没必要变化 lambda_b
    lambda_b = 0.08
    # gamma, path loss component. Usually take 4 as value.

    p = 0.008
    L = p*lambda_m / lambda_b


    X_END = p*max(lambda_m)/lambda_b


    gamma = 3.7
    thetha_dB = 3.0 # unit dB
    theta = np.power(10, thetha_dB/10)

    FIG_DST = "/Users/qsong/Documents/Communication_Letter_MRC"
    FIG_NAME = 'slot_packet_loss_rate_mpr_gamma_{0}_theta_{1}.eps'.format(3, thetha_dB)


     # Analytical results for SC macro diversity
    slot_p_f_rx_div_8 = sgam.bs_rx_div_op(lambda_m, lambda_b, gamma, p, thetha_dB, 8, False)
    # Analytical results for MRC macro diversity
    slot_p_f_mrc_div = []
    for l in L:
        slot_p_f_mrc_div.append(mrc_curve_fitting.lt2plr(thetha_dB, l, gamma, pure=False, itf_mean=True))
    slot_p_f_mrc_div = np.array(slot_p_f_mrc_div)


    slot_p_f_rx_div ={}
    slot_p_f_rx_div['30'] = sgam.bs_rx_div_op(lambda_m, lambda_b, 3.0, p, thetha_dB, 8, False, True)
    slot_p_f_rx_div['33'] = sgam.bs_rx_div_op(lambda_m, lambda_b, 3.3, p, thetha_dB, 8, False, True)
    slot_p_f_rx_div['35'] = sgam.bs_rx_div_op(lambda_m, lambda_b, 3.5, p, thetha_dB, 8, False, True)
    slot_p_f_rx_div['37'] = sgam.bs_rx_div_op(lambda_m, lambda_b, 3.7, p, thetha_dB, 8, False, True)
    slot_p_f_rx_div['40'] = sgam.bs_rx_div_op(lambda_m, lambda_b, 4.0, p, thetha_dB, 8, False, True)
    slot_p_f_rx_div['60'] = sgam.bs_rx_div_op(lambda_m, lambda_b, 6.0, p, thetha_dB, 8, False, True)

    # Analytical results for MRC macro diversity
    slot_p_f_mrc_div = {}
    slot_p_f_mrc_div['30'] = []
    slot_p_f_mrc_div['33'] = []
    slot_p_f_mrc_div['35'] = []
    slot_p_f_mrc_div['37'] = []
    slot_p_f_mrc_div['40'] = []
    slot_p_f_mrc_div['60'] = []
    for l in L:
        slot_p_f_mrc_div['30'].append(mrc_curve_fitting.lt2plr(thetha_dB, l, gamma=3.0, pure=False, itf_mean=True))
        slot_p_f_mrc_div['33'].append(mrc_curve_fitting.lt2plr(thetha_dB, l, gamma=3.3, pure=False, itf_mean=True))
        slot_p_f_mrc_div['35'].append(mrc_curve_fitting.lt2plr(thetha_dB, l, gamma=3.5, pure=False, itf_mean=True))
        slot_p_f_mrc_div['37'].append(mrc_curve_fitting.lt2plr(thetha_dB, l, gamma=3.7, pure=False, itf_mean=True))
        slot_p_f_mrc_div['40'].append(mrc_curve_fitting.lt2plr(thetha_dB, l, gamma=4.0, pure=False, itf_mean=True))
        slot_p_f_mrc_div['60'].append(mrc_curve_fitting.lt2plr(thetha_dB, l, gamma=6.0, pure=False, itf_mean=True))

    slot_p_f_mrc_div['30'] = np.array(slot_p_f_mrc_div['30'])
    slot_p_f_mrc_div['33'] = np.array(slot_p_f_mrc_div['33'])
    slot_p_f_mrc_div['35'] = np.array(slot_p_f_mrc_div['35'])
    slot_p_f_mrc_div['37'] = np.array(slot_p_f_mrc_div['37'])
    slot_p_f_mrc_div['40'] = np.array(slot_p_f_mrc_div['40'])
    slot_p_f_mrc_div['60'] = np.array(slot_p_f_mrc_div['60'])


    # Simulation results from Xavier
    sim_intensity = np.array([0.02,	0.04,	0.06,	0.08,	0.1,	0.15,	0.2,	0.3,	0.4,	0.5])
    sim_plr_sc_divers_slot ={}
    sim_plr_sc_divers_slot['35'] = np.array([0,	0.00034,	0.00192,	0.00862,	0.02042,	0.0721,	    0.14202,	0.27964,	0.3832,	    0.47192])
    sim_plr_sc_divers_slot['37'] = np.array([0,	0.00034,	0.00212,	0.0075,	    0.01836,	0.06554,	0.13026,	0.25732,	0.35842,	0.44508])
    sim_plr_sc_divers_slot['40'] = np.array([0,	0.00022,	0.00154,	0.00586,	0.01602,	0.05926,	0.11746,	0.23758,	0.33132,	0.41822])


    sim_plr_sc_divers_slot_semi_ci ={}
    sim_plr_sc_divers_slot_semi_ci['35'] = np.array([0,	0.00016267,	0.00041053,	0.00095842,	0.0010907,	0.0025557,	0.0030719,	0.0041154,	0.0043109,	0.0046833])
    sim_plr_sc_divers_slot_semi_ci['37'] = np.array([0,	0.00017186,	0.00036622,	0.00071048,	0.0010486,	0.0021721,	0.0028462,	0.0034212,	0.0035966,	0.0045431])
    sim_plr_sc_divers_slot_semi_ci['40'] = np.array([0,	0.00014968,	0.00034268,	0.00071643,	0.0010259,	0.0020778,	0.0032544,	0.0042423,	0.0040931,	0.0040349])




    sim_plr_mrc_divers_slot ={}
    sim_plr_mrc_divers_slot['35'] = np.array([0,	0,	0,	0,	0,	0.00094,	0.0061,	 0.05834,	0.16264,	0.27686])
    sim_plr_mrc_divers_slot['37'] = np.array([0,	0,	0,	0,	0,	0.00078,	0.00852,	0.06576,	0.16552,	0.27324])
    sim_plr_mrc_divers_slot['40'] = np.array([0,	0,	0,	0,	2.00E-05,	0.00158,	0.00996,	0.07194,	0.1682,	0.2691])

    sim_plr_mrc_divers_slot_semi_ci ={}
    sim_plr_mrc_divers_slot_semi_ci['35'] = np.array([0,	0,	0,	0,	0,	0.00026244,	0.00076264,	0.0022645,	0.0035379,	0.004126])
    sim_plr_mrc_divers_slot_semi_ci['37'] = np.array([0,	0,	0,	0,	0,	0.00023709,	0.00074633,	0.0021731,	0.0032274,	0.0033458])
    sim_plr_mrc_divers_slot_semi_ci['40'] = np.array([0,	0,	0,	0,	3.88e-05,	0.0003041,	0.0010503,	0.0018966,	0.0033876,	0.0035224])



    fig, axes = plt.subplots(1, 1, figsize=FIGSIZE, sharey=False)

    axes.set_yscale("log")

    axes.plot(
        L,
        slot_p_f_rx_div['35'],
        color='r',  marker='', linestyle='-', linewidth=LINEWIDTH, label=r"SC Diversity,ANA, $\gamma=3.5$"
    )
    # axes.plot(
    #     L,
    #     slot_p_f_rx_div['60'],
    #     color='r',  marker='', linestyle='-.', linewidth=LINEWIDTH, label=r"SC Diversity,ANA, $\gamma=6.0$"
    # )
    axes.plot(
        L,
        slot_p_f_rx_div['40'],
        color='r',  marker='', linestyle='--', linewidth=LINEWIDTH, label=r"SC Diversity,ANA, $\gamma=4.0$"
    )

    axes.plot(
        L,
        slot_p_f_mrc_div['35'],
        color='m',  marker='', linestyle='-', linewidth=LINEWIDTH, label="MRC Diversity,ANA, $\gamma=3.5$"
    )

    # axes.plot(
    #     L,
    #     slot_p_f_mrc_div['60'],
    #     color='m',  marker='', linestyle='-.', linewidth=LINEWIDTH, label="MRC Diversity,ANA, $\gamma=6.0$"
    # )

    axes.plot(
        L,
        slot_p_f_mrc_div['40'],
        color='m',  marker='', linestyle='--', linewidth=LINEWIDTH, label="MRC Diversity,ANA, $\gamma=4.0$"
    )

    # axes.errorbar(
    #     sim_intensity,
    #     sim_plr_sc_divers_slot['60'],
    #     yerr=[sim_plr_sc_divers_slot_semi_ci['60'],  sim_plr_sc_divers_slot_semi_ci['60']],
    #     fmt='*',
    #     mfc='none',
    #     ecolor='r',
    #     capthick=2,
    #     label="SC Diversity,SIM,$\gamma=6.0$"
    # )

    axes.errorbar(
        sim_intensity,
        sim_plr_sc_divers_slot['40'],
        yerr=[sim_plr_sc_divers_slot_semi_ci['40'],  sim_plr_sc_divers_slot_semi_ci['40']],
        fmt='s',
        mfc='none',
        ecolor='r',
        capthick=2,
        label="SC Diversity,SIM,$\gamma=4.0$"
    )

    axes.errorbar(
        sim_intensity,
        sim_plr_sc_divers_slot['35'],
        yerr=[sim_plr_sc_divers_slot_semi_ci['35'],  sim_plr_sc_divers_slot_semi_ci['35']],
        fmt='o',
        mfc='none',
        ecolor='r',
        capthick=2,
        label="SC Diversity,SIM,$\gamma=3.5$"
    )

    axes.errorbar(
        sim_intensity,
        sim_plr_mrc_divers_slot['40'],
        yerr=[sim_plr_mrc_divers_slot_semi_ci['40'],  sim_plr_mrc_divers_slot_semi_ci['40']],
        fmt='v',
        mfc='none',
        ecolor='b',
        capthick=2,
        label=r"MRC Diversity,SIM,$\gamma=4.0$"
    )

    axes.errorbar(
        sim_intensity,
        sim_plr_mrc_divers_slot['35'],
        yerr=[sim_plr_mrc_divers_slot_semi_ci['35'],  sim_plr_mrc_divers_slot_semi_ci['35']],
        fmt='^',
        mfc='none',
        ecolor='b',
        capthick=2,
        label=r"MRC Diversity,SIM,$\gamma=3.5$"
    )

    # axes.errorbar(
    #     sim_intensity,
    #     sim_plr_mrc_divers_slot['60'],
    #     yerr=[sim_plr_mrc_divers_slot_semi_ci['60'],  sim_plr_mrc_divers_slot_semi_ci['60']],
    #     fmt='>',
    #     mfc='none',
    #     ecolor='b',
    #     capthick=2,
    #     label=r"MRC Diversity,SIM,$\gamma=6.0$"
    # )



    axes.grid()
    # axes.set_xticks(np.arange(X_START, X_END, X_STEP))
    axes.set_xlabel(r"Normalized Load")
    axes.set_ylabel(r"Packet loss rate")
    axes.axis([X_START, X_END, Y_START, Y_END])


    # plt.legend(bbox_to_anchor=(0.119, 0.01, 0.79, 1), loc=1, ncol=3, mode="expand", bbox_transform=plt.gcf().transFigure)
    plt.legend(loc='best')
    plt.savefig(os.path.join(FIG_DST, FIG_NAME), format='eps', dpi=300)

    plt.show()
