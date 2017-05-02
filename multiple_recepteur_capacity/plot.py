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

from analytical_model import sgam

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
    X_END = 0.03
    X_STEP = 0.002
    Y_START = 1e-3
    Y_END = 0.6
    Y_STEP = 0.1

    MAX_TRANS = 1
    LINEWIDTH = 3
    #logs_dir_1 = glob.glob(os.path.join("..",  LOG_DIR, SUB_DIR, CASE_DIR, SUB_CASE_DIR, "*.csv"))
    bs_rx_div_sigma_0 = glob.glob(os.path.join(
        "..",  LOG_DIR, SUB_DIR, "fading_shadowing", "BS_RX_DIVERS", "p_0.04_bs_0.004_R_40",
        "sigma_0dB", "*.csv")
    )

    bs_rx_div_sigma_8 = glob.glob(os.path.join(
        "..",  LOG_DIR, SUB_DIR, "fading_shadowing", "BS_RX_DIVERS", "p_0.04_bs_0.004_R_40",
        "sigma_8dB", "*.csv")
    )

    bs_nst_att_sigma_8 = glob.glob(os.path.join(
        "..",  LOG_DIR, SUB_DIR, "fading_shadowing", "BS_NST_ATT", "p_0.01_bs_0.05_R_25",
        "sigma_8dB_50per", "*.csv")
    )

    bs_nst_att_sigma_0 = glob.glob(os.path.join(
        "..",  LOG_DIR, SUB_DIR, "fading_shadowing", "BS_NST_ATT", "p_0.002_bs_0.01_R_25",
        "sigma_0dB_50per", "*.csv")
    )


    sim_intensity_0dB, sim_plr_list_0dB = sgam.sim_data_process(bs_rx_div_sigma_0)
    sim_intensity_8dB, sim_plr_list_8dB = sgam.sim_data_process(bs_rx_div_sigma_8)
    sim_intensity_nst_8dB, sim_plr_list_nst_8dB = sgam.sim_data_process(bs_nst_att_sigma_8)
    sim_intensity_nst_0dB, sim_plr_list_nst_0dB = sgam.sim_data_process(bs_nst_att_sigma_0)

    # sim_intensity_4, sim_plr_list_4 = sim_data_process(logs_dir_4)
    # sim_intensity_5, sim_plr_list_5 = sim_data_process(bs_nst_att_8dB)

    sim_plr_0 = [element[1] for element in sim_plr_list_0dB]
    sim_plr_lower_0 = [element[1]-element[0] for element in sim_plr_list_0dB]
    sim_plr_upper_0 = [element[2]-element[1] for element in sim_plr_list_0dB]


    sim_plr_8 = [element[1] for element in sim_plr_list_8dB]
    sim_plr_lower_8 = [element[1]-element[0] for element in sim_plr_list_8dB]
    sim_plr_upper_8 = [element[2]-element[1] for element in sim_plr_list_8dB]


    sim_plr_nst_8 = [element[1] for element in sim_plr_list_nst_8dB]
    sim_plr_lower_nst_8 = [element[1]-element[0] for element in sim_plr_list_nst_8dB]
    sim_plr_upper_nst_8 = [element[2]-element[1] for element in sim_plr_list_nst_8dB]
    #
    sim_plr_nst_0 = [element[1] for element in sim_plr_list_nst_0dB]
    sim_plr_lower_nst_0 = [element[1]-element[0] for element in sim_plr_list_nst_0dB]
    sim_plr_upper_nst_0 = [element[2]-element[1] for element in sim_plr_list_nst_0dB]
    #
    # sim_plr_5 = [element[1] for element in sim_plr_list_5]
    # sim_plr_lower_5 = [element[1]-element[0] for element in sim_plr_list_5]
    # sim_plr_upper_5 = [element[2]-element[1] for element in sim_plr_list_5]

    # 生成 lambda_m 的 ndarray
    lambda_m = np.linspace(0, X_END, 1000)
    # 原则上我们让 基站的密度是个肯定值，毕竟这个东西投资大，没必要变化 lambda_b
    lambda_b = 0.004
    # gamma, path loss component. Usually take 4 as value.
    gamma = 4.0
    p = 0.04
    thetha_dB = 3.0 # unit dB

    X_END = p*max(lambda_m)/lambda_b


    # Define p_f_2 as the outage probability over infinite plane
    p_f_rx_div_0 = sgam.bs_rx_div_op(lambda_m, lambda_b, gamma, p, thetha_dB, 0)
    p_f_bs_nst_att_0 = sgam.bs_nearest_atch_op(lambda_m, lambda_b, gamma, p, thetha_dB, 0)
    p_f_bs_nst_att_8 = sgam.bs_nearest_atch_op(lambda_m, lambda_b, gamma, p, thetha_dB, 8)
    p_f_bs_bst_att_8 = sgam.bs_best_atch_op(lambda_m, lambda_b, gamma, p, thetha_dB, 8)

    fig, axes = plt.subplots(1, 1, figsize=FIGSIZE, sharey=False)

    axes.set_yscale("log")
    print zip(p*lambda_m/lambda_b, p_f_bs_nst_att_8)
    axes.plot(
        p*lambda_m/lambda_b,
        p_f_rx_div_0,
        color='r',  marker='', linestyle='-.', linewidth=LINEWIDTH, label="Diversity,ANA"
    )
    axes.plot(
        p*lambda_m/lambda_b,
        p_f_bs_nst_att_8,
        color='g',  marker='', linestyle='-', linewidth=LINEWIDTH, label="Nearest,ANA,8dB shadowing"
    )

    axes.plot(
        p*lambda_m/lambda_b,
        p_f_bs_bst_att_8,
        color='k',  marker='', linestyle='-', linewidth=LINEWIDTH, label="Best,ANA,8dB shadowing"
    )

    axes.plot(
        p*lambda_m/lambda_b,
        p_f_bs_nst_att_0,
        color='b',  marker='', linestyle='--', linewidth=LINEWIDTH, label="Nearest,ANA,no shadowing"
    )


    axes.errorbar(
        p*sim_intensity_0dB/lambda_b,
        sim_plr_0,
        yerr=[sim_plr_lower_0, sim_plr_upper_0],
        fmt='*',
        mfc='r',
        ecolor='r',
        capthick=2,
        label="Diversity SIM,no shadowing"
    )
    axes.errorbar(
        p*sim_intensity_8dB/lambda_b,
        sim_plr_8,
        yerr=[sim_plr_lower_8, sim_plr_upper_8],
        fmt='s',
        mfc='r',
        ecolor='k',
        capthick=2,
        label="Diversity SIM,8dB shadowing"
    )
    axes.errorbar(
        sim_intensity_nst_0dB/5.0,
        sim_plr_nst_0,
        yerr=[sim_plr_lower_nst_0, sim_plr_upper_nst_0],
        fmt='o',
        ecolor='b',
        mfc='b',
        capthick=2,
        label="Nearest SIM,no shadowing"
    )
    axes.errorbar(
        sim_intensity_nst_8dB/5.0,
        sim_plr_nst_8,
        yerr=[sim_plr_lower_nst_8, sim_plr_upper_nst_8],
        fmt='d',
        mfc='g',
        ecolor='g',
        capthick=2,
        label="Nearest SIM,8dB shadowing"
    )

    axes.grid()
    # axes.set_xticks(np.arange(X_START, X_END, X_STEP))
    axes.set_xlabel(r"Normalized Load")
    axes.set_ylabel(r"Packet loss rate")
    axes.axis([X_START, X_END, Y_START, Y_END])


    # plt.legend(bbox_to_anchor=(0.119, 0.01, 0.79, 1), loc=1, ncol=3, mode="expand", bbox_transform=plt.gcf().transFigure)
    plt.legend(loc='best')
    plt.savefig('new_packet_loss_rate_mpr.eps', format='eps', dpi=300)

    plt.show()
