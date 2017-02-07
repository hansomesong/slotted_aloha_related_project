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


params = {
    'legend.fontsize': 15,
    'xtick.labelsize': 20,
    'ytick.labelsize': 20,
    'axes.labelsize': 20,
}
plt.rcParams.update(params)

LOG_DIR = 'logs'
SUB_DIR = 'multiple_reception'
CASE_DIR = 'fading'
SUB_CASE_DIR = "bs_0.01_p_0.002_R_40|100"
DATA_FOLDED = '.'
FIGSIZE = (15, 8)

A_P_DECREMENT ="analytical, power decrement"
A_P_IDENTIC = "analytical, identical power"
A_P_INCREMENT = "analytical, power increment"
S_P_IDENTIC = "simulation, identical power"
S_P_INCREMENT = "simulation, power increment"
S_P_DECREMENT = "simulation, power decrement"

X_START = 0.0
X_END = 0.03
X_STEP = 0.002
Y_START = 1e-3
Y_END = 0.6
Y_STEP = 0.1

MAX_TRANS = 1
LINEWIDTH = 2

SCALE = ["log", "linear"]









if __name__ == '__main__':
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
        "..",  LOG_DIR, SUB_DIR, "fading_shadowing", "BS_NST_ATT", "p_0.04_bs_0.004_R_40",
        "sigma_8dB_50per", "*.csv")
    )

    bs_nst_att_sigma_0 = glob.glob(os.path.join(
        "..",  LOG_DIR, SUB_DIR, "fading_shadowing", "BS_NST_ATT", "p_0.04_bs_0.004_R_40",
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
    lambda_m = np.linspace(0, X_END, 100)
    # 原则上我们让 基站的密度是个肯定值，毕竟这个东西投资大，没必要变化 lambda_b
    lambda_b = 0.004
    # gamma, path loss component. Usually take 4 as value.
    gamma = 4.0
    p = 0.04
    thetha_dB = 3.0 # unit dB
    theta = np.power(10, 3.0/10)
    mu_shadow = 0.0
    # shadowing effect 的标准差，单位分贝
    sigma_dB = 12
    BETA = np.log(10.0)/10.0
    # 需要对标准差做个转换，因为对数基底的不同
    sigma_G = BETA*sigma_dB
    sigma_X = 2.0*sigma_G/gamma

    constant_A = gamma_f(1+2.0/gamma)*gamma_f(1-2.0/gamma)*theta**(2.0/gamma)

    # Define p_f_2 as the outage probability over infinite plane
    p_f_rx_div_0 = sgam.bs_rx_div_op(lambda_m, lambda_b, gamma, p, thetha_dB, 0)
    p_f_bs_nst_att_0 = sgam.bs_nearest_atch_op(lambda_m, lambda_b, gamma, p, thetha_dB, 0)
    p_f_bs_nst_att_8 = sgam.bs_nearest_atch_op(lambda_m, lambda_b, gamma, p, thetha_dB, 8)


    fig, axes = plt.subplots(1, 2, figsize=FIGSIZE, sharey=False)

    for i in range(len(axes)):
        axes[i].set_yscale(SCALE[i])
        axes[i].plot(lambda_m, p_f_rx_div_0, color='k',  marker='', linestyle='--', linewidth=LINEWIDTH, label="BS_RX_DIVERS ANA")
        axes[i].plot(lambda_m, p_f_bs_nst_att_8, color='g',  marker='', linestyle='-', linewidth=LINEWIDTH, label="BS_NST_ATT, 8dB")
        axes[i].plot(lambda_m, p_f_bs_nst_att_0, color='g',  marker='', linestyle='--', linewidth=LINEWIDTH, label="BS_NST_ATT, 0dB")


        axes[i].errorbar(
            sim_intensity_0dB,
            sim_plr_0,
            yerr=[sim_plr_lower_0, sim_plr_upper_0],
            fmt='*',
            ecolor='r',
            capthick=2,
            label="BS_RX_DIVERS SIM,0dB")
        axes[i].errorbar(
            sim_intensity_8dB,
            sim_plr_8,
            yerr=[sim_plr_lower_8, sim_plr_upper_8],
            fmt='d',
            ecolor='m',
            capthick=2,
            label="BS_RX_DIVERS SIM,8dB"
        )
        axes[i].errorbar(
            sim_intensity_nst_0dB,
            sim_plr_nst_0,
            yerr=[sim_plr_lower_nst_0, sim_plr_upper_nst_0],
            fmt='d',
            ecolor='m',
            capthick=2,
            label="BS_NST_ATT SIM,0dB"
        )
        axes[i].errorbar(
            sim_intensity_nst_8dB,
            sim_plr_nst_8,
            yerr=[sim_plr_lower_nst_8, sim_plr_upper_nst_8],
            fmt='d',
            ecolor='m',
            capthick=2,
            label="BS_NST_ATT SIM,8dB"
        )

        axes[i].grid()
        axes[i].axis([X_START, X_END, Y_START, Y_END])
        # axes[i].set_xticks(np.arange(X_START, X_END, X_STEP))
        axes[i].set_xlabel(r"Spatial density of devices")
        axes[i].set_ylabel(r"Packet loss rate")


    plt.legend(bbox_to_anchor=(0.119, 0.01, 0.79, 1), loc=1, ncol=4, mode="expand", bbox_transform=plt.gcf().transFigure)
    plt.savefig('new_packet_loss_rate_mpr.eps', format='eps', dpi=300)

    plt.show()
