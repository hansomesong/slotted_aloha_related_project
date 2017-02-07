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
# from analytical_model.sgam import bs_nearest_atch_op, bs_rx_div_op

params = {
    'legend.fontsize': 15,
    'xtick.labelsize': 20,
    'ytick.labelsize': 20,
    'axes.labelsize': 30,
}
plt.rcParams.update(params)

LOG_DIR = 'logs'
SUB_DIR = 'multiple_reception'
CASE_DIR = 'fading'
SUB_CASE_DIR = "bs_0.01_p_0.002_R_40|100"
DATA_FOLDED = '.'
FIGSIZE = (15, 10)

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
Y_END = 0.5
Y_STEP = 0.1

MAX_TRANS = 1
LINEWIDTH = 2

SCALE = ["log", "linear"]

if __name__ == '__main__':


    # 生成 lambda_m 的 ndarray
    lambda_m = np.linspace(0, X_END, 190)
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
    pure_p_f_rx_div_8 = sgam.bs_rx_div_op(lambda_m, lambda_b, gamma, p, thetha_dB, 8, True)
    slot_p_f_rx_div_8 = sgam.bs_rx_div_op(lambda_m, lambda_b, gamma, p, thetha_dB, 8, False)
    pure_p_f_rx_div_0 = sgam.bs_rx_div_op(lambda_m, lambda_b, gamma, p, thetha_dB, 0, True)
    slot_p_f_rx_div_0 = sgam.bs_rx_div_op(lambda_m, lambda_b, gamma, p, thetha_dB, 0, False)


    pure_p_f_bs_nst_att_8 = sgam.bs_nearest_atch_op(lambda_m, lambda_b, gamma, p, thetha_dB, 8, True)
    slot_p_f_bs_nst_att_8 = sgam.bs_nearest_atch_op(lambda_m, lambda_b, gamma, p, thetha_dB, 8, False)
    pure_p_f_bs_nst_att_0 = sgam.bs_nearest_atch_op(lambda_m, lambda_b, gamma, p, thetha_dB, 0, True)
    slot_p_f_bs_nst_att_0 = sgam.bs_nearest_atch_op(lambda_m, lambda_b, gamma, p, thetha_dB, 0, False)

    bs_rx_div_sigma_0_pure = glob.glob(os.path.join(
        "..",  LOG_DIR, SUB_DIR, "fading_shadowing", "BS_RX_DIVERS", "p_0.04_bs_0.004_R_40",
        "pure_aloha", "sigma_0dB", "*.csv")
    )
    sim_intensity_0dB, sim_plr_list_0dB = sgam.sim_data_process(bs_rx_div_sigma_0_pure)
    sim_plr_divers_0_pure = [element[1] for element in sim_plr_list_0dB]
    sim_plr_lower_divers_0_pure = [element[1]-element[0] for element in sim_plr_list_0dB]
    sim_plr_upper_divers_0_pure = [element[2]-element[1] for element in sim_plr_list_0dB]
    fig, axes = plt.subplots(1, 2, figsize=FIGSIZE, sharey=False)

    bs_nst_att_sigma_0_pure = glob.glob(os.path.join(
        "..",  LOG_DIR, SUB_DIR, "fading_shadowing", "pure_aloha", "p_0.04_bs_0.004_R_40",
        "BS_NST_ATT",  "sigma_0dB_50per", "*.csv")
    )
    sim_intensity_nst_0dB, sim_plr_list_nst_0dB = sgam.sim_data_process(bs_nst_att_sigma_0_pure)
    sim_plr_nst_0_pure = [element[1] for element in sim_plr_list_nst_0dB]
    sim_plr_lower_nst_0_pure = [element[1]-element[0] for element in sim_plr_list_nst_0dB]
    sim_plr_upper_nst_0_pure = [element[2]-element[1] for element in sim_plr_list_nst_0dB]
    fig, axes = plt.subplots(1, 2, figsize=FIGSIZE, sharey=False)

    for i in range(len(axes)):
        axes[i].set_yscale(SCALE[i])
        axes[i].plot(
            lambda_m,
            pure_p_f_rx_div_8,
            color='r',  marker='', linestyle='--', linewidth=LINEWIDTH, label="BS_RX_DIVERS,pure,mean_itf"
        )
        axes[i].errorbar(
            sim_intensity_0dB,
            sim_plr_divers_0_pure,
            yerr=[sim_plr_lower_divers_0_pure, sim_plr_upper_divers_0_pure],
            fmt='*',
            ecolor='r',
            capthick=2,
            label="BS_RX_DIVERS,SIM,0dB")
        axes[i].plot(
            lambda_m,
            sgam.bs_rx_div_op(lambda_m, lambda_b, gamma, p, thetha_dB, 8, True, False),
            color='r',  marker='', linestyle='-.', linewidth=LINEWIDTH, label="BS_RX_DIVERS,pure,max_itf"
        )
        axes[i].errorbar(
            sim_intensity_nst_0dB,
            sim_plr_nst_0_pure,
            yerr=[sim_plr_lower_nst_0_pure, sim_plr_upper_nst_0_pure],
            fmt='o',
            ecolor='b',
            capthick=2,
            label="BS_NST_ATT,SIM,0dB")
        axes[i].plot(
            lambda_m,
            slot_p_f_rx_div_8,
            color='r',  marker='', linestyle='-', linewidth=LINEWIDTH, label="BS_RX_DIVERS,slot"
        )
        axes[i].plot(
            lambda_m,
            pure_p_f_bs_nst_att_8,
            color='g',  marker='', linestyle='--', linewidth=LINEWIDTH, label="BS_NST_ATT,pure,mean_itf,8dB"
        )
        axes[i].plot(
            lambda_m,
            sgam.bs_nearest_atch_op(lambda_m, lambda_b, gamma, p, thetha_dB, 8, True, False),
            color='g',  marker='', linestyle='-.', linewidth=LINEWIDTH, label="BS_NST_ATT,pure,max_itf,8dB"
        )
        axes[i].plot(
            lambda_m,
            slot_p_f_bs_nst_att_8,
            color='g',  marker='', linestyle='-', linewidth=LINEWIDTH, label="BS_NST_ATT,slot,8dB"
        )
        axes[i].plot(
            lambda_m,
            pure_p_f_bs_nst_att_0,
            color='b',  marker='', linestyle='--', linewidth=LINEWIDTH, label="BS_NST_ATT,pure,mean_itf,0dB"
        )
        axes[i].plot(
            lambda_m,
            sgam.bs_nearest_atch_op(lambda_m, lambda_b, gamma, p, thetha_dB, 0, True, False),
            color='b',  marker='', linestyle='-.', linewidth=LINEWIDTH, label="BS_NST_ATT,pure,max_itf,0dB"
        )
        axes[i].plot(
            lambda_m,
            slot_p_f_bs_nst_att_0,
            color='b',  marker='', linestyle='-', linewidth=LINEWIDTH, label="BS_NST_ATT,slot,0dB"
        )
        axes[i].grid()
        axes[i].axis([X_START, X_END, Y_START, Y_END])
        axes[i].set_xlabel(r"Device Spatial Density")
        axes[i].set_ylabel("Packet Loss Rate")


    print sgam.bs_rx_div_op(lambda_m, lambda_b, gamma, p, thetha_dB, 8, True, False),

    plt.legend(bbox_to_anchor=(0.119, 0.02, 0.79, 1), loc=1, ncol=3, mode="expand", bbox_transform=plt.gcf().transFigure)
    plt.savefig('pure_slot_packet_loss_rate_mpr.eps', format='eps', dpi=300)

    plt.show()

