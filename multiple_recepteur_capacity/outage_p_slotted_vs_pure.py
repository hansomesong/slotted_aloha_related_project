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
FIGSIZE = (15, 6)

A_P_DECREMENT ="analytical, power decrement"
A_P_IDENTIC = "analytical, identical power"
A_P_INCREMENT = "analytical, power increment"
S_P_IDENTIC = "simulation, identical power"
S_P_INCREMENT = "simulation, power increment"
S_P_DECREMENT = "simulation, power decrement"

X_START = 0.0
X_END = 2.0
X_STEP = 0.1
Y_START = 1e-3
Y_END = 1.0
Y_STEP = 0.1

MAX_TRANS = 1
LINEWIDTH = 2

SCALE = ["linear", "log"]

if __name__ == '__main__':


    # 生成 lambda_m 的 ndarray
    lambda_m = np.linspace(0, X_END, 190)
    # 原则上我们让 基站的密度是个肯定值，毕竟这个东西投资大，没必要变化 lambda_b
    lambda_b = 0.05
    # gamma, path loss component. Usually take 4 as value.
    gamma = 4.0
    p = 0.01
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

    fig, axes = plt.subplots(1, 1, figsize=FIGSIZE, sharey=False)

    axes.set_yscale(SCALE[0])
    axes.plot(
        lambda_m,
        slot_p_f_rx_div_8/pure_p_f_rx_div_8,
        color='r',  marker='', linestyle='-', linewidth=LINEWIDTH, label="OP Ratio,BS_RX_DIVERS,pure"
    )
    axes.plot(
        lambda_m,
        slot_p_f_bs_nst_att_0/pure_p_f_bs_nst_att_0,
        color='g',  marker='', linestyle='--', linewidth=LINEWIDTH, label="OP Ratio,BS_NST_ATT,0dB"
    )
    axes.plot(
        lambda_m,
        slot_p_f_bs_nst_att_8/pure_p_f_bs_nst_att_8,
        color='b',  marker='', linestyle='-', linewidth=LINEWIDTH, label="OP Ratio,BS_NST_ATT,8dB"
    )
    axes.grid()
    axes.axis([X_START, X_END, Y_START, Y_END])
    axes.set_xticks(np.arange(X_START, X_END, X_STEP))
    axes.set_title("Outage Probability Ratio between slotted and Pure Aloha")
    axes.set_xlabel(r"path loss exponent")
    axes.legend(loc='best', numpoints=2)

    print slot_p_f_bs_nst_att_8/pure_p_f_bs_nst_att_8,
    plt.savefig('op_ratio_slotted_vs_pure.eps', format='eps', dpi=300)

    plt.show()

