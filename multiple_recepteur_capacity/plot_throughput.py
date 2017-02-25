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
FIGSIZE = (15, 8)

A_P_DECREMENT ="analytical, power decrement"
A_P_IDENTIC = "analytical, identical power"
A_P_INCREMENT = "analytical, power increment"
S_P_IDENTIC = "simulation, identical power"
S_P_INCREMENT = "simulation, power increment"
S_P_DECREMENT = "simulation, power decrement"

X_START = 0.0
X_END = 1.0
X_STEP = 0.002
Y_START = 1e-3
Y_END = 0.5
Y_STEP = 0.1

MAX_TRANS = 1
LINEWIDTH = 2

SCALE = ["log", "linear"]

if __name__ == '__main__':


    # 生成 lambda_m 的 ndarray
    lambda_m = np.linspace(0, X_END, 100)
    # 原则上我们让 基站的密度是个肯定值，毕竟这个东西投资大，没必要变化 lambda_b
    lambda_b = 0.008
    # gamma, path loss component. Usually take 4 as value.
    gamma = 4.0
    p = 0.08
    thetha_dB = 3.0 # unit dB

    Y_END = p*X_END

    # Define p_f_2 as the outage probability over infinite plane
    pure_p_f_rx_div_8 = sgam.bs_rx_div_thrpt(lambda_m, lambda_b, gamma, p, thetha_dB, 8, True)
    slot_p_f_rx_div_8 = sgam.bs_rx_div_thrpt(lambda_m, lambda_b, gamma, p, thetha_dB, 8, False)
    pure_p_f_rx_div_0 = sgam.bs_rx_div_thrpt(lambda_m, lambda_b, gamma, p, thetha_dB, 0, True)
    slot_p_f_rx_div_0 = sgam.bs_rx_div_thrpt(lambda_m, lambda_b, gamma, p, thetha_dB, 0, False)


    pure_p_f_bs_nst_att_8 = sgam.bs_nearest_atch_thrpt(lambda_m, lambda_b, gamma, p, thetha_dB, 8, True)
    slot_p_f_bs_nst_att_8 = sgam.bs_nearest_atch_thrpt(lambda_m, lambda_b, gamma, p, thetha_dB, 8, False)
    pure_p_f_bs_nst_att_0 = sgam.bs_nearest_atch_thrpt(lambda_m, lambda_b, gamma, p, thetha_dB, 0, True)
    slot_p_f_bs_nst_att_0 = sgam.bs_nearest_atch_thrpt(lambda_m, lambda_b, gamma, p, thetha_dB, 0, False)

    fig, axes = plt.subplots(1, 1, figsize=FIGSIZE, sharey=False)

    axes.set_yscale("linear")
    axes.plot(
        p*lambda_m/lambda_b,
        pure_p_f_rx_div_8,
        color='r',  marker='', linestyle='--', linewidth=LINEWIDTH, label="BS_RX_DIVERS,pure,mean_itf"
    )
    axes.plot(
        p*lambda_m/lambda_b,
        slot_p_f_rx_div_8,
        color='r',  marker='', linestyle='-', linewidth=LINEWIDTH, label="BS_RX_DIVERS,slot"
    )
    axes.plot(
        p*lambda_m/lambda_b,
        pure_p_f_bs_nst_att_8,
        color='g',  marker='', linestyle='--', linewidth=LINEWIDTH, label="BS_NST_ATT,pure,mean_itf,8dB"
    )

    # axes.plot(
    #     lambda_m,
    #     sgam.bs_nearest_atch_op(lambda_m, lambda_b, gamma, p, thetha_dB, 8, True, False),
    #     color='g',  marker='', linestyle='-.', linewidth=LINEWIDTH, label="BS_NST_ATT,pure,max_itf,8dB"
    # )
    axes.plot(
        p*lambda_m/lambda_b,
        slot_p_f_bs_nst_att_8,
        color='g',  marker='', linestyle='-', linewidth=LINEWIDTH, label="BS_NST_ATT,slot,8dB"
    )
    axes.plot(
        p*lambda_m/lambda_b,
        pure_p_f_bs_nst_att_0,
        color='b',  marker='', linestyle='--', linewidth=LINEWIDTH, label="BS_NST_ATT,pure,mean_itf,0dB"
    )
    # axes.plot(
    #     lambda_m,
    #     sgam.bs_nearest_atch_op(lambda_m, lambda_b, gamma, p, thetha_dB, 0, True, False),
    #     color='b',  marker='', linestyle='-.', linewidth=LINEWIDTH, label="BS_NST_ATT,pure,max_itf,0dB"
    # )
    axes.plot(
        p*lambda_m/lambda_b,
        slot_p_f_bs_nst_att_0,
        color='b',  marker='', linestyle='-', linewidth=LINEWIDTH, label="BS_NST_ATT,slot,0dB"
    )
    axes.grid()
    # axes.axis([X_START, X_END, Y_START, Y_END])
    axes.set_xlabel(r"Device Spatial Density")
    axes.set_ylabel("Packet Loss Rate")


    print 1-sgam.bs_rx_div_op(lambda_m, lambda_b, gamma, p, thetha_dB, 8, True, False),
    print lambda_m
    print (1-sgam.bs_rx_div_op(lambda_m, lambda_b, gamma, p, thetha_dB, 8, True, False))*lambda_m
    plt.legend(bbox_to_anchor=(0.119, 0.02, 0.79, 1), loc=1, ncol=3, mode="expand", bbox_transform=plt.gcf().transFigure)
    plt.savefig('pure_slot_throughput_mpr.eps', format='eps', dpi=300)

    plt.show()

