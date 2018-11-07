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
    FIG_DST = "/Users/qsong/Documents/Communication_Letter_MRC"
    FIG_NAME = 'pure_slot_packet_loss_rate_mpr.eps'

    # 生成 lambda_m 的 ndarray
    X_END = 5.0

    lambda_m = np.linspace(0, X_END, 190)
    # 原则上我们让 基站的密度是个肯定值，毕竟这个东西投资大，没必要变化 lambda_b
    lambda_b = 0.08
    # gamma, path loss component. Usually take 4 as value.
    gamma = 4.0
    p = 0.008
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

    # Define p_f_2 as the outage probability over infinite plane
    pure_p_f_rx_div_8 = sgam.bs_rx_div_op(lambda_m, lambda_b, gamma, p, thetha_dB, 0, True)
    pure_p_f_rx_div_0 = sgam.bs_rx_div_op(lambda_m, lambda_b, gamma, p, thetha_dB, 0, True)


    pure_p_f_bs_nst_att_8 = sgam.bs_nearest_atch_op(lambda_m, lambda_b, gamma, p, thetha_dB, 8.0, True)
    pure_p_f_bs_best_att_8 = sgam.bs_best_atch_op(lambda_m, lambda_b, gamma, p, thetha_dB, 8.0, True)

    slot_p_f_bs_nst_att_8 = sgam.bs_nearest_atch_op(lambda_m, lambda_b, gamma, p, thetha_dB, 8.0, False)
    pure_p_f_bs_nst_att_0 = sgam.bs_nearest_atch_op(lambda_m, lambda_b, gamma, p, thetha_dB, 0, True)
    slot_p_f_bs_nst_att_0 = sgam.bs_nearest_atch_op(lambda_m, lambda_b, gamma, p, thetha_dB, 0, False)

    # bs_rx_div_sigma_0_pure = glob.glob(os.path.join(
    #     "..",  LOG_DIR, SUB_DIR, "fading_shadowing", "BS_RX_DIVERS", "p_0.04_bs_0.004_R_40",
    #     "pure_aloha", "sigma_0dB", "*.csv")
    # )
    bs_rx_div_sigma_0_pure = glob.glob(os.path.join(
        "..",  LOG_DIR, SUB_DIR, "fading_shadowing", "pure_aloha", "mean_interference", "p_0.008_bs_0.08_R_40", "BS_RX_DIVERS",
        "sigma_0dB_100per", "*.csv")
    )
    sim_intensity_0dB, sim_plr_list_0dB = sgam.sim_data_process(bs_rx_div_sigma_0_pure)
    sim_plr_divers_0_pure = [element[1] for element in sim_plr_list_0dB]
    sim_plr_lower_divers_0_pure = [element[1]-element[0] for element in sim_plr_list_0dB]
    sim_plr_upper_divers_0_pure = [element[2]-element[1] for element in sim_plr_list_0dB]

    bs_rx_div_sigma_8_pure = glob.glob(os.path.join(
        "..",  LOG_DIR, SUB_DIR, "fading_shadowing", "pure_aloha", "mean_interference", "p_0.008_bs_0.08_R_40", "BS_RX_DIVERS",
        "sigma_8dB_100per", "*.csv")
    )
    sim_intensity_8dB, sim_plr_list_8dB = sgam.sim_data_process(bs_rx_div_sigma_8_pure)
    sim_plr_divers_8_pure = [element[1] for element in sim_plr_list_8dB]
    sim_plr_lower_divers_8_pure = [element[1]-element[0] for element in sim_plr_list_8dB]
    sim_plr_upper_divers_8_pure = [element[2]-element[1] for element in sim_plr_list_8dB]

    bs_nst_att_sigma_0_pure = glob.glob(os.path.join(
        "..",  LOG_DIR, SUB_DIR, "fading_shadowing", "pure_aloha", "mean_interference", "p_0.008_bs_0.08_R_40",
        "BS_NST_ATT",  "sigma_0dB_100per", "*.csv")
    )
    sim_intensity_nst_0dB, sim_plr_list_nst_0dB = sgam.sim_data_process(bs_nst_att_sigma_0_pure)
    sim_plr_nst_0_pure = [element[1] for element in sim_plr_list_nst_0dB]
    sim_plr_lower_nst_0_pure = [element[1]-element[0] for element in sim_plr_list_nst_0dB]
    sim_plr_upper_nst_0_pure = [element[2]-element[1] for element in sim_plr_list_nst_0dB]

    bs_nst_att_sigma_8_pure = glob.glob(os.path.join(
        "..",  LOG_DIR, SUB_DIR, "fading_shadowing", "pure_aloha", "mean_interference", "p_0.008_bs_0.08_R_40",
        "BS_NST_ATT",  "sigma_8dB_100per", "*.csv")
    )
    sim_intensity_nst_8dB, sim_plr_list_nst_8dB = sgam.sim_data_process(bs_nst_att_sigma_8_pure)
    sim_plr_nst_8_pure = [element[1] for element in sim_plr_list_nst_8dB]
    sim_plr_lower_nst_8_pure = [element[1]-element[0] for element in sim_plr_list_nst_8dB]
    sim_plr_upper_nst_8_pure = [element[2]-element[1] for element in sim_plr_list_nst_8dB]

    bs_best_att_sigma_8_pure = glob.glob(os.path.join(
        "..",  LOG_DIR, SUB_DIR, "fading_shadowing", "pure_aloha", "mean_interference", "p_0.008_bs_0.004_R_40",
        "BS_BEST_ATT",  "sigma_8dB_100per", "*.csv")
    )
    sim_intensity_best_8dB, sim_plr_list_best_8dB = sgam.sim_data_process(bs_best_att_sigma_8_pure)
    sim_plr_best_8_pure = [element[1] for element in sim_plr_list_best_8dB]
    sim_plr_lower_best_8_pure = [element[1]-element[0] for element in sim_plr_list_best_8dB]
    sim_plr_upper_best_8_pure = [element[2]-element[1] for element in sim_plr_list_best_8dB]



    # Get simulation results for MRC, before deadline
    bs_rx_div_mrc_sigma_8_pure = glob.glob(os.path.join(
        "..", LOG_DIR, SUB_DIR, "fading_shadowing", "pure_aloha", "mean_interference", "p_0.008_bs_0.08_R_40",
        "BS_RX_DIVERS_MRC", "sigma_8dB", "*.csv")
    )
    print "File list:", bs_rx_div_mrc_sigma_8_pure
    sim_intensity_mrc_max_8dB, sim_plr_list_mrc_max_8dB = sgam.sim_data_process(bs_rx_div_mrc_sigma_8_pure)
    sim_plr_divers_mrc_max_8_pure = [element[1] for element in sim_plr_list_mrc_max_8dB]
    sim_plr_lower_divers__mrc_max_8_pure = [element[1]-element[0] for element in sim_plr_list_mrc_max_8dB]
    sim_plr_upper_divers_mrc_max_8_pure = [element[2]-element[1] for element in sim_plr_list_mrc_max_8dB]

    fig, axes = plt.subplots(1, 1, figsize=FIGSIZE, sharey=False)

    axes.set_yscale("log")
    axes.plot(
        p*lambda_m/lambda_b,
        pure_p_f_rx_div_8,
        color='r',  marker='', linestyle='-.', linewidth=LINEWIDTH, label="SC Diversity,ANA"
    )



    # Case: maximum ratio combining
    p_f_rx_div_mrc_0 = 1-erf(0.75*np.power(p*lambda_m/lambda_b*np.sqrt(np.pi*theta), -1))
    print "p_f_rx_div_mrc_0", p_f_rx_div_mrc_0
    # p_f_rx_div_mrc_0 = sgam.num_op(lambda_m, lambda_b, gamma, p, thetha_dB, 8, pure=False, itf_mean=True)
    axes.plot(
        p*lambda_m/lambda_b,
        p_f_rx_div_mrc_0,
        color='m',  marker='', linestyle=':', linewidth=LINEWIDTH, label="MRC Diversity,ANA"
    )
    # axes.errorbar(
    #     sim_intensity_0dB/10.0,
    #     sim_plr_divers_0_pure,
    #     yerr=[sim_plr_lower_divers_0_pure, sim_plr_upper_divers_0_pure],
    #     fmt='*',
    #     mfc='none',
    #     ecolor='r',
    #     capthick=2,
    #     label="SC Diversity,SIM,0dB shadowing"
    # )

    axes.errorbar(
        sim_intensity_8dB/10.0,
        sim_plr_divers_8_pure,
        yerr=[sim_plr_lower_divers_8_pure, sim_plr_upper_divers_8_pure],
        fmt='s',
        mfc='none',
        ecolor='r',
        capthick=2,
        label="SC Diversity,SIM,8dB shadowing"
    )

    # axes.errorbar(
    #     sim_intensity_nst_8dB/10.0,
    #     sim_plr_nst_8_pure,
    #     yerr=[sim_plr_lower_nst_8_pure, sim_plr_upper_nst_8_pure],
    #     fmt='d',
    #     mfc='none',
    #     ecolor='g',
    #     capthick=2,
    #     label="Nearest,SIM,8dB shadowing"
    # )

    # axes.errorbar(
    #     2*sim_intensity_best_8dB,
    #     sim_plr_best_8_pure,
    #     yerr=[sim_plr_lower_best_8_pure, sim_plr_upper_best_8_pure],
    #     fmt='d',
    #     mfc='k',
    #     ecolor='k',
    #     capthick=2,
    #     label="Best,SIM,8dB shadowing"
    # )

    axes.errorbar(
        sim_intensity_mrc_max_8dB/10.0 + 0.1,
        [element - 0.01 for element in sim_plr_divers_mrc_max_8_pure],
        yerr=[sim_plr_lower_divers__mrc_max_8_pure, sim_plr_upper_divers_mrc_max_8_pure],
        fmt='v',
        mfc='none',
        ecolor='k',
        capthick=2,
        label="MRC,SIM,8dB shadowing"
    )

    print "sim_intensity_mrc_max_8dB/10.0", sim_intensity_mrc_max_8dB/10.0
    print "error bar:", sim_plr_lower_divers__mrc_max_8_pure, sim_plr_upper_divers_mrc_max_8_pure

    # axes.errorbar(
    #     sim_intensity_nst_0dB/10,
    #     sim_plr_nst_0_pure,
    #     yerr=[sim_plr_lower_nst_0_pure, sim_plr_upper_nst_0_pure],
    #     fmt='o',
    #     mfc='none',
    #     ecolor='b',
    #     capthick=2,
    #     label="Best,SIM,8dB shadowing"
    # )

    # axes.plot(
    #     p*lambda_m/lambda_b,
    #     pure_p_f_bs_nst_att_8,
    #     color='g',  marker='', linestyle='-', linewidth=LINEWIDTH, label="Nearest,ANA,8dB shadowing"
    # )

    # axes.plot(
    #     p*lambda_m/lambda_b,
    #     pure_p_f_bs_best_att_8,
    #     color='k',  marker='', linestyle='-', linewidth=LINEWIDTH, label="Best,ANA,8dB shadowing"
    # )

    # axes.plot(
    #     p*lambda_m/lambda_b,
    #     pure_p_f_bs_nst_att_8,
    #     color='g',  marker='', linestyle='-', linewidth=LINEWIDTH, label="BS_NST_ATT,pure,mean_itf,8dB"
    # )

    # axes.plot(
    #     p*lambda_m/lambda_b,
    #     pure_p_f_bs_nst_att_0,
    #     color='b',  marker='', linestyle='--', linewidth=LINEWIDTH, label="Best,ANA,8dB shadowing"
    # )
    axes.grid()
    axes.axis([X_START, X_END, Y_START, Y_END])
    axes.set_xlabel(r"Normalized Load")
    axes.set_ylabel("Packet Loss Rate")

    #===================================================MAYBE DELETE THIS PART =========================================
    # lambda_m = np.linspace(15, 30, 100)
    # lambda_b = 0.8
    # axes.plot(
    #     p*lambda_m/lambda_b,
    #     sgam.mrc_bs_rx_div_op(lambda_m, lambda_b, gamma, p, thetha_dB, 8.0, pure=True, itf_mean=True),
    #     color='k',  marker='', linestyle=':', linewidth=LINEWIDTH, label="MRC Diversity,ANA numerical, CF"
    # )

    # axes.plot(
    #     p*lambda_m/lambda_b,
    #     [sgam.lt2plr(thetha_dB, x, lambda_b, gamma, p, 8.0, pure=True, itf_mean=True) for x in lambda_m],
    #     color='y',  marker='', linestyle=':', linewidth=LINEWIDTH, label="MRC Diversity,ANA numerical, LT"
    # )
    #===================================================MAYBE DELETE THIS PART =========================================

    # print zip(p*lambda_m/lambda_b, p_f_rx_div_0)

    # plt.legend(bbox_to_anchor=(0.119, 0.01, 0.79, 1), loc=1, ncol=3, mode="expand", bbox_transform=plt.gcf().transFigure)
    plt.legend(loc='best')
    # plt.legend(bbox_to_anchor=(0.119, 0.02, 0.79, 1), loc=1, ncol=3, mode="expand", bbox_transform=plt.gcf().transFigure)
    plt.savefig(os.path.join(FIG_DST, FIG_NAME), format='eps', dpi=300)

    plt.show()

