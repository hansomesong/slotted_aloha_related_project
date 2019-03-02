# -*- coding: utf-8 -*-
# Update 07-25-2018:
# Xavier suggested to resubmit our work refused by IEEE WCL to journal of Telecommunication system.
# Xavier hope that I can modify some figures about curve fitting techniques.
# The destination folder of generated figures is updated.
# We only care about the comparison between simulation results and the curve obtained by curve fitting with respect to
# simulation results.

# Update 09-15-2018:

# Update 02-08-2019:
# Xavier suggests to discuss the effect of different BS configurations in MRC, i.e. compare the performance difference
# when 2, 3, 4 best surrounding BS are involved into MRC.

__author__ = 'qsong'

import json
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
    'legend.numpoints': 2
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
    Y_START = 5e-3
    Y_END = 1.0
    Y_STEP = 0.1

    MAX_TRANS = 1
    LINEWIDTH = 2

    SCALE = ["log", "linear"]
    FIG_DST = "/Users/qsong/Dropbox/Qipeng's PhD Doc/MRC-V4-09192018"
    FIG_NAME = 'mrc_with_different_bs_settings.pdf'

    SIM_LOG_DIR = '/Users/qsong/Documents/slotted_aloha_related_project/logs/SimuApril23'


    SIM_CONFIG_PATH = os.path.join(SIM_LOG_DIR, 'sim_config.json')
    with open(SIM_CONFIG_PATH) as json_file:
        sim_config_dict = json.load(json_file)

    SKIP_NB_AVG = sim_config_dict['avg']
    SKIP_NB_CI = sim_config_dict['ci']
    # The number of entries in a table.
    CASE_NB = sim_config_dict['case_nb']

    # The dict for csv file row index for basic-based simulation, i.e., dependence interference, more close to reality
    B_SIM_MAP = sim_config_dict['dept_itf_case']
    MRC = {}
    MRC['B_SIM_all'] = B_SIM_MAP['pure_mrc_with_shadowing_avg_itf'] - 1
    MRC['B_SIM_2ru'] = B_SIM_MAP["pure_mrc_with_shadowing_avg_itf_2ru"] - 1
    MRC['B_SIM_4ru'] = B_SIM_MAP["pure_mrc_with_shadowing_avg_itf_4ru"] - 1
    MRC['B_SIM_6ru'] = B_SIM_MAP["pure_mrc_with_shadowing_avg_itf_6ru"] - 1

    file_name = "/Users/qsong/Documents/slotted_aloha_related_project/logs/SimuApril23/AlohaMultg40s8t3b20.csv"

    file_name = "/Users/qsong/Documents/slotted_aloha_related_project/logs/SimuApril23/AlohaMultg33s8t3b20-1000eNB.csv"

    plr_df = \
        pd.read_csv(file_name, sep=";", decimal=',', skiprows=SKIP_NB_AVG, nrows=CASE_NB).dropna(axis=1, how='all')

    sim_intensity = [float(x.replace(',', '.')) for x in plr_df.columns[1:].values.tolist()]

    fig, axes = plt.subplots(1, 1, figsize=FIGSIZE)

    axes.plot(
        sim_intensity,
        plr_df.loc[41, plr_df.columns != 'Charge'].values,
        marker='D',
        linewidth=LINEWIDTH,
        label="1 RU",
        linestyle='-'
    )

    axes.plot(
        sim_intensity,
        plr_df.loc[MRC['B_SIM_2ru'], plr_df.columns != 'Charge'].values,
        marker='o',
        linewidth=LINEWIDTH,
        label="2 RUs",
        linestyle='--'
    )


    axes.plot(
        sim_intensity,
        plr_df.loc[MRC['B_SIM_6ru'], plr_df.columns != 'Charge'].values,
        marker='+',
        linewidth=LINEWIDTH,
        label="6 RUs",
        linestyle='-.'
    )

    axes.plot(
        sim_intensity,
        plr_df.loc[60, plr_df.columns != 'Charge'].values,
        marker='d',
        linewidth=LINEWIDTH,
        label="20 RUs",
        linestyle=':'
    )

    axes.plot(
        sim_intensity,
        plr_df.loc[MRC['B_SIM_all'], plr_df.columns != 'Charge'].values,
        marker='s',
        linewidth=LINEWIDTH,
        label="all RUs",
        linestyle='-'
    )

    axes.grid()
    axes.axis([X_START, X_END, Y_START, Y_END])
    axes.set_xlabel(r"Normalized Load")
    axes.set_ylabel("Packet Loss Rate")

    axes.set_yscale("log")

    MARCO_DIV_TYPE = ['MRC']

    # plot curves for sc
    lines_styles = ['-', '-.', '--', ':', '-']
    mrc_marker_styles = ['s', 'o', 'v', '+', 'd', 'D']


    # plt.legend(bbox_to_anchor=(0.119, 0.01, 0.79, 1), loc=1, ncol=3, mode="expand", bbox_transform=plt.gcf().transFigure)
    plt.legend(loc='best')
    plt.tight_layout(rect=[0, 0.0, 1, 0.95])

    # plt.legend(bbox_to_anchor=(0.119, 0.02, 0.79, 1), loc=1, ncol=3, mode="expand", bbox_transform=plt.gcf().transFigure)
    plt.savefig(os.path.join(FIG_DST, FIG_NAME), format='pdf', dpi=300)
    plt.show()

