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

from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
from mpl_toolkits.axes_grid.inset_locator import inset_axes

import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
# from analytical_model.sgam import bs_nearest_atch_op, bs_rx_div_op

params = {
    'legend.fontsize': 20,
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
FIGSIZE = (12, 10)

X_START = 0.0
X_END = .1
X_STEP = 0.002
Y_START = 1e-3
Y_END = 0.5
Y_STEP = 0.1

MAX_TRANS = 1
LINEWIDTH = 3
MARKEVERY= 200
MARKER_SIZE = 11

SCALE = ["log", "linear"]

if __name__ == '__main__':


    # 生成 lambda_m 的 ndarray
    lambda_m = np.linspace(0, X_END, 1000)
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
    pure_max_p_f_bs_nst_att_8 = sgam.bs_nearest_atch_thrpt(lambda_m, lambda_b, gamma, p, thetha_dB, 8, True, False)
    slot_p_f_bs_nst_att_8 = sgam.bs_nearest_atch_thrpt(lambda_m, lambda_b, gamma, p, thetha_dB, 8, False)
    pure_p_f_bs_nst_att_0 = sgam.bs_nearest_atch_thrpt(lambda_m, lambda_b, gamma, p, thetha_dB, 0, True)
    slot_p_f_bs_nst_att_0 = sgam.bs_nearest_atch_thrpt(lambda_m, lambda_b, gamma, p, thetha_dB, 0, False)



#=============================================Best, 3 cases==================================================
    slotted_0dB_best_bs_plot_d = {
        "x": p*lambda_m/lambda_b,
        "y": (1-sgam.bs_best_atch_op(lambda_m, lambda_b, gamma, p, thetha_dB, 0, False))*p*lambda_m/lambda_b,
        "color" : 'k',
        "marker" : 'v',
        "markersize" :  MARKER_SIZE,
        "markevery" :  MARKEVERY,
        "linestyle" : '-',
        "linewidth" : LINEWIDTH,
        "label": "Best,slotted"
    }

    # slotted_8dB_best_bs_plot_d = {
    #     "x": p*lambda_m/lambda_b,
    #     "y": (1-sgam.bs_best_atch_op(lambda_m, lambda_b, gamma, p, thetha_dB, 8, False))*p*lambda_m/lambda_b,
    #     "color" : 'k',
    #     "marker" : "*",
    #     "markersize" :  MARKER_SIZE,
    #     "markevery" :  MARKEVERY,
    #     "linestyle" : '--',
    #     "linewidth" : LINEWIDTH,
    #     "label": "Best,slotted,8dB"
    # }

    best_bs_plot_d = {
        "x": p*lambda_m/lambda_b,
        "y": (1-sgam.bs_best_atch_op(lambda_m, lambda_b, gamma, p, thetha_dB, 0, True, False))*p*lambda_m/lambda_b,
        "color" : 'k',
        "marker" : '^',
        "markersize" :  MARKER_SIZE,
        "markevery" :  MARKEVERY,
        "linestyle" : '-.',
        "linewidth" : LINEWIDTH,
        "label": "Best,pure ALOHA,\nmax.interference"
    }
#=============================================Macro Diversity, 3 cases==================================================
    diver_slot = {
        "x": p*lambda_m/lambda_b,
        "y": sgam.bs_rx_div_thrpt(lambda_m, lambda_b, gamma, p, thetha_dB, 0, False),
        "color" : 'r',
        "marker" : 's',
        "markersize" :  MARKER_SIZE,
        "markevery" :  MARKEVERY,
        "linestyle" : '-',
        "linewidth" : LINEWIDTH,
        "label": "Diversity,slot"
    }

    diver_pure_avg_itf = {
        "x": p*lambda_m/lambda_b,
        "y": sgam.bs_rx_div_thrpt(lambda_m, lambda_b, gamma, p, thetha_dB, 0, True),
        "color" : 'r',
        "marker" : 'o',
        "markersize" :  MARKER_SIZE,
        "markevery" :  MARKEVERY,
        "linestyle" : '--',
        "linewidth" : LINEWIDTH,
        "label": "Diversity,pure ALOHA,\n avg.interference"
    }

    diver_pure_max_itf = {
        "x": p*lambda_m/lambda_b,
        "y": sgam.bs_rx_div_thrpt(lambda_m, lambda_b, gamma, p, thetha_dB, 0, True, False),
        "color" : 'r',
        "marker" : 'd',
        "markersize":  MARKER_SIZE,
        "markevery" :  MARKEVERY,
        "linestyle" : '-.',
        "linewidth" : LINEWIDTH,
        "label": "Diversity,pure ALOHA, \n max.interference,0dB"
    }
#=============================================Macro Diversity, 3 cases==================================================

    line_d = [
        slotted_0dB_best_bs_plot_d,
        best_bs_plot_d,
        diver_slot,
        diver_pure_avg_itf,
        diver_pure_max_itf
    ]

    fig, axes = plt.subplots(1, 1, figsize=FIGSIZE, sharey=False)

    axes.set_yscale("linear")

    for element in line_d:
        axes.plot(
            element["x"],
            element["y"],
            color=element["color"],
            marker=element["marker"],
            markevery=element["markevery"],
            markersize=element["markersize"],
            linestyle=element["linestyle"],
            linewidth=element["linewidth"],
            label=element["label"]
        )

    axes.legend(loc=4)
    axes.grid()
    # axes.axis([X_START, X_END, Y_START, Y_END])
    axes.set_xlabel(r"Normalized Load")
    axes.set_ylabel("Spatial throughput")

    # axins = zoomed_inset_axes(axes, 1, loc=4) # zoom-factor: 2.5, location: upper-left
    # change MARKEVERY from 100 to 5
    # axins = inset_axes(axes, 3, 4, loc=2, bbox_to_anchor=(0.1, 0.75), bbox_transform=axes.figure.transFigure) # no zoom
    # axins = inset_axes(axes, 3, 4, loc=1) # no zoom
    #
    # MARKEVERY = 10
    # for element in line_d:
    #     axins.plot(
    #         element["x"],
    #         element["y"],
    #         color=element["color"],
    #         marker=element["marker"],
    #         markevery=MARKEVERY,
    #         linestyle=element["linestyle"],
    #         linewidth=element["linewidth"],
    #         label=element["label"]
    #     )
    # axins.set_xlim(0.05, 0.25)
    # axins.grid()
    # axins.set_ylim(0.06, 0.18)
    # axins.set_xticks([0.05, 0.1, 0.15, 0.2, 0.25])
    # axins.patch.set_facecolor('#909090')
    #
    # mark_inset(axes, axins, loc1=2, loc2=4, fc="none", ec="0.5")

    plt.savefig('pure_slot_throughput_mpr.eps', format='eps', dpi=300)

    plt.show()

