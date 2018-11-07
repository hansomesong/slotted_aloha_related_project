# -*- coding: utf-8 -*-
# This script is used for generating figures used in our second communication letter
# The idea is to generate three figures for gamma = 3.3, 4.0, 4.5
# Each figure contains three curves respectively for case: slotted ALOHA, pure ALOHA mean
# pure ALOHA max

__author__ = 'qsong'

import matplotlib.pyplot as plt
import os
import numpy as np
import matplotlib.ticker as ticker
import numpy as np

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
    'legend.fontsize': 10,
    'xtick.labelsize': 15,
    'ytick.labelsize': 15,
    'axes.labelsize': 15,
    'legend.numpoints': 2,
    'legend.handlelength': 3
}
plt.rcParams.update(params)

if __name__ == '__main__':

    # 生成 lambda_m 的 ndarray
    X_END = 5.0

    lambda_m = np.linspace(0, X_END, 190)
    # 原则上我们让 基站的密度是个肯定值，毕竟这个东西投资大，没必要变化 lambda_b
    lambda_b = 0.08
    # gamma, path loss component. Usually take 4 as value.
    p = 0.008

    L = p*lambda_m/lambda_b

    thetha_dB = 3.0 # unit dB
    theta = np.power(10.0, thetha_dB/10)

    mu_shadow = 0.0
    X_END = p*max(lambda_m)/lambda_b

    X_START = 0.0
    X_STEP = 0.002
    Y_START = 1e-2
    Y_END = 1.0
    Y_STEP = 0.1

    MAX_TRANS = 1
    LINEWIDTH = 3

    SCALE = ["log", "linear"]
    FIG_DST = "/Users/qsong/Documents/Communication_Letter_MRC_New_Organisation"


    SIM_LOG_DIR = '/Users/qsong/Documents/slotted_aloha_related_project/logs/SimuNov23PlusdePoints'

    # each element indicates that whether ALOHA type is pure. i.e. False => slotted, True  => pure
    ALOHA_TYPES = [False, True, True]
    ITF_MEANS = [True, True, False]
    COLORS = ['b', 'r', 'g', 'm', 'k', 'y']
    LINES_STYLES = ['-', ':', '--', '-.', '-']
    MRC_MRAKER_STYLES = ['s', 'o', 'd', '*', '+', 'D']
    SC_MARKER_STYLES = ['v', '^', '>', '<', 'h', 'p']

    CURVE_LABELS = [r'S-ALOHA', r'P-ALOHA,AVG', r'P-ALOHA,MAX']

    SC_SHOW = False # Whether the curves about SC should be present in the figure.
    SIM_SC_SHOW = False
    MRC_SHOW = True # Whether the curves about MRC should be present in the figure.
    SIM_MRC_SHOW = True
    FIT_MRC_SHOW = False # Whether the fitted curves about MRC should be present in the figure.
    SIM_FIT_MRC_SHOW = False # Whether the fitted curves about MRC should be present in the figure.



    gammas = [3.0, 3.3, 4.0, 4.2, 4.5, 4.7, 5.0, 6.0]
    for i, gamma in enumerate(gammas):
        gamma_label = str(int(gammas[i]*10))
        FIG_NAME = 'mrc_packet_loss_rate_gamma_{0}_theta_{1}.eps'.format(gamma_label, thetha_dB)

        # fig, axes = plt.subplots(1, 1, figsize=FIGSIZE, sharey=False)
        fig, axes = plt.subplots(1, 1)

        axes.set_yscale("log")
        axes.grid()
        axes.axis([X_START, X_END, Y_START, Y_END])
        axes.set_xlabel(r"Normalized Load")
        axes.set_ylabel("Packet Loss Rate")
        axes.set_title(r'$\gamma={0}$'.format(gamma))

        # Each element is a dict corresponding to an ALOHA type
        sim_plr_mrc_divers_group = []
        sim_plr_mrc_divers_semi_ci_group =[]
        sim_plr_sc_divers_group =[]
        sim_plr_sc_divers_semi_ci_group = []
        sim_intensity_group = []
        sim_fit_p_f_mrc_div_group = []
        p_f_mrc_div_group = []
        p_f_rx_div_group = []
        fit_p_f_mrc_div_group = []

        sim_fit_p_f_mrc_div_group = []
        sim_fit_log = '/Users/qsong/Documents/slotted_aloha_related_project/analytical_model/sim_fit_result_theta_3.csv'

        for j in range(3):
            PURE, ITF_MEAN = ALOHA_TYPES[j], ITF_MEANS[j]
            sim_intensity, sim_plr_sc_divers, sim_plr_sc_divers_semi_ci, sim_plr_mrc_divers, sim_plr_mrc_divers_semi_ci = \
            mrc_curve_fitting.sim_parser(SIM_LOG_DIR, PURE, ITF_MEAN)

            sim_intensity_group.append(sim_intensity)
            sim_plr_mrc_divers_group.append(sim_plr_mrc_divers)
            sim_plr_mrc_divers_semi_ci_group.append(sim_plr_mrc_divers_semi_ci)
            sim_plr_sc_divers_group.append(sim_plr_sc_divers)
            sim_plr_sc_divers_semi_ci_group.append(sim_plr_sc_divers_semi_ci)

            p_f_rx_div, p_f_mrc_div, fit_p_f_mrc_div, empirical_p_f_mrc_div \
                = mrc_curve_fitting.sc_mrc_anayltical_parser(lambda_m, lambda_b, p, thetha_dB, gammas, PURE, ITF_MEAN)
            p_f_mrc_div_group.append(p_f_mrc_div)
            p_f_rx_div_group.append(p_f_rx_div)
            fit_p_f_mrc_div_group.append(fit_p_f_mrc_div)


            sim_fit_p_f_mrc_div ={}
            sim_fit_p_f_mrc_div[gamma_label] = mrc_curve_fitting.sim_fitted_function(
                    sim_fit_log, thetha_dB, L, gamma, PURE, ITF_MEAN
            )


            # plot curves for sc
            # Iterate for ALOHA type
            if SC_SHOW:
                CURVE_LABEL = CURVE_LABELS[j] + r",SC,ANA"
                axes.plot(
                L,
                p_f_rx_div_group[j][gamma_label],
                color=COLORS[j],  marker='', linestyle=LINES_STYLES[j], linewidth=LINEWIDTH, label=CURVE_LABEL
                )
            # plot curves for MRC
            if MRC_SHOW:
                CURVE_LABEL = CURVE_LABELS[j] + r",ANA"
                axes.plot(
                    L,
                    p_f_mrc_div_group[j][gamma_label],
                    color=COLORS[j],  marker='', linestyle=LINES_STYLES[j], linewidth=LINEWIDTH, label=CURVE_LABEL
                )

            if FIT_MRC_SHOW:
                CURVE_LABEL = CURVE_LABELS[j] + r",FIT_ANA"
                axes.plot(
                    L,
                    fit_p_f_mrc_div_group[j][gamma_label],
                    color=COLORS[j],  marker='', linestyle=LINES_STYLES[j], linewidth=LINEWIDTH, label=CURVE_LABEL
                )

            if SIM_FIT_MRC_SHOW:
                CURVE_LABEL = CURVE_LABELS[j] + r",FIT_SIM"
                axes.plot(
                    L,
                    sim_fit_p_f_mrc_div_group[j][gamma_label],
                    color=COLORS[j],  marker='', linestyle=LINES_STYLES[j], linewidth=LINEWIDTH, label=CURVE_LABEL
                )

            if SIM_SC_SHOW:
                axes.errorbar(
                    sim_intensity_group[j][gamma_label],
                    sim_plr_sc_divers_group[j][gamma_label],
                    yerr=[sim_plr_sc_divers_semi_ci_group[j][gamma_label],  sim_plr_sc_divers_semi_ci_group[j][gamma_label]],
                    fmt=SC_MARKER_STYLES[i],
                    mfc='none',
                    ecolor=COLORS[j],
                    capthick=2,
                    label="SC,SIM"
                )

            if SIM_MRC_SHOW:
                curve_label = CURVE_LABELS[j] + r",SIM"
                axes.errorbar(
                    sim_intensity_group[j][gamma_label],
                    sim_plr_mrc_divers_group[j][gamma_label],
                    yerr=[sim_plr_mrc_divers_semi_ci_group[j][gamma_label],  sim_plr_mrc_divers_semi_ci_group[j][gamma_label]],
                    fmt=MRC_MRAKER_STYLES[j],
                    mfc='none',
                    ecolor=COLORS[j],
                    capthick=2,
                    label= curve_label
                )

            plt.legend(loc='best')
            fig.tight_layout()

        # plt.savefig(os.path.join(FIG_DST, FIG_NAME), format='eps', bbox_inches='tight', dpi=300)
        # Without parameter setting bbox_inches='tight', the generated eps figure, I observe that y-lable is
        # is disappreaded and x-lable is cut off. Two possible causes:
        # 1) I should call fig.tight_layout at the end
        # 2) bbox_inches='tight' has effect when saving figures while fig.tight_layout() has no effect.
        # plt.savefig(os.path.join(FIG_DST, FIG_NAME), format='eps', dpi=300)





    # lines = axes.get_lines()
    # legend1 = plt.legend([lines[i] for i in [0,1,2]], ["algo1", "algo2", "algo3"], loc=1)
    # legend2 = plt.legend([lines[i] for i in [0,3,6]], parameters, loc=4)
    # axes.add_artist(legend1)
    # axes.add_artist(legend2)

    # plt.legend(bbox_to_anchor=(0.119, 0.01, 0.79, 1), loc=1, ncol=3, mode="expand", bbox_transform=plt.gcf().transFigure)
    # plt.legend(bbox_to_anchor=(0.119, 0.02, 0.79, 1), loc=1, ncol=3, mode="expand", bbox_transform=plt.gcf().transFigure)
    # plt.savefig(os.path.join(FIG_DST, FIG_NAME), format='eps', dpi=300)

    plt.show()
