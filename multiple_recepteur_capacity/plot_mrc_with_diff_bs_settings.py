# -*- coding: utf-8 -*-
# This script is used to generate figures used in our second IEEE communication letter
# The idea is to generate three figures for gamma = 3.3, 4.0, 4.5
# Each figure contains three curves respectively for case: slotted ALOHA, pure ALOHA mean
# pure ALOHA max

# 2018-01-20
# Update: our submission to IEEE communication letter has been rejected.
# Now we aim at submitting to IEEE wireless communication letter

# 2018-04-17
# Update: Bad news! After nearly two months review, our submission to IEEE WCL has also been rejected
# but can be resubmitted after revision.
# Xavier has conducted a new series of simulation, with a 500X500 km area.
# Let's see whether in this case the bias between simulation and analytical is still obvious.

# 2019-02-08
# Update: We submit our work to Springer Telecommunication system.
# Now We received a major revision.
# To respond to first comment of reviewer 1, Xavier suggest to illustrate the performance with different BS settings.

# Sur simulation les numéros correspondent à

#% 1 - Slotte avec shadow et macrodiv SC, ex nbPcktOK
#% 2 - Slotte sans shadow et macrodiv SC, ex nbPcktOKSsShadow
#% 3 - Slotte avec shadow sur meilleur, ex nbPcktOKBest
#% 4 - Slotte sans shadow sur meilleur, ex nbPcktOKSsShadowBest
#% 5 - Pur interf moyenne avec shadow et macrodiv SC, ex nbPcktOKPureAvgAvecShad
#% 6 - Pur interf moyenne sans shadow et macrodiv SC, ex nbPcktOKPureAvgSsShad
#% 7 - Pur interf moyenne avec shadow sur meilleur, ex nbPcktOKPureAvgAvecShadBest
#% 8 - Pur interf moyenne  sans shadow sur meilleur, ex nbPcktOKPureAvgSsShadBest
#% 9 - Pur interf max avec shadow et macrodiv SC,  ex nbPcktOKPureMaxAvecShad
#% 10 - Pur interf max sans shadow sur plus proche, ex nbPcktOKPureMaxSsShad
#% 11 - Pur interf max avec shadow et macrodiv SC, ex nbPcktOKPureMaxAvecShadBest
#% 12 - Pur interf max  sans shadow sur meilleur, ex nbPcktOKPureMaxSsShadBest
#% 13 - slotte avec shadow et macrodiv MRC,  ex nbPcktMrcOK
#% 14 - Pur interf moyenne avec shadow et macrodiv MRC, ex nbPcktMrcOKPureAvg
#% 15 - Pur interf max avec shadow et macrodiv MRC, ex nbPcktMrcOKPureMax
#% 16 - Slotte avec shadow et macrodiv SC mais interferences independantes, nbPckScIIDInter
#% 17 - Slotte avec shadow et macrodiv MRC mais interferences independantes,nbPckMrcIIDInter
#% 18 - Pure interf moyenne avec shadow et macro SC mais interferences independantes
#% 19 - Pure interf moyenne avec shadow et macro MRC mais interferences independantes,
#% 20 - Pure interf max avec shadow et macro SC mais interferences independantes,
#% 21 - Pure interf max avec shadow et macro MRCmais interferences independantes,

# Les numéros 22 à 41 correspondent à MRC limité en slotted Aloha, 42 à 61 à MRC limité en pure Aloha avg Inteférence,
# 62 à 81 à MRC limité en pure Aloha max Int. Chaque ligne correspond à un nombre de branches. Par exemple 47 correspond
# à MRC où on somme les 47-42+1=6 meilleurs SIR en pure Aloha avg Interf.
# Tu remarques que les cas 22, 42 et  62 correspondent en fait au SC (on ne prend que le meilleur SIR,
# si cela marche, c'est bon ; en revanche si le meilleur est inférieur au seuil, tous les autres sont également plus petits).

__author__ = 'qsong'

import matplotlib.pyplot as plt
import os
from simulator_config import SIMULATOR_PATH
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
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'axes.labelsize': 15,
    'legend.numpoints': 2,
    'legend.handlelength': 3
}
plt.rcParams.update(params)

from matplotlib import rc
rc('font', **{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

if __name__ == '__main__':

    # 生成 lambda_m 的 ndarray
    X_END = 4.0

    # Let the spatial density of BS is constant value and make the density of devices vary
    lambda_b = 0.08
    lambda_m = np.linspace(1, X_END, 200)

    p = 0.008
    # L is the normalized load
    L = p*lambda_m/lambda_b

    thetha_dB = 3.0 # unit dB
    theta = np.power(10.0, thetha_dB/10)

    mu_shadow = 0.0
    L_END = p * max(lambda_m) / lambda_b

    X_START = p * min(lambda_m) / lambda_b
    X_STEP = 0.002
    Y_START = 1e-2
    Y_END = 1.0
    Y_STEP = 0.1

    MAX_TRANS = 1
    LINEWIDTH = 2
    FIGSIZE = (15, 5)

    SCALE = ["log", "linear"]
    FIG_DST = os.path.join(SIMULATOR_PATH, 'figures', 'wireless_communication_letter')

    if not os.path.exists(FIG_DST):
        exit("The folder in which the generated figure will be saved does not exist!"
             "The given path is {0}: "
             "Please double check and modify the path in your PC!!".format(FIG_DST))

    print FIG_DST

    # Simulation data used to generate figures in submission to IEEE communication letter
    # is: SimuNov23PlusdePoints. The simulation data is generated by Xavier's MATLab programs.
    # SIM_LOG_DIR = os.path.join(SIMULATOR_PATH, 'logs', 'SimuJan20')

    # Update: 2018-04-16, generate new figures with new simulation results.
    SIM_LOG_DIR = os.path.join(SIMULATOR_PATH, 'logs', 'SimuApril23')


    if not os.path.exists(SIM_LOG_DIR):
        exit("The folder in which the simulation data is saved does not exist!"
             "The given folder is: {0}".format(SIM_LOG_DIR))

    # each element indicates that whether ALOHA type is pure. i.e. False => slotted, True  => pure
    ALOHA_TYPES = [False, True, True]
    ITF_MEANS = [True, True, False]
    # Whether the interference at BS is independent.
    ITF_INDEPENDENCES = [True, True, True]
    # ITF_INDEPENDENCES = [False, False, False]

    COLORS = ['b', 'r', 'g', 'm', 'k', 'y']
    LINES_STYLES = ['-', ':', '--', '-.', '-']
    MRC_MRAKER_STYLES = ['s', 'o', 'd', '*', '+', 'D']
    SC_MARKER_STYLES = ['v', '^', '>', '<', 'h', 'p']

    CURVE_LABELS = [r'S-ALOHA', r'Pa-ALOHA', r'Pm-ALOHA']

    # B-SIM approach (basic simulation),  a set of interferers is considered
    # for each transmitted packet and the interferences are computed.
    # Hence, the interference on each RU is computed without any independence assumption.
    B_SIM = True
    # In the M-SIM approach (model-based simulation)
    # we consider independent interferences for different RUs and
    # we thus should get results close the ones given by the theoretical analysis.
    M_SIM = True

    SC_SHOW = False # Whether the curves about SC should be present in the figure.
    SIM_SC_SHOW = False
    MRC_SHOW = True # Whether the curves about MRC should be present in the figure.
    SIM_MRC_SHOW = True
    FIT_MRC_SHOW = False # Whether the fitted curves about MRC should be present in the figure.
    SIM_FIT_MRC_SHOW = False # Whether the fitted curves about MRC should be present in the figure.







    # gammas = [3.0, 3.3, 4.0, 4.2, 4.5, 4.7, 5.0, 6.0]
    gammas = [3.3, 4.0, 4.5]

    fig, axes = plt.subplots(1, 3, figsize=FIGSIZE, gridspec_kw = {'wspace':0.08, 'hspace':0}, sharey=True)

    FIG_NAME = 'mrc_packet_loss_rate_gamma_3in1.eps'
    axes[0].set_ylabel("Packet Loss Rate")

    for i, gamma in enumerate(gammas):
        gamma_label = str(int(gammas[i]*10))

        axes[i].set_yscale("log")
        axes[i].grid('on', linestyle='--')
        axes[i].axis([X_START, L_END, Y_START, Y_END])
        axes[i].set_xlabel(r"Normalized Load")
        # axes[i].set_ylabel("Packet Loss Rate")
        # axes[i].set_title(r'$\gamma={0}$'.format(gamma))

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
        sim_fit_log = os.path.join(SIMULATOR_PATH, 'analytical_model', 'sim_fit_result_theta_3.csv')
        if not os.path.exists(sim_fit_log):
            exit("The folder in which the fit data according to simulation data"
                 " is saved does not exist!"
                 "The given folder is: {0}".format(sim_fit_log))

        for j in range(3):
            PURE, ITF_MEAN, ITF_INDEPENDENCE = ALOHA_TYPES[j], ITF_MEANS[j], ITF_INDEPENDENCES[j]
            sim_intensity, sim_plr_sc_divers, sim_plr_sc_divers_semi_ci, sim_plr_mrc_divers, sim_plr_mrc_divers_semi_ci = \
            mrc_curve_fitting.sim_parser(SIM_LOG_DIR, PURE, ITF_MEAN, ITF_INDEPENDENCE)

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
                axes[i].plot(
                L,
                p_f_rx_div_group[j][gamma_label],
                color=COLORS[j],  marker='', linestyle=LINES_STYLES[j], linewidth=LINEWIDTH, label=CURVE_LABEL
                )
            # plot curves for MRC
            if MRC_SHOW:
                CURVE_LABEL = CURVE_LABELS[j] + r",ANA"
                curve_ana_mrc, = axes[i].plot(
                    L,
                    p_f_mrc_div_group[j][gamma_label],
                    color=COLORS[j],  marker='', linestyle=LINES_STYLES[j], linewidth=LINEWIDTH, label=CURVE_LABEL
                )

            if FIT_MRC_SHOW:
                CURVE_LABEL = CURVE_LABELS[j] + r",FIT_ANA"
                curve_ana_fitmrc, = axes[i].plot(
                    L,
                    fit_p_f_mrc_div_group[j][gamma_label],
                    color=COLORS[j],  marker='', linestyle=LINES_STYLES[j], linewidth=LINEWIDTH, label=CURVE_LABEL
                )

            if SIM_FIT_MRC_SHOW:
                CURVE_LABEL = CURVE_LABELS[j] + r",FIT_SIM"
                axes[i].plot(
                    L,
                    sim_fit_p_f_mrc_div_group[j][gamma_label],
                    color=COLORS[j],  marker='', linestyle=LINES_STYLES[j], linewidth=LINEWIDTH, label=CURVE_LABEL
                )

            if SIM_SC_SHOW:
                axes[i].errorbar(
                    sim_intensity_group[j][gamma_label],
                    sim_plr_sc_divers_group[j][gamma_label],
                    yerr=[sim_plr_sc_divers_semi_ci_group[j][gamma_label],  sim_plr_sc_divers_semi_ci_group[j][gamma_label]],
                    fmt=SC_MARKER_STYLES[i],
                    mfc='none',
                    ecolor=COLORS[j],
                    capthick=2,
                    label="SC,SIM"
                )

            if SIM_MRC_SHOW and M_SIM:
                # simulation curve corresponding to model-based simulation. A ideal but non-realistic model
                curve_label = CURVE_LABELS[j] + r",M-SIM"
                axes[i].errorbar(
                    sim_intensity_group[j][gamma_label],
                    sim_plr_mrc_divers_group[j]['M_SIM'][gamma_label],
                    yerr=[
                        sim_plr_mrc_divers_semi_ci_group[j]['M_SIM'][gamma_label],
                        sim_plr_mrc_divers_semi_ci_group[j]['M_SIM'][gamma_label]
                    ],
                    fmt=MRC_MRAKER_STYLES[j],
                    mfc=COLORS[j], # I guess 'mfc' is short for 'marker fulfilled color'?
                    markeredgecolor=COLORS[j],
                    ecolor=COLORS[j],
                    capthick=2,
                    label= curve_label
                )

            if SIM_MRC_SHOW and B_SIM:
                # simulation curve corresponding to basic model simulation
                curve_label = CURVE_LABELS[j] + r",R-SIM"
                axes[i].errorbar(
                    sim_intensity_group[j][gamma_label],
                    sim_plr_mrc_divers_group[j]['B_SIM'][gamma_label],
                    yerr=[
                        sim_plr_mrc_divers_semi_ci_group[j]['B_SIM'][gamma_label],
                        sim_plr_mrc_divers_semi_ci_group[j]['B_SIM'][gamma_label]
                    ],
                    fmt=MRC_MRAKER_STYLES[j],
                    mfc='none',
                    markeredgecolor = COLORS[j],
                    ecolor=COLORS[j],
                    capthick=2,
                    label= curve_label
                )

            # plt.legend(loc='best')
            # fig.tight_layout()

    axes[0].set_title(r'(a)$\gamma=3.3$')
    axes[1].set_title(r'(b)$\gamma=4.0$')
    axes[2].set_title(r'(c)$\gamma=4.5$')
    handles, labels = axes[0].get_legend_handles_labels()
    handles_ordered = [handles[0], handles[3], handles[4], handles[1], handles[5], handles[6], handles[2], handles[7], handles[8]]
    labels_ordered = [labels[0], labels[3], labels[4], labels[1], labels[5], labels[6], labels[2], labels[7], labels[8]]
    plt.legend(
        handles_ordered,
        labels_ordered,
        bbox_to_anchor=(0.119, 0.96, 0.785, .102),
        loc='upper center',
        ncol=3,
        mode="expand",
        bbox_transform=plt.gcf().transFigure
    )
    # plt.savefig(os.path.join(FIG_DST, FIG_NAME), format='eps', bbox_inches='tight', dpi=300)
    # Without parameter setting bbox_inches='tight', the generated eps figure, I observe that y-lable is
    # is disappreaded and x-lable is cut off. Two possible causes:
    # 1) I should call fig.tight_layout at the end
    # 2) bbox_inches='tight' has effect when saving figures while fig.tight_layout() has no effect.
    # fig.tight_layout()
    # fig.savefig(os.path.join(FIG_DST, FIG_NAME), bbox_inches='tight', format='eps', dpi=300)





    # lines = axes[i].get_lines()
    # legend1 = plt.legend([lines[i] for i in [0,1,2]], ["algo1", "algo2", "algo3"], loc=1)
    # legend2 = plt.legend([lines[i] for i in [0,3,6]], parameters, loc=4)
    # axes[i].add_artist(legend1)
    # axes[i].add_artist(legend2)

    # plt.legend(bbox_to_anchor=(0.119, 0.01, 0.79, 1), loc=1, ncol=3, mode="expand", bbox_transform=plt.gcf().transFigure)
    # plt.legend(bbox_to_anchor=(0.119, 0.02, 0.79, 1), loc=1, ncol=3, mode="expand", bbox_transform=plt.gcf().transFigure)
    # plt.savefig(os.path.join(FIG_DST, FIG_NAME), format='eps', dpi=300)


    # The following loop is used to calculate the max load under different gamma
    for gamma in gammas:
        print "gamma", gamma, "theta", thetha_dB
        print "sc max load", sgam.div_max_load(gamma, thetha_dB, p_max=0.1, pure=True, itf_mean=True)
    plt.show()

