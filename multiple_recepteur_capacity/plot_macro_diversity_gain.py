# -*- coding: utf-8 -*-
# This script is used for plotting the evolution of macro reception diversity with respect to maximum allowed probability
# loss rate
__author__ = 'qsong'


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


# To avoid using type-3 fonts, which are not allowed by some journal.
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
# from analytical_model.sgam import bs_nearest_atch_op, bs_rx_div_op

params = {
    'legend.fontsize': 20,
    'xtick.labelsize': 20,
    'ytick.labelsize': 20,
    'axes.labelsize': 20,
    'text.usetex' : True
}
plt.rcParams.update(params)


FIGSIZE = (12, 10)
LINEWIDTH = 3
MARKER_SIZE = 8
MARKEVERY  = 50

if __name__ == "__main__":
    GAMMA = 4
    SIGMA_DB = 8
    P_F_MAX = np.linspace(1e-3, 1e-1, 1000)
    macro_div_gain2best = sgam.macro_div_gain(p_f_max=P_F_MAX, target="BEST", sigma_dB=SIGMA_DB, gamma=GAMMA, trans_rep_nb=1, max_trans_nb=1)
    macro_div_gain2nearest = sgam.macro_div_gain(p_f_max=P_F_MAX, target="NEAREST", sigma_dB=SIGMA_DB, gamma=GAMMA, trans_rep_nb=1, max_trans_nb=1)

    m_macro_div_gain2best = sgam.macro_div_gain(p_f_max=P_F_MAX, target="BEST", sigma_dB=SIGMA_DB, gamma=GAMMA, trans_rep_nb=3, max_trans_nb=3)
    m_macro_div_gain2nearest = sgam.macro_div_gain(p_f_max=P_F_MAX, target="NEAREST", sigma_dB=SIGMA_DB, gamma=GAMMA, trans_rep_nb=3, max_trans_nb=3)

    fig, axes = plt.subplots(1, 1, figsize=FIGSIZE, sharey=False)
    print sgam.macro_div_gain(p_f_max=0.1, target="BEST", sigma_dB=SIGMA_DB, gamma=GAMMA, trans_rep_nb=1, max_trans_nb=1)
    print sgam.macro_div_gain(p_f_max=0.001, target="BEST", sigma_dB=SIGMA_DB, gamma=GAMMA, trans_rep_nb=3, max_trans_nb=3)

    axes.set_yscale("linear")
    axes.set_xscale('log')

    axes.plot(
        P_F_MAX,
        m_macro_div_gain2best,
        color='g',
        marker='o',
        markersize=MARKER_SIZE,
        markevery=MARKEVERY,
        linestyle='--',
        linewidth=LINEWIDTH,
        label=r"Macro Diversity Gain, $N_{rep} = 3, N_{max\_trans} = 3$"
    )

    axes.plot(
        P_F_MAX,
        sgam.macro_div_gain(p_f_max=P_F_MAX, target="BEST", sigma_dB=SIGMA_DB, gamma=GAMMA, trans_rep_nb=4, max_trans_nb=4),
        color='r',
        marker='*',
        markersize=MARKER_SIZE,
        markevery=MARKEVERY,
        linestyle='-.',
        linewidth=LINEWIDTH,
        label=r"Macro Diversity Gain, $N_{rep} = 4, N_{max\_trans} = 4$"
    )

    axes.plot(
        P_F_MAX,
        sgam.macro_div_gain(p_f_max=P_F_MAX, target="BEST", sigma_dB=SIGMA_DB, gamma=GAMMA, trans_rep_nb=5, max_trans_nb=5),
        color='b',
        marker='s',
        markersize=MARKER_SIZE,
        markevery=MARKEVERY,
        linestyle='-',
        linewidth=LINEWIDTH,
        label=r"Macro Diversity Gain, $N_{rep} = 5, N_{max\_trans} = 5$"
    )

    axes.set_xlabel(r"Packet Loss Rate")
    axes.set_ylabel(r"Macro Diversity Gain against Best BS Attach" )
    axes.grid()
    # axes.plot(
    #     P_F_MAX,
    #     macro_div_gain2nearest,
    #     color='b',
    #     marker='o',
    #     markersize=MARKER_SIZE,
    #     markevery=MARKEVERY,
    #     linestyle='-',
    #     linewidth=LINEWIDTH,
    #     label="Macro Diversity Gain against Best Nearest Attach"
    # )
    plt.legend(loc='best')
    plt.savefig('macro_div_gain_retrans.eps', format='eps', dpi=300)

    plt.show()

