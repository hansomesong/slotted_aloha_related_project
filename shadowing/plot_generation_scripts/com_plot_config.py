__author__ = 'qsong'

import matplotlib.pyplot as plt
import os

params = {
    'legend.fontsize': 17,
    'xtick.labelsize': 15,
    'ytick.labelsize': 15,
    'axes.labelsize': 20,
}
plt.rcParams.update(params)

# mpl.rcParams['text.usetex'] = True
# mpl.rcParams.update({'figure.autolayout': True})
FIGSIZE = (10, 10)
FIG_DEST = '/Users/qsong/Documents/Final_Thesis/phd-thesis-template-2.2.2/Chapter4/Figures'
PROJECT_DST = '/Users/qsong/Documents/slotted_aloha_related_project'


A_P_IDENTIC = "$v=1.0$, $\sigma=1.0$"
A_P_INCREMENT = "$v=2.0$, $\sigma=1.0$"
A_P_DECREMENT ="$v=0.5$, $\sigma=1.0$"

P_A_P_IDENTIC = "$v=1.0$, \n$\sigma=0.0$"
P_A_P_INCREMENT = "$v=2.0$, \n$\sigma=0.0$"
P_A_P_DECREMENT = "$v=0.5$, \n$\sigma=0.0$"

MAX_TRANS = 5
LOG_DIR = os.path.join(PROJECT_DST, 'logs')
SUB_DIR = 'shadowing'
ANA_DIR = 'analytical'

SETTING_0 ={
    'MAX_TRANS': 5,
    'LOG_DIR':  'logs',
    'SUB_DIR':  'perfect',
    'ANA_DIR':  'analytical',
    'BACKOFF': 36,
    'MU_FADING': 0.0,
    'SIGMA_S':  0.0,
    'LINEWIDTH': 2,
    'LM_PAIR':  [(1, 2), (1, 1), (2, 1)],
    'LINE_COLOR': ['r', 'b', 'g'],
    'LINE_TYPE': ['--', '-', '-.'],
    'MARKER_T': ['', '', ''],
    'MARKER_EVERY': [20, 20, 20],
    'LABEL':  ["$v=0.5$, $\sigma=0.0$", "$v=1.0$, $\sigma=0.0$", "$v=2.0$, $\sigma=0.0$"],
}

SETTING_3 ={
    'MAX_TRANS': 5,
    'LOG_DIR':  'logs',
    'SUB_DIR':  'shadowing',
    'ANA_DIR':  'analytical',
    'BACKOFF': 36,
    'MU_FADING': 0.0,
    'SIGMA_S':  3.0,
    'LINEWIDTH': 2,
    'LM_PAIR':  [(1, 2), (1, 1), (2, 1)],
    'LINE_COLOR': ['r', 'b', 'g'],
    'LINE_TYPE': ['--', '-', '-.'],
    'MARKER_T': ['v', 'o', '^'],
    'MARKER_EVERY': [20, 10, 10],
    'LABEL':  ["$v=0.5$, $\sigma=3.0$", "$v=1.0$, $\sigma=3.0$", "$v=2.0$, $\sigma=3.0$"],
}
#


SETTING = [SETTING_0, SETTING_3]

MU_FADING = 0.0
SIGMA_S = 1.0
LINEWIDTH = 2
LM_Pair = [(1, 1), (2, 1), (1, 2)]

METRICS = ['Packet Loss Rate', 'Throughput', 'Energy Efficiency', 'Expected Nb. of Transmissions']
Y_SCALE = ['log', 'linear', 'linear', 'linear']
SUBFIG_ENUM = ['a', 'b', 'c', 'd']


