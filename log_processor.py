# -*- coding: utf-8 -*-
__author__ = 'qsong'

import csv
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as st
import os
import pandas as pd
import glob

LOG_DIR = 'logs'

SUB_DIR = 'fading_shadowing_diff_ipc'
CASE_DIR = "case_0dB"
SUBSUB_DIR = "backoff_50"
SUBSUBSUB_DIR = "l_1_m_1_sigma_s_1"

all_logs = glob.glob(os.path.join(LOG_DIR, SUB_DIR, CASE_DIR, SUBSUB_DIR, SUBSUBSUB_DIR, "*.csv"))
all_logs = [
    os.path.join(
        LOG_DIR,
        "simd=5000_N=400_threshold=0.0dB_l=1_m=1_backoff=2_alpha=0.6_mufading=1.0_mushadowing=0.0_sigmashadowing=1.0_tmp=20160823143844.csv")
]
for csv_file in all_logs:

    csv_df = pd.read_csv(csv_file, sep=',', header=None)

    plr = csv_df.values[:, -2]

    alpha = csv_df.values[:, -1][0]

    ci =st.t.interval(0.95, len(plr)-1, loc=np.mean(plr), scale=st.sem(plr))

    print alpha, np.mean(plr), ci