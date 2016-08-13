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
SUB_DIR = 'shadowing'
SUBSUB_DIR = "backoff_50"
SUBSUBSUB_DIR = "l_1_m_1_sigma_s_1"

all_logs = glob.glob(os.path.join(LOG_DIR, SUB_DIR, SUBSUB_DIR, SUBSUBSUB_DIR, "*.csv"))
# all_logs = [os.path.join(LOG_DIR,
#                  "simd=5000_N=500_threshold=3.0_l=1_m=1_backoff=50_start=0.4_end=0.4_simstep=0.02_sigmas=1.0_20160813014646.csv")
# ]
for csv_file in all_logs:

    csv_df = pd.read_csv(csv_file, sep=',', header=None)

    plr = csv_df.values[:, -2]

    alpha = csv_df.values[:, -1][0]

    ci =st.t.interval(0.95, len(plr)-1, loc=np.mean(plr), scale=st.sem(plr))

    print alpha, np.mean(plr), ci