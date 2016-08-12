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
SUB_DIR = 'fading_shadowing'

all_logs = glob.glob(os.path.join(LOG_DIR, SUB_DIR, "*.csv"))

f_name = "simd=5000_N=500_threshold=3.0_l=1_m=1_backoff=150_start=0.8_end=0.8_simstep=0.02_20160812040142.csv"

for csv_file in all_logs:

    csv_df = pd.read_csv(csv_file, sep=',', header=None)

    plr = csv_df.values[:, -2]

    alpha = csv_df.values[:, -1][0]

    ci =st.t.interval(0.95, len(plr)-1, loc=np.mean(plr), scale=st.sem(plr))

    print alpha, np.mean(plr), ci