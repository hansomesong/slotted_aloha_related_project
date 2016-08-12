# -*- coding: utf-8 -*-
__author__ = 'qsong'

import csv
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as st
import os
import pandas as pd
LOG_DIR = 'logs'
f_name = "simd=5000_N=500_threshold=3.0_l=1_m=1_backoff=150_start=0.44_end=0.44_simstep=0.02_032136.csv"

csv_file = os.path.join(
    LOG_DIR,
    f_name
)

csv_df = pd.read_csv(csv_file, sep=',', header=None)

plr = csv_df.values[:, -2]

ci =st.t.interval(0.95, len(plr)-1, loc=np.mean(plr), scale=st.sem(plr))

print np.mean(plr), ci