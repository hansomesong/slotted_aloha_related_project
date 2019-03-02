# -*- coding: utf-8 -*-
__author__ = 'qsong'

import pandas as pd
import numpy as np
import glob
import os

MAX_TRANS = 5
LOG_DIR = 'logs'
SUB_DIR = 'shadowing'
ANA_DIR = 'analytical'
SINR_THLD = [3.0]


def sim_data_process(sinr_thrsld, l, m, sigma_shadowing, backoff):
    SUBSUB_DIR = "backoff_{0}".format(backoff)
    CASE_DIR = 'case_{0}dB'.format(int(sinr_thrsld))
    POWER_DIR = "l_{0}_m_{1}_sigma_s_{2}".format(l, m, sigma_shadowing)
    logs_dir = glob.glob(os.path.join('.', LOG_DIR, SUB_DIR, CASE_DIR, SUBSUB_DIR, POWER_DIR, "*.csv"))
    #TODO: to be deleted
    print os.path.join('.', LOG_DIR, SUB_DIR, CASE_DIR, SUBSUB_DIR, POWER_DIR, "*.csv"), logs_dir

    sim_intensity = []
    sim_plr =  []
    sim_thrpt = []
    sim_ee = []
    sim_avg_nb = []
    for csv_file in logs_dir:
        csv_df = pd.read_csv(csv_file, sep=',', header=None)
        plr = csv_df.values[:, -2]
        alpha = csv_df.values[:, -1][0]
        ee_l  = csv_df.values[:, 6]
        prob_vector_df = csv_df.iloc[:, MAX_TRANS+2:-1]
        avg_nb = prob_vector_df.apply(lambda x: sum(x.values[0:-1]), axis=1).mean()
        sim_avg_nb.append(avg_nb)

        sim_intensity.append(alpha)
        # plt = 0.0
        if len(plr) > 10:
            plt = (np.sum(plr) - np.max(plr) - np.min(plr))/(len(plr)-2.0)
            ee = (np.sum(ee_l) - np.max(ee_l) - np.min(ee_l))/(len(ee_l)-2.0)
        else:
            plt = np.sum(plr)/len(plr)
            ee =  np.sum(ee_l)/len(ee_l)

        sim_plr.append(plt)
        sim_thrpt.append(alpha*(1-plt))
        sim_ee.append(ee)

        # print POWER_DIR, alpha, plt, alpha*(1-plt), ee, avg_nb
    print  l, m, sigma_shadowing, backoff
    print sim_intensity
    print sim_plr
    return sim_intensity, sim_plr, sim_thrpt, sim_ee, sim_avg_nb

if __name__ == "__main__":

    sinr_thrsld, l, m, sigma_shadowing, backoff = 3.0, 1, 2, 1, 36
    sim_data_process(sinr_thrsld, l, m, sigma_shadowing, backoff)

    table_labels = ["charges", "p_ana", "p_sim", "ratio"]

    p_sim = np.array([0.01064, 0.0156, 0.02401, 0.03549, 0.05356, 0.0852, 0.1215, 0.1756, 0.2287, 0.2992, 0.3782])
    p_ana = np.array([0.008151, 0.01301, 0.02068, 0.03284, 0.05213, 0.08222, 0.1269, 0.1868, 0.2577, 0.3322, 0.4044])

    p_ana_l_2_m_1 = np.array([0.001083, 0.002361, 0.005206, 0.01144, 0.02437, 0.04409, 0.08453, 0.1313, 0.1838, 0.2381, 0.2914])
    p_sim_l_2_m_1 = np.array([0.002117, 0.004768, 0.007510, 0.01650, 0.02920, 0.05796, 0.08329, 0.1207, 0.1650, 0.2064, 0.2682])

    p_ana_l_1_m_2 = np.array([0.02233, 0.02862, 0.03623, 0.04536, 0.05616, 0.06883, 0.08353, 0.1004, 0.1196, 0.1411, 0.1652])
    p_sim_l_1_m_2 = np.array([0.02202, 0.02787, 0.03596, 0.04406, 0.05361, 0.06954, 0.08443, 0.09885, 0.1197, 0.1373, 0.1594])

    print 10*np.log10(p_sim/p_ana)

    print 10*np.log10(p_sim_l_2_m_1/p_ana_l_2_m_1)

    print 10*np.log10(p_sim_l_1_m_2/p_ana_l_1_m_2)
