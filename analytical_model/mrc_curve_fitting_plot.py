# -*- coding: utf-8 -*-
# The curve fit function is here...
__author__ = 'qsong'

from analytical_model import mrc_curve_fitting
import numpy as np
import pandas as pd
import os
import glob
import re
import matplotlib
import matplotlib.pyplot as plt
from scipy.optimize import leastsq

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

params = {
    'legend.fontsize': 15,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'axes.labelsize': 10,
    'legend.numpoints': 1
}
plt.rcParams.update(params)


K_labels = ['K_slot',  'K_pure_mean',  'K_pure_max']
B_labels = ['B_slot',  'B_pure_mean',  'B_pure_max']

pure = [False, True, True]
itf_mean = [True, True, False]
title = ['Slotted ALOHA', 'Pure ALOHA, mean', 'Pure ALOHA, max']

# The default ploynomial fitting is 3 order.
POLY_FIT_DEGREE = 3

if __name__ == "__main__":


    # Normally, this flag is set as false
    ana_need_process_again = False
    sim_need_process_again = True

    thetha_dB = 3.0
    ANA_FIT_CSV_NAME = './ana_fit_result_theta_{0}.csv'.format(int(thetha_dB))
    FIG_DST = "/Users/qsong/Documents/Communication_Letter_MRC_curve_fitting_20180721"
    if not os.path.exists(FIG_DST):
        os.makedirs(FIG_DST)

    FIG_NAME1 = 'ana_curve_fit_K_B.eps'
    FIG_NAME2 = 'sim_curve_fit_K_B.eps'


    x_interval = [[0.15, 0.3], [0.12, 0.22], [0.08, 0.15]]

    ana_gammas = np.arange(3.0, 6.1, 0.1)
    # I have tried that for any ALOHA type, B is always the same.
    # K, B depends only on gamam (instead of thetha_T)
    ana_column_name = [
        'gamma',
        'K_slot', 'B_slot', 'slot_goodness',
        'K_pure_mean', 'B_pure_mean', 'pure_mean_goodness',
        'K_pure_max', 'B_pure_max', 'pure_max_goodness'
    ]
    #
    if not os.path.isfile(ANA_FIT_CSV_NAME) or os.stat(ANA_FIT_CSV_NAME).st_size == 0 or ana_need_process_again:
        ana_table_content = []
        for gamma in ana_gammas:
            ana_tmp = [gamma]
            for i in range(3):
                loads, plrs, X_ref, Y_ref, K, B, fit_func, goodness = mrc_curve_fitting.curve_fitting(
                    thetha_dB, gamma, pure[i], itf_mean[i], plr_min=0.001, plr_max=0.1
                )
                print "ALOHA type:", title[i], "\t\tgamma: ", gamma, ",K and B: ", K, ', ', B, "\tgoodness: ", goodness
                ana_tmp.append(K)
                # when the loop for aloha type finishes, add the value of B at the end
                ana_tmp.append(B)
                ana_tmp.append(goodness)
            # Add the obtained list, [K_slot, K_pure_mean, K_pure_max, B] to the table content
            ana_table_content.append(ana_tmp)

        # convert table_content into a pandas dataframe
        ana_table_content = pd.DataFrame(ana_table_content, columns=ana_column_name)

        print ana_table_content
        ana_table_content.to_csv(ANA_FIT_CSV_NAME, index=False, sep=',')

    ana_table_content = 0
    if os.path.isfile(ANA_FIT_CSV_NAME) and os.stat(ANA_FIT_CSV_NAME).st_size != 0:
        ana_table_content = pd.read_csv(ANA_FIT_CSV_NAME)
        print "ana_table_content"
        print ana_table_content


        def func(para, x):
            k, b = para
            return k*x + b
            ## 偏差函数：x,y都是列表:这里的x,y更上面的Xi,Yi中是一一对应的
        def error(para, x, y):
            return func(para, x) - y


        gammas = np.array(ana_table_content.gamma.tolist())
        print "gammas", gammas


        fig, axes = plt.subplots(3, 2)

        for i in range(3):
            K_list = np.array(ana_table_content.loc[:, K_labels[i]].tolist())
            B_list = np.array(ana_table_content.loc[:, B_labels[i]].tolist())
            # print K_list
            # print B_list
            # It seems that leastsq function only accept np array! with list there will be error!! (not running error)
            K_fit_para1 = np.polyfit(gammas, K_list, POLY_FIT_DEGREE)
            B_fit_para1 = np.polyfit(gammas, B_list, POLY_FIT_DEGREE)
            print "Polyfit K: ", K_fit_para1
            print "Polyfit B: ", B_fit_para1
            print K_labels[i], B_labels[i]


            K_fit1 = 0
            B_fit1 = 0
            # since index i is used at outter level loop. Use other index for this level loop
            for j, coeff in enumerate(K_fit_para1):
                # print j, coeff
                K_fit1 += coeff*np.power(gammas, 3-j)
            for j, coeff in enumerate(B_fit_para1):
                B_fit1 += coeff*np.power(gammas, 3-j)


            Y_LABELS =['K', 'B']

            TITLES = ['Curve fitting for parameter K', 'Curve fitting for parameter B']

            K_fit_std1 = np.sqrt(np.sum(np.power((K_list - K_fit1)/K_fit1, 2.0))/K_list.size)
            K_fit_std2 = np.sqrt(np.sum(np.power((K_list - K_fit1)/K_list, 2.0))/K_list.size)

            B_fit_std1 = np.sqrt(np.sum(np.power((B_list - B_fit1)/B_list, 2.0))/B_list.size)
            B_fit_std2 = np.sqrt(np.sum(np.power((B_list - B_fit1)/B_fit1, 2.0))/B_list.size)
            print np.power((B_list - B_fit1)/B_list, 2.0)
            print "K estimator goodness: ", K_fit_std1, K_fit_std2, "and B estimator goodness: ", B_fit_std1, B_fit_std2


            axes[i, 0].scatter(gammas, K_list, c='r')
            # ax.plot(gammas, K_fit, 'b')
            axes[i, 0].plot(gammas, K_fit1, 'g')
            axes[i, 0].set_xlabel(r"gamma")
            axes[i, 0].set_ylabel(K_labels[i])
            axes[i, 0].grid()

            axes[i, 1].scatter(gammas, B_list, c='r')
            # ax.plot(gammas, K_fit, 'b')
            axes[i, 1].plot(gammas, B_fit1, 'g')
            axes[i, 1].set_xlabel(r"gamma")
            axes[i, 1].set_ylabel(B_labels[i])
            axes[i, 1].grid()

        plt.savefig(os.path.join(FIG_DST, FIG_NAME1), format='eps', dpi=300)



            # for j, ax in enumerate(axes[i]):
            #     if j == 0:
            #         ax.scatter(gammas, K_list, c='r')
            #         # ax.plot(gammas, K_fit, 'b')
            #         ax.plot(gammas, K_fit1, 'g')
            #         ax.set_xlabel(r"gamma")
            #         ax.set_ylabel(Y_LABELS[j])
            #     if j == 1:
            #         ax.scatter(gammas, B_list, c='r')
            #         # ax.plot(gammas, B_fit, 'b')
            #         ax.plot(gammas, B_fit1, 'g')
            #         ax.set_xlabel(r"gamma")
            #         ax.set_ylabel(Y_LABELS[j])




            # print "K fit error: ", np.array(K_list) - np.array(K_fit)
            # print "B fit error: ", np.array(B_list) - np.array(B_fit)

    SIM_LOG_DIR = '/Users/qsong/Documents/slotted_aloha_related_project/logs/SimuNov23PlusdePoints'
    thetha_dB = 3
    SIM_FIT_CSV_NAME = './sim_fit_result_theta_{0}.csv'.format(int(thetha_dB))

    sim_table_content = []

    if not os.path.isfile(SIM_FIT_CSV_NAME) or os.stat(SIM_FIT_CSV_NAME).st_size == 0 or sim_need_process_again:

        sim_column_name = [
            'gamma',
            'K_slot', 'B_slot', 'slot_goodness',
            'K_pure_mean', 'B_pure_mean', 'pure_mean_goodness',
            'K_pure_max', 'B_pure_max', 'pure_max_goodness'
        ]
        for filename in glob.iglob(os.path.join(SIM_LOG_DIR, 'Aloha*.csv')):
            sim_tmp = []
            pattern = re.compile(r"AlohaMultg(\d+)s(\d+)t(\d+)")
            simu_setting = re.findall(pattern, filename)[0]
            label = simu_setting[0]
            gamma = float(simu_setting[0])/10.0
            sigma_dB = float(simu_setting[1])
            theta_dB = float(simu_setting[2])
            # print gamma, sigma_dB, theta_dB
            sim_tmp.append(gamma)
            for i in range(3):
                loads, plrs, X_ref, Y_ref, K, B, fit_func, goodness =\
                    mrc_curve_fitting.sim_curve_fitting(
                        filename, thetha_dB, sigma_dB, gamma, pure[i], itf_mean[i], plr_min=0.01, plr_max=0.1
                    )
                sim_tmp.append(K)
                sim_tmp.append(B)
                sim_tmp.append(goodness)
                print "ALOHA type:", title[i], "\t\tgamma: ", gamma, ",K and B: ", K, ', ', B, "\tgoodness: ", goodness
            sim_table_content.append(sim_tmp)
            # convert table_content into a pandas dataframe

        sim_table_content = pd.DataFrame(sim_table_content, columns=sim_column_name).sort_values(by=['gamma'])
        # print sim_table_content
        sim_table_content.to_csv(SIM_FIT_CSV_NAME, index=False, sep=',')


    if os.path.isfile(SIM_FIT_CSV_NAME) and os.stat(SIM_FIT_CSV_NAME).st_size != 0:
        sim_table_content = pd.read_csv(SIM_FIT_CSV_NAME)
        print "sim_table_content"
        print sim_table_content

        # Filter gamma.
        # We sometimes need a narrower range of gamma than that present in simulation.
        # For example, simulation resulats may cover a range [3.0, 6.0], but we only care about [3.3, 5.0]
        # filter_mask = sim_table_content['gamma'] < 6.0
        gammas = np.array(sim_table_content.gamma.tolist())


        fig, axes = plt.subplots(nrows=3, ncols=2, sharex=True)
        # fig.tight_layout()

        #
        k_y_lims = [[1.6, 1.9], [2.25, 2.40], [2.9, 3.5]]


        for i in range(3):
            print title[i]
            sim_K_list = np.array(sim_table_content[K_labels[i]].tolist())
            sim_B_list = np.array(sim_table_content[B_labels[i]].tolist())
            sim_K_fit_para1 = np.polyfit(gammas, sim_K_list, POLY_FIT_DEGREE)
            sim_B_fit_para1 = np.polyfit(gammas, sim_B_list, POLY_FIT_DEGREE)

            sim_K_fit1 = 0
            sim_B_fit1 = 0
            # since index i is used at outter level loop. Use other index for this level loop
            for j, coeff in enumerate(sim_K_fit_para1):
                # print j, coeff
                sim_K_fit1 += coeff*np.power(gammas, 3-j)
            for j, coeff in enumerate(sim_B_fit_para1):
                sim_B_fit1 += coeff*np.power(gammas, 3-j)

            sim_table_content.plot.scatter(ax=axes[i, 0], x='gamma', y=K_labels[i])
            sim_table_content.plot.scatter(ax=axes[i, 1], x='gamma', y=B_labels[i])
            axes[i, 0].plot(gammas, sim_K_fit1, label=K_labels[i])
            axes[i, 1].plot(gammas, sim_B_fit1, label=B_labels[i])
            axes[i, 1].set_xlim([3.0, 5.0])
            axes[i, 0].set_ylim(k_y_lims[i])
            axes[i, 0].grid()
            axes[i, 1].grid()

            K_fit_std1 = np.sqrt(np.sum(np.power((sim_K_list - sim_K_fit1)/sim_K_list, 2.0))/sim_K_list.size)
            B_fit_std1 = np.sqrt(np.sum(np.power((sim_B_list - sim_B_fit1)/sim_B_list, 2.0))/sim_B_list.size)

            print "Goodness of K", sim_K_fit_para1, " from simulation estimation: ", K_fit_std1
            print "Goodness of B", sim_B_fit_para1, " from simulation estimation: ", B_fit_std1


        plt.tight_layout(rect=[0, 0.03, 1, 0.95])
        plt.savefig(os.path.join(FIG_DST, FIG_NAME2), format='eps', dpi=300)

        # sim_fit_log = '/Users/qsong/Documents/slotted_aloha_related_project/analytical_model/sim_fit_result_theta_3.csv'
        # load = np.linspace(0.1, 0.5, 20)
        # mrc_curve_fitting.sim_fitted_function(sim_fit_log, thetha_dB, load, 3.5, True, True)

    # Do the second curve fitting for K and B

    filename = '/Users/qsong/Documents/slotted_aloha_related_project/logs/SimuNov23PlusdePoints/AlohaMultg60s8t3.csv'
    pattern = re.compile(r"AlohaMultg(\d+)s(\d+)t(\d+)")
    simu_setting = re.findall(pattern, filename)[0]
    label = simu_setting[0]
    gamma = float(simu_setting[0])/10.0
    sigma_dB = float(simu_setting[1])
    theta_dB = float(simu_setting[2])
    loads, plrs, X_ref, Y_ref, K, B, fit_func, goodness =\
                    mrc_curve_fitting.sim_curve_fitting(
                        filename, theta_dB, sigma_dB, gamma, pure[1], itf_mean[1], plr_min=0.001, plr_max=0.1
                    )
    # fig, axes = plt.subplots(1, 1, figsize=(8, 6))
    #
    # axes.scatter(X_ref, Y_ref)
    # axes.plot(X_ref, K*X_ref + B)
    #
    # print gamma, K, B, goodness
    #
    # plt.show()
    plt.show()








