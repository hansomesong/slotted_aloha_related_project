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

if __name__ == "__main__":
    LINEWIDTH = 2.0

    gamma = 4.0
    thetha_dB = 3.0
    sigma_dB = 8.0

    p_max_lb = 0.01
    p_max_ub = 0.1
    p_max_nb = 100
    p_max_vector = np.linspace(p_max_lb, p_max_ub, p_max_nb)

    slot_l_b_max_vector = sgam.best_bs_max_load(gamma, thetha_dB, sigma_dB, p_max_vector, pure=False, itf_mean=True)
    pure_max_l_b_max_vector = sgam.best_bs_max_load(gamma, thetha_dB, sigma_dB, p_max_vector, pure=True, itf_mean=False)
    pure_max_l_m_max_vector = sgam.div_max_load(gamma, thetha_dB, p_max_vector, pure=True, itf_mean=False)



    print
    fig, axes = plt.subplots(1, 1)

    axes.plot(
            p_max_vector,
            pure_max_l_m_max_vector,
            color='r',  marker='', linestyle='-', linewidth=LINEWIDTH, label="Diversity, pure-ALOHA, max.itf, theta=3.0 dB"
    )

    axes.plot(
            p_max_vector,
            pure_max_l_b_max_vector,
            color='b',  marker='', linestyle='-', linewidth=LINEWIDTH, label="Best ATT, pure-ALOHA, max.itf, theta=3.0 dB"
    )

    axes.plot(
            p_max_vector,
            slot_l_b_max_vector,
            color='g',  marker='', linestyle='-', linewidth=LINEWIDTH, label="Best ATT, slotted, theta=3.0 dB"
    )

    thetha_dB = 6.0
    slot_l_b_max_vector = sgam.best_bs_max_load(gamma, thetha_dB, sigma_dB, p_max_vector, pure=False, itf_mean=True)
    pure_max_l_b_max_vector = sgam.best_bs_max_load(gamma, thetha_dB, sigma_dB, p_max_vector, pure=True, itf_mean=False)
    pure_max_l_m_max_vector = sgam.div_max_load(gamma, thetha_dB, p_max_vector, pure=True, itf_mean=False)

    print sgam.bs_best_atch_op(lambda_m=3.0, lambda_b=0.08, gamma=4, p=0.008, thetha_dB=3, sigma_dB=0, pure=False, itf_mean=True)
    print 0.2*(1-sgam.bs_best_atch_op(lambda_m=3.0, lambda_b=0.08, gamma=4, p=0.008, thetha_dB=3, sigma_dB=0, pure=False, itf_mean=True))

    print sgam.bs_rx_div_op(lambda_m=3.0, lambda_b=0.08, gamma=4, p=0.008, thetha_dB=3, sigma_dB=0, pure=True, itf_mean=False)
    print 0.2*(1-sgam.bs_rx_div_op(lambda_m=3.0, lambda_b=0.08, gamma=4, p=0.008, thetha_dB=3, sigma_dB=0, pure=True, itf_mean=False))

    axes.plot(
            p_max_vector,
            pure_max_l_m_max_vector,
            color='r',  marker='', linestyle='--', linewidth=LINEWIDTH, label="Diversity, pure-ALOHA, max.itf, theta=6.0 dB"
    )

    axes.plot(
            p_max_vector,
            pure_max_l_b_max_vector,
            color='b',  marker='', linestyle='--', linewidth=LINEWIDTH, label="Best ATT, pure-ALOHA, max.itf, theta=6.0 dB"
    )

    axes.plot(
            p_max_vector,
            slot_l_b_max_vector,
            color='g',  marker='', linestyle='--', linewidth=LINEWIDTH, label="Best ATT, slotted, theta=6.0 dB"
    )
    axes.set_xticks(np.arange(0.01, 0.1, 0.01))
    axes.set_yticks(np.arange(0.0, 0.1, 0.01))

    axes.axis([0.01, 0.1, 0, 0.1])
    plt.legend(loc='best')


    axes.grid()


    plt.show()





