# -*- coding: utf-8 -*-
# The name of this script, sga is short for "stochastic geometry analytical model".
# A short name is convinient for import usage elsewhere
__author__ = 'qsong'

import numpy as np
import pandas as pd
from scipy.special import gamma as gamma_f
import scipy.stats as st


def bs_nearest_atch_op(lambda_m, lambda_b, gamma, p, thetha_dB, sigma_dB, pure=False, itf_mean=True):
    """
    This method is used to calculate the outage probability for the case where each device should attach to the nearest
    BS(i.e., Base Station) before transmission. The path-loss function is in form of simple law and has a singularity
    at the origin. The channel randomness effect we take into account fading and large-scale shadowing.
    :param lambda_m:    numpy array, spatial device density array
    :param lambda_b:    scalar, spatial BS density
    :param gamma:       scalar,  path-loss exponent, if gamma=0, no shadowing effect
    :param p:           scalar, the probability to transmit one message
    :param thetha_dB:   scalar, capture effect SINR threshold, unit dB!
    :param sigma_dB:    scalar, standard error of shadowing effect, unit dB
    :param pure:        boolean, False refers to slotted Aloha otherwise pure Aloha (i.e., non-slotted Aloha)
    :return: numpy array, the outage probability as function of lambda_m
    """
    THETA = np.power(10, thetha_dB/10)
    BETA = np.log(10.0)/10.0
    sigma = BETA*sigma_dB
    sigma_X = 2.0*sigma/gamma
    if not pure:
        k = np.pi*gamma_f(1-2.0/gamma)
    else:
        if itf_mean:
            factor = 2*gamma/(2+gamma)
        else:
            factor = 2.0
        k = np.pi*gamma_f(1-2.0/gamma)*factor

    # Fractional moment of fading effect
    fm_fading = gamma_f(1+2.0/gamma)
    # Fractional moment of shadowing effect
    fm_shadowing = np.exp(0.5*sigma_X**2)

    constant_A = p*lambda_m*k*fm_fading*fm_shadowing*np.power(THETA, 2.0/gamma)
    constant_B = constant_A/(lambda_b*np.pi)
    B_power = np.power(1+np.pi*sigma_X**2/8.0, -0.5)

    return 1.0 - 1.0/(1+np.power(constant_B, B_power))

def bs_rx_div_op(lambda_m, lambda_b, gamma, p, thetha_dB, sigma_dB, pure=False, itf_mean=True):
    """
    This method is used to calculate the outage probability for the case where each device employs a 'fire and forget'
    emission strategy. All the BSes are the potential receiver for the considered device.
    :param lambda_m:    numpy array, spatial device density array
    :param lambda_b:    scalar, spatial BS density
    :param gamma:       scalar,  path-loss exponent
    :param p:           scalar, the probability to transmit one message
    :param thetha_dB:   scalar, capture effect SINR threshold, unit dB!
    :param sigma_dB:    scalar, standard error of shadowing effect, unit dB
    :param pure:        boolean, False refers to slotted Aloha otherwise pure Aloha (i.e., non-slotted Aloha)
    :return: numpy array, the outage probability as function of lambda_m
    """
    THETA = np.power(10, thetha_dB/10)
    BETA = np.log(10.0)/10.0
    sigma = BETA*sigma_dB
    sigma_X = 2.0*sigma/gamma
    variant_1 = p*lambda_m/lambda_b
    if not pure:
        k = np.pi*gamma_f(1-2.0/gamma)
    else:
        if itf_mean:
            factor = 2*gamma/(2+gamma)
        else:
            factor = 2.0
        k = np.pi*gamma_f(1-2.0/gamma)*factor

    # Fractional moment of fading effect
    fm_fading = gamma_f(1+2.0/gamma)
    # Fractional moment of shadowing effect
    fm_shadowing = np.exp(0.5*sigma_X**2)

    # In this case, shadowing effect related term is not present in the constant_A...
    constant_A = p*lambda_m*k*fm_fading*np.power(THETA, 2.0/gamma)
    constant_B = constant_A/(lambda_b*np.pi)
    # B_power = np.power(1+np.pi*sigma_X**2/8.0, -0.5)
    # return np.exp(1.0/(1+np.power(constant_B, B_power))-np.exp(0.5*(2.0*sigma/gamma)**2)/constant_B)
    return np.exp(-1.0/constant_B)

def bs_nearest_atch_thrpt(lambda_m, lambda_b, gamma, p, thetha_dB, sigma_dB, pure=False, itf_mean=True):
    """
    This method is used to calculate the outage probability for the case where each device should attach to the nearest
    BS(i.e., Base Station) before transmission. The path-loss function is in form of simple law and has a singularity
    at the origin. The channel randomness effect we take into account fading and large-scale shadowing.
    :param lambda_m:    numpy array, spatial device density array
    :param lambda_b:    scalar, spatial BS density
    :param gamma:       scalar,  path-loss exponent, if gamma=0, no shadowing effect
    :param p:           scalar, the probability to transmit one message
    :param thetha_dB:   scalar, capture effect SINR threshold, unit dB!
    :param sigma_dB:    scalar, standard error of shadowing effect, unit dB
    :param pure:        boolean, False refers to slotted Aloha otherwise pure Aloha (i.e., non-slotted Aloha)
    :return: numpy array, the outage probability as function of lambda_m
    """
    THETA = np.power(10, thetha_dB/10)
    BETA = np.log(10.0)/10.0
    sigma = BETA*sigma_dB
    sigma_X = 2.0*sigma/gamma
    if not pure:
        k = np.pi*gamma_f(1-2.0/gamma)
    else:
        if itf_mean:
            factor = 2*gamma/(2+gamma)
        else:
            factor = 2.0
        k = np.pi*gamma_f(1-2.0/gamma)*factor

    # Fractional moment of fading effect
    fm_fading = gamma_f(1+2.0/gamma)
    # Fractional moment of shadowing effect
    fm_shadowing = np.exp(0.5*sigma_X**2)

    constant_A = p*lambda_m*k*fm_fading*fm_shadowing*np.power(THETA, 2.0/gamma)
    constant_B = constant_A/(lambda_b*np.pi)
    B_power = np.power(1+np.pi*sigma_X**2/8.0, -0.5)

    return p*lambda_m*(1.0/(1+np.power(constant_B, B_power)))/lambda_b


def bs_rx_div_thrpt(lambda_m, lambda_b, gamma, p, thetha_dB, sigma_dB, pure=False, itf_mean=True):

    return p*lambda_m*(1-bs_rx_div_op(lambda_m, lambda_b, gamma, p, thetha_dB, sigma_dB, pure, itf_mean))/lambda_b

def sim_data_process(folder_dir):
    sim_intensity = []
    sim_plr = []
    for csv_file in folder_dir:
        csv_df = pd.read_csv(csv_file, sep=',', skiprows=[0], header=None, dtype=np.float64)
        # Hard code is not good...
        plr_df = csv_df.values[1:, -2]
        # sort the obtained to remove the max and min
        plr_df.sort()
        # remove the max and min
        # plr_df = plr_df[1:-1]
        avg_plr = plr_df.mean()
        # Calculates the standard error of the mean (or standard error of measurement)
        # Here we have to be careful to keep all y values positive:
        ci_min_plr = max(avg_plr-1.96*st.sem(plr_df), 1e-7)
        ci_max_plr = avg_plr + 1.96*st.sem(plr_df)
        alpha = csv_df.values[:, 1][0]
        sim_intensity.append(alpha)
        sim_plr.append([ci_min_plr, avg_plr, ci_max_plr])
    return sim_intensity, sim_plr




