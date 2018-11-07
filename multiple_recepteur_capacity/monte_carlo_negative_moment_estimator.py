__author__ = 'qsong'


import numpy as np
import pandas as pd
from scipy.special import gamma as gamma_f
from scipy.special import erf as erf
import scipy.stats as st
import scipy.integrate as integrate


def deploy_nodes(width, alpha):
    """
    :param width:                   scalar, the radius of simulation area,
    :param alpha:                   scalar, the spatial MTC device density,
    :param intensity_bs:            scalar, the spatial BS density
    :return:
    """
    # The involved device number will be one sample from a spatial PPP over a finit disk area
    # Generate the needed device nb and base station in this simulation
    AREA_SURFACE = np.pi*np.power(width, 2)
    device_nb = int(np.random.poisson(alpha*AREA_SURFACE, 1))
    # Uniformelly distribute MTC devices using polar coordinates.
    # We assume that always one device at the origin
    device_rho = np.concatenate(([0.0], width*np.sqrt(np.random.uniform(0, 1, device_nb-1))))
    device_arguments = np.random.uniform(-np.pi, np.pi, device_nb)
    coordinates_devices_array = zip(device_rho, device_arguments)

    return device_nb, coordinates_devices_array