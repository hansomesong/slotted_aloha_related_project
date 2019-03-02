__author__ = 'qsong'

import numpy as np
from scipy.special import gamma as gamma_f
from scipy.special import erf as erf

from analytical_model import sgam
import matplotlib.pyplot as plt




gamma = 4.0

r = 1.0

thetha_dB = 3.0
threshold = np.power(10, thetha_dB/10)
intensity_b = 0.08
p = 0.008
intensity_m = 0.01
sigma_dB = 8.0
BETA = np.log(10)/10
sigma = sigma_dB * BETA
sigma_X = 2.0/gamma*sigma
length = 30000
chi = np.random.lognormal(0, sigma, length)
A = gamma_f(1+2.0/gamma) * np.exp(0.5*np.power(sigma_X, 2)) * np.power(threshold, 2.0/gamma)
K = np.pi * gamma_f(1-2.0/gamma)



B = A*K*p*intensity_m/(np.pi*intensity_b)


def psr(x, chi):
    # x, the distance
    # chi, np.power()
    return np.exp(-p*intensity_m * A * K * x**2 * (np.exp(-chi))**(2.0/gamma))

def pdf(x):
    """
    pdf is actually the probability density function of a RV of normal distribution. Its standard deviation error is
    noted as \sigma_X and mean as \log(B), where B = A * K * p * \lambda_m / (pi * \lambda_b ),
    A = gamma( 1 + 2.0/\gamma) * np.exp(0.5 \sigma_X ^ 2) * \theta ^{ 2.0/\gamma},
    K depends on slotted or pure ALOHA. In slotted ALOHA, k = pi * gamma(1-2.0/\gamma)
    :param x:           a vector
    :return:
    """
    tmp1 = np.exp(-0.5 * np.power((x - np.log(B)) / sigma_X, 2))
    tmp2 = np.power(np.sqrt(2 * np.pi) * sigma_X, -1.0)
    return tmp2*tmp1

def logistical(x):
    return np.power(1+np.exp(x),  -1)

def logistical2(x):
    return np.power(1+np.exp(-x),  -1)

def q(x):
    return 0.5 + 0.5*erf(x/1.595769/np.sqrt(2))

def logist(x):
    return np.power(1+np.exp(-x),  -1)


def repara_logist(x):
    gamma_3 = np.sqrt(1+np.pi*sigma_X**2/8)
    return np.power(1 + np.exp(-x/gamma_3), -1)

def repara_logist2(x):
    gamma_3 = -np.power(1+np.pi*sigma_X**2/8, -0.5)
    return np.power(1 + x**gamma_3, -1)

f_chi = np.power(B*np.power(chi, -2.0/gamma) + 1, -1.0)

p_f_bs_nst_att_8 = sgam.bs_nearest_atch_op2(intensity_m, intensity_b, gamma, p, thetha_dB, 8)

x = np.linspace(-100, 100, 100000)




y = logistical(x)*pdf(x)

plt.plot(x, np.power(logistical(x), 2), 'b')
plt.plot(x, logistical2(x), 'r')

gamma_2 = np.sqrt(1+np.pi*sigma_X**2/8)

print "K", K, "B", B
print "p_{f,n}=", 1-np.mean(f_chi), p_f_bs_nst_att_8
print "p_{f,n}, trap, ",  1-np.trapz(y, x)
print "p_{f,n}, trap, ",  np.trapz(logistical2(x)*pdf(x), x)
print "p_{f,n}, trap, logis-normal-app ",  np.trapz(q(x)*pdf(x), x)
print "p_{f,n}, trap, logis-normal-app-further ",  0.5 + 0.5*erf(np.log(B)*np.sqrt(np.pi)/(4*gamma_2))
print "p_{f,n}, repara",  repara_logist(np.log(B))
print "p_{f,n}, repara2",  repara_logist2(B)
print "p_{f,n}, lower bound",  1-np.power(1+np.exp(np.log(B) + 0.5*sigma_X), -1)

print "p_{f,n}, upper bound",  1-np.power(1+np.exp(np.log(B) - 0.5*sigma_X), -1)

print "p_{s}(r)",  np.mean(psr(1, chi))
print "p_{f}(best)",  1 - np.power(B+1, -1)


print "p_{f}(best)",  1 - np.power(B+1, -1)

x = np.linspace(0, 2, 100)

y1 = np.trapz((1+np.exp(x))**-2*pdf(x), x)
y2 = np.trapz((1+3*np.exp(x))**-1*pdf(x), x)


print sgam.bs_rx_div_op(1.8, intensity_b, gamma, p, thetha_dB, 8, True)













