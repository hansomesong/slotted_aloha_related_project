# -*- coding: utf-8 -*-
__author__ = 'qsong'

import numpy as np
import matplotlib.pyplot as plt
import scipy
import scipy.stats
from scipy.stats import norm


intensity = 0.8
alpha = intensity
vector_p = [1.,          0.55066884,  0.,          0.,          0.        ]
P = vector_p
intensities = [intensity*p for p in vector_p]

# print sum(intensities)

mu = 0.0
sigma = .0001
beta = 0.1*np.log(10.0)
BETA = 0.1*np.log(10.0)

v = 2.0
K = 4

NB = 40000
MAX_TRANS = 5
k = 0
SINR = 3.0

# Y = [
#     sum(
#         [v*(m-k)*sum(
#             [np.random.lognormal(mu, np.sqrt(2)*beta*sigma, np.random.poisson(intensity))]
#         )
#         for m, inensity in enumerate(intensities)]
#     )
#     for i in range(NB)
# ]
# Y = []
#
# for i in range(NB):
#     ln_rv = 0.0
#     for m, intensity in enumerate(intensities):
#         ln_rv += v**(m-k)*sum(np.random.lognormal(mu, np.sqrt(2)*beta*sigma, np.random.poisson(intensity)))
#
#     Y.append(ln_rv)
# print Y

Y = [sum([v**(m-k)*sum(np.random.lognormal(mu, np.sqrt(2)*beta*sigma, np.random.poisson(intensity))) for m, intensity in enumerate(intensities)]) for i in range(NB)]

# print Y

Y_plus = [element for element in Y if element > 0]
print Y_plus

mean_Y_plus = np.mean(Y_plus)
mean_Y_2_plus = np.mean([x**2 for x in Y_plus])

mu_Y_tmp = 2*np.log(mean_Y_plus) - 0.5*np.log(mean_Y_2_plus)
var_Y_tmp = np.log(mean_Y_2_plus) - 2*np.log(mean_Y_plus)

p_zero_inf = np.exp(-1.0*alpha*sum(P))
mean_Y = (1+BETA**2*sigma**2)*sum([np.power(v, m-k)*alpha*P[m] for m in range(K)])
var_Y = np.exp(np.power(2*BETA*sigma, 2))*sum([np.power(v, 2*(m-k))*alpha*P[m] for m in range(K)])
mean_Y_2 = np.power(mean_Y, 2) + var_Y

mean_Y_plus = mean_Y/(1-p_zero_inf)
mean_Y_2_plus = mean_Y_2/(1-p_zero_inf)


mu_eq = 2*np.log(mean_Y_plus) - 0.5*np.log(mean_Y_2_plus)
sigma_eq = np.log(mean_Y_2_plus) - 2*np.log(mean_Y_plus)


plt.figure(1)
log_z = np.log(Y_plus)
log_z_mean = np.mean(log_z)
log_z_std = np.std(log_z)
count, bins, ignored = plt.hist(log_z, 50, normed=True, color='b', label="Histogram of log(Z) with G = 0.8")
xmin, xmax = plt.xlim()
x = np.linspace(xmin, xmax, 1000)
mu, std = scipy.stats.norm.fit(log_z)
pdf = scipy.stats.norm.pdf(x, mu, std)

ERROR = (-BETA*SINR - mu_eq)/(np.sqrt(2*sigma_eq**2))


print mu_Y_tmp, var_Y_tmp, (1 - norm.cdf(-beta*SINR, mu_Y_tmp, np.sqrt(var_Y_tmp)))*(1-p_zero_inf), 0.5*(1-scipy.special.erf(ERROR))*(1-p_zero_inf)
print mu, std**2
print mu_eq, sigma_eq, (1 - norm.cdf(-beta*SINR, mu_eq, np.sqrt(sigma_eq)))*(1-p_zero_inf), 0.5*(1-scipy.special.erf(ERROR))*(1-p_zero_inf)

plt.plot(x, pdf, linewidth=2, color='r', label=r"pdf of log(Z) with Mean:{:0.3f}, STD:{:1.3f}".format(log_z_mean, log_z_std))
plt.legend(loc='upper left')
plt.show()







