# -*- coding: utf-8 -*-
__author__ = 'qsong'

import numpy as np
import matplotlib.pyplot as plt
import scipy
import scipy.stats

intensity = 0.4
vector_p = [1.0, 0.60, 0.36, 0.22, 0.13]

intensities = [intensity*p for p in vector_p]

print sum(intensities)

mu = 0.0
sigma = 1.0
beta = 0.1*np.log(10.0)
rho = 1.0
v = 1.0
K = 4

# First order moment of Y
fm_Y = intensity*rho*np.exp(0.5*beta**2*sigma**2)*sum([v**k*p for k, p in enumerate(vector_p)])
fm_Y_plus = fm_Y/(1-np.exp(-sum(intensities)))
# Second order moment of Y
x = 0
for m in range(K):
    for n in range(m+1, K+1):
        x +=v**(m+n)*vector_p[m]*vector_p[n]


sm_Y = 2*(rho*intensity)**2*np.e**(beta**2*sigma**2)*x \
       + intensity*rho**2*np.exp(beta**2*sigma**2)*\
         sum([v**(2*k)*vector_p[k]*(np.exp(beta**2*sigma**2) + intensity*vector_p[k]) for k in range(K)])
sm_Y_plus = sm_Y/(1-np.exp(-sum(intensities)))

mu_eq = 2*np.log(fm_Y_plus) - 0.5*np.log(sm_Y_plus)
sigma_eq = np.sqrt(np.log(sm_Y_plus) - 2*np.log(fm_Y_plus))
print fm_Y
print sm_Y

print mu_eq
print sigma_eq


plt.figure(1)
mu, sigma = 0.0, 1.0 # mean and standard deviation
intensity = sum(vector_p)


nb = 4000000
z1 = [sum(np.random.lognormal(mu, sigma, max([np.random.poisson(intensity), 1]))) for i in range(nb)]
log_z = np.log(z1)
log_z_mean = np.mean(log_z)
log_z_std = np.std(log_z)
print "MIN:", min(z1)
print "MAX:", max(z1)
print "mean:", log_z_mean
print "stand error:", log_z_std

print "0:", sum([1 for x in z1 if x == 0])
# count, bins, ignored = plt.hist(z, 100, normed=True, align='mid')
# sns.distplot(z1, bins=[-0.1+0.05*x for x in range(600)]).set(xlim=(0, 30))
# sns.distplot(z1, 600, hist=True)


weights = np.ones_like(np.array(log_z))/float(len(z1))
count, bins, ignored = plt.hist(log_z, 100, normed=False, color='w')
plt.title(r"Histogram of $ln^{Z}$, with packet arrival rate 0.8", fontsize=30)

xmin, xmax = plt.xlim()
x = np.linspace(xmin, xmax, 1000)
mu, std = mu_eq, sigma_eq
pdf = scipy.stats.norm.pdf(x, mu, std)
plt.plot(x, pdf, linewidth=2, color='r', label=r"pdf of log(Z) with Mean:{:0.3f}, STD:{:1.3f}".format(mu_eq, sigma_eq))
plt.legend(loc='upper left')

# weights=weights,

plt.xlabel("Value of Z")
plt.ylabel("Frequency")
plt.show()


