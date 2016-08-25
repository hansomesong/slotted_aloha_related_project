__author__ = 'qsong'


import matplotlib.pyplot as plt
import numpy as np
import scipy.special

params = {
    'legend.fontsize': 20,
    'xtick.labelsize': 'large',
    'ytick.labelsize': 'large',
    'axes.labelsize': 30,
}
plt.rcParams.update(params)

print scipy.special.erf(0)
x = np.linspace(-10, 10, 1000)
y_1 = [1.0/(1 + np.exp(-e)) for e in x]
y_2 =[ 0.5+ 0.5*scipy.special.erf(np.sqrt(np.pi)/4.0*e) for e in x]

print 1.0/(1 + np.exp(-3)) - 0.5 - 0.5*scipy.special.erf(np.sqrt(np.pi)/4.0*3)

plt.figure(1, figsize=(12, 10))

plt.plot(x, y_1, 'r--', label=r"$y = \frac{1}{1+e^{-x}}$")
plt.plot(x, y_2, 'b-', label=r"$y = 0.5 + 0.5erf(\frac{\sqrt{\pi}}{4}x)$")

plt.legend(loc='lower right', numpoints=2)
plt.ylabel(r"y")
plt.xlabel(r"x")
plt.savefig('erf_approx.eps', format='eps', dpi=300)
plt.show()