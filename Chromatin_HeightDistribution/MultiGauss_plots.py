''' This script plots the results obtained from the Matlab script, in which the Chromatin height
distribution is fitted using a GaussianMixture model.'''


import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns


# Defining the gaussian function
def gaussian(x,mu,sigma):
    return 1/(sigma*(2*np.pi)**0.5)*np.exp(-0.5*((x-mu)/sigma)**2)

# loading data of the chromatin height distribution
data = np.loadtxt('20231117_ChromatinDistribution_points.txt', skiprows=1) # load data
print(data)
#  loading file containing the parameters defining the 5 gaussians obtained from the Matlab script
fit = np.loadtxt('five_gaussians.dat', delimiter=',')
# prefactor
H = fit[:, 0]
# means
mu = fit[:, 1]
# sigmas
sig = fit[:, 2]

# Plotting
fig, ax = plt.subplots()
ax.set_xlabel("Chromatin height (nm)")
sns.histplot(x=data, bins=100, stat='density', element="step", alpha=.5, ax=ax, color='darkgray')
# ax.hist(data, 100, density=True)
x = np.linspace(min(data), max(data), 200)
# Initializing the cumulative distribution of the n gaussians
y_tot = np.zeros(np.shape(x))
num_gauss = 5
for i in range(num_gauss):
    y = H[i] * gaussian(x, mu[i], sig[i] ** 0.5)
    y_tot = y_tot + y
    ax.plot(x, y, linewidth=2, linestyle='--')
plt.plot(x, y_tot, linewidth=2, linestyle='-', c='dimgray')
plt.tight_layout()
# saving the figure
# plt.savefig("MultiGauss_fit.pdf", dpi=500)
plt.show()
