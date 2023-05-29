############################################################
# imports

import numpy as np
import matplotlib.pyplot as plt
from demotools import source_counts

freq = '230'
BHMF = 'KH'

############################################################
# make a source counts plot

# create a grid in (theta,sigma)
logtheta_min = -1.0
logtheta_max = 2.0
logsigma_min = -3.0
logsigma_max = 0.0
theta1d = 10.0**np.linspace(logtheta_min,logtheta_max,100)
sigma1d = 10.0**np.linspace(logsigma_min,logsigma_max,100)
theta, sigma = np.meshgrid(theta1d,sigma1d,indexing='ij')

# compute the source counts on the grid
counts = source_counts(theta,sigma,freq=freq,counts_type='full',BHMF_type='KH')
logtheta = np.log10(theta)
logsigma = np.log10(sigma*1000.0)
logcounts = np.log10(counts)

# mask below some minimum source counts
logcounts_min = -2.0
mask = (logcounts <= logcounts_min)
counts_plot = np.ma.array(logcounts,mask=mask)

fig = plt.figure(figsize=(4,4))
ax = fig.add_axes([0.1,0.1,0.8,0.8])

# make a colorized contour plot
cf = ax.contourf(logtheta,logsigma,counts_plot,255,cmap='YlOrRd',vmin=logcounts_min)

# remove white lines between contour levels
for c in cf.collections:
    c.set_edgecolor("face")

# add contours
levels_neg = np.arange(logcounts_min,0.0)
levels_pos = np.arange(1.0,10.0)
cs1 = ax.contour(logtheta,logsigma,logcounts,levels=levels_neg,colors='k',linewidths=1,linestyles='--')
cs2 = ax.contour(logtheta,logsigma,logcounts,levels=[0.0],colors='k',linewidths=2,linestyles='-')
cs3 = ax.contour(logtheta,logsigma,logcounts,levels=levels_pos,colors='k',linewidths=1,linestyles='-')

# remaining plot items
ax.set_xlim(-1,2)
ax.set_ylim(0,3)
ax.set_xlabel(r'$\log(\theta/\mu\rm{as})$')
ax.set_ylabel(r'$\log(S_{\nu}/\rm{mJy})$')

# save
plt.savefig('source_counts_'+freq+'GHz.png',dpi=300,bbox_inches='tight')
plt.close()
