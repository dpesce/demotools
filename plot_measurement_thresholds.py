##################################################
# imports etc.

import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
import glob
import os
from demotools import source_counts

# output directory to save plots in
outdir = './survey'

##################################################
# read in survey data tables

# identify relevant directories
dirlist = np.sort(glob.glob('./survey/weather_inst*'))

# initialize arrays
F0_dum, ring_diam_dum, Delta_d, Delta_alpha, Delta_argb2, Delta_ampb1I, Delta_argb1I, Delta_ampb1, Delta_argb1, Delta_ampb2 = np.loadtxt(dirlist[0]+'/paramtable.txt',unpack=True)
F0_arr = np.sort(np.unique(F0_dum))
ring_diam_arr = np.sort(np.unique(ring_diam_dum))
Delta_d = Delta_d.reshape((len(F0_arr),len(ring_diam_arr))) * 0.0
Delta_alpha = Delta_alpha.reshape((len(F0_arr),len(ring_diam_arr))) * 0.0
Delta_argb2 = Delta_argb2.reshape((len(F0_arr),len(ring_diam_arr))) * 0.0
Delta_ampb1I = Delta_ampb1I.reshape((len(F0_arr),len(ring_diam_arr))) * 0.0
Delta_argb1I = Delta_argb1I.reshape((len(F0_arr),len(ring_diam_arr))) * 0.0
Delta_ampb1 = Delta_ampb1.reshape((len(F0_arr),len(ring_diam_arr))) * 0.0
Delta_argb1 = Delta_argb1.reshape((len(F0_arr),len(ring_diam_arr))) * 0.0
Delta_ampb2 = Delta_ampb2.reshape((len(F0_arr),len(ring_diam_arr))) * 0.0
F0_dum, ring_diam_dum, argb2, ampb1I, argb1I, ampb1, argb1, ampb2 = np.loadtxt(dirlist[0]+'/ring_paramtable.txt',unpack=True)
argb2 = argb2.reshape((len(F0_arr),len(ring_diam_arr))) * 0.0
ampb1I = ampb1I.reshape((len(F0_arr),len(ring_diam_arr))) * 0.0
argb1I = argb1I.reshape((len(F0_arr),len(ring_diam_arr))) * 0.0
ampb1 = ampb1.reshape((len(F0_arr),len(ring_diam_arr))) * 0.0
argb1 = argb1.reshape((len(F0_arr),len(ring_diam_arr))) * 0.0
ampb2 = ampb2.reshape((len(F0_arr),len(ring_diam_arr))) * 0.0

# loop through directories and read in files
count = 0.0
for indir in dirlist:

    filename_inputs = indir+'/paramtable.txt'
    filename_meas = indir+'/ring_paramtable.txt'
    if os.path.exists(filename_inputs) & os.path.exists(filename_meas):
        
        F0_dum, ring_diam_dum, Delta_d_here, Delta_alpha_here, Delta_argb2_here, Delta_ampb1I_here, Delta_argb1I_here, Delta_ampb1_here, Delta_argb1_here, Delta_ampb2_here = np.loadtxt(filename_inputs,unpack=True)
        Delta_d += Delta_d_here.reshape((len(F0_arr),len(ring_diam_arr)))
        Delta_alpha += Delta_alpha_here.reshape((len(F0_arr),len(ring_diam_arr)))
        Delta_argb2 += Delta_argb2_here.reshape((len(F0_arr),len(ring_diam_arr)))
        Delta_ampb1I += Delta_ampb1I_here.reshape((len(F0_arr),len(ring_diam_arr)))
        Delta_argb1I += Delta_argb1I_here.reshape((len(F0_arr),len(ring_diam_arr)))
        Delta_ampb1 += Delta_ampb1_here.reshape((len(F0_arr),len(ring_diam_arr)))
        Delta_argb1 += Delta_argb1_here.reshape((len(F0_arr),len(ring_diam_arr)))
        Delta_ampb2 += Delta_ampb2_here.reshape((len(F0_arr),len(ring_diam_arr)))

        F0_dum, ring_diam_dum, argb2_here, ampb1I_here, argb1I_here, ampb1_here, argb1_here, ampb2_here = np.loadtxt(filename_meas,unpack=True)
        argb2 += argb2_here.reshape((len(F0_arr),len(ring_diam_arr)))
        ampb1I += ampb1I_here.reshape((len(F0_arr),len(ring_diam_arr)))
        argb1I += argb1I_here.reshape((len(F0_arr),len(ring_diam_arr)))
        ampb1 += ampb1_here.reshape((len(F0_arr),len(ring_diam_arr)))
        argb1 += argb1_here.reshape((len(F0_arr),len(ring_diam_arr)))
        ampb2 += ampb2_here.reshape((len(F0_arr),len(ring_diam_arr)))

        count += 1.0

# average over all instantiations
Delta_d /= count
Delta_alpha /= count
Delta_argb2 /= count
Delta_ampb1I /= count
Delta_argb1I /= count
Delta_ampb1 /= count
Delta_argb1 /= count
Delta_ampb2 /= count
argb2 /= count
ampb1I /= count
argb1I /= count
ampb1 /= count
argb1 /= count
ampb2 /= count

##################################################
# make source counts plot

# create a grid in (theta,sigma)
logtheta_min = -1.0
logtheta_max = 2.0
logsigma_min = -3.0
logsigma_max = 1.0
theta1d = 10.0**np.linspace(logtheta_min,logtheta_max,100)
sigma1d = 10.0**np.linspace(logsigma_min,logsigma_max,100)
theta, sigma = np.meshgrid(theta1d,sigma1d,indexing='ij')

# compute the source counts on the grid
counts = source_counts(theta,sigma,counts_type='full',BHMF_type='KH')
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

# add contours
levels_neg = np.arange(logcounts_min,0.0)
levels_pos = np.arange(1.0,10.0)
cs1 = ax.contour(logtheta,logsigma,logcounts,levels=levels_neg,colors='k',linewidths=1,linestyles='--')
cs2 = ax.contour(logtheta,logsigma,logcounts,levels=[0.0],colors='k',linewidths=2,linestyles='-')
cs3 = ax.contour(logtheta,logsigma,logcounts,levels=levels_pos,colors='k',linewidths=1,linestyles='-')

##################################################
# mass

# get current measurements mappable
mappable = np.copy(Delta_d)
for i in range(len(F0_arr)):
    mappable[i,:] /= ring_diam_arr

# plot and retrieve the 20% contour
cs = ax.contour(np.log10(ring_diam_arr),np.log10(F0_arr*1000.0),mappable,levels=[0.2],colors='white',alpha=0.0)
paths = cs.get_paths()

max_counts = 0.0
for ipath, path in enumerate(paths):
    points = path.vertices
    x = points[:,0]
    y = points[:,1]
    ax.plot(x,y,'w-')
    ax.plot(x,y,'r--')

    # determine maximum source counts along the contour
    test_theta = 10.0**x
    test_sigma = (10.0**y)/1000.0
    contour_counts = source_counts(test_theta,test_sigma,counts_type='full',BHMF_type='KH')
    if np.max(contour_counts) > max_counts:
        max_counts = np.max(contour_counts)

    # export contour
    with open(outdir+'/mass_contour_'+str(ipath).zfill(3)+'.txt','w') as outtable:
        for ix in range(len(x)):
            strhere = ''
            strhere += str(x[ix]).ljust(24)
            strhere += str(y[ix]) + '\n'
            outtable.write(strhere)

count_here = max_counts
print('Maximum source counts for mass measurement: '+str(count_here))

##################################################
# spin

# get current measurements mappables
mappable_here1 = np.copy(Delta_argb2)
mappable_here2 = np.copy(Delta_ampb2) / (0.3 / np.sqrt(2.0))
mappable_here3 = np.copy(Delta_ampb1I) / (0.5 / np.sqrt(2.0))
mappable_here4 = np.copy(Delta_ampb1) / (0.05 / np.sqrt(2.0))
mappable_here5 = np.copy(Delta_argb1)

cs1 = ax.contour(np.log10(ring_diam_arr),np.log10(F0_arr*1000.0),mappable_here1,levels=[0.2],colors='white',linewidths=[0])
paths = cs1.get_paths()
maxlength = 0
for path in paths:
    points = path.vertices
    if len(points) >= maxlength:
        maxlength = len(points)
        x1 = points[:,0]
        y1 = points[:,1]

cs2 = ax.contour(np.log10(ring_diam_arr),np.log10(F0_arr*1000.0),mappable_here2,levels=[0.2],colors='white',linewidths=[0])
paths = cs2.get_paths()
maxlength = 0
for path in paths:
    points = path.vertices
    if len(points) >= maxlength:
        maxlength = len(points)
        x2 = points[:,0]
        y2 = points[:,1]

cs3 = ax.contour(np.log10(ring_diam_arr),np.log10(F0_arr*1000.0),mappable_here3,levels=[0.2],colors='white',linewidths=[0])
paths = cs3.get_paths()
maxlength = 0
for path in paths:
    points = path.vertices
    if len(points) >= maxlength:
        maxlength = len(points)
        x3 = points[:,0]
        y3 = points[:,1]

cs4 = ax.contour(np.log10(ring_diam_arr),np.log10(F0_arr*1000.0),mappable_here4,levels=[0.2],colors='white',linewidths=[0])
paths = cs4.get_paths()
maxlength = 0
for path in paths:
    points = path.vertices
    if len(points) >= maxlength:
        maxlength = len(points)
        x4 = points[:,0]
        y4 = points[:,1]

cs5 = ax.contour(np.log10(ring_diam_arr),np.log10(F0_arr*1000.0),mappable_here5,levels=[0.2],colors='white',linewidths=[0])
paths = cs5.get_paths()
maxlength = 0
for path in paths:
    points = path.vertices
    if len(points) >= maxlength:
        maxlength = len(points)
        x5 = points[:,0]
        y5 = points[:,1]

if ((len(x1) >= 2) & (len(x2) >= 2) & (len(x3) >= 2) & (len(x4) >= 2) & (len(x5) >= 2)):
    
    f1 = interpolate.interp1d(x1, y1, fill_value='extrapolate')
    f2 = interpolate.interp1d(x2, y2, fill_value='extrapolate')
    f3 = interpolate.interp1d(x3, y3, fill_value='extrapolate')
    f4 = interpolate.interp1d(x4, y4, fill_value='extrapolate')
    f5 = interpolate.interp1d(x5, y5, fill_value='extrapolate')

    dumx = np.linspace(-1.0,2.0,100)
    dumy1 = f1(dumx)
    dumy2 = f2(dumx)
    dumy3 = f3(dumx)
    dumy4 = f4(dumx)
    dumy5 = f5(dumx)
    dumy = np.zeros_like(dumx)
    for i in range(len(dumx)):
        yhere = np.max([dumy1[i],dumy2[i],dumy3[i],dumy4[i],dumy5[i]])
        dumy[i] = yhere
    ax.plot(dumx,dumy,linestyle='-',color='white')
    ax.plot(dumx,dumy,'g--')

    # export contour
    with open(outdir+'/spin_contour.txt','w') as outtable:
        for ix in range(len(dumx)):
            strhere = ''
            strhere += str(dumx[ix]).ljust(24)
            strhere += str(dumy[ix]) + '\n'
            outtable.write(strhere)

    # determine maximum source counts (per steradian) along the contour
    test_theta = 10.0**dumx
    test_sigma = (10.0**dumy)/1000.0
    contour_counts = source_counts(test_theta,test_sigma,counts_type='full',BHMF_type='KH')
    count_here = np.max(contour_counts)
    print('Maximum source counts for spin measurement: '+str(count_here))

##################################################
# shadow

# get current measurements mappable
mappable = np.copy(Delta_alpha)
for i in range(len(F0_arr)):
    mappable[i,:] /= ((2.0/3.0)*ring_diam_arr)

# plot and retrieve the 20% contour
cs = ax.contour(np.log10(ring_diam_arr),np.log10(F0_arr*1000.0),mappable,levels=[0.2],colors='white',alpha=0.0)
paths = cs.get_paths()

max_counts = 0.0
for ipath, path in enumerate(paths):
    points = path.vertices
    x = points[:,0]
    y = points[:,1]
    ind = np.argsort(x)
    ax.plot(x[ind],y[ind],'w-')
    ax.plot(x[ind],y[ind],'b--')

    # determine maximum source counts along the contour
    test_theta = 10.0**x
    test_sigma = (10.0**y)/1000.0
    contour_counts = source_counts(test_theta,test_sigma,counts_type='full',BHMF_type='KH')
    if np.max(contour_counts) > max_counts:
        max_counts = np.max(contour_counts)

    # export contour
    with open(outdir+'/shadow_contour_'+str(ipath).zfill(3)+'.txt','w') as outtable:
        for ix in range(len(x)):
            strhere = ''
            strhere += str(x[ix]).ljust(24)
            strhere += str(y[ix]) + '\n'
            outtable.write(strhere)

count_here = max_counts
print('Maximum source counts for shadow measurement: '+str(count_here))

##################################################
# remaining plot items

# label contours
xtext = -0.95
ytext = 2.42
ax.text(xtext,ytext,r'$N=1$',fontsize=10,ha='left',va='bottom')
ytext = 1.72
ax.text(xtext,ytext,r'$N=10$',fontsize=10,ha='left',va='bottom')
ytext = 1.03
ax.text(xtext,ytext,r'$N=10^2$',fontsize=10,ha='left',va='bottom')
ytext = 0.32
ax.text(xtext,ytext,r'$N=10^3$',fontsize=10,ha='left',va='bottom')

ax.set_xlim(-1,2)
ax.set_ylim(0,4)

ax.set_xlabel(r'$\log(\theta/\mu\rm{as})$')
ax.set_ylabel(r'$\log(S_{\nu}/\rm{mJy})$')

plt.savefig(outdir+'/source_counts.png',dpi=300,bbox_inches='tight')
plt.close()
