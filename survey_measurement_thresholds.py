##################################################
# imports etc.

import numpy as np
import ngehtsim.obs.obs_generator as og
import ngEHTforecast.fisher as fp
import ehtim as eh
import ehtim.parloop as ploop
import os

# set random number seed
np.random.seed(12345)

##################################################
# source counts survey settings

# flux density grid, in Jy
F0_arr = 10.0**np.linspace(-3.0,1.0,40)

# ring diameter grid, in microarcseconds
ring_diam_arr = 10.0**np.linspace(-1.0,2.0,30)

# number of weather instantiations to use per grid point
Nweather = 20

##################################################
# array and observation settings

# coordinates of source
RA = 12.513729
DEC = 12.391123

# telescopes in the array
sites = ['ALMA',
         'IRAM',
         'LMT',
         'NOEMA',
         'SMA',
         'SMT']

# default settings
settings = {'RA': RA,
            'DEC': DEC,
            'sites': sites,
            'frequency': 230.0,
            'month': 'Apr',
            'day': 11,
            'year': 2017,
            't_start': 0.0,
            'dt': 24.0,
            't_int': 600.0,
            't_rest': 1200.0,
            'fringe_finder': ['fringegroups', [5.0, 10.0]],
            'ttype': 'fast',
            'weather': 'random'}

##################################################
# miscellaneous functions

# extract marginalized uncertainties using FF
def marginalized_uncertainties(ff,obs,p,i,verbosity=0,**kwargs) :
    C = ff.fisher_covar(obs,p,verbosity=verbosity,**kwargs)
    M = np.zeros((ff.size-1,ff.size-1))
    v = np.zeros(ff.size-1)
    ilist = np.arange(ff.size)
    ini = ilist[ilist!=i]
    for j2,j in enumerate(ini) :
        for k2,k in enumerate(ini) :
            M[j2,k2] = C[j][k]
        v[j2] = C[i][j]
    N = C[i][i]
    Minv = fp.fisher_forecast._invert_matrix(M)
    mN = N - np.matmul(v,np.matmul(Minv,v))
    Sig_marg = np.sqrt(2.0/mN)
    return Sig_marg

##################################################
# loop through weather instantiations

for ii in range(Nweather):

    # make directory to hold files for this weather instantiation
    indir = './survey/weather_inst_'+str(ii).zfill(4)
    os.makedirs(indir,exist_ok=True)

    # initialize arrays to contain parameter values
    ampb1I = np.zeros((len(F0_arr),len(ring_diam_arr)))
    argb1I = np.zeros((len(F0_arr),len(ring_diam_arr)))
    ampb1 = np.zeros((len(F0_arr),len(ring_diam_arr)))
    argb1 = np.zeros((len(F0_arr),len(ring_diam_arr)))
    ampb2 = np.zeros((len(F0_arr),len(ring_diam_arr)))
    argb2 = np.zeros((len(F0_arr),len(ring_diam_arr)))
    Delta_d = np.zeros((len(F0_arr),len(ring_diam_arr)))
    Delta_alpha = np.zeros((len(F0_arr),len(ring_diam_arr)))
    Delta_ampb1I = np.zeros((len(F0_arr),len(ring_diam_arr)))
    Delta_argb1I = np.zeros((len(F0_arr),len(ring_diam_arr)))
    Delta_ampb1 = np.zeros((len(F0_arr),len(ring_diam_arr)))
    Delta_argb1 = np.zeros((len(F0_arr),len(ring_diam_arr)))
    Delta_ampb2 = np.zeros((len(F0_arr),len(ring_diam_arr)))
    Delta_argb2 = np.zeros((len(F0_arr),len(ring_diam_arr)))

    # loop through grid
    for iF0, F0 in enumerate(F0_arr):
        for iring, ring_diam in enumerate(ring_diam_arr):

            print('Running progress:  Weather ('+str(ii+1)+'/'+str(Nweather)+'), Flux values ('+str(iF0+1)+'/'+str(len(F0_arr))+'), Diameter values ('+str(iring+1)+'/'+str(len(ring_diam_arr))+')')

            ##################################################
            # set up m-ring model

            d = ring_diam*eh.RADPERUAS
            alpha = d / 3.0
            stretch = 1.0
            stretch_PA = 0.0
            beta_list = np.random.uniform(-0.5,0.5,size=1) + (1.0j)*np.random.uniform(-0.5,0.5,size=1)

            # randomly select ring parameters within GRMHD-motivated range
            bpol0 = np.random.uniform(-0.1,0.1) + (1.0j)*np.random.uniform(-0.1,0.1)
            bpoln1 = np.random.uniform(-0.05,0.05) + (1.0j)*np.random.uniform(-0.05,0.05)
            bpol1 = np.random.uniform(-0.05,0.05) + (1.0j)*np.random.uniform(-0.05,0.05)
            bpoln2 = np.random.uniform(-0.1,0.1) + (1.0j)*np.random.uniform(-0.1,0.1)
            bpol2 = np.random.uniform(-0.3,0.3) + (1.0j)*np.random.uniform(-0.3,0.3)
            beta_list_pol = [bpoln2,bpoln1,bpol0,bpol1,bpol2]

            mod = eh.model.Model()
            mod = mod.add_stretched_thick_mring(F0=F0,
                                                d=d,
                                                alpha=alpha,
                                                x0=0.0,
                                                y0=0.0,
                                                beta_list=beta_list,
                                                beta_list_pol=beta_list_pol,
                                                stretch=stretch,
                                                stretch_PA=stretch_PA)

            mod.ra = RA
            mod.dec = DEC

            # record ring values
            argb2[iF0,iring] = np.angle(bpol2)
            ampb1I[iF0,iring] = np.abs(beta_list[0])
            argb1I[iF0,iring] = np.angle(beta_list[0])
            ampb1[iF0,iring] = np.abs(bpol1)
            argb1[iF0,iring] = np.angle(bpol1)
            ampb2[iF0,iring] = np.abs(bpol2)

            ##################################################
            # generate the observation

            # use new random seed
            settings.update({'random_seed': (101*ii) + 102})

            # initialize the obs_generator
            obsgen = og.obs_generator(settings)

            # generate obs
            with ploop.HiddenPrints():
                obs = obsgen.make_obs(mod,
                                      addnoise=True,
                                      addgains=False,
                                      flagday=True,
                                      flagwind=True)

            ##################################################
            # estimate measurement uncertainties

            if len(obs.data) > 0:

                # initialize FF object
                m = 1
                mp = 2
                mc = 0
                stokes = 'full'
                ffpg = fp.ff_models.FF_thick_mring(m,mp,mc,stokes=stokes)
                ffg = fp.ff_complex_gains.FF_complex_gains(ffpg)

                # parameter list
                params = mod.params[0]
                p = list()
                p.append(params['F0'])
                p.append(params['d']/eh.RADPERUAS)
                p.append(params['alpha']/eh.RADPERUAS)
                p.append(params['beta_list'])
                p.append(params['beta_list_pol'])

                pcheck = list()
                for i in range(len(p)):
                    if (type(p[i]) != np.ndarray):
                        pcheck.append(p[i])
                    else:
                        for item in p[i]:
                            pcheck.append(np.abs(item))
                            pcheck.append(np.angle(item))

                # set priors
                ffg.add_gaussian_prior(0,1e-10)
                ffg.add_gaussian_prior(1,d/eh.RADPERUAS)
                ffg.add_gaussian_prior(2,alpha/eh.RADPERUAS)
                ffg.set_gain_epochs(scans=True)
                ffg.set_gain_amplitude_prior(0.1)
                ffg.set_gain_phase_prior(3.0)
                ffg.set_gain_ratio_amplitude_prior(0.01)
                ffg.set_gain_ratio_phase_prior(0.01)

                # compute marginalized uncertainties
                Delta_d[iF0,iring] = marginalized_uncertainties(ffg,obs,pcheck,1)
                Delta_alpha[iF0,iring] = marginalized_uncertainties(ffg,obs,pcheck,2)
                Delta_ampb1I[iF0,iring] = marginalized_uncertainties(ffg,obs,pcheck,3)
                Delta_argb1I[iF0,iring] = marginalized_uncertainties(ffg,obs,pcheck,4)
                Delta_ampb1[iF0,iring] = marginalized_uncertainties(ffg,obs,pcheck,11)
                Delta_argb1[iF0,iring] = marginalized_uncertainties(ffg,obs,pcheck,12)
                Delta_ampb2[iF0,iring] = marginalized_uncertainties(ffg,obs,pcheck,13)
                Delta_argb2[iF0,iring] = marginalized_uncertainties(ffg,obs,pcheck,14)

            else:

                Delta_d[iF0,iring] = np.nan
                Delta_alpha[iF0,iring] = np.nan
                Delta_ampb1I[iF0,iring] = np.nan
                Delta_argb1I[iF0,iring] = np.nan
                Delta_ampb1[iF0,iring] = np.nan
                Delta_argb1[iF0,iring] = np.nan
                Delta_ampb2[iF0,iring] = np.nan
                Delta_argb2[iF0,iring] = np.nan

    ##################################################
    # write tables

    # table containing the input values
    with open(indir+'/ring_paramtable.txt','w') as outtable:
        for iF0, F0 in enumerate(F0_arr):
            for iring, ring_diam in enumerate(ring_diam_arr):
                outstr = ''
                outstr += str(np.round(F0,10)).ljust(24)
                outstr += str(np.round(ring_diam,10)).ljust(24)
                outstr += str(np.round(argb2[iF0,iring],10)).ljust(24)
                outstr += str(np.round(ampb1I[iF0,iring],10)).ljust(24)
                outstr += str(np.round(argb1I[iF0,iring],10)).ljust(24)
                outstr += str(np.round(ampb1[iF0,iring],10)).ljust(24)
                outstr += str(np.round(argb1[iF0,iring],10)).ljust(24)
                outstr += str(np.round(ampb2[iF0,iring],10)) + '\n'
                outtable.write(outstr)

    # table containing the estimated measurement uncertainties
    with open(indir+'/paramtable.txt','w') as outtable:
        for iF0, F0 in enumerate(F0_arr):
            for iring, ring_diam in enumerate(ring_diam_arr):
                outstr = ''
                outstr += str(np.round(F0,10)).ljust(24)
                outstr += str(np.round(ring_diam,10)).ljust(24)
                outstr += str(np.round(Delta_d[iF0,iring],10)).ljust(24)
                outstr += str(np.round(Delta_alpha[iF0,iring],10)).ljust(24)
                outstr += str(np.round(Delta_argb2[iF0,iring],10)).ljust(24)
                outstr += str(np.round(Delta_ampb1I[iF0,iring],10)).ljust(24)
                outstr += str(np.round(Delta_argb1I[iF0,iring],10)).ljust(24)
                outstr += str(np.round(Delta_ampb1[iF0,iring],10)).ljust(24)
                outstr += str(np.round(Delta_argb1[iF0,iring],10)).ljust(24)
                outstr += str(np.round(Delta_ampb2[iF0,iring],10)) + '\n'
                outtable.write(outstr)

##################################################
