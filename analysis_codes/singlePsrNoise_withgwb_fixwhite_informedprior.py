#!/usr/bin/env python
# coding: utf-8

import os
import sys
import numpy as np
from enterprise import constants as const
from enterprise.pulsar import Pulsar
from enterprise.signals import signal_base
from enterprise.signals import white_signals
from enterprise.signals import gp_signals
from enterprise.signals import parameter
from enterprise.signals import utils
from enterprise.signals import selections
from enterprise.signals import gp_priors
from enterprise.signals import deterministic_signals
from enterprise.signals import gp_bases
from enterprise_extensions.chromatic.solar_wind import solar_wind, createfourierdesignmatrix_solar_dm
# from enterprise_extensions import blocks
from enterprise_extensions import hypermodel
# from PTMCMCSampler.PTMCMCSampler import PTSampler as ptmcmc

from ppta_dr3_models import *
    
# Read in pulsar name, data directory, and chain number
psrname = sys.argv[1]
chainnum = sys.argv[2]
dir = sys.argv[3]
datadir = os.path.abspath("./data/" + str(dir) + '/')
reddm_only = False

chain_label = ''
if len(sys.argv) > 4:
    chain_label = '_'+str(sys.argv[4])
    
# load in pulsar data
pfile = os.path.join(datadir, psrname + ".par")
tfile = os.path.join(datadir, psrname + ".tim")
psr = Pulsar(pfile, tfile, ephem='DE440')


if 'all' in dir:
    dir = 'all'
if 'dr2' in dir:
    dir = 'dr2'
if 'uwl' in dir:
    dir = 'uwl'

dirnoise = os.path.abspath("./data/")
noisefiles = sorted(glob.glob(dirnoise + '/' + str(dir) + '/noiseFiles_maxlike/{}*.json'.format(psrname)))

noisedict = {}
### DO UWL NOISE FIRST, TO RE-NAME ECORR_ALL
for nf in noisefiles:
    if 'dr2' in nf:
        continue
    with open(nf, 'r') as f:
        noisedict.update(json.load(f))
    if 'n_earth' in noisedict.keys():
            noisedict['{}_n_earth'.format(os.path.basename(nf).split('_')[0])] = noisedict['n_earth']
            del(noisedict['n_earth'])
keys = np.array([str(key) for key in noisedict.keys()]).squeeze()
for key in keys:
    if "basis_ecorr_all" in key and not 'uwl' in key and not 'dr2' in key:
       new_key = key.replace("basis_ecorr_all", "basis_ecorr_all_uwl")
       noisedict[new_key] = noisedict[key]


               
### NOW DO DR2 NOISE, TO RE-NAME ECORR_ALL
for nf in noisefiles:
    if 'uwl' in nf:
        continue
    with open(nf, 'r') as f:
        noisedict.update(json.load(f))
    if 'n_earth' in noisedict.keys():
        del(noisedict['n_earth']) 
keys = np.array([str(key) for key in noisedict.keys()]).squeeze()
for key in keys:
    if "basis_ecorr_all" in key and not 'uwl' in key and not 'dr2' in key:
       new_key = key.replace("basis_ecorr_all", "basis_ecorr_all_dr2")
       noisedict[new_key] = noisedict[key]

#NOW FUSS OVER GROUP TERM ECORR NAMING:
for key in keys:
    if 'ecorr' not in key:
        continue
    for grouppsr, grouplist_ in psr_groupecorr_dict_dict[dir].items():
        
        if grouppsr not in key:
            
            continue
        for group_ in grouplist_:
            if group_ not in key:
                continue
            group_done = False
            if f'{grouppsr}_basis_ecorr_{group_}' in key:
                new_key = key.replace(f'basis_ecorr_{group_}', f'basis_ecorr_group_{group_}')
                noisedict[new_key] = noisedict[key]
                print(f'renamed {key} to {new_key}')
                del noisedict[key]
                group_done = True
            #check if it's already been converted to new name
            elif f'{grouppsr}_basis_ecorr_group_{group_}' in key:
                group_done = True
            #check if group isn't added in - e.g. it is not in the `all` noise files
            if not group_done:
                print('COULD NOT FIND {} TERM FOR {} (KEY IS {})'.format(group_, grouppsr, key))
                raise ValueError
                print('looking in uwl noise file for {} ecorr term'.format(group_))
                uwl_gnf = sorted(glob.glob(dirnoise+f'/uwl/noiseFiles_maxlike/{grouppsr}*.json'))[0]
                with open(uwl_gnf, 'r') as uwl_nf:
                    #uwl noisedict
                    uwl_nd = json.load(uwl_nf)
                    keys__ = np.array([str(key) for key in uwl_nd.keys()]).squeeze()
                    for k in keys__:
                        if f'{grouppsr}_basis_ecorr_{group_}' in k or f'{grouppsr}_basis_ecorr_group_{group_}' in k:
                            noisedict[k.replace(f'basis_ecorr_{group_}', f'basis_ecorr_group_{group_}')] = uwl_nd[k]
    
print(noisedict)
    
if not dir=='ppta15':
    # Set up group noise
    psr_groupnoiselist_dict = {psr: None for psr in psrnames}
    psr_groupnoiselist_dict.update(psr_groupnoise_dict_dict[dir])
    by_group_dict = {key: selections.Selection(sel_by_group_factory(item)._sel_by_group) for key, item in psr_groupnoiselist_dict.items()}

    # Set up group-selected ecorrs
    psr_groupecorrlist_dict = {psr: None for psr in psrnames}
    psr_groupecorrlist_dict.update(psr_groupecorr_dict_dict[dir])
    by_group_ecorr_dict = {key: selections.Selection(sel_by_group_factory(item)._sel_by_group) for key, item in psr_groupecorrlist_dict.items()}


if dir == 'all':
    if len(glob.glob(dirnoise + '/all/noiseFiles/3sig/*p0015s.json')) > 0:
        noisedict_p0015 = {}
        noisefiles_p0015 = sorted(glob.glob(dirnoise+'/all/noiseFiles/3sig/*p0015s.json'))
        for nf in noisefiles_p0015:
            with open(nf, 'r') as f:
                noisedict_p0015.update(json.load(f))
                psrname__ = nf.split('_')[0].split('/')[-1]
                print('noise dict pulsar name {}'.format(psrname__))
                if "n_earth" in noisedict_p0015.keys():
                    noisedict_p0015[psrname__ + '_n_earth'] = noisedict_p0015['n_earth']
                del noisedict_p0015['n_earth']
        noisedict_p9985 = {}
        noisefiles_p9985 = sorted(glob.glob(dirnoise+'/all/noiseFiles/3sig/*p9985s.json'))
        for nf in noisefiles_p9985:
            with open(nf, 'r') as f:
                noisedict_p9985.update(json.load(f))
                psrname__ = nf.split('_')[0].split('/')[-1]
                if 'n_earth' in noisedict_p9985.keys():
                    noisedict_p9985[psrname__ + '_n_earth'] = noisedict_p9985['n_earth']
                del noisedict_p9985['n_earth']

    else:
        # take 3-sigma bounds from dr2 only
        print("Couldn't find 3-sigma noisefiles!")
        dir_3sig = os.path.abspath("./data/")
        noisefiles_p9985_dr2  = sorted(glob.glob(dir_3sig + '/dr2/noiseFiles/3sig/*p9985s.json'))
        noisefiles_p0015_dr2  = sorted(glob.glob(dir_3sig + '/dr2/noiseFiles/3sig/*p0015s.json'))
        noisefiles_p9985_uwl  = sorted(glob.glob(dir_3sig + '/uwl/noiseFiles/3sig/*p9985s.json'))
        noisefiles_p0015_uwl  = sorted(glob.glob(dir_3sig + '/uwl/noiseFiles/3sig/*p0015s.json'))    
        #print(noisefiles)
        
        noisedict_p9985_dr2 = {}
        for nf in noisefiles_p9985_dr2:
            with open(nf, 'r') as f:
                #print(nf)
                noisedict_p9985_dr2.update(json.load(f))
        noisedict_p0015_dr2 = {}
        for nf in noisefiles_p0015_dr2:
            with open(nf, 'r') as f:

                #print(nf)
                noisedict_p0015_dr2.update(json.load(f))
        noisedict_p9985_uwl = {}
        for nf in noisefiles_p9985:
            with open(nf, 'r') as f:
                #print(nf)
                noisedict_p9985.update(json.load(f))
        noisedict_p0015_uwl = {}
        for nf in noisefiles_p0015:
            with open(nf, 'r') as f:
                #print(nf)
                noisedict_p0015.update(json.load(f))

        #go through each noise term in DR2 and add it into the overall p0015 dict
        #if the noise term is also in UWL, add in the minimum of the DR2 or UWL values
        noisedict_p0015 = {}
        for key, item in noisedict_p0015_dr2.items():
            if key in noisedict_p0015_uwl:
                val_ = np.min(noisedict_p0015_uwl[key], item)
            else:
                val_ = item
            noisedict_p0015[key] = val_
        #add in any leftover terms that weren't in DR2 dict:
        for key, item in noisedict_p0015_uwl.items():
            #if they're not in the overall dict then they weren't in DR2
            if not key in noisedict_p0015:
                noisedict_p0015[key] = item
                
        #go through each noise term in DR2 and add it into the overall 3sig_max dict
        #if the noise term is also in UWL, add in the max of the DR2 or UWL values
        noisedict_p9985 = {}
        for key, item in noisedict_p9985_dr2.items():
            if key in noisedict_p9985_uwl:
                val_ = np.max(noisedict_p9985_uwl[key], item)
            else:
                val_ = item
            noisedict_p9985[key] = val_
        #add in any leftover terms that weren't in DR2 dict:
        for key, item in noisedict_p9985_uwl.items():
            #if they're not in the overall dict then they weren't in DR2
            if not key in noisedict_p9985:
                noisedict_p9985[key] = item
                
else:
    dir_3sig = os.path.abspath("./data/")
    noisefiles_p9985  = sorted(glob.glob(dir_3sig + '/' + str(dir) + '/noiseFiles/3sig/*p9985s.json'))
    noisefiles_p0015  = sorted(glob.glob(dir_3sig + '/' + str(dir) + '/noiseFiles/3sig/*p0015s.json'))
    noisedict_p0015 = {}
    for nf in noisefiles_p0015:
        with open(nf, 'r') as f:
            #print(nf)
            noisedict_p0015.update(json.load(f))
        if 'n_earth' in noisedict_p0015.keys():
            noisedict_p0015['{}_n_earth'.format(os.path.basename(nf).split('_')[0])] = noisedict_p0015['n_earth']
            del(noisedict_p0015['n_earth'])
    noisedict_p9985 = {}
    for nf in noisefiles_p9985:
        with open(nf, 'r') as f:
            #print(nf)
            noisedict_p9985.update(json.load(f))
        if 'n_earth' in noisedict_p9985.keys():
            noisedict_p9985['{}_n_earth'.format(os.path.basename(nf).split('_')[0])] = noisedict_p9985['n_earth']
            del(noisedict_p9985['n_earth'])

"""
Choose whether to marginalise or sample the timing model
"""
tm = gp_signals.MarginalizingTimingModel(use_svd=True)

"""
Define white noise model
"""
# EFAC "MeasurementNoise" can add equad, but only t2equad - we want the tnequad
efac_prior = parameter.Uniform(0.01, 10.0)
wn = white_signals.MeasurementNoise(
    efac=efac_prior,
    selection=by_backend)

# EQUAD - TempoNest definition: sigma = sqrt((efac*sigma_0)**2 + (tnequad)**2)
log10_equad_prior = parameter.Constant()#Uniform(-10, -5)
wn += white_signals.TNEquadNoise(
    log10_tnequad=log10_equad_prior,
    selection=by_backend)

# ECORR - we will swap to "white_signals.EcorrKernelNoise" later
if not dir=='ppta15':
    log10_ecorr_prior = parameter.Constant()#Uniform(-10, -5)
    wn += gp_signals.EcorrBasisModel(
        log10_ecorr=log10_ecorr_prior,
        selection=ecorr_selection)
    if psr.name in psr_groupecorr_dict_dict[dir].keys():
        wn += gp_signals.EcorrBasisModel(log10_ecorr = log10_ecorr_prior,
                                         selection = by_group_ecorr_dict[psr.name],
                                         name='basis_ecorr_group')
    if dir == 'uwl' or dir == 'dr2':
        wn += gp_signals.EcorrBasisModel(
            log10_ecorr=log10_ecorr_prior,
            selection=no_selection, name='basis_ecorr_all')
    elif dir == 'all':
        wn += gp_signals.EcorrBasisModel(
            log10_ecorr=log10_ecorr_prior,
            selection=global_ecorr_selection, name='basis_ecorr_all')

"""
Define red noise model
"""
log10_A_prior = parameter.Uniform(-20, -11)
gamma_prior = parameter.Uniform(0, 7)
# # powerlaw
rn_model = gp_priors.powerlaw(log10_A=log10_A_prior,
                              gamma=gamma_prior)

priors = {}

rn_model, rn_lgA_prior, rn_gam_prior, priors = get_informed_rednoise_priors(psr, 'red_noise', noisedict_p0015, noisedict_p9985, priors)
# broken powerlaw
#lfb_prior = parameter.Uniform(-10, -7)  # log10 transition frequency
#kappa_prior = 0.1  # smoothness of transition (Default = 0.1)
#delta_prior = 0  # slope for frequencies > f_break
#rn_model = gp_priors.broken_powerlaw(log10_A=log10_A_prior,
#                                     gamma=gamma_prior,
#                                     delta=delta_prior,
#                                     log10_fb=lfb_prior,
#                                     kappa=kappa_prior)
Tspan = psr.toas.max() - psr.toas.min()  # seconds
max_cadence = 240  # days
red_components = int(Tspan / (max_cadence*86400))
print("Using {} red noise components".format(red_components))
rn = gp_signals.FourierBasisGP(rn_model, components=red_components,
                               selection=no_selection, name='red_noise')

gwb_model = gp_priors.powerlaw(log10_A=log10_A_prior,
                              gamma=4.33)
print("Using {} gwb components".format(red_components))
gwb = gp_signals.FourierBasisGP(gwb_model, components=red_components,
                               selection=no_selection, name='gwb')


hf_model, hf_lgA_prior, hf_gam_prior, priors = get_informed_rednoise_priors(psr, 'hf_noise', noisedict_p0015, noisedict_p9985, priors)
    
max_cadence = 30  # days
hf_components = int(Tspan / (max_cadence*86400))
print("Using {} hf achromatic noise components".format(hf_components))
hf = gp_signals.FourierBasisGP(hf_model, components=hf_components,
                            selection=no_selection, name='hf_noise')


band_model, band_lgA_prior, band_gam_prior, priors = get_informed_rednoise_priors(psr, 'band_noise_low', noisedict_p0015, noisedict_p9985, priors)

bandmid_model, bandmid_lgA_prior, bandmid_gam_prior, priors = get_informed_rednoise_priors(psr, 'band_noise_mid', noisedict_p0015, noisedict_p9985, priors)
bandhigh_model, bandhigh_lgA_prior, bandhigh_gam_prior, priors = get_informed_rednoise_priors(psr, 'band_noise_high', noisedict_p0015, noisedict_p9985, priors)

max_cadence = 60  # days
band_components = int(Tspan / (max_cadence*86400))
bn = gp_signals.FourierBasisGP(band_model, components=band_components,
                               selection=low_freq, name='band_noise_low')

bn_mid = gp_signals.FourierBasisGP(bandmid_model, components=band_components,
                                   selection=mid_freq, name='band_noise_mid')

bn_high = gp_signals.FourierBasisGP(bandhigh_model, components=band_components,
                                    selection=high_freq, name='band_noise_high')

"""
Define system noise model
"""
if not dir == 'ppta15':
    max_cadence = 30  # days
    gn_components = int(Tspan / (max_cadence * 86400.0))
    gn_lgA_prior_min = np.inf
    gn_lgA_prior_max = -np.inf
    gn_gamma_prior_min = np.inf
    gn_gamma_prior_max = -np.inf        
    if psr_groupnoiselist_dict[psr.name] is None:
        gn = None
    else:
        for group in psr_groupnoiselist_dict[psr.name]:

            gn_model__, gn_lgA_prior, gn_gam_prior, gn_lgA_prior_min_, gn_lgA_prior_max_, gn_gamma_prior_min_, gn_gamma_prior_max_, priors = get_informed_rednoise_priors(psr, 'group_noise_' + group, noisedict_p0015, noisedict_p9985, priors, return_priorvals = True)
            if gn_lgA_prior_min_ < gn_lgA_prior_min:
                gn_lgA_prior_min = gn_lgA_prior_min_
            if gn_lgA_prior_max_ > gn_lgA_prior_max:
                gn_lgA_prior_max = gn_lgA_prior_max_
            if gn_gamma_prior_min_ < gn_gamma_prior_min:
                gn_gamma_prior_min = gn_gamma_prior_min_
            if gn_gamma_prior_max_ > gn_gamma_prior_max:
                gn_gamma_prior_max = gn_gamma_prior_max_


        gn_lgA_prior_mins = [-18, gn_lgA_prior_min]
        gn_lgA_prior_maxs = [-11, gn_lgA_prior_max]
        gn_gamma_prior_mins = [0, gn_gamma_prior_min]
        gn_gamma_prior_maxs = [7, gn_gamma_prior_max]
        #for group noise, we look over all priors and set a global one based on the overall max and mins of the individual priors
        gn_lgA_prior = parameter.Uniform(np.max(gn_lgA_prior_mins), np.min(gn_lgA_prior_maxs))
        gn_gamma_prior = parameter.Uniform(np.max(gn_gamma_prior_mins), np.min(gn_gamma_prior_maxs))
        gn_model = gp_priors.powerlaw(log10_A=gn_lgA_prior,
                                      gamma=gn_gamma_prior)

        gn = FourierBasisGP_ppta(gn_model, fmax=1/(30*86400),
                                 selection = by_group_dict[psr.name],
                                 name='group_noise')


"""
Define DM noise model
"""

dm_model, dm_lgA_prior, dm_gam_prior, priors = get_informed_rednoise_priors(psr, 'dm_gp', noisedict_p0015, noisedict_p9985, priors)

Tspan = psr.toas.max() - psr.toas.min()  # seconds
max_cadence = 60  # days
dm_components = int(Tspan / (max_cadence*86400))
print("Using {} DM components".format(dm_components))
dm_basis = gp_bases.createfourierdesignmatrix_dm(nmodes=dm_components)
dm = gp_signals.BasisGP(dm_model, dm_basis, name='dm_gp')

"""
Define chromatic noise model
"""

chrom_model, chrom_lgA_prior, chrom_gam_prior, priors = get_informed_rednoise_priors(psr, 'chrom_gp', noisedict_p0015, noisedict_p9985, priors)

idx = 4  # Define freq^-idx scaling
max_cadence = 240  # days
chrom_components = int(Tspan / (max_cadence*86400))
print("Using {} Chrom components".format(chrom_components))
chrom_basis = gp_bases.createfourierdesignmatrix_chromatic(nmodes=chrom_components,
                                                           idx=idx)
chrom = gp_signals.BasisGP(chrom_model, chrom_basis, name='chrom_gp')

"""
DM annual
"""
log10_Amp_dm1yr = parameter.Uniform(-10, -2)
phase_dm1yr = parameter.Uniform(0, 2*np.pi)

wf = chrom_yearly_sinusoid(log10_Amp=log10_Amp_dm1yr,
                           phase=phase_dm1yr, idx=2)

dm1yr = deterministic_signals.Deterministic(wf, name="dm1yr")

"""
define solar wind model
"""

sw, priors = get_informed_nearth_priors(psr, noisedict_p0015, noisedict_p9985, priors)

"""
DM Gaussian
"""
log10_Amp = parameter.Uniform(-10, -2)
log10_sigma_gauss = parameter.Uniform(0, 3)
epoch_gauss = parameter.Uniform(53800, 54000)

wf = dm_gaussian(log10_Amp=log10_Amp, epoch=epoch_gauss, log10_sigma=log10_sigma_gauss, idx=2)

dmgauss = deterministic_signals.Deterministic(wf, name="dmgauss") 

"""
Define total model by summing all components
"""

"""
Define total model by summing all components
"""
# define model
s = tm + wn + rn

add_gwb = True
if 'gwb' in chain_label:
    add_gwb = True
elif add_gwb:
    chain_label += '_gwb'
if add_gwb:
    s += gwb

chain_label += '_fixwhite_informedprior'

add_jump = False
if add_jump:
    s += achrom_JUMP

do_solar_wind = True
if do_solar_wind:
    s += sw
    
if not dir == 'ppta15':
    s += dm 


"""
Add special noise model components for some pulsars
"""
# Define exponential dip parameters for 0437, 1643, and 1713
if psr.name == 'J1713+0747' and psr.toas.min() < 57500*86400:
    expdip = True
    num_dips = 2
    idx = [parameter.Uniform(1.0, 3.0), parameter.Uniform(0.0, 2.0)]
    tmin = [54650, 57400]  # centred 54750 and 57510
    tmax = [54850, 57600]
elif psr.name == 'J0437-4715' and psr.toas.min() < 57100*86400:
    expdip = True
    num_dips = 1
    idx = [parameter.Uniform(-1.0, 2.0)]
    tmin = [57000]
    tmax = [57200]
elif psr.name == 'J1643-1224' and psr.toas.min() < 57100*86400:
    expdip = True
    num_dips = 1
    idx = [parameter.Uniform(-2.0, 0.0)]
    tmin = [57000]
    tmax = [57200]
elif psr.name == 'J2145-0750' and psr.toas.min() < 56450*86400  and psr.toas.max() > 56300*86400:
    expdip = True
    num_dips = 1
    idx = [parameter.Uniform(-2.0, 2.0)]
    tmin = [56250]#, psr.toas.min()/86400]
    tmax = [56450]#, psr.toas.max()/86400]
else:
    expdip = False

# Add exponential dips to model
if expdip:
    name = ['dmexp_{0}'.format(ii+1) for ii in range(num_dips)]

    for idip in range(num_dips):
        t0_exp = parameter.Uniform(tmin[idip], tmax[idip])
        log10_Amp_exp = parameter.Uniform(-10, -2)
        log10_tau_exp = parameter.Uniform(0, 2.5)
        # Define chromatic exponential decay waveform
        wf = chrom_exp_decay(log10_Amp=log10_Amp_exp,
                             t0=t0_exp, log10_tau=log10_tau_exp,
                             sign_param=-1.0, idx=idx[idip])
        expdip = deterministic_signals.Deterministic(wf, name=name[idip])
        s += expdip
        #if model_comp: s4 += expdip

# Annual DM for J0613
if psr.name == 'J0613-0200':
    s += dm1yr
    #if model_comp: s4 += dm1yr

# Gaussian DM for J1603
if psr.name == 'J1603-7202' and psr.toas.min() < 57500*86400:
    s += dmgauss
    #if model_comp: s4 += dmgauss

# Chromatic noise for several pulsars in Goncharov+ or Lentati+ (1600 and 1643)
if psr.name in ['J0437-4715', 'J0613-0200', 'J1017-7156', 'J1045-4509', 'J1600-3053', 'J1643-1224', 'J1939+2134']:
    s += chrom
    #if model_comp: s4 += chrom

# Excess low-frequency band noise for several pulsars in Goncharov+ or Lentati+
if psr.name in ['J0437-4715', 'J0613-0200', 'J1017-7156', 'J1045-4509', 'J1600-3053', 'J1643-1224', 'J1713+0747', 'J1909-3744', 'J1939+2134']:
    s += bn 
    #if model_comp: s4 += bn

if psr.name in ['J0437-4715']:
    s += bn_mid
    #if model_comp: s4 += bn_mid

# add in group noise
if psr.name in psr_groupnoise_dict_dict[dir].keys():
    print(f"{psr.name} IN GROUP NOISE DICT")
    s += gn
    #if model_comp: s4 += gn
if psr.name in psr_groupecorr_dict_dict[dir].keys():
    ecg = gp_signals.EcorrBasisModel(log10_ecorr = log10_ecorr_prior,
                                     selection = by_group_ecorr_dict[psr.name],
                                     name = 'basis_ecorr_group')
    s += ecg
    #if model_comp: s4 += ecg

# Add some high-frequency (up to 1/30 days) achromatic process
if psr.name in ['J0437-4715', 'J1017-7156', 'J1022+1001', 'J1600-3053', 'J1713+0747', 'J1744-1134', 'J1909-3744', 'J2241-5236']:
    s += hf
    #if model_comp: s4 += hf

if psr.name in ['J0437-4715']:
    s += bn_high
    #if model_comp: s4 += bn_high


#model.append(s(psr))
#if model_comp: model2.append(s4(psr))


"""
Set up your PTA likelihood and your sampler
"""
# set up PTA
nmodels = 1
pta = dict.fromkeys(np.arange(nmodels))
pta[0] = signal_base.PTA(s(psr))
pta[0].set_default_params(noisedict)
hyper_model = hypermodel.HyperModel(pta)

# set initial parameters drawn from prior
x0 = hyper_model.initial_sample()
ndim = len(x0)

params = hyper_model.param_names
chaindir = datadir + '/noiseFiles/chains/'
print(chaindir)
chainfile = sorted(glob.glob(chaindir + psrname + '*_chain.npy'))[0]

informed_sample = True
if informed_sample:
    
    for i in range(0, ndim):
        print("Parameter {}, prior sample {}".format(params[i], x0[i]))
        
        if psrname not in chainfile: continue
        
        #chainfile = glob.glob(chaindir + psrname + '*.npy')[0] #assuming only one chain file per pulsar
        parsfile = chainfile.replace('_chain.npy', '_pars.npy')
        pars = np.load(parsfile)#, dtype =  np.unicode_)
        chain = np.load(chainfile)
        #chain = np.loadtxt(chaindir + psrname + '_999/chain_1.txt')
        #pars = np.loadtxt(chaindir + psrname + '_999/pars.txt'
        
        ind_sample = np.random.randint(0, len(chain[:, 0]))
        ind_par = np.argwhere([params[i] == p for p in pars])
        if len(ind_par) == 0:
            print("Parameter {} not in chain - continuing.".format(params[i]))
            continue
        else:
            print("Grabbing initial sample for parameter {} from noisemodelling chain {}".format(params[i], chainfile))
            sample = chain[ind_sample, ind_par].squeeze()
            print(f'{params[i]} = {sample}')
            x0[i] = sample
            del chain
            print("Parameter {}, posterior sample {}".format(params[i], x0[i]))
        
        
        
        
# sampler for N steps
N = int(1e6)

# output directory:
if reddm_only:
    outdir = datadir + '/chains/singlePsrNoise_reddm/' + psrname + "_" + chainnum
else:
    if not dir == 'ppta15':
        outdir = datadir + '/chains/singlePsrNoise' + chain_label +'/' + psrname + "_" + chainnum
    else:
        outdir = datadir + '/chains/singlePsrNoise' + chain_label + '_gwb' +'/' + psrname + "_" + chainnum

# Use PTMCMC sampler. We can easily update this to use e.g. Bilby instead
sampler = hyper_model.setup_sampler(outdir=outdir, resume=False)

# Print parameter names and write to a file
print(hyper_model.param_names)
filename = outdir + "/pars.txt"
if os.path.exists(filename):
    os.remove(filename)
with open(filename, "a") as f:
    for par in hyper_model.param_names:
        f.write(par + '\n')

# print components to file:
filename = outdir + "/components.txt"
if os.path.exists(filename):
    os.remove(filename)
with open(filename, "a") as f:
    f.write('Red: {}\n'.format(red_components))
    f.write('DM: {}\n'.format(dm_components))
    f.write('Band: {}\n'.format(band_components))
    f.write('Chrom: {}\n'.format(chrom_components))
    f.write('hf: {}\n'.format(hf_components))

# copy current file to output directory
os.system('cp {0} {1}/{0}'.format(sys.argv[0], outdir))

# Sample! The sampler parameters can be left as default.
sampler.sample(x0, N, SCAMweight=30, AMweight=15, DEweight=50)
