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
log10_A_GW_prior = parameter.LinearExp(-20, -11)
gamma_prior = parameter.Uniform(0, 7)
# # powerlaw
rn_model = gp_priors.powerlaw(log10_A=log10_A_prior,
                              gamma=gamma_prior)
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

gwb_model = gp_priors.powerlaw(log10_A=log10_A_GW_prior,
                              gamma=4.33)
print("Using {} gwb components".format(red_components))
gwb = gp_signals.FourierBasisGP(gwb_model, components=red_components,
                               selection=no_selection, name='gwb')

max_cadence = 30  # days
hf_components = int(Tspan / (max_cadence*86400))
print("Using {} hf achromatic noise components".format(hf_components))
hf = gp_signals.FourierBasisGP(rn_model, components=hf_components,
                               selection=no_selection, name='hf_noise')

max_cadence = 60  # days
band_components = int(Tspan / (max_cadence*86400))
bn = gp_signals.FourierBasisGP(rn_model, components=band_components,
                               selection=low_freq, name='band_noise')

bn_mid = gp_signals.FourierBasisGP(rn_model, components=band_components,
                                    selection=mid_freq, name='band_noise_mid')

bn_high = gp_signals.FourierBasisGP(rn_model, components=band_components,
                                    selection=high_freq, name='band_noise_high')

# rn = blocks.red_noise_block(components=30)  # The easy way

"""
Define DM noise model
"""
log10_A_dm_prior = parameter.Uniform(-20, -11)
gamma_dm_prior = parameter.Uniform(0, 7)
dm_model = gp_priors.powerlaw(log10_A=log10_A_dm_prior,
                              gamma=gamma_dm_prior)
# # broken powerlaw
# lfb_dm_prior = parameter.Uniform(-10, -7)  # log10 transition frequency
# kappa_dm_prior = 0.1  # smoothness of transition (Default = 0.1)
# delta_dm_prior = 0  # slope for frequencies > f_break
# dm_model = gp_priors.broken_powerlaw(log10_A=log10_A_dm_prior,
#                                      gamma=gamma_dm_prior,
#                                      delta=delta_dm_prior,
#                                      log10_fb=lfb_dm_prior,
#                                      kappa=kappa_dm_prior)
# components = 120
Tspan = psr.toas.max() - psr.toas.min()  # seconds
max_cadence = 60  # days
dm_components = int(Tspan / (max_cadence*86400))
print("Using {} DM components".format(dm_components))
dm_basis = gp_bases.createfourierdesignmatrix_dm(nmodes=dm_components)
dm = gp_signals.BasisGP(dm_model, dm_basis, name='dm_gp')
# dm = blocks.dm_noise_block(components=30)  # The easy way

"""
Define system noise model
"""

if not dir == 'ppta15':
    max_cadence = 30 #days
    log10_A_gn_prior = parameter.Uniform(-20, -11)
    gamma_gn_prior = parameter.Uniform(0, 7)
    gn_model = gp_priors.powerlaw(log10_A = log10_A_gn_prior, gamma = gamma_gn_prior)
    gn = FourierBasisGP_ppta(gn_model, components=None, fmax=1/(max_cadence*86400),
                             selection = by_group_dict[psr.name],
                             name = 'group_noise')
"""
Define chromatic noise model
"""
log10_A_chrom_prior = parameter.Uniform(-20, -11)
gamma_chrom_prior = parameter.Uniform(0, 7)
chrom_model = utils.powerlaw(log10_A=log10_A_chrom_prior,
                             gamma=gamma_chrom_prior)
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
DM Gaussian
"""
log10_Amp = parameter.Uniform(-10, -2)
log10_sigma_gauss = parameter.Uniform(0, 3)
epoch_gauss = parameter.Uniform(53800, 54000)

wf = dm_gaussian(log10_Amp=log10_Amp, epoch=epoch_gauss, log10_sigma=log10_sigma_gauss, idx=2)

dmgauss = deterministic_signals.Deterministic(wf, name="dmgauss")

"""
SOLAR WIND. By default only models 30 components, regarless of Tspan
"""
n_earth = parameter.Uniform(0, 20)('n_earth')
deter_sw = solar_wind(n_earth=n_earth)
mean_sw = deterministic_signals.Deterministic(deter_sw, name='n_earth')

Tspan = psr.toas.max() - psr.toas.min()
max_cadence = 60
sw_components = int(Tspan / (max_cadence*86400))

log10_A_sw = parameter.Uniform(-10, 1)
gamma_sw = parameter.Uniform(-4, 4)
sw_prior = utils.powerlaw(log10_A=log10_A_sw, gamma=gamma_sw)
sw_basis = createfourierdesignmatrix_solar_dm(nmodes=sw_components, Tspan=Tspan)

sw = mean_sw + gp_signals.BasisGP(sw_prior, sw_basis, name='gp_sw')

"""
Achrom JUMP
"""

jump_search = False
achrom_JUMP = get_achrom_jump(psr.toas, jump_search = jump_search)

# if jump_search == True:
#     MJD_JUMP_boundaries = np.arange(50000,60000,182.625)
#     start_ind = np.squeeze(np.argwhere(psr.toas.min()/86400.0 > MJD_JUMP_boundaries)[-1])
#     end_ind = np.squeeze(np.argwhere(psr.toas.max()/86400.0 < MJD_JUMP_boundaries)[0])
#     print(start_ind, end_ind)
#     log10_Amp_JUMP = parameter.Uniform(-10,-6)
#     t0_JUMP = parameter.Uniform(psr.toas.min()/86400.0, psr.toas.max()/86400.0)
#     signpar_JUMP = parameter.Uniform(-1, 1)

#     wf = step_achrom_jump(log10_Amp=log10_Amp_JUMP, sign_param=signpar_JUMP, t0=t0_JUMP)
#     jump_wfs = [wf]
#     jump_labels = ['ALL']
#     for _ind in range(start_ind, end_ind):
#         _tmin = MJD_JUMP_boundaries[_ind]
#         _tmax = MJD_JUMP_boundaries[_ind + 1]
#         log10_Amp_JUMP_ = parameter.Uniform(-10,-6)
#         t0_JUMP_ = parameter.Uniform(_tmin - 30.0, _tmax + 30.0) #enable some overlap between JUMP search blocks
#         signpar_JUMP_ = parameter.Uniform(-1, 1)
#         jump_wfs.append(step_achrom_jump(log10_Amp=log10_Amp_JUMP_, sign_param=signpar_JUMP_, t0=t0_JUMP_))
#         jump_labels.append('{}'.format(_ind))
#     achrom_JUMP = deterministic_signals.Deterministic(jump_wfs[0], name=f'MJD_JUMP_{jump_labels[0]}')
#     for i, jump_wf_ in enumerate(jump_wfs[1:]):
#         print(i+1, jump_labels[i+1])
#         achrom_JUMP += deterministic_signals.Deterministic(jump_wf_, name=f'MJD_JUMP_{jump_labels[i+1]}')
# else:
#     if dir == 'uwl':
#         log10_Amp_JUMP = parameter.Uniform(-10,-6)
#         t0_JUMP_1 = parameter.Uniform(58256, 59200)
#         #t0_JUMP_2 = parameter.Uniform(58925, 59200)
#         t0_JUMP_2 = parameter.Uniform(59200, 59645)
#         #t0_JUMP_58900 = parameter.Uniform(58800, 59000
#         signpar_JUMP = parameter.Uniform(-1, 1)
#         wf = step_achrom_jump(log10_Amp=log10_Amp_JUMP, sign_param=signpar_JUMP, t0=t0_JUMP_1)
#         achrom_JUMP = deterministic_signals.Deterministic(wf, name='MJD_JUMP_1')
#         wf = step_achrom_jump(log10_Amp=log10_Amp_JUMP, sign_param=signpar_JUMP, t0=t0_JUMP_2)
#         achrom_JUMP += deterministic_signals.Deterministic(wf, name='MJD_JUMP_2')
#     elif dir == 'dr2':
#         log10_Amp_JUMP = parameter.Uniform(-10,-6)
#         t0_JUMP_1 = parameter.Uniform(53040, 59144)
#         signpar_JUMP = parameter.Uniform(-1, 1)
#         wf = step_achrom_jump(log10_Amp=log10_Amp_JUMP, sign_param=signpar_JUMP, t0=t0_JUMP_1)
#         achrom_JUMP = deterministic_signals.Deterministic(wf, name='MJD_JUMP_1')
#     else:
#         log10_Amp_JUMP = parameter.Uniform(-10,-6)
#         t0_JUMP_1 = parameter.Uniform(58256, 59200)
#         t0_JUMP_2 = parameter.Uniform(59200, 59645)
#         t0_JUMP_3 = parameter.Uniform(53040, 58256)
#         signpar_JUMP = parameter.Uniform(-1, 1)
#         wf = step_achrom_jump(log10_Amp=log10_Amp_JUMP, sign_param=signpar_JUMP, t0=t0_JUMP_1)
#         achrom_JUMP = deterministic_signals.Deterministic(wf, name='MJD_JUMP_1')
#         wf = step_achrom_jump(log10_Amp=log10_Amp_JUMP, sign_param=signpar_JUMP, t0=t0_JUMP_2)
#         achrom_JUMP += deterministic_signals.Deterministic(wf, name='MJD_JUMP_2')
#         wf = step_achrom_jump(log10_Amp=log10_Amp_JUMP, sign_param=signpar_JUMP, t0=t0_JUMP_3)
#         achrom_JUMP += deterministic_signals.Deterministic(wf, name='MJD_JUMP_3')
        


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

chain_label += 'lineargw_fixwhite'

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
if psr.name == 'J1713+0747' and psr.toas.min() < 57500*86400 and psr.toas.max() > 54650*86400:
    expdip = True
    num_dips = 2
    idx = [parameter.Uniform(1.0, 3.0), parameter.Uniform(0.0, 2.0)]
    tmin = [54650, 57400]  # centred 54750 and 57510
    tmax = [54850, 57600]
elif psr.name == 'J0437-4715' and psr.toas.min() < 57100*86400  and psr.toas.max() > 57000*86400:
    expdip = True
    num_dips = 1
    idx = [parameter.Uniform(-1.0, 2.0), parameter.Uniform(-7.0, 7.0)]
    tmin = [57000]#, psr.toas.min()/86400]
    tmax = [57200]#, psr.toas.max()/86400]
elif psr.name == 'J1643-1224' and psr.toas.min() < 57100*86400 and psr.toas.max() > 57000*86400:
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
        log10_tau_exp = parameter.Uniform(0, 4)
        # Define chromatic exponential decay waveform
        wf = chrom_exp_decay(log10_Amp=log10_Amp_exp,
                             t0=t0_exp, log10_tau=log10_tau_exp,
                             sign_param=-1.0, idx=idx[idip])
        expdip = deterministic_signals.Deterministic(wf, name=name[idip])
        s += expdip

# Annual DM for J0613
if psr.name == 'J0613-0200':
    s += dm1yr

# Gaussian DM for J1603
if psr.name == 'J1603-7202' and psr.toas.min() < 57500*86400:
    s += dmgauss

if not reddm_only and not dir=='ppta15':

    # Chromatic noise for several pulsars in Goncharov+ or Lentati+ (1600 and 1643)
    if psr.name in ['J0437-4715', 'J0613-0200', 'J1017-7156', 'J1045-4509', 'J1600-3053', 'J1643-1224', 'J1939+2134']:
        s += chrom

    # Excess low-frequency band noise for several pulsars in Goncharov+ or Lentati+
    # cut out 1045 due to insignificant results
    if psr.name in ['J0437-4715', 'J0613-0200', 'J1017-7156', 'J1045-4509', 'J1600-3053', 'J1643-1224', 'J1713+0747', 'J1909-3744', 'J1939+2134']:
        s += bn 

    if psr.name in ['J0437-4715']:
        s += bn_mid

    # add in group noise
    if psr.name in psr_groupnoise_dict_dict[dir].keys():
        print(f"{psr.name} IN GROUP NOISE DICT")
        s += gn

    # Add some high-frequency (up to 1/30 days) achromatic process
    if psr.name in ['J0437-4715', 'J1017-7156', 'J1022+1001', 'J1600-3053', 'J1713+0747', 'J1744-1134', 'J1909-3744', 'J2241-5236']:
        s += hf

    if psr.name in ['J0437-4715']:
        s += bn_high

elif not reddm_only:
    if psr.name in ['J0437-4715', 'J1017-7156', 'J1022+1001', 'J1600-3053', 'J1713+0747', 'J1744-1134', 'J1909-3744', 'J2241-5236']:
        s += hf

    if psr.name in ['J0437-4715']:
        s += bn_high


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

informed_sample = False
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
