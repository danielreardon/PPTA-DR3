#!/usr/bin/env python
# coding: utf-8

import os
import sys
import glob
import json
import numpy as np
import subprocess
from enterprise import constants as const
from enterprise.pulsar import Pulsar
from enterprise.signals import signal_base
from enterprise.signals import white_signals
from enterprise.signals import gp_signals
from enterprise.signals import parameter
from enterprise.signals import selections
from enterprise.signals import gp_priors
from enterprise.signals import utils
from enterprise.signals import deterministic_signals
from enterprise.signals import gp_bases
from enterprise_extensions.blocks import common_red_noise_block
# from PTMCMCSampler.PTMCMCSampler import PTSampler as ptmcmc
from enterprise_extensions import hypermodel

from ppta_dr3_models import *
from ppta_dr3_utils import *

"""
USAGE:
python3 commonNoise.py <chainnum> <dir> <crn_name> 

chain_num is the chain index, used for launching multiple identical samplers
dir is either "all", "dr2", or "uwl" for different sub-sections of the DR3 dataset

Select one or more of the <crn_names> listed in the menu below
To use multiple CRN signals, separate them with a comma
e.g.
pl_nocorr_freegam,freespec_dipole
would add both a free-gamma CRN and a dipole-correlated free spectrum

To do model comparison between CRN models, separate with a hyphen:
e.g.
pl_nocorr_freegam-pl_hd_freegam
This would do model comparison for a powerlaw free-gamma uncorrelated CRN vs powerlaw HD-correlated free gamma CRN

This can also be combined with commas so that multiple models can be compared
e.g.
pl_nocorr_freegam,pl_hd_fixgam-pl_hd_freegam

models can be added by addiing an entry into the crn_model_dict in ppta_dr3_models
"""

"""
CRN SELECTIONS ON THE MENU:
crn_model_dict = {
'pl_nocorr_freegam'  :  #powerlaw, no correlation, free spectral index
'pl_nocorr_fixgam'   :  #powerlaw, no correlation, fixed 13/3 spectral index
'bpl_nocorr_freegam' :  #broken power law, no correlation
'freespec_nocorr'    :  #free-spectrum no correlation
'pl_orf_bins'        :  #power law red noise, ORF sampled in angular bins
'pl_orf_singlebin'   :  #power law red noise, ORF sampled in a singl angular bin
'pl_orf_singlebin_fixA':  #power law red noise, ORF sampled in a singl angular bin with fixed log10_A
'pl_orf_spline'      :  #power law red noise, spline-interpolated ORF bins    
'pl_hd_fixgam'       :  #power law red noise, fixed spectral index, Hellings-Downs correlations
'pl_hd_freegam'      :  #power law red noise, free spectral index, Hellings-Downs correlations,
'pl_hdnoauto_fixgam' :  #power law red noise, fixed spectral index, Hellings-Downs correlations excluding auto-correlations
'freespec_hd'        :  #free spectrum red noise, Hellings-Downs correlations
'pl_dipole_freegam'  :  #power law red noise, free spectral index, Dipole correlations,
'freespec_dipole'    :  #free spectrum red noise, Dipole correlations
'pl_monopole_freegam':  #power law red noise, free spectral index, Monopole correlations,
'freespec_monopole'  :  #free spectrum red noise, Monopole correlations
}

"""

# Read in data directory, and chain number
chainnum = sys.argv[1]
dir = sys.argv[2]
crn_name = sys.argv[3]
if "freespec_hdcorr" in crn_name:
    crn_name = crn_name.replace("freespec_hdcorr", "freespec_hd")


if '_slice' in crn_name:
    #then set up the slice index and slice length
    # e.g. _slice1_yr9
    islice = int(crn_name.split('_slice')[-1].split('_yr')[0])
    yearslice = int(crn_name.split('_slice')[-1].split('_yr')[-1])
    crn_name = crn_name.replace('_slice{}_yr{}'.format(islice, yearslice),'')
else:
    islice = None
    yearslice = None

if '_quad_' in crn_name:
    #add a generic quadratic form to the toas
    do_quad = True
else:
    do_quad = False

#take SSE as input argument
ephem = sys.argv[4]

model_comp, CRN_Name = parse_crn_name(crn_name)


bayesephem_1 = False
bayesephem_2 = False
try:
    #parse in bayesephem arguments
    _BE_str_list = sys.argv[5].split('-')
    if 'BE' == _BE_str_list[0]:
        bayesephem_1 = True
    if len(_BE_str_list) > 1:
        if 'BE' == _BE_str_list[1]:
            bayesephem_2 = True
except Exception as e:
    print('here1', e)
    try:
        if 'BE' == sys.argv[5]:
            bayesephem_1 = True
            bayesephem_2 = True
    except Exception as e:
        print('here2', e)


if model_comp: chain_label = '_' + '_'.join(CRN_Name[0]) + '_v_' + '_'.join(CRN_Name[1])
#set up output label for chain
else: chain_label = '_' + '_'.join(CRN_Name[0])


datadir = os.path.abspath("./data/" + str(dir))

#hard-coded bug-fix
if 'all' in dir:
    dir = 'all'
if 'dr2' in dir:
    dir = 'dr2'
if 'uwl' in dir:
    dir = 'uwl'

#find the noise files
dirnoise = os.path.abspath("./data/")
noisefiles = sorted(glob.glob(dirnoise + '/' + str(dir) + '/noiseFiles_maxlike/*.json'))

#set up the noise dictionaries

noisedict = update_tidy_noisedict(noisefiles, dir, dirnoise)

#for key_, val_ in noisedict.items():
#    print(key_, val_)
#print(noisedict)

#now load in the 3-sigma bounds
noisedict_p0015, noisedict_p9985 = get_3sig_noisedict(dir, dirnoise)


# Ignore J1824 and J1741 because they're rubbish
#psrnames = [psrname]
psrnames = ['J0030+0451',
            'J0125-2327',
            'J0437-4715',
            'J0613-0200',
            'J0614-3329',
            'J0711-6830',
            'J0900-3144',
            'J1017-7156',
            'J1022+1001',
            'J1024-0719',
            'J1045-4509',
            'J1125-6014',
            'J1446-4701',
            'J1545-4550',
            'J1600-3053',
            'J1603-7202',
            'J1643-1224',
            'J1713+0747',
            'J1730-2304',
            'J1744-1134',
            'J1832-0836',
            'J1857+0943',
            'J1902-5105',
            'J1909-3744',
            'J1933-6211',
            'J1939+2134',
            'J2124-3358',
            'J2129-5721',
            'J2145-0750',
            'J2241-5236']


#set up the single-pair correlation inference
if np.any(['singlebin' in i for i in  CRN_Name[0]]):
    ipair = 1
    pairs = {}
    for i in range(0, len(psrnames)):
        for j in range(0, len(psrnames)):
            if j >= i:
                continue
            pairs[str(ipair)] = [psrnames[i], psrnames[j]]
            ipair += 1

    psrnames = pairs[str(1 + (int(chainnum) - 1) % 435)]
    #435 = N(N-1)/2 for 30 psrs
    print("USING A SINGLE BIN FOR PULSAR PAIR {}".format('_'.join(psrnames)))



#read in par and tim files
parfiles = []
timfiles = []
# if dir == 'dr2' or dir == 'uwl':
#     suffix = '_dmo'
# else:
suffix = ''
for psrname in psrnames:
    parfiles.append(datadir + '/' + psrname + suffix + '.par')
    timfiles.append(datadir + '/' + psrname + suffix + '.tim')

#set up time-slicing
if islice is not None:
    # Do time slicing:
    min_mjd = 53040 + (islice-1)*365.2425
    max_mjd = min_mjd + yearslice*365.2425
    print("Cutting data from {} to {}".format(min_mjd, max_mjd))
    chain_label += '_{}-{}'.format(round(min_mjd), round(max_mjd))
else:
    min_mjd, max_mjd = (None, None)

psrs = construct_psrs(parfiles, timfiles, islice, yearslice, chainnum, ephem, min_mjd=min_mjd, max_mjd=max_mjd)
        
if len(psrs) == 0:
    raise ValueError("No pulsars left after cutting")

tspan, fundamental_freq = get_tspan_fundamental_freq(psrs)
print('TSpan is {:.3f} days'.format(tspan / 86400.0))
print('Fundamental frequency is {:.3g} Hz'.format(fundamental_freq))


#set up the group noise and group-specific ecorr selections
if not dir=='ppta15':
    # Set up group noise
    psr_groupnoiselist_dict = {psr: None for psr in psrnames}
    psr_groupnoiselist_dict.update(psr_groupnoise_dict_dict[dir])
    psr_groupnoiselist_dict = get_groups_in_toas(psrs, psr_groupnoiselist_dict)
    by_group_dict = {key: selections.Selection(sel_by_group_factory(item)._sel_by_group) for key, item in psr_groupnoiselist_dict.items()}

    # Set up group-selected ecorrs
    psr_groupecorrlist_dict = {psr: None for psr in psrnames}
    psr_groupecorrlist_dict.update(psr_groupecorr_dict_dict[dir])
    psr_groupecorrlist_dict = get_groups_in_toas(psrs, psr_groupecorrlist_dict) 
    by_group_ecorr_dict = {key: selections.Selection(sel_by_group_factory(item)._sel_by_group) for key, item in psr_groupecorrlist_dict.items()}


        
"""
Choose whether to marginalise or sample the timing model
"""
tm = gp_signals.MarginalizingTimingModel(use_svd=True)

"""
Define white noise model
"""
# EFAC "MeasurementNoise" can add equad, but only t2equad - we want the tnequad
efac_prior = parameter.Constant()
wn = white_signals.MeasurementNoise(
    efac=efac_prior,
    selection=by_backend)

# EQUAD - TempoNest definition: sigma = sqrt((efac*sigma_0)**2 + (tnequad)**2)
log10_equad_prior = parameter.Constant()
wn += white_signals.TNEquadNoise(
    log10_tnequad=log10_equad_prior,
    selection=by_backend)

# ECORR - we will swap to "white_signals.EcorrKernelNoise" later
if not dir=='ppta15':
    log10_ecorr_prior = parameter.Constant()
    wn += gp_signals.EcorrBasisModel(
        log10_ecorr=log10_ecorr_prior,
        selection=ecorr_selection)
    if dir == 'uwl' or dir == 'dr2':
        wn += gp_signals.EcorrBasisModel(
            log10_ecorr=log10_ecorr_prior,
            selection=no_selection, name='basis_ecorr_all')
    if dir == 'all':
        wn += gp_signals.EcorrBasisModel(
            log10_ecorr=log10_ecorr_prior,
            selection=global_ecorr_selection, name='basis_ecorr_all')

"""
Set up Common Noise
"""

crn_model_dict = get_crn_model_dict(tspan)
if model_comp:
    crns_1 = crn_model_dict[CRN_Name[0][0]]
    if len(CRN_Name[0]) > 1:
        for crn_name_ in CRN_Name[0][1:]:
            crns_1 += crn_model_dict[crn_name_]

    #for model_ind, crn in enumerate(CRN_Name[1]):
    crns_2 = crn_model_dict[CRN_Name[1][0]]
    if len(CRN_Name[1]) > 1:
        for crn_name_ in CRN_Name[1][1:]:
            crns_2 += crn_model_dict[crn_name_]

    crn = [crns_1, crns_2]
    crn_names = ["-".join([i_ for _, i_ in CRN_Name[0].items()]), "-".join([i_ for _, i_ in CRN_Name[1].items()])]
    crn_name = "_v_".join(crn_names)

else:
    crn_name_ = CRN_Name[0][0].replace('_quad_', '')
    crn = crn_model_dict[crn_name_]
    if len(CRN_Name[0]) > 1:
        for crn_name_ in CRN_Name[0][1:]:
            crn += crn_model_dict[crn_name_.replace('_quad_', '')]

"""
combine global model parameters before adding pulsar-specific
"""
if model_comp:
    s = tm + wn + crn[0]
    s3 = tm + wn + crn[1]
else:
    s = tm + wn + crn

if bayesephem_1 or bayesephem_2:
    beph = deterministic_signals.PhysicalEphemerisSignal(use_epoch_toas=True,
                                                         model='setIII_1980', sat_orb_elements = True)
    if bayesephem_1:
        print("Adding Bayesephem to model 1") 
        s += beph
    if bayesephem_2:
        print("Adding Bayesephem to model 2")
        s3 += beph

model = []
model2 = []
priors = {}

components_dict = {key: [] for key in ['red', 'dm', 'band', 'chrom', 'hf']}

for psr in psrs:
    """
    Define red noise model
    """
    rn_model, rn_lgA_prior, rn_gam_prior, priors = get_informed_rednoise_priors(psr, 'red_noise', noisedict_p0015, noisedict_p9985, priors)

    Tspan = psr.toas.max() - psr.toas.min()  # seconds
    max_cadence = 240.0  # days
    red_components = int(Tspan / (max_cadence*86400))
    components_dict['red'].append(red_components)
    print("Using {} red noise components".format(red_components))
    rn = gp_signals.FourierBasisGP(rn_model, components=red_components,
                                   selection=no_selection, name='red_noise')



    hf_model, hf_lgA_prior, hf_gam_prior, priors = get_informed_rednoise_priors(psr, 'hf_noise', noisedict_p0015, noisedict_p9985, priors)
    
    max_cadence = 30  # days
    hf_components = int(Tspan / (max_cadence*86400))
    components_dict['hf'].append(hf_components)
    print("Using {} hf achromatic noise components".format(hf_components))
    hf = gp_signals.FourierBasisGP(hf_model, components=hf_components,
                                   selection=no_selection, name='hf_noise')


    band_model, band_lgA_prior, band_gam_prior, priors = get_informed_rednoise_priors(psr, 'band_noise_low', noisedict_p0015, noisedict_p9985, priors)

    bandmid_model, bandmid_lgA_prior, bandmid_gam_prior, priors = get_informed_rednoise_priors(psr, 'band_noise_mid', noisedict_p0015, noisedict_p9985, priors)
    bandhigh_model, bandhigh_lgA_prior, bandhigh_gam_prior, priors = get_informed_rednoise_priors(psr, 'band_noise_high', noisedict_p0015, noisedict_p9985, priors)
    
    max_cadence = 60  # days
    band_components = int(Tspan / (max_cadence*86400))
    components_dict['band'].append(band_components)
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
        #setting up group noise priors using the supmax and supmin of all group noise bounds
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
    components_dict['dm'].append(dm_components)
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
    components_dict['chrom'].append(chrom_components)
    print("Using {} Chrom components".format(chrom_components))
    chrom_basis = gp_bases.createfourierdesignmatrix_chromatic(nmodes=chrom_components,
                                                               idx=idx)
    chrom = gp_signals.BasisGP(chrom_model, chrom_basis, name='chrom_gp')

    """achromatic quadratic"""
    #coefficient of quadratic term
    log10_a_quad = parameter.Uniform(-24, -16)
    sgn_a = parameter.Uniform(-1.0, 1.0)
    #coefficient of linear term
    log10_b_quad = parameter.Uniform(-16, -9)
    sgn_b = parameter.Uniform(-1.0, 1.0)
    #coefficient of DC term
    log10_c_quad = parameter.Uniform(-12, -8)
    sgn_c = parameter.Uniform(-1.0, 1.0)

    wf = achrom_tm_quadratic(a=log10_a_quad, b=log10_b_quad, c=log10_c_quad, sgn_a=sgn_a, sgn_b=sgn_b, sgn_c=sgn_c)
    tm_quad = deterministic_signals.Deterministic(wf, name="tm_quad")
    
    
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
    Gaussian 20cm - J1600-3053
    """
    epoch_gauss_20cm = parameter.Uniform(57385, 57785)
    wf = gaussian_20cm(log10_Amp=log10_Amp, epoch=epoch_gauss_20cm, log10_sigma=log10_sigma_gauss, nu1=1000, nu2=2000)
    gauss_20cm = deterministic_signals.Deterministic(wf, name="gauss_20cm")

    """
    Chromatic-Gaussian Gaussian - J1600-3053
    """
    nu0_chrom_gaussian = parameter.Uniform(1000, 2000)
    log10_sigma_nu_gauss = parameter.Uniform(-1, 6)
    wf = gaussian_chrom_gaussian(log10_Amp=log10_Amp, epoch=epoch_gauss_20cm, log10_sigma=log10_sigma_gauss, nu_0=nu0_chrom_gaussian, log10_sigma_nu=log10_sigma_nu_gauss)
    gauss_chrom_gauss = deterministic_signals.Deterministic(wf, name="gauss_chrom_gauss")
    
    """
    Define total model by summing all components
    """
    # define model, add DM variations
    s2 = s + rn  + dm + sw
    if model_comp: s4 = s3 + rn  + dm + sw
    
    """
    Add special noise model components for some pulsars
    """

    if do_quad:
        #add achromatic quadratic term if you want to explicitly marginalise over F1
        s2 += tm_quad
        if model_comp: s4 += tm_quad
    
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
            s2 += expdip
            if model_comp: s4 += expdip
    
    # Annual DM for J0613
    if psr.name == 'J0613-0200':
        s2 += dm1yr
        if model_comp: s4 += dm1yr

    # Gaussian DM for J1603
    if psr.name == 'J1603-7202' and psr.toas.min() < 57500*86400:
        s2 += dmgauss
        if model_comp: s4 += dmgauss

    # 20CM Gaussian bump for J1600:
    bump_1600 = False
    if psr.name == "J1600-3053" and psr.toas.min() < 57585*86400 and bump_1600:
        s2 += gauss_20cm
        if model_comp: s4 += gauss_20cm

    # chromatic-gaussian Gaussian bump for J1600:
    bump_1600_gauss = False
    if psr.name == "J1600-3053" and psr.toas.min() < 57585*86400 and bump_1600_gauss:
        s2 += gauss_chrom_gauss
        if model_comp: s4 += gauss_chrom_gauss
    
    # Chromatic noise for several pulsars in Goncharov+ or Lentati+ (1600 and 1643)
    if psr.name in ['J0437-4715', 'J0613-0200', 'J1017-7156', 'J1045-4509', 'J1600-3053', 'J1643-1224', 'J1939+2134']:
        s2 += chrom
        if model_comp: s4 += chrom
   
    # Excess low-frequency band noise for several pulsars in Goncharov+ or Lentati+
    if psr.name in ['J0437-4715', 'J0613-0200', 'J1017-7156', 'J1045-4509', 'J1600-3053', 'J1643-1224', 'J1713+0747', 'J1909-3744', 'J1939+2134']:
        s2 += bn 
        if model_comp: s4 += bn
    
    if psr.name in ['J0437-4715']:
        s2 += bn_mid
        if model_comp: s4 += bn_mid
    
    # add in group noise
    if  psr_groupnoiselist_dict[psr.name] is not  None:
        print(f"{psr.name} IN GROUP NOISE DICT")
        s2 += gn
        if model_comp: s4 += gn
    if  psr_groupecorrlist_dict[psr.name] is not None:
        ecg = gp_signals.EcorrBasisModel(log10_ecorr = log10_ecorr_prior,
                                         selection = by_group_ecorr_dict[psr.name],
                                         name = 'basis_ecorr_group')
        s2 += ecg
        if model_comp: s4 += ecg
    
    # Add some high-frequency (up to 1/30 days) achromatic process
    if psr.name in ['J0437-4715', 'J1017-7156', 'J1022+1001', 'J1600-3053', 'J1713+0747', 'J1744-1134', 'J1909-3744', 'J2241-5236']:
        s2 += hf
        if model_comp: s4 += hf
    
    if psr.name in ['J0437-4715']:
        s2 += bn_high
        if model_comp: s4 += bn_high
    
        
    model.append(s2(psr))
    if model_comp: model2.append(s4(psr))
    
"""
Set up your PTA likelihood and your sampler
"""
# set up PTA
nmodels = 2 if model_comp else 1
pta = dict.fromkeys(np.arange(nmodels))
pta[0] = signal_base.PTA(model)
pta[0].set_default_params(noisedict)
if model_comp:
    pta[1] = signal_base.PTA(model2)
    pta[1].set_default_params(noisedict)


hyper_model = hypermodel.HyperModel(pta)



params = hyper_model.param_names
#print(params)
chaindir = datadir + '/noiseFiles/chains/'
print(chaindir + '*chain.npy')
chainfiles = sorted([sorted(glob.glob(chaindir + '{}*chain.npy'.format(psr_.name)))[0] for psr_ in psrs])
print('Found {} chain files'.format(len(chainfiles)))

print("")
print("Priors:")
print(priors)
print("")

#setting up initial samples by sampling from single-pulsar noise runs: speed up burn-in
x0 = get_initial_samples_chain(chainfiles, hyper_model, priors, informed_sample = True)
ndim = len(x0)

# sampler for N steps
N = int(1e6)

# output directory:
print(datadir, crn_name, chainnum)

savename = 'commonNoise' + chain_label + '_' + ephem
if yearslice is not None:
    slicedir = f"slice_{yearslice}yr"
else:
    slicedir = ""
if bayesephem_1:
    savename += '_BE'
if bayesephem_2:
    savename += '_BE2'
if not 'J0437-4715' in psrnames and not 'singlebin' in chain_label:
    savename += '_no0437'
if bump_1600:
    savename += '_bump1600'
    
outdir = datadir + '/chains/' + slicedir + '/' + savename + '/' + chainnum
print('saving to {}'.format(outdir))

# Use PTMCMC sampler. We can easily update this to use e.g. Bilby instead
sampler = hyper_model.setup_sampler(outdir=outdir, resume=False)

# Print parameter names and write to a file
filename = outdir + "/pars.txt"
if os.path.exists(filename):
    os.remove(filename)
with open(filename, "a") as f:
    for par in hyper_model.param_names:
        f.write(par + '\n')
filename = outdir + "/custom_priors.json"
with open(filename, 'w') as f:
    json.dump(priors, f)

# print components to file:
filename = outdir + "/components.json"
if os.path.exists(filename):
    os.remove(filename)
with open(filename, "w") as f:
    json.dump(components_dict, f)
    # f.write('Red: {}\n'.format(red_components))
    # f.write('DM: {}\n'.format(dm_components))
    # f.write('Band: {}\n'.format(band_components))
    # f.write('Chrom: {}\n'.format(chrom_components))
    # f.write('hf: {}\n'.format(hf_components))

# copy current file to output directory
os.system('cp {0} {1}/{0}'.format(os.path.basename(sys.argv[0]), outdir))

# Sample! The sampler parameters can be left as default.
sampler.sample(x0, N, SCAMweight=30, AMweight=15, DEweight=50)
