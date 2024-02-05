#!/usr/bin/env python
# coding: utf-8

import os
import sys
import glob
import json
import types
import numpy as np
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
from enterprise_extensions.blocks import common_red_noise_block as crn_block
# from PTMCMCSampler.PTMCMCSampler import PTSampler as ptmcmc
from enterprise_extensions import hypermodel, model_orfs
from enterprise_extensions.frequentist.optimal_statistic import OptimalStatistic as OS
from enterprise_extensions.chromatic.solar_wind import solar_wind, createfourierdesignmatrix_solar_dm
from enterprise.signals.parameter import function

@function
def createfourierdesignmatrix_red_ppta(
    toas, nmodes=30, Tspan=None, logf=False, fmin=None, fmax=None, pshift=False, modes=None, pseed=None
):
    """
    Construct fourier design matrix from eq 11 of Lentati et al, 2013
    :param toas: vector of time series in seconds
    :param nmodes: number of fourier coefficients to use
    :param freq: option to output frequencies
    :param Tspan: option to some other Tspan
    :param logf: use log frequency spacing
    :param fmin: lower sampling frequency
    :param fmax: upper sampling frequency
    :param pshift: option to add random phase shift
    :param pseed: option to provide phase shift seed
    :param modes: option to provide explicit list or array of
                  sampling frequencies
    :return: F: fourier design matrix
    :return: f: Sampling frequencies
    """

    T = Tspan if Tspan is not None else toas.max() - toas.min()

    # define sampling frequencies
    if modes is not None:
        nmodes = len(modes)
        f = modes
    elif fmin is None and fmax is None and not logf:
        # make sure partially overlapping sets of modes
        # have identical frequencies
        f = 1.0 * np.arange(1, nmodes + 1) / T
    else:
        # more general case

        if fmin is None:
            fmin = 1 / T

        if fmax is None:
            fmax = nmodes / T

        if nmodes is None:
            nmodes = int(T * fmax)
            print('HERE! Using {} components'.format(nmodes))

        if logf:
            f = np.logspace(np.log10(fmin), np.log10(fmax), nmodes)
        else:
            f = np.linspace(fmin, fmax, nmodes)

    # if requested, add random phase shift to basis functions
    if pshift or pseed is not None:
        if pseed is not None:
            # use the first toa to make a different seed for every pulsar
            seed = int(toas[0] / 17) + int(pseed)
            np.random.seed(seed)

        ranphase = np.random.uniform(0.0, 2 * np.pi, nmodes)
    else:
        ranphase = np.zeros(nmodes)

    Ffreqs = np.repeat(f, 2)

    N = len(toas)
    F = np.zeros((N, 2 * nmodes))

    # The sine/cosine modes
    F[:, ::2] = np.sin(2 * np.pi * toas[:, None] * f[None, :] + ranphase[None, :])
    F[:, 1::2] = np.cos(2 * np.pi * toas[:, None] * f[None, :] + ranphase[None, :])

    return F, Ffreqs

@signal_base.function
def singlebin_orf(pos1, pos2, param):
    '''
    used for inferring the correlation for a single pair of pulsars.
    param is the correlation value for the pair, passed in as a Uniform distr.
    '''
    if np.all(pos1 == pos2):
        return 1
    else:
        return param

@signal_base.function
def zero_diag_crn(pos1, pos2):
    """
    Off-diagonal uncorrelated CRN correlation function (i.e. correlation = 0)
    Explicitly sets cross-correlation terms to 0. Auto terms are 1 (do not run with additional CURN term)
    """
    if np.all(pos1 == pos2):
        return 1e-20
    else:
        return 1e-20
    
    
def FourierBasisGP_ppta(
    spectrum,
    coefficients=False,
    combine=True,
    components=None,
    selection=selections.Selection(selections.no_selection),
    Tspan=None,
    modes=None,
    name="red_noise",
    pshift=False,
    pseed=None,
    fmax=None
):
    """Convenience function to return a BasisGP class with a
    fourier basis."""

    #now feeding fmax as kwarg to createfourierdesignmatrix_red:
    basis = createfourierdesignmatrix_red_ppta(nmodes=components, Tspan=Tspan, modes=modes, pshift=pshift, pseed=pseed, fmax=fmax)
    BaseClass = gp_signals.BasisGP(spectrum, basis, coefficients=coefficients, combine=combine, selection=selection, name=name)

    class FourierBasisGP(BaseClass):
        signal_type = "basis"
        signal_name = "red noise"
        signal_id = name

    return FourierBasisGP

@signal_base.function
def achrom_tm_quadratic(toas, a=-22, b=-12, c=-10, sgn_a=0.5, sgn_b=0.5, sgn_c=0.5):
    """
    general achromatic quadratic waveform
    designed to marginalise over F0 and F1
    """
    t1 = (10**a)
    t1 *= np.sign(sgn_a)
    t2 = (10**b)
    t2 *= np.sign(sgn_b)
    t3 = (10**c)
    t3 *= np.sign(sgn_c)
    quad = t1 + t2 + t3
    wf = quad
    return wf


@signal_base.function
def step_achrom_jump(toas, log10_Amp=-7, sign_param=0.0, t0=54000):
    """
    Achromatic JUMP in ToAs.
    :param t0: time of JUMP [MJD]
    :param log10_Amp: amplitude of JUMP
    :param sign_param: sign of waveform
    :return wf: delay time-series [s]
    """
    t0 *= const.day
    wf = 10**log10_Amp * np.heaviside(toas - t0, 1)
    
    return np.sign(sign_param) * wf

@signal_base.function
def chrom_exp_decay(toas, freqs, log10_Amp=-7, sign_param=-1.0,
                    t0=54000, log10_tau=1.7, idx=2, ref_freq=1400):
    """
    Chromatic exponential-dip delay term in TOAs.
    :param t0: time of exponential minimum [MJD]
    :param tau: 1/e time of exponential [s]
    :param log10_Amp: amplitude of dip
    :param sign_param: sign of waveform
    :param idx: index of chromatic dependence
    :param ref_freq: reference frequency in MHz
    :return wf: delay time-series [s]
    """
    t0 *= const.day
    tau = 10**log10_tau * const.day
    ind = np.where(toas > t0)[0]
    wf = 10**log10_Amp * np.heaviside(toas - t0, 1)
    wf[ind] *= np.exp(- (toas[ind] - t0) / tau)

    return np.sign(sign_param) * wf * (ref_freq / freqs) ** idx


@signal_base.function
def chrom_yearly_sinusoid(toas, freqs, log10_Amp=-7, phase=0, idx=2):
    """
    Chromatic annual sinusoid.
    :param log10_Amp: amplitude of sinusoid
    :param phase: initial phase of sinusoid
    :param idx: index of chromatic dependence
    :return wf: delay time-series [s]
    """

    wf = 10**log10_Amp * np.sin(2 * np.pi * const.fyr * toas + phase)
    return wf * (1400 / freqs) ** idx

@signal_base.function
def dm_gaussian(toas, freqs, log10_Amp=-7, epoch=53910, log10_sigma=2, idx=2):
    """
    DM Gaussian.
    :param log10_Amp: amplitude of gaussian
    :param epoch: epoch of DM event centre
    :param log10_sigma: width of DM event
    :return wf: delay time-series [s]
    """

    wf = 10**log10_Amp * np.exp( - (toas - epoch*86400)**2 / (2* (10**log10_sigma * 86400)**2) )
    return wf * (1400 / freqs) ** 2

@signal_base.function
def gaussian_20cm(toas, freqs, log10_Amp=-7, epoch=57585, log10_sigma=2, nu1 = 1000, nu2 = 2000):
    """
    20cm wavelength Gaussian bump.
    :param log10_Amp: amplitude of gaussian
    :param epoch: epoch of event centre
    :param log10_sigma: width of event
    :return wf: delay time-series [s]
    """

    wf = 10**log10_Amp * np.exp( - (toas - epoch*86400)**2 / (2* (10**log10_sigma * 86400)**2) )
    #modify only toas between nu1 and nu2
    tophat = ((freqs > nu1 ) * (freqs < nu2)).astype(int)
    return wf * tophat

@signal_base.function
def gaussian_chrom_gaussian(toas, freqs, log10_Amp=-7, epoch=57585, log10_sigma=2, nu_0=1600, log10_sigma_nu=2):
    """
    chromatic-gaussian envelope Gaussian temporal bump.
    :param log10_Amp: amplitude of gaussian
    :param epoch: epoch of event centre
    :param log10_sigma: width of event
    :param nu_0: central frequency of chromatic gaussain envelope
    :param log10_sigma_nu: width of chromatic gaussian envelope
    :return wf: delay time-series [s]
    """

    wf = 10**log10_Amp * np.exp( - (toas - epoch*86400)**2 / (2* (10**log10_sigma * 86400)**2) )
    #set up a normalised chromatic gaussian envelope
    chrom_g = np.exp(-(freqs - nu_0)**2.0 / (2.0 * (10**log10_sigma_nu)**2.0))
    #modify the ToAs with the envelope
    return wf * chrom_g

def dm_annual_signal(idx=2, name='dm_s1yr'):
    """
    Returns chromatic annual signal (i.e. TOA advance):
    :param idx:
        index of radio frequency dependence (i.e. DM is 2). If this is set
        to 'vary' then the index will vary from 1 - 6
    :param name: Name of signal
    :return dm1yr:
        chromatic annual waveform.
    """
    log10_Amp_dm1yr = parameter.Uniform(-10, -2)
    phase_dm1yr = parameter.Uniform(0, 2*np.pi)

    wf = chrom_yearly_sinusoid(log10_Amp=log10_Amp_dm1yr,
                               phase=phase_dm1yr, idx=idx)
    dm1yr = deterministic_signals.Deterministic(wf, name=name)

    return dm1yr

def get_groups_in_toas(psrs, groupsel_dict):
    """take in groupnoise list dict and return same dict without groups that are missing in the toas"""
    for psr in psrs:
        if len(psr.toas) == 0:
            continue
        
        grouplist = groupsel_dict[psr.name]
        if grouplist is None:
            continue
        grouplist_new = []
        for i_, group in enumerate(grouplist):
            #ignore groups with less than 10 ToAs
            print('HERE', psr.name, group, np.sum(psr.flags["group"] == group))
            if np.sum(psr.flags["group"] == group) <= 10:
                continue
            if group in np.unique(psr.flags["group"]):
                grouplist_new.append(group)
                
        if len(grouplist_new) == 0:
            grouplist_new = None
        groupsel_dict[psr.name] = grouplist_new
    return(groupsel_dict)
    #for psrname_, grouplist in d.items():
        
        


def delete_empty_keys(d):
    """
    function to delete any empty keys in a group flag dictionary
    """
    delkeys = []
    for key in d.keys():
        if np.sum(d[key]) == 0:  # all False
            delkeys.append(key)
    for key in delkeys:
        del d[key]
    return d


def sel_by_group(flags, flagvals = None):
    """Selection function to split by PPTA frequency band under group flag"""
    if flagvals is None:
        
        flagvals = np.unique(flags["group"])
    else:
        if isinstance(flagvals, list) or isinstance(flagvals, np.ndarray):
            flagvals = np.unique(flagvals)
        else:
            raise ValueError('kwarg flagvals should be list or np.ndarray, but is {}'.format(type(flagvals)))
    
    return delete_empty_keys({val: flags["group"] == val for val in flagvals})

        
class sel_by_group_factory:
    """
    class that is instantiated with a list of group flags.
    Used for making group and band noise selections
    """
    def __init__(self, flagvals = None):
        if isinstance(flagvals, list) or isinstance(flagvals, np.ndarray):
            self.flagvals = np.unique(flagvals)
        elif flagvals == None:
            self.flagvals = flagvals
        else:
            raise ValueError(f'flagvals should be None or array-like, not {type(self.flagvals)}')
        
    def _sel_by_group(self, flags):
        if self.flagvals is None:
            #return all flags
            return delete_empty_keys({val: [val in g for g in flags["group"]] for val in np.unique(flags["group"])})
        else:
            return delete_empty_keys({val: np.array([val in g for g in flags["group"]]) for val in self.flagvals})
    

        
class get_toa_portion_factory:
    """selection function factory to split data into portion from MJD 'mjd_min' to MJD 'mjd_max'"""

    def __init__(self, mjd_min = None, mjd_max = None):
        #we turn single values of mjd_min, mjd_max into an array for easier handling later on
        if isinstance(mjd_min, int) or isinstance(mjd_min, float):
            self.mjd_min = np.floor([mjd_min])
        if isinstance(mjd_max, int) or isinstance(mjd_max, float):
            self.mjd_max = np.ceil([mjd_max])

        if mjd_min is None:
            if mjd_max is None:
                self.mjd_min = [None]
                self.mjd_max = [None]
            else:
                self.mjd_min = [None for m in mjd_max]
                
        elif mjd_max is None:
            self.mjd_max = [None for m in mjd_min]
            
        #if it's a list, ensure same length and rint them
        if isinstance(mjd_min, list) or isinstance(mjd_min, np.ndarray):
            if len(mjd_min) != len(mjd_max):
                raise ValueError('mjd_min (length {}) differs in length from mjd_max (length {})'.format(len(mjd_min), len(mjd_max)))
            self.mjd_min = np.array(mjd_min)
            self.mjd_max = np.array(mjd_max)
        
            
    def get_toa_portion_(self, toas):
        """Selection function to split by data segment"""

        #if either min or max not specified, default to toa min/max
        self.mjd_min = np.floor([m if m is not None else toas.min()/86400.0 for m in self.mjd_min])
        self.mjd_max = np.ceil([m if m is not None else toas.max()/86400.0 for m in self.mjd_max])
        return delete_empty_keys(dict(zip(["{}-{}".format(int(_mjd_min), int(_mjd_max)) for _mjd_min, _mjd_max in zip(self.mjd_min, self.mjd_max)],
                        [(toas >= _mjd_min) * (toas < _mjd_max) for _mjd_min, _mjd_max in zip(self.mjd_min, self.mjd_max)])))
    

def low_frequencies(freqs):
    """Selection for obs frequencies <=960MHz"""
    return delete_empty_keys(dict(zip(['low'], [freqs <= 960])))

def mid_frequencies(freqs):
    """Selection for obs 960MHz < frequencies < 2048MHz"""
    return delete_empty_keys(dict(zip(['mid'], [(freqs > 960) * (freqs < 2048)])))

def high_frequencies(freqs):
    """Selection for obs frequencies >=2048MHz"""
    return delete_empty_keys(dict(zip(['high'], [freqs >= 2048])))

def high_mid_frequencies(freqs):
    """Selection for obs 960MHz < frequencies - to replicate 10-20cm paired band noise"""
    return delete_empty_keys(dict(zip(['mid'], [freqs > 960])))

def global_ecorr(backend_flags):
    """Selection for only UWL data"""
    d = dict(zip(['dr2', 'uwl'],[~np.array(['UWL' in val for val in backend_flags]), np.array(['UWL' in val for val in backend_flags])]))
    delkeys = []
    for key in d.keys():
        if np.sum(d[key]) == 0:  # all False
            delkeys.append(key)
    for key in delkeys:
        del d[key]
    return d

def uwl_all(backend_flags):
    """Selection for only UWL data"""
    return delete_empty_keys({'uwl': np.array(['UWL' in val for val in backend_flags])})



def band_split(freqs, backend_flags):
    """ Selection for splitting the band in 3"""
    arr = np.array(['UWL' in val for val in backend_flags]).squeeze()
    d =  dict(zip(['40CM', '20CM', '10CM', '40CM_uwl', '20CM_uwl', '10CM_uwl'],
                  [(freqs < 960) * ~arr, (960 < freqs) * (freqs < 2048) * ~arr, (2048 < freqs) * (freqs < 4032) * ~arr,
                   (freqs < 960) *  arr, (960 < freqs) * (freqs < 2048) *  arr, (2048 < freqs) * (freqs < 4032) *  arr]))
    delkeys = []
    for key in d.keys():
        if np.sum(d[key]) == 0:  # all False
            delkeys.append(key)
    for key in delkeys:
        del d[key]
    return delete_empty_keys(d)


#def uwl_all(backend_flags):
#    """Selection for only UWL data"""
#    return {'uwl': ['UWL' in val for val in backend_flags]}
def not_uwl_all(backend_flags):
    """Selection for all not-UWL data"""
    return delete_empty_keys({'uwl': np.array(['UWL' not in val for val in backend_flags])})



#@function
def get_achrom_jump(toas, jump_search = False):
    """
    function to return achrom_JUMP signal.
    Can either be tailored for UWL or DR2 portions, or do a jump search in 243-day windows
    """
    if jump_search == True:
        MJD_JUMP_boundaries = np.arange(50000,60000,182.625)
        start_ind = np.squeeze(np.argwhere(toas.min()/86400.0 > MJD_JUMP_boundaries)[-1])
        end_ind = np.squeeze(np.argwhere(toas.max()/86400.0 < MJD_JUMP_boundaries)[0])
        print(start_ind, end_ind)
        log10_Amp_JUMP = parameter.Uniform(-10,-6)
        t0_JUMP = parameter.Uniform(toas.min()/86400.0, toas.max()/86400.0)
        signpar_JUMP = parameter.Uniform(-1, 1)
        
        wf = step_achrom_jump(log10_Amp=log10_Amp_JUMP, sign_param=signpar_JUMP, t0=t0_JUMP)
        jump_wfs = [wf]
        jump_labels = ['ALL']
        for _ind in range(start_ind, end_ind):
            _tmin = MJD_JUMP_boundaries[_ind]
            _tmax = MJD_JUMP_boundaries[_ind + 1]
            log10_Amp_JUMP_ = parameter.Uniform(-10,-6)
            t0_JUMP_ = parameter.Uniform(_tmin - 30.0, _tmax + 30.0) #enable some overlap between JUMP search blocks
            signpar_JUMP_ = parameter.Uniform(-1, 1)
            jump_wfs.append(step_achrom_jump(log10_Amp=log10_Amp_JUMP_, sign_param=signpar_JUMP_, t0=t0_JUMP_))
            jump_labels.append('{}'.format(_ind))
            achrom_JUMP = deterministic_signals.Deterministic(jump_wfs[0], name=f'MJD_JUMP_{jump_labels[0]}')
        for i, jump_wf_ in enumerate(jump_wfs[1:]):
            print(i+1, jump_labels[i+1])
            achrom_JUMP += deterministic_signals.Deterministic(jump_wf_, name=f'MJD_JUMP_{jump_labels[i+1]}')
    else:
        if dir == 'uwl':
            log10_Amp_JUMP = parameter.Uniform(-10,-6)
            t0_JUMP_1 = parameter.Uniform(58256, 59200)
            t0_JUMP_2 = parameter.Uniform(59200, 59645)
            signpar_JUMP = parameter.Uniform(-1, 1)
            wf = step_achrom_jump(log10_Amp=log10_Amp_JUMP, sign_param=signpar_JUMP, t0=t0_JUMP_1)
            achrom_JUMP = deterministic_signals.Deterministic(wf, name='MJD_JUMP_1')
            wf = step_achrom_jump(log10_Amp=log10_Amp_JUMP, sign_param=signpar_JUMP, t0=t0_JUMP_2)
            achrom_JUMP += deterministic_signals.Deterministic(wf, name='MJD_JUMP_2')
        elif dir == 'dr2':
            log10_Amp_JUMP = parameter.Uniform(-10,-6)
            t0_JUMP_1 = parameter.Uniform(53040, 59144)
            signpar_JUMP = parameter.Uniform(-1, 1)
            wf = step_achrom_jump(log10_Amp=log10_Amp_JUMP, sign_param=signpar_JUMP, t0=t0_JUMP_1)
            achrom_JUMP = deterministic_signals.Deterministic(wf, name='MJD_JUMP_1')
        else:
            log10_Amp_JUMP = parameter.Uniform(-10,-6)
            t0_JUMP_1 = parameter.Uniform(58256, 59200)
            t0_JUMP_2 = parameter.Uniform(59200, 59645)
            t0_JUMP_3 = parameter.Uniform(53040, 58256)
            signpar_JUMP = parameter.Uniform(-1, 1)
            wf = step_achrom_jump(log10_Amp=log10_Amp_JUMP, sign_param=signpar_JUMP, t0=t0_JUMP_1)
            achrom_JUMP = deterministic_signals.Deterministic(wf, name='MJD_JUMP_1')
            wf = step_achrom_jump(log10_Amp=log10_Amp_JUMP, sign_param=signpar_JUMP, t0=t0_JUMP_2)
            achrom_JUMP += deterministic_signals.Deterministic(wf, name='MJD_JUMP_2')
            wf = step_achrom_jump(log10_Amp=log10_Amp_JUMP, sign_param=signpar_JUMP, t0=t0_JUMP_3)
            achrom_JUMP += deterministic_signals.Deterministic(wf, name='MJD_JUMP_3')
    return achrom_JUMP



def red_noise_block(psd='powerlaw', prior='log-uniform', Tspan=None,
                    components=30, gamma_val=None, coefficients=False,
                    select=None, modes=None, wgts=None, combine=True,
                    break_flat=False, break_flat_fq=None,
                    logmin=None, logmax=None, dropout=False, k_threshold=0.5,
                    logrhomin=None, logrhomax=None):
    """
    MODIFIED TO TAKE RHOMIN RHOMAX PRIORS AS KWARG FOR FREE-SPECTRUM ANALYSIS.
    1824-2452A was hitting default prior boundary for the lowest frequency bins
    Returns red noise model:
        Red noise modeled as a power-law with 30 sampling frequencies
    :param psd:
        PSD function [e.g. powerlaw (default), turnover, spectrum, tprocess]
    :param prior:
        Prior on log10_A. Default if "log-uniform". Use "uniform" for
        upper limits.
    :param Tspan:
        Sets frequency sampling f_i = i / Tspan. Default will
        use overall time span for indivicual pulsar.
    :param components:
        Number of frequencies in sampling of red noise
    :param gamma_val:
        If given, this is the fixed slope of the power-law for
        powerlaw, turnover, or tprocess red noise
    :param coefficients: include latent coefficients in GP model?
    :param dropout: Use a dropout analysis for intrinsic red noise models.
        Currently only supports power law option.
    :param k_threshold: Threshold for dropout analysis.
    """
    # red noise parameters that are common
    if psd in ['powerlaw', 'powerlaw_genmodes', 'turnover',
               'tprocess', 'tprocess_adapt']:
        # parameters shared by PSD functions
        if logmin is not None and logmax is not None:
            if prior == 'uniform':
                log10_A = parameter.LinearExp(logmin, logmax)
            elif prior == 'log-uniform':
                log10_A = parameter.Uniform(logmin, logmax)
        else:
            if prior == 'uniform':
                log10_A = parameter.LinearExp(-18, -11)
            elif prior == 'log-uniform' and gamma_val is not None:
                if np.abs(gamma_val - 4.33) < 0.1:
                    log10_A = parameter.Uniform(-18, -11)
                else:
                    log10_A = parameter.Uniform(-18, -11)
            else:
                log10_A = parameter.Uniform(-18, -11)

        if gamma_val is not None:
            gamma = parameter.Constant(gamma_val)
        else:
            gamma = parameter.Uniform(0, 7)

        # different PSD function parameters
        if psd == 'powerlaw' and dropout:
            k_drop = parameter.Uniform(0, 1)
            pl = drop.dropout_powerlaw(log10_A=log10_A, gamma=gamma,
                                       dropout_psr='all', k_drop=k_drop,
                                       k_threshold=k_threshold)
        elif psd == 'powerlaw':
            pl = utils.powerlaw(log10_A=log10_A, gamma=gamma)
        elif psd == 'powerlaw_genmodes':
            pl = gp_priors.powerlaw_genmodes(log10_A=log10_A, gamma=gamma, wgts=wgts)
        elif psd == 'turnover':
            kappa = parameter.Uniform(0, 7)
            lf0 = parameter.Uniform(-9, -7)
            pl = utils.turnover(log10_A=log10_A, gamma=gamma,
                                lf0=lf0, kappa=kappa)
        elif psd == 'tprocess':
            df = 2
            alphas = gp_priors.InvGamma(df/2, df/2, size=components)
            pl = gp_priors.t_process(log10_A=log10_A, gamma=gamma, alphas=alphas)
        elif psd == 'tprocess_adapt':
            df = 2
            alpha_adapt = gp_priors.InvGamma(df/2, df/2, size=1)
            nfreq = parameter.Uniform(-0.5, 10-0.5)
            pl = gp_priors.t_process_adapt(log10_A=log10_A, gamma=gamma,
                                     alphas_adapt=alpha_adapt, nfreq=nfreq)

    if psd == 'spectrum':
        if logrhomin is None:
            logrhomin = -10
        if logrhomax is None:
            logrhomax = -4
        if prior == 'uniform':
            log10_rho = parameter.LinearExp(logrhomin, logrhomax, size=components)
        elif prior == 'log-uniform':
            log10_rho = parameter.Uniform(logrhomin, logrhomax, size=components)

        pl = gp_priors.free_spectrum(log10_rho=log10_rho)

    if select == 'backend':
        # define selection by observing backend
        selection = selections.Selection(selections.by_backend)
    elif select == 'band' or select == 'band+':
        # define selection by observing band
        selection = selections.Selection(selections.by_band)
    elif isinstance(select, types.FunctionType):
        selection = selections.Selection(select)
    else:
        # define no selection
        selection = selections.Selection(selections.no_selection)

    if break_flat:
        log10_A_flat = parameter.Uniform(-18, -11)
        gamma_flat = parameter.Constant(0)
        pl_flat = utils.powerlaw(log10_A=log10_A_flat, gamma=gamma_flat)

        freqs = 1.0 * np.arange(1, components+1) / Tspan
        components_low = sum(f < break_flat_fq for f in freqs)
        if components_low < 1.5:
            components_low = 2

        rn = gp_signals.FourierBasisGP(pl, components=components_low,
                                       Tspan=Tspan, coefficients=coefficients,
                                       combine=combine, selection=selection)

        rn_flat = gp_signals.FourierBasisGP(pl_flat,
                                            modes=freqs[components_low:],
                                            coefficients=coefficients,
                                            selection=selection,
                                            combine=combine,
                                            name='red_noise_hf')
        rn = rn + rn_flat
    else:
        rn = gp_signals.FourierBasisGP(pl, components=components,
                                       Tspan=Tspan,
                                       combine=combine,
                                       coefficients=coefficients,
                                       selection=selection,
                                       modes=modes)

    if select == 'band+':  # Add the common component as well
        rn = rn + gp_signals.FourierBasisGP(pl, components=components,
                                            Tspan=Tspan, combine=combine,
                                            coefficients=coefficients)

    return rn


def common_red_noise_block(psd='powerlaw', prior='log-uniform',
                           Tspan=None, components=30, combine=True,
                           log10_A_val=None, gamma_val=None, delta_val=None,
                           logmin=None, logmax=None, select = None,
                           orf=None, orf_ifreq=0, leg_lmax=5,
                           name='gw', coefficients=False,
                           pshift=False, pseed=None):
    """
    MODIFIED TO TAKE `SELECT` AS A KWARG
    
    Returns common red noise model:

        1. Red noise modeled with user defined PSD with
        30 sampling frequencies. Available PSDs are
        ['powerlaw', 'turnover' 'spectrum']

    :param psd:
        PSD to use for common red noise signal. Available options
        are ['powerlaw', 'turnover' 'spectrum', 'broken_powerlaw']
    :param prior:
        Prior on log10_A. Default if "log-uniform". Use "uniform" for
        upper limits.
    :param Tspan:
        Sets frequency sampling f_i = i / Tspan. Default will
        use overall time span for individual pulsar.
    :param log10_A_val:
        Value of log10_A parameter for fixed amplitude analyses.
    :param gamma_val:
        Value of spectral index for power-law and turnover
        models. By default spectral index is varied of range [0,7]
    :param delta_val:
        Value of spectral index for high frequencies in broken power-law
        and turnover models. By default spectral index is varied in range [0,7].\
    :param logmin:
        Specify the lower bound of the prior on the amplitude for all psd but 'spectrum'.
        If psd=='spectrum', then this specifies the lower prior on log10_rho_gw
    :param logmax:
        Specify the lower bound of the prior on the amplitude for all psd but 'spectrum'.
        If psd=='spectrum', then this specifies the lower prior on log10_rho_gw
    :param orf:
        String representing which overlap reduction function to use.
        By default we do not use any spatial correlations. Permitted
        values are ['hd', 'dipole', 'monopole'].
    :param orf_ifreq:
        Frequency bin at which to start the Hellings & Downs function with
        numbering beginning at 0. Currently only works with freq_hd orf.
    :param leg_lmax:
        Maximum multipole of a Legendre polynomial series representation
        of the overlap reduction function [default=5]
    :param pshift:
        Option to use a random phase shift in design matrix. For testing the
        null hypothesis.
    :param pseed:
        Option to provide a seed for the random phase shift.
    :param name: Name of common red process

    """

    
    orfs = {'crn': None, 'hd': model_orfs.hd_orf(),
            'gw_monopole': model_orfs.gw_monopole_orf(),
            'gw_dipole': model_orfs.gw_dipole_orf(),
            'st': model_orfs.st_orf(),
            'gt': model_orfs.gt_orf(tau=parameter.Uniform(-1.5, 1.5)('tau')),
            'dipole': model_orfs.dipole_orf(),
            'monopole': model_orfs.monopole_orf(),
            'param_hd': model_orfs.param_hd_orf(a=parameter.Uniform(-1.5, 3.0)('gw_orf_param0'),
                                                b=parameter.Uniform(-1.0, 0.5)('gw_orf_param1'),
                                                c=parameter.Uniform(-1.0, 1.0)('gw_orf_param2')),
            'spline_orf': model_orfs.spline_orf(params=parameter.Uniform(-0.9, 0.9, size=7)('gw_orf_spline')),
            'bin_orf': model_orfs.bin_orf(params=parameter.Uniform(-1.0, 1.0, size=7)('gw_orf_bin')),
            'single_bin_orf': singlebin_orf(param=parameter.Uniform(-1.0, 1.0)('gw_orf_bin')),
            'zero_diag_crn': zero_diag_crn(),
            'zero_diag_hd': model_orfs.zero_diag_hd(),
            'zero_diag_bin_orf': model_orfs.zero_diag_bin_orf(params=parameter.Uniform(
                -1.0, 1.0, size=7)('gw_orf_bin_zero_diag')),
            'freq_hd': model_orfs.freq_hd(params=[components, orf_ifreq]),
            'legendre_orf': model_orfs.legendre_orf(params=parameter.Uniform(
                -1.0, 1.0, size=leg_lmax+1)('gw_orf_legendre')),
            'zero_diag_legendre_orf': model_orfs.zero_diag_legendre_orf(params=parameter.Uniform(
                -1.0, 1.0, size=leg_lmax+1)('gw_orf_legendre_zero_diag'))}

    # common red noise parameters
    if psd in ['powerlaw', 'turnover', 'turnover_knee', 'broken_powerlaw']:
        amp_name = '{}_log10_A'.format(name)
        if log10_A_val is not None:
            log10_Agw = parameter.Constant(log10_A_val)(amp_name)

        elif logmin is not None and logmax is not None:
            if prior == 'uniform':
                log10_Agw = parameter.LinearExp(logmin, logmax)(amp_name)
            elif prior == 'log-uniform' and gamma_val is not None:
                if np.abs(gamma_val - 4.33) < 0.1:
                    log10_Agw = parameter.Uniform(logmin, logmax)(amp_name)
                else:
                    log10_Agw = parameter.Uniform(logmin, logmax)(amp_name)
            else:
                log10_Agw = parameter.Uniform(logmin, logmax)(amp_name)

        else:
            if prior == 'uniform':
                log10_Agw = parameter.LinearExp(-18, -11)(amp_name)
            elif prior == 'log-uniform' and gamma_val is not None:
                if np.abs(gamma_val - 4.33) < 0.1:
                    log10_Agw = parameter.Uniform(-18, -14)(amp_name)
                else:
                    log10_Agw = parameter.Uniform(-18, -11)(amp_name)
            else:
                log10_Agw = parameter.Uniform(-18, -11)(amp_name)

        gam_name = '{}_gamma'.format(name)
        if gamma_val is not None:
            gamma_gw = parameter.Constant(gamma_val)(gam_name)
        else:
            gamma_gw = parameter.Uniform(0, 7)(gam_name)

        # common red noise PSD
        if psd == 'powerlaw':
            cpl = utils.powerlaw(log10_A=log10_Agw, gamma=gamma_gw)
        elif psd == 'broken_powerlaw':
            delta_name = '{}_delta'.format(name)
            kappa_name = '{}_kappa'.format(name)
            log10_fb_name = '{}_log10_fb'.format(name)
            kappa_gw = parameter.Uniform(0.01, 0.5)(kappa_name)
            log10_fb_gw = parameter.Uniform(-10, -7)(log10_fb_name)

            if delta_val is not None:
                delta_gw = parameter.Constant(delta_val)(delta_name)
            else:
                delta_gw = parameter.Uniform(0, 7)(delta_name)
            cpl = gp_priors.broken_powerlaw(log10_A=log10_Agw,
                                      gamma=gamma_gw,
                                      delta=delta_gw,
                                      log10_fb=log10_fb_gw,
                                      kappa=kappa_gw)
        elif psd == 'turnover':
            kappa_name = '{}_kappa'.format(name)
            lf0_name = '{}_log10_fbend'.format(name)
            kappa_gw = parameter.Uniform(0, 7)(kappa_name)
            lf0_gw = parameter.Uniform(-9, -7)(lf0_name)
            cpl = utils.turnover(log10_A=log10_Agw, gamma=gamma_gw,
                                 lf0=lf0_gw, kappa=kappa_gw)
        elif psd == 'turnover_knee':
            kappa_name = '{}_kappa'.format(name)
            lfb_name = '{}_log10_fbend'.format(name)
            delta_name = '{}_delta'.format(name)
            lfk_name = '{}_log10_fknee'.format(name)
            kappa_gw = parameter.Uniform(0, 7)(kappa_name)
            lfb_gw = parameter.Uniform(-9.3, -8)(lfb_name)
            delta_gw = parameter.Uniform(-2, 0)(delta_name)
            lfk_gw = parameter.Uniform(-8, -7)(lfk_name)
            cpl = gp_priors.turnover_knee(log10_A=log10_Agw, gamma=gamma_gw,
                                    lfb=lfb_gw, lfk=lfk_gw,
                                    kappa=kappa_gw, delta=delta_gw)

    if psd == 'spectrum':
        rho_name = '{}_log10_rho'.format(name)

        # checking if priors specified, otherwise give default values
        if logmin is None:
            logmin = -9
        if logmax is None:
            logmax = -4

        if prior == 'uniform':
            log10_rho_gw = parameter.LinearExp(logmin, logmax,
                                               size=components)(rho_name)
        elif prior == 'log-uniform':
            log10_rho_gw = parameter.Uniform(logmin, logmax, size=components)(rho_name)

        cpl = gp_priors.free_spectrum(log10_rho=log10_rho_gw)

    if select == 'backend':
        # define selection by observing backend
        selection = selections.Selection(selections.by_backend)
    elif select == 'band' or select == 'band+':
        # define selection by observing band
        selection = selections.Selection(selections.by_band)
    elif isinstance(select, types.FunctionType):
        selection = selections.Selection(select)
    else:
        # define no selection
        selection = selections.Selection(selections.no_selection)
        
    if orf is None:
        crn = gp_signals.FourierBasisGP(cpl, coefficients=coefficients,
                                        components=components, Tspan=Tspan,
                                        name=name, pshift=pshift, pseed=pseed,
                                        selection=selection)
    elif orf in orfs.keys():
        if orf == 'crn':
            crn = gp_signals.FourierBasisGP(cpl, coefficients=coefficients,
                                            components=components, Tspan=Tspan,
                                            name=name, pshift=pshift, pseed=pseed,
                                            selection=selection)
        else:
            crn = gp_signals.FourierBasisCommonGP(cpl, orfs[orf],
                                                  components=components,
                                                  Tspan=Tspan,
                                                  name=name, pshift=pshift,
                                                  pseed=pseed)
    elif isinstance(orf, types.FunctionType):
        crn = gp_signals.FourierBasisCommonGP(cpl, orf,
                                              components=components,
                                              Tspan=Tspan,
                                              name=name, pshift=pshift,
                                              pseed=pseed)
    else:
        raise ValueError('ORF {} not recognized'.format(orf))

    return crn


