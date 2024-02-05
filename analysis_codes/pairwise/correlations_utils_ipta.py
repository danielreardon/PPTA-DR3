#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 16 21:01:48 2023

@author: dreardon
"""

import numpy as np
from scipy.signal import correlate

def acor(arr):
    arr -= np.mean(arr)
    auto_correlation = correlate(arr, arr, mode='full')
    auto_correlation = auto_correlation[auto_correlation.size//2:]
    auto_correlation /= auto_correlation[0]
    indices = np.where(auto_correlation<0.5)[0]
    if len(indices)>0:
        return indices[0]
    else:
        return 0

def hd_orf(theta):
    """Hellings & Downs spatial correlation function."""
    omc2 = (1 - np.cos(theta * np.pi/180)) / 2
    orf = (1/2) - (1/4) * omc2 + (3/2) * omc2 * np.log(omc2)
    orf[theta == 0] = 0.5
    return orf

def monopole(theta):
    """Monopole spatial correlation function."""
    return np.ones(np.shape(theta))

def dipole(theta):
    """Dipole spatial correlation function."""
    return np.cos(theta * np.pi/180)

def plot_violin(ax, pos, chain, width=22.5, colour='darkgreen', alpha=1, edgecolour=None, linewidth=0, ylabel=r'Correlation coefficient, $\Gamma$'):

    ax.set_xlabel(r'Sky separation angle, $\zeta$ (degrees)')
    ax.set_ylabel(ylabel)

    violin_dict = ax.violinplot(chain, positions = [pos], widths=width, showextrema=False, ) #, showextrema = True, showmeans = False, showmedians = False)
    for violin_body in violin_dict['bodies']:
        violin_body.set_alpha(alpha)
        violin_body.set_facecolor(colour)
        violin_body.set_edgecolor(edgecolour)
        violin_body.set_linewidth(linewidth)

    return ax

def get_pairs(psrnames):

    ipair = 1
    pairs = {}

    for i in range(0, len(psrnames)):
        for j in range(0, len(psrnames)):
            if j >= i:
                continue
            pairs[str(ipair)] = [psrnames[i], psrnames[j]]
            ipair += 1

    return pairs

def get_psrnames():
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
    return psrnames













