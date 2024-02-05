#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  4 14:53:00 2019

@author: dreardon

Runs basic white, red, and DM noise model for all pulsars in datadir
"""

import os
import json
import sys
import glob
import numpy as np
import subprocess
from enterprise_extensions import model_utils
import matplotlib.pyplot as plt
import corner

from matplotlib import rc
import matplotlib
import os
os.environ["PATH"] += os.pathsep + '/Library/TeX/texbin/'
rc('text', usetex=True)
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
#matplotlib.rcParams.update({'font.size': 26})


noise_model = 'singlePsrNoise_AZ_jumpsearch'
dir = sys.argv[1]
cpars = ['dm', 'red', 'gwb']
if dir == 'all':
    datadir = os.path.abspath("./data/[up]*")
else:
    datadir = os.path.abspath("./data/" + str(dir))

nchain = 8 #for PPTA15,  5 is 440 dr3 noise models, 6 is 440 basic noise, 7 is 421 basic noise, 8 is 438 basic noise

psrs = np.loadtxt('./psrs.dat', dtype = 'str')
if 'dr2' in dir:
    psrs = np.loadtxt('./dr2_psrs.dat', dtype = 'str')
    

#chainfiles = sorted(glob.glob(datadir + "/chains/{}/".format(noise_model) + '/J*_' + str(nchain) + '/chain_1.txt'))
njump = int(sys.argv[2])
number = 5#len(psrs)
print(number)
cmap = plt.get_cmap('jet')

cmap = plt.get_cmap('tab20')#gist_rainbow')
colors = [cmap(i) for i in np.linspace(0, 0.8, number)]
#colors = [cmap(i) for i in np.linspace(0, 1, number)]

first_chain_overall = True
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(2.0*3.3, 2.0*1.5*3.3))
psrnames = []
psr_chains = []
psr_pars = []

for psrnum, psrname in enumerate(psrs):
    #print('{}/{}_chain.npy'.format(datadir + "/chains/{}/".format(noise_model), psrname))
    if os.path.exists('{}/{}_chain.npy'.format(datadir + "/chains/{}/".format(noise_model), psrname)):
        outdir = datadir + "/chains/{}/".format(noise_model) + psrname + '_1'
        if not os.path.exists(outdir):
            continue
        chain_total = np.load('{}/{}_chain.npy'.format(datadir + "/chains/{}/".format(noise_model), psrname), allow_pickle = True)
        if chain_total.size > 1:
            pars = np.loadtxt(outdir + '/pars.txt', dtype=np.unicode_)
            psr_chains.append(chain_total)
            psr_pars.append(pars)
            print(psrname)
        #psrnames.append(psrname)
        continue
    else:
        max_nchain = 1  # number of chains to read
        
        first_chain = True
        chain_total = None
        for i in range(1, max_nchain+1): 
            outdir = datadir + "/chains/{}/".format(noise_model) + psrname + '_' + str(i)
            chainfile = outdir + '/chain_1.txt'
            print(chainfile)
            nlines = int(subprocess.check_output(f'cat {chainfile} | wc -l', shell = True, text = True))
            if nlines < 2500:
                print(nlines, "dead")
                continue
            if not os.path.exists(chainfile):
                print(chainfile, 'non essisto')
                continue
            if os.path.getsize(chainfile) == 0:
                continue
            print("loading {}".format(chainfile))
            chain_i = np.loadtxt(chainfile).squeeze()
            pars = np.loadtxt(outdir + '/pars.txt', dtype=np.unicode_)
            psr_pars.append(pars)
            
            #pp = model_utils.PostProcessing(chain_i, pars)
            #pp.plot_trace()
            #plt.savefig(chainfile.replace('.txt', '_{}_{}_trace.png'.format(psrname, dir)))
            #plt.close()
            print("Made {}".format(chainfile.replace('.txt', '_{}_{}_trace.png'.format(psrname, dir))))
            
            # Burn larger of first 25% or 2500
            burn = int(max([0.25*chain_i.shape[0], 1300]))
            chain = chain_i[burn:, :]
            if chain.size == 0:
                continue

            if first_chain:
                chain_total = chain
                first_chain = False
            else:
                chain_total = np.concatenate((chain_total, chain)).squeeze()
        psr_chains.append(chain_total)
        np.save('{}/{}_chain.npy'.format(datadir + "/chains/{}/".format(noise_model), psrname), chain_total)

psrnum = 0
MJD_JUMP_boundaries = np.arange(50000,60000,182.625)
i = 0
good_psrnum = 0

for psrname, pars, chain_total in zip(psrs, psr_pars, psr_chains):
    if chain_total is None:
        continue
    
    #JUMP_t0_pars = [p for p in pars if ('JUMP' in p ) and ('t0' in p)]
    #njumps = [int(p.replace(f'{psr}_MJD_JUMP_', '').replace('_t0', '')) for p in JUMP_t0_pars]

    

    if psrname in ['J1909-3744', 'J0437-4715', 'J1125-6014', 'J1713+0747', 'J1744-1134']:
        psrnames.append(psrname)
        color = colors[good_psrnum]
        linewidth = 2.5
        alpha = 0.8
        good_psrnum += 1
        label = psrname
        if '0437' in psrname:
            zorder += 500
    else:
        color = '0.8'
        alpha = 0.6
        linewidth = None
        label = None
        zorder = 1
    i+=1
    ind = np.argwhere(['JUMP_{}_t0'.format(njump) in p for p in pars]).squeeze()
    linestyle = '-'
    #print(psrname, chain_total.size)
    chain_corner = chain_total[:, ind]

    p, bins, patches = ax1.hist(chain_corner, bins=100, range=(MJD_JUMP_boundaries[njump], MJD_JUMP_boundaries[njump+1]), density=True, color=color, histtype='step', zorder=psrnum, linestyle = linestyle, linewidth = linewidth, alpha = alpha, label = label)
    ratio_ind = np.argmin(np.abs(bins - 58925))
    ratio = p[ratio_ind] / np.median(p[(bins < 58850)[:-1]])
    print(psrname, ratio)
    ind = np.argwhere(['JUMP_{}_log10_Amp'.format(njump) in p for p in pars]).squeeze()
    good_inds = chain_corner > 58850
    chain_corner2 = chain_total[:, ind].squeeze()#[good_inds]
    p2, bins2, patches2 = ax2.hist(chain_corner2, bins=100, range=(-9,-6), density=True, color=color, histtype='step', zorder=psrnum, linestyle = linestyle, linewidth = linewidth, alpha = alpha, label = label)
    

    indmax = np.argmax(p2)
    centres = bins2[0:-1] + np.diff(bins2)/2.0
    print('this is the jump amp:', centres[indmax])

    if first_chain_overall:
        p[np.isnan(p)] = 1e-30
        p_total = (p + 1e-20)
        p2[np.isnan(p2)] = 1e-30
        p_total2 = (p2 + 1e-20)
        first_chain_overall = False
        
        print('P TOTAL', len(p_total))
    else:
        p[np.isnan(p)] = 1e-30
        p_total *= (p + 1e-20)
        p2[np.isnan(p2)] = 1e-30
        p_total2 *= (p2 + 1e-20)
    psrnum += 1

#print(np.sum(np.isnan(p_total)), 'NANs in p_total')

#psrnames.append('Total')
#plt.xlabel(r'log$_{10}$ A$_{\rm GW}$')
plt.sca(ax1)
plt.stairs(p_total / np.nansum(p_total) / np.nanmean(np.diff(bins)), bins, color='k', zorder=0, linewidth=2, label = 'Total')
plt.xlabel(r'JUMP epoch', fontsize = 14)
plt.ylabel('Probability density', fontsize = 14)
plt.legend(loc='upper left', fontsize=12)
plt.xlim(MJD_JUMP_boundaries[njump], MJD_JUMP_boundaries[njump+1])#[-20, -12])
plt.sca(ax2)
plt.stairs(p_total2 / np.nansum(p_total2) / np.nanmean(np.diff(bins2)), bins2, color='k', zorder=0, linewidth=2, label = 'Total')
plt.xlabel(r'$\log_{{10}} A_\mathrm{{JUMP}}$', fontsize = 14)
plt.ylabel('Probability density', fontsize = 14)
#plt.legend(loc='upper left', fontsize=12)
plt.xlim(-9,-6)
#plt.xlim(MJD_JUMP_boundaries[njump], MJD_JUMP_boundaries[njump+1])#[-20, -12])
plt.tight_layout()

indmax = np.argmax(p_total2)
centres = bins2[0:-1] + np.diff(bins2)/2.0
print('this is the jump amp:', centres[indmax])

ind = np.argmax(p_total)
centres = bins[0:-1] + np.diff(bins)
print(centres[ind])

plt.savefig('{}_factorlike_JUMP_{}.pdf'.format(dir, njump))
plt.savefig('{}_factorlike_JUMP_{}.png'.format(dir, njump))
plt.close()    
