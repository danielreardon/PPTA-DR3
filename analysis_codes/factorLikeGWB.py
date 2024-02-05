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

from scipy.signal import correlate

def acor(arr):
    arr -= np.mean(arr)
    auto_correlation = correlate(arr, arr, mode='full')
    auto_correlation = auto_correlation[auto_correlation.size//2:]
    auto_correlation /= auto_correlation[0]
    indices = np.where(auto_correlation<0.5)[0]
    if len(indices)>0:
        if indices[0] == 0:
            indices[0] = 1
        return indices[0]
    else:
        return 1





noise_model = 'singlePsrNoise_sw_nesw0_gwb_fixwhite_fixwhite'
dir = sys.argv[1]
cpars = ['dm', 'red', 'gwb']
#if dir == 'all':
#    datadir = os.path.abspath("./data/[up]*")
#else:
datadir = os.path.abspath("./data/" + str(dir))

nchain = 8 #for PPTA15,  5 is 440 dr3 noise models, 6 is 440 basic noise, 7 is 421 basic noise, 8 is 438 basic noise

psrs = np.loadtxt('./psrs.dat', dtype = 'str')

#chainfiles = sorted(glob.glob(datadir + "/chains/{}/".format(noise_model) + '/J*_' + str(nchain) + '/chain_1.txt'))

number = len(psrs)
number = 6
print(number)
cmap = plt.get_cmap('tab20')#gist_rainbow')
colors = [cmap(i) for i in np.linspace(0, 0.8, number)]

first_chain_overall = True
plt.figure(figsize=(2.0*3.3, 1.5*3.3))
psrnames = []
psr_chains = []
psr_ps = []
psrnum_ = 0
colors = {i: f'C{i}' for i in range(number)}
colors[3] = 'goldenrod' #manually over-ride red for colorblind-friendliness

for psrnum, psrname in enumerate(psrs):
    max_nchain = 30  # number of chains to read
    first_chain = True
    chain_total = None
    if psrname in ["J1741+1351", "J1824-2452A"]:
        continue
    pname = datadir + "/chains/{}/".format(noise_model)+f'/{psrname}_log10A_GW_p.npy'
    print(pname)
    
    if os.path.exists(pname):#:+'dfadsfafa'):
        #pass
        max_nchain = 0
        p = np.load(pname, allow_pickle = True)
        if first_chain_overall:
            p_total = (p + 1e-20)
            first_chain_overall = False
        else:
            p_total *= (p + 1e-20)
        bins = np.linspace(-20, -11, 111)
        centres = bins[0:-1] + np.diff(bins)/2.0

        psrnames.append(psrname)
        if psrname in ['J1909-3744', 'J0437-4715', 'J2145-0750']:#, 'J1939+2134']:
            
            color = colors[psrnum_]
            label = psrname
            #lw = 2.0
            lw = 1.8
            psrnum_ += 1
            alpha = 1.0
            zorder = psrnum + 8000
        elif psrname in ['J1713+0747', 'J1744-1134', 'J1603-7202']:
            color = colors[psrnum_]
            label = psrname
            lw = 1.8
            psrnum_ += 1
            alpha = 1.0
            zorder = psrnum  + 7990
        else:
            color = '0.8'
            label = None
            lw = 0.8
            alpha = 0.6
            zorder = 1

        
        plt.stairs(p / np.sum(p) / np.mean(np.diff(bins)), bins, alpha=alpha, color=color, zorder=zorder, lw = lw, label = label)
        #chain_total = np.load(NAME)
        #max_nchain = 1

    for i in range(1, max_nchain+1): 
        outdir = datadir + "/chains/{}/".format(noise_model) + psrname + '_' + str(i)
        chainfile = outdir + '/chain_1.txt'
        print(chainfile)
        nlines = int(subprocess.check_output(f'cat {chainfile} | wc -l', shell = True, text = True))
        print(nlines)
        if nlines < 10000:
            print(nlines, "dead")
            continue
        if not os.path.exists(chainfile):
            continue
        if os.path.getsize(chainfile) == 0:
            continue
        print("loading {}".format(chainfile))
        chain_i = np.loadtxt(chainfile).squeeze()
        pars = np.loadtxt(outdir + '/pars.txt', dtype=np.unicode_)
        ind = np.argwhere(['gw' in p for p in pars]).squeeze()
        #pp = model_utils.PostProcessing(chain_i, pars)
        #pp.plot_trace()
        #plt.savefig(chainfile.replace('.txt', '_{}_{}_trace.png'.format(psrname, dir)))
        #plt.close()
        print("Made {}".format(chainfile.replace('.txt', '_{}_{}_trace.png'.format(psrname, dir))))
        
        # Burn larger of first 25% or 2500
        burn = int(max([0.25*chain_i.shape[0], 10000]))
        #chain = chain_i[burn:, :]
        if chain_i.size == 0:
            continue
        thin = acor(chain_i[:, ind])
        print('autocorrelation length = {}'.format(thin))
        chain = chain_i[burn::thin, :]
        if first_chain:
            chain_total = chain
            first_chain = False
        else:
            chain_total = np.concatenate((chain_total, chain)).squeeze()
        #psr_chains.append(chain_total)
        chain_i = None
        del(chain_i)

    if chain_total is None:
        continue

    chain = None

    del(chain)

    psrnames.append(psrname)
    if psrname in ['J1909-3744', 'J0437-4715', 'J2145-0750']:#, 'J1939+2134']:
        
        color = colors[psrnum_]
        label = psrname
        #lw = 2.0
        lw = 1.2
        psrnum_ += 1
    elif psrname in ['J1713+0747', 'J1744-1134', 'J1603-7202']:
        color = colors[psrnum_]
        label = psrname
        lw = 2.0
        psrnum_ += 1
    else:
        color = '0.8'
        label = None
        lw = 0.8
    

    #if you didn't load in the chain:
    if max_nchain > 1 or len(chain_total.shape) > 1:
        ind = np.argwhere(['gw' in p for p in pars]).squeeze()
        chain_corner = chain_total[:, ind]
    else:
        #the chain is actually chain_corner
        chain_corner = chain_total
    # if "1939" not in psrname:
    #     ind = np.argwhere(['gw' in p for p in pars]).squeeze()
    #     chain_corner = chain_total[:, ind]
    # else:
    #     ind = np.argwhere(['red_noise_log10_A' in p for p in pars]).squeeze()
    #     ind_gamma = np.argwhere(['red_noise_gamma' in p for p in pars]).squeeze()
    #     rowinds = np.argwhere(np.abs(chain_total[:, ind_gamma] - 13.0 / 3.0) < 0.5).squeeze()
    #     chain_corner = chain_total[:, ind]
    #     chain_corner = chain_corner[rowinds, None]

    #if "1713" in psrname:
    #    lw = 3.5
    #    color = 'r'
    # if "1741" in psrname:
    #     lw = 3.5
    #     color = 'r'
    # elif 'J1824' in psrname:
    #     lw = 3.5
    #     color = 'b'
    # elif 'J1939' in psrname:
    #     lw = 3.5
    #     color = 'g'
    #else:

    color = color

    print(chain_corner.shape)
    p, bins, patches = plt.hist(chain_corner, bins=110, range=(-20, -11), density=True, alpha=0.6, color=color, histtype='step', zorder=psrnum, lw = lw, label = label)
    pname = datadir + "/chains/{}/".format(noise_model)+f'/{psrname}_log10A_GW_p.npy'
    if True:#not os.path.exists(pname):
        np.save(pname, p)
    psr_ps.append(p)
    ind = np.argmax(p)
    centres = bins[0:-1] + np.diff(bins)/2.0
    print(centres[ind])

    if first_chain_overall:
        p_total = (p + 1e-20)
        first_chain_overall = False
    else:
        p_total *= (p + 1e-20)

    #if not(os.path.exists(pname)):
        
    #    np.save(NAME, chain_corner)
    chain_corner = None
    del(chain_corner)
    chain_total = None
    del(chain_total)




prior_density = 1.0 / np.abs(-11- (-20))
avg_density = np.mean((p_total / np.sum(p_total) / np.mean(np.diff(bins)))[centres < -16.5])
sd_bf = prior_density / avg_density
print(np.log(sd_bf))
print(np.log10(sd_bf))
print(sd_bf)
print('asadfasfafad')
print(avg_density)
#sys.exit()
plt.axhline(prior_density, ls = '--', color = 'darkgreen', label = "Prior", zorder = 9999)
plt.stairs(p_total / np.sum(p_total) / np.mean(np.diff(bins)), bins, color='k', zorder=10000, linewidth=2, label = 'Total')
psrnames.append('Total')
plt.xlabel(r'$\log_{{10}} A_{{\mathrm{{GW}}}}$', fontsize = 14)
plt.ylabel('Probability density', fontsize = 14)
plt.legend(fontsize = 12, loc = 'center left')#loc='upper left', fontsize=13)
plt.xlim([-20, -11])
plt.ylim([10**-9, 10.0])

plt.yscale('log')
#plt.ylim(1e-10, 20)
plt.tight_layout()




ind = np.argmax(p_total)
centres = bins[0:-1] + np.diff(bins)/2.0
overall_best=centres[ind]
print(len(centres))
print(len(p_total))
print(centres[ind])
print('{}_factorlike_GWB.pdf'.format(dir))
plt.savefig('{}_factorlike_GWB.pdf'.format(dir), bbox_inches = 'tight')
plt.savefig('{}_factorlike_GWB.png'.format(dir), bbox_inches = 'tight', dpi = 300)
plt.close()    

sys.exit()
psr_ratios = np.array([p[ind] / np.median(p[centres <-16]) for p in psr_ps])



order = np.argsort(psr_ratios)
fig, ax = plt.figure(figsize = (4.0*3.3, 1.5*3.3))
ax.plot(np.array(psr_ratios)[order], 'o')
ax.set_xticks(np.arange(0, len(psr_ratios), 1))
ax.set_xticklabels(np.array(psrnames)[order], rotation = 45.0, ha = 'right')
ax.set_xlim(-1, len(psr_ratios))
ax.set_xlabel('Pulsar')
ax.set_ylabel(r'$p(\log_{{10}} A_\mathrm{{GW}} = -14.6) / \langle p(\log_{{10}} A_\mathrm{{GW}} < -16) \rangle$')
ax.set_yscale('log')
ax.grid(linestyle = ':')

