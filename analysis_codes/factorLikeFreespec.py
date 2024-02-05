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
from scipy.signal import correlate
from scipy.stats import gaussian_kde

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


def plot_violin_psr(fig, ax, psrname, chain, freqs, pars, outdir = 'noisefiles_freespec/', noise_model = '', fc = 'C0'):

    cpars = [f'{psrname}_red_noise_log10_rho_{i}' for i in range(60)]
    print(cpars)
    indices = get_par_indices(pars, cpars)
    print(indices)
    corner_pars = [p.replace('{}_'.format(psrname), '') for p in pars[indices]]
    chain = chain[:,:-4]
    chain = chain[:, indices]
    print(chain.shape)
    df = np.median(np.diff(freqs))
    #fig, ax = plt.subplots(1,1)
    
    for pl, orb_freq in planet_frequencies.items():
        ax.axvline(orb_freq, ls = '--', c = 'C0', alpha = 0.5)
        ax.text(1.0*orb_freq, 1E-5, f'{pl}', ha = 'center', va = 'center')
    
    ax.set_xlabel('Frequency (Hz)', fontsize = 18)
    ax.set_ylabel(r'$\rho$ (s)', fontsize =  18)
    ax.set_xscale('log')
    ax.set_yscale('log')
    
    violin_dict = ax.violinplot(10.0**(1.0*chain), positions = freqs, widths = df, showextrema = False, showmeans = False, showmedians = False)
    for violin_body in violin_dict['bodies']:
        violin_body.set_alpha(0.08)
        
        violin_body.set_facecolor("None")
        violin_body.set_edgecolor('k')
        #violin_body.set_facecolor(fc)
        violin_body.set_linewidth(0)

    return fig,ax,outdir


def power_spec(freqs, A, gamma, fc):

    psd = A**2.0/(12.0 * np.pi**2.0) * (freqs / fc ) ** -gamma * fc**-3.0

    return psd

def freespec_psd_var(var, tspan):
    return (np.sqrt(var / tspan ))

Tspan = (59645 - 53040)*86400
def rho_to_cp(var, tspan = Tspan):
    return (12.0 * np.pi**2.0 * tspan * (10.0**var)**2.0)**0.5

os.environ["PATH"] += os.pathsep + '/Library/TeX/texbin/'
rc('text', usetex=True)
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
#matplotlib.rcParams.update({'font.size': 26})


noise_model = 'singlePsrNoise_sw_nesw0_fixwhite_rn_freespec_fixwhite'
dir = sys.argv[1]
cpars = ['dm', 'red', 'gwb']
#if dir == 'all':
#    datadir = os.path.abspath("./data/{*")
#else:
datadir = os.path.abspath("./data/" + str(dir))

nchain = 8 #for PPTA15,  5 is 440 dr3 noise models, 6 is 440 basic noise, 7 is 421 basic noise, 8 is 438 basic noise

psrs = np.loadtxt('./psrs.dat', dtype = 'str')

#chainfiles = sorted(glob.glob(datadir + "/chains/{}/".format(noise_model) + '/J*_' + str(nchain) + '/chain_1.txt'))

number = len(psrs)
print(number)
cmap = plt.get_cmap('gist_rainbow')
colors = [cmap(i) for i in np.linspace(0, 1, number)]

first_chain_overall = True
plt.figure(figsize=(12, 12))
psrnames = []
psr_chains = []
psr_pars = []


components = 30
freqs = 1.0*np.arange(1.0, components + 1, 1.0)/Tspan
df = np.median(np.diff(freqs))
print(freqs)
print(freqs[7])

for psrnum, psrname in enumerate(psrs):
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
        max_nchain = 30  # number of chains to read
        
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
            burn = int(max([0.25*chain_i.shape[0], 20000]))
            print(chain_i.shape[0], burn)
            gw_log10A_ind = np.squeeze(np.argwhere(['rho_0' in p for p in pars]))
            thin = acor(chain_i[:, gw_log10A_ind])
            #thin = 1
            print('autocorrelation length = {}'.format(thin))
            chain = chain_i[burn::thin, :]
            print(chain.shape)
            
            #chain = chain_i[burn:, :]
            if chain.size == 0:
                continue

            if first_chain:
                chain_total = chain
                first_chain = False
            else:
                chain_total = np.concatenate((chain_total, chain)).squeeze()
            chain_i = None
            del(chain_i)
            chain = None
            del(chain)
        psr_chains.append(chain_total)
        np.save('{}/{}_chain.npy'.format(datadir + "/chains/{}/".format(noise_model), psrname), chain_total)
        chain_total = None
        del(chain_total)

psrnum = 0

#freespec_ps = np.zeros((len(psrs), components))


rho_min = -9
rho_max = -3

fig, ax = plt.subplots(1,1,figsize = (4.0*3.3, 2.0*3.3))



log10_rho_total = {}
log10_rho_psr = {}

log10_rho_total_array = np.zeros((256, components))
for psrname, pars, chain_total in zip(psrs, psr_pars, psr_chains):
    if psrname in ["J1824-2452A"]:
        continue
    bws = [0.01 for i in np.arange(0, components, 1)]
    bws[0] = 0.1
    bws[1] = 0.1
    bws[2] = 0.08
    bws[3] = 0.07
    for freq_ind in np.arange(0, components, 1):
        
        freespec_data_ind = np.argwhere(pars == f"{psrname}_red_noise_log10_rho_{freq_ind}")
        freespec_data = chain_total[:, freespec_data_ind]
        freespec_data_reflect_pos = freespec_data.copy()
        freespec_data_reflect_neg = freespec_data.copy()
        freespec_data_reflect_pos[:, 0] = 2*rho_max - freespec_data_reflect_pos[:, 0]
        freespec_data_reflect_neg[:, 0] = 2*rho_min - freespec_data_reflect_neg[:, 0]
        freespec_data_reflect = np.concatenate((freespec_data_reflect_neg, freespec_data, freespec_data_reflect_pos))
        freespec_data_reflect_pos = None
        freespec_data_reflect_neg = None
        del(freespec_data_reflect_pos)
        del(freespec_data_reflect_neg)
        rvs = freespec_data_reflect.squeeze()
        print(psrname, freq_ind)
        #print(rvs.shape)
        kde = gaussian_kde(rvs, bw_method = bws[freq_ind])
        #x_flat = np.array(freqs[freq_ind])
        y_flat = np.linspace(rho_min, rho_max, 256)

        #x, y = np.meshgrid(x_flat, y_flat)
        #grid_coords = np.append(x.reshape(-1, 1), y.reshape(-1,1),axis = 1)
        z = kde(y_flat)
        #z = z.reshape(256,1)
        z = z.squeeze()
        z /= (np.sum(z)*np.mean(np.diff(y_flat)))
        log10_rho_total_array[:, freq_ind] = z
        try: log10_rho_psr[psrname].append(z)
        except KeyError: log10_rho_psr[psrname] = [z]
        try:
            if psrname == "J1643-1224" and freq_ind == 26:
                log10_rho_total[f"rho_{freq_ind}"] *= np.ones_like(z)/(rho_max - rho_min)

            else:
                log10_rho_total[f"rho_{freq_ind}"] *= z
                #pass
        except KeyError: log10_rho_total[f"rho_{freq_ind}"] = z

        draws = np.random.choice(y_flat, p=z / np.sum(z), size=100000)
        #draws = rho_to_cp(draws)
        # violin_dict = ax.violinplot((freespec_data.squeeze()), positions = np.array([freqs[freq_ind] + df/2.0]), widths = df, showextrema = False, showmeans = False, showmedians = False, bw_method = 0.1)
        # for violin_body in violin_dict['bodies']:
        #     violin_body.set_alpha(0.02)
        
        #     violin_body.set_facecolor('None')
        #     #violin_body.set_edgecolor('0.5')
        #     violin_body.set_edgecolor('k')
        #     #violin_body.set_facecolor(fc)
        #     violin_body.set_linewidth(0.8)

for freq_ind in np.arange(0, components, 1):
    draws = np.random.choice(y_flat, p = log10_rho_total[f"rho_{freq_ind}"] / np.sum(log10_rho_total[f"rho_{freq_ind}"]), size = 100000)
    print(draws.shape)
    #draws = rho_to_cp(draws)
    
    violin_dict = ax.violinplot((1.0*draws), positions = [freqs[freq_ind] ], widths = df, showextrema = False, showmeans = False, showmedians = False)
    for violin_body in violin_dict['bodies']:
        violin_body.set_alpha(1)
        
        violin_body.set_facecolor('darkgreen')
        violin_body.set_edgecolor(None)
        #violin_body.set_facecolor(fc)
        violin_body.set_linewidth(0)
        
np.save('freespec_factorised_posterior.npy', log10_rho_total_array)
np.save('freespec_factorised_rhos.npy', y_flat)

        
freespec_pars = np.loadtxt('data/all/chains/commonNoise_freespec_hd_DE440/1/pars.txt', dtype = np.unicode_)
rho_inds = np.argwhere(['gw_freespec_hd_log10_rho' in f for f in freespec_pars])

freespec_chain = np.load('freespec_chain_hd.npy')
chain_rho = freespec_chain[:, rho_inds]
for freq_ind in np.arange(0, components, 1):

    draws = chain_rho[:, freq_ind]
    #draws = rho_to_cp(draws)
    violin_dict = ax.violinplot((1.0*draws), positions = [freqs[freq_ind]], widths = df, showextrema = False, showmeans = False, showmedians = False, bw_method=0.08)
    for violin_body in violin_dict['bodies']:
        violin_body.set_alpha(1.0)
        
        violin_body.set_facecolor('None')
        violin_body.set_edgecolor('goldenrod')
        #violin_body.set_facecolor(fc)
        violin_body.set_linewidth(0.8)


freespec_pars = np.loadtxt('data/all/chains/commonNoise_freespec_nocorr_basicprior_DE440/400/pars.txt', dtype = np.unicode_)
rho_inds = np.argwhere(['gw_freespec_nocorr_log10_rho' in f for f in freespec_pars])
freespec_chain = np.load('freespec_chain_bp.npy')
chain_rho = freespec_chain[:, rho_inds]
# for freq_ind in np.arange(0, components, 1):

#     draws = chain_rho[:, freq_ind]
#     violin_dict = ax.violinplot((1.0*draws), positions = [freqs[freq_ind]], widths = df, showextrema = False, showmeans = False, showmedians = False)
#     #violin_dict = ax_main.violinplot((1.0*draws), positions = np.log10([freqs[freq_ind]]), widths = np.log10(df)*np.log10(freqs[freq_ind]), showextrema = False, showmeans = False, showmedians = False)
#     for violin_body in violin_dict['bodies']:
#         violin_body.set_alpha(0.5)
        
#         violin_body.set_facecolor('None')
#         violin_body.set_edgecolor('mediumblue')
#         #violin_body.set_facecolor(fc)
#         violin_body.set_linewidth(0.8)

freespec_pars = np.loadtxt('data/all/chains/commonNoise_freespec_nocorr_DE440/1/pars.txt', dtype = np.unicode_)
rho_inds = np.argwhere(['gw_freespec_nocorr_log10_rho' in f for f in freespec_pars])
freespec_chain = np.load('freespec_chain.npy')
chain_rho = freespec_chain[:, rho_inds]
for freq_ind in np.arange(0, components, 1):

    draws = chain_rho[:, freq_ind]
    #draws = rho_to_cp(draws)
    violin_dict = ax.violinplot((1.0*draws), positions = [freqs[freq_ind]], widths = df, showextrema = False, showmeans = False, showmedians = False)
    #violin_dict = ax_main.violinplot((1.0*draws), positions = np.log10([freqs[freq_ind]]), widths = np.log10(df)*np.log10(freqs[freq_ind]), showextrema = False, showmeans = False, showmedians = False)
    for violin_body in violin_dict['bodies']:
        violin_body.set_alpha(1.0)
        
        violin_body.set_facecolor('None')
        violin_body.set_edgecolor('k')
        #violin_body.set_facecolor(fc)
        violin_body.set_linewidth(0.8)


freqs_plot = np.logspace(-9, -7, 10000)
psd = power_spec(freqs_plot, 10.0**-14.5, 3.87, (1.0/365.25/86400.0))
#print((1.0 / 365.25 / 24.0 / 3600.0))
print(psd[0])
rho = freespec_psd_var(psd, Tspan)
#rho += 2e-9
#ax.plot(freqs_plot, rho_to_cp(np.log10(rho)), c = 'k', lw = 2.0)
ax.set_xlabel('Frequency (Hz)', fontsize = 18)

#ax.set_ylabel(r'$\rho$ (s)')
#ax.set_ylabel(r'$\log_{{10}}(\rho / \mathrm{{s}})$', fontsize = 18)
#ax.set_ylabel(r'
ax.set_xscale('log')
#ax.set_yscale('log')
#ax.set_xlim(9e-10, 5.5e-8)
ax.set_xlim(freqs[0]-df/2.0, freqs[components-1])

ax.set_ylim(-9, -5.5)
ax.grid(which='both', ls=':', alpha=0.5, zorder=1)
fig.savefig('test2.png', dpi = 300, bbox_inches = 'tight')
fig.savefig('common_freespec_factorised.png', dpi = 300, bbox_inches = 'tight')
fig.savefig('common_freespec_factorised.pdf', bbox_inches = 'tight')


# MJD_JUMP_boundaries = np.arange(50000,60000,182.625)
# for psrname, pars, chain_total in zip(psrs, psr_pars, psr_chains):
#     if chain_total is None:
#         continue
    
#     #JUMP_t0_pars = [p for p in pars if ('JUMP' in p ) and ('t0' in p)]
#     #njumps = [int(p.replace(f'{psr}_MJD_JUMP_', '').replace('_t0', '')) for p in JUMP_t0_pars]
#     psrnames.append(psrname)
#     color = colors[psrnum]
    
#     ind = np.argwhere(['JUMP_{}_t0'.format(njump) in p for p in pars]).squeeze()
#     print(psrname, chain_total.size)
#     chain_corner = chain_total[:, ind]
#     p, bins, patches = plt.hist(chain_corner, bins=100, range=(MJD_JUMP_boundaries[njump], MJD_JUMP_boundaries[njump+1]), density=True, alpha=0.6, color=color, histtype='step', zorder=psrnum)

#     ind = np.argmax(p)
#     centres = bins[0:-1] + np.diff(bins)
#     print(centres[ind])

#     if first_chain_overall:
#         p[np.isnan(p)] = 1e-30
#         p_total = (p + 1e-20)
#         first_chain_overall = False
#     else:
#         p[np.isnan(p)] = 1e-30
#         p_total *= (p + 1e-20)
#     psrnum += 1

# print(np.sum(np.isnan(p_total)), 'NANs in p_total')
# plt.stairs(p_total / np.nansum(p_total) / np.nanmean(np.diff(bins)), bins, color='k', zorder=0, linewidth=2)
# psrnames.append('Total')
# #plt.xlabel(r'log$_{10}$ A$_{\rm GW}$')
# plt.xlabel(r'JUMP epoch')
# plt.ylabel('Probability density')
# plt.legend(psrnames, loc='upper left', fontsize=13)
# plt.xlim(MJD_JUMP_boundaries[njump], MJD_JUMP_boundaries[njump+1])#[-20, -12])
# plt.tight_layout()

# ind = np.argmax(p_total)
# centres = bins[0:-1] + np.diff(bins)
# print(centres[ind])

# plt.savefig('{}_factorlike_JUMP_{}.pdf'.format(dir, njump))
# plt.savefig('{}_factorlike_JUMP_{}.png'.format(dir, njump))
# plt.close()    
