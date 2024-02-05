#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 16 21:01:00 2023

@author: dreardon
"""

import json
import numpy as np
import matplotlib.pyplot as plt
from correlations_utils_ipta import get_psrnames, get_pairs, hd_orf, dipole, \
    plot_violin

"""
Define Matplotlib settings
"""
from matplotlib import rc
import matplotlib
rc('text', usetex=True)
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
matplotlib.rcParams.update({'font.size': 14})

figsize=(2.0*3.3,1.5*3.3)

"""
Load pulsar names, pair names, positions and correlation data
"""
psrnames = get_psrnames()
pairs = get_pairs(psrnames)
with open('positions.json', 'r') as f:
    positions = json.load(f)
string_scram = 'scrambles_02_ppta'  # 10,000 skies with match threshold <= 0.2
scrambles = np.load('{}.npz'.format(string_scram))

# Choose the kernel function and bandwidth method:

bw = 0.274
# bw = 'silverman'
# kernel = 'epa'
kernel = 'gaussian'

corr_hd_all = {}  # fixed amplitude
corr_hd_all_reweight = {}  # reweighted amplitude
from scipy.interpolate import interp1d
mxcorr = 0.99
y_flat = np.linspace(-mxcorr, mxcorr, 256)
for i in range(1, 436): # 435 pulsar pairs
    pair = pairs[str(1 + (int(i) - 1) % 435)]
    pdf_1d_reweight = np.load("corr_chains/{}_corr_hd_reweight_more1713_{}_{}.npz".format(i, kernel, bw))
    pdf_1d = np.load("corr_chains/{}_corr_hd_fixA_{}_{}.npz".format(i, kernel, bw))
    x_data_reweight = pdf_1d_reweight["arr_1"].squeeze()
    y_data_reweight = pdf_1d_reweight["arr_0"].squeeze()
    x_data = pdf_1d["arr_1"].squeeze()
    y_data = pdf_1d["arr_0"].squeeze()

    f = interp1d(x_data, y_data, kind='linear', assume_sorted=True)
    corr_hd_all['_'.join(pair)] = f(y_flat)
    f = interp1d(x_data_reweight, y_data_reweight, kind='linear', assume_sorted=True)
    corr_hd_all_reweight['_'.join(pair)] = f(y_flat)


# Define sky separations
nseps = 8
vals = np.linspace(0, 180, nseps+1)
halfdiff = np.mean(np.diff(vals))/2
vals = (vals - halfdiff)[1:]
seps = {}
for ibin in range(1, nseps+1):
    seps["bin_{}".format(ibin)] = vals[ibin-1]

"""
Loop through Nscramble sky scrambles
"""
Nscramble = 1
plot = True
ts_scrambles = []
ts_scrambles_reweight = []

for ns in range(0, Nscramble):

    npos = len(psrnames)
    new_pos = {}
    new_pos_arr = []
    for nump, psr in enumerate(psrnames):
        # get new sky position
        theta = scrambles['thetas'][ns, nump]
        phi = scrambles['phis'][ns, nump]
        new_pos[psr] = np.array([np.cos(phi)*np.sin(theta),
                                 np.sin(phi)*np.sin(theta),
                                 np.cos(theta)])
        new_pos_arr.append([np.cos(phi)*np.sin(theta),
                            np.sin(phi)*np.sin(theta),
                            np.cos(theta)])

    # Use original data if not scrambling
    pos = positions if Nscramble == 1 else new_pos

    # Look through pairs and make plots
    kde = {}
    pdf = {}
    ptot = {}

    orf_bins_total = {}
    orf_bins_total_reweight = {}
    likelihood_hd = {}
    likelihood_hd_reweight = {}
    likelihood_curn = {}
    likelihood_mono = {}
    likelihood_curn_reweight = {}
    likelihood_hd_global = 0
    likelihood_hd_global_reweight = 0
    likelihood_mono_global = 0
    likelihood_dipole_global = 0
    likelihood_curn_global = 0
    likelihood_curn_global_reweight = 0
    likelihood_mono_global_reweight = 0
    likelihood_dipole_global_reweight = 0

    bf_hd = {}
    bf_hd_reweight = {}
    n_bins = {}
    angsep_array = {}
    orf_val_array = {}
    n_tot = 0

    # y_flat = np.linspace(-1, 1, 256)
    null_prob = np.ones(np.shape(y_flat))
    null_prob /= np.sum(null_prob)
    null_prob /= np.mean(np.diff(y_flat))

    """
    Loop through all pulsar pairs
    """

    for i in range(1, 436): # 435 pulsar pairs

        pair = pairs[str(1 + (int(i) - 1) % 435)]

        psr1 = pair[0]
        psr2 = pair[1]

        corr_hd = np.array(corr_hd_all['_'.join(pair)])
        corr_hd_reweight = np.array(corr_hd_all_reweight['_'.join(pair)])

        # normalise
        corr_hd /= np.sum(corr_hd)
        corr_hd /= np.mean(np.diff(y_flat))

        f = interp1d(y_flat, corr_hd, kind='linear')

        # normalise
        corr_hd_reweight /= np.sum(corr_hd_reweight)
        corr_hd_reweight /= np.mean(np.diff(y_flat))

        f_reweight = interp1d(y_flat, corr_hd_reweight, kind='linear')

        # calculate angular separation and ORF values
        pos1 = pos[psr1]
        pos2 = pos[psr2]
        angsep = np.arccos(np.dot(pos1, pos2)) * 180/np.pi
        orf_val = hd_orf(np.array([angsep]))[0]

        angsep_array['_'.join(pair)] = angsep
        orf_val_array['_'.join(pair)] = orf_val

        dipole_val = dipole(np.array([angsep]))[0]
        if dipole_val >= 0.99:
            dipole_val = 0.99
        if dipole_val <= -0.99:
            dipole_val = -0.99

        # Append to likelihoods for psr1
        try:
            likelihood_hd[psr1] += np.log(f(orf_val))
            likelihood_mono[psr1] += np.log(f(0.99))
            likelihood_curn[psr1] += np.log(f(0))
            likelihood_hd_reweight[psr1] += np.log(f_reweight(orf_val))
            likelihood_curn_reweight[psr1] += np.log(f_reweight(0))
        except KeyError:
            likelihood_hd[psr1] = np.log(f(orf_val))
            likelihood_mono[psr1] = np.log(f(0.99))
            likelihood_curn[psr1] = np.log(f(0))
            likelihood_hd_reweight[psr1] = np.log(f_reweight(orf_val))
            likelihood_curn_reweight[psr1] = np.log(f_reweight(0))

        # Append to likelihoods for psr2
        try:
            likelihood_hd[psr2] += np.log(f(orf_val))
            likelihood_mono[psr2] += np.log(f(0.99))
            likelihood_curn[psr2] += np.log(f(0))
            likelihood_hd_reweight[psr2] += np.log(f_reweight(orf_val))
            likelihood_curn_reweight[psr2] += np.log(f_reweight(0))
        except KeyError:
            likelihood_hd[psr2] = np.log(f(orf_val))
            likelihood_mono[psr2] = np.log(f(0.99))
            likelihood_curn[psr2] = np.log(f(0))
            likelihood_hd_reweight[psr2] = np.log(f_reweight(orf_val))
            likelihood_curn_reweight[psr2] = np.log(f_reweight(0))

        ibin = np.argmin(np.abs(vals - angsep)) + 1
        try:
            orf_bins_total["bin_{}".format(ibin)] *= corr_hd
            orf_bins_total_reweight["bin_{}".format(ibin)] *= corr_hd_reweight
            n_bins["bin_{}".format(ibin)] +=1
        except KeyError:
            orf_bins_total["bin_{}".format(ibin)] = corr_hd
            orf_bins_total_reweight["bin_{}".format(ibin)] = corr_hd_reweight
            n_bins["bin_{}".format(ibin)] = 1
        try:
            likelihood_hd_global += np.log(f(orf_val))
            likelihood_curn_global += np.log(f(0))
            likelihood_mono_global += np.log(f(0.99))
            likelihood_dipole_global += np.log(f(dipole_val))
            likelihood_hd_global_reweight += np.log(f_reweight(orf_val))
            likelihood_curn_global_reweight += np.log(f_reweight(0))
            likelihood_mono_global_reweight += np.log(f_reweight(0.99))
            likelihood_dipole_global_reweight += np.log(f_reweight(dipole_val))
        except KeyError:
            likelihood_hd_global = np.log(f(orf_val))
            likelihood_curn_global = np.log(f(0))
            likelihood_mono_global = np.log(f(0.99))
            likelihood_dipole_global = np.log(f(dipole_val))
            likelihood_hd_global_reweight = np.log(f_reweight(orf_val))
            likelihood_curn_global_reweight = np.log(f_reweight(0))
            likelihood_mono_global_reweight = np.log(f_reweight(0.99))
            likelihood_dipole_global_reweight = np.log(f_reweight(dipole_val))

    if plot:
        figsize2=(2.4*3.3,1.5*3.3)
        fig, ax1 = plt.subplots(1, 1, figsize=figsize2)
        ax2 = ax1.twinx()

        values = []
        for bini in range(1, nseps+1):
            try:
                values.append(n_bins['bin_{}'.format(bini)])
            except KeyError:
                n_bins['bin_{}'.format(bini)] = 0
                values.append(n_bins['bin_{}'.format(bini)])
        edges = np.append(np.array(list(seps.values()) - halfdiff), 180)

        ax2.fill_between(np.insert(np.array(list(seps.values()) + halfdiff),0, 0),
                         np.insert(np.array(values),0,values[0]),
                         y2=-np.ones(len(values)+1),
                         step="pre", alpha=0.2, color='k')

        for k in orf_bins_total.keys():
            orf_bins_total[k] /= np.sum(orf_bins_total[k])
            orf_bins_total[k] /= np.mean(np.diff(y_flat))
            draws = np.random.choice(y_flat, p=orf_bins_total_reweight[k]/np.sum(orf_bins_total_reweight[k]), size=100000)
            ax = plot_violin(ax1, seps[k], draws, width=np.floor(np.mean(np.diff(vals)))-2, colour='darkgreen', alpha=0.9, edgecolour='darkgreen', linewidth=0)

            draws = np.random.choice(y_flat, p=orf_bins_total[k]/np.sum(orf_bins_total[k]), size=100000)
            ax = plot_violin(ax1, seps[k], draws, width=np.floor(np.mean(np.diff(vals)))-2, alpha=1, colour='none', edgecolour='k', linewidth=2)


            ibin = int(k.split('_')[-1])

        theta = np.linspace(0, 180, 1000)
        orf = hd_orf(theta)
        ax1.plot(theta, orf, color='k', linewidth=3, linestyle='--')
        ax1.set_ylim([-1, 1])
        ax2.set_ylim([0, round(max(n_bins.values())+10, -1)])
        ax2.set_ylabel('Number of pulsar pairs')
        plt.xlim([0, 180])
        plt.tight_layout()
        ax1.set_zorder(ax2.get_zorder()+1)
        ax1.set_frame_on(False)
        plt.savefig('corr_plots/corr_total_{}_{}.png'.format(bw, kernel))
        plt.savefig('corr_plots/corr_total_{}_{}.pdf'.format(bw, kernel))
        plt.show()
        plt.close()

    ts = (likelihood_hd_global - likelihood_curn_global)
    ts_reweight = (likelihood_hd_global_reweight - likelihood_curn_global_reweight)

    print(" ")
    print('HD test statistic')
    print(ns, ts)
    print(ns, ts_reweight)

    if ns >= 1:
        ts_scrambles.append(ts)
        ts_scrambles_reweight.append(ts_reweight)

if ns >= 1:
    np.save("likelihood_ratios_{}_{}_{}.npy".format(bw, string_scram, kernel), ts_scrambles)
    np.save("likelihood_ratios_{}_{}_{}_reweight.npy".format(bw, string_scram, kernel), ts_scrambles_reweight)
ts_scrambles = np.load("likelihood_ratios_{}_{}_{}.npy".format(bw, string_scram, kernel))
ts_scrambles_reweight = np.load("likelihood_ratios_{}_{}_{}_reweight.npy".format(bw, string_scram, kernel))


plt.figure(figsize=figsize)
plt.hist(np.log10(np.exp(ts_scrambles_reweight)), bins=16, range=(-3.5, 3.5), alpha=0.9, color='darkgreen', log=False, density=False)
plt.hist(np.log10(np.exp(ts_scrambles)), bins=16, range=(-3.5, 3.5), alpha=1, facecolor='None', linestyle='-', lw=2, ec='k', log=False, density=False)
yl=plt.ylim()
plt.plot([np.log10(np.exp(ts)), np.log10(np.exp(ts))], yl, 'k--', linewidth=2)
plt.plot([np.log10(np.exp(ts_reweight)), np.log10(np.exp(ts_reweight))], yl, 'darkgreen', linestyle='--', linewidth=2)

plt.ylim(yl)
plt.ylabel(r'Number of randomized skies, $N_{\rm sky}$')
plt.xlabel(r'$\log_{10} \Delta\mathcal{L}^{\rm HD}_{\rm CRN}$')
plt.savefig('corr_plots/scrambles_{}_{}.png'.format(bw, kernel))
plt.savefig('corr_plots/scrambles_{}_{}.pdf'.format(bw, kernel))
plt.show()

print(len(np.array(ts_scrambles_reweight)[(ts_scrambles_reweight > ts_reweight)]) / len(np.array(ts_scrambles_reweight)) * 100)
print(len(np.array(ts_scrambles)[(ts_scrambles > ts)]) / len(np.array(ts_scrambles)) * 100)





