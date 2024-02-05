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
from enterprise_extensions import model_utils
import matplotlib.pyplot as plt
import corner
from matplotlib import rc
import matplotlib
from astropy.time import Time
from scipy.stats import gaussian_kde
os.environ["PATH"] += os.pathsep + '/Library/TeX/texbin/'
rc('text', usetex=True)
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
#matplotlib.rcParams.update({'font.size': 12})

def make_noise_files(psrname, chain, pars, outdir='noiseFiles/'):
    x = {}
    for ct, par in enumerate(pars):
        x[par] = np.median(chain[:, ct])

    os.system('mkdir -p {}'.format(outdir))
    with open(outdir + '/{}_noise.json'.format(psrname), 'w') as fout:
        json.dump(x, fout, sort_keys=True, indent=4, separators=(',', ': '))


def make_noise_files_maxlike(psrname, chain, pars, outdir='noiseFiles_maxlike/'):
    x = {}
    for ct, par in enumerate(pars):
        ind = np.argmax(chain[:, -4])  # find index of maximum posterior
        x[par] = chain[ind, ct] 

    os.system('mkdir -p {}'.format(outdir))
    with open(outdir + '/{}_noise.json'.format(psrname), 'w') as fout:
        json.dump(x, fout, sort_keys=True, indent=4, separators=(',', ': '))


def plot_corner(chain, pars, cpars, noise_model='', outdir='noiseFiles/'):

    chain = chain[:, :-4]
    if cpars is not None:

        all_indices = []

        for cp in cpars:
            indices = [cp in p for p in pars]
            all_indices.append(indices)
        indices = np.any(all_indices, axis=0)

        corner_pars = pars[indices]
        corner_file_label = noise_model+'_' + '_'.join(cpars)
    else:
        # plot all pars
        indices = np.array([True for p in pars])
        corner_file_label = noise_model

    chain_corner = chain[:, indices]
    fig = corner.corner(chain_corner, bins=30, labels=corner_pars,
                        quantiles=(0.16, 0.84), show_titles=True)
    for ax in fig.axes:
        xlab = ax.get_xlabel()
        ylab = ax.get_ylabel()
        ti = ax.get_title()
        ax.set_title(ti, fontsize=9)
        ax.set_xlabel(xlab, fontsize=9)
        ax.set_ylabel(ylab, fontsize=9)
    os.system('mkdir -p {}'.format(outdir))
    figsavename = outdir + '/corner' + corner_file_label + '.png'
    print(figsavename)
    plt.savefig(figsavename, dpi=300, bbox_inches='tight')

#dir = sys.argv[1]
dir = 'all'


sliceyrs = [6,9]
uplims_6yr = [True, True, True, True, True, True, False, False, False, False, False, False, False]
uplims_9yr = [True, True, True, False, False, True, False, False, False, False]

uplims_slice = {6: uplims_6yr, 9: uplims_9yr}

fig, (ax1, ax2) = plt.subplots(2,1,figsize = (2.0*3.3, 3.0*3.3), sharex = True, sharey = True)
for sliceyr, ax in zip(sliceyrs, [ax1,ax2]):
    
        noise_models = [i.replace('data/all/chains/slice_{}yr/'.format(sliceyr), '').replace('/', '') for i in sorted(glob.glob('data/all/chains/slice_{}yr/commonNoise_pl_nocorr_freegam_5????-5????_DE440'.format(sliceyr)))]
        #print(noise_models)
        #print("HEHFEDHAVFSDFIASDNVSADFA")
        log10_A_list = []
        z_list = []
        #noise_model = sys.argv[2]#'commonNoise_pl_nocorr_fixgam_freespec_monopole'#commonNoise_pl_nocorr_fixgam_pl_hdnoauto_fixgam'
        #cpars = None
        cpars = ['gw']
        if len(sys.argv) > 3:
            if ',' in sys.argv[3]:
                cpars = [p.strip() for p in sys.argv[3].split(',')]
            else:
                cpars = [p.strip() for p in sys.argv[3].split(' ')]
        datadir = os.path.abspath("./data/" + str(dir))
        max_nchain = 100  # number of chains to read
        slice_ind = 0
        for noise_model in noise_models:
            #print(noise_model)
            first_chain = True
            chain_total = None
            log10_A = None
            chain = None
            chain_i = None
            #print(f'data/all/chains/slice_{sliceyr}yr/{noise_model}_gw_log10_A.npy')
            if os.path.exists(f'data/all/chains/slice_{sliceyr}yr/{noise_model}_gw_log10_A.npy'):
                #print('found existing processed chain')
                #print(f'data/all/chains/slice_{sliceyr}yr/{noise_model}_gw_log10_A.npy')
                
                log10_A = np.load(f'data/all/chains/slice_{sliceyr}yr/{noise_model}_gw_log10_A.npy', allow_pickle = True)
                
                #print("THIS IS LIOGA [0]",  log10_A[0])
                #print(log10_A.shape)
                #print(log10_A[0])
                max_nchain = 1
            iii=0
            for i in range(1, max_nchain+1):
                pass
            #     if log10_A is not None:
            #         break
            #     if os.path.exists(f'data/all/chains/slice_{}yr{noise_model}_gw_log10_A.npy'.format(sliceyr)):
            #         #pass
            #         break
            #     outdir = datadir + "/chains/slice_{}yr/".format(sliceyr) + noise_model + '/' + str(i) 
            #     chainfile = outdir + '/chain_1.txt'
            #     if not os.path.exists(chainfile):
            #         continue
            #     if os.path.getsize(chainfile) == 0:
            #         continue

            #     print('loading {}...'.format(chainfile))
            #     try:

            #         chain_i = np.loadtxt(chainfile).squeeze()
            #         print('chain loaded')
            #     except ValueError as e:
            #         print(e)
            #         continue
            #     pars = np.loadtxt(outdir + '/pars.txt', dtype=np.unicode_)

            #     #pp = model_utils.PostProcessing(chain_i, pars)
            #     #pp.plot_trace()
            #     #plt.savefig(chainfile.replace('.txt', '.png'))
            #     #plt.close()
            #     print('saved {}'.format(chainfile.replace('.txt', '.png')))

            #     # Burn larger of 25% or first 25000
            # #    burn = int(max([0.25*chain_i.shape[0], 25000]))
            #     nburnmin = 15000
            #     burn = int(max([0.25*chain_i.shape[0], nburnmin]))
            #     print(burn)

            #     gw_log10A_ind = np.squeeze(np.argwhere(np.logical_and(['gw' in p for p in pars], ['log10_A' in p for p in pars])))
            #     gw_gamma_ind = np.squeeze(np.argwhere(np.logical_and(['gw' in p for p in pars], ['gamma' in p for p in pars])))
            #     print('slicing chain')
            #     chain = chain_i[burn:, :]
            #     smbhb_inds =  np.abs(chain[:, gw_gamma_ind] - 4.3333) < 0.1
            #     chain = chain[smbhb_inds]
            #     chain = chain[:, gw_log10A_ind]
            #     print('done slicing')
            #     if chain.size == 0:
            #         print('fewer than {} samples. Continuing'.format(nburnmin))
            #         continue
            #     print(chain.shape)

            #     #if cpars is not None:
            #     #    plot_corner(chain, pars, cpars, outdir=outdir)
            #         #print('Made corner plot for {} using chain {}'.format(psrname, str(i)))
            #     iii += 1

            #     if first_chain:
            #         chain_total = chain
            #         first_chain = False
            #     else:
            #         #if chain.shape[1] != chain_total.shape[1]:
            #         #    continue
            #         print('concatenating and deleting chains')
            #         chain_total = np.concatenate((chain_total, chain)).squeeze()
            #         print(chain_total.shape)
            #         chain_i = None
            #         del(chain_i)
            #         chain = None
            #         del(chain)
            #         print('done deleting')
            #     if iii > 100: break
            logA_max = -11
            logA_min = -18
            y_flat = np.linspace(logA_min, logA_max, 256)
            log10_A_reflect_pos = log10_A.copy()
            log10_A_reflect_neg = log10_A.copy()
            log10_A_reflect_pos = 2*logA_max - log10_A_reflect_pos
            log10_A_reflect_neg = 2*logA_min - log10_A_reflect_neg
            log10_A_reflect = np.concatenate((log10_A_reflect_neg, log10_A, log10_A_reflect_pos))
            log10_A_reflect_pos = None
            log10_A_reflect_neg = None
            del(log10_A_reflect_pos)
            del(log10_A_reflect_neg)
            rvs = log10_A_reflect.squeeze()
            #print(rvs.shape)
            kde = gaussian_kde(rvs, bw_method = 0.01)
            
            y_flat = np.linspace(logA_min, logA_max, 256)
            
            #x, y = np.meshgrid(x_flat, y_flat)
            #grid_coords = np.append(x.reshape(-1, 1), y.reshape(-1,1),axis = 1)
            z = kde(y_flat)
            #z = z.reshape(256,1)
            z = z.squeeze()
            z /= (np.sum(z)*np.mean(np.diff(y_flat)))
            #print(z.shape)

            if uplims_slice[sliceyr][slice_ind]:
                z *= (10.0**(y_flat) / np.abs(10.0**(-11.0) - 10.0**(-18.0)))
                #print('LOG A UNIFORM: {:.3f}'.format(y_flat[np.argmax(z)]))
            # else:
            #     print('LOG A LOGUNIFORM: {:.3f}'.format(y_flat[np.argmax(z)]))
            try:
            
                if not os.path.exists(f'data/all/chains/slice_{sliceyr}yr/{noise_model}_gw_log10_A.npy'):
                    log10_A = chain_total#[gw_log10A_ind]
                    np.save(f'data/all/chains/slice_{sliceyr}yr/{noise_model}_gw_log10_A.npy', log10_A)
                z_list.append(z)

                log10_A_list.append(log10_A)
                
            except NameError as e:
                print(e)
                z_list.append(z)
                log10_A_list.append(np.zeros(1000) - 17.5)
            slice_ind +=1
            #print(iii)
            #print(np.median(log10_A_list[0]))
            #if chain_total is None or log10_A is None:
            #    log10_A_list.append(np.zeros(1000) - 17.5)
            # max_nchain = 100


        yearslice = sliceyr
        times = []
        for islice in range(1, len(log10_A_list) + 1):

            min_mjd = 53040 + (islice-1)*365.2425
            max_mjd = min_mjd + yearslice*365.2425
            midslice = (max_mjd + min_mjd)/ 2.0
            times.append(midslice)
        #print("HFEHRHERRRRRRRRRRR")
        #for l in log10_A_list:
            #print(l)
            
        #p975s =
        p95s = []
        p05s = []
        p025s = []
        p975s = []
        p_commons = []
        print(sliceyr)
        for _ind_ in range(len(z_list)):
            z = z_list[_ind_]
            log10_A = log10_A_list[_ind_]
            draws = np.random.choice(y_flat, p=z / np.sum(z), size=100000)
            t = Time(times[_ind_], format = 'mjd').decimalyear
            p95 = np.percentile(draws, 95.0)
            p975 = np.percentile(draws, 97.5)
            p05 = np.percentile(draws, 5.0 )
            p025 = np.percentile(draws, 2.5 )
            p16 = np.percentile(draws, 15)
            p84 = np.percentile(draws, 84)
            bf = (1.0 / np.abs(-11 - -18)) / np.nanmedian(z[y_flat <= -16.5])
            print(times[_ind_], p95, uplims_slice[sliceyr][_ind_], bf, 10.0**np.median(draws), 10.0**p16, 10.0**p84, np.sum(z), np.sum(z[y_flat <= -14.7])/np.sum(z),np.sum(z[y_flat > -14.7])/np.sum(z))
            p95s.append( p95)
            p975s.append(p975)
            p05s.append( p05)
            p025s.append(p025)
            p_commons.append(np.sum(z[y_flat < -14.7]))


            

            #print('Time\t0.05\t0.5\t0.95')

            if uplims_slice[sliceyr][_ind_]:
                color_ = 'darkorange'
            else:
                color_ = 'mediumblue'

            #print('{:.2f}\t{:.4f}\t{:.4f}\t{:.4f}'.format(t, p05, med_, p95))
            violin_dict = ax.violinplot(draws, positions = [Time(times[_ind_], format = 'mjd').decimalyear], widths = yearslice*1 * 0.5, showextrema = False, showmeans = False, showmedians = False)#, quantiles = [[0.05, 0.5, 0.95] for l in log10_A_list])
            for violin_body, _ in zip(violin_dict['bodies'], log10_A_list):
                #violin_body.set_alpha(0.3 + 0.7*(1.0 - (np.percentile(nearth_, 99.0) - np.percentile(nearth_, 1.0))/20.0))
                violin_body.set_facecolor(color_)
                violin_body.set_edgecolor(None)
                violin_body.set_linewidth(0)
                
            if uplims_slice[sliceyr][_ind_]:
                ax.annotate("", xy=(t, -18.02), xycoords = 'data', xytext=(t, p95), textcoords='data', arrowprops=dict(arrowstyle="->", color = 'k', alpha = 0.4, linewidth=1.4))
                #ax.arrow(t, p95, 0, -np.abs(p95 - -18), color = 'k', alpha = 0.4, width = 0.1, head_width = 100, head_length =0.08, length_includes_head=True)
                #ax.plot([t,t], [-18, p95], 'k-', alpha = 0.4)
            else:
                ax.plot([t,t], [p025, p975], 'k-', lw = 1.4, alpha = 0.4)
            

        
        #widths = np.abs(elats)/2.0
        #print(len(log10_A_list), len(times))


        #for t,p05,p95 in zip(times, p05s, p95s):
            
            #ax.plot([t,t], [p05,p95], 'k-', alpha=0.4)
        ax.set_ylabel(r'$\log_{{10}} A^{{\mathrm{{CRN}}}}_{{13/3}}$', fontsize = 14)
        ax.grid(linestyle = ':')
        ax.axhline(-15.0, ls = '--', c = 'k', alpha = 0.6, lw = 1.1)

        xlims = ax.get_xlim()
        
        ax.fill_between((2004,2023), -14.69 - 0.05, -14.69 + 0.05, color = 'k', alpha = 0.2)
        ax.axhline(-14.69, ls = '-', c = 'k', alpha = 0.8, lw = 2.0)
        ax.text(2017.0,-14,f"{sliceyr} yr slices", ha = 'left', va='center', fontsize = 12)
        if ax == ax2:
            ax.set_xlabel('Year', fontsize = 14)
        ax.set_xlim(2005,2022)

            
plt.subplots_adjust(hspace = 0.05)
fig.savefig('log10_A_time_uniformul2.png', dpi = 300, bbox_inches = 'tight')
fig.savefig('log10_A_time_uniformul2.pdf', bbox_inches = 'tight')

    
#plt.errorbar(times, [np.nanmedian(i) for i in log10_A_list], [np.nanstd(i) for i in log10_A_list])
#plt.savefig('tmp.png', dpi = 300)
    
# # Now, save noise files
# try:
#     #make_noise_files(psrname, chain_total, pars, outdir=datadir+'/noiseFiles/')
#     #make_noise_files_maxlike(psrname, chain_total, pars, outdir=datadir+'/noiseFiles_maxlike/')
#     print(chain_total.shape)
#     #print('Wrote noise files for {}'.format(psrname))
#     if cpars is not None:
#         plot_corner(chain_total, pars, cpars, noise_model=noise_model, outdir=datadir+'/corners/')
#         print('Made corner plot')
# except Exception as e:
#     print(e)
    # print('No valid chains for {}'.format(psrname))
