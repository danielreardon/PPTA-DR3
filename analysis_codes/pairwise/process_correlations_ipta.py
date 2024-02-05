#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 16 21:01:00 2023

@author: dreardon
"""

import numpy as np
from KDEpy import FFTKDE

numerator = 0
numerator_68 = 0
denominator = 0
denominator_68 = 0
for i in range(1, 436): # 435 pulsar pairs

    # Find index corresponding to the correlation corefficient
    ind = -6  # corr coeff
    ind_amp = -7  # amplitude

    chain_name = "corr_chains/{}_fixA.npy".format(i)
    chain = np.load(chain_name)

    if "fixA" in chain_name:
        chain_1d = chain[:, ind].squeeze()

        s = np.std(chain_1d)
        s68 = (np.percentile(chain_1d, q=84) - np.percentile(chain_1d, q=16)) / 2
        numerator += 1/(s**2)
        numerator_68 += 1/(s68**2)
        denominator += 1/(s**4)
        denominator_68 += 1/(s68**4)

        chain_1d_reflect_corrs_pos = chain_1d.copy()
        chain_1d_reflect_corrs_pos = 2 - chain_1d_reflect_corrs_pos
        chain_1d_reflect_corrs_neg = chain_1d.copy()
        chain_1d_reflect_corrs_neg = -2 - chain_1d_reflect_corrs_neg

        chain_mirror = np.concatenate((chain_1d_reflect_corrs_neg, chain_1d, chain_1d_reflect_corrs_pos))

        # Choose bandwidth method
        for bw in ['silverman', 0.274]:

            # Choose kernel
            for kernel in ['gaussian','epa']:

                data = chain_mirror.squeeze()
                grid_points = 4096  # Grid points in each dimension

                kde = FFTKDE(kernel=kernel, bw=bw)
                grid, points = kde.fit(data).evaluate(grid_points)
                # The grid is of shape (obs, dims), points are of shape (obs, 1)
                y = np.unique(grid)
                # y is corr coeff
                z = points.reshape(grid_points)

                z[y<=-1] = 0  # Set the KDE to zero outside of the domain
                z[y>=1] = 0  # Set the KDE to zero outside of the domain
                z = z.squeeze() * 3  # multiply the kde to get integral of 1

                indy = np.argwhere((y>=-1)*(y<=1)).squeeze()

                y = y[indy]
                z = z[indy]

                np.savez("corr_chains/{}_corr_hd_fixA_{}_{}.npz".format(i, kernel, bw), z, y)
                print("Saved:", "corr_chains/{}_corr_hd_fixA_{}_{}.npz".format(i, kernel, bw))
                print("")

    else:
        chain_2d = chain[:, [ind_amp, ind]].squeeze()

        # estimate number of samples of interest, for bandwidth computation
        nsamps = len(np.argwhere( (chain_2d[:, 0] <= -14.69 + 0.05) *  (chain_2d[:, 0] >= -14.69 - 0.05))) / 0.68

        chain_2d_reflect_corrs_pos = chain_2d.copy()
        chain_2d_reflect_corrs_pos[:, 1] = 2*(1) - chain_2d_reflect_corrs_pos[:, 1]
        chain_2d_reflect_corrs_neg = chain_2d.copy()
        chain_2d_reflect_corrs_neg[:, 1] = 2*(-1) - chain_2d_reflect_corrs_neg[:, 1]

        chain_mirror_corr = np.concatenate((chain_2d_reflect_corrs_neg, chain_2d, chain_2d_reflect_corrs_pos))

        chain_2d_reflect_amp_pos = chain_mirror_corr.copy()
        chain_2d_reflect_amp_pos[:, 0] = 2*(-14) - chain_2d_reflect_amp_pos[:, 0]
        chain_2d_reflect_amp_neg = chain_mirror_corr.copy()
        chain_2d_reflect_amp_neg[:, 0] = 2*(-18) - chain_2d_reflect_amp_neg[:, 0]

        chain_mirror = np.concatenate((chain_2d_reflect_amp_neg, chain_mirror_corr, chain_2d_reflect_amp_pos))

        data = chain_mirror.squeeze()
        grid_points = 4096  # Grid points in each dimension
        # Bandwidth is the standard deviation of the Gaussian reweighting distribution.
        #     A different bandwidth/kernel can be used for the correlation
        kde = FFTKDE(kernel='gaussian', bw=0.05)
        grid, points = kde.fit(data).evaluate((grid_points, grid_points))
        # The grid is of shape (obs, dims), points are of shape (obs, 1)
        x, y = np.unique(grid[:, 0]), np.unique(grid[:, 1])
        # x is amplitudes
        # y is corr coeff
        z = points.reshape(grid_points, grid_points).T

        z[y<=-1, :] = 0  # Set the KDE to zero outside of the domain
        z[y>=1, :] = 0  # Set the KDE to zero outside of the domain
        z[:, x<=-18] = 0  # Set the KDE to zero outside of the domain
        z[:, x>=-14] = 0  # Set the KDE to zero outside of the domain
        z = z.squeeze() * 9  # multiply the kde to get integral of 1

        indx = np.argwhere((x>=-18)*(x<=-14)).squeeze()
        indy = np.argwhere((y>=-1)*(y<=1)).squeeze()

        y = y[indy]
        x = x[indx]
        z = z[indy, :][:, indx]

        for bw in ['silverman', 0.274]:

            for kernel in ['gaussian','epa']:

                nsamp = 100

                samples = np.array([])

                if nsamp > 1:
                    for _ in range(0, nsamp):
                        # draw an amplitude
                        amp = np.random.normal(loc=-14.69, scale=0.05)
                        indx = np.argmin(np.abs(x - amp)).squeeze()

                        pdf = z[:, indx]
                        pdf /= np.sum(pdf)

                        # Make sure number of samples matches original data, for bw computation
                        samps = np.random.choice(y, size=int(np.ceil(nsamps/nsamp)), p=pdf)
                        samples = np.concatenate((samples, samps))

                    samples_pos = samples.copy()
                    samples_pos = 2*(1) - samples_pos
                    samples_neg = samples.copy()
                    samples_neg = 2*(-1) - samples_neg

                    samples = np.concatenate((samples_neg, samples, samples_pos))

                    data = samples.squeeze()
                    grid_points = 4096  # Grid points in each dimension
                    kde = FFTKDE(kernel=kernel, bw=bw)
                    grid, points = kde.fit(data).evaluate(grid_points)
                    # The grid is of shape (obs, dims), points are of shape (obs, 1)
                    y2 = np.unique(grid)

                    z2 = points.reshape(grid_points).T

                    z2[y2<=-1] = 0  # Set the KDE to zero outside of the domain
                    z2[y2>=1] = 0  # Set the KDE to zero outside of the domain
                    z2 = z2.squeeze() * 3  # multiply the kde to get integral of ~1

                    indy = np.argwhere((y2>=-1)*(y2<=1)).squeeze()
                    y2 = y2[indy]

                    corr_hd2 = z2[indy] / np.sum(z2[indy])

                corr_hd2 /= np.mean(np.diff(y2))

                np.savez("corr_chains/{}_corr_hd_reweight_more1713_{}_{}.npz".format(i, kernel, bw), corr_hd2, y2)
                print("Saved:", "corr_chains/{}_corr_hd_reweight_more1713_{}_{}.npz".format(i, kernel, bw))
                print("")


if "fixA" in chain_name:
    neff = numerator**2 / denominator
    neff_68 = numerator_68**2 / denominator_68
    print("Number of effective pulsar pairs, using standard deviation = {}".format(neff))
    print("Number of effective pulsar pairs, using 68% confidence = {}".format(neff_68))





