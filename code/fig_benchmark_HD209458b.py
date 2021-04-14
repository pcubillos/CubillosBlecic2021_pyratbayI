from copy import copy

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage.filters import gaussian_filter1d as gaussf

import pyratbay as pb
import pyratbay.atmosphere as pa
import pyratbay.constants as pc
import mc3


# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Read transmission data:
mm2017 = pb.run(
    'mcmc_transmission_HD209458b_MM2017.cfg', run_step='init', no_logfile=True)
mcmc_mm2017 = np.load(mm2017.ret.mcmcfile)
bestp_mm2017 = mcmc_mm2017['bestp']

pyrat = pb.run(
    'mcmc_transmission_HD209458b_P2019.cfg', run_step='init', no_logfile=True)

bf = pyrat.ret.mcmcfile.replace('.npz', '_bestfit_spectrum.dat')
mcmc_p2019 = np.load(pyrat.ret.mcmcfile)
wn, best_spec_p2019 = pb.io.read_spectrum(bf)

bestp = np.copy(pyrat.ret.params)
for pname, val in zip(mcmc_mm2017['pnames'], mcmc_mm2017['bestp']):
    if pname in mcmc_p2019['pnames']:
        bestp[list(mcmc_p2019['pnames']).index(pname)] = val
bestp[list(mcmc_p2019['pnames']).index('log(CO)')] = -12.0
bestp[list(mcmc_p2019['pnames']).index('log(CO2)')] = -12.0
pyrat.atm.refpressure = mm2017.atm.refpressure
best_spec_mm2017, _ = pyrat.eval(bestp)

# Temperature posteriors:
mm2017_posterior, _, _ = mc3.utils.burn(mcmc_mm2017)
ifree = mm2017.ret.pstep[mm2017.ret.itemp] > 0
mm2017_tpost = pa.temperature_posterior(
    mm2017_posterior[:,mm2017.ret.itemp], mm2017.atm.tmodel,
    mm2017.ret.params[mm2017.ret.itemp],
    ifree, mm2017.atm.press)

p2019_posterior, _, _ = mc3.utils.burn(mcmc_p2019)
ifree = pyrat.ret.pstep[pyrat.ret.itemp] > 0
p2019_tpost = pa.temperature_posterior(
    p2019_posterior[:,pyrat.ret.itemp], pyrat.atm.tmodel,
    pyrat.ret.params[pyrat.ret.itemp],
    ifree, pyrat.atm.press)


nbins = 16
posterior1 = p2019_posterior[:,0:6]
nsamples1, npars1 = np.shape(posterior1)
pnames1 = mcmc_p2019['texnames'][0:6]
pnames1 = [pname.replace('log_{10}', 'log') for pname in pnames1]
ranges1 = [
    [-6, 2], [-6, 2], [-2, 2], [0.01, 1], [0.01, 1], [800,1650],]

hist1 = []
xran1, yran1, lmax1 = [], [], []
for irow in range(1, npars1):
    for icol in range(irow):
        ran = None
        if ranges1[icol] is not None:
            ran = [ranges1[icol], ranges1[irow]]
        h, x, y = np.histogram2d(
            posterior1[:,icol], posterior1[:,irow], bins=nbins, range=ran)
        hist1.append(h.T)
        xran1.append(x)
        yran1.append(y)
        lmax1.append(np.amax(h)+1)


texnames_mm2017 = [
    pname.replace('Ray','ray').replace('(f_','(X_')
    for pname in mcmc_mm2017['texnames']]
texnames_p2019 = [
    pname.replace('Ray','ray').replace('(f_','(X_')
    for pname in mcmc_p2019['texnames']]
texnames_mm2017[13] = '$\\log_{10}(f_{\\rm ray})$'
texnames_p2019[15]  = '$\\log_{10}(f_{\\rm ray})$'

posterior2 = p2019_posterior[:,6:]
nsamples2, npars2 = np.shape(posterior2)
pnames2 = [
    pname.replace('log_{10}','log').replace('planet','p')
    for pname in copy(texnames_p2019[6:])]

ranges2 = [
    [1.33, 1.39], [-7, -2], [-9, -2.9], [-6, -2], [-10, -2], [-10, -2],
    [-10, -2], [-10, -2], [-10,-2], [0, 8], [-20, 2], [-6,-2], [0.1, 0.9]]

hist2 = []
xran2, yran2, lmax2 = [], [], []
for irow in range(1, npars2):
    for icol in range(irow):
        ran = None
        if ranges2[icol] is not None:
            ran = [ranges2[icol], ranges2[irow]]
        h, x, y = np.histogram2d(posterior2[:,icol], posterior2[:,irow],
            bins=nbins, range=ran)
        hist2.append(h.T)
        xran2.append(x)
        yran2.append(y)
        lmax2.append(np.amax(h)+1)


p2019_median = [
    -1.85, -4.22, 0.05, 0.59, 0.35, 949.0,
    1.359, -4.92, -6.46, -4.66, -8.60, -8.56, -8.66, -8.37, -9.8,
    4.57, -14.82, -4.47, 0.52] 

p2019_lo = [
    1.28, 1.22, 1.07, 0.2, 0.17, 109,
    0.0, 0.57, 0.64, 0.3, 2.22, 2.28, 2.21, 2.38, 1.45,
    0.74, 3.45, 0.48, 0.07]

p2019_hi = [
    1.66, 1.8, 1.29, 0.27, 0.38, 252,
    0.0, 0.83, 0.84, 0.39, 2.2, 2.22, 2.25, 2.57, 1.57,
    0.58, 4.79, 0.52, 0.06]

mm2017_median = [
    -1.46, -4.03, 0.62, 0.67, 0.57, 1071.0,
    1.359, -5.13, -20, -5.24, -7.84, -6.03, -6.35, -20, -20,
    4.44, -15.03, -4.34, 0.47]

mm2017_lo = [
    1.83, 1.38, 1.35, 0.18, 0.26, 161,
    0.0, 0.57, 0, 0.27, 1.44, 1.88, 2.40, 0, 0,
    0.95, 3.36, 0.53, 0.0]

mm2017_hi = [
    1.68, 2.05, 0.97, 0.21, 0.28, 149,
    0.0, 0.73, 0, 0.36, 1.5, 0.46, 0.85, 0, 0,
    0.63, 4.65, 0.6, 0.0]


# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Plot:

wl = 1.0/(pyrat.spec.wn*pc.um)
bandwl = 1.0/(pyrat.obs.bandwn*pc.um)

sigma = 12.0
fs = 9.75
lw = 1.0
nlevels = 20
margin = 0.003
nb = 14

nhist = npars1 + npars2
margin2 = 0.01
ymargin = 0.03
nx = 7
ny = 3

rect1 = [0.073, 0.8, 0.615, 0.19]
rect2 = [0.76, 0.8, 0.23, 0.19]
rect3 = [0.67, 0.51, 1.02, 0.76]
rect4 = [0.115, 0.328, 0.84, 0.823]
rect5 = [0.02, 0.03, 0.99, 0.28]

palette = copy(plt.cm.viridis_r)
palette.set_under(color='w')
palette.set_bad(color='w')
p2019_col = pb.plots.alphatize('r', 0.6, 'darkred')


fig = plt.figure(18, (8.5, 11.0))
plt.clf()
ax = plt.axes(rect1)
ax.plot(wl, gaussf(best_spec_mm2017, sigma)/pc.percent, lw=lw,
    c='royalblue', label='pyratbay fit to MM2017')
ax.plot(wl, gaussf(best_spec_p2019, sigma)/pc.percent, lw=lw,
    c='orange', label='pyratbay fit to P2019')
ax.errorbar(
    bandwl, pyrat.obs.data/pc.percent, pyrat.obs.uncert/pc.percent,
    fmt='o', alpha=0.7, ms=2.5, color=p2019_col,
    elinewidth=lw, capthick=lw, zorder=3, label='data')

ax.tick_params(labelsize=fs-1, direction='in', which='both')
ax.tick_params(length=0, which='minor')
ax.set_xscale('log')
plt.gca().xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
ax.set_xticks(pyrat.inputs.logxticks)
ax.set_xlim(0.29, 5.5)
ax.set_ylabel(r'$(R_{\rm p}/R_{\rm s})^2$  (%)', fontsize=fs)
ax.set_xlabel(r'Wavelength (um)', fontsize=fs)
ax.legend(loc='upper right', fontsize=fs-1)

ax2 = plt.axes(rect2)
pb.plots.temperature(
    mm2017.atm.press, profiles=mm2017_tpost[0:1], labels=['MM2017'],
    bounds=mm2017_tpost[1:3], ax=ax2, theme='blue', fs=fs,
    colors=['mediumblue'])
pb.plots.temperature(
    pyrat.atm.press, profiles=p2019_tpost[0:1], labels=['P2019'],
    bounds=p2019_tpost[1:3], ax=ax2, theme='orange', fs=fs,
    colors=['darkorange'], alpha=[0.45, 0.25])
ax2.tick_params(labelsize=fs-1, direction='in', which='both')
ax2.tick_params(length=0, which='minor')
ax2.set_xlim(800, 1900)

axes = np.tile(None, (npars2, npars2))
k = 0 # Histogram index
for irow in range(1, npars2):
    for icol in range(irow):
        h = npars2*irow + icol + 1  # Subplot index
        ax = mc3.plots.subplotter(rect4, margin, h, npars2)
        ax.tick_params(labelsize=fs-3, direction='in')
        if icol == 0:
            ax.set_ylabel(pnames2[irow], size=fs-2.5,
                ha='right', va='center', rotation='horizontal')
        else:
            ax.get_yaxis().set_visible(False)
        if irow == npars2-1:
            ax.set_xlabel(pnames2[icol], size=fs-2.5)
            plt.setp(ax.xaxis.get_majorticklabels(), rotation=90)
        else:
            ax.get_xaxis().set_visible(False)
        cont = ax.contourf(hist2[k], cmap=palette, vmin=1, origin='lower',
            levels=[0]+list(np.linspace(1,lmax2[k], nlevels)),
            extent=(xran2[k][0], xran2[k][-1], yran2[k][0], yran2[k][-1]))
        for c in cont.collections:
            c.set_edgecolor("face")
        ax.set_xlim(ranges2[icol])
        ax.set_ylim(ranges2[irow])
        k += 1

axes = np.tile(None, (npars1, npars1))
k = 0 # Histogram index
for irow in range(1, npars1):
    for icol in range(irow):
        h = npars1*irow + icol + 1  # Subplot index
        ax = mc3.plots.subplotter(rect3, margin, h, npars1)
        ax.tick_params(labelsize=fs-3, direction='in')
        if icol == 0:
            ax.set_ylabel(pnames1[irow], size=fs-2)
        else:
            ax.get_yaxis().set_visible(False)
        if irow == npars1-1:
            ax.set_xlabel(pnames1[icol], size=fs-2)
            plt.setp(ax.xaxis.get_majorticklabels(), rotation=90)
        else:
            ax.get_xaxis().set_visible(False)
        cont = ax.contourf(hist1[k], cmap=palette, vmin=1, origin='lower',
            levels=[0]+list(np.linspace(1,lmax1[k], nlevels)),
            extent=(xran1[k][0], xran1[k][-1], yran1[k][0], yran1[k][-1]))
        for c in cont.collections:
            c.set_edgecolor("face")
        ax.set_xlim(ranges1[icol])
        ax.set_ylim(ranges1[irow])
        k += 1

bounds = np.linspace(0, 1.0, nlevels)
norm = matplotlib.colors.BoundaryNorm(bounds, palette.N)
ax2 = plt.axes([0.58, 0.57, 0.015, 0.15])
cb = matplotlib.colorbar.ColorbarBase(ax2, cmap=palette, norm=norm,
    spacing='proportional', boundaries=bounds, format='%.1f')
cb.set_label("P2019 Posterior Density", fontsize=fs-1)
cb.ax.yaxis.set_ticks_position('left')
cb.ax.yaxis.set_label_position('left')
cb.ax.tick_params(labelsize=fs-2)
cb.set_ticks(np.linspace(0, 1, 5))
for c in ax2.collections:
    c.set_edgecolor("face")


axes = [
    mc3.plots.subplotter(rect5, margin2, j+1, nx=nx, ny=ny, ymargin=ymargin)
    for j in range(nhist)]
idx = np.array([i for i in range(nhist) if i not in [13,14,18]])
zaxes = mc3.plots.histogram(
    mm2017_posterior,
    pnames=texnames_mm2017, ranges=np.array(ranges1+ranges2)[idx],
    quantile=0.683,
    axes=np.array(axes)[idx], fs=fs-2.0, lw=1.2, yscale=None, theme='blue')

axes = mc3.plots.histogram(
    p2019_posterior,
    pnames=texnames_p2019, ranges=ranges1+ranges2, quantile=0.683,
    axes=axes, fs=fs-2.0, lw=1.2, yscale=None, theme='orange')
axes[8].set_xticks([-8, -6, -4])

for ax, med, lo, hi in zip(axes, mm2017_median, mm2017_lo, mm2017_hi):
    y = 0.2*ax.get_ylim()[1]
    ax.errorbar([med], [y], xerr=[[lo],[hi]], color='navy', marker='D', ms=3)
for ax, med, lo, hi in zip(axes, p2019_median, p2019_lo, p2019_hi):
    ax.tick_params(length=0, axis='y')
    plt.setp(ax.xaxis.get_majorticklabels(), rotation=0)
    y = 0.1*ax.get_ylim()[1]
    ax.errorbar([med], [y], xerr=[[lo],[hi]], color=p2019_col, marker='D', ms=3)
axes[0].legend(
    ['pyratbay fit to MM2017', 'pyratbay fit to P2019', None, None,
     'MacDonald & Madhusudhan (2017)', 'Pinhas et al. (2019)'],
    loc=(5.4, -3.5), fontsize=fs-2)
leg = axes[0].get_legend()
handles = leg.legendHandles
axes[0].legend(
    [(handles[0],handles[2]), (handles[1],handles[3]), handles[4], handles[5]],
    ['pyratbay fit to MM2017', 'pyratbay fit to P2019',
     'MacDonald & Madhusudhan (2017)', 'Pinhas et al. (2019)'], 
    loc=(5.4, -2.95), fontsize=fs-2)

plt.savefig('../plots/pyratbay-madhu_HD209458b_comparison.pdf')
plt.savefig('../plots/pyratbay-madhu_HD209458b_comparison.png', dpi=300)
