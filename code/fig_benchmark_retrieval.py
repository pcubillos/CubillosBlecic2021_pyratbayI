import os
import sys

import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter, ScalarFormatter, MultipleLocator
import numpy as np
from scipy.ndimage.filters import gaussian_filter1d as gaussf
import scipy.interpolate as si

import pyratbay as pb
import pyratbay.constants as pc
import mc3.plots as mp
import mc3.utils as mu
import mc3.stats as ms

sys.path.append('../code')
import taurex_ariel_sim as tas


def adjust_spines(ax, spines):
    for loc, spine in ax.spines.items():
        if loc in spines:
            spine.set_position(('outward', 0))
        else:
            spine.set_color('none')  # don't draw spine
    if 'left' in spines:
        ax.yaxis.set_ticks_position('left')
    else:
        ax.yaxis.set_ticks([])
    if 'bottom' in spines:
        ax.xaxis.set_ticks_position('bottom')
    else:
        ax.xaxis.set_ticks([])


def posterior_pt(posterior, tmodel, tpars, ifree, pressure):
    nlayers = len(pressure)
    u, uind, uinv = np.unique(posterior[:,0], return_index=True,
        return_inverse=True)
    nsamples = len(u)

    # Evaluate posterior PT profiles:
    profiles = np.zeros((nsamples, nlayers), np.double)
    for i in range(nsamples):
        tpars[ifree] = posterior[uind[i]]
        profiles[i] = tmodel(tpars)
    # Get percentiles (for 1,2-sigma boundaries and median):
    low1   = np.zeros(nlayers, np.double)
    low2   = np.zeros(nlayers, np.double)
    median = np.zeros(nlayers, np.double)
    high1  = np.zeros(nlayers, np.double)
    high2  = np.zeros(nlayers, np.double)
    for i in range(nlayers):
        tpost = profiles[uinv,i]
        low2[i]   = np.percentile(tpost,  2.275)
        low1[i]   = np.percentile(tpost, 15.865)
        median[i] = np.percentile(tpost, 50.000)
        high1[i]  = np.percentile(tpost, 84.135)
        high2[i]  = np.percentile(tpost, 97.725)
    return np.array([median, low1, low2, high1, high2])


# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Read transmission data:
transmission = [
    'mcmc_transmission_TRAPPIST-1b.cfg',
    'mcmc_transmission_GJ1214b.cfg',
    'mcmc_transmission_GJ436b.cfg',
    'mcmc_transmission_WASP-107b.cfg',
    'mcmc_transmission_HAT-P-11b.cfg',
    'mcmc_transmission_HAT-P-26b.cfg',
    'mcmc_transmission_WASP-39b.cfg',
    'mcmc_transmission_HD189733b.cfg',
    'mcmc_transmission_HD209458b.cfg',
    'mcmc_transmission_KELT-11b.cfg',
    'mcmc_transmission_55Cnce.cfg',
    'mcmc_transmission_WASP-76b.cfg',
    ]
nplanets = len(transmission)
tnames = [mc_file.split('_')[-1].split('.')[0]
    for mc_file in transmission]

planet_data = np.loadtxt(
    '../inputs/retrieval_benchmark_transmission.txt', dtype=str)
planets = list(planet_data[:,0])
ifit = [6, 4, 7, 8, 9, 10, 11]
params = np.array(planet_data[:,ifit], np.double)

nobs = 52  # Number of Ariel observing filters
tr_npars = 7
tr_data = np.zeros((nplanets, nobs))
tr_uncert = np.zeros((nplanets, nobs))

tr_bestfit = [None for _ in range(nplanets)]
tr_true = np.zeros((nplanets, tr_npars))
tr_best = np.zeros((nplanets, tr_npars))
tr_posterior = np.zeros((nplanets, tr_npars), object)

for i in range(nplanets):
    pyrat = pb.run(transmission[i], init=True, no_logfile=True)
    tr_data[i] = pyrat.obs.data
    tr_uncert[i] = pyrat.obs.uncert
    with np.load(pyrat.ret.mcmcfile) as mcmc:
        tr_texnames = mcmc['texnames']
        tr_best[i] = mcmc['bestp']
        post, zchain, zmask = mu.burn(mcmc)
    tr_bestfit[i], _ = pyrat.eval(tr_best[i])
    for p in range(tr_npars):
        tr_posterior[i,p] = post[:,p]
    j = planets.index(tnames[i])
    tr_true[i] = params[j]


tr_PDF  = [None for _ in range(nplanets)]
tr_Xpdf = [None for _ in range(nplanets)]
tr_HPDmin = np.zeros((nplanets, tr_npars))
quantile = 0.683
for i in range(nplanets):
    tr_PDF[i], tr_Xpdf[i] = [], []
    for p in range(tr_npars):
        pdf, xpdf, tr_HPDmin[i,p] = ms.cred_region(tr_posterior[i,p], quantile)
        tr_PDF[i].append(pdf)
        tr_Xpdf[i].append(xpdf)

# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Read emission data:
emission = [
    'mcmc_emission_WASP-80b.cfg',
    'mcmc_emission_TrES-1b.cfg',
    'mcmc_emission_WASP-6b.cfg',
    'mcmc_emission_HAT-P-1b.cfg',
    'mcmc_emission_WASP-43b.cfg',
    'mcmc_emission_WASP-94Ab.cfg',
    'mcmc_emission_WASP-17b.cfg',
    'mcmc_emission_WASP-31b.cfg',
    'mcmc_emission_HD149026b.cfg',
    'mcmc_emission_WASP-14b.cfg',
    'mcmc_emission_WASP-121b.cfg',
    'mcmc_emission_WASP-12b.cfg',
    ]
nplanets = len(emission)
enames = [mc_file.split('_')[-1].split('.')[0]
    for mc_file in emission]

em_npars = 8
nlayers = 22
em_true = np.zeros((nplanets, em_npars))
em_best = np.zeros((nplanets, em_npars))
em_posterior = np.zeros((nplanets, em_npars), object)
tp_post = np.zeros((nplanets, 5, nlayers))
tp_best = np.zeros((nplanets, nlayers))
planet_data = np.loadtxt(
    '../inputs/retrieval_benchmark_emission.txt', dtype=str)
planets = list(planet_data[:,0])
ifit = [7, 8, 6, 4, 9, 10, 11, 12]
params = np.array(planet_data[:,ifit], np.double)

em_data = np.zeros((nplanets, nobs))
em_uncert = np.zeros((nplanets, nobs))
em_bestfit = [None for _ in range(nplanets)]

for i, mc_file in enumerate(emission):
    pyrat = pb.run(mc_file, init=True, no_logfile=True)
    em_data[i] = pyrat.obs.data
    em_uncert[i] = pyrat.obs.uncert

    with np.load(pyrat.ret.mcmcfile) as mcmc:
        ifree = mcmc['ifree']
        em_texnames = mcmc['texnames'][ifree]
        em_best[i] = mcmc['bestp'][ifree]
        post, zchain, zmask = mu.burn(mcmc)
        rplanet = mcmc['bestp'][pyrat.ret.irad] * pc.rearth
        tp_best[i] = pyrat.atm.tmodel(mcmc['bestp'][pyrat.ret.itemp])
    for p in range(em_npars):
        em_posterior[i,p] = post[:,p]
    em_true[i] = params[i]
    bestp = pyrat.ret.params
    bestp[ifree] = em_best[i]
    em_bestfit[i], _ = pyrat.eval(bestp)
    em_bestfit[i] *= (rplanet/pyrat.phy.rstar)**2 / pyrat.spec.starflux

    #ifree = pyrat.ret.pstep[pyrat.ret.itemp] > 0
    #itemp_free = np.arange(np.sum(ifree))
    itemp_free = [itemp in ifree for itemp in pyrat.ret.itemp]
    ntemp_free = np.sum(itemp_free)
    tp_post[i] = posterior_pt(
        post[:,0:ntemp_free],
        pyrat.atm.tmodel,
        pyrat.ret.params[pyrat.ret.itemp],
        itemp_free,
        pyrat.atm.press)


PDF  = [None for _ in range(nplanets)]
Xpdf = [None for _ in range(nplanets)]
HPDmin = np.zeros((nplanets, em_npars))
for i in range(nplanets):
    PDF[i], Xpdf[i] = [], []
    for p in range(em_npars):
        pdf, xpdf, HPDmin[i,p] = ms.cred_region(em_posterior[i,p], quantile)
        PDF[i].append(pdf)
        Xpdf[i].append(xpdf)

# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Plot spectra:

wl = 1.0/(pyrat.spec.wn*pc.um)
bandwl = 1.0/(pyrat.obs.bandwn*pc.um)

margin = 0.05
sigma = 8.0
fs = 10
lw = 1.0
logxticks = [0.5, 1.0, 2.0, 3.0, 5.0 , 8.0]
rect1 = [0.075, 0.54, 0.985, 0.99]
rect2 = [0.075, 0.04, 0.985, 0.49]

fig = plt.figure(8, (8.2, 10))
plt.clf()
axes = []
for i in range(nplanets):
    ax = mp.subplotter(rect1, margin, i+1, nx=3, ny=4, ymargin=0.025)
    axes.append(ax)
    if tr_bestfit[i] is not None:
        plt.plot(wl, gaussf(tr_bestfit[i], sigma)/pc.percent, lw=lw, c='orange')
    plt.errorbar(bandwl, tr_data[i]/pc.percent, tr_uncert[i]/pc.percent,
        fmt='o', alpha=0.7, ms=1.5, color='b', ecolor='0.3',
        elinewidth=lw, capthick=lw, zorder=3)
    yran = ax.get_ylim()
    ax.tick_params(labelsize=fs-1, direction='in', which='both')
    ax.set_xscale('log')
    plt.gca().xaxis.set_minor_formatter(NullFormatter())
    ax.get_xaxis().set_major_formatter(ScalarFormatter())
    ax.set_xticks(logxticks)
    plt.xlim(np.amin(wl), np.amax(wl))
    plt.text(0.03, 0.86, tnames[i], fontsize=fs-1, ha='left',
        transform=ax.transAxes)
    if i%3 == 0:
        plt.ylabel(r'$(R_{\rm p}/R_{\rm s})^2$  (%)', fontsize=fs)
    if i >= 9:
        ax.set_xlabel(r'Wavelength (um)', fontsize=fs)
for i in range(nplanets):
    ax = mp.subplotter(rect2, margin, i+1, nx=3, ny=4, ymargin=0.025)
    axes.append(ax)
    if em_bestfit[i] is not None:
        plt.plot(wl, gaussf(em_bestfit[i], sigma)/pc.ppt,lw=lw, c='orange')
    plt.errorbar(bandwl, em_data[i]/pc.ppt, em_uncert[i]/pc.ppt,
        fmt='o', alpha=0.7, ms=1.5, color='b', ecolor='0.3',
        elinewidth=lw, capthick=lw, zorder=3)
    yran = ax.get_ylim()
    ax.set_xscale('log')
    plt.gca().xaxis.set_minor_formatter(NullFormatter())
    ax.get_xaxis().set_major_formatter(ScalarFormatter())
    ax.set_xticks(logxticks)
    ax.tick_params(labelsize=fs-1, direction='in', which='both')
    plt.xlim(np.amin(wl), np.amax(wl))
    plt.text(0.03, 0.86, enames[i], fontsize=fs-1, ha='left',
        transform=ax.transAxes)
    if i%3 == 0:
        ax.set_ylabel(r'$F_{\rm p}/F_{\rm s}$ (ppt)', fontsize=fs)
    if i >= 9:
        ax.set_xlabel(r'Wavelength (um)', fontsize=fs)
axes[3].yaxis.set_major_locator(MultipleLocator(0.2))
axes[4].yaxis.set_major_locator(MultipleLocator(0.02))
axes[5].yaxis.set_major_locator(MultipleLocator(0.1))
axes[12].yaxis.set_major_locator(MultipleLocator(0.5))
axes[13].yaxis.set_major_locator(MultipleLocator(0.5))
axes[14].yaxis.set_major_locator(MultipleLocator(1.0))
axes[15].yaxis.set_major_locator(MultipleLocator(0.5))
axes[19].yaxis.set_major_locator(MultipleLocator(1.0))
axes[20].yaxis.set_major_locator(MultipleLocator(0.25))
plt.savefig('../plots/pyratbay-taurex_ret_spectra_comparision.pdf')
# Lower temperature+higher radius(p0) flattens IR bands
# (at constant abundances)




# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Transmission posterior plot:

hkw = {'histtype':'step', 'lw':1.25, 'density':False}
bkw1 = {'color':'crimson', 'dashes':(5,1), 'lw':1.5}
bkw2 = {'color':'green', 'dashes':(), 'lw':1.5}
fs = 10
tmin, tmax = 0.4, 0.6


width = np.array([
#   TRA1  G1214  G436  W107  H11   H26  W39  HD1  HD2  K11   55C   W76
    [ 80,  150,  180,  70,   355,  300, 110, 860, 270, 420, 2390, 1000], # T
    [0.09, 0.12, 0.06, 0.18, 0.18, 0.9, 0.2, 0.3, 0.2, 0.8,  1.0,  0.9], # R
    [1.1,  0.6,  0.7,  0.4,  2.0,  0.0, 0.7, 3.2, 0.7, 0.7,  6.0,  1.5], # H2O
    [0.7,  0.5,  0.7,  0.3,  1.7,  1.7, 0.0, 0.0, 0.7, 0.4,  0.0,  0.0], # CH4
    [0.0,  0.0,  0.0,  1.0,  0.0,  6.0, 0.0, 5.0, 2.0, 1.0,  0.0,  4.0], # CO
    [0.0,  1.7,  5.0,  2.2,  0.0,  0.0, 0.0, 7.0, 1.2, 0.0,  0.0,  0.0], # CO2
    [3.0,  0.8,  0.0,  1.6,  0.0,  2.5, 0.0, 3.5, 1.7, 0.4,  0.0,  2.5], # cl
    ])
ranges = np.tile(None, (tr_npars, nplanets, 2))

for p in range(tr_npars):
    Tmax, Tmin = np.amax(tr_true[:,p]), np.amin(tr_true[:,p])
    if 2 <= p <= 5:
        Tmin = np.amin(tr_true[:,p][tr_true[:,p]>-12])
    if p == 6:
        Tmax = np.amax(tr_true[:,p][tr_true[:,p]<2.0])
    a = (tmax-tmin)/(Tmax-Tmin)
    b = tmax - a*Tmax
    t = a*tr_true[:,p] + b
    ranges[p,:,0] = tr_true[:,p] - t*width[p]
    ranges[p,:,1] = tr_true[:,p] + (1-t)*width[p]
    if p == 6:
        ranges[p,:,1] = tr_true[:,p] - t*width[p]
        ranges[p,:,0] = tr_true[:,p] + (1-t)*width[p]

ranges[2,width[2]==0] = -12.3, -1.0
ranges[3,width[3]==0] = -12.3, -1.0
ranges[4,width[4]==0] = -12.3, -1.0
ranges[5,width[5]==0] = -12.3, -1.0
ranges[4,3] = -12.3, -1.0
ranges[5,7] = -12.3, -1.0
idx = 0, 2, 4, 6, 10
ranges[6,idx] = 2.1, -3.5

fig = plt.figure(10, (8.2, 10))
plt.clf()
axes = np.empty((tr_npars, nplanets), object)
plt.subplots_adjust(0.07, 0.05, 0.99, 0.99, hspace=0.07, wspace=0.0)
for p in range(tr_npars):
    for i in range(nplanets):
        ecol = 'k' if width[p,i] == 0 else 'b'
        axes[p,i] = ax = fig.add_subplot(tr_npars, 12, i+1 + 12*p)
        vals, bins, h = ax.hist(tr_posterior[i,p], bins=25, zorder=0,
            orientation='horizontal', edgecolor=ecol, **hkw)
        vals = np.r_[0, vals, 0]
        bins = np.r_[bins[0] - (bins[1]-bins[0]), bins]
        f = si.interp1d(bins+0.5*(bins[1]-bins[0]), vals, kind='nearest')
        ax.fill_betweenx(tr_Xpdf[i][p], 0, f(tr_Xpdf[i][p]),
            where=tr_PDF[i][p]>=tr_HPDmin[i,p], facecolor='0.75',
            edgecolor='none',
            interpolate=False, zorder=-2)
        xtop = np.amax(vals)
        ax.axhline(tr_true[i,p], xmax=0.5, **bkw2)
        ax.axhline(tr_best[i,p], xmax=0.5, **bkw1)
        ax.set_xlim(0, 2.2*xtop)
        ax.set_xticklabels([])
        ax.set_xticks([])
        adjust_spines(ax, ['left', 'bottom'])
        ax.tick_params(labelsize=fs-2, direction='in')
        if width[p,i] == 0 and p != tr_npars-1:
            ax.yaxis.set_major_locator(MultipleLocator(3))
        if i == 0:
            ax.set_ylabel(tr_texnames[p], fontsize=fs)
        if ranges[p,i,0] is not None:
            ax.set_ylim(ranges[p,i])
        if p == tr_npars - 1:
            ax.set_xlabel(tnames[i], fontsize=fs-3)

axes[1,1].yaxis.set_major_locator(MultipleLocator(0.03))
axes[1,3].set_yticks(np.arange(10.25, axes[1,3].get_ylim()[1], 0.04))
axes[2,0].yaxis.set_major_locator(MultipleLocator(0.3))
axes[3,5].yaxis.set_major_locator(MultipleLocator(0.4))
axes[4,3].yaxis.set_major_locator(MultipleLocator(3))
axes[5,7].yaxis.set_major_locator(MultipleLocator(3))
plt.savefig('../plots/pyratbay-taurex_ret_transmission_comparision.pdf')


# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Emission posterior plot:

# TauREx temperature profiles:
tas.setup_cache()

tau_temp = np.zeros((nplanets, nlayers))
sys_params = np.array(planet_data[:,1:7], np.double)
for i,planet in enumerate(enames):
    j = planets.index(planet)
    rstar, tstar, kmag, rplanet, mplanet, tplanet = sys_params[j]
    tm, wn_tau, depth = tas.simulate_ariel(
        rplanet*pc.rearth, mplanet*pc.mearth, params[j,2], params[j,4:],
        100, rstar, tstar, tpars=params[j,0:2], mode='eclipse')
    tau_temp[i] = tm.temperatureProfile
tau_press = tm.pressureProfile * pc.pascal/pc.bar


rect = [0.07, 0.05, 0.99, 0.38]
margin = 0.02

width = np.array([
#    W80   Tr1   W6    H1   W43  W94  W17  W31  149  W14  W121  W12
    [150,  195,  80,   300, 110, 900, 900, 300, 400, 400, 2340, 2340], # kap
    [0.12, 0.06, 0.15, 0.9, 0.2, 0.3, 0.3, 0.2, 0.8, 0.8,  1.0,  1.0], # gam
    [0.6,  0.7,  0.4,  0.0, 0.7, 3.2, 3.2, 0.7, 0.7, 0.7,  6.0,  6.0], # Tirr
    [1.8,  1.6,  3.0,  1.3, 0.8, 2.5, 5.0, 3.5, 1.2, 2.2,  3.0,  2.8], # R
    [5.5,  4.5,  4.2,  0.0, 3.0, 4.0, 3.0, 6.0, 3.8, 5.0,  4.0,  5.8], # H2O
    [4.5,  0.0,  0.0,  0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0,  0.0,  0.0], # CH4
    [1.0,  1.0,  0.0,  0.0, 5.0, 1.0, 2.5, 0.0, 0.0, 1.0,  1.0,  1.0], # CO
    [1.0,  1.0,  6.0,  0.0, 3.0, 0.0, 0.0, 0.0, 1.0, 0.0,  0.0,  1.0], # CO2
    ])
ranges = np.tile(None, (em_npars, nplanets,2))

for p in range(em_npars):
    Tmax, Tmin = np.amax(em_true[:,p]), np.amin(em_true[:,p])
    if 2 <= p <= 5:
        Tmin = np.amin(em_true[:,p][em_true[:,p]>-12])
    if p == 6:
        Tmax = np.amax(em_true[:,p][em_true[:,p]<2.0])
    a = (tmax-tmin)/(Tmax-Tmin)
    b = tmax - a*Tmax
    t = a*em_true[:,p] + b
    ranges[p,:,0] = em_true[:,p] - t*width[p]
    ranges[p,:,1] = em_true[:,p] + (1-t)*width[p]
    if p == 6:
        ranges[p,:,0] = em_true[:,p] - t*width[p]
        ranges[p,:,1] = em_true[:,p] + (1-t)*width[p]

for i in range(nplanets):
    if enames[i] not in ['WASP-17b']:
        ranges[4,i] = -5.8, -0.99
ranges[4,:] = -7, -0.99
for p in range(4, em_npars):
    ranges[p,width[p]==0] = -12.3, -0.99
ranges[6,width[6]==1] = -12.3, -0.99
ranges[7,width[7]==1] = -12.3, -0.99
ranges[5,6] = -5.3, -3.9
ranges[6,4] = -4.2, -0.95
ranges[6,6] = -4.5, -1.9
ranges[7,4] = -7.2, -4.9

fig = plt.figure(11, (8.2, 10))
plt.clf()
axes = np.empty((em_npars, nplanets), object)
plt.subplots_adjust(0.07, 0.05, 0.99, 0.99, hspace=0.07, wspace=0.0)
for p in range(3,em_npars):
    for i in range(nplanets):
        ecol = 'k' if width[p,i] == 0 else 'b'
        axes[p,i] = ax = fig.add_subplot(em_npars, 12, i+1 + 12*(p-3))
        vals, bins, h = ax.hist(
            em_posterior[i,p], bins=25, range=None, zorder=0,
            orientation='horizontal', edgecolor=ecol, **hkw)
        vals = np.r_[0, vals, 0]
        bins = np.r_[bins[0] - (bins[1]-bins[0]), bins]
        f = si.interp1d(bins+0.5*(bins[1]-bins[0]), vals, kind='nearest')
        ax.fill_betweenx(
            Xpdf[i][p], 0, f(Xpdf[i][p]), where=PDF[i][p]>=HPDmin[i,p],
            facecolor='0.75', edgecolor='none', interpolate=False, zorder=-2)
        xtop = np.amax(vals)
        ax.axhline(em_true[i,p], xmax=0.5, **bkw2)
        ax.axhline(em_best[i,p], xmax=0.5, **bkw1)
        ax.set_xlim(0, 2.2*xtop)
        ax.set_xticks([])
        adjust_spines(ax, ['left', 'bottom'])
        ax.tick_params(labelsize=fs-2, direction='in')
        if ranges[p,i,0] == -12.3:
            ax.yaxis.set_major_locator(MultipleLocator(3))
        if i == 0:
            ax.set_ylabel(em_texnames[p], fontsize=fs)
        if ranges[p,i,0] is not None:
            ax.set_ylim(ranges[p,i])
        if p == em_npars - 1:
            ax.set_xlabel(enames[i], fontsize=fs-3)
for i in range(nplanets):
    ax = mp.subplotter(rect, margin, i+1, nx=6, ny=2)
    ax.fill_betweenx(
        pyrat.atm.press/pc.bar, tp_post[i,2], tp_post[i,4],
        facecolor='royalblue', edgecolor='none', alpha=0.6)
    ax.fill_betweenx(
        pyrat.atm.press/pc.bar, tp_post[i,1], tp_post[i,3],
        facecolor='royalblue', edgecolor='none', alpha=0.8)
    plt.plot(tp_post[i,0], pyrat.atm.press/pc.bar,'navy',lw=1.5,label='Median')
    ax.plot(tau_temp[i], tau_press, c='limegreen', lw=1.5)
    plt.plot(tp_best[i], pyrat.atm.press/pc.bar, **bkw1)
    ax.set_ylim(100, 1e-5)
    ax.set_xlim(200, 3150)
    ax.set_yscale('log')
    if i % 6:
        ax.set_yticklabels([])
    ax.tick_params(labelsize=fs-2, direction='in')
    plt.text(
        0.98, 0.92, enames[i], fontsize=fs-3, ha='right',
        transform=ax.transAxes)
    if i == 0:
        plt.text(
            -0.42, -0.4, 'Pressure (bar)', fontsize=fs,
            rotation='vertical', transform=ax.transAxes)
        plt.text(
            2.8, -1.35, 'Temperature (K)', fontsize=fs, transform=ax.transAxes)
plt.savefig('../plots/pyratbay-taurex_ret_emission_comparision.pdf')

