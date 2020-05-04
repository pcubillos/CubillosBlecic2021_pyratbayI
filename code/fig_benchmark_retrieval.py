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

sys.path.append('code')
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

planet_data = np.loadtxt('retrieval_validation_transmission.txt', dtype=str)
planets = list(planet_data[:,0])
ifit = [6, 4, 7, 8, 9, 10, 11]
params = np.array(planet_data[:,ifit], np.double)

tnpars = 7
tdata, tuncert = [], []
tbest = [None for _ in range(nplanets)]
ttrue_vals = np.zeros((nplanets, tnpars))
tbest_vals = np.zeros((nplanets, tnpars))
tlo_vals   = np.zeros((nplanets, tnpars))
thi_vals   = np.zeros((nplanets, tnpars))
tposterior = np.zeros((nplanets, tnpars), object)

for i in range(nplanets):
    pyrat = pb.run(transmission[i], init=True, no_logfile=True)
    tdata.append(pyrat.obs.data)
    tuncert.append(pyrat.obs.uncert)
    bf = pyrat.ret.mcmcfile.replace('.npz', '_bestfit_spectrum.dat')
    if not os.path.exists(bf):
        continue
    if transmission[i] == 'mcmc_transmission_GJ436b.cfg':
        continue

    mcmc = np.load(pyrat.ret.mcmcfile)
    wn, best_spec = pb.io.read_spectrum(bf)
    tbest[i] = best_spec
    post, zchain, zmask = mu.burn(mcmc)
    for p in range(tnpars):
        tposterior[i,p] = post[:,p]
    tbest_vals[i] = mcmc['bestp']
    tlo_vals[i] = mcmc['CRlo']
    thi_vals[i] = mcmc['CRhi']
    texnames = mcmc['texnames']
    j = planets.index(tnames[i])
    ttrue_vals[i] = params[j]
    mcmc.close()


PDF  = [None for _ in range(nplanets)]
Xpdf = [None for _ in range(nplanets)]
HPDmin = np.zeros((nplanets, tnpars))
quantile = 0.683
for i in range(nplanets):
    if tbest[i] is None:
        continue
    PDF[i], Xpdf[i] = [], []
    for p in range(tnpars):
        pdf, xpdf, HPDmin[i,p] = ms.cred_region(tposterior[i,p], quantile)
        PDF[i].append(pdf)
        Xpdf[i].append(xpdf)

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
enames = [mc_file.split('_')[-1].split('.')[0]
    for mc_file in emission]

edata, euncert = [], []
ebest = []
for mc_file in emission:
    pyrat = pb.run(mc_file, init=True, no_logfile=True)
    edata.append(pyrat.obs.data)
    euncert.append(pyrat.obs.uncert)
    bf = pyrat.ret.mcmcfile.replace('.npz', '_bestfit_spectrum.dat')
    try:
        wn, best_spec = pb.io.read_spectrum(bf.replace('retri', 'old_retri'))
        #wn, best_spec = pb.io.read_spectrum(bf)
        with np.load(pyrat.ret.mcmcfile.replace('retri', 'old_retri')) as d:
            rplanet = d['bestp'][5] * pc.rearth
        best_spec *= (rplanet/pyrat.phy.rstar)**2 / pyrat.spec.starflux
    except:
        best_spec = None
    ebest.append(best_spec)

wl = 1.0/(wn*pc.um)
bandwl = 1.0/(pyrat.obs.bandwn*pc.um)


# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Plot spectra:
rect = [0.075, 0.05, 0.985, 0.99]
margin = 0.05
sigma = 8.0
fs = 10
lw = 1.0
ms = 1.5
logxticks = [0.5, 1.0, 2.0, 3.0, 5.0 , 8.0]

fig = plt.figure(9, (8.2, 10))
plt.clf()
axes = []
for i in range(nplanets):
    ax = mp.subplotter(rect, margin, i+1, nx=3, ny=8, ymargin=0.025)
    axes.append(ax)
    if tbest[i] is not None:
        plt.plot(wl, gaussf(tbest[i], sigma)/pc.percent, lw=lw, c='orange')
    plt.errorbar(bandwl, tdata[i]/pc.percent, tuncert[i]/pc.percent, fmt='o',
        alpha=0.8,
        ms=ms, color='b', ecolor='0.4', elinewidth=lw, capthick=lw, zorder=3)
    yran = ax.get_ylim()
    ax.set_xscale('log')
    plt.gca().xaxis.set_minor_formatter(NullFormatter())
    ax.get_xaxis().set_major_formatter(ScalarFormatter())
    ax.set_xticks(logxticks)
    ax.tick_params(labelsize=fs-1)
    plt.xlim(np.amin(wl), np.amax(wl))
    plt.text(0.01, 0.86, tnames[i], fontsize=fs-1, ha='left',
        transform=ax.transAxes)
    if i%3 == 0:
        plt.ylabel(r'$(R_{\rm p}/R_{\rm s})^2$  (%)', fontsize=fs)
for i in range(nplanets):
    ax = mp.subplotter(rect, margin, i+1+nplanets, nx=3, ny=8, ymargin=0.025)
    if ebest[i] is not None:
        plt.plot(wl, gaussf(ebest[i], sigma)/pc.ppt,lw=lw, c='orange')
    plt.errorbar(bandwl, edata[i]/pc.ppt, euncert[i]/pc.ppt, fmt='o', alpha=0.8,
        ms=ms, color='b', ecolor='0.4', elinewidth=lw, capthick=lw, zorder=3)
    yran = ax.get_ylim()
    ax.set_xscale('log')
    plt.gca().xaxis.set_minor_formatter(NullFormatter())
    ax.get_xaxis().set_major_formatter(ScalarFormatter())
    ax.set_xticks(logxticks)
    ax.tick_params(labelsize=fs-1)
    plt.xlim(np.amin(wl), np.amax(wl))
    plt.text(0.01, 0.86, enames[i], fontsize=fs-1, ha='left',
        transform=ax.transAxes)
    if i%3 == 0:
        ax.set_ylabel(r'$(F_{\rm p}/F_{\rm s})^2$ (ppt)', fontsize=fs)
    if i >= 9:
        ax.set_xlabel(r'Wavelength (um)', fontsize=fs)
axes[3].yaxis.set_major_locator(MultipleLocator(0.2))
axes[4].yaxis.set_major_locator(MultipleLocator(0.02))
plt.savefig('../plots/pyratbay-taurex_ret_spectra_comparision.pdf')
# Lower temperature+higher radius(p0) flattens IR bands
# (at constant abundances)




# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Transmission posterior plot:

hkw = {'histtype':'step', 'lw':1.5, 'density':False}
bkw1 = {'color':'crimson', 'dashes':(5,1), 'lw':1.5}
bkw2 = {'color':'green', 'dashes':(), 'lw':1.5}
fs = 10
lw = 2.0
tmin, tmax = 0.4, 0.6


width = np.array([
#   TRA1  G1214  G436  W107  H11   H26  W39  HD1  HD2  K11   55C   W76
    [ 80,  150,  180,  70,   355,  300, 110, 860, 270, 420, 2390, 1000], # T
    [0.09, 0.12, 0.06, 0.15, 0.18,  0.9, 0.2, 0.3, 0.2, 0.8,  1.0,  0.9], # R
    [1.1,  0.6,  0.7,  0.4,  2.0,  0.0, 0.7, 3.2, 0.7, 0.7,  6.0,  1.5], # H2O
    [0.7,  0.5,  0.7,  0.3,  1.7,  1.7, 0.0, 0.0, 0.7, 0.4,  0.0,  0.0], # CH4
    [0.0,  0.0,  0.0,  1.0,  0.0,  6.0, 0.0, 5.0, 2.0, 1.0,  0.0,  4.0], # CO
    [0.0,  1.7,  5.0,  3.0,  0.0,  0.0, 0.0, 7.0, 1.2, 0.0,  0.0,  0.0], # CO2
    [3.0,  0.8,  0.0,  1.6,  0.0,  2.5, 0.0, 3.5, 1.7, 0.4,  0.0,  2.5], # cl
    ])
ranges = np.tile(None, (tnpars, nplanets, 2))
igood = np.array([best is not None for best in tbest])

for p in range(tnpars):
    Tmax, Tmin = np.amax(ttrue_vals[igood,p]), np.amin(ttrue_vals[igood,p])
    if 2 <= p <= 5:
        Tmin = np.amin(ttrue_vals[igood,p][ttrue_vals[igood,p]>-12])
    if p == 6:
        Tmax = np.amax(ttrue_vals[igood,p][ttrue_vals[igood,p]<2.0])
    a = (tmax-tmin)/(Tmax-Tmin)
    b = tmax - a*Tmax
    t = a*ttrue_vals[igood,p] + b
    ranges[p,igood,0] = ttrue_vals[igood,p] - t*width[p,igood]
    ranges[p,igood,1] = ttrue_vals[igood,p] + (1-t)*width[p,igood]
    if p == 6:
        ranges[p,igood,1] = ttrue_vals[igood,p] - t*width[p,igood]
        ranges[p,igood,0] = ttrue_vals[igood,p] + (1-t)*width[p,igood]

ranges[2,width[2]==0] = -12.3, -1.0
ranges[3,width[3]==0] = -12.3, -1.0
ranges[4,width[4]==0] = -12.3, -1.0
ranges[5,width[5]==0] = -12.3, -1.0
ranges[4,3] = -12.3, -1.0
idx = 0, 2, 4, 6, 10
ranges[6,idx] = 2.1, -3.5

fig = plt.figure(10, (8.2, 10))
plt.clf()
axes = np.empty((tnpars, nplanets), object)
plt.subplots_adjust(0.07, 0.05, 0.99, 0.99, hspace=0.07, wspace=0.0)
for p in range(tnpars):
    for i in range(nplanets):
        if not igood[i]:
            continue
        ecol = 'k' if width[p,i] == 0 else 'b'
        axes[p,i] = ax = fig.add_subplot(tnpars, 12, i+1 + 12*p)
        vals, bins, h = ax.hist(tposterior[i,p], bins=25, range=None, zorder=0,
            orientation='horizontal', edgecolor=ecol, **hkw)
        vals = np.r_[0, vals, 0]
        bins = np.r_[bins[0] - (bins[1]-bins[0]), bins]
        f = si.interp1d(bins+0.5*(bins[1]-bins[0]), vals, kind='nearest')
        ax.fill_betweenx(Xpdf[i][p], 0, f(Xpdf[i][p]),
            where=PDF[i][p]>=HPDmin[i,p], facecolor='0.75', edgecolor='none',
            interpolate=False, zorder=-2)
        xtop = np.amax(vals)
        ax.axhline(ttrue_vals[i,p], xmax=0.5, **bkw2)
        ax.axhline(tbest_vals[i,p], xmax=0.5, **bkw1)
        ax.set_xlim(0, 2.2*xtop)
        ax.set_xticklabels([])
        ax.set_xticks([])
        adjust_spines(ax, ['left', 'bottom'])
        ax.tick_params(labelsize=fs-2, direction='in')
        if width[p,i] == 0 and p != tnpars-1:
            ax.yaxis.set_major_locator(MultipleLocator(3))
        if i == 0:
            ax.set_ylabel(texnames[p], fontsize=fs)
        if ranges[p,i,0] is not None:
            ax.set_ylim(ranges[p,i])
        if p == tnpars - 1:
            ax.set_xlabel(tnames[i], fontsize=fs-3)


axes[1,1].yaxis.set_major_locator(MultipleLocator(0.03))
axes[1,3].set_yticks(np.arange(10.25, axes[1,3].get_ylim()[1], 0.04))
axes[2,0].yaxis.set_major_locator(MultipleLocator(0.3))
axes[3,5].yaxis.set_major_locator(MultipleLocator(0.4))
axes[5,8].yaxis.set_major_locator(MultipleLocator(0.3))
plt.savefig('../plots/pyratbay-taurex_ret_transmission_comparision.pdf')


# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Emission posterior plot:
isort = {
    'WASP-80b':  0,
    'TrES-1b':   1,
    'WASP-6b':   2,
    'HAT-P-1b':  3,
    'WASP-43b':  4,
    'WASP-94Ab': 5,
    'WASP-17b':  6,
    'WASP-31b':  7,
    'HD149026b': 8,
    'WASP-14b':  9,
    'WASP-12b': 10,
    'WASP-121b':11,
    }

mcfiles = sorted([
    mcmc_file for mcmc_file in os.listdir('.')
    if mcmc_file.startswith('MCMC_emission')
    if mcmc_file.endswith('.npz')])

nplanets = len(mcfiles)
npars = 8
nlayers = 22
names = np.zeros(nplanets, object)
argsort = np.zeros(nplanets, int)
true_vals = np.zeros((npars, nplanets))
best_vals = np.zeros((npars, nplanets))
lo_vals = np.zeros((npars, nplanets))
hi_vals = np.zeros((npars, nplanets))
posterior = np.zeros((npars, nplanets), object)
tpost = np.zeros((nplanets, 5, nlayers))
best_pt = np.zeros((nplanets, nlayers))
planet_data = np.loadtxt('retrieval_validation_emission.txt', dtype=str)
planets = list(planet_data[:,0])
ifit = [7, 8, 6, 4, 9, 10, 11, 12]
params = np.array(planet_data[:,ifit], np.double)

for i in range(nplanets):
    mcfile = mcfiles[i]
    mcmc = np.load(mcfile)
    if 'bestp' not in mcmc:
        mcmc.close()
        continue

    logname = mcfile[0:mcfile.find('.')]
    names[i] = logname.split('_')[-1]
    post, zchain, zmask = mu.burn(mcmc)
    for p in range(npars):
        posterior[p,i] = post[:,p]
    ifree = mcmc['ifree']
    best_vals[:,i] = mcmc['bestp'][ifree]
    lo_vals[:,i] = mcmc['CRlo'][ifree]
    hi_vals[:,i] = mcmc['CRhi'][ifree]
    texnames = mcmc['texnames'][ifree]
    j = planets.index(names[i])
    true_vals[:,i] = params[j]
    argsort[i] = isort[names[i]]

    pyrat = pb.run(f'mcmc_emission_{names[i]}.cfg', init=True, no_logfile=True)
    ifree = pyrat.ret.pstep[pyrat.ret.itemp] > 0
    itemp = np.arange(np.sum(ifree))
    tpost[i] = posterior_pt(post[:,itemp], pyrat.atm.tmodel,
          pyrat.ret.params[pyrat.ret.itemp], ifree, pyrat.atm.press)
    best_pt[i] = pyrat.atm.tmodel(mcmc['bestp'][pyrat.ret.itemp])
    mcmc.close()

# Filter out unfinished runs (sort by temp):
mask = np.where(best_vals[0]!=0)[0]
ind = np.argsort(argsort[mask])
names = names[mask][ind]
tpost = tpost[mask][ind]
best_pt = best_pt[mask][ind]
true_vals = true_vals[:,mask][:,ind]
best_vals = best_vals[:,mask][:,ind]
lo_vals = lo_vals[:,mask][:,ind]
hi_vals = hi_vals[:,mask][:,ind]
posterior = posterior[:,mask][:,ind]
nplanets = len(names)

PDF, Xpdf = [], []
HPDmin = np.zeros((npars, nplanets))
quantile = 0.683
for p in range(npars):
    PDF.append([])
    Xpdf.append([])
    for i in range(nplanets):
        pdf, xpdf, HPDmin[p,i] = ms.cred_region(posterior[p,i], quantile)
        PDF[p].append(pdf)
        Xpdf[p].append(xpdf)

# TauREx temperature profiles:
tr_temp = np.zeros((nplanets, nlayers))
sys_params = np.array(planet_data[:,1:7], np.double)
for i,planet in enumerate(names):
    j = planets.index(planet)
    rstar, tstar, kmag, rplanet, mplanet, tplanet = sys_params[j]
    tm, wn_tau, depth = tas.simulate_ariel(
        rplanet*pc.rearth, mplanet*pc.mearth, params[j,2], params[j,4:],
        100, rstar, tstar, tpars=params[j,0:2], mode='eclipse')
    tr_temp[i] = tm.temperatureProfile
tr_press = tm.pressureProfile * pc.pascal/pc.bar


hkw = {'histtype':'step', 'lw':1.5, 'density':False}
bkw1 = {'color':'crimson', 'dashes':(5,1), 'lw':1.5}
bkw2 = {'color':'green', 'dashes':(), 'lw':1.5}

rect = [0.07, 0.05, 0.99, 0.38]
margin = 0.02
fs = 10
lw = 2.0
tmin, tmax = 0.4, 0.6

width = np.array([
#    W80   Tr1   W6    H1   W43  W94* W17  W31  149* W14  W12   W121*
    [150,  195,  80,   300, 110, 900, 900, 300, 400, 400, 2340, 2340], # kap
    [0.12, 0.06, 0.15, 0.9, 0.2, 0.3, 0.3, 0.2, 0.8, 0.8,  1.0,  1.0], # gam
    [0.6,  0.7,  0.4,  0.0, 0.7, 3.2, 3.2, 0.7, 0.7, 0.7,  6.0,  6.0], # bet
    [5.0,  14,   5.0,  5.5, 0.9, 4.5, 5.5, 5.5, 2.2, 2.2,  4.0,  5.0], # R
    [5.5,  4.5,  4.2,  0.0, 3.0, 4.0, 3.0, 6.0, 3.8, 5.0,  4.0,  5.3], # H2O
    [4.5,  0.0,  0.0,  0.0, 0.0, 0.0, 3.0, 0.0, 0.0, 0.0,  0.0,  0.0], # CH4
    [1.0,  1.0,  0.0,  0.0, 5.0, 1.0, 4.0, 0.0, 0.0, 1.0,  1.0,  1.0], # CO
    [1.0,  1.0,  6.0,  0.0, 3.0, 0.0, 0.0, 0.0, 1.0, 0.0,  1.0,  0.0], # CO2
    ])
ranges = np.tile(None, (npars, nplanets,2))
for p in range(npars):
    Tmax, Tmin = np.amax(true_vals[p]), np.amin(true_vals[p])
    if 2 <= p <= 5:
        Tmin = np.amin(true_vals[p][true_vals[p]>-12])
    if p == 6:
        Tmax = np.amax(true_vals[p][true_vals[p]<2.0])
    a = (tmax-tmin)/(Tmax-Tmin)
    b = tmax - a*Tmax
    t = a*true_vals[p] + b
    ranges[p,:,0] = true_vals[p] - t*width[p]
    ranges[p,:,1] = true_vals[p] + (1-t)*width[p]
    if p == 6:
        ranges[p,:,0] = true_vals[p] - t*width[p]
        ranges[p,:,1] = true_vals[p] + (1-t)*width[p]

ranges[4,:] = -7.0, -1.0
for i in range(4, npars):
    ranges[i,width[i]==0] = -12.3, -1.0
ranges[6,width[6]==1] = -12.3, -1.0
ranges[7,width[7]==1] = -12.3, -1.0


fig = plt.figure(11, (8.2, 10))
plt.clf()
axes = np.empty((npars, nplanets), object)
plt.subplots_adjust(0.07, 0.05, 0.99, 0.99, hspace=0.07, wspace=0.0)
for p in range(3,npars):
    for i in range(nplanets):
        ecol = 'k' if width[p,i] == 0 else 'b'
        axes[p,i] = ax = fig.add_subplot(npars, 12, i+1 + 12*(p-3))
        vals, bins, h = ax.hist(posterior[p,i], bins=25, range=None, zorder=0,
            orientation='horizontal', edgecolor=ecol, **hkw)
        vals = np.r_[0, vals, 0]
        bins = np.r_[bins[0] - (bins[1]-bins[0]), bins]
        f = si.interp1d(bins+0.5*(bins[1]-bins[0]), vals, kind='nearest')
        ax.fill_betweenx(Xpdf[p][i], 0, f(Xpdf[p][i]),
            where=PDF[p][i]>=HPDmin[p,i], facecolor='0.75', edgecolor='none',
            interpolate=False, zorder=-2)
        xtop = np.amax(vals)
        ax.axhline(true_vals[p,i], xmax=0.5, **bkw2)
        ax.axhline(best_vals[p,i], xmax=0.5, **bkw1)
        ax.set_xlim(0, 2.2*xtop)
        ax.set_xticklabels([])
        ax.set_xticks([])
        adjust_spines(ax, ['left', 'bottom'])
        ax.tick_params(labelsize=fs-2, direction='in')
        if width[p,i] == 0 and p != npars-1:
            ax.yaxis.set_major_locator(MultipleLocator(3))
        if i == 0:
            ax.set_ylabel(texnames[p], fontsize=fs)
        if ranges[p,i,0] is not None:
            ax.set_ylim(ranges[p,i])
        if p == npars - 1:
            ax.set_xlabel(names[i], fontsize=fs-3)
for i in range(nplanets):
    ax = mp.subplotter(rect, margin, i+1, nx=6, ny=2)
    ax.fill_betweenx(pyrat.atm.press/pc.bar, tpost[i,2], tpost[i,4],
        facecolor='royalblue', edgecolor='none', alpha=0.6)
    ax.fill_betweenx(pyrat.atm.press/pc.bar, tpost[i,1], tpost[i,3],
        facecolor='royalblue', edgecolor='none', alpha=0.8)
    plt.plot(tpost[i,0], pyrat.atm.press/pc.bar, 'navy', lw=1.5, label='Median')
    ax.plot(tr_temp[i], tr_press, c='limegreen', lw=1.5)
    plt.plot(best_pt[i], pyrat.atm.press/pc.bar, **bkw1)
    ax.set_ylim(100, 1e-5)
    ax.set_xlim(200, 3150)
    ax.set_yscale('log')
    if i % 6:
        ax.set_yticklabels([])
    ax.tick_params(labelsize=fs-2, direction='in')
    plt.text(0.98,0.92,names[i],fontsize=fs-3,ha='right',transform=ax.transAxes)
    if i == 0:
        plt.text(-0.42, -0.4, 'Pressure (bar)', fontsize=fs,
            rotation='vertical', transform=ax.transAxes)
        plt.text(2.8, -1.35, 'Temperature (K)', fontsize=fs,
            transform=ax.transAxes)
plt.savefig('../plots/pyratbay-taurex_ret_emission_comparision.pdf')

