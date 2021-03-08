import os

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from matplotlib.ticker import MultipleLocator, NullFormatter, ScalarFormatter
from scipy.ndimage.filters import gaussian_filter1d as gaussf

import pyratbay.atmosphere as pa
import pyratbay.constants as pc
import pyratbay.io as io
import mc3


def roll(array, indices):
    """roll rows in indices to the end of the array"""
    first = indices[0]
    last = np.shape(array)[-1] - 1
    nroll = last - indices[-1]

    if np.ndim(array) == 1:
        array[first:] = np.roll(array[first:], nroll, axis=-1)
    if np.ndim(array) == 2:
        array[:,first:] = np.roll(array[:,first:], nroll, axis=-1)
    return array


cfiles = [
    'run_GJ-0436b_tsiaras/mcmc_GJ-0436b_tsiaras_w00000-mc.cfg',
    'run_HATP-03b_tsiaras/mcmc_HATP-03b_tsiaras_w00000-mc.cfg',
    'run_HATP-17b_tsiaras/mcmc_HATP-17b_tsiaras_w00000-mc.cfg',
    'run_HATP-18b_tsiaras/mcmc_HATP-18b_tsiaras_w00000-mc.cfg',
    'run_HATP-32b_damiano/mcmc_HATP-32b_damiano_w00000-mc.cfg',
    'run_HATP-32b_tsiaras/mcmc_HATP-32b_tsiaras_w00000-mc.cfg',
    'run_HATP-38b_bruno/mcmc_HATP-38b_bruno_w00000-mc.cfg',

    'run_HATP-38b_tsiaras/mcmc_HATP-38b_tsiaras_w00000-mc.cfg',
    'run_HATP-41b_tsiaras/mcmc_HATP-41b_tsiaras_w00000-mc.cfg',
    'run_HD-149026b_tsiaras/mcmc_HD-149026b_tsiaras_w00000-mc.cfg',
    'run_K2-018b_benneke/mcmc_K2-018b_benneke_wmdm00-mc.cfg',
    'run_K2-018b_tsiaras/mcmc_K2-018b_tsiaras_w00000-mc.cfg',
    'run_WASP-029b_tsiaras/mcmc_WASP-029b_tsiaras_w00000-mc.cfg',
    'run_WASP-043b_stevenson/mcmc_WASP-043b_stevenson_wmdm00-mc.cfg',
    'run_WASP-043b_tsiaras/mcmc_WASP-043b_tsiaras_w00000-mc.cfg',

    'run_WASP-063b_kilpatrick/mcmc_WASP-063b_kilpatrick_w000h0-mc.cfg',
    'run_WASP-063b_tsiaras/mcmc_WASP-063b_tsiaras_w000h0-mc.cfg',
    'run_WASP-067b_bruno/mcmc_WASP-067b_bruno_w00000-mc.cfg',
    'run_WASP-067b_tsiaras/mcmc_WASP-067b_tsiaras_w00000-mc.cfg',
    'run_WASP-069b_tsiaras/mcmc_WASP-069b_tsiaras_w00000-mc.cfg',
    'run_WASP-074b_tsiaras/mcmc_WASP-074b_tsiaras_w00000-mc.cfg',
    'run_WASP-080b_tsiaras/mcmc_WASP-080b_tsiaras_w00000-mc.cfg',
    'run_WASP-101b_tsiaras/mcmc_WASP-101b_tsiaras_w00000-mc.cfg',

    'run_WASP-101b_wakeford/mcmc_WASP-101b_wakeford_w00000-mc.cfg',
    'run_WASP-103b_kreidberg/mcmc_WASP-103b_kreidberg_wmdm00-mc.cfg',
    'run_XO-1b_tsiaras/mcmc_XO-1b_tsiaras_w00000-mc.cfg',
    ]
ntargets = len(cfiles)


data = [None for _ in range(ntargets)]
uncert = [None for _ in range(ntargets)]
wl = [None for _ in range(ntargets)]
best_spec = [None for _ in range(ntargets)]
band_wl = [None for _ in range(ntargets)]
band_model = [None for _ in range(ntargets)]
names = [None for _ in range(ntargets)]
bestp = [None for _ in range(ntargets)]
texnames = [None for _ in range(ntargets)]
posteriors = [None for _ in range(ntargets)]
low1 = [None for _ in range(ntargets)]
low2 = [None for _ in range(ntargets)]
high1 = [None for _ in range(ntargets)]
high2 = [None for _ in range(ntargets)]
temp_ranges = [None for _ in range(ntargets)]

# Mass fractions (Asplund et al. 2009):
# X + Y + Z = mu, with:
# X = N_H  * m_H
# Y = N_He * m_He
# Z = sum(N_Z*m_z) = mu - X - Y
Zsun = 0.0134
Xsun = 0.7381

# Indices for species in atmosphere:
units, specs, press, temp, q, rad = io.read_atm(
    'run_setup/isothermal_1500K_uniform.atm')
iHe  = np.where(specs == 'He')[0][0]
iH2  = np.where(specs == 'H2')[0][0]
iH   = np.where(specs == 'H')[0][0]
iH2O = np.where(specs == 'H2O')[0][0]
iCH4 = np.where(specs == 'CH4')[0][0]
iHCN = np.where(specs == 'HCN')[0][0]

def ZX(params, pyrat):
    """Compute metals mass fraction relative to solar."""
    q2 = pa.qscale(
        pyrat.atm.qbase, pyrat.mol.name, pyrat.atm.molmodel,
        pyrat.atm.molfree, params[pyrat.ret.imol],
        pyrat.atm.bulk, iscale=pyrat.atm.ifree, ibulk=pyrat.atm.ibulk,
        bratio=pyrat.atm.bulkratio, invsrat=pyrat.atm.invsrat)[0]
    mu = pa.mean_weight(q2, pyrat.mol.name, mass=pyrat.mol.mass)[0]
    Y = pyrat.mol.mass[iHe] * q2[iHe]
    X = pyrat.mol.mass[iH]*(
           2*q2[iH2] + q2[iH] + 2*q2[iH2O] + q2[iHCN] + 4*q2[iCH4])
    Z = mu - X - Y
    return np.log10(Z/X / (Zsun/Xsun))


for i in range(ntargets):
    folder, cfile = os.path.split(cfiles[i])
    print(folder)
    _, planet, dset = folder.split('_')
    names[i] = f'{planet} {dset.capitalize()}'
    pickle_file = [f for f in os.listdir(folder) if f.endswith('pickle')][0]
    pyrat = io.load_pyrat(f'{folder}/{pickle_file}')
    with np.load(pyrat.ret.mcmcfile) as mcmc:
        post, zchain, zmask = mc3.utils.burn(mcmc)
        band_model[i] = mcmc['best_model']
        bestp[i] = mcmc['bestp']
        texnames[i] = mcmc['texnames']
        ifree = mcmc['ifree']

    wl[i] = 1.0/(pyrat.spec.wn*pc.um)
    band_wl[i] = 1.0/(pyrat.obs.bandwn*pc.um)
    data[i] = pyrat.obs.data
    uncert[i] = pyrat.obs.uncert
    best_spec[i], _ = pyrat.eval(bestp[i])
    temp_ranges[i] = (pyrat.ret.pmin[0], pyrat.ret.pmax[0])
    posteriors[i] = pyrat.ret.posterior = post
    low1[i]  = pyrat.ret.spec_low1
    low2[i]  = pyrat.ret.spec_low2
    high1[i] = pyrat.ret.spec_high1
    high2[i] = pyrat.ret.spec_high2

    u, uind, uinv = np.unique(post[:,0], return_index=True, return_inverse=True)
    nunique = np.size(u)
    # Get mean molecular mass:
    ZX_ratio = np.zeros(nunique, np.double)
    params = pyrat.ret.params
    imol = np.in1d(ifree, pyrat.ret.imol)
    for k in range(nunique):
        params[imol] = post[uind[k],imol]
        ZX_ratio[k] = ZX(params, pyrat)
    iN2 = pyrat.ret.pnames.index('log(N2)')
    # Metal mass fraction relative to solar values, i.e., [Z/X]:
    posteriors[i][:,iN2] = ZX_ratio[uinv]
    bestp[i][iN2] = ZX(bestp[i], pyrat)
    texnames[i][iN2] = r"$[\rm Z/\rm X]$"

    posteriors[i] = roll(posteriors[i], indices=pyrat.ret.imol[:-1])
    bestp[i] = roll(bestp[i], indices=pyrat.ret.imol[:-1])
    texnames[i] = roll(texnames[i], indices=pyrat.ret.imol[:-1])


themes = {
    '$\\log_{10}(f_{\\rm H2O})$': 'blue',
    '$\\log_{10}(f_{\\rm CO})$': 'red',
    '$\\log_{10}(f_{\\rm CO2})$': 'black',
    '$\\log_{10}(f_{\\rm CH4})$': 'green',
    '$\\log_{10}(f_{\\rm HCN})$': 'orange',
    }

line1 = mlines.Line2D(
    [], [], color='blue', linewidth=2.0, label=r'$X_{\rm i} = X_{\rm H2O}$')
line2 = mlines.Line2D(
    [], [], color='red',  linewidth=2.0, label=r'$X_{\rm i} = X_{\rm CO}$')
line3 = mlines.Line2D(
    [], [], color='black', linewidth=2.0, label=r'$X_{\rm i} = X_{\rm CO2}$')
line4 = mlines.Line2D(
    [], [], color='green', linewidth=2.0, label=r'$X_{\rm i} = X_{\rm CH4}$')
line5 = mlines.Line2D(
    [], [], color='orange', linewidth=2.0, label=r'$X_{\rm i} = X_{\rm HCN}$')


rect  = [0.075, 0.045, 0.985, 0.95]
rect2 = [0.36,  0.045, 0.985, 0.95]
margin1 = 0.05
margin2 = 0.01
ymargin = 0.035
sigma = 25.0
fs = 9.0
lw = 1.0
ms = 4.0
nhist = 5
logxticks = [1.0, 1.4, 2.0, 3.0, 5.0]
ranges = [(100, 3000), None, (-6,3), (-6,2), (-12,0)]

sticks = [None] * ntargets
sticks[ 3] = 0.03
sticks[ 9] = 0.005
sticks[13] = 0.02

new_page = np.array([0, 7, 15, 23, 26])
ny = 8

for i in range(ntargets):
    page = sum(i >= new_page) - 1
    pos = i - new_page[page]
    if i in new_page:
        fig = plt.figure(20+page, (8.2, 10.8))
        plt.clf()
    # The spectra:
    ax1 = mc3.plots.subplotter(
        rect, margin1, 1+3*pos, nx=3, ny=ny, ymargin=ymargin)
    if low2[i] is not None:
        ax1.fill_between(
            wl[i], gaussf(low2[i], sigma)/pc.percent,
            gaussf(high2[i],sigma)/pc.percent,
            facecolor="gold", edgecolor="none", alpha=1.0)
        ax1.fill_between(
            wl[i], gaussf(low1[i], sigma)/pc.percent,
            gaussf(high1[i],sigma)/pc.percent,
            facecolor="orange", edgecolor="none", alpha=1.0)
    ax1.plot(
        wl[i], gaussf(best_spec[i], sigma)/pc.percent,
        lw=lw, c='red', label='Bestfit model')
    ax1.plot(
        band_wl[i], band_model[i]/pc.percent, "o",
        color='orangered', mew=0.25, mec="k", ms=ms)
    ax1.errorbar(
        band_wl[i], data[i]/pc.percent, uncert[i]/pc.percent,
        fmt='o', alpha=0.7, ms=ms, mew=0.25, color='b',
        elinewidth=lw, capthick=lw, zorder=3, label='Data')
    ax1.tick_params(labelsize=fs-1, direction='in', which='both')
    if ax1.get_xlim()[1] > 3:
        ax1.set_xscale('log')
        plt.gca().xaxis.set_minor_formatter(NullFormatter())
        ax1.get_xaxis().set_major_formatter(ScalarFormatter())
        ax1.set_xticks(logxticks)
        ax1.set_xlim(1.0, 5.5)
    else:
        ax1.set_xlim(1.0, 2.0)
    ax1.text(
        0.997, 0.89, names[i], fontsize=fs-1, ha='right',
        transform=ax1.transAxes)
    if sticks[i] is not None:
        ax1.yaxis.set_major_locator(MultipleLocator(sticks[i]))
    plt.ylabel(r'$(R_{\rm p}/R_{\rm s})^2$ (%)', fontsize=fs)
    # The posteriors:
    axes = [
        mc3.plots.subplotter(
            rect2, margin2, nhist*pos+k+1, nx=nhist, ny=ny, ymargin=ymargin)
        for k in range(nhist)]
    ranges[0] = temp_ranges[i]
    mc3.plots.histogram(
        posteriors[i][:,:nhist-1], pnames=texnames[i], bestp=bestp[i],
        ranges=ranges, quantile=0.683, axes=axes, fs=fs, lw=1.2,
        yscale=None, theme='blue')
    for j in range(nhist-1, len(bestp[i])):
        mc3.plots.histogram(
            posteriors[i][:,j], pnames=['$\\log_{10}(X_{\\rm i})$'],
            bestp=bestp[i][j:j+1], ranges=[ranges[-1]], quantile=0.683,
            fs=fs, lw=1.2, axes=[axes[-1]], yscale=None,
            theme=themes[texnames[i][j]])
    for ax in axes:
        ax.tick_params(labelsize=fs-2, direction='in', which='both')
    if i in new_page:
        ax1.legend(
            bbox_to_anchor=(1.0, 1.3), ncol=2, fontsize=fs,
            loc='upper right', borderaxespad=0.0)
        axes[-1].legend(
            handles=[line1, line2, line3, line4, line5],
            bbox_to_anchor=(1.0, 1.55), ncol=3, fontsize=fs, borderpad=0.3,
            loc='upper right', borderaxespad=0., handlelength=1.5)
    if i+1 in new_page:
        ax1.set_xlabel('Wavelength (um)', fontsize=fs)
        plt.savefig(f'plots/WFC3_sample_spectra_histograms_0{page+1}.pdf')
        plt.savefig(
            f'plots/WFC3_sample_spectra_histograms_0{page+1}.png', dpi=300)
    else:
        for ax in axes:
            ax.set_xlabel('')
