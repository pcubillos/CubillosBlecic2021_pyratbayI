import os

import numpy as np
import matplotlib.pyplot as plt
plt.ioff()

import pyratbay as pb
import pyratbay.atmosphere as pa
import pyratbay.io as io
import pyratbay.tools as pt
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
    # Mass fractions (Asplund et al. 2009):
    Zsun = 0.0134
    Xsun = 0.7381

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


def main():
    cfile = 'run_HATP-41b_tsiaras/mcmc_HATP-41b_tsiaras_w00000-mc.cfg'
    folder, cfile = os.path.split(cfile)
    _, planet, dset = folder.split('_')
    with pt.cd(folder):
        pyrat = pb.run(cfile, init=True, no_logfile=True)

    with np.load(pyrat.ret.mcmcfile) as mcmc:
        post, zchain, zmask = mc3.utils.burn(mcmc)
        texnames = mcmc['texnames']
        bestp = mcmc['bestp']
        imol = np.in1d(mcmc['ifree'], pyrat.ret.imol)

    posteriors = pyrat.ret.posterior = post
    u, uind, uinv = np.unique(
        posteriors[:,0], return_index=True, return_inverse=True)
    nunique = np.size(u)
    # Get mean molecular mass:
    ZX_ratio = np.zeros(nunique, np.double)
    params = pyrat.ret.params
    for k in range(nunique):
        params[imol] = posteriors[uind[k],imol]
        ZX_ratio[k] = ZX(params, pyrat)
    iN2 = pyrat.ret.pnames.index('log(N2)')
    # Metal mass fraction relative to solar values, i.e., [Z/X]:
    posteriors[:,iN2] = ZX_ratio[uinv]
    bestp[iN2] = ZX(bestp, pyrat)
    texnames[iN2] = r"$[Z/X]$"

    posteriors = roll(posteriors, indices=pyrat.ret.imol[:-1])
    bestp = roll(bestp, indices=pyrat.ret.imol[:-1])
    texnames = roll(texnames, indices=pyrat.ret.imol[:-1])

    texnames[4] = texnames[4].replace('f', 'X')
    rect = [0.105, 0.11, 0.99, 0.99]
    ranges = [
        (100, 3000),
        (1.49, 1.68),
        (-4.5, 2.6),
        (-5.5, 2.5),
        (-7., -0.45)]


    plt.figure(316, (6, 6))
    plt.clf()
    axes, cax = mc3.plots.pairwise(
        posteriors, pnames=texnames, ranges=ranges,
        fignum=8, rect=rect, nbins=20)
    ax = axes[1][2]
    plt.text(0.15, 0.125, 'cloudy', rotation=45, transform=ax.transAxes)
    plt.text(0.70, 0.25, 'heavy', rotation=45, transform=ax.transAxes)
    plt.text(0.75, 0.67, 'deplet', rotation=90, transform=ax.transAxes)
    ax = axes[3][3]
    plt.text(0.05, 0.30, 'cloudy', rotation=305, transform=ax.transAxes)
    plt.text(0.55, 0.75, 'heavy', rotation=0, transform=ax.transAxes)
    plt.text(0.65, 0.10, 'deplet', rotation=0, transform=ax.transAxes)

    plt.savefig(f'plots/WFC3_{planet}_{dset}_correlation.pdf')
    plt.savefig(f'plots/WFC3_{planet}_{dset}_correlation.png', dpi=300)


if __name__ == '__main__':
    main()
