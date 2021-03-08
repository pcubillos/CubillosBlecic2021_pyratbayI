import os

import numpy as np
import matplotlib.pyplot as plt

import pyratbay as pb
import pyratbay.atmosphere as pa
import pyratbay.io as io
import pyratbay.tools as pt
import mc3


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
    Zsun = 0.0134  # Asplund et al. 2009
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
    # Only those with ADI > 3:
    cfiles = [
        'run_HATP-18b_tsiaras/mcmc_HATP-18b_tsiaras_w00000-0c.cfg',
        'run_HATP-32b_tsiaras/mcmc_HATP-32b_tsiaras_w00000-0c.cfg',
        'run_HATP-38b_bruno/mcmc_HATP-38b_bruno_w00000-0c.cfg',
        'run_HATP-41b_tsiaras/mcmc_HATP-41b_tsiaras_w00000-0c.cfg',
        'run_K2-018b_benneke/mcmc_K2-018b_benneke_w00000-0c.cfg',

        'run_WASP-043b_stevenson/mcmc_WASP-043b_stevenson_w00000-0c.cfg',
        'run_WASP-063b_kilpatrick/mcmc_WASP-063b_kilpatrick_w000h0-0c.cfg',
        'run_WASP-069b_tsiaras/mcmc_WASP-069b_tsiaras_w00000-0c.cfg',
        'run_WASP-103b_kreidberg/mcmc_WASP-103b_kreidberg_w00000-0c.cfg',
        'run_XO-1b_tsiaras/mcmc_XO-1b_tsiaras_w00000-0c.cfg',
        ]
    ntargets = len(cfiles)

    names = np.array([None for _ in range(ntargets)])
    posteriors = [None for _ in range(ntargets)]
    teq = np.zeros(ntargets)

    for i in range(ntargets):
        folder, cfile = os.path.split(cfiles[i])
        _, planet, dset = folder.split('_')
        name = f'{planet} {dset.capitalize()}'
        names[i] = name.replace('-0','-')
        with pt.cd(folder):
            pyrat = pb.run(cfile, init=True, no_logfile=True)
            teq[i], _ = pa.equilibrium_temp(
                pyrat.phy.tstar, pyrat.phy.rstar, pyrat.phy.smaxis)
            with np.load(pyrat.ret.mcmcfile) as mcmc:
                posteriors[i], zchain, zmask = mc3.utils.burn(mcmc)


    tposteriors = np.array([post[:,0] for post in posteriors])
    ci1 = np.zeros((2,ntargets))
    for i in range(ntargets):
        pdf, xpdf, HPDmin = mc3.stats.cred_region(tposteriors[i])
        ci1[:,i] = -np.amin(xpdf[pdf>HPDmin]), np.amax(xpdf[pdf>HPDmin])


    # Print some stats:
    decimals = [0, 3, 2, 2, 2, 2]
    split = [None, None, -4, None, None, None, None, -2, -4, None]

    table = []
    for j,posterior in enumerate(posteriors):
        texts = []
        for i in range(len(posterior.T)):
            dec = decimals[i]
            text = ''
            pdf, xpdf, hpd_min = mc3.stats.cred_region(posterior[:,i])
            if i == 2 and split[j] is not None:
                x = xpdf[xpdf < split[j]]
                p = pdf[xpdf < split[j]]
                mode = x[np.argmax(p)]
                lo =   mode - np.amin(x[p>hpd_min])
                hi = -(mode - np.amax(x[p>hpd_min]))
                text = f'${mode:.{dec}f}_{{-{lo:.{dec}f}}}^{{+{hi:.{dec}f}}}$ and '
                pdf = pdf[xpdf > split[j]]
                xpdf = xpdf[xpdf > split[j]]
            mode = xpdf[np.argmax(pdf)]
            lo =   mode - np.amin(xpdf[pdf>hpd_min])
            hi = -(mode - np.amax(xpdf[pdf>hpd_min]))
            text += f'${mode:.{dec}f}_{{-{lo:.{dec}f}}}^{{+{hi:.{dec}f}}}$'
            texts.append(text)
        text = "  &  ".join(texts)
        table.append(text)
    table_text = " \\\\\n".join(table)
    print(table_text)

    # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # Temperatures:
    lw = 2.0
    fs = 14
    y = np.zeros(ntargets)

    plt.figure(22, (6.5,4.75))
    plt.clf()
    plt.subplots_adjust(0.12, 0.1, 0.99, 0.99)
    ax = plt.subplot(111)
    plt.errorbar(
        teq, y, yerr=ci1, capsize=5.0,
        capthick=lw, lw=lw, ls='none', color='orangered', label='This work')
    plt.plot(teq, teq, "o", color='mediumblue')
    plt.ylim(ymin=0)
    ax.set_xlim(200, 2560)
    ax.tick_params(labelsize=fs-2, direction='in', which='both')
    for xt, name in zip(teq, names):
        name = name.split()[0]
        plt.text(xt-61, xt+45, name.split()[0],
            va='bottom', rotation=90, fontsize=fs-4)
    ax.set_xlabel('Equilibrium temperature (K)', fontsize=fs)
    ax.set_ylabel('Retrieved temperature (K)', fontsize=fs)
    plt.savefig('plots/WFC3_sample_temperature.pdf')

    # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # H2O abundances:
    idx_H2O = 2
    H2O_posts = [post[:,idx_H2O] for post in posteriors]

    nx = 1
    ny = ntargets
    rect = [0.1, 0.06, 0.95, 0.99]
    margin = 0.0
    ranges = [(-12, 0) for _ in range(ntargets)]
    pnames = ['' for _ in range(ntargets)]
    fs = 13

    plt.figure(32, (6,8))
    plt.clf()
    axes = [
        mc3.plots.subplotter(rect, margin, i+1, nx, ny)
        for i in range(ntargets)]

    for i in range(ntargets):
        mc3.plots.histogram(
            H2O_posts[i], pnames=pnames[i:i+1], quantile=0.683,
            ranges=ranges[i:i+1], nbins=60,
            axes=axes[i:i+1], lw=2.0, fs=fs-1, theme='blue')
        plt.text(0.02, 0.75, names[i], transform=axes[i].transAxes)
        if i < ntargets - 1:
            axes[i].set_xticklabels([])
        axes[i].set_yticks([])
    axes[i].set_xlabel('$\\log_{10}(X_{\\rm H2O})$', fontsize=fs)
    axes[4].set_ylabel('Marginal posterior density' + 15*' ', fontsize=fs)
    plt.setp(axes[i].xaxis.get_majorticklabels(), rotation=0)
    plt.savefig('plots/WFC3_sample_H2O_abundance.pdf')


if __name__ == '__main__':
    main()
