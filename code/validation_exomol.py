import sys
import pickle

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker

import pyratbay as pb
import pyratbay.constants as pc


with open('../inputs/taurex/H2O_pokazatel.R10000.TauREx.pickle', "rb") as f:
    xs_H2O = pickle.load(f)['xsecarr']
with open('../inputs/taurex/CO_Li2015.R10000.TauREx.pickle', "rb") as f:
    xs_CO  = pickle.load(f)['xsecarr']
with open('../inputs/taurex/CO2.R10000.TauREx.pickle', "rb") as f:
    d = pickle.load(f)
    xs_CO2 = d['xsecarr']
    wl_exomol = 1.0/(d['wno']*pc.um)
    press = d['p']
    temp  = d['t']

itemp = 9
ipressures = [3, 12, 21]
fs = 11
xlim = 1.0, 10.0
xticks = [1, 2, 3, 4, 6, 10]

# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# H2O:
pyrat = pb.run('H2O_spectrum.cfg')
imol = list(pyrat.mol.name).index('H2O')
wl_pyrat = 1.0/(pyrat.spec.wn*pc.um)


plt.figure(20, (8, 7))
plt.clf()
ax = plt.subplot(111)
plt.subplots_adjust(0.1, 0.08, 0.97, 0.98, wspace=0.2, hspace=0.15)
for i,ipress in enumerate(ipressures):
    for j in [0,1]:
        ax = plt.subplot(3, 2, 2*i+1+j)
        plt.plot(wl_pyrat, pyrat.ex.ec[ipress]/pyrat.atm.d[ipress,imol], c='b',
            label='pyratbay')
        plt.plot(wl_exomol, xs_H2O[ipress,itemp], alpha=0.7, c='orange',
            label='exomol')
        plt.yscale('log')
        if j == 0:
            plt.text(0.03, 0.9, f'T={temp[itemp]} K', fontsize=fs,
                     transform=ax.transAxes)
            plt.text(0.03, 0.8, f'p={press[ipress]:.1e} bar', fontsize=fs,
                     transform=ax.transAxes)
            if i == 1:
                plt.ylabel(r'Extinction coefficient (cm$^{2}$ molec$^{-1}$)')
        if i == 0 and j == 0:
            plt.legend(loc='upper right')
        if i == 2:
            plt.xlabel('Wavelength (um)')
        if j == 0:
            plt.xscale('log')
            plt.xlim(xlim)
            plt.ylim(1e-27, 1e-14)
            ax.get_xaxis().set_major_formatter(ticker.ScalarFormatter())
            ax.set_xticks(xticks)
            plt.gca().xaxis.set_minor_formatter(ticker.NullFormatter())
        else:
            plt.xticks(np.arange(1.4, 1.44, 0.01))
            plt.xlim(1.4, 1.43)
            plt.ylim(1e-25, 1e-18)
plt.savefig('../plots/H2O_pokazatel_comparison_pyratbay-exomol.pdf')


# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# CO:
pyrat = pb.run('CO_spectrum.cfg')
imol = list(pyrat.mol.name).index('CO')
# Consider only the first isotope:
pyrat.iso.ratio[1:] = 1e-30
pyrat.run()


plt.figure(30, (8, 7))
plt.clf()
ax = plt.subplot(111)
plt.subplots_adjust(0.1, 0.08, 0.97, 0.98, wspace=0.2, hspace=0.15)
for i,ipress in enumerate(ipressures):
    for j in [0,1]:
        ax = plt.subplot(3,2, 2*i+1+j)
        plt.plot(wl_pyrat, pyrat.ex.ec[ipress]/pyrat.atm.d[ipress,imol], c='b',
            label='pyratbay')
        plt.plot(wl_exomol, xs_CO[ipress,itemp], alpha=0.7, c='orange',
            label='exomol')
        plt.yscale('log')
        if j == 0:
            plt.text(0.03, 0.9, f'T={temp[itemp]} K', fontsize=fs,
                     transform=ax.transAxes)
            plt.text(0.03, 0.8, f'p={press[ipress]:.1e} bar', fontsize=fs,
                     transform=ax.transAxes)
            if i == 1:
                plt.ylabel(r'Extinction coefficient (cm$^{2}$ molec$^{-1}$)')
        if i == 0 and j == 0:
            plt.legend(loc='upper right')
        if i == 2:
            plt.xlabel('Wavelength (um)')
        if j == 0:
            plt.xscale('log')
            plt.xlim(xlim)
            plt.ylim(1e-45, 1e-10)
            ax.get_xaxis().set_major_formatter(ticker.ScalarFormatter())
            ax.set_xticks(xticks)
            plt.gca().xaxis.set_minor_formatter(ticker.NullFormatter())
        else:
            plt.xlim(1.177, 1.23)
            plt.ylim(1e-34, 3e-24)
plt.savefig('../plots/CO_HITEMP_comparison_pyratbay-exomol.pdf')


# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# CO2:
pyrat = pb.run('CO2_spectrum.cfg')
imol = list(pyrat.mol.name).index('CO2')
# Consider only the first isotope:
pyrat.iso.ratio[1:] = 1e-30
pyrat.run()


plt.figure(40, (8, 7))
plt.clf()
ax = plt.subplot(111)
plt.subplots_adjust(0.1, 0.08, 0.97, 0.98, wspace=0.2, hspace=0.15)
for i,ipress in enumerate(ipressures):
    for j in [0,1]:
        ax = plt.subplot(3, 2, 2*i+1+j)
        plt.plot(wl_pyrat, pyrat.ex.ec[ipress]/pyrat.atm.d[ipress,imol], c='b',
            label='pyratbay')
        plt.plot(wl_exomol, xs_CO2[ipress,itemp], alpha=0.7, c='orange',
            label='exomol')
        plt.yscale('log')
        if j == 0:
            plt.text(0.03, 0.9, f'T={temp[itemp]} K', fontsize=fs,
                     transform=ax.transAxes)
            plt.text(0.03, 0.8, f'p={press[ipress]:.1e} bar', fontsize=fs,
                     transform=ax.transAxes)
            if i == 1:
                plt.ylabel(r'Extinction coefficient (cm$^{2}$ molec$^{-1}$)')
        if i == 0 and j == 0:
            plt.legend(loc='lower right')
        if i == 2:
            plt.xlabel('Wavelength (um)')
        if j == 0:
            plt.xscale('log')
            plt.xlim(xlim)
            plt.ylim(1e-45, 1e-10)
            ax.get_xaxis().set_major_formatter(ticker.ScalarFormatter())
            ax.set_xticks(xticks)
            plt.gca().xaxis.set_minor_formatter(ticker.NullFormatter())
        else:
            plt.xlim(2.6, 2.7)
            plt.ylim(1e-25, 1e-18)
plt.savefig('../plots/CO2_HITEMP_comparison_pyratbay-exomol.pdf')
