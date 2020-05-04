import os
import pickle

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker

import pyratbay as pb
import pyratbay.io as io
import pyratbay.constants as pc
import pyratbay.atmosphere as pa


HOME = os.path.abspath('..')

# Read ExoMol cross-section opacity data:
xs_H2O, press_exomol, temp_exomol, wn_exomol, _ = io.import_exomol_xs(
    f'{HOME}/inputs/taurex/H2O_pokazatel.R10000.TauREx.pickle')
xs_CO2 = io.import_exomol_xs(
    f'{HOME}/inputs/taurex/CO2.R10000.TauREx.pickle', read_all=False)
xs_CO  = io.import_exomol_xs(
    f'{HOME}/inputs/taurex/CO_Li2015.R10000.TauREx.pickle', read_all=False)

wl_exomol = 1.0/(wn_exomol*pc.um)
xs_H2O = xs_H2O[3::3]
xs_CO  = xs_CO [3::3]
xs_CO2 = xs_CO2[3::3]
press_exomol = press_exomol[3::3]

# Create an atmospheric file at same selected pressures:
nlayers = 7
press = pa.pressure('1e-4 bar', '1e2 bar', nlayers)
temp  = np.tile(1000.0, nlayers)
species    = 'H2    He     H2O   CH4   CO    CO2'.split()
abundances = [0.86, 0.14, 1e-4, 1e-4, 1e-4, 1e-4]
atmfile = 'exomol_opacity_benchmark.atm'

q_atm = pa.uniform(press, temp, species, abundances, atmfile=atmfile)

# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
pyrat = pb.run('H2O_spectrum.cfg', init=True)
wl_pyrat = 1.0/(pyrat.spec.wn*pc.um)
pyrat.run()
ec_H2O = pyrat.ex.ec

pyrat = pb.run('CO_spectrum.cfg', init=True)
pyrat.iso.ratio[1:] = 1e-30  # Consider only the first isotope
pyrat.run()
ec_CO = pyrat.ex.ec

pyrat = pb.run('CO2_spectrum.cfg', init=True)
pyrat.iso.ratio[1:] = 1e-30  # Consider only the first isotope
pyrat.run()
ec_CO2 = pyrat.ex.ec

# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
ipress = 3
itemp = 9

iH2O = species.index('H2O')
iCO  = species.index('CO')
iCO2 = species.index('CO2')
ec = [ec_H2O[ipress]/pyrat.atm.d[ipress,iH2O],
      ec_CO [ipress]/pyrat.atm.d[ipress,iCO],
      ec_CO2[ipress]/pyrat.atm.d[ipress,iCO2]]
xs = xs_H2O[ipress,itemp], xs_CO[ipress,itemp], xs_CO2[ipress,itemp]

fs = 11
xlim = 1.0, 10.0
xticks = [1, 2, 3, 4, 6, 10]
xran = [(1.4, 1.432), (1.177, 1.232), (2.6, 2.73)]
yran1 = [(3e-27, 5e-14), (1e-48, 1e-15), (1e-33, 1e-15)]
yran2 = [(1e-25, 1e-18), (8e-34, 3e-24), (1e-25, 1e-18)]
labels = ['H2O/POKAZATEL', 'CO/Li', 'CO2/HITEMP']

plt.figure(50, (8, 7))
plt.clf()
ax = plt.subplot(111)
plt.subplots_adjust(0.1, 0.08, 0.97, 0.98, wspace=0.2, hspace=0.15)
for i in range(len(xs)):
    for j in [0,1]:
        ax = plt.subplot(3, 2, 2*i+1+j)
        plt.plot(wl_pyrat, ec[i], c='b', label='pyratbay')
        plt.plot(wl_exomol, xs[i], c='orange', label='exomol', alpha=0.7)
        plt.yscale('log')
        if i == 0 and j == 0:
            plt.legend(loc='upper right')
        if i == 1 and j == 0:
            plt.ylabel(r'Opacity (cm$^{2}$ molec$^{-1}$)', fontsize=fs)
        if i == 0 and  j == 1:
            plt.text(0.03, 0.9, f'$T$={temp_exomol[itemp]:.0f} K',
                     fontsize=fs, transform=ax.transAxes)
            plt.text(0.03, 0.8, f'$p$={press_exomol[ipress]/pc.bar:.1f} bar',
                     fontsize=fs, transform=ax.transAxes)
            plt.xticks(np.arange(1.4, 1.44, 0.01))
        if i == 2:
            plt.xlabel('Wavelength (um)', fontsize=fs)
        if j == 0:
            plt.text(0.03, 0.9, labels[i], fontsize=fs, transform=ax.transAxes)
            plt.xscale('log')
            plt.xlim(1, 10)
            plt.ylim(yran1[i])
            ax.get_xaxis().set_major_formatter(ticker.ScalarFormatter())
            ax.set_xticks(xticks)
            plt.gca().xaxis.set_minor_formatter(ticker.NullFormatter())
        else:
            plt.xlim(xran[i])
            plt.ylim(yran2[i])
plt.savefig('../plots/pyratbay-exomol_opacity_comparison.pdf')
plt.savefig('../plots/pyratbay-exomol_opacity_comparison.png',
    dpi=300)

