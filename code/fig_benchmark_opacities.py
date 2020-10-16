import os

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker

import pyratbay as pb
import pyratbay.io as io
import pyratbay.constants as pc
import pyratbay.atmosphere as pa


HOME = os.path.abspath('..')

# Read ExoMol cross-section opacity data:
xs_H2O, press_exomol, temp_exomol, wn_exomol, _ = io.import_xs(
    f'{HOME}/inputs/taurex/1H2-16O__POKAZATEL__R15000_0.3-50mu.xsec.TauREx.h5',
    source='exomol')
xs_CO2 = io.import_xs(
    f'{HOME}/inputs/taurex/12C-16O2__UCL-4000.R15000_0.3-50mu.xsec.TauREx.h5',
    source='exomol', read_all=False)
xs_CO = io.import_xs(
    f'{HOME}/inputs/taurex/C-O-NatAbund__Li2015.R15000_0.3-50mu.xsec.TauREx.h5',
    source='exomol', read_all=False)

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
ipress = 3
itemp = 9

pyrat = pb.run('H2O_spectrum.cfg', init=True)
wl_pyrat = 1.0/(pyrat.spec.wn*pc.um)
ec_H2O = pyrat.get_ec(ipress)[0][0]

pyrat = pb.run('CO_spectrum.cfg', init=True)
ec_CO = pyrat.get_ec(ipress)[0][0]

pyrat = pb.run('CO2_spectrum.cfg', init=True)
pyrat.iso.ratio[1:] = 1e-30  # Consider only the first isotope
ec_CO2 = pyrat.get_ec(ipress)[0][0]

# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

iH2O = species.index('H2O')
iCO  = species.index('CO')
iCO2 = species.index('CO2')
ec = [ec_H2O/pyrat.atm.d[ipress,iH2O],
      ec_CO /pyrat.atm.d[ipress,iCO],
      ec_CO2/pyrat.atm.d[ipress,iCO2]]
xs = xs_H2O[ipress,itemp], xs_CO[ipress,itemp], xs_CO2[ipress,itemp]

fs = 10
lw = 1.25
xlim = 1.0, 10.0
xticks = [1, 2, 3, 4, 6, 10]
xran = [(1.4, 1.432), (1.177, 1.232), (2.6, 2.73)]
xran = [(1.409, 1.432), (1.177, 1.232), (2.62, 2.685)]
yran1 = [(3e-27, 5e-14), (1e-48, 1e-15), (1e-33, 1e-15)]
yran2 = [(1e-25, 1e-18), (8e-34, 3e-24), (1e-25, 1e-18)]
labels = ['H2O/ExoMol-POKAZATEL', 'CO/HITEMP-Li', 'CO2/Exomol-UCL4000']

plt.figure(150, (8, 7))
plt.clf()
ax = plt.subplot(111)
plt.subplots_adjust(0.1, 0.08, 0.99, 0.99, wspace=0.2, hspace=0.15)
for i in range(len(xs)):
    for j in [0,1]:
        ax = plt.subplot(3, 2, 2*i+1+j)
        plt.plot(wl_pyrat, ec[i], c='b', lw=lw, label='pyratbay')
        plt.plot(wl_exomol, xs[i], c='orange', lw=lw, label='exomol', alpha=0.7)
        plt.yscale('log')
        if i == 0 and j == 0:
            plt.legend(loc='upper right', fontsize=fs)
        if i == 1 and j == 0:
            plt.ylabel(
                r'Opacity cross section (cm$^{2}$ molec$^{-1}$)', fontsize=fs+1)
        if i == 0 and  j == 0:
            plt.text(0.03, 0.8, f'$T$ = {temp_exomol[itemp]:.0f} K',
                     fontsize=fs, transform=ax.transAxes)
            plt.text(0.03, 0.7, f'$p$ = {press_exomol[ipress]/pc.bar:.1f} bar',
                     fontsize=fs, transform=ax.transAxes)
        if i == 0 and  j == 1:
            plt.xticks(np.arange(1.4, 1.44, 0.01))
        if i == 2:
            plt.xlabel('Wavelength (um)', fontsize=fs+1)
        if j == 0:
            plt.text(0.03, 0.9, labels[i], fontsize=fs, transform=ax.transAxes)
            plt.xscale('log')
            plt.ylim(yran1[i])
            ax.get_xaxis().set_major_formatter(ticker.ScalarFormatter())
            ax.set_xticks(xticks)
            plt.gca().xaxis.set_minor_formatter(ticker.NullFormatter())
            plt.xlim(1.0, 10.0)
        else:
            plt.xlim(xran[i])
            plt.ylim(yran2[i])
plt.savefig('../plots/pyratbay-exomol_opacity_comparison.pdf')
plt.savefig('../plots/pyratbay-exomol_opacity_comparison.png',
    dpi=300)

