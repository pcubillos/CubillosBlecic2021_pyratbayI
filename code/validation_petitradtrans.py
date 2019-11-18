import numpy as np
from scipy.ndimage.filters import gaussian_filter1d as gaussf

import pyratbay as pb
import pyratbay.constants as pc

from petitRADTRANS import Radtrans
from petitRADTRANS import nat_cst as nc


def ticker(ax, xaxis=True, yaxis=True):
    if yaxis:
        ax2 = ax.twinx()
        ax2.tick_params(axis='y', direction='in', which='both')
        ax2.set_yscale(ax.get_yscale())
        ax2.set_yticklabels([])
        ax2.set_yticks(ax.get_yticks())
        ax2.set_ylim(ax.get_ylim())
    if xaxis:
        ax2 = ax.twiny()
        ax2.tick_params(axis='x', direction='in', which='both')
        ax2.set_xscale(ax.get_xscale())
        ax2.set_xticklabels([])
        ax2.set_xticks(ax.get_xticks())
        ax2.set_xlim(ax.get_xlim())

# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Preamble [post, since I need to run the code below ...]
iwn = atmosphere.cia_h2h2_lambda/pc.um <= 20.0
pb.io.write_cs('CIA_H2H2_molliere.dat',
    np.fliplr(atmosphere.cia_h2h2_alpha_grid[iwn].T),
    ['H2','H2'],
    atmosphere.cia_h2h2_temp,
    1/np.flipud(atmosphere.cia_h2h2_lambda[iwn]),
    header='# P Molliere data\n')

pb.io.write_atm('petit.atm', pyrat.atm.press, pyrat.atm.temp, pyrat.mol.name,
    pyrat.atm.q, 'bar', '# Petit profile\n',
    radius=atmosphere.radius, runits='km')

# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#atm_model = pb.run('petit_atmosphere.cfg')

pyrat = pb.run('petit_spectrum_no_lbl.cfg')
pyrat_wl = 1/(pyrat.spec.wn*pc.um)

# Transmission spectrum with no LBL opacities:
atmosphere = Radtrans(
    rayleigh_species = ['H2'],
    continuum_opacities = ['H2-H2'],
    wlen_bords_micron = [0.3, 10.0])
pressures = pyrat.atm.press/pc.bar
atmosphere.setup_opa_structure(pressures)
temperature = pyrat.atm.temp
abundances = {}
i = list(pyrat.mol.name).index('H2')
abundances['H2'] = pyrat.mol.mass[i]/pyrat.atm.mm * pyrat.atm.q[:,i]
i = list(pyrat.mol.name).index('He')
abundances['He'] = pyrat.mol.mass[i]/pyrat.atm.mm * pyrat.atm.q[:,i]
MMW  = pyrat.atm.mm
P0   = pyrat.atm.refpressure/pc.bar
R_pl = pyrat.phy.rplanet
gravity = pyrat.phy.gplanet
atmosphere.calc_transm(temperature, abundances, gravity, MMW,
    R_pl=R_pl, P0_bar=P0)
petit_wl = nc.c/atmosphere.freq/1e-4


pyrat_no_lbl = np.sqrt(pyrat.spec.spectrum) * pyrat.phy.rstar / pc.rjup
petit_no_lbl = atmosphere.transm_rad/pc.rjup


# Transmission spectrum:
pyrat = pb.run('petit_spectrum.cfg')

atmosphere = Radtrans(
    #line_species = ['H2O', 'CO_all_iso', 'Na', 'K'],
    line_species = ['CO_all_iso', 'Na', 'K'],
    rayleigh_species = ['H2'],
    continuum_opacities = ['H2-H2'],
    wlen_bords_micron = [0.3, 10.0])
pressures = pyrat.atm.press/pc.bar
atmosphere.setup_opa_structure(pressures)
temperature = pyrat.atm.temp
abundances = {}
i = list(pyrat.mol.name).index('H2')
abundances['H2'] = pyrat.mol.mass[i]/pyrat.atm.mm * pyrat.atm.q[:,i]
i = list(pyrat.mol.name).index('He')
abundances['He'] = pyrat.mol.mass[i]/pyrat.atm.mm * pyrat.atm.q[:,i]
i = list(pyrat.mol.name).index('H2O')
abundances['H2O'] = pyrat.mol.mass[i]/pyrat.atm.mm * pyrat.atm.q[:,i]
i = list(pyrat.mol.name).index('CO')
abundances['CO_all_iso'] = pyrat.mol.mass[i]/pyrat.atm.mm * pyrat.atm.q[:,i]
i = list(pyrat.mol.name).index('Na')
abundances['Na'] = pyrat.mol.mass[i]/pyrat.atm.mm * pyrat.atm.q[:,i]
i = list(pyrat.mol.name).index('K')
abundances['K'] = pyrat.mol.mass[i]/pyrat.atm.mm * pyrat.atm.q[:,i]
MMW  = pyrat.atm.mm
P0   = pyrat.atm.refpressure/pc.bar
R_pl = pyrat.phy.rplanet
gravity = pyrat.phy.gplanet
atmosphere.calc_transm(temperature, abundances, gravity, MMW,
    R_pl=R_pl, P0_bar=P0)

pyrat_trans = np.sqrt(pyrat.spec.spectrum) * pyrat.phy.rstar / pc.rjup
petit_trans = atmosphere.transm_rad/pc.rjup

fs = 12
plt.figure(11)
plt.clf()
ax = plt.subplot(111)
plt.plot(petit_wl, petit_no_lbl, c='navy')
plt.plot(pyrat_wl, pyrat_no_lbl, c='darkorange', alpha=0.8)
plt.plot(petit_wl, petit_trans, c='b', label='petitRADTRANS')
plt.plot(pyrat_wl, gaussf(pyrat_trans,7), c='orange', alpha=0.8,
    label='pyratbay')
plt.xscale('log')
plt.xlabel('Wavelength (microns)', fontsize=fs)
plt.ylabel(r'Transit radius ($\rm R_{Jup}$)', fontsize=fs)
plt.xlim(0.3, 10.0)
ax.tick_params(labelsize=fs-1)
ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
ax.set_xticks([0.3, 0.6, 1, 2, 3, 6, 10])
plt.gca().xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
plt.legend(loc='lower right', fontsize=fs)
plt.tight_layout()
ticker(ax)
plt.savefig('../plots/petit_pyratbay_comparison_transmission.pdf')


# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Emission:

kappa_IR = 0.01
gamma = 0.4
T_int = 200.
T_equ = 1500.
temp2 = nc.guillot_global(pressures, kappa_IR, gamma, gravity, T_int, T_equ)

atmosphere.calc_flux(temp2, abundances, gravity, MMW)
petit_emission = atmosphere.flux*pc.c

pyrat.od.path = 'eclipse'
pyrat.run(temp=temp2, radius=atmosphere.radius)
pyrat_emission = gaussf(pyrat.spec.spectrum, 5)

plt.figure(21)
plt.clf()
ax = plt.subplot(111)
plt.plot(petit_wl, petit_emission, c='b', label='petitRADTRANS')
plt.plot(pyrat_wl, pyrat_emission, c='orange', alpha=0.8, label='pyratbay')
plt.xscale('log')
plt.xlabel('Wavelength (microns)', fontsize=fs)
plt.ylabel(r'Planet flux $F_\nu$ (erg cm$^{-2}$ s$^{-1}$ cm)', fontsize=fs)
plt.xlim(0.3, 10.0)
ax.tick_params(labelsize=fs-1)
ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
ax.set_xticks([0.3, 0.6, 1, 2, 3, 6, 10])
plt.gca().xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
plt.legend(loc='lower right', fontsize=fs)
plt.tight_layout()
ticker(ax)
plt.savefig('../plots/petit_pyratbay_comparison_emission.pdf')


# Two-panel plot:
fs = 14
plt.figure(30, (7,9))
plt.clf()
ax = plt.subplot(211)
plt.plot(petit_wl, petit_no_lbl, c='navy')
plt.plot(pyrat_wl, pyrat_no_lbl, c='darkorange', alpha=0.8)
plt.plot(petit_wl, petit_trans, c='b', label='petitRADTRANS')
plt.plot(pyrat_wl, gaussf(pyrat_trans,7), c='orange', alpha=0.8,
    label='pyratbay')
plt.xscale('log')
plt.xlabel('Wavelength (microns)', fontsize=fs)
plt.ylabel(r'Transit radius ($\rm R_{Jup}$)', fontsize=fs)
plt.xlim(0.3, 10.0)
ax.tick_params(labelsize=fs-1)
ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
ax.set_xticks([0.3, 0.6, 1, 2, 3, 6, 10])
plt.gca().xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
plt.legend(loc='lower right', fontsize=fs)
plt.tight_layout()
ticker(ax)
ax = plt.subplot(212)
plt.plot(petit_wl, petit_emission, c='b')
plt.plot(pyrat_wl, pyrat_emission, c='orange', alpha=0.8)
plt.xscale('log')
plt.xlabel('Wavelength (microns)', fontsize=fs)
plt.ylabel(r'Planet flux $F_\nu$ (erg cm$^{-2}$ s$^{-1}$ cm)', fontsize=fs)
plt.xlim(0.3, 10.0)
ax.tick_params(labelsize=fs-1)
ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
ax.set_xticks([0.3, 0.6, 1, 2, 3, 6, 10])
plt.gca().xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
plt.tight_layout()
ticker(ax)
plt.savefig('../plots/petit_pyratbay_comparison.pdf')
plt.savefig('../plots/petit_pyratbay_comparison.ps')


