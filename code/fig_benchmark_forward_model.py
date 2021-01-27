import os
import sys

import numpy as np
from scipy.ndimage.filters import gaussian_filter1d as gaussf
import matplotlib
import matplotlib.pyplot as plt

import pyratbay as pb
import pyratbay.io as io
import pyratbay.atmosphere as pa
import pyratbay.constants as pc

if '__file__' not in locals():
    __file__ = os.getcwd() + '/'
sys.path.append(
    os.path.realpath(os.path.dirname(__file__) + '/../petitRADTRANS'))
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
# Transmission model with only H2 Rayleigh and CIA opacity:

nlayers = 41
T0 = 1300.0
pressure = pa.pressure('1e-8 bar', '100 bar', nlayers)
temperature = pa.tmodels.Isothermal(nlayers)(T0)

species = 'H2  He  Na  K  H2O  CH4  CO  CO2'.split()
q = [0.85, 0.1492, 3.0e-06, 1.0e-06, 4.0e-04, 1.0e-04, 3.0e-04, 1.0e-07]
Q = pa.uniform(pressure, temperature, species, q)
mu = pa.mean_weight(Q, species)

iH2  = species.index('H2')
iHe  = species.index('He')
iH2O = species.index('H2O')
iCO  = species.index('CO')
iCH4 = species.index('CH4')
iNa  = species.index('Na')
iK   = species.index('K')

molecs, mass, diameter = io.read_molecs(f'{pc.ROOT}/inputs/molecules.dat')
H2_mass  = mass[molecs=='H2'][0]
He_mass  = mass[molecs=='He'][0]
H2O_mass = mass[molecs=='H2O'][0]
CO_mass  = mass[molecs=='CO'][0]
CH4_mass = mass[molecs=='CH4'][0]
Na_mass  = mass[molecs=='Na'][0]
K_mass   = mass[molecs=='K'][0]

# Reference pressure-radius point (bar and CGS units):
P0_bar = 0.1
R_pl = 1.0 * pc.rjup
gravity = 1487.2

# Transmission spectrum with no LBL opacities:
atmosphere = Radtrans(
    rayleigh_species = ['H2'],
    continuum_opacities = ['H2-H2'],
    wlen_bords_micron = [0.3, 10.0])
atmosphere.setup_opa_structure(pressure/pc.bar)

# petit uses abundances in mass mixing fraction units:
abundances = {}
abundances['H2'] = H2_mass/mu * Q[:,iH2]
abundances['He'] = He_mass/mu * Q[:,iHe]
atmosphere.calc_transm(temperature, abundances, gravity, mu,
    R_pl=R_pl, P0_bar=P0_bar)
petit_wl = nc.c/atmosphere.freq/1e-4
petit_no_lbl = atmosphere.transm_rad/pc.rjup

# Create atmospheric file with same values as petitRADTRANS atmosphere:
io.write_atm('petit.atm', pressure, temperature, species, Q,
    header='# Petit profile\n', radius=atmosphere.radius,
    punits='bar', runits='km')

# Run pyratbay model (using the atmospheric file above as input):
pyrat = pb.run('petit_spectrum_CIA-Rayleigh.cfg')
pyrat_wl_no_lbl = 1/(pyrat.spec.wn*pc.um)
# Transmission radius in Rjup units:
pyrat_no_lbl = np.sqrt(pyrat.spec.spectrum) * pyrat.phy.rstar / pc.rjup


# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Now, transmission model with line transitions:

atmosphere = Radtrans(
    line_species = ['H2O', 'CO_all_iso'],
    rayleigh_species = ['H2'],
    continuum_opacities = ['H2-H2'],
    wlen_bords_micron = [0.3, 10.0])
atmosphere.setup_opa_structure(pressure/pc.bar)
abundances = {}
abundances['H2']         = H2_mass /mu * Q[:,iH2]
abundances['He']         = He_mass /mu * Q[:,iHe]
abundances['H2O']        = H2O_mass/mu * Q[:,iH2O]
abundances['CO_all_iso'] = CO_mass /mu * Q[:,iCO]
atmosphere.calc_transm(temperature, abundances, gravity, mu,
    R_pl=R_pl, P0_bar=P0_bar)
petit_trans = atmosphere.transm_rad/pc.rjup

pyrat = pb.run('petit_spectrum_H2O-CO.cfg')
pyrat_wl = 1/(pyrat.spec.wn*pc.um)
pyrat_trans = np.sqrt(pyrat.spec.spectrum) * pyrat.phy.rstar / pc.rjup


# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Emission models:

# This one requires a non-isothermal profile:
kappa_IR = 0.01
gamma = 0.4
T_int = 200.0
T_equ = 1500.0
temp2 = nc.guillot_global(
    pressure/pc.bar, kappa_IR, gamma, gravity, T_int, T_equ)

atmosphere.calc_flux(temp2, abundances, gravity, mu)
petit_emission = atmosphere.flux*pc.c

pyrat.od.path = 'eclipse'
pyrat.run(temp=temp2, radius=atmosphere.radius)
pyrat_emission = np.copy(pyrat.spec.spectrum)


# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Additional calculation with methane opacity:
atmosphere = Radtrans(
    line_species = ['CH4', 'Na_lor_cut', 'K_lor_cut'],
    rayleigh_species = ['H2'],
    continuum_opacities = ['H2-H2'],
    wlen_bords_micron = [0.3, 10.0])
atmosphere.setup_opa_structure(pressure/pc.bar)
abundances = {}
abundances['H2']  = H2_mass /mu * Q[:,iH2]
abundances['He']  = He_mass /mu * Q[:,iHe]
abundances['CH4'] = CH4_mass/mu * Q[:,iCH4]
abundances['Na_lor_cut'] = Na_mass /mu * Q[:,iNa]
abundances['K_lor_cut']  = K_mass  /mu * Q[:,iK]
# Transmission model:
atmosphere.calc_transm(temperature, abundances, gravity, mu,
    R_pl=R_pl, P0_bar=P0_bar)
petit_wl = nc.c/atmosphere.freq/1e-4
petit_methane = atmosphere.transm_rad/pc.rjup
# And emission model:
atmosphere.calc_flux(temp2, abundances, gravity, mu)
petit_methane_emission = atmosphere.flux*pc.c

# Run pyratbay model (emission and transmission):
pyrat = pb.run('petit_spectrum_CH4-Na-K.cfg')
pyrat_wl = 1/(pyrat.spec.wn*pc.um)
pyrat_methane = np.sqrt(pyrat.spec.spectrum) * pyrat.phy.rstar / pc.rjup
pyrat.od.path = 'eclipse'
pyrat.run(temp=temp2, radius=atmosphere.radius)
pyrat_methane_emission = np.copy(pyrat.spec.spectrum)

# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Now, plot everything:

# T = 1359, 1723
planck_hot  = pb.spectrum.bbflux(pyrat.spec.wn, np.amax(temp2))
planck_cold = pb.spectrum.bbflux(pyrat.spec.wn, np.amin(temp2))

fs = 14
xticks = [0.3, 0.6, 1, 2, 3, 6, 10]
xran = 0.5, 10.0
yrant = 0.983, 1.046
yrane = -2000, 96000

def boliler_plate(plt, ax):
    plt.xscale('log')
    plt.xlabel('Wavelength (microns)', fontsize=fs)
    ax.tick_params(labelsize=fs-1)
    ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax.set_xticks(xticks)
    plt.gca().xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    plt.xlim(0.5, 10)
    ticker(ax)


plt.figure(325, (12.5,8))
plt.clf()
plt.subplots_adjust(0.085, 0.07, 0.98, 0.99, hspace=0.18, wspace=0.22)
ax = plt.subplot(221)
plt.text(0.02, 0.92, 'H2O, CO, H2 Ray, H2-H2 CIA', fontsize=fs,
    transform=ax.transAxes)
plt.plot(pyrat_wl_no_lbl, pyrat_no_lbl, c='navy')
plt.plot(petit_wl, petit_no_lbl, c='darkorange', alpha=0.75)
plt.plot(1e4/pyrat.spec.wn, gaussf(pyrat_trans,12), c='blue', label='pyratbay')
plt.plot(petit_wl, petit_trans, c='orange', label='petitRADTRANS', alpha=0.8)
plt.ylabel(r'Transit radius ($\rm R_{Jup}$)', fontsize=fs)
plt.ylim(yrant)
plt.legend(loc='lower right', fontsize=fs)
boliler_plate(plt, ax)
ax = plt.subplot(223)
plt.plot(pyrat_wl, planck_hot, c='0.5')
plt.plot(pyrat_wl, planck_cold, c='0.5')
plt.plot(pyrat_wl, gaussf(pyrat_emission, 10.0), c='blue')
plt.plot(petit_wl, petit_emission, c='orange', alpha=0.75)
plt.ylabel(r'Planet flux (erg s$^{-1}$ cm$^{-2}$ cm)', fontsize=fs)
plt.ylim(yrane)
boliler_plate(plt, ax)
ax = plt.axes([0.139, 0.295, 0.08, 0.18])
plt.plot(temp2, pressure/pc.bar, c='red', lw=2)
ax.set_yscale('log')
ax.set_ylim(np.amax(pressure/pc.bar), np.amin(pressure/pc.bar))
ax.set_xlabel('Temperature (K)', fontsize=fs-3)
ax.set_ylabel('Pressure (bar)', fontsize=fs-3)
ax = plt.subplot(222)
plt.text(0.02, 0.92, 'CH4, Na, K, H2 Ray, H2-H2 CIA', fontsize=fs,
    transform=ax.transAxes)
ax.plot(pyrat_wl, gaussf(pyrat_methane, 8.0), c='blue')
ax.plot(petit_wl, petit_methane, c='orange', alpha=0.75)
ax.set_ylabel(r'Transit radius ($\rm R_{Jup}$)', fontsize=fs)
plt.ylim(yrant)
boliler_plate(plt, ax)
ax = plt.subplot(224)
plt.plot(pyrat_wl, planck_hot, c='0.5')
plt.plot(pyrat_wl, planck_cold, c='0.5')
plt.plot(pyrat_wl, gaussf(pyrat_methane_emission, 8.0), c='blue')
plt.plot(petit_wl, petit_methane_emission, c='orange', alpha=0.75)
plt.ylabel(r'Planet flux (erg s$^{-1}$ cm$^{-2}$ cm)', fontsize=fs)
plt.ylim(yrane)
boliler_plate(plt, ax)
plt.savefig('../plots/pyratbay-petit_spectrum_comparison.pdf')
plt.savefig('../plots/pyratbay-petit_spectrum_comparison.png', dpi=300)
