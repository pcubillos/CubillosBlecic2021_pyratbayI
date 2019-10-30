import numpy as np
from scipy.ndimage.filters import gaussian_filter1d as gaussf

import pyratbay as pb
import pyratbay.constants as pc

from petitRADTRANS import Radtrans
from petitRADTRANS import nat_cst as nc


pyrat = pb.run('petit_spectrum.cfg')

atmosphere = Radtrans(
    #line_species = ['H2O', 'CO_all_iso', 'CO2', 'Na', 'K'],
    line_species = ['H2O', 'CO_all_iso', 'Na', 'K'],
    rayleigh_species = ['H2'],
    continuum_opacities = ['H2-H2', 'H2-He'],
    wlen_bords_micron = [0.3, 15.0])

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
#i = list(pyrat.mol.name).index('CO2')
#abundances['CO2'] = 0.00001 * np.ones_like(temperature)
#abundances['CH4'] = 0.000001 * np.ones_like(temperature)

MMW = pyrat.atm.mm

P0   = pyrat.atm.refpressure/pc.bar
R_pl = pyrat.phy.rplanet
gravity = pyrat.phy.gplanet

atmosphere.calc_transm(temperature, abundances, gravity, MMW,
    R_pl=R_pl, P0_bar=P0)

# Plot:
spec = np.sqrt(pyrat.spec.spectrum) * pyrat.phy.rstar / pc.rjup

plt.figure(10)
plt.clf()
plt.plot(nc.c/atmosphere.freq/1e-4, atmosphere.transm_rad/pc.rjup,
    c='b', label='petitRADTRANS')
plt.plot(1e4/pyrat.spec.wn, 1.004*gaussf(spec,8), c='orange', alpha=0.8,
    label='pyratbay')
plt.xscale('log')
plt.xlabel('Wavelength (microns)')
plt.ylabel(r'Transit radius ($\rm R_{Jup}$)')
plt.xlim(0.4, 10.0)
plt.legend(loc='lower right')
plt.savefig('../plots/petit_pyratbay_comparison.pdf')


plt.figure(11)
plt.clf()
plt.plot(nc.c/atmosphere.freq/1e-4, (atmosphere.transm_rad/pyrat.phy.rstar)**2,
    c='b', label='petitRADTRANS')
plt.plot(1e4/pyrat.spec.wn, 1.008*gaussf(pyrat.spec.spectrum,8),
    c='orange', alpha=0.8, label='pyratbay')
plt.xscale('log')
plt.xlabel('Wavelength (microns)')
plt.ylabel(r'Transit radius ($\rm R_{Jup}$)')
plt.xlim(0.4, 10.0)
plt.legend(loc='lower right')


