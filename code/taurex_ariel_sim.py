import os
from itertools import islice

import numpy as np

import taurex.log
taurex.log.disableLogging()
from taurex.cache import OpacityCache, CIACache
from taurex.planet import Planet
from taurex.stellar import BlackbodyStar
from taurex.chemistry import TaurexChemistry, ConstantGas
from taurex.model import TransmissionModel, EmissionModel
from taurex.temperature import Guillot2010, Isothermal
from taurex.contributions import CIAContribution, AbsorptionContribution, \
    RayleighContribution, SimpleCloudsContribution

import pyratbay.io as io
import pyratbay.tools as pt
import pyratbay.constants as pc


def parse_list(array, n, fmt='.8f'):
    """Pretty parsing list of numbers into text"""
    iter_values = iter(array)
    s = []
    while True:
        values = list(islice(iter_values, n))
        if len(values) == 0:
            break
        s += ['  '.join(f'{val:{fmt}}' for val in values)]
    return s


def simulate_ariel(rplanet, mplanet, tplanet, qplanet, pcloud,
        rstar, tstar, roffset=0.0, mode='transit', tpars=None,
        kmag=None, noise_floor=0.00002, band_trans=None, band_wn=None):
    """
    Generate a TauREx synthetic model
    If kmag is not None, also generate synthetic ARIEL observations.
    There is plenty of hardcoded stuff down here.
    """
    # Prep up TauREx:
    inputs_dir = '../inputs/taurex/'
    OpacityCache().clear_cache()
    OpacityCache().set_opacity_path(inputs_dir)
    CIACache().set_cia_path(inputs_dir)

    kflux_HD209 = 1780045.0
    kmag_HD209 = 6.31
    # Pressure boundaries (pascal):
    pmin = 1.0e-5 * pc.bar/pc.pascal
    pmax = 1.0e+2 * pc.bar/pc.pascal
    nlayers = 22

    if tpars is None:
        kappa_irr = 0.01
        kappa_v1 = 0.005
        kappa_v2 = 0.1
    else:
        kappa_irr = 10**tpars[0]
        kappa_v1 = 10**(tpars[1]+tpars[0])
        kappa_v2 = 0.1

    if mode == 'transit':
        taurex_model = TransmissionModel
        tmodel = Isothermal(T=tplanet)
    elif mode == 'eclipse':
        taurex_model = EmissionModel
        tmodel = Guillot2010(
            T_irr=tplanet, kappa_irr=kappa_irr,
            kappa_v1=kappa_v1, kappa_v2=kappa_v2, alpha=0.0)

    tr_rplanet = rplanet/pc.rjup + roffset
    planet = Planet(planet_radius=tr_rplanet, planet_mass=mplanet/pc.mjup)
    star = BlackbodyStar(temperature=tstar, radius=rstar)

    chemistry = TaurexChemistry(fill_gases=['H2', 'He'], ratio=0.172)
    chemistry.addGas(ConstantGas('H2O', mix_ratio=10**qplanet[0]))
    chemistry.addGas(ConstantGas('CH4', mix_ratio=10**qplanet[1]))
    chemistry.addGas(ConstantGas('CO',  mix_ratio=10**qplanet[2]))
    chemistry.addGas(ConstantGas('CO2', mix_ratio=10**qplanet[3]))

    tm = taurex_model(
        planet=planet, temperature_profile=tmodel, chemistry=chemistry,
        star=star,
        atm_min_pressure=pmin, atm_max_pressure=pmax, nlayers=nlayers)

    tm.add_contribution(AbsorptionContribution())
    tm.add_contribution(CIAContribution(cia_pairs=['H2-H2']))
    tm.add_contribution(RayleighContribution())
    tm.add_contribution(
        SimpleCloudsContribution(clouds_pressure=pcloud*pc.bar/pc.pascal))
    tm.build()
    wn_tau, depth, tau, _ = tm.model()
    wl_tau = 1/(wn_tau*pc.um)

    if kmag is None:
        return tm, wn_tau, depth

    # Stellar flux in erg s-1 cm-2 cm:
    starflux = 1e7 * tm.star.sed / wn_tau**2
    # Scale flux according to Kmag:
    kband_idx = np.where(wl_tau<2.2)[0][0]
    kband_flux = starflux[kband_idx]
    starflux = starflux/kband_flux * kflux_HD209 * 10**(0.4*(kmag_HD209-kmag))

    # Sum flux in bands, rather than average:
    sflux = []
    for btrans, bwn in zip(band_trans, band_wn):
        resampled, wnidx = pt.resample(btrans, bwn, wn_tau)
        resampled /= np.amax(resampled)
        sflux.append(np.trapz(starflux[wnidx]*resampled, wn_tau[wnidx]))
    sflux = np.array(sflux)
    bin_depth = np.array(pt.band_integrate(depth, wn_tau, band_trans, band_wn))

    # Poisson noise for transmission spectroscopy:
    if mode == 'transit':
        Fout = sflux
        Fin  = sflux * (1-bin_depth)
    elif mode == 'eclipse':
        Fout = sflux * (1 + bin_depth)
        Fin  = sflux
    # This ignores integration times:
    sigma_out = np.sqrt(Fout)
    sigma_in  = np.sqrt(Fin)
    poisson_noise = np.sqrt((sigma_in/Fout)**2 + (Fin*sigma_out/Fout**2)**2)
    # Flatten it a bit:
    poisson_noise = 0.2*poisson_noise + 0.8*np.amin(poisson_noise)
    # Add noise floor in quadrature:
    noise = np.sqrt(poisson_noise**2+noise_floor**2)

    return tm, wn_tau, depth, bin_depth, bin_wl, noise



if __name__ == "__main__":
    # Prep up ARIEL filters:
    filter_dir = '../inputs/filters/'
    filters = [f'{filter_dir}{ffile}' for ffile in os.listdir(filter_dir)
               if ffile.startswith('ARIEL')]
    bin_wl = [ffile.split('_')[-1].split('um')[0] for ffile in filters]
    wn_sort = np.argsort(np.array(bin_wl, float))
    filters = [filters[idx] for idx in wn_sort]
    band_trans, band_wn = [], []
    for filter_file in filters:
        bwn, btrans = io.read_spectrum(filter_file)
        band_trans.append(btrans)
        band_wn.append(bwn)


    # Generate synthetic ARIEL transmission observations with TauREx:
    planet_data = np.loadtxt('retrieval_validation_transmission.txt', dtype=str)
    planets = planet_data[:,0]
    sys_params = np.array(planet_data[:,1:7], np.double)
    qplanets   = np.array(planet_data[:,7:11], np.double)
    pclouds    = np.array(planet_data[:,11], np.double)
    noise_floor = 10**np.array(planet_data[:,12], np.double)
    nplanets = len(planets)

    sim_file = open('ariel_taurex_synthetic_transmission.dat', 'w')
    for i in range(nplanets):
        rstar, tstar, kmag, rplanet, mplanet, tplanet = sys_params[i]
        rplanet *= pc.rearth
        mplanet *= pc.mearth
        qplanet = qplanets[i]
        pcloud = 10**pclouds[i]
        tm, wn_tau, depth, bin_depth, bin_wl, uncert = simulate_ariel(
            rplanet, mplanet, tplanet, qplanet, pcloud, rstar, tstar,
            mode='transit', noise_floor=noise_floor[i], kmag=kmag,
            band_trans=band_trans, band_wn=band_wn)

        data   = '    ' + '\n    '.join(parse_list(bin_depth*100, 5))
        uncert = '    ' + '\n    '.join(parse_list(uncert*100, 5))

        sim_file.write(f'\n{planets[i]}:\n')
        sim_file.write('dunits = percent\ndata =\n')
        sim_file.write(data)
        sim_file.write('\nuncert =\n')
        sim_file.write(uncert)
        sim_file.write('\n'
            f'rstar = {rstar:.3f} rsun\n'
            f'tstar = {tstar:.1f}\n'
            f'rplanet = {rplanet/pc.rearth:.2f} rearth\n'
            f'mplanet = {mplanet/pc.mearth:.2f} mearth\n\n')
    sim_file.close()


    # Generate synthetic ARIEL emission obseravtions with TauREx:
    planet_data = np.loadtxt('retrieval_validation_emission.txt', dtype=str)
    planets = planet_data[:,0]
    sys_params = np.array(planet_data[:,1:7], np.double)
    tpars      = np.array(planet_data[:,7:9], np.double)
    qplanets   = np.array(planet_data[:,9:13], np.double)
    nplanets = len(planets)
    pcloud = 100.0
    noise_floor = 2e-5

    sim_file = open('ariel_taurex_synthetic_emission.dat', 'w')
    for i in range(nplanets):
        rstar, tstar, kmag, rplanet, mplanet, tplanet = sys_params[i]
        rplanet *= pc.rearth
        mplanet *= pc.mearth
        qplanet = qplanets[i]
        tm, wn_tau, depth, bin_depth, bin_wl, uncert = simulate_ariel(
            rplanet, mplanet, tplanet, qplanet, pcloud, rstar, tstar,
            tpars=tpars[i], mode='eclipse', noise_floor=noise_floor,
            kmag=kmag, band_trans=band_trans, band_wn=band_wn)

        data = '    ' + '\n    '.join(parse_list(bin_depth/pc.ppm, 8, '7.2f'))
        uncert = '    ' + '\n    '.join(parse_list(uncert/pc.ppm, 8, '7.2f'))

        sim_file.write(f'\n{planets[i]}:\n')
        sim_file.write('dunits = ppm\ndata =\n')
        sim_file.write(data)
        sim_file.write('\nuncert =\n')
        sim_file.write(uncert)
        sim_file.write('\n'
            f'rstar = {rstar:.3f} rsun\n'
            f'tstar = {tstar:.1f}\n'
            f'rplanet = {rplanet/pc.rearth:.2f} rearth\n'
            f'mplanet = {mplanet/pc.mearth:.2f} mearth\n\n')
    sim_file.close()
