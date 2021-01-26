import os

import pyratbay.io as io
import pyratbay.atmosphere as pa


# Reformat TauREx cross section into pyratbay format:
taurex_xs_files = [
    'inputs/taurex/H2O.R15000.TauREx.pickle',
    'inputs/taurex/CH4.R15000.TauREx.pickle',
    'inputs/taurex/CO.R15000.TauREx.pickle',
    'inputs/taurex/CO2.R15000.TauREx.pickle',
   ]
for tr_file in taurex_xs_files:
    mol = os.path.basename(tr_file).split('.')[0]
    pb_file = f'run_setup/taurex_{mol}_opacity_R15000.npz'
    xs, press, temp, wn, specs = io.import_xs(
        tr_file, source='taurex', ofile=pb_file)

# Generate an atmospheric profile at TauREx pressure levels:
T0 = 1480.0
nlayers = len(press)
temperature = pa.tmodels.Isothermal(nlayers)(T0)
species = 'H2     He      H2O      CH4      CO       CO2'.split()
abund   = [0.853, 0.1467, 1.0e-04, 1.0e-04, 1.0e-04, 1.0e-07]
q = pa.uniform(press, temperature, species, abund)
io.write_atm(
    'run_benchmark_retrieval/taurex.atm',
    press, temperature, species, q,
    header='# TauREx profile\n', punits='bar')

