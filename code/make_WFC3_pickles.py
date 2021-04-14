#! /usr/bin/env python

import os
import numpy as np

import pyratbay as pb
import pyratbay.io as io
import pyratbay.tools as pt

import mc3


cfiles = [
    'run_GJ-0436b_tsiaras/mcmc_GJ-0436b_tsiaras_w00000-mc.cfg',
    'run_HATP-03b_tsiaras/mcmc_HATP-03b_tsiaras_w00000-mc.cfg',
    'run_HATP-17b_tsiaras/mcmc_HATP-17b_tsiaras_w00000-mc.cfg',
    'run_HATP-18b_tsiaras/mcmc_HATP-18b_tsiaras_w00000-mc.cfg',
    'run_HATP-32b_damiano/mcmc_HATP-32b_damiano_w00000-mc.cfg',
    'run_HATP-32b_tsiaras/mcmc_HATP-32b_tsiaras_w00000-mc.cfg',
    'run_HATP-38b_bruno/mcmc_HATP-38b_bruno_w00000-mc.cfg',

    'run_HATP-38b_tsiaras/mcmc_HATP-38b_tsiaras_w00000-mc.cfg',
    'run_HATP-41b_tsiaras/mcmc_HATP-41b_tsiaras_w00000-mc.cfg',
    'run_HD-149026b_tsiaras/mcmc_HD-149026b_tsiaras_w00000-mc.cfg',
    'run_K2-018b_benneke/mcmc_K2-018b_benneke_wmdm00-mc.cfg',
    'run_K2-018b_tsiaras/mcmc_K2-018b_tsiaras_w00000-mc.cfg',
    'run_WASP-029b_tsiaras/mcmc_WASP-029b_tsiaras_w00000-mc.cfg',
    'run_WASP-043b_stevenson/mcmc_WASP-043b_stevenson_wmdm00-mc.cfg',

    'run_WASP-043b_tsiaras/mcmc_WASP-043b_tsiaras_w00000-mc.cfg',
    'run_WASP-063b_kilpatrick/mcmc_WASP-063b_kilpatrick_w000h0-mc.cfg',
    'run_WASP-063b_tsiaras/mcmc_WASP-063b_tsiaras_w000h0-mc.cfg',
    'run_WASP-067b_bruno/mcmc_WASP-067b_bruno_w00000-mc.cfg',
    'run_WASP-067b_tsiaras/mcmc_WASP-067b_tsiaras_w00000-mc.cfg',
    'run_WASP-069b_tsiaras/mcmc_WASP-069b_tsiaras_w00000-mc.cfg',
    'run_WASP-074b_tsiaras/mcmc_WASP-074b_tsiaras_w00000-mc.cfg',

    'run_WASP-080b_tsiaras/mcmc_WASP-080b_tsiaras_w00000-mc.cfg',
    'run_WASP-101b_tsiaras/mcmc_WASP-101b_tsiaras_w00000-mc.cfg',
    'run_WASP-101b_wakeford/mcmc_WASP-101b_wakeford_w00000-mc.cfg',
    'run_WASP-103b_kreidberg/mcmc_WASP-103b_kreidberg_wmdm00-mc.cfg',
    'run_XO-1b_tsiaras/mcmc_XO-1b_tsiaras_w00000-mc.cfg',
    ]
ntargets = len(cfiles)

for cfile in cfiles:
    folder, cfg = os.path.split(cfile)
    _, planet, dset = folder.split('_')
    with pt.cd(folder):
        pyrat = pb.run(cfg, run_step='init', no_logfile=True)
        with np.load(pyrat.ret.mcmcfile) as mcmc:
            post, zchain, zmask = mc3.utils.burn(mcmc)
            accept_rate = mcmc['acceptance_rate']
        pyrat.ret.posterior = post
        # I want to evaluate less than ~1e5 unique samples:
        nmax = int(1e5 / (6 * 0.01 * accept_rate))
        pyrat.percentile_spectrum(nmax)
        pfile = pyrat.ret.mcmcfile.replace('.npz', '.pickle')
        io.save_pyrat(pyrat, pfile)

