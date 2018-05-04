#! /usr/bin/env python
import sys
import numpy as np
import scipy.optimize as so
import ConfigParser as configparser

sys.path.append("../pyratbay/modules/MCcubed")
import MCcubed.fit as fit


# Flat model:
def sfit(params, data=None, uncert=None, res=False):
  if res:
    return (data-params[0])/uncert
  return np.tile(params[0], len(data))

# Data files:
cfile = [
  "../run02_HATP11b/mcmc_hatp11b_wm-00000-c.cfg",
  "../run03_HATP32b/mcmc_hatp32b_wm-00000-c.cfg",
  "../run04_HATP38b/mcmc_hatp38b_wm-00000-c.cfg",
  "../run05_WASP043b/mcmc_wasp043b_wm-00000-c.cfg",
  "../run06_WASP063b/mcmc_wasp063b_wm-00000-c.cfg",
  "../run07_WASP067b/mcmc_wasp067b_wm-00000-c.cfg",
  "../run08_WASP101b/mcmc_wasp101b_wm-00000-c.cfg",
  "../run09_WASP107b/mcmc_wasp107b_wm-00000-c.cfg",
  ]
nplanets = len(cfile)


print("Flat-curve fit to transmission data:\n"
      "  Nfree = 1\n\n"
      "Planet      chi-square  red-chisq      BIC\n"
      "------------------------------------------")
for i in np.arange(nplanets):
  # Read data:
  config = configparser.SafeConfigParser()
  config.read(cfile[i])
  name = config.get('pyrat','logfile').split("_")[1]
  data   = np.array(config.get('pyrat','data').split(), np.double)
  uncert = np.array(config.get('pyrat','uncert').split(), np.double)

  indparams = [data]
  params    = np.array([np.mean(data)])
  stepsize  = np.array([1.0])
  chisq, bpars, bmodel, lsfit = fit.modelfit(params, sfit, data, uncert,
                                             indparams, stepsize, lm=True)
  dof = len(data) - len(params)
  print("{:10s}  {:10.3f}  {:9.3f}  {:7.3f}".
        format(name, chisq, chisq/dof, chisq+len(params)*np.log(len(data))))
