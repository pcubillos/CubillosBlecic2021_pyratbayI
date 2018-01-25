#! /usr/bin/env python
import sys
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import scipy.constants   as sc
from scipy.ndimage.filters import gaussian_filter1d as gaussf
#plt.ioff()
plt.ion()

# Edit the path to the Pyrat-Bay package if necessary:
sys.path.append("pyratbay")
import pyratbay as pb
import pyratbay.constants  as pc
#import pyratbay.tools      as pt
import pyratbay.wine       as pw
#import pyratbay.atmosphere as pa
import pyratbay.starspec as ps


# Load best-fitting models:
paths = ["run02_HATP11b/",  "run03_HATP32b/",  "run04_HATP38b/",
         "run05_WASP043b/", "run06_WASP063b/", "run07_WASP067b/",
         "run08_WASP101b/", "run09_WASP107b/"]
bestm =["MCMC_HAT-P-11b_1.0-5.5um_w0-00000-c_bestfit_spectrum.dat",
        "MCMC_HAT-P-32b_1.0-5.5um_w0-00000-c_bestfit_spectrum.dat",
        "MCMC_HAT-P-38b_1.0-5.5um_w0-00000-c_bestfit_spectrum.dat",
        "MCMC_WASP-043b_1.0-5.5um_w0-00000-c_bestfit_spectrum.dat",
        "MCMC_WASP-063b_1.0-5.5um_w0-000h0-c_bestfit_spectrum.dat",
        "MCMC_WASP-067b_1.0-5.5um_w0-00000-c_bestfit_spectrum.dat",
        "MCMC_WASP-101b_1.0-5.5um_wm-00000-0_bestfit_spectrum.dat",
        "MCMC_WASP-107b_1.0-5.5um_wm-00000-c_bestfit_spectrum.dat",
]

nplanets = len(paths)
spec = []
for j in np.arange(nplanets):
  try:
    wn, fl = ps.readpyrat(paths[j] + bestm[j])
    spec.append(fl)
  except:
    spec.append(None)
wl = 1.0/(wn*pc.um)


# Load data:
with open("code/filter_info.txt", "r") as f:
  lines = f.readlines()

data    = []
uncert  = []
filters = []
name    = []
i = 0
for k in np.arange(nplanets):
  name.append(lines[i].strip())
  while not lines[i].startswith("data"):
    i += 1
    continue
  info = lines[i]
  i += 1
  while not lines[i].startswith("uncert"):
    info += lines[i]
    i += 1
  data.append(np.asarray(info.split()[2:], np.double))
  info = lines[i]
  i += 1
  while not lines[i].startswith("filter"):
    info += lines[i]
    i += 1
  uncert.append(np.asarray(info.split()[2:], np.double))
  info = lines[i]
  i += 1
  while not lines[i].strip() == "":
    info += lines[i]
    i += 1
  filters.append(info.split()[2:])
  i += 1

bandwl = []
bandfl = []

for j in np.arange(nplanets):
  bwl = np.zeros(len(data[j]))
  bfl = np.zeros(len(data[j]))
  for i in np.arange(len(data[j])):
    fwn, ftr = pw.readfilter("code/"+filters[j][i])
    nifilter, wnidx, dummy = pw.resample(wn, fwn, ftr)
    bwl[i] = 1.0/(np.sum(fwn*ftr)/np.sum(ftr)*pc.um)
    if spec[j] is not None:
      bfl[i] = pw.bandintegrate(spec[j], wn, wnidx, nifilter)
  bandwl.append(bwl)
  bandfl.append(bfl)

# Plot it:
lw = 1.25
ms = 5
fs = 14
sigma = 2.0
xticks = [1.0, 1.3, 1.6, 2.0, 3.0, 4.0, 5.0]
fscale = 100

plt.figure(0, (7, 9))
plt.clf()
plt.subplots_adjust(0.12, 0.06, 0.97, 0.97, hspace=0.04, wspace=0.04)
for j in np.arange(nplanets):
  ax = plt.subplot(8,1,1+j)
  if spec[j] is not None:
    plt.plot(wl, fscale*gaussf(spec[j], sigma), color='orange', lw=lw)
    plt.plot(bandwl[j], fscale*bandfl[j], "o", color='orange', mew=0.15,
             mec="0.25", ms=ms)
  plt.errorbar(bandwl[j], fscale*data[j], uncert[j]*fscale, fmt="ob",
               label="Data", ms=ms, mew=0, elinewidth=lw, capthick=lw,
               capsize=0, zorder=3)
  plt.xlim(1.0, 5.5)
  plt.xscale('log')
  ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
  ax.set_xticks(xticks)
  if j < nplanets-1:
    ax.set_xticklabels([])
  ax.locator_params(axis='y', nbins=7)
  plt.text(0.9, 0.8, name[j], horizontalalignment='center',
           transform=ax.transAxes, fontsize=fs-2)
  # 100 ppm:
  bot, top = ax.get_ylim()
  print(top-bot)
  dh = 50.0 * fscale/1e6
  plt.errorbar(1.04, 0.9*top +0.1*bot-dh, dh, fmt="k")
  if j == 4:
    plt.ylabel(r"${\rm Transit\ depth,}\ (R_{\rm p}/R_{\rm s})^2\ (\%)$",
               fontsize=fs)
plt.xlabel(r"$\rm Wavelength\ \ (um)$", fontsize=fs)
plt.savefig("plots/transmission.ps")


