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
sys.path.append("../pyratbay")
import pyratbay as pb
import pyratbay.constants  as pc
import pyratbay.tools      as pt
import pyratbay.wine       as pw
import pyratbay.starspec   as ps


# Load best-fitting models:
paths = ["../run02_HATP11b/",  "../run03_HATP32b/",  "../run04_HATP38b/",
         "../run05_WASP043b/", "../run06_WASP063b/", "../run07_WASP067b/",
         "../run08_WASP101b/", "../run09_WASP107b/"]
bestm =["MCMC_HAT-P-11b_1.0-5.5um_w0-00000-c_bestfit_spectrum.dat",
        "MCMC_HAT-P-32b_1.0-5.5um_w0-00000-c_bestfit_spectrum.dat",
        "MCMC_HAT-P-38b_1.0-5.5um_w0-00000-c_bestfit_spectrum.dat",
        "MCMC_WASP-043b_1.0-5.5um_w0-00000-c_bestfit_spectrum.dat",
        "MCMC_WASP-063b_1.0-5.5um_w0-000h0-c_bestfit_spectrum.dat",
        "MCMC_WASP-067b_1.0-5.5um_w0-00000-c_bestfit_spectrum.dat",
        "MCMC_WASP-101b_1.0-5.5um_wm-00000-0_bestfit_spectrum.dat",
        "MCMC_WASP-107b_1.0-5.5um_wm-00000-0_bestfit_spectrum.dat",
       ]

mcmcs = ["mcmc_hatp11b_w0-00000-c.cfg",
         "mcmc_hatp32b_w0-00000-c.cfg",
         "mcmc_hatp38b_w0-00000-c.cfg",
         "mcmc_wasp043b_w0-00000-c.cfg",
         "mcmc_wasp063b_w0-000h0-c.cfg",
         "mcmc_wasp067b_w0-00000-c.cfg",
         "mcmc_wasp101b_wm-00000-0.cfg",
         "mcmc_wasp107b_wm-00000-0.cfg",
        ]

nplanets = len(paths)
spec = []
for j in np.arange(nplanets):
  wn, fl = ps.readpyrat(paths[j] + bestm[j])
  spec.append(fl)
wl = 1.0/(wn*pc.um)

# Load data:
with open("filter_info.txt", "r") as f:
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

# Compute posterior spectra percentiles:
low2  = [None] * nplanets
low1  = [None] * nplanets
high1 = [None] * nplanets
high2 = [None] * nplanets
bmodel = [None] * nplanets
spec  = [None] * nplanets
#for i in np.arange(nplanets):
#  dummy, dummy, low2[i], low1[i], high1[i], high2[i] = \
#     pt.specpercent("../" + paths[i] + mcmcs[i], nmax=100000)

for i in np.arange(nplanets):
  dummy, bmodel[i], low2[i], low1[i], high1[i], high2[i] = \
    np.loadtxt(paths[i]+bestm[i].replace("bestfit_spectrum", "percentiles"),
               unpack=True)

bandwl = []
bandfl = []
for j in np.arange(nplanets):
  bwl = np.zeros(len(data[j]))
  bfl = np.zeros(len(data[j]))
  for i in np.arange(len(data[j])):
    fwn, ftr = pw.readfilter(filters[j][i])
    nifilter, wnidx, dummy = pw.resample(wn, fwn, ftr)
    bwl[i] = 1.0/(np.sum(fwn*ftr)/np.sum(ftr)*pc.um)
    bfl[i] = pw.bandintegrate(bmodel[j], wn, wnidx, nifilter)
  bandwl.append(bwl)
  bandfl.append(bfl)

# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Plot it:
lw = 1.25
ms = 4.5
fs = 14
sigma = 4.0
xticks = [1.0, 1.3, 1.6, 2.0, 3.0, 4.0, 5.0]
f = 100

yran = [None]*8
yran[0] = [f*0.00328,f*0.00372]  # H11
yran[1] = [f*0.02240,f*0.0244]   # H32
yran[2] = [f*0.00830,f*0.0096]   # H38
yran[3] = [f*0.02522,f*0.02568]  # W43
#yran[3] = [f*0.02485,f*0.0257]
yran[4] = [f*0.00588,f*0.00685]  # W63
yran[5] = [f*0.025  ,f*0.0275]   # W67
yran[6] = [f*0.0115 ,f*0.01215]  # W101
yran[7] = [f*0.0203 ,f*0.0211]   # W107

plt.figure(0, (8.5, 7))
plt.clf()
plt.subplots_adjust(0.1, 0.1, 0.99, 0.97, hspace=0.08, wspace=0.13)
for j in np.arange(nplanets):
  ax = plt.subplot(4,2,1+j)
  # Percentiles:
  ax.fill_between(wl, f*gaussf(low2[j],sigma), f*gaussf(high2[j],sigma),
    facecolor="gold", edgecolor="none", alpha=1.0)
  ax.fill_between(wl, f*gaussf(low1[j],sigma), f*gaussf(high1[j],sigma),
    facecolor="orange", edgecolor="none", alpha=1.0)
  # Best fit and data:
  plt.plot(wl, f*gaussf(bmodel[j], sigma), color='r', lw=lw)
  plt.plot(bandwl[j], f*bandfl[j], "o", color='orange', mew=0.15,
           mec="0.25", ms=ms)
  plt.errorbar(bandwl[j], f*data[j], uncert[j]*f, fmt="ob",
               label="Data", ms=ms, mew=0, elinewidth=lw, capthick=lw,
               capsize=0, zorder=3)
  plt.xlim(1.0, 5.5)
  plt.xscale('log')
  ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
  ax.set_xticks(xticks)
  if j/2 == 3:
    plt.xlabel(r"$\rm Wavelength\ \ (um)$", fontsize=fs)
  else:
    ax.set_xticklabels([])
  if j == 4:
    plt.ylabel(25*" "+r"${\rm Transit\ depth,}\ (R_{\rm p}/R_{\rm s})^2\ (\%)$",
               fontsize=fs)
  plt.ylim(yran[j])
  ax.locator_params(axis='y', nbins=7)
  plt.text(0.96, 0.85, name[j], horizontalalignment='right',
           transform=ax.transAxes, fontsize=fs-2)
  # 100 ppm:
  bot, top = ax.get_ylim()
  dh = 50.0 * f/1e6
  plt.errorbar(1.04, 0.95*top +0.05*bot-dh, dh, fmt="k", lw=lw, capthick=lw)
  # Mean errorbar:
  #plt.errorbar(1.04, 0.95*top +0.05*bot-dh, dh, fmt="k", lw=lw, capthick=lw)
  # Equilibrium scale height:
  #plt.errorbar(1.04, 0.95*top +0.05*bot-dh, dh, fmt="k", lw=lw, capthick=lw)
plt.savefig("../plots/transmission.ps")


