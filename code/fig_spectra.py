#! /usr/bin/env python
import sys
import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import scipy.constants   as sc
import ConfigParser as configparser
from scipy.ndimage.filters import gaussian_filter1d as gaussf
plt.ioff()

# Edit the path to the Pyrat-Bay package if necessary:
sys.path.append("../pyratbay")
import pyratbay as pb
import pyratbay.constants  as pc
import pyratbay.wine       as pw
import pyratbay.atmosphere as pa

# Load best-fitting models:
cfiles = [
   "../run02_HATP11b/mcmc_hatp11b_w0-00000-c.cfg",
   "../run03_HATP32b/mcmc_hatp32b_w0-00000-c.cfg",
   "../run04_HATP38b/mcmc_hatp38b_w0-00000-c.cfg",
   "../run05_WASP043b/mcmc_wasp043b_w0-00000-c.cfg",
   "../run06_WASP063b/mcmc_wasp063b_w0-000h0-c.cfg",
   "../run07_WASP067b/mcmc_wasp067b_wm-00000-0.cfg",
   "../run08_WASP101b/mcmc_wasp101b_wm-00000-0.cfg",
   "../run09_WASP107b/mcmc_wasp107b_w0-00000-c.cfg",
   ]

nplanets = len(cfiles)

# Get posterior spectra percentiles:
low2    = [None] * nplanets
low1    = [None] * nplanets
high1   = [None] * nplanets
high2   = [None] * nplanets
bmodel  = [None] * nplanets
data    = [None] * nplanets
uncert  = [None] * nplanets
filters = [None] * nplanets
name    = [None] * nplanets
bfpars  = [None] * nplanets
Rp = np.zeros(nplanets)
Rs = np.zeros(nplanets)
Mp = np.zeros(nplanets)
g  = np.zeros(nplanets)
Tmean  = np.zeros(nplanets)
mumean = np.zeros(nplanets)
Hm1  = np.zeros(nplanets)
Hm2  = np.zeros(nplanets)

for j in np.arange(nplanets):
  path = os.path.dirname(cfiles[j]) + "/"
  config = configparser.SafeConfigParser()
  config.read(cfiles[j])
  log = config.get('pyrat','logfile').strip('.log')
  name[j] = log.split("_")[1]
  # Model spectra:
  wl, bmodel[j], low2[j], low1[j], high1[j], high2[j] = \
    np.loadtxt(path + log + "_percentiles.dat", unpack=True)
  # Data:
  data[j]    = np.array(config.get('pyrat', 'data').split(), float)
  uncert[j]  = np.array(config.get('pyrat', 'uncert').split(), float)
  filters[j] = config.get('pyrat', 'filter').split()
  # Physics:
  rs = config.get('pyrat', 'rstar').split()
  Rs[j] = float(rs[0]) * getattr(pc, rs[1])
  mp = config.get('pyrat','mplanet').split()
  Mp[j] = float(mp[0]) * getattr(pc, mp[1]) # g
  Rp[j] = np.mean(np.sqrt(data[j])*Rs[j])   # cm
  g[j]  = (1e3*sc.G) * Mp[j] / Rp[j]**2     # cm s-2
  # Best-fitting parameters:
  bfpars[j] = np.load(path + log+".npz")["bestp"]
  # Get mean molecular mass:
  burn = (config.getint('pyrat','burnin')
         *config.getint('pyrat','nchains')
         /config.getint('pyrat','thinning'))
  post = np.load(path + log+".npz")["Z"][burn:]
  # Open Pyrat Object:
  nolog = open("dummy.log", "w")  # Avoid overwriting the log file
  pyrat = pb.pyrat.init(cfiles[j], log=nolog)
  spec  = pyrat.mol.name
  q     = pyrat.atm.q[0:1,:]
  Tmean[j] = np.mean(post[:,pyrat.ret.itemp])
  T, uind, uinv = np.unique(post[:,pyrat.ret.itemp], return_index=True,
                            return_inverse=True)
  nunique = np.size(T)
  mu = np.zeros(nunique, np.double)
  H  = np.zeros(nunique, np.double)
  params = pyrat.ret.params
  ifree = np.where(pyrat.ret.stepsize>0)[0]
  imol = np.intersect1d(ifree, pyrat.ret.iabund)
  pmol = post[:, np.in1d(ifree, pyrat.ret.iabund)]
  for k in np.arange(nunique):
    idx = uind[k]  # Sample index
    params[imol] = pmol[idx]
    q2 = pa.qscale(q[0:1,:], spec, params[pyrat.ret.iabund], pyrat.ret.molscale,
          bulk=pyrat.ret.bulk, iscale=pyrat.ret.iscale,
          ibulk=pyrat.ret.ibulk, bratio=pyrat.ret.bulkratio[0:1,:],
          invsrat=pyrat.ret.invsrat[0:1])
    mu[k] = pa.meanweight(q2, spec)
    H[k] = 100 * sc.k * T[k] / (mu[k]*(g[j]/100.0)) # cm
    if k % (nunique/10) == 0:
      print("{}/{}".format(k, nunique))
  mumean[j] = np.mean(mu[uinv]) * sc.atomic_mass
  Hm1[j] = np.mean(H[uinv])/sc.atomic_mass
  Hm2[j] = 100 * sc.k * Tmean[j] / (mumean[j]*(g[j]/100.0))

molscale = config.get('pyrat', 'molscale').split()
wn = 1.0/(wl*pc.um)

# Filter transmission:
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


# Equilibrium temperature:
Teq = np.array([880, 1790, 1080, 1440, 1500, 1040, 1560, 770], float)
# Solar mean molecular mass:
musun = 2.3 * sc.atomic_mass # kg
# Equilibrium scale height at solar abundance:
H = 100 * sc.k * Teq / (musun*(g/100.0)) # cm
dHeq = ((Rp+H)/Rs)**2 - (Rp/Rs)**2  # Scale height in depth units
# Best-fitting temperature:
Tbf = np.array(bfpars)[:,0]
# Best-fitting mean molecular mass:
spec, press, temp, q0 = pa.readatm("../run01/isothermal_1500K_uniform.atm")
mubf = np.zeros(nplanets)
for j in np.arange(nplanets):
  q = pa.qscale(np.copy(q0), spec, bfpars[j][2:11], molscale, ["H2","He"])
  mubf[j] = pa.meanweight(q, spec)[0] * sc.atomic_mass
# Best-fit scale height:
Hf = 100 * sc.k * Tbf / (mubf*(g/100.0)) # cm
dHf = ((Rp+Hf)/Rs)**2 - (Rp/Rs)**2  # Scale height in depth units
# Mean mean molecular mass:
dHm = ((Rp+Hm1)/Rs)**2 - (Rp/Rs)**2



# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Plot it:
lw = 1.2
ms = 4.5
fs = 13
cs = 1.3
sigma = 4.0
ye = np.arange(1.02, 1.5, 0.02)
xticks = [1.0, 1.3, 1.6, 2.0, 3.0, 4.0, 5.0]
f = 100

fmt = ["limegreen", "k", "0.6", "crimson"]

yran = [
  [f*0.00328,f*0.00372],  # H11
  [f*0.02240,f*0.0244],   # H32
  [f*0.00830,f*0.0093],   # H38
  [f*0.02522,f*0.02568],  # W43
  #[f*0.02485,f*0.0257],   # W43
  [f*0.00588,f*0.00685],  # W63
  [f*0.025  ,f*0.0275],   # W67
  [f*0.0115 ,f*0.01215],  # W101
  [f*0.0203 ,f*0.0211],   # W107
]

plt.figure(0, (8.5, 7))
plt.clf()
plt.subplots_adjust(0.1, 0.1, 0.99, 0.97, hspace=0.08, wspace=0.15)
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
  plt.gca().xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
  if j/2 == 3:
    plt.xlabel(r"$\rm Wavelength\ \ (um)$", fontsize=fs)
  else:
    ax.set_xticklabels([])
  if j == 4:
    plt.ylabel(25*" "+r"${\rm Transit\ depth,}\ (R_{\rm p}/R_{\rm s})^2\ (\%)$",
               fontsize=fs)
  ax.set_ylim(yran[j])
  ax.locator_params(axis='y', nbins=7)
  plt.text(0.96, 0.85, name[j], horizontalalignment='right',
           transform=ax.transAxes, fontsize=fs-2)
  # 100 ppm, Mean errorbar, Equilibrium scale height, Mean scale height:
  dh = np.array([100/1e6, np.median(uncert[j]), dHeq[j], dHm[j]]) * f/2
  for i in np.arange(len(dh)):
    bot, top = ax.get_ylim()
    dy = 0.96*top + 0.04*bot - dh[i]
    plt.errorbar(ye[i], dy, dh[i], fmt=fmt[i], lw=lw, capthick=lw, capsize=cs)
  ax2 = ax.twinx()
  ax2.tick_params(right='on', direction='in')
  ax2.set_yticklabels([])
  ax2.set_ylim(yran[j])
  ax2.locator_params(axis='y', nbins=7)
plt.savefig("../plots/transmission.ps")
