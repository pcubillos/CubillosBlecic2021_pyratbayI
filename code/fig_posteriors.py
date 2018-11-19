#! /usr/bin/env python

import sys
import os
import copy
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from scipy.ndimage.filters import gaussian_filter1d as gaussf
import ConfigParser as configparser
plt.ioff()

sys.path.append("../pyratbay")
import pyratbay as pb
import pyratbay.constants  as pc
import pyratbay.atmosphere as pa

sys.path.append('../pyratbay/modules/MCcubed/')
import MCcubed.plots as mp


def allocateaxes(allfreep, pnames, hrect, margin):
  nfree = len(allfreep)
  nruns = len(pnames)
  # Allocate histograms:
  allaxes = []
  axes = []
  for j in np.arange(nruns):
    axes.append([])
    allaxes.append([])
    for k in np.arange(nfree):
      ax = mp.subplotter(hrect, margin, 1+j*nfree+k, nfree, nruns,
                         ymargin=3*margin)
      if allfreep[k] in pnames[j]:
        axes[j].append(ax)
      else:
        axclear(ax)
      allaxes[j].append(ax)
  return allaxes, axes


def axclear(ax):
  ax.set_xticks([])
  ax.set_yticks([])
  ax.get_xaxis().set_visible(False)
  ax.get_yaxis().set_visible(False)
  plt.setp(ax.spines.values(), color='0.7')
  plt.setp([ax.get_xticklines(), ax.get_yticklines()], color='r')


cfiles = [
  # HATP11b
  ["../run02_HATP11b/mcmc_hatp11b_w0-00000-c.cfg",
   "../run02_HATP11b/mcmc_hatp11b_wm-00000-0.cfg",
   "../run02_HATP11b/mcmc_hatp11b_wm-00000-c.cfg"],
  # HATP32b
  ["../run03_HATP32b/mcmc_hatp32b_w0-00000-c.cfg",
   "../run03_HATP32b/mcmc_hatp32b_wm-00000-0.cfg",
   "../run03_HATP32b/mcmc_hatp32b_wm-00000-c.cfg"],
  # HATP38b
  ["../run04_HATP38b/mcmc_hatp38b_w0-00000-c.cfg",
   "../run04_HATP38b/mcmc_hatp38b_wm-00000-0.cfg",
   "../run04_HATP38b/mcmc_hatp38b_wm-00000-c.cfg"],
  # WASP043b
 ["../run05_WASP043b/mcmc_wasp043b_w0-00000-c.cfg",
  "../run05_WASP043b/mcmc_wasp043b_wm-00000-0.cfg",
  "../run05_WASP043b/mcmc_wasp043b_wm-00000-c.cfg",
  "../run05_WASP043b/mcmc_wasp043b_w0-000h0-c.cfg"],
  # WASP063b
 ["../run06_WASP063b/mcmc_wasp063b_w0-00000-c.cfg",
  "../run06_WASP063b/mcmc_wasp063b_wm-00000-0.cfg",
  "../run06_WASP063b/mcmc_wasp063b_wm-00000-c.cfg",
  "../run06_WASP063b/mcmc_wasp063b_w0-000h0-c.cfg",
  "../run06_WASP063b/mcmc_wasp063b_wm-000h0-0.cfg",
  "../run06_WASP063b/mcmc_wasp063b_wm-000h0-c.cfg"],
  # WASP067b
 ["../run07_WASP067b/mcmc_wasp067b_w0-00000-c.cfg",
  "../run07_WASP067b/mcmc_wasp067b_wm-00000-0.cfg",
  "../run07_WASP067b/mcmc_wasp067b_wm-00000-c.cfg"],
  # WASP101b
 ["../run08_WASP101b/mcmc_wasp101b_w0-00000-c.cfg",
  "../run08_WASP101b/mcmc_wasp101b_wm-00000-0.cfg",
  "../run08_WASP101b/mcmc_wasp101b_wm-00000-c.cfg"],
  # WASP107b
 ["../run09_WASP107b/mcmc_wasp107b_w0-00000-c.cfg",
  "../run09_WASP107b/mcmc_wasp107b_wm-00000-0.cfg",
  "../run09_WASP107b/mcmc_wasp107b_wm-00000-c.cfg"],
 ]

# Metallicity (by mass):
# X + Y + Z = 1.0
# X = N_H  * m_H
# Y = N_He * m_He
# Z = sum(N_Z*m_z) = mu - X - Y
Zsun = 0.0134
Xsun = 0.7381
# Indices for species in atmosphere:
spec, p, t, q = pa.readatm("../run01/isothermal_1500K_uniform.atm")
iHe   = np.where(spec == 'He')[0][0]
iH2   = np.where(spec == 'H2')[0][0]
iH    = np.where(spec == 'H')[0][0]
iH2O  = np.where(spec == 'H2O')[0][0]
iCH4  = np.where(spec == 'CH4')[0][0]
iHCN  = np.where(spec == 'HCN')[0][0]
iNH3  = np.where(spec == 'NH3')[0][0]
iC2H2 = np.where(spec == 'C2H2')[0][0]
iC2H4 = np.where(spec == 'C2H4')[0][0]

# This index sets the planet:
i = 2

nruns = len(cfiles[i])
nparams = 12
freep = np.zeros(nparams, bool)
posterior, pnames = [], []

# Read MCMC results:
for j in np.arange(nruns):
  path = os.path.dirname(cfiles[i][j]) + "/"
  config = configparser.SafeConfigParser()
  config.read([cfiles[i][j]])
  # Burned-in posterior sample:
  rootname = config.get('pyrat','logfile').replace('.log','')
  planet = rootname.split("_")[1]
  savefile = rootname + '.npz'
  burn = (config.getint('pyrat','burnin') * config.getint('pyrat','nchains')
         / config.getint('pyrat','thinning'))
  post = np.load(path+savefile)["Z"]
  posterior.append(post[burn:])
  # Open Pyrat Object:
  log = open("dummy.log", "w")  # Avoid overwriting the log file
  pyrat = pb.pyrat.init(cfiles[i][j], log=log)
  spec  = pyrat.mol.name
  q     = pyrat.atm.q
  freep |= pyrat.ret.stepsize > 0

  posterior[j][:,pyrat.ret.irad] *= pc.km/pc.rjup
  parname = np.asarray(pyrat.ret.parname, "|S100")
  pnames.append(parname[pyrat.ret.stepsize>0])
  ifree = np.where(pyrat.ret.stepsize>0)[0]
  if 'N2' in np.array(pyrat.ret.molscale)[np.in1d(pyrat.ret.iabund, ifree)]:
    # Unique posterior values:
    u, uind, uinv = np.unique(posterior[j][:,0], return_index=True,
                              return_inverse=True)
    nunique = np.size(u)
    # Get mean molecular mass:
    mu = np.zeros(nunique, np.double)
    X  = np.zeros(nunique, np.double)
    Z  = np.zeros(nunique, np.double)
    params = pyrat.ret.params
    imol = np.intersect1d(ifree, pyrat.ret.iabund)
    for k in np.arange(nunique):
      idx = uind[k]  # Sample index
      params[imol] = posterior[j][idx, np.in1d(ifree, pyrat.ret.iabund)]
      q2 = pa.qscale(q, spec, params[pyrat.ret.iabund],
            pyrat.ret.molscale, bulk=pyrat.ret.bulk, iscale=pyrat.ret.iscale,
            ibulk=pyrat.ret.ibulk, bratio=pyrat.ret.bulkratio,
            invsrat=pyrat.ret.invsrat)
      mu[k] = pa.meanweight(q2, spec)[0]
      Y    =    q2[0,iHe]  *pyrat.mol.mass[iHe]
      X[k] = (2*q2[0,iH2]  *pyrat.mol.mass[iH]
             +  q2[0,iH]   *pyrat.mol.mass[iH]
             +4*q2[0,iH2O] *pyrat.mol.mass[iH]
             +  q2[0,iHCN] *pyrat.mol.mass[iH]
             +3*q2[0,iCH4] *pyrat.mol.mass[iH]
             +2*q2[0,iC2H2]*pyrat.mol.mass[iH]
             +4*q2[0,iC2H4]*pyrat.mol.mass[iH])
      Z[k] = mu[k] - X[k] - Y  # X + Y + Z = mu
      if (k+1) % (nunique/10) == 0:
        print("{:3.0f}% done.".format(k*100.0/nunique))
    mmm = mu[uinv]
    # Metallicity by mass with respect to solar values, i.e., [M/H]:
    MH = np.log10(Z/X / (Zsun/Xsun))
    MH = MH[uinv]

    iN2 = np.where(ifree==pyrat.ret.iabund[pyrat.ret.molscale.index('N2')])
    iN2 = iN2[0][0]
    posterior[j][:,iN2] = MH

# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

nfree = np.sum(freep)
allfreep = parname[freep]
# Adjust parameter names:
for j in np.arange(nruns):
  for k in np.arange(len(pnames[j])):
    if pnames[j][k].find("Radius")>=0:
      pnames[j][k] = r"$R_0\ (R_{\rm Jup})$"
    if pnames[j][k].find("H2O")>=0:
      pnames[j][k] = r"$\log(X_{\rm H2O})$"
    if pnames[j][k].find("HCN")>=0:
      pnames[j][k] = r"$\log(X_{\rm HCN})$"
    if pnames[j][k].find("N2")>=0:
      pnames[j][k] = r"$[\rm M/\rm H]$"
    if pnames[j][k].find("top")>=0:
      pnames[j][k] = r"$\log(p_{\rm top})\ (\rm bar)$"

for k in np.arange(nfree):
  if allfreep[k].find("Radius")>=0:
    allfreep[k] = r"$R_0\ (R_{\rm Jup})$"
  if allfreep[k].find("H2O")>=0:
    allfreep[k] = r"$\log(X_{\rm H2O})$"
  if allfreep[k].find("HCN")>=0:
    allfreep[k] = r"$\log(X_{\rm HCN})$"
  if allfreep[k].find("N2")>=0:
    allfreep[k] = r"$[\rm M/\rm H]$"
  if allfreep[k].find("top")>=0:
    allfreep[k] = r"$\log(p_{\rm top})\ (\rm bar)$"

# Swap positions of HCN and [H/M]:
if i in [3,4]:
  hcnsort = [0,1,2,4,3,5]
  allfreep = allfreep[hcnsort]
  for j in np.arange(nruns):
    if r"$[\rm M/\rm H]$" in pnames[j] and r"$\log(X_{\rm HCN})$" in pnames[j]:
      tmp = np.copy(posterior[j][:,3])
      posterior[j][:,3] = np.copy(posterior[j][:,4])
      posterior[j][:,4] = tmp


ranges = [
 # Temp        R0              H2O        HCN     [M/H]     pc
 [[ 300, 1000], [0.415,0.445], [-6.5, 0],         [-2,  4.0], [-5,2]],   # H11
 [[ 300, 2750], [1.68, 1.86],  [-6.2, 0],         [-2,  3.0], [-6,2]],   # H32
 [[ 300, 1500], [0.79, 0.88],  [-9.0, 0],         [-2,  3.0], [-6,2]],   # H38
 [[ 300, 1500], [1.048,1.062], [-6.5, 0], [-9,0], [-2,  3.5], [-5.5,2]], # W43
 [[ 300, 3000], [1.32, 1.48],  [-7.0, 0], [-7,0], [-2, 3.5],  [-6.5,2]], # W63
 [[ 300, 3000], [1.22, 1.42],  [-9.0, 0],         [-1.8,3.8], [-8,2]],   # W67
 [[ 300, 3000], [1.12, 1.43],  [-9.0, 0],         [-2,  4.0], [-8,2]],   # W101
 [[ 300, 1000], [0.915, 0.947],[-6.5, 0],         [ 2,  4.0], [-6,2]],   # W107
]

ticks = [
  # HAT11b
  [[300, 500, 700, 900], [0.42, 0.43, 0.44], [-6, -4, -2, 0.0],
   [-1, 1, 3], [-4, -2, 0, 2]],
  # HAT32b
  [[500, 1000, 1500, 2000, 2500], [1.7, 1.75, 1.8, 1.85], [-6,-4,-2,0],
   [-1, 1, 3], [-4, -2, 0, 2]],
  # HAT38b
  [[300, 600, 900, 1200, 1500], [0.8, 0.82, 0.84, 0.86, 0.88], [-8,-6,-4,-2,0],
   [-1, 0, 1, 2], [-6, -4, -2, 0, 2]],
  # WASP043
  [[ 300, 700, 1100, 1500], [1.05, 1.055, 1.06], [-6, -3, 0],
   [-9, -6, -3, -0], [-1, 1, 3], [-4, -2, 0, 2]],
  # WASP063
  [[ 300, 1200, 2100, 3000], [1.35, 1.4, 1.45], [-6, -3, 0],
   [-6, -3, -0], [-1, 1, 3], [-5.5, -3, -0.5, 2]],
  # WASP067
  [[300, 1200, 2100, 3000], [1.25, 1.3, 1.35, 1.4], [-8, -6, -4, -2, 0.0],
   [-1, 1, 3], [-7, -4, -1, 2]],
  # WASP101
  [[300, 1200, 2100, 3000], [1.2, 1.3, 1.4], [-8, -6, -4, -2, 0.0],
   [-1, 1, 3], [-7, -4, -1, 2]],
  # WASP107
  [[300, 500, 700, 900], [0.92, 0.93, 0.94], [-6, -4, -2, 0],
   [2.5, 3.0, 3.5], [-6, -4, -2, 0, 2]],
]

# Modify ranges and ticks for pairwise plots:
pranges = copy.deepcopy(ranges)
pticks  = copy.deepcopy(ticks)
if i == 5:
  pranges[i][1] = 1.325, 1.417
  pticks [i][1] = 1.34, 1.36, 1.38, 1.4
elif i == 6:
  pranges[i][1] = 1.355, 1.4
  pticks [i][1] = 1.36, 1.37, 1.38, 1.39, 1.4
elif i == 7:
  pranges[i][2] = -0.7, -0.1
  pticks [i][2] = -0.6, -0.4, -0.2


# OK for 2 planets per page:
#hrect = [0.06, 0.23, 0.52, 0.87]
#prect = [0.60, 0.2,  0.97, 0.9]
#margin = 0.015

# 3 planets per page:
#        left  bottom right top
hrect = [0.06, 0.33,  0.52, 0.82]
prect = [0.63, 0.3,   0.97, 0.85]
margin = 0.01
fs = 9

if   nruns == 4:
  hrect[1], hrect[3] = 0.27, 0.88
elif nruns >= 5:
  hrect[1], hrect[3] = 0.21, 0.94

# Index of model to select for pairwise:
jbest = [0, 0, 0, 0, 3, 1, 1, 0]


plt.figure(0, (8.5, 5))
plt.clf()
# Background:
ax0 = plt.axes([0.0, 0.0, 1.0, 1.0])
ax0.set_xlim(0,1)
ax0.set_ylim(0,1)
ax0.set_axis_off()
ax0.text(hrect[0], hrect[3]+margin, planet, fontsize=fs+2, va='bottom')
ax0.text(0.005, 0.5*(hrect[1]+hrect[3]), r"$\rm Posterior\ density$",
         rotation=90, va='center', fontsize=fs+1)
# Allocate histograms:
allaxes, axes = allocateaxes(allfreep, pnames, hrect, margin)

pos = allaxes[jbest[i]][nfree-1].get_position().get_points()
bot, top = pos[0][1], pos[1][1]
# lelo leup riup, rilo
verts = [(hrect[2]+margin, bot-margin)] + [(hrect[2]+margin, top+margin)] + \
        [(prect[0], prect[3])] + [(prect[0], prect[1])]
poly = matplotlib.patches.Polygon(verts, facecolor='0.85',
            edgecolor='none')
ax0.add_patch(poly)
hverts = [(hrect[0]-margin,bot-margin), (hrect[0]-margin,top+margin),
          (hrect[2]+margin,top+margin), (hrect[2]+margin,bot-margin)]
ax0.add_patch(matplotlib.patches.Polygon(hverts, facecolor='0.85',
            edgecolor='none', zorder=5))
# The posteriors:
for j in np.arange(nruns):
  ifree = np.in1d(allfreep, pnames[j])
  mp.histogram(posterior[j], pnames[j], percentile=0.683,
      lw=1.25, fs=fs, axes=axes[j], ranges=np.array(ranges[i])[ifree])
  for k in np.arange(nfree):
    allaxes[j][k].set_yticklabels([])
    if j == nruns-1:
      allaxes[j][k].set_xlabel(allfreep[k])
    else:
      allaxes[j][k].set_xticklabels([])
      allaxes[j][k].set_xlabel("")
    allaxes[j][k].set_xlim(ranges[i][k])
ifree = np.in1d(allfreep, pnames[jbest[i]])
axp = mp.pairwise(posterior[jbest[i]], pnames[jbest[i]], fs=fs,
         rect=prect, margin=1.5*margin, ranges=np.array(pranges[i])[ifree])
# Set ticks:
for j in np.arange(nruns):
  for k in np.arange(nfree):
    allaxes[j][k].set_xticks(ticks[i][k])
    if j == nruns-1:
      allaxes[j][k].get_xaxis().set_visible(True)
      allaxes[j][k].tick_params(labelsize=fs-1)
      plt.setp(allaxes[j][k].xaxis.get_majorticklabels(), rotation=90)
  allaxes[j][0].set_ylabel(r"$\rm M{:d}$".format(j+1), fontsize=fs)
nbfree = len(axes[jbest[i]])
for k in np.arange(nbfree-1):
    for l in np.arange(k, nbfree-1):
      axp[k][l].set_xticks(pticks[i][np.where(ifree)[0][k]])
      axp[k][l].set_yticks(pticks[i][np.where(ifree)[0][l+1]])
plt.savefig("../plots/posterior_{:s}.ps".format(planet))
plt.savefig("../plots/posterior_{:s}.pdf".format(planet))
