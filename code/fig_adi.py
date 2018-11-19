#! /usr/bin/env python

import sys
import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import ConfigParser as configparser
plt.ioff()


# Data files:
folders = [
  "../run02_HATP11b/",
  "../run03_HATP32b/",
  "../run04_HATP38b/",
  "../run05_WASP043b/",
  "../run06_WASP063b/",
  "../run07_WASP067b/",
  "../run08_WASP101b/",
  "../run09_WASP107b/",
  ]
nplanets = len(folders)

# Tsiaras 2018:                                                #dummy
snr_t = np.array([7.62, 14.32, 5.47, 7.34, 12.22, 5.87, 14.03, 1.0])
adi_t = np.array([6.61, 16.44, 0.67, 1.93,  0.00, 0.27,  0.00, -1.0])


# Read flatfit results:
with open("flatfit.dat", "r") as f:
  lines = f.readlines()
for i in np.arange(len(lines)):
  if lines[i].startswith("---"):
    break
# Extract values:
name = np.zeros(nplanets, "|S15")
flatchisq = np.zeros(nplanets)
flatbic   = np.zeros(nplanets)
for j in np.arange(nplanets):
  info = lines[i+j+1].split()
  name[j], flatchisq[j], flatbic[j] = info[0], info[2], info[3]


bic      = np.tile(np.inf, nplanets)
snr_pb   = np.zeros(nplanets)
# Read PB MCMC results:
cwd = os.getcwd() + "/"
for j in np.arange(nplanets):
  files = os.listdir(cwd + folders[j])
  cfiles = []
  for f in sorted(files)[::-1]:
    if f.startswith("mcmc") and f.endswith(".cfg"):
      config = configparser.SafeConfigParser()
      config.read(cwd + folders[j] + f)
      data  = np.array(config.get('pyrat', 'data').split(), float)
      error = np.array(config.get('pyrat', 'uncert').split(), float)
      if os.path.exists(cwd + folders[j] + config.get('pyrat','logfile')):
        cfiles.append(f)

  nruns = len(cfiles)
  for i in np.arange(nruns):
    config = configparser.SafeConfigParser()
    config.read(cwd + folders[j] + cfiles[i])
    # Extract info:
    d = np.load(cwd + folders[j] +
                config.get('pyrat','logfile').replace("log","npz"))
    try:
      if d['BIC'] < bic[j]:
        bic[j]      = d['BIC']
        snr_pb[j]   = np.median(data/error)
    except:
      pass
    d.close()

adi_pb = 0.5*np.clip(flatbic-bic, 0, np.inf)


# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Plot

fs = 12
ms = 6
dy1 = np.tile(1.04, nplanets)
dy2 = np.tile(1.04, nplanets)
dy1[1] = dy2[2] = dy2[6] = 0.9

label = ["HAT-P-11b", "HAT-P-32b", "HAT-P-38b", "WASP-43b",
         "WASP-63b",  "WASP-67b",  "WASP-101b", "WASP-107b"]

col = ['navy', 'orangered']

dadi = adi_pb - adi_t
dadi[-1] = 0.0
isort = np.argsort(dadi)

lab = np.asarray(label)[isort]
adit  = adi_t [isort]
adipb = adi_pb[isort]


plt.figure(1, (6.5, 4.25))
plt.clf()
plt.subplots_adjust(0.06, 0.12, 0.98, 0.97)
ax = plt.subplot(111)
for i in np.arange(nplanets):
  if adit[i] >= 0:
    tsi = plt.plot(adit[i],i, "D", color="w", ms=ms+1, mec=col[0], mew=1.5,
                   label='Tsiaras et al. (2018)')
  pbi = plt.plot(adipb[i], i, "o", color=col[1], ms=ms, label='This work')
  plt.text(-1.28, i, lab[i], ha='left', va='center', fontsize=fs)
ax.set_xscale('symlog')
plt.xlim(-1.32, 60)
plt.ylim(-0.5, 9.2)
plt.gca().xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
ax.set_xticks([0,1,3, 10, 20, 50])
ax.set_yticklabels([])
ax.set_yticks(np.arange(nplanets))
ax.tick_params(labelsize=fs-1)

ax.legend(handles=tsi+pbi, loc='upper right', fontsize=fs)
ax.set_xlabel(r"ADI = 0.5$\,\Delta$BIC", fontsize=fs)
ax.set_ylabel("Planet", fontsize=fs)

plt.axvline( 0, dashes=(7,2,2,2), lw=1.5, color="0.8", zorder=-10)
plt.axvline( 3, ls="--", lw=1.5, color="0.8", zorder=-10)
plt.axvline(11, ls="--", lw=1.5, color="0.8", zorder=-10)
plt.text( 3, 0, r"3$\sigma$", ha='right', va='center', rotation=90,
         color="0.7", fontsize=fs)
plt.text(11, 0, r"5$\sigma$", ha='right', va='center', rotation=90,
         color="0.7", fontsize=fs)
plt.savefig("../plots/adi.ps")
plt.savefig("../plots/adi.pdf")

