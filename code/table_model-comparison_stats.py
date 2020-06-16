#! /usr/bin/env python
import os
import numpy as np

import pyratbay as pb
from pyratbay.tools import cd
import mc3

"""
This horrible piece of Python script, takes info from the MCMC
config and log files and makes a summary table.
"""

def sfit(params, data=None, uncert=None, res=False):
    """Fit a flat curve to a dataset."""
    if res:
        return (data-params[0])/uncert
    return np.tile(params[0], len(data))

# Data files:
targets = [
    'run_GJ-0436b_tsiaras',
    'run_HATP-03b_tsiaras',
    'run_HATP-17b_tsiaras',
    'run_HATP-18b_tsiaras',
    'run_HATP-32b_damiano',
    'run_HATP-32b_tsiaras',
    'run_HATP-38b_bruno',
    'run_HATP-38b_tsiaras',
    'run_HATP-41b_tsiaras',
    'run_HD-149026b_tsiaras',
    'run_K2-018b_benneke',
    'run_K2-018b_tsiaras',
    'run_WASP-029b_tsiaras',
    'run_WASP-043b_stevenson',
    'run_WASP-043b_tsiaras',
    'run_WASP-063b_kilpatrick',
    'run_WASP-063b_tsiaras',
    'run_WASP-067b_bruno',
    'run_WASP-067b_tsiaras',
    'run_WASP-069b_tsiaras',
    'run_WASP-074b_tsiaras',
    'run_WASP-080b_tsiaras',
    'run_WASP-101b_tsiaras',
    'run_WASP-101b_wakeford',
    'run_WASP-103b_kreidberg',
    'run_XO-1b_tsiaras',
    ]

ntargets = len(targets)

# Extract values:
names = [None] * ntargets
flat_chisq = np.zeros(ntargets)
flat_rchisq = np.zeros(ntargets)
flat_bic = np.zeros(ntargets)

for i in range(ntargets):
    cfile = [f for f in os.listdir(targets[i]) if f.endswith(".cfg")][0]
    with cd(targets[i]):
        pyrat = pb.run(cfile, init=True, no_logfile=True)

    data   = pyrat.obs.data
    uncert = pyrat.obs.uncert
    indparams = [data]
    params = np.array([np.mean(data)])
    pstep = np.array([1.0])
    mc3_fit = mc3.fit(data, uncert, sfit, params, indparams, pstep)

    dof = len(data) - len(params)
    names[i] = cfile.split('_')[1] + ' ' + cfile.split('_')[2]
    flat_chisq[i] = mc3_fit['best_chisq']
    flat_rchisq[i] = mc3_fit['best_chisq']/dof
    flat_bic[i] = mc3_fit['best_chisq'] + len(params) * np.log(len(data))

with open("stats/flat_fit.dat", "w") as f:
    f.write(
        "Flat-curve fit to transmission data (Nfree = 1):\n\n"
        "Planet                chi-square  red-chisq      BIC\n"
        "----------------------------------------------------\n")
    for i in range(ntargets):
        f.write(f"{names[i]:20s}  {flat_chisq[i]:10.3f}  "
                f"{flat_rchisq[i]:9.3f}  {flat_bic[i]:7.3f}\n")



table = ''
for j in range(ntargets):
    cfiles = [f for f in os.listdir(targets[j]) if f.endswith(".cfg")]
    nruns = len(cfiles)
    bic = np.zeros(nruns+1)
    chisq = np.zeros(nruns+1)
    red_chisq = np.zeros(nruns+1)
    nfree = np.zeros(nruns+1, int)
    idx   = np.zeros(nruns+1, int)
    pars  = np.zeros(nruns+1, "|U240")
    for i in np.arange(nruns):
        with cd(targets[j]):
            pyrat = pb.run(cfiles[i], True, True)
        with np.load(pyrat.ret.mcmcfile) as mcmc:
            bic[i] = mcmc['BIC']
            chisq[i] = mcmc['best_chisq']
            red_chisq[i] = mcmc['red_chisq']

        nfree[i] = np.sum(pyrat.ret.pstep > 0)

        key = cfiles[i][-14:-4]
        pars[i] = "$T$, $R\sb{\\rm p}$, {\\water}"
        if 'log(CO)' in pyrat.ret.pnames:
            pars[i] += ", CO"
            idx[i] += 3
        if 'log(CO2)' in pyrat.ret.pnames:
            pars[i] += ", {\\carbdiox}"
            idx[i] += 3
        if 'log(CH4)' in pyrat.ret.pnames:
            pars[i] += ", {\\methane}"
            idx[i] += 3
        if 'log(HCN)' in pyrat.ret.pnames:
            pars[i] += ", HCN"
            idx[i] += 3
        if 'log(N2)' in pyrat.ret.pnames:
            pars[i] += ", [M/H]"
            idx[i] += 2
        if 'log(p_top)' in pyrat.ret.pnames:
            pars[i] += ", {\\pcloud}"
            idx[i] += 1

    ndata = len(pyrat.obs.data)
    red_chisq[-1] = flat_rchisq[j]
    bic     [-1] = flat_bic[j]
    nfree   [-1] = 1
    pars    [-1] = "Flat fit"
    idx[nruns] = nruns + 1
    delta_bic = bic - np.amin(bic)
    # Uncertainty in reduced chi squared:
    chi_unc = np.sqrt(2.0/(ndata-nfree))

    file_name = '_'.join(names[j].split())
    for i in range(nruns+1):
        pars[i] = (
            f"M{idx[i]%(nruns+1)} & {pars[i]:50s} & {nfree[i]} & "
            f"${delta_bic[i]:4.1f}$ & "
            f"${red_chisq[i]:.2f}\\pm{chi_unc[i]:.2f}$ \\\\ \n")

    table += (
        f"\multicolumn{{5}}{{l}}{{ {names[j]} }}\\\\ \n" +
        "\\hline\n"
        "ID & Free parameters & $N\sb{\\rm free}$ & $\\Delta$BIC & \\redchisq \\\\ \n"
        "\\hline \n")
    for i in range(nruns+1):
        table += pars[np.argsort(idx)[i]]
    table += "\\hline\n\\\\ \n"

with open(f"stats/sample_stats.tex", "w") as f:
    f.write(
        "{\\renewcommand{\\arraystretch}{1.0}\n"
        "\\setlength{\\tabcolsep}{3pt}\n"
        "\\begin{table}\n"
        "\\centering\n" +
        "\\caption{Model comparison statistics.}\n" +
        "\\label{table:model_comparison}\n" +
        "\\begin{tabular}{@{\\extracolsep{\\fill}}c@{\\hskip 6pt}lccc}\n"
        )
    f.write(table)
    f.write(
        "\\hline\n"
        "\\end{tabular}\n"
        "\\end{table}\n"
        "\\normalsize\n"
        "}\n")


