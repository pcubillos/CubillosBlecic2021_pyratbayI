import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
plt.ioff()


# Taken from stats/sample_stats.tex:
pyrat = [
    'run_GJ-0436b_tsiaras/MCMC_GJ-0436b_tsiaras_1.0-2.0um_w00000-0c.npz',
    'run_HATP-03b_tsiaras/MCMC_HATP-03b_tsiaras_1.0-2.0um_w00000-0c.npz',
    'run_HATP-17b_tsiaras/MCMC_HATP-17b_tsiaras_1.0-2.0um_w00000-0c.npz',
    'run_HATP-18b_tsiaras/MCMC_HATP-18b_tsiaras_1.0-2.0um_w00000-m0.npz',
    'run_HATP-32b_tsiaras/MCMC_HATP-32b_tsiaras_1.0-2.0um_w00000-0c.npz',
    'run_HATP-38b_tsiaras/MCMC_HATP-38b_tsiaras_1.0-2.0um_w00000-0c.npz',
    'run_HATP-41b_tsiaras/MCMC_HATP-41b_tsiaras_1.0-2.0um_w00000-0c.npz',
    'run_HD-149026b_tsiaras/MCMC_HD-149026b_tsiaras_1.0-2.0um_w00000-0c.npz',
    'run_K2-018b_tsiaras/MCMC_K2-018b_tsiaras_1.0-2.0um_w00000-m0.npz',
    'run_WASP-029b_tsiaras/MCMC_WASP-029b_tsiaras_1.0-2.0um_w00000-0c.npz',
    'run_WASP-043b_tsiaras/MCMC_WASP-043b_tsiaras_1.0-2.0um_w00000-0c.npz',
    'run_WASP-063b_tsiaras/MCMC_WASP-063b_tsiaras_1.0-2.0um_w00000-0c.npz',
    'run_WASP-067b_tsiaras/MCMC_WASP-067b_tsiaras_1.0-2.0um_w00000-0c.npz',
    'run_WASP-069b_tsiaras/MCMC_WASP-069b_tsiaras_1.0-2.0um_w00000-m0.npz',
    'run_WASP-074b_tsiaras/MCMC_WASP-074b_tsiaras_1.0-2.0um_w00000-m0.npz',
    'run_WASP-080b_tsiaras/MCMC_WASP-080b_tsiaras_1.0-2.0um_w00000-0c.npz',
    'run_WASP-101b_tsiaras/MCMC_WASP-101b_tsiaras_1.0-2.0um_w00000-m0.npz',
    None,
    'run_XO-1b_tsiaras/MCMC_XO-1b_tsiaras_1.0-2.0um_w00000-0c.npz',
    ]

other = [
    None,
    None,
    None,
    None,
    'run_HATP-32b_damiano/MCMC_HATP-32b_damiano_1.0-2.0um_w00000-0c.npz',
    'run_HATP-38b_bruno/MCMC_HATP-38b_bruno_1.0-2.0um_w00000-0c.npz',
    None,
    None,
    'run_K2-018b_benneke/MCMC_K2-018b_benneke_1.0-5.5um_w00000-0c.npz',
    None,
    'run_WASP-043b_stevenson/MCMC_WASP-043b_stevenson_1.0-5.5um_w00000-0c.npz',
    'run_WASP-063b_kilpatrick/MCMC_WASP-063b_kilpatrick_1.0-2.0um_w000h0-0c.npz',
    'run_WASP-067b_bruno/MCMC_WASP-067b_bruno_1.0-2.0um_w00000-0c.npz',
    None,
    None,
    None,
    'run_WASP-101b_wakeford/MCMC_WASP-101b_wakeford_1.0-2.0um_w00000-m0.npz',
    'run_WASP-103b_kreidberg/MCMC_WASP-103b_kreidberg_1.0-5.5um_w00000-0c.npz',
    None,
    ]


def main():
    nplanets = len(pyrat)
    names = [None] * nplanets
    for i in range(nplanets):
        cfg = pyrat[i] if pyrat[i] is not None else other[i]
        names[i] = os.path.split(cfg)[0].split('_')[1]

    # Read flatfit results:
    with open("stats/flat_fit.dat", "r") as f:
        lines = f.readlines()
    for i in range(len(lines)):
        if lines[i].startswith("---"):
            break

    # Extract values:
    flat_bic = np.zeros([nplanets,2])
    for j in range(i+1, len(lines)):
        info = lines[j].split()
        idx = list(names).index(info[0])
        ts = int(info[1] != 'tsiaras')
        flat_bic[idx, ts] = info[4]

    pyrat_bic = np.tile(-np.inf, nplanets)
    for i in range(nplanets):
        if pyrat[i] is None:
            continue
        with np.load(pyrat[i]) as mcmc:
            pyrat_bic[i] = mcmc['BIC']

    other_bic = np.tile(-np.inf, nplanets)
    for i in range(nplanets):
        if other[i] is None:
            continue
        with np.load(other[i]) as mcmc:
            other_bic[i] = mcmc['BIC']

    adi_pb = 0.5*np.clip(flat_bic[:,0]-pyrat_bic, 0, np.inf)
    adi_other = 0.5*np.clip(flat_bic[:,1]-other_bic, 0, np.inf)
    adi_ts = np.array([
        0.0,  0.0, 0.28,  5.71, 16.44, 0.67, 7.29, 0.0, 4.7, 1.25,
        1.93, 0.0, 0.27, 13.3,   0.0,  1.16, 0.0, np.inf,  3.15])


    # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # The plot:
    fs = 12
    ms = 6
    dy1 = np.tile(1.04, nplanets)
    dy2 = np.tile(1.04, nplanets)
    dy1[1] = dy2[2] = dy2[6] = 0.9
    col = ['navy', 'orangered', 'limegreen']


    plt.figure(1, (6.5, 6))
    plt.clf()
    plt.subplots_adjust(0.06, 0.08, 0.98, 0.97)
    ax = plt.subplot(111)
    for i in np.arange(nplanets):
        pbi = plt.plot(
            adi_pb[i], i, "o", color=col[1], ms=ms,
            label='This work / Tsiaras data')
        tsi = plt.plot(
            adi_ts[i], i, "D", color=col[0], ms=ms+1, mfc='none',
            mew=1.5, label='Tsiaras et al. (2018, 2019)')
        oi = plt.plot(
            adi_other[i], i, "+", color=col[2], ms=ms+1, mfc='none',
            mew=1.5, label="This work / Other's data")
        plt.text(-1.28, i, names[i], ha='left', va='center', fontsize=fs)
    ax.set_xscale('symlog')
    plt.xlim(-1.32, 25)
    plt.ylim(18.5, -2)
    plt.gca().xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax.set_xticks([0, 1, 3, 10, 20])
    ax.set_yticklabels([])
    ax.tick_params(labelsize=fs-1, direction='in', left=False)
    ax.tick_params(labelsize=fs)

    ax.legend(handles=tsi+pbi+oi, loc='upper right', fontsize=fs)
    ax.set_xlabel("ADI", fontsize=fs+3)
    ax.set_ylabel("Planet", fontsize=fs+3)

    plt.axvline( 0, dashes=(7,2,2,2), lw=1.5, color="0.8", zorder=-10)
    plt.axvline( 3, ls="--", lw=1.5, color="0.8", zorder=-10)
    plt.axvline(11, ls="--", lw=1.5, color="0.8", zorder=-10)
    plt.text(
        2.9, 17, r"3$\sigma$", ha='right', va='center', rotation=90,
        color="0.7", fontsize=fs)
    plt.text(
        11, 17, r"5$\sigma$", ha='right', va='center', rotation=90,
        color="0.7", fontsize=fs)
    for i in range(10):
        plt.axhline(2*i+0.5, c='0.95', lw=1.0)
    plt.savefig("plots/WFC3_sample_adi.pdf")


if __name__ == '__main__':
    main()
