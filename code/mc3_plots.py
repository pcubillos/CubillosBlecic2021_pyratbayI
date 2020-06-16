# Copyright (c) 2015-2019 Patricio Cubillos and contributors.
# MC3 is open-source software under the MIT license (see LICENSE).

import os
import sys

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy.interpolate as si

from mc3 import utils as mu
from mc3 import stats as ms


if int(np.__version__.split('.')[1]) >= 15:
    histkeys = {'density':False}
else:
    histkeys = {'normed':False}

themes = {
    'blue':{'edgecolor':'blue', 'facecolor':'royalblue', 'color':'navy'},
    'red': {'edgecolor':'crimson', 'facecolor':'orangered', 'color':'darkred'},
    'black':{'edgecolor':'0.3', 'facecolor':'0.3', 'color':'black'},
    'green':{'edgecolor':'forestgreen', 'facecolor':'limegreen', 'color':'darkgreen'},
    'orange':{'edgecolor':'darkorange', 'facecolor':'gold', 'color':'darkgoldenrod'},
    }


def histogram(posterior, pnames=None, thinning=1, fignum=1100,
              savefile=None, bestp=None, quantile=None, pdf=None,
              xpdf=None, ranges=None, axes=None, lw=2.0, fs=11,
              yscale='absolute', theme='blue',
              # Deprecated: Remove by 2020-07-01
              percentile=None):
    """
    Plot parameter marginal posterior distributions

    Parameters
    ----------
    posterior: 1D or 2D float ndarray
        An MCMC posterior sampling with dimension [nsamples] or
        [nsamples, nparameters].
    pnames: Iterable (strings)
        Label names for parameters.
    thinning: Integer
        Thinning factor for plotting (plot every thinning-th value).
    fignum: Integer
        The figure number.
    savefile: Boolean
        If not None, name of file to save the plot.
    bestp: 1D float ndarray
        If not None, plot the best-fitting values for each parameter
        given by bestp.
    quantile: Float
        If not None, plot the quantile- highest posterior density region
        of the distribution.  For example, set quantile=0.68 for a 68% HPD.
    pdf: 1D float ndarray or list of ndarrays
        A smoothed PDF of the distribution for each parameter.
    xpdf: 1D float ndarray or list of ndarrays
        The X coordinates of the PDFs.
    ranges: List of 2-element arrays
        List with custom (lower,upper) x-ranges for each parameter.
        Leave None for default, e.g., ranges=[(1.0,2.0), None, (0, 1000)].
    axes: List of matplotlib.axes
        If not None, plot histograms in the currently existing axes.
    lw: Float
        Linewidth of the histogram contour.
    fs: Float
        Font size for texts.

    Deprecated Parameters
    ---------------------
    percentile: Float
        Deprecated. Use quantile instead.

    Returns
    -------
    axes: 1D list of matplotlib.axes.Axes
        List of axes containing the marginal posterior distributions.
    """
    if theme is None:
        theme = 'blue'
    theme = themes[theme]

    if np.ndim(posterior) == 1:
        posterior = np.expand_dims(posterior, axis=1)
    nsamples, npars = np.shape(posterior)

    if pdf is None:
        pdf  = [None]*npars
        xpdf = [None]*npars
    if not isinstance(pdf, list):  # Put single arrays into list
        pdf  = [pdf]
        xpdf = [xpdf]
    # Histogram keywords depending whether one wants the HPD or not:
    hkw = {'edgecolor':'navy', 'color':'b'}
    # Bestfit keywords:
    bkw = {'zorder':2, 'color':'orange'}
    if quantile is not None:
        hkw = {'histtype':'step', 'lw':lw}
        bkw = {'zorder':-1}
    hkw.update(histkeys)

    if yscale != 'absolute':
        hkw.update({'density':True})

    # Set default parameter names:
    if pnames is None:
        pnames = mu.default_parnames(npars)

    # Xranges:
    if ranges is None:
        ranges = np.repeat(None, npars)

    # Set number of rows:
    nrows, ncolumns, npanels = 4, 3, 12
    npages = int(1 + (npars-1)/npanels)

    if axes is None:
        figs = []
        axes = []
        for j in range(npages):
            fig = plt.figure(fignum+j, figsize=(8.5, 11.0))
            figs.append(fig)
            fig.clf()
            fig.subplots_adjust(left=0.1, right=0.97, bottom=0.08, top=0.98,
                hspace=0.5, wspace=0.1)
            for ipar in range(np.amin([npanels, npars-npanels*j])):
                ax = fig.add_subplot(nrows, ncolumns, ipar+1)
                axes.append(ax)
                if ipar%ncolumns == 0:
                    ax.set_ylabel(r"$N$ samples", fontsize=fs)
                else:
                    ax.set_yticklabels([])
    else:
        npages = 1  # Assume there's only one page
        figs = [axes[0].get_figure()]
        for ax in axes:
            ax.set_yticklabels([])

    maxylim = 0
    for ipar in range(npars):
        ax = axes[ipar]
        ax.tick_params(labelsize=fs-1)
        plt.setp(ax.xaxis.get_majorticklabels(), rotation=90)
        ax.set_xlabel(pnames[ipar], size=fs)
        vals, bins, h = ax.hist(posterior[0::thinning,ipar], bins=25,
            range=ranges[ipar], zorder=0, ec=theme['edgecolor'], **hkw)
        # Plot HPD region:
        if quantile is not None:
            PDF, Xpdf, HPDmin = ms.cred_region(posterior[:,ipar], quantile,
                                               pdf[ipar], xpdf[ipar])
            vals = np.r_[0, vals, 0]
            bins = np.r_[bins[0] - (bins[1]-bins[0]), bins]
            # Interpolate xpdf into the histogram:
            f = si.interp1d(bins+0.5*(bins[1]-bins[0]), vals, kind='nearest')
            # Plot the HPD region as shaded areas:
            if ranges[ipar] is not None:
                xran = np.argwhere((Xpdf>ranges[ipar][0])
                                 & (Xpdf<ranges[ipar][1]))
                Xpdf = Xpdf[np.amin(xran):np.amax(xran)]
                PDF  = PDF [np.amin(xran):np.amax(xran)]
            ax.fill_between(Xpdf, 0, f(Xpdf), where=PDF>=HPDmin,
                facecolor=theme['facecolor'], edgecolor='none',
                interpolate=False, zorder=-2, alpha=0.4)
        if bestp is not None:
            ax.axvline(bestp[ipar], dashes=(7,4), lw=1.25,
                color=theme['color'], **bkw)
        maxylim = np.amax((maxylim, ax.get_ylim()[1]))
        if ranges[ipar] is not None:
            ax.set_xlim(
                np.clip(ax.get_xlim(), ranges[ipar][0], ranges[ipar][1]))

    for ax in axes:
        if yscale == 'absolute':
            ax.set_ylim(0, maxylim)

    if savefile is not None:
        for page, fig in enumerate(figs):
            if npages > 1:
                sf = os.path.splitext(savefile)
                fig.savefig("{:s}_page{:02d}{:s}".format(sf[0], page, sf[1]),
                            bbox_inches='tight')
            else:
                fig.savefig(savefile, bbox_inches='tight')

    return axes

