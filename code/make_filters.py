#! /usr/bin/env python
import os

import numpy as np

import pyratbay.tools as pt
import pyratbay.spectrum as ps


def read_obs_file(obs_file, out_path='.', make_filters=True):
    """
    Read observation file.  Generate top-hat filters if requested.
    """
    # SPITZER IRAC filters:
    irac_channel = {3.6: '1', 4.5: '2', 5.8: '3', 8.0: '4'}

    path, data_file = os.path.split(obs_file)
    data_set = data_file.split("_")[0:2]
    data_set = '_'.join(data_set)

    wl, widths, rprs, rprs_err, filters = [], [], [], [], []
    for line in open(obs_file, "r"):
        if line.strip() == "" or line.strip().startswith("#"):
            continue
        elif line.strip().startswith("@"):
            inst = line.strip()[1:]
            continue
        wave, width, data, uncert = np.asarray(line.split(","), float)

        if inst.startswith("irac"):
            channel = irac_channel[wave]
            ffile = f"{{ROOT}}/inputs/filters/spitzer_irac{channel}_sa.dat"
        elif inst.startswith("mips"):
            ffile = "{ROOT}/inputs/filters/spitzer_mips24.dat"
        else:
            ffile = f"{out_path}/{data_set}_{inst}_{wave:5.3f}um.dat"
            if make_filters:
                margin = 0.1 * width
                dlambda = width / 500.0
                dummy = ps.tophat(
                    wave, width, margin=margin, dlambda=dlambda, ffile=ffile)
        wl.append(wave)
        widths.append(width)
        rprs.append(data)
        rprs_err.append(uncert)
        filters.append(ffile)
    depth, depth_err = pt.radius_to_depth(rprs, rprs_err)

    return data_set, wl, width, depth, depth_err, filters


if __name__ == "__main__":
    data_dir = 'inputs/data/'
    pfiles = [
        f"{data_dir}{pfile}"
        for pfile in os.listdir(data_dir)
        if pfile.endswith(".csv")]
    for pfile in pfiles:
        dummy = read_obs_file(pfile, out_path='inputs/filters')

