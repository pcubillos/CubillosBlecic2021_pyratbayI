#! /usr/bin/env python

import os
import numpy as np
import pyratbay.tools as pt


if __name__ == "__main__":
    # Gather files:
    dpath = "../inputs/data/"
    fdata = sorted(f for f in os.listdir(dpath) if f.endswith(".csv"))

    # Prep up SPITZER filters:
    irac = ["{ROOT}/inputs/filters/spitzer_irac1_sa.dat",
            "{ROOT}/inputs/filters/spitzer_irac2_sa.dat",
            "{ROOT}/inputs/filters/spitzer_irac3_sa.dat",
            "{ROOT}/inputs/filters/spitzer_irac4_sa.dat"]
    swave = np.array([3.6, 4.5, 5.8, 8.0, 24])

    # Read data:
    for pfile in fdata:
        planet = pfile.split("_")[0]

        data    = "data =   "
        error   = "uncert = "
        filters = "filters =\n"
        k = 0

        for line in open(dpath + pfile, "r"):
            if line.strip() == "" or line.strip().startswith("#"):
                continue
            elif line.strip().startswith("@"):
                inst = line.strip()[1:]
                continue

            wl, width, rprs, runc = line.split(",")
            wl    = float(wl)
            width = float(width)
            rprs  = float(rprs)
            runc  = float(runc)
            if inst.startswith("irac"):
                k = np.argmin(np.abs(swave-wl))
                ffile = irac[k]
            elif inst.startswith("mips"):
                ffile = "{ROOT}/inputs/filters/spitzer_mips24.dat"
            else:
                ffile = f"../inputs/filters/{planet}_{inst}_{wl:5.3f}um.dat"
            if inst.startswith("stis") or inst.startswith("wfc3") \
               or inst.startswith("acs") or inst.startswith("nicmos"):
                w,t = pt.tophat(wl, width, 0.1*width, width/500.0, ffile=ffile)
            depth, depth_err = pt.radius_to_depth(rprs, runc)
            data  += "{:.8f}  ".format(depth)
            error += "{:.8f}  ".format(depth_err)
            if (k+1) % 5 == 0 and k != 0:
                data  += "\n         "
                error += "\n         "
            filters += "         {:s}\n".format(ffile)
            k += 1

        print(planet)
        print(data)
        print(error)
        print(filters)
