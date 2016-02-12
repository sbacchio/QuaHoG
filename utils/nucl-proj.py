from __future__ import print_function
import re
import argparse
import numpy as np
import h5py
import sys
import time

NS = 4
I = complex(0,1)

one = np.eye(NS, dtype=complex)

gx = np.array([[ 0, 0, 0,-I],
               [ 0, 0,-I, 0],
               [ 0, I, 0, 0],
               [ I, 0, 0, 0]])

gy = np.array([[ 0, 0, 0,-1],
               [ 0, 0, 1, 0],
               [ 0, 1, 0, 0],
               [-1, 0, 0, 0]])

gz = np.array([[ 0, 0,-I, 0],
               [ 0, 0, 0, I],
               [ I, 0, 0, 0],
               [ 0,-I, 0, 0]])

gt = np.array([[ 0, 0, 1, 0],
               [ 0, 0, 0, 1],
               [ 1, 0, 0, 0],
               [ 0, 1, 0, 0]])

g5 = gx.dot(gy.dot(gz.dot(gt)))

twist = (one + I*g5)/np.sqrt(2), (one - I*g5)/np.sqrt(2)
gproj = (one - gt)/4, (one + gt)/4

def get_dset_names(fname):
    # Returns the list of dataset names in fname of 1-1 ppm and pmm
    # nucleons.
    names = []
    with h5py.File(fname, "r") as fp:
        fp.visititems(lambda x,t: names.append(x) if type(t) is h5py.Dataset else None)
    # Get the datasets ending, with 1-1 and ppm or pmm in the name
    names = [n for n in names if ("ppm" in n or "pmm" in n) and "1-1" in n]
    # Now get the unique datasets before the momentum squared groups
    names = [n for n in names if "arr" in n]
    names = ["/".join(n.split("/")[:-2]) for n in names]
    return list(set(names))

def get_dsets(fname, grp):
    # Returns all the datasets under grp
    with h5py.File(fname, "r") as fp:
        msqs = list(fp[grp].keys())
        dset = {msq: {} for msq in msqs}
        for msq in msqs:
            arr = np.array(fp[grp][msq]['arr'])
            mom = np.array(fp[grp][msq]['mvec'])
            dset[msq] = {'arr': arr, 'mom': mom}
    return dset

def project(arr, flav):
    # Twist depends on flavor (proton or neutron)
    tw = twist[["ppm", "pmm"].index(flav)]
    shape = arr.shape
    T = shape[0]
    arr = arr.reshape([-1, NS, NS])
    fwd = np.array(list(map(lambda x: x.dot(tw.dot(gproj[0].dot(tw))).trace(), arr)))
    bwd = np.array(list(map(lambda x: x.dot(tw.dot(gproj[1].dot(tw))).trace(), arr)))
    # back to original shape
    fwd = fwd.reshape(shape[:-2])
    bwd = bwd.reshape(shape[:-2])
    # reverse "bwd"
    return {"fwd": fwd, "bwd": bwd}
        
def init_proj_h5file(outname, fname):
    # Initialize outname based on fname. Take root and attibs from
    # fname
    with h5py.File(fname, "r") as fp:
        root = list(fp.keys())[0]
        nsrc = list(fp[root].keys())[0]
        attrs = {k: v for k,v in fp[root][nsrc].attrs.items()}
    with h5py.File(outname, "w") as fp:
        top = root + "/" + nsrc
        grp = fp.create_group(top)
        for k,v in attrs.items():
            grp.attrs[k] = np.array(v,'S')
    return

def write_dset(fname, grp_name, arr, mvec):
    with h5py.File(fname, "a") as fp:
        grp = fp.create_group(grp_name)
        grp.create_dataset('arr', arr.shape, dtype=arr.dtype, data=arr)
        grp.create_dataset('mvec', mvec.shape, dtype=mvec.dtype, data=mvec)
    return        

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("FNAME", type=str, help="input file name")
    parser.add_argument("-o", "--output", metavar="F", type=str, default=None,
                        help="output file name (default: FNAME.proj)")
    args = parser.parse_args()
    fname = args.FNAME
    output = args.output    
    if output is None:
        output = fname + ".proj"
    names = get_dset_names(fname)
    data = {n: get_dsets(fname, n) for n in names}
    def flav(n):
        return "ppm" if "ppm" in n else "pmm"
    proj = {n: {msq: project(v['arr'], flav(n)) for msq,v in data[n].items()} for n in names}
    init_proj_h5file(output, fname)
    for n in names:
        for msq,v in sorted(proj[n].items()):
            fwd = proj[n][msq]['fwd']
            bwd = proj[n][msq]['bwd']
            mvec = data[n][msq]['mom']
            grp_name = n + "/fwd/" + msq
            write_dset(output, grp_name, fwd, mvec)
            grp_name = n + "/bwd/" + msq
            write_dset(output, grp_name, bwd, mvec)
    return 0

if __name__ == "__main__":
    sys.exit(main())
