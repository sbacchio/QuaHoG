from __future__ import print_function
import re
import getopt
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
    fwd = np.array(list(map(lambda x: tw.dot(x.dot(gproj[0].dot(tw))), arr))).trace(1,2)
    bwd = np.array(list(map(lambda x: tw.dot(x.dot(gproj[1].dot(tw))), arr))).trace(1,2)
    # back to original shape
    fwd = fwd.reshape(shape[:-2])
    bwd = bwd.reshape(shape[:-2])
    # reverse "bwd"
    t0 = np.arange(T, dtype=int)
    tm = (T - t0) % T
    bwd = bwd[tm,:]
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

class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg

def main(argv=None):
    # defaults
    output = None
    if argv is None:
        argv = sys.argv
        usage =  " Usage: %s [OPTIONS] FNAME" % argv[0]
    try:
        try:
            opts, args = getopt.getopt(argv[1:], "ho:", ["help", "output="])
        except getopt.GetoptError as msg:
            raise Usage(msg)

        for o, a in opts:
            if o in ["-h", "--help"]:
                print((usage), file=sys.stderr)
                print((" Options:"), file=sys.stderr)
                print(("  -o, --output=F\t: Write to file F (default: FNAME.proj)"), file=sys.stderr)
                print(("  -h, --help\t\t: This help message"), file=sys.stderr)
                print 
                return 2
            elif o in ["-o", "--output"]:
                output = a
            else:
                print((" %s: ignoring unhandled option" % o), file=sys.stderr)

        if len(args) != 1:
            raise Usage(usage)

    except Usage as err:
        print((err.msg), file=sys.stderr)
        print((" for help use --help"), file=sys.stderr)
        return 2

    fname = args[0]
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
