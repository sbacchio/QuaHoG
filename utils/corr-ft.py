from __future__ import print_function
import getopt
import numpy as np
import h5py
import sys
import time

ND = 4

def get_moms(cut, dims):
    # Returns an array of shape [:, 4], of all momenta with msq =
    # \sum_i(p_i**2), i=x,y,z for which msq <= cut. [:,0] are the
    # msq values, while [:,1:] are the px, py and pz values.
    Lx,Ly,Lz = dims
    v = np.mgrid[-Lx/2:Lx/2, -Ly/2:Ly/2, -Lz/2:Lz/2]
    r = v[0,:,:,:]**2 + v[1,:,:,:]**2 + v[2,:,:,:]**2
    idx = r <= cut
    v = v[:,idx].reshape([3,-1])
    mv = np.array([(v**2).sum(axis=0),v[0,:],v[1,:],v[2,:]])
    mv = np.array(sorted(mv.T, key=lambda x: tuple(x)), int)
    return mv

def get_origin(fname):
    with h5py.File(fname, "r") as fp:
        origin = fp.attrs['Origin']
        # infer whether this file is a three-point function file by
        # inspecting whether it has an attribute "Sink-source
        # separation". If so, set time-component of origin to zero,
        # since it is already shifted
        if "Sink-source separation" in list(fp.attrs.keys()):
            origin[0] = 0
    return origin

def get_dset_names(fname):
    # Returns the full list of dataset names in fname
    with h5py.File(fname, "r") as fp:
        names = []
        fp.visititems(lambda x,t: names.append(x) if type(t) is h5py.Dataset else None)
    return names

def get_dset(fname, name):
    # Returns the dataset with name "name" of fname
    with h5py.File(fname, "r") as fp:
        dset = np.array(fp[name])
    return dset

def get_attrs(fname):
    # Returns HDF5 file attributes
    with h5py.File(fname, "r") as fp:
        d = dict(fp.attrs)
    return d

def init_h5file(fname, root):
    with h5py.File(fname, "w") as fp:
        fp.create_group(root)
    return        

def write_dset(fname, grp_name, arr):
    with h5py.File(fname, "a") as fp:
        grp = fp.create_group(grp_name)
        grp.create_dataset('arr', arr['arr'].shape, dtype=arr['arr'].dtype, data=arr['arr'])
        grp.create_dataset('mvec', arr['mom'].shape, dtype=arr['mom'].dtype, data=arr['mom'])
    return        

class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg

def main(argv=None):
    # defaults
    output = None
    root = "/"
    max_msq = 24
    if argv is None:
        argv = sys.argv
        usage =  " Usage: %s [OPTIONS] FNAME" % argv[0]        
    try:
        try:
            opts, args = getopt.getopt(argv[1:], "ho:m:r:", ["help", "output=", "root=", "max-msq="])
        except getopt.GetoptError as msg:
            raise Usage(msg)

        for o, a in opts:
            if o in ["-h", "--help"]:
                print((usage), file=sys.stderr)
                print((" Options:"), file=sys.stderr)
                print(("  -o, --output=F\t: Write to file F (default: FNAME.mom)"), file=sys.stderr)
                print(("  -m, --max-msq=Q\t: Write out to max momentum squared Q (default: %d)" % max_msq), file=sys.stderr)
                print(("  -r, --root=R\t\t: HDF5 file top group (default: %s)" % root), file=sys.stderr)
                print(("  -h, --help\t\t: This help message"), file=sys.stderr)
                print 
                return 2
            elif o in ["-o", "--output"]:
                output = a
            elif o in ["-r", "--root"]:
                root = a
            elif o in ["-m", "--max-msq"]:
                max_msq = a
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
        output = fname + ".mom"
        
    origin = get_origin(fname)
    datasets = get_dset_names(fname)
    mom_corr = {d: {} for d in datasets}
    for d in datasets:
        corr = get_dset(fname, d)
        dims = np.array(corr.shape[:-1])
        corr = corr.reshape([np.prod(dims),2])
        corr = corr[:,0] + complex(0,1)*corr[:,1]
        corr = corr.reshape(dims)
        for i in range(ND):
            sh = origin[i]
            corr = np.roll(corr, axis=i, shift=-sh)
        coft = np.fft.fftn(corr, s=dims[1:ND], axes=(1,2,3))
        trail_dims = coft.shape[ND:]
        if trail_dims == ():
            trail_dims = (1,)
        coft = coft.reshape(tuple(dims[:ND]) + trail_dims)
        mvecs = get_moms(max_msq, dims[1:ND])
        for msq in sorted(set(mvecs[:,0])):            
            mv = mvecs[mvecs[:,0] == msq, 1:]
            mom_corr[d][msq] = {
                "arr": np.zeros((dims[0], len(mv),) + trail_dims, complex),
                "mom": mv}
            for t in range(dims[0]):
                for i,m in enumerate(mv):
                    p = m
                    x = (p[0]+dims[1])%dims[1]
                    y = (p[1]+dims[2])%dims[2]
                    z = (p[2]+dims[3])%dims[3]
                    mom_corr[d][msq]['arr'][t,i,:] = coft[t,x,y,z,:]
            # If the trailing dimensions are (1,), squeeze them out
            if trail_dims == (1,):
                shape = mom_corr[d][msq]['arr'].shape
                mom_corr[d][msq]['arr'] = mom_corr[d][msq]['arr'].reshape(shape[:2])

    a = get_attrs(fname)
    # Origin will be used as the group under root
    pos = a.pop("Origin")
    pos = "sx%02dsy%02dsz%02dst%02d" % (pos[1], pos[2], pos[3], pos[0])
    top = root + "/" + pos
    # Pop the index order, it is irrelevant
    a.pop("IndexOrder")
    # If there is a "projector" attrib. put as next level group
    if "Projector" in a.keys():
        val = a["Projector"].decode()
        top = top + "/" + val
    # If there is a "flavor" attrib. put as next level group
    if "Flavor" in a.keys():
        val = a["Flavor"].decode()
        top = top + "/" + val
    # Initialize the file, writing the top-level group
    init_h5file(output, top)
    # Loop over dataset names writing the FTed arrays
    for d in datasets:
        grp = top + "/" + "/".join(d.split("/")[:-1])
        for msq in mom_corr[d]:
            subgrp = grp + "/msq%04d" % msq
            write_dset(output, subgrp, mom_corr[d][msq])
    return 0

if __name__ == "__main__":
    sys.exit(main())
