from __future__ import print_function
import argparse
import numpy as np
import h5py
import sys
import time
import re

ND = 4
NS = 4

def get_moms(cut):
    # Returns an array of shape [:, 4], of all momenta with msq =
    # \sum_i(p_i**2), i=x,y,z for which msq <= cut. [:,0] are the
    # msq values, while [:,1:] are the px, py and pz values.
    v = np.mgrid[-cut:cut, -cut:cut, -cut:cut]
    r = v[0,:,:,:]**2 + v[1,:,:,:]**2 + v[2,:,:,:]**2
    idx = r <= cut
    v = v[:,idx].reshape([3,-1])
    mv = np.array([(v**2).sum(axis=0),v[0,:],v[1,:],v[2,:]])
    mv = np.array(sorted(mv.T, key=lambda x: tuple(x)), int)
    return mv

def get_first_line(fname):
    # Get first line of fname. Quit if the file does not appear to be
    # ascii
    with open(fname, "rb") as fp:
        line = fp.readline()
        try:
            line = line.decode()
        except Exception:
            print(" %s: does not appear to be a text file" % fname)
            sys.exit(1)                
    return line

def get_data_type(fname):
    line = get_first_line(fname)
    # Infer data type from file name. Check first line in file for
    # consistency
    line = line.split()
    if "thrp" in fname:
        errline = " %s: file name suggests three point function but file content fails sanity check" % fname
        assert len(line) == 7, errline
        assert re.match("[0-9]*", line[0]) != None, errline
        assert re.match("[-+][0-9]* [-+][0-9]* [-+][0-9]*", " ".join(line[1:4])) != None, errline                    
        assert re.match("[-+0-9\.e]]* *[-+0-9\.e]]*", " ".join(line[4:6])) != None, errline                    
        assert re.match("=[a-z:0-9]*=", line[6]) != None, errline                    
        return "thrp"
    if "mesons" in fname:
        errline = " %s: file name suggests meson two point function but file content fails sanity check" % fname        
        assert len(line) == 8, errline
        assert re.match("[0-9]*", line[0]) != None, errline
        assert re.match("[-+][0-9]* [-+][0-9]* [-+][0-9]*", " ".join(line[1:4])) != None, errline                    
        assert re.match("[-+0-9\.e]]* *[-+0-9\.e]]*", " ".join(line[4:6])) != None, errline                    
        assert re.match("=[a-z:0-9]*=", line[6]) != None, errline                    
        assert re.match("[ud][pn]", line[7]) != None, errline                    
        return "meson"
    if "nucleons" in fname:
        errline = " %s: file name suggests nucleon two point function but file content fails sanity check" % fname        
        assert len(line) == 14, errline
        assert re.match("[0-9]*", line[0]) != None, errline
        assert re.match("[-+][0-9]* [-+][0-9]* [-+][0-9]*", " ".join(line[1:4])) != None, errline                    
        assert re.match("[-+0-9\.e]]* *[-+0-9\.e]]*", " ".join(line[ 4: 6])) != None, errline                    
        assert re.match("[-+0-9\.e]]* *[-+0-9\.e]]*", " ".join(line[ 6: 8])) != None, errline                    
        assert re.match("[-+0-9\.e]]* *[-+0-9\.e]]*", " ".join(line[ 8:10])) != None, errline                    
        assert re.match("[-+0-9\.e]]* *[-+0-9\.e]]*", " ".join(line[10:12])) != None, errline                    
        assert re.match("p[pm]m", line[12]) != None, errline
        assert re.match("[12]-[12]", line[13]) != None, errline                    
        return "nucleon"
    if "nproj" in fname:
        errline = " %s: file name suggests projected nucleon two point function but file content fails sanity check" % fname        
        assert len(line) == 12, errline
        assert re.match("[0-9]*", line[0]) != None, errline
        assert re.match("[-+][0-9]* [-+][0-9]* [-+][0-9]*", " ".join(line[1:4])) != None, errline                    
        assert re.match("[-+0-9\.e]]* *[-+0-9\.e]]*", " ".join(line[ 4: 6])) != None, errline                    
        assert re.match("[-+0-9\.e]]* *[-+0-9\.e]]*", " ".join(line[ 6: 8])) != None, errline                    
        assert re.match("[-+0-9\.e]]* *[-+0-9\.e]]*", " ".join(line[ 8:10])) != None, errline                    
        assert re.match("[-+0-9\.e]]* *[-+0-9\.e]]*", " ".join(line[10:12])) != None, errline                    
        return "proj-nucl"

def nproj_toh5(fname, outname, root):
    data = {}
    def parse_line(line):
        x = line.split()
        if x == []:
            return None
        t = int(x[0])
        mom = tuple(map(int, x[1:4]))
        reim = np.array(list(map(float, x[4:]))).reshape([-1,2])
        p0,p1,n0,n1 = [complex(*x) for x in reim]
        return t,mom,p0,p1,n0,n1
    with open(fname, "r") as fp:
        for line in fp.readlines():
            x = parse_line(line)
            if x is not None:
                t,mom,prfw,prbw,nufw,nubw = x
                data[t,mom,"fwd","ppm"] = prfw
                data[t,mom,"bwd","ppm"] = prbw
                data[t,mom,"fwd","pmm"] = nufw
                data[t,mom,"bwd","pmm"] = nubw
            else:
                break
    with h5py.File(outname, "w") as fp:
        top = fp.require_group(root)
        T = len(set([x[0] for x in data.keys()]))
        moms = np.array([x[1] for x in data.keys()])
        max_qsq = max((moms**2).sum(axis=1))
        mvec = get_moms(max_qsq)
        for flv in ["ppm", "pmm"]:
            for d in ["fwd", "bwd"]:
                for msq in sorted(set(mvec[:,0])):
                    grp_name = "%s/1-1/%s/msq%04d/" % (flv, d, msq)
                    grp = top.create_group(grp_name)
                    mv = mvec[mvec[:,0] == msq, 1:]
                    arr = np.empty([T,len(mv)], complex)
                    mom = np.empty([len(mv),3], int)
                    for t in range(T):
                        for i,m in enumerate(mv):
                            arr[t,i] = data[t,tuple(m),d,flv]
                            mom[i,:] = m
                    grp.create_dataset('arr', arr.shape, dtype=arr.dtype, data=arr)
                    grp.create_dataset('mvec', mom.shape, dtype=mom.dtype, data=mom)
    return

def nucleon_toh5(fname, outname, root):
    data = {}
    def parse_line(line):
        x = line.split()
        if x == []:
            return None
        t = int(x[0])
        mom = tuple(map(int, x[1:4]))
        reim = np.array(list(map(float, x[4:12]))).reshape([-1,2])
        reim = np.array([complex(*x) for x in reim]).reshape([ND])
        flv = x[12]
        intrp = x[13]
        return t,mom,reim,flv,intrp
    with open(fname, "r") as fp:
        for line in fp.readlines():
            x = parse_line(line)
            if x is not None:
                t,mom,reim,flv,intrp = x
                data.setdefault((t,mom,flv,intrp), []).append(reim)
            else:
                break
    with h5py.File(outname, "w") as fp:
        top = fp.require_group(root)
        T = len(set([x[0] for x in data.keys()]))
        moms = np.array([x[1] for x in data.keys()])
        max_qsq = max((moms**2).sum(axis=1))
        mvec = get_moms(max_qsq)
        flvs = list(set([x[2] for x in data.keys()]))
        intrps = list(set([x[3] for x in data.keys()]))
        for flv in flvs:
            for intrp in intrps:
                for msq in sorted(set(mvec[:,0])):
                    grp_name = "%s/%s/msq%04d/" % (flv, intrp, msq)
                    grp = top.create_group(grp_name)
                    mv = mvec[mvec[:,0] == msq, 1:]
                    arr = np.empty([T,len(mv),ND,ND], complex)
                    mom = np.empty([len(mv),3], int)
                    for t in range(T):
                        for i,m in enumerate(mv):
                            arr[t,i,:,:] = np.array(data[t,tuple(m),flv,intrp])
                            mom[i,:] = m
                    grp.create_dataset('arr', arr.shape, dtype=arr.dtype, data=arr)
                    grp.create_dataset('mvec', mom.shape, dtype=mom.dtype, data=mom)
    return

def meson_toh5(fname, outname, root):
    data = {}
    def parse_line(line):
        x = line.split()
        if x == []:
            return None
        t = int(x[0])
        mom = tuple(map(int, x[1:4]))
        reim = complex(float(x[4]), float(x[5]))
        gam = x[6]
        flv = x[7]
        return t,mom,reim,gam,flv
    with open(fname, "r") as fp:
        for line in fp.readlines():
            x = parse_line(line)
            if x is not None:
                t,mom,reim,gam,flv = x
                data[t,mom,gam,flv] = reim
            else:
                break
    with h5py.File(outname, "w") as fp:
        top = fp.require_group(root)
        T = len(set([x[0] for x in data.keys()]))
        moms = np.array([x[1] for x in data.keys()])
        max_qsq = max((moms**2).sum(axis=1))
        mvec = get_moms(max_qsq)
        gams = list(set([x[2] for x in data.keys()]))
        flvs = list(set([x[3] for x in data.keys()]))
        for flv in flvs:
            for gam in gams:
                for msq in sorted(set(mvec[:,0])):
                    grp_name = "%s/%s/msq%04d/" % (flv, gam, msq)
                    grp = top.create_group(grp_name)
                    mv = mvec[mvec[:,0] == msq, 1:]
                    arr = np.empty([T,len(mv)], complex)
                    mom = np.empty([len(mv),3], int)
                    for t in range(T):
                        for i,m in enumerate(mv):
                            arr[t,i] = np.array(data[t,tuple(m),gam,flv])
                            mom[i,:] = m
                    grp.create_dataset('arr', arr.shape, dtype=arr.dtype, data=arr)
                    grp.create_dataset('mvec', mom.shape, dtype=mom.dtype, data=mom)
    return

def thrp_toh5(fname, outname, root):
    data = {}
    def parse_line(line):
        x = line.split()
        if x == []:
            return None
        t = int(x[0])
        mom = tuple(map(int, x[1:4]))
        reim = complex(float(x[4]), float(x[5]))
        gam = x[6]
        return t,mom,reim,gam
    with open(fname, "r") as fp:
        for line in fp.readlines():
            x = parse_line(line)
            if x is not None:
                t,mom,reim,gam = x
                data[t,mom,gam] = reim
            else:
                break
    with h5py.File(outname, "w") as fp:
        top = fp.require_group(root)
        T = len(set([x[0] for x in data.keys()]))
        moms = np.array([x[1] for x in data.keys()])
        max_qsq = max((moms**2).sum(axis=1))
        mvec = get_moms(max_qsq)
        gams = list(set([x[2] for x in data.keys()]))
        for gam in gams:
            for msq in sorted(set(mvec[:,0])):
                grp_name = "%s/msq%04d/" % (gam, msq)
                grp = top.create_group(grp_name)
                mv = mvec[mvec[:,0] == msq, 1:]
                arr = np.empty([T,len(mv)], complex)
                mom = np.empty([len(mv),3], int)
                for t in range(T):
                    for i,m in enumerate(mv):
                        arr[t,i] = np.array(data[t,tuple(m),gam])
                        mom[i,:] = m
                grp.create_dataset('arr', arr.shape, dtype=arr.dtype, data=arr)
                grp.create_dataset('mvec', mom.shape, dtype=mom.dtype, data=mom)
    return

def main():
    parser = argparse.ArgumentParser(epilog=" Note: this script infers the file structure from the input file, so use with care ")
    parser.add_argument("FNAME", type=str)
    parser.add_argument("-o", "--output", metavar="F", type=str,
                        help="output file name (default: FNAME.h5)")
    parser.add_argument("-r", "--root", metavar="R", default="/", type=str,
                        help="root of HDF5 file. Allows multiple groups, e.g. via /A/B/C/... (default: %s)" % "/")
    args = parser.parse_args()    
    fname = args.FNAME
    output = args.output
    root = args.root
    if output is None:
        output = fname + ".h5"

    ### What kind of file are we parsing?
    dt = get_data_type(fname)
    print(" %s file detected" % dt)
    if dt == "proj-nucl":
        nproj_toh5(fname, output, root)
    if dt == "nucleon":
        nucleon_toh5(fname, output, root)
    if dt == "meson":
        meson_toh5(fname, output, root)
    if dt == "thrp":
        thrp_toh5(fname, output, root)
    return 0

if __name__ == "__main__":
    sys.exit(main())
