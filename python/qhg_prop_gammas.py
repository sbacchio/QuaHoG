import os
import sympy as sy
from gammas import *

fname = os.path.basename(__file__)
fname = fname.replace('.py','.h')

with open('_'+fname+'_','r') as fp:
    template = fp.read()

NC = 3
NS = 4

one = sy.eye(NS)
gamma_5 = gamma_x*gamma_y*gamma_z*gamma_t
gammas = [
    # Scalar
    dict(name = "1", doc = "1", mat = sy.Matrix.eye(NS)),
    # Vector
    dict(name = "g0", doc = "\gamma_0", mat = gamma_0),
    dict(name = "gx", doc = "\gamma_x", mat = gamma_x),
    dict(name = "gy", doc = "\gamma_y", mat = gamma_y),
    dict(name = "gz", doc = "\gamma_z", mat = gamma_z),
    dict(name = "gt", doc = "\gamma_t", mat = gamma_t),
    # Axial
    dict(name = "g5", doc = "\gamma_5", mat = gamma_5),
    dict(name = "g5g0", doc = "\gamma_5\gamma_0", mat = gamma_5*gamma_0),
    dict(name = "g5gx", doc = "\gamma_5\gamma_x", mat = gamma_5*gamma_x),
    dict(name = "g5gy", doc = "\gamma_5\gamma_y", mat = gamma_5*gamma_y),
    dict(name = "g5gz", doc = "\gamma_5\gamma_z", mat = gamma_5*gamma_z),
    dict(name = "g5gt", doc = "\gamma_5\gamma_t", mat = gamma_5*gamma_t),
    # Di-quark
    dict(name = "C", doc = "C", mat = gamma_t*gamma_y),
    dict(name = "Cg5", doc = "C\gamma_5", mat = gamma_t*gamma_y*gamma_5),
    dict(name = "Cgx", doc = "C\gamma_x", mat = gamma_t*gamma_y*gamma_x),
    dict(name = "Cgy", doc = "C\gamma_y", mat = gamma_t*gamma_y*gamma_y),
    dict(name = "Cgz", doc = "C\gamma_z", mat = gamma_t*gamma_y*gamma_z),
    dict(name = "Cgt", doc = "C\gamma_t", mat = gamma_t*gamma_y*gamma_t),
    # Tensor
    dict(name = "g5sitx", doc = "\gamma_5\sigma_{t,x} = -0.5\gamma_5[\gamma_t, \gamma_x]",
         mat = - 0.5*gamma_5*(gamma_t*gamma_x - gamma_x*gamma_t)),
    dict(name = "g5sity", doc = "\gamma_5\sigma_{t,y} = -0.5\gamma_5[\gamma_t, \gamma_y]",
         mat = - 0.5*gamma_5*(gamma_t*gamma_y - gamma_y*gamma_t)),
    dict(name = "g5sitz", doc = "\gamma_5\sigma_{t,z} = -0.5\gamma_5[\gamma_t, \gamma_z]",
         mat = - 0.5*gamma_5*(gamma_t*gamma_z - gamma_z*gamma_t)),
    dict(name = "g5si0x", doc = "\gamma_5\sigma_{0,x} = -0.5\gamma_5[\gamma_0, \gamma_x]",
         mat = - 0.5*gamma_5*(gamma_0*gamma_x - gamma_x*gamma_0)),
    dict(name = "g5si0y", doc = "\gamma_5\sigma_{0,y} = -0.5\gamma_5[\gamma_0, \gamma_y]",
         mat = - 0.5*gamma_5*(gamma_0*gamma_y - gamma_y*gamma_0)),
    dict(name = "g5si0z", doc = "\gamma_5\sigma_{0,z} = -0.5\gamma_5[\gamma_0, \gamma_z]",
         mat = - 0.5*gamma_5*(gamma_0*gamma_z - gamma_z*gamma_0)),
    dict(name = "g5sixy", doc = "\gamma_5\sigma_{x,y} = -0.5\gamma_5[\gamma_x, \gamma_y]",
         mat = - 0.5*gamma_5*(gamma_x*gamma_y - gamma_y*gamma_x)),
    dict(name = "g5sixz", doc = "\gamma_5\sigma_{x,z} = -0.5\gamma_5[\gamma_x, \gamma_z]",
         mat = - 0.5*gamma_5*(gamma_x*gamma_z - gamma_z*gamma_x)),
    dict(name = "g5siyz", doc = "\gamma_5\sigma_{y,z} = -0.5\gamma_5[\gamma_y, \gamma_z]",
         mat = - 0.5*gamma_5*(gamma_y*gamma_z - gamma_z*gamma_y)),
    #
    dict(name = "g5sixt", doc = "\gamma_5\sigma_{x,t} = -0.5\gamma_5[\gamma_x, \gamma_t]",
         mat = - 0.5*gamma_5*(gamma_x*gamma_t - gamma_t*gamma_x)),
    dict(name = "g5siyt", doc = "\gamma_5\sigma_{y,t} = -0.5\gamma_5[\gamma_y, \gamma_t]",
         mat = - 0.5*gamma_5*(gamma_y*gamma_t - gamma_t*gamma_y)),
    dict(name = "g5sizt", doc = "\gamma_5\sigma_{z,t} = -0.5\gamma_5[\gamma_z, \gamma_t]",
         mat = - 0.5*gamma_5*(gamma_z*gamma_t - gamma_t*gamma_z)),
    dict(name = "g5six0", doc = "\gamma_5\sigma_{x,0} = -0.5\gamma_5[\gamma_x, \gamma_0]",
         mat = - 0.5*gamma_5*(gamma_x*gamma_0 - gamma_0*gamma_x)),
    dict(name = "g5siy0", doc = "\gamma_5\sigma_{y,0} = -0.5\gamma_5[\gamma_y, \gamma_0]",
         mat = - 0.5*gamma_5*(gamma_y*gamma_0 - gamma_0*gamma_y)),
    dict(name = "g5siz0", doc = "\gamma_5\sigma_{z,0} = -0.5\gamma_5[\gamma_z, \gamma_0]",
         mat = - 0.5*gamma_5*(gamma_z*gamma_0 - gamma_0*gamma_z)),
    dict(name = "g5siyx", doc = "\gamma_5\sigma_{y,x} = -0.5\gamma_5[\gamma_y, \gamma_x]",
         mat = - 0.5*gamma_5*(gamma_y*gamma_x - gamma_x*gamma_y)),
    dict(name = "g5sizx", doc = "\gamma_5\sigma_{z,x} = -0.5\gamma_5[\gamma_z, \gamma_x]",
         mat = - 0.5*gamma_5*(gamma_z*gamma_x - gamma_x*gamma_z)),
    dict(name = "g5sizy", doc = "\gamma_5\sigma_{z,y} = -0.5\gamma_5[\gamma_z, \gamma_y]",
         mat = - 0.5*gamma_5*(gamma_z*gamma_y - gamma_y*gamma_z)),    
    # Projectors
    dict(name = "1pg0", doc = "(1+\gamma_0)", mat = (one+gamma_0)),
    dict(name = "1pgx", doc = "(1+\gamma_x)", mat = (one+gamma_x)),
    dict(name = "1pgy", doc = "(1+\gamma_y)", mat = (one+gamma_y)),
    dict(name = "1pgz", doc = "(1+\gamma_z)", mat = (one+gamma_z)),
    dict(name = "1pgt", doc = "(1+\gamma_t)", mat = (one+gamma_t)),
    #
    dict(name = "1mg0", doc = "(1-\gamma_0)", mat = (one-gamma_0)),
    dict(name = "1mgx", doc = "(1-\gamma_x)", mat = (one-gamma_x)),
    dict(name = "1mgy", doc = "(1-\gamma_y)", mat = (one-gamma_y)),
    dict(name = "1mgz", doc = "(1-\gamma_z)", mat = (one-gamma_z)),
    dict(name = "1mgt", doc = "(1-\gamma_t)", mat = (one-gamma_t)),
]

code = ""
for elem in gammas:
    name = elem["name"]
    doc = elem["doc"]
    mat = elem["mat"]
    for i,side in enumerate(["left", "right"]):
        code += "/* multiply prop by %s from the %s */\n" % (doc if side == "left" else "\bar{%s}" % doc, side)
        code += "static void\n"
        nn = (name, "G") if i == 0 else ("G", name)
        code += "prop_%s_%s(_Complex double out[NS*NC][NS*NC], _Complex double in[NS*NC][NS*NC])\n{\n" % nn
        #
        p = sy.zeros(NS,NS)
        for c0 in range(NC):
            for c1 in range(NC):
                for s0 in range(NS):
                    for s1 in range(NS):
                        p[s0,s1] = sy.Symbol('in[%2d][%2d]' % (s0*NC+c0,s1*NC+c1))
                q = mat*p if i == 0 else p*gamma_0*mat.H*gamma_0
                for s0 in range(NS):
                    for s1 in range(NS):
                        s = str(q[s0,s1])
                        s = '+' + s if s[0] != '-' else s
                        code += ("  out[%2d][%2d] = %s;\n" % (s0*NC+c0, s1*NC+c1, s))
                code += "\n"
        code += "  return;\n}\n\n"
                    
code = template.replace('__PYTHON_SEGMENT_GOES_HERE__',code)
print(code)
