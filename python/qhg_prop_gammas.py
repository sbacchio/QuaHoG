import os
import sympy as sy
from gammas import *

fname = os.path.basename(__file__)
fname = fname.replace('.py','.h')

with open('_'+fname+'_','r') as fp:
    template = fp.read()

NC = 3
NS = 4

gamma_5 = gamma_x*gamma_y*gamma_z*gamma_t
gammas = [dict(name = "gx", doc = "\gamma_x", mat = gamma_x),
          dict(name = "gy", doc = "\gamma_y", mat = gamma_y),
          dict(name = "gz", doc = "\gamma_z", mat = gamma_z),
          dict(name = "gt", doc = "\gamma_t", mat = gamma_t),
          dict(name = "g5", doc = "\gamma_5", mat = gamma_5),
          dict(name = "g5gx", doc = "\gamma_5\gamma_x", mat = gamma_5*gamma_x),
          dict(name = "g5gy", doc = "\gamma_5\gamma_y", mat = gamma_5*gamma_y),
          dict(name = "g5gz", doc = "\gamma_5\gamma_z", mat = gamma_5*gamma_z),
          dict(name = "g5gt", doc = "\gamma_5\gamma_t", mat = gamma_5*gamma_t),
          dict(name = "C", doc = "C", mat = gamma_t*gamma_y),
          dict(name = "Cg5", doc = "C\gamma_5", mat = gamma_t*gamma_y*gamma_5),
          dict(name = "Cgx", doc = "C\gamma_x", mat = gamma_t*gamma_y*gamma_x),
          dict(name = "Cgy", doc = "C\gamma_y", mat = gamma_t*gamma_y*gamma_y),
          dict(name = "Cgz", doc = "C\gamma_z", mat = gamma_t*gamma_y*gamma_z),
          dict(name = "Cgt", doc = "C\gamma_t", mat = gamma_t*gamma_y*gamma_t),]

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
