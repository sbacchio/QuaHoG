import os
import sympy as sy
from gammas import *

fname = os.path.basename(__file__)
fname = fname.replace('.py','.h')

with open('_'+fname+'_','r') as fp:
    template = fp.read()

NC = 3
NS = 4

U = sy.zeros(NC,NC)
for c0 in range(NC):
    for c1 in range(NC):
        U[c0,c1] = sy.Symbol('u[%2d][%2d]' % (c0,c1))
        
code = ""
for i,side in enumerate(["left", "right"]):
    for j,u in enumerate(["U","D"]):
        code += "/* multiply prop by %s from the %s */\n" % ("U" if u == "U" else "U^\dagger", side)
        code += "static void\n"
        nn = (u, "G") if i == 0 else ("G", u)
        code += "prop_mul_su3_%s_%s(_Complex double out[NS*NC][NS*NC], _Complex double in[NS*NC][NS*NC], _Complex double u[NC][NC])\n{\n" % nn
        #
        p = sy.zeros(NC,NC)
        for s0 in range(NS):
            for s1 in range(NS):
                for c0 in range(NC):
                    for c1 in range(NC):
                        p[c0,c1] = sy.Symbol('in[%2d][%2d]' % (s0*NC+c0,s1*NC+c1))
                if u == "U":
                    q = p*U if side == "right" else U*p
                if u == "D":
                    q = p*U.H if side == "right" else U.H*p
                for c0 in range(NC):
                    for c1 in range(NC):
                        s = str(q[c0,c1])
                        s = '+' + s if s[0] != '-' else s
                        code += ("  out[%2d][%2d] = %s;\n" % (s0*NC+c0, s1*NC+c1, s)).replace("conjugate", "conj")
                code += "\n"
        code += "  return;\n}\n\n"
                    
code = template.replace('__PYTHON_SEGMENT_GOES_HERE__',code)
print(code)
