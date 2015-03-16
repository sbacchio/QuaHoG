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
twist = (one+I*gamma_5)/sy.sqrt(2)
G0 = (one+gamma_0)/4.0
G4 = (one+gamma_0)/4.0*I*gamma_5*gamma_x
G5 = (one+gamma_0)/4.0*I*gamma_5*gamma_y
G6 = (one+gamma_0)/4.0*I*gamma_5*gamma_z
G3 = (G4 + G5 + G6)
twdoc = ("(1+i\gamma_5)/\sqrt{2}",)*2
gammas = [
    # Projectors without twist
    dict(name = "P0", doc = "\Gamma_0", mat = G0),
    dict(name = "P4", doc = "\Gamma_x", mat = G4),
    dict(name = "P5", doc = "\Gamma_y", mat = G5),
    dict(name = "P6", doc = "\Gamma_z", mat = G6),
    dict(name = "P3", doc = "\sum_{k=x,y,z}\Gamma_k", mat = G3),
    # Projectors with twist
    dict(name = "TP0T", doc = "%s\Gamma_0%s" % twdoc, mat = twist*G0*twist),
    dict(name = "TP4T", doc = "%s\Gamma_x%s" % twdoc, mat = twist*G4*twist),
    dict(name = "TP5T", doc = "%s\Gamma_y%s" % twdoc, mat = twist*G5*twist),
    dict(name = "TP6T", doc = "%s\Gamma_z%s" % twdoc, mat = twist*G6*twist),
    dict(name = "TP3T", doc = "%s\sum_{k=x,y,z}\Gamma_k%s" % twdoc, mat = twist*G3*twist),]

code = ""
for elem in gammas:
    name = elem["name"]
    doc = elem["doc"]
    mat = elem["mat"]
    mat.simplify()
    for i,side in enumerate(["left", "right"]):
        code += "/* multiply prop by %s from the %s */\n" % (doc, side)
        code += "static void\n"
        nn = (name, "G") if i == 0 else ("G", name)
        code += "prop_proj_%s_%s(_Complex double out[NS*NC][NS*NC], _Complex double in[NS*NC][NS*NC])\n{\n" % nn
        #
        p = sy.zeros(NS,NS)
        for c0 in range(NC):
            for c1 in range(NC):
                for s0 in range(NS):
                    for s1 in range(NS):
                        p[s0,s1] = sy.Symbol('in[%2d][%2d]' % (s0*NC+c0,s1*NC+c1))
                q = mat*p if i == 0 else p*mat
                for s0 in range(NS):
                    for s1 in range(NS):
                        s = str(q[s0,s1])
                        s = '+' + s if s[0] != '-' else s
                        code += ("  out[%2d][%2d] = %s;\n" % (s0*NC+c0, s1*NC+c1, s))
                code += "\n"
        code += "  return;\n}\n\n"
                    
code = template.replace('__PYTHON_SEGMENT_GOES_HERE__',code)
print(code)
