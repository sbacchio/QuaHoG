import os
import sympy as sy
import itertools
import textwrap

fname = os.path.basename(__file__)
fname = fname.replace('.py','.h')

with open('_'+fname+'_','r') as fp:
    template = fp.read()

NC = 3
NS = 4
I = sy.I
eps = [dict(idx = (0,1,2), sign = +1),
       dict(idx = (2,0,1), sign = +1),
       dict(idx = (1,2,0), sign = +1),
       dict(idx = (2,1,0), sign = -1),
       dict(idx = (0,2,1), sign = -1),
       dict(idx = (1,0,2), sign = -1)]

code = ""
A = sy.zeros(NS*NC,NS*NC)
B = sy.zeros(NS*NC,NS*NC)
for c0 in range(NC):
    for c1 in range(NC):
        for s0 in range(NS):
            for s1 in range(NS):
                A[s0*NC+c0, s1*NC+c1] = sy.Symbol('A[%2d][%2d]' % (s0*NC+c0, s1*NC+c1))
                B[s0*NC+c0, s1*NC+c1] = sy.Symbol('B[%2d][%2d]' % (s0*NC+c0, s1*NC+c1))                

ctype = "_Complex double"
for s0 in range(NS):
    for s1 in range(s0+1, NS):
        doc = 'C_{\mu\\nu}^{a\'a} = \sum_{\kappa}\epsilon^{abc}\epsilon^{a\'b\'c\'} A_{%s%s}^{bb\'} B_{%s%s}^{cc\'}'
        C = sy.Matrix(NS*NC, NS*NC, lambda i, j: 0+0*I)
        xx = [""]*4
        xx[s0] = '\kappa'
        xx[s1] = '\kappa'
        xx[xx.index("")] = '\mu'
        xx[xx.index("")] = '\\nu'        
        code += "/*\n %s\n */\n" % (doc % tuple(xx))
        code += "static void\n"
        code += ("prop_contract_%d%d(%s C[NC*NS][NC*NS], %s A[NC*NS][NC*NS], %s B[NC*NS][NC*NS])\n{\n" 
                 % (s0, s1, ctype, ctype, ctype))
        for col0 in eps:
            a0, b0, c0 = col0["idx"]
            sign0 = col0["sign"] 
            for col1 in eps:
                a1, b1, c1 = col1["idx"]
                sign1 = col1["sign"] 
                for mu in range(NS):
                    for nu in range(NS):
                        for ku in range(NS):
                            idx = [-1]*4
                            idx[s0] = ku
                            idx[s1] = ku
                            idx[idx.index(-1)] = mu
                            idx[idx.index(-1)] = nu
                            csp0 = a0+mu*NC
                            csp1 = a1+nu*NC
                            C[csp0, csp1] += A[b0+idx[0]*NC, b1+idx[1]*NC] * B[c0+idx[2]*NC, c1+idx[3]*NC] * sign0 * sign1
        lines = [str(C[i, j].expand()) for i,j in itertools.product(range(NS*NC), range(NS*NC))]
        for i,l in enumerate(lines):
            lines[i] = '+ ' + l if l[0] != '-' else '- ' + l[1:]
        lines = ["\n  C[%2d][%2d] = \n\t%s;\n" % (i, j, "\n\t".join(textwrap.wrap(lines[i*NS*NC+j], 88)))
                 for i,j in itertools.product(range(NS*NC), range(NS*NC))]
        code += "".join(lines)
        code += "\n  return;\n}\n\n"

code = template.replace('__PYTHON_SEGMENT_GOES_HERE__',code)
print(code)
        
