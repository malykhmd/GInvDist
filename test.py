from ginv.monom import *
from ginv.poly import *
from ginv.gb import *

import json
import os

def init(variables, order = Monom.TOPdeglex):
  #variables = ['x1', 'x2', 'x3', 'x4', 'x5']
  Monom.init(variables)
  Monom.variables = variables.copy()
  Monom.zero = Monom(0 for _ in Monom.variables)
  Monom.cmp = order
  for i in range(len(Monom.variables)):
    p = Poly()
    p.append([Monom(0 if l != i else 1 for l in range(len(Monom.variables))), 1])
    globals()[Monom.variables[i]] = p


init(['x1', 'x2', 'x3', 'x4', 'x5'])
eqs = [
  x1 + x2 + x3 + x4 + x5,
  x1*x2 + x1*x5 + x2*x3 + x3*x4 + x4*x5,
  x1*x2*x3 + x1*x2*x5 + x1*x4*x5 + x2*x3*x4 + x3*x4*x5,
  x1*x2*x3*x4 + x1*x2*x3*x5 + x1*x2*x4*x5 + x1*x3*x4*x5 + x2*x3*x4*x5,
  x1*x2*x3*x4*x5 - 1
]
#print(eqs)
#_eqs = []
#for eq in eqs:
  #_eqs.append(sympy.va(eq))
#print(_eqs)
G = GB()
G.algorithm2(eqs, output=True)
print(G)
print(", ".join(str(g.lm()) for g in G))
print("crit1 =", G.crit1, "crit2 =", G.crit2)
print("time %.2f" % G.time) 