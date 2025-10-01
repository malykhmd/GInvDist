# GInvDist
```
>>> from ginv.monom import *
>>> from ginv.poly import *
>>> from ginv.gb import *
>>> variables = ['x' ,'y']
>>> Monom.init(variables)
>>> Monom.variables = variables.copy()
>>> Monom.zero = Monom(0 for _ in Monom.variables)
>>> p_x = Poly()
>>> p_x.append([Monom(0), 1])
>>> p_y = Poly()
>>> p_y.append([Monom(1), 1])
>>> globals()[Monom.variables[0]] = p_x
>>> globals()[Monom.variables[1]] = p_y
>>> eq = [x**2+y*x-2,  x**3+x*y-4]
>>> G = GB()
>>> G.algorithm2(eq)
>>> print(G)
y**2*-1 + x*4 + y + -6,
x*y + x*2 + y + -2,
x**2 + x*-2 + y*-1
```
