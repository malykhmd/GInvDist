import sympy
#from sympy import S, Basic, gcd, igcd_lehmer, Integer
from sympy import S, Integer, gcd, cancel, expand
from coefs import Coef
if __name__ == '__main__':
  from monom import *
else:
  try:
    from .monom import *
  except:
    from monom import *

class PolyCoef(list):
  def __init__(self, *args):
    if args:
      assert len(args) == 1
      if isinstance(args[0], Monom):
        self.append([args[0], Integer(1)])
      elif isinstance(args[0], PolyCoef):
        for m, k in args[0]:
          self.append([m, k])
      else:
        assert args[0]
        if isinstance(args[0], int):
          self.append([Monom(), Coef(args[0])])
        else:
          self.append([Monom(), args[0]])
    assert self.assertValid()

  def __str__(self):
    if not self:
      return '0'
    else:
      r = []
      for m, c in self:
        if m.degree() == 0 and m.position() < 0:
          r.append(str(c))
        elif c == S.One:
          r.append(str(m))
        else:
          r.append(f"{m}*{c}")
      return " + ".join(r)

  def lm(self):
    assert self
    return self[0][0]

  def lc(self):
    assert self
    return self[0][1]

  def find(self, m): # возвращает себя, если моном м делится на старший терм
    assert self
    return self if m.divisible(self[0][0]) else None

  def cmp(self, other):
    return self.lm().cmp(other.lm())

  def __pos__(self):
    return self

  def __neg__(self):
    p = self.__class__()
    for m, k in self:
      p.append([m, -k])
    return p

  def __radd__(self, other):
    return self+other

  def __add__(self, other):
    if other == S.Zero:
      return self
    if isinstance(other, int): 
      other = self.__class__(Integer(other))
    else:
      other = self.__class__(other)
    p = self.__class__()
    i, j, iend, jend = 0, 0, len(self), len(other)
    while i < iend and j < jend:
      c = self[i][0].cmp(other[j][0])
      if c > 0:
        p.append(self[i])
        i += 1
      elif c < 0:
        p.append(other[j])
        j += 1
      else:
        k = self[i][1] + other[j][1]
        if k:
          p.append([self[i][0], k])
        i += 1
        j += 1
    while i < iend:
      p.append(self[i])
      i += 1
    while j < jend:
      p.append(other[j])
      j += 1
    assert p.assertValid()
    return p

  def __rsub__(self, other):
    return self-other

  def __sub__(self, other):
    if other == S.Zero:
      return self
    if isinstance(other, int): 
      other = self.__class__(Integer(other))
    else:
      other = self.__class__(other)
    p = self.__class__()
    i, j, iend, jend = 0, 0, len(self), len(other)
    while i < iend and j < jend:
      c = self[i][0].cmp(other[j][0])
      if c > 0:
        p.append(self[i])
        i += 1
      elif c < 0:
        p.append([other[j][0], -other[j][1]])
        j += 1
      else:
        k = self[i][1] - other[j][1]
        if k:
          p.append([self[i][0], k])
        i += 1
        j += 1
    while i < iend:
      p.append(self[i])
      i += 1
    while j < jend:
      p.append([other[j][0], -other[j][1]])
      j += 1
    assert p.assertValid()
    return p

  def __rmul__(self, other):
    return self*other

  def __mul__(self, other):
    if other == S.One:
      return self
    if isinstance(other, int): 
      other = self.__class__(Integer(other))
    else:
      other = self.__class__(other)
    t = {}
    for m1, k1 in self:
      for m2, k2 in other:
        m, k = m1*m2, k1*k2
        if m in t:
          t[m] += k
        else:
          t[m] = k
    p = self.__class__()
    for m, k in sorted(t.items(), reverse=True):
      if k:
        p.append([m, k])
    assert p.assertValid()
    return p

  def __pow__(self, other):
    assert other >= 0
    if other == 0:
      p = self.__class__(Integer(1))
    else:
      p = self
      for i in range(1, other):
        p = p*self
    return p

  def prolong(self, var):  # вернуть полином, домноженный на переменную
    p = self.__class__()
    for m, k in self:
      p.append([m.prolong(var), k])
    assert p.assertValid()
    return p

  def diff(self, dv):
    p = self.__class__()
    for m, k in self:
      k = sympy.diff(k, dv)
      if k:
        p.append([m, k])
    assert p.assertValid()
    return p

  def mult(self, m1, k1):
    assert k1 != S.Zero
    p = self.__class__()
    for m2, k2 in self:
      p.append([m1*m2, k1*k2])
    assert p.assertValid()
    return p

  def simplify(self, output=False):
    assert self
    g = self[0][1]
    for monom in self[1:]:
      if g == S.One:
        return
      g = gcd(g, monom[1])
    if output: print('simplify polinom', str(self), ' by ', g)
    for monom in self:
      monom[1] = monom[1]/g # Change?
      if not monom[1].is_integer:
        monom[1] = cancel(monom[1]) # упрощает символьные выражения (сокращает дробь)
    assert self.assertValid()

  def reduction(self, i, other): # делим i-тый моном на другой полином
    assert id(self) != id(other) # проверям что полиномы не равны
    assert 0 <= i < len(self) # проверям что есть i-тый моном
    assert self[i][0].divisible(other.lm()) #проверяем, что i-тый моном делится на старший моном другого полинома
    k = gcd(self[i][1], other.lc()) # считаем коэф. равный НОДу от коеф. i-того монома и старшего коэф. другого полинома
    monom, k1, k2 = self[i][0]/other.lm(), other.lc()/k, -self[i][1]/k  
    # cчитаем итоговый моном, а также старший коэф. другого полинома и отрицательный i-тый коэф., делённые на НОД
    if not k1.is_integer: # если k1 дробь - сокращаем
        k1 = cancel(k1)
    if not k2.is_integer: # если k2 дробь - сокращаем
        k2 = cancel(k2)
    for m in self[:i]: # домножаем все коэф. до i-того на k1
      m[1] *= k1
    del self[i] # удаляем текущий i-тый моном
    j = 1
    while i < len(self) and j < len(other): # пока не кончилисиь мономы в этом или другом полиноме
      c = self[i][0].cmp(other[j][0]*monom) # сравниваем порядок (степени) i-того монома и домноженного j-того монома другого полинома
      if c > 0: # если i-тый моном больше
        self[i][1] *= k1 # домножаем коэфициент i-того монома
        i += 1 # переходим на следующий моном
      elif c < 0: #
        self.insert(i, [other[j][0]*monom, other[j][1]*k2]) # если меньше, вставляем домноженый j-тый моном другого полинома
        i += 1 
        j += 1 # переходим на следующий моном в другом полиноме
      else: # если порядок одинаков, т.е. мономы одинаковы
        k = self[i][1]*k1 + other[j][1]*k2 # домножаем и суммируем коэф.
        if not k.is_integer: # если есть скобки
          k = expand(k) # раскрываем 
        if k == S.Zero: # если коэф. нулевой
          del self[i] # удаляем моном
        else: # иначе 
          self[i][1] = k # задаём расчитанный коэф.
          i += 1 # переходим к правому моному 
        j += 1 # также в переходим и в другом полиноме

    while i < len(self): # если остались не обработанные мономы
      self[i][1] *= k1 # домножаем их
      i += 1
    while j < len(other): # если остались не обработанные мономы в другом полиноме
      self.append([other[j][0]*monom, other[j][1]*k2]) # домножаем и добавляем их
      j += 1
    assert self.assertValid() # проверяем правильность вычислений

  def redRat(self, i, other):
    assert id(self) != id(other)
    assert 0 <= i < len(self)
    assert self[i][0].divisible(other[0][0])
    m2, k2 = self[i][0]/other[0][0], -self[i][1]/other[0][1]
    if not k2.is_integer:
        k2 = cancel(k2)
    del self[i]
    j = 1
    while i < len(self) and j < len(other):
      c = self[i][0].cmp(other[j][0]*m2)
      if c > 0:
        i += 1
      elif c < 0:
        self.insert(i, [other[j][0]*m2, other[j][1]*k2])
        i += 1
        j += 1
      else:
        k = self[i][1] + other[j][1]*k2
        if not k.is_integer:
          k = expand(k)
        if k == S.Zero:
          del self[i]
        else:
          self[i][1] = k
          i += 1
        j += 1
    while j < len(other):
      self.append([other[j][0]*m2, other[j][1]*k2])
      j += 1
    assert self.assertValid()

  def reduce(self, other): # редукция (деление) одного полинома на другой
    assert isinstance(other, Poly) # проверяем, что делим на полином
    assert self and other # проверяем что оба полинома верно корректны
    i, iend, lm = 0, len(self), other[0][0] 
    while i < iend and self[i][0].cmp(lm) > 0: # определяем первый моном меньший(по порядку) чем старший моном делителя
      i += 1 
    if i < iend and self[i][0].cmp(lm) == 0: # редуцируем i-тый моном, если он делится на старший член другого полинома
      self.reduction(i, other)
      assert self.assertValid() # проверяем что результат корректен

  def NFhead(self, pset): # если голова делится на полином из pset, редуцирует голову на подходящий полином, в цикле
    while self:
      p = pset.find(self[0][0])
      if p:
        self.reduction(0, p)
      else:
        break

  def NFtail(self, pset): # если хвост делится на полином из pset, редуцирует хвост на подходящий полином, в цикле
    assert self
    i = 1
    while i < len(self):
      p = pset.find(self[i][0])
      if p:
        self.reduction(i, p)
      else:
        i += 1

  def NF(self, pset):
    i = 0
    while i < len(self):
      p = pset.find(self[i][0])
      if p:
        self.redRat(i, p)
      else:
        i += 1
    return self

  @staticmethod
  def S(self, other):
    assert self and other
    assert self.lm().position() == other.lm().position()
    m, k = self[0][0].lcm(other[0][0]), gcd(self[0][1], other[0][1])
    #assert self[0][1]%k == S.Zero and other[0][1]%k == S.Zero
    return self.mult(m/self[0][0], other[0][1]/k) - other.mult(m/other[0][0], self[0][1]/k)

  def assertValid(self):
    for i in range(len(self)-1):
      if self[i][1] == S.Zero or self[i][0].cmp(self[i+1][0]) <= 0:
        print(self[i], self[i+1])
        return False
    return not self or self[-1][1]

class PolyDiff(Poly):
  __var = None
  __fun = None

  @staticmethod
  def df(*args):
    assert args[0] in PolyDiff.__fun
    m, i = [0 for i in Monom._Monom__var], 1
    while i < len(args):
        assert args[i] in PolyDiff.__var
        v = PolyDiff.__var.index(args[i])
        if i+1 == len(args) or not isinstance(args[i+1], int):
            m[v] += 1
            i += 1
        else:
            m[v] += args[i+1]
            i += 2
    return PolyDiff(Monom(m, pos=PolyDiff.__fun.index(args[0])))

  @staticmethod
  def diff2poly(a): # перевод из симольного выражения SymPy в дифференциальный полином
    if not a.args:
      return PolyDiff(a) # если а не содержит аргументов, переводим в дифф. полином
    else:
      if a.func == sympy.Add: # для функции сложения
        r = PolyDiff.diff2poly(a.args[0]) # переводим первый аргумент в дифф. полином
        for s in a.args[1:]:
          r = r + PolyDiff.diff2poly(s) # и прибавляем все последующие, также переводя их в дифф. полином
      elif a.func == sympy.Mul: # для функции умножения
        r = PolyDiff.diff2poly(a.args[0]) # переводим первый аргумент в дифф. полином
        for s in a.args[1:]:
          r = r * PolyDiff.diff2poly(s) # и домножаем на все последующие, также переводя их в дифф. полином
      elif a.func == sympy.Pow: # для функции степени
        r = PolyDiff.diff2poly(a.args[0])**a.args[1] # (первый аргумент в дифф. полином)^(значение второго аргумента)
      elif repr(a.func) in Monom._Monom__fun: # обработка специальных функций, заданных для мономов)
        r = PolyDiff(Monom(pos=Monom._Monom__fun.index(repr(a.func))))
      else:
        assert a.func == sympy.Derivative # последний вариант - производная
        assert repr(a.args[0].func) in Monom._Monom__fun # проверяем, что под производной стоит обрабатываемая функция
        m = [0 for _ in Monom._Monom__var] # задаём список нулей по количеству объявленных переменных
        for v, d in a.args[1:]:
          m[Monom._Monom__var.index(repr(v))] = d # заполняем список соответствующими значениями для переменных
        r = PolyDiff(Monom(m, pos=Monom._Monom__fun.index(repr(a.args[0].func)))) # переводим получаемый моном в дифф. полином
      return r       

  @staticmethod
  def init():
    assert len(Monom._Monom__fun) >= 0
    PolyDiff.__var = sympy.symbols(Monom._Monom__var, real=True)
    PolyDiff.__fun = sympy.symbols(Monom._Monom__fun, cls=sympy.Function, real=True)
    return PolyDiff.__var, PolyDiff.__fun

  def __init__(self, *args):
    super().__init__(*args)

  def __str__(self):
    if not self:
      return '0'
    else:
      r = []
      for m, c in self:
        if not m:
          r.append(f"{c}")
        else:
          if c == S.One:
            r.append(f"{m.df()}")
          else:
            k = f"{c}"
            if k.find("+") < 0 and k.find("-") < 0:
              r.append(f"{m.df()}*{c}")
            else:
              r.append(f"{m.df()}*({c})")
      return " + ".join(r)

  def prolong(self, var):
    p1, p2 = PolyDiff(), PolyDiff()
    for m, k in self:
      if m:
        p1.append([m.prolong(var), k])
    for m, k in self:
      k = sympy.diff(k, PolyDiff.__var[var])
      if k:
        p2.append([m, k])
    return p1 + p2

  def NFhead(self, pset):
    while self:
      p = pset.find(self[0][0])
      if not p:
        break
      else:
        m2 = self[0][0]/p[0][0]
        for df in m2.expand():
          p = p.prolong(df)
        self.reduction(0, p)

  def NFtail(self, pset):
    assert self
    i = 1
    while i < len(self):
      p = pset.find(self[i][0])
      if not p:
        i += 1
      else:
        m2 = self[i][0]/p[0][0]
        for df in m2.expand():
          p = p.prolong(df)
        self.reduction(i, p)

  def NF(self, pset):
    i = 0
    while i < len(self):
      p = pset.find(self[i][0])
      if not p:
        i += 1
      else:
        m2 = self[i][0]/p[0][0]
        for df in m2.expand():
          p = p.prolong(df)
        self.redRat(i, p)
    return self

if __name__ == '__main__':
  sympy.init_printing()

  var = ['a', 'b', 'c', 'd', 'e', 'f']
  fun = ['u', 'v', 'w']
  Monom.init(var, fun)

  Monom.cmp = Monom.POTlex

  for i, g in enumerate(var):
    globals()[g] = Poly(Monom(i))
  for i, g in enumerate(fun):
    globals()[g] = Poly(Monom(pos=i))

  print(repr(Poly()))
  print(Poly())

  alpha, beta, tau = sympy.symbols('alpha, beta, tau', real=True)
  print(alpha)

  print(repr(Poly(21546)))
  print(Poly(21546))

  print(repr(Poly(tau**4*21546)))
  print(Poly(tau**4*21546))


  h = u*c**5 + f*a*b*tau*1236537 + d*alpha
  print(repr(h))

  g = Poly(h)
  print(g)

  print(g != Poly(h))

  print(h.lm(), h.lc())

  print(g.prolong(4))

  print(repr(g.diff(tau)))
  print(g.diff(tau))

  g.reduction(0, h)
  print(g)

  h *= (4*b**4 + f*tau + beta)**11
  print(h)

  h.NFtail(4*b**4 + f*tau + beta)
  print(h)

  h.NFhead(4*b**4 + f*tau + beta)
  print(h)

  h = (4*b**4 + a*b*c**2 + f*tau + beta)
  h.reduce(a*b*c**2 + f*tau)
  print(h)

  h = (4*b**4 + a*b*c**2 + f*tau + beta)
  h.reduce(a*b*c + f*tau)
  print(h)

  df = PolyDiff.df
  var, fun = PolyDiff.init()
  a, b, c, d, e, f = var
  u, v, w = fun

  print(f"{w}")
  print(f"{4*df(u, b, a, b, 3)!r}")

  h = (4*a*b + c*(c + b**2 + a) + a*c)*(3*df(u, b, e, f, 3) + 4*df(v, b, a) + d)
  print(f"{h!r}")
  print(f"{h}")

  g = h.prolong(2).prolong(3)
  print(f"{g!r}")
  print(f"{g}")

  h.NFhead(4*df(v, a, b)*d)
  print(f"{h!r}")
  print(f"{h}")

  h.NFtail(df(u, b))
  print(f"{h!r}")
