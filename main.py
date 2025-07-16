from memory_profiler import memory_usage

from ginv.monom import Monom
from ginv.poly import Poly
from ginv.gb import GB
from ginv.ginv import *
import json
import os
#import sympy

VERY_QUICK = ['quadfor2', 'sparse5', 'hunecke', 'solotarev', 'chandra4', 'quadgrid', 'lorentz', 'liu', 'hemmecke', 'boon', 'chandra5', 'caprasse', 'issac97', 'hcyclic5', 'redcyc5', 'cyclic5', 'extcyc4', 'chemequ', 'uteshev_bikker', 'chandra6', 'geneig']
QUICK = ['chemequs', 'vermeer', 'camera1s', 'reimer4', 'redeco7', 'tangents', 'cassou', 'butcher', 'eco7', 'cohn2', 'dessin1', 'des18_3', 'hcyclic6', 'noon5', 'katsura6', 'cyclic6', 'butcher8', 'redcyc6', 'cpdm5', 'extcyc5']
MEDIUM = ['noon6', 'reimer5', 'kotsireas', 'assur44']


def init(variables, order = Monom.TOPdeglex):
  #variables = ['x1', 'x2', 'x3', 'x4', 'x5']
  Monom.init(variables)
  Monom.variables = variables.copy()
  Monom.zero = Monom(0 for _ in Monom.variables)
  Monom.cmp = order 
  Poly.cmp = order
  for i in range(len(Monom.variables)):
    p = Poly()
    p.append([Monom(0 if l != i else 1 for l in range(len(Monom.variables))), 1])
    globals()[Monom.variables[i]] = p


def test0():
  init(['a', 'b', 'c'], Monom.POTlex)

  print(b**3 - c**2 + a - 1)
  G = GB()
  G.algorithm1([a**3 - b**2 + c - 1, b**3 - c**2 + a - 1, - a**2 + c**3  + b - 1])
  print(G)
  print(", ".join(str(g.lm()) for g in G))
  print("crit1 =", G.crit1, "crit2 =", G.crit2)
  print("time %.2f" % G.time)


def test1():
  init(['x1', 'x2', 'x3', 'x4', 'x5'], Monom.TOPlex)
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
  G.algorithm11(eqs, output=True)
  print(G)
  print(", ".join(str(g.lm()) for g in G))
  print("crit1 =", G.crit1, "crit2 =", G.crit2)
  print("time %.2f" % G.time)

def test2():
  variables = [ 'u0', 'u1', 'u2', 'u3', 'u4', 'u5', 'u6' ]
  init(variables, Monom.POTdeglex)
  equations = [ 'u0**2+2*u1**2+2*u2**2+2*u3**2+2*u4**2+2*u5**2+2*u6**2-u0',
                '2*u0*u1+2*u1*u2-u1+2*u2*u3+2*u3*u4+2*u4*u5+2*u5*u6',
                '2*u0*u2+u1**2+2*u1*u3+2*u2*u4-u2+2*u3*u5+2*u4*u6', 
                '2*u0*u3+2*u1*u2+2*u1*u4+2*u2*u5+2*u3*u6-u3',
                '2*u0*u4+2*u1*u3+2*u1*u5+u2**2+2*u2*u6-u4',
                '2*u0*u5+2*u1*u4+2*u1*u6+2*u2*u3-u5',
                '2*u1+2*u2+2*u3+2*u4+2*u5+2*u6+u0-1'
    ]
  equations = [ 'u0**2-u0+2*u1**2+2*u2**2+2*u3**2+2*u4**2+2*u5**2+2*u6**2',
                '2*u0*u1+2*u1*u2-u1+2*u2*u3+2*u3*u4+2*u4*u5+2*u5*u6',
                '2*u0*u2+u1**2+2*u1*u3+2*u2*u4-u2+2*u3*u5+2*u4*u6', 
                '2*u0*u3+2*u1*u2+2*u1*u4+2*u2*u5+2*u3*u6-u3',
                '2*u0*u4+2*u1*u3+2*u1*u5+u2**2+2*u2*u6-u4',
                '2*u0*u5+2*u1*u4+2*u1*u6+2*u2*u3-u5',
                'u0+2*u1+2*u2+2*u3+2*u4+2*u5+2*u6-1'
    ]
  eqs = []
  for eq in equations:
    eqs.append(eval(eq))
  G = GB()
  G.algorithm2(eqs, output=True)#
  print(G)
  file = open('out3.txt', 'w')
  file.write(str(G))
  file.close()
  print(", ".join(str(g.lm()) for g in G))
  print("crit1 =", G.crit1, "crit2 =", G.crit2)
  print("time %.2f" % G.time)

def receiving_json(test_name, json_data, out=False):
  try:
    size = json_data["dimension"]
    if out:
      print(f"Тест для {test_name}")    
    variables = json_data["variables"]
    init(variables)
    #Monom.init(variables, Monom.TOPdeglex)
    
    equations = json_data['equations']
    eqs = []
    for eq in equations:
      eqs.append(eval(eq.replace('^', '**')))
    G = GB()
    G.algorithm2(eqs)
    leads = ', '.join(str(g.lm()) for g in G)
    data = {
        'time' : str(G.time),
        'dimension' : size,
        'crit1' : str(G.crit1),
        'crit2' : str(G.crit2),
        'leads' : leads,
        'basis' : str(G),}
    if out:
      print(G)
      print("crit1 =", G.crit1, "crit2 =", G.crit2)
      print("time %.2f" % G.time)
    return data
    #array = [poly_int.to_monom("TOP", "deglex", monom.variable(i, len(variables_list), -1)) for i in range(len(variables_list))]
  except Exception as e:
    print(f"Error processing JSON {test_name}: {e}")
    return -1
  

def test_with_memory(test):
  print(test, 'start...')
  res = []
  def f():
    res.append(receiving_json(test, json.load(open('json/'+test+'.json')), True))
  mem_usage = memory_usage(f)
  res[-1]['avr memory'] = sum(mem_usage) / len(mem_usage)
  res[-1]['max memory'] = max(mem_usage)
  print(res[-1])
  with open('resultsM/' + test + '.json', 'w') as file:
    json.dump(res[-1], file)
  print(test, 'complete!')

def test_json():
  tests = os.listdir('json') 
  print(tests)
  print(tests[0])
  res= []
  for test in tests:
    res.append(receiving_json(test.split('.')[0], json.load(open('json/'+test))))
  print(res)
  print('Sum', sum(res) + res.count(-1))
  print('Max', max(res))
  print('Min', min(res))
  print('Unreaded', res.count(-1))

def test_select(tests):
  print('Selected:', tests)
  for test in tests:
    print(test)
    test_with_memory(test)

#test2()

def check():
  file1 = open('out.txt', 'r').read().strip()
  file2 = open('out1.txt', 'r').read().strip()
  file3 = open('out2.txt', 'r').read().strip()
  file4 = open('out3.txt', 'r').read().strip()
  print(len(file1))
  print(len(file2))
  print(len(file3))
  print(len(file4))
  for i in range(len(file3)):
    if file3[i] != file4[i]:
      print(file3[i-50:i], ' !! ', file3[i:i+1000])
      print('__'*100)
      print(file4[i-50:i], ' !! ', file4[i:i+1000])
      break
  #print('file1 == file2', file1[:30000] == file2[:30000])
  #print('file2 == file3', file2[:1500] == file3[:1500])
  #print('file3 == file4', file3[:1500] == file4[:1500])
  #print('file4 == file1', file4[:2000] == file1[:2000])




def collect_results():
  import pandas as pd
  all_data = pd.DataFrame({
    'name' : [],
    'time' : [],
    'avr memory' : [],
    'reduce' : [],
    'dimension' : []
                           })
  files = os.listdir('resultsM')
  for file_name in files:
    data = json.load(open('resultsM/'+file_name))
    all_data.loc[len(all_data)] = [ file_name, data['time'], data['avr memory'], int(data['crit1']) + int(data['crit2']), data['dimension']]

  print(all_data)
  all_data.to_csv('results1.csv', sep=';')


#test1()
#crit1 = 280 crit2 = 401
#collect_results()


#collect_results()

def reverse_order(eq):
  print(eq)
  

#receiving_json('assur44', json.load(open('json/assur44.json')), True)
#check()


# POT position of a term
# -> TOP term of position 

if __name__ == '__main__':
  test_select(VERY_QUICK)
  test_select(QUICK)
  #test_select(MEDIUM)
  #tests = os.listdir('json')
  #tests.remove(VERY_QUICK)
  #tests.remove(QUICK)
  #tests.remove(MEDIUM)
  #test_select(tests)
  collect_results()
  