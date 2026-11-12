#!/usr/bin/env python3
"""Independent permutation control of the new run compiler."""
import importlib.util,itertools,random,json
from pathlib import Path
P=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('family',P/'families/recurrence_bridge.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
r=random.Random(277)
checks=0
for q in range(1,7):
 for repeat in range(12):
  A=[[0]*q for _ in range(q)]
  for i in range(q):
   for j in range(i+1,q): A[i][j]=r.randrange(2);A[j][i]=1-A[i][j]
  sizes=[0]*q
  for _ in range(r.randrange(1,8)): sizes[r.randrange(q)]+=1
  labels=[(i,j) for i,s in enumerate(sizes) for j in range(s)]
  def edge(a,b):return a[1]<b[1] if a[0]==b[0] else A[a[0]][b[0]]
  brute=sum(all(edge(a,b) for a,b in zip(p,p[1:])) for p in itertools.permutations(labels))
  got=m.compiled_H(A,sizes)
  assert got==brute,(q,sizes,A,got,brute)
  checks+=1
# Transitive quotient: unique path for all module sizes, multinomial total q**n.
for q in range(1,7):
 A=[[int(i<j) for j in range(q)] for i in range(q)]
 for n in range(7):
  assert sum(__import__('math').factorial(n)//__import__('math').prod(__import__('math').factorial(s) for s in z)*m.compiled_H(A,z) for z in m.compositions(n,q))==q**n
# Independent regular aggregate via differential recurrence y'= (1+y),
# or explicit coefficients of (d+1-d exp(x))^{-1} using convolution.
import math
for q in (3,5,7,11):
 d=(q-1)//2; C=[1]
 for n in range(1,31):
  C.append(d*sum(math.comb(n,k)*C[n-k] for k in range(1,n+1)))
  assert q*C[n]//d==m.regular_exact(q,n)
print(json.dumps({'independent_permutation_checks':checks,'transitive_aggregate_checks':42,'regular_convolution_checks':120,'status':'PASS'},indent=2))
