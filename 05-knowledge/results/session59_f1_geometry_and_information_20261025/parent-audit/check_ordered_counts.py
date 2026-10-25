from itertools import product,combinations
from math import comb

def points(k,q):
 out=[]
 for j in range(k):
  out += [(0,)*j+(1,)+a for a in product(range(q),repeat=k-j-1)]
 return out

def formula(k,q):
 t=q-1
 num=(1+3*t*t+t**3)**k-3*(1+t*t)**k-3*t*(1+t+t*t)**k+t*(2*t+3)*(1+t)**k+2-2*t*t
 if num%t**3: raise ValueError('not integral')
 return num//t**3
for q in (2,3,5):
 for k in range(1,5 if q<5 else 4):
  pp=points(k,q)
  brute=6*sum(all(sum(v[i]!=0 for v in triple)!=1 for i in range(k)) for triple in combinations(pp,3))
  f=formula(k,q)
  if brute!=f: raise ValueError((q,k,brute,f))
  print(q,k,len(pp),brute,'PASS')
