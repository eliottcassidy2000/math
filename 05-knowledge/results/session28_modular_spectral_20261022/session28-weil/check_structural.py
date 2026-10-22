#!/usr/bin/env python3
"""Exact structural checks. Standard-library arithmetic only."""
from fractions import Fraction as F
from itertools import permutations, combinations
from math import factorial
import json

def require(condition):
 if not condition: raise RuntimeError("Exact certificate failed")

def run():
 out={"hessian":[],"cayley":[],"moments":[]}
 for n in range(2,8):
  edges=list(combinations(range(n),2)); idx={e:i for i,e in enumerate(edges)}; m=len(edges)
  H=[[0]*m for _ in edges]
  # Independently sum all labelled Hamiltonian path polynomials.
  for p in permutations(range(n)):
   terms=[(idx[tuple(sorted((a,b)))],1 if a<b else -1) for a,b in zip(p,p[1:])]
   for (i,s),(j,t) in combinations(terms,2): H[i][j]+=s*t;H[j][i]+=s*t
  B=[[int(v==a)-int(v==b) for a,b in edges] for v in range(n)]
  c=F(factorial(n-2),2**(n-2))
  errors=0
  for i in range(m):
   for j in range(m):
    dot=sum(B[v][i]*B[v][j] for v in range(n))
    predicted=c*(2*int(i==j)-dot)
    errors += F(H[i][j],2**(n-1)) != predicted
  require(errors==0)
  out['hessian'].append({'n':n,'permutations':factorial(n),'entries':m*m,'c':str(c),'negative_eigenvalue':str(c*(2-n)),'negative_multiplicity':n-1 if n>2 else 0,'positive_eigenvalue':str(2*c),'positive_multiplicity':(n-1)*(n-2)//2,'errors':errors})
 for x,y in [(-1,1),(1,-1),(2,0),(0,2)]:
  denom=x*x+(y+1)**2; num=x*x+(y-1)**2
  require(denom-num==4*y)
  out['cayley'].append({'x':x,'y':y,'norm_squared':str(F(num,denom)),'denominator_minus_numerator':denom-num})
 for d in range(1,9):
  vals=[(2*m+1)**d for m in range(4)]
  determinants=[vals[j]*vals[j+2]-vals[j+1]**2 for j in range(2)]
  require(all(x<0 for x in determinants))
  out['moments'].append({'degree':d,'even_moments':vals,'hankel2':determinants})
 out['fixed63']={'a0_a1_a2':[2,4,6],'ordinary_hankel':-4,'even_subsequence':[2,6,10],'even_hankel':-16}
 # Positive controls: genuine positive spectral measure, and Jordan double pole.
 vals=[2*1**m+3*2**m for m in range(7)]
 require(all(vals[j]*vals[j+2]-vals[j+1]**2>=0 for j in range(5)))
 out['positive_spectral_control']={'moments':vals,'hankel2':[vals[j]*vals[j+2]-vals[j+1]**2 for j in range(5)]}
 # Jordan J=[[1,1],[0,1]], u=[1,0],v=[1,1]: u^T J^m v=m+1.
 J=[[1,1],[0,1]];v=[1,1];seq=[]
 for m in range(8):seq.append(v[0]);v=[sum(J[i][j]*v[j] for j in range(2)) for i in range(2)]
 require(seq==list(range(1,9)))
 out['jordan_control']=seq
 return out
if __name__=='__main__':print(json.dumps(run(),indent=2,sort_keys=True))
