#!/usr/bin/env python3
"""Small independent exact positivity control for sign-even Walsh coordinates."""
from itertools import combinations, permutations
from fractions import Fraction
import json

def need(x,s):
 if not x: raise RuntimeError(s)

def setup(n):
 edges=list(combinations(range(n),2)); index={e:i for i,e in enumerate(edges)}
 actions=[]
 for p in permutations(range(n)):
  dest=[];rev=0
  for k,(a,b) in enumerate(edges):
   u,v=p[a],p[b];dest.append(index[tuple(sorted((u,v)))])
   if u>v:rev|=1<<k
  actions.append((dest,rev))
 return edges,actions

def act(mask,action,tournament=False):
 dest,rev=action; value=mask^rev if tournament else mask
 return sum(((value>>k)&1)<<j for k,j in enumerate(dest))

def rank(A):
 A=[[Fraction(x) for x in row] for row in A];r=0
 for j in range(len(A[0])):
  piv=next((i for i in range(r,len(A)) if A[i][j]),None)
  if piv is None:continue
  A[r],A[piv]=A[piv],A[r];a=A[r][j];A[r]=[x/a for x in A[r]]
  for i in range(len(A)):
   if i!=r and A[i][j]:
    a=A[i][j];A[i]=[x-a*y for x,y in zip(A[i],A[r])]
  r+=1
 return r

def main():
 rows=[]
 for n in range(1,5):
  edges,actions=setup(n);m=len(edges);full=(1<<m)-1
  ts=sorted({min(act(t,a,True) for a in actions) for t in range(1<<m)})
  gs=sorted({min(act(g,a) for a in actions) for g in range(1<<m)})
  even=[];orbits=[]
  for g in gs:
   if any(act(g,a)==g and (g&a[1]).bit_count()%2 for a in actions):continue
   orbit={}
   for a in actions:
    h=act(g,a);s=(-1)**((g&a[1]).bit_count()%2)
    need(h not in orbit or orbit[h]==s,'orbit sign not well-defined')
    orbit[h]=s
   even.append(g);orbits.append(orbit)
  M=[[sum(s*((-1)**((h&(~t)&full).bit_count())) for h,s in orbit.items()) for orbit in orbits] for t in ts]
  need(len(ts)==len(even),'type count discrepancy')
  need(rank(M)==len(ts),'Walsh type matrix singular')
  reverse=[ts.index(min(act(t^full,a,True) for a in actions)) for t in ts]
  diag=[(-1)**g.bit_count() for g in even]
  for i in range(len(ts)):
   for j in range(len(even)):need(M[reverse[i]][j]==M[i][j]*diag[j],'converse transform failure')
  need(any(v<0 for v in diag)==(n==4),'wrong first negative converse coordinate')
  rows.append({'n':n,'tournament_representatives':ts,'sign_even_representatives':even,'walsh_rows':M,'converse_permutation':reverse,'converse_fourier_diagonal':diag})
 last=rows[-1];need(7 in last['sign_even_representatives'],'star3 missing')
 star=last['sign_even_representatives'].index(7)
 need(last['converse_fourier_diagonal'][star]==-1,'star converse not negative')
 # Exact inverse and cone inequalities on every vertex and integer mixture.
 for coeff in __import__('itertools').product(range(3),repeat=4):
  m=[sum(coeff[i]*last['walsh_rows'][i][j] for i in range(4)) for j in range(4)]
  mass,a,b,c=m
  inv=[Fraction(3*mass+c+a,8),Fraction(mass-c+b,8),Fraction(3*mass+c-a,8),Fraction(mass-c-b,8)]
  need(inv==list(coeff),'explicit inverse mismatch')
  need(abs(a)<=3*mass+c and abs(b)<=mass-c,'positive mixture outside cone')
 # No diagonal change of sign affects a diagonal linear operator.
 for signs in range(1<<len(last['walsh_rows'])):
  ss=[(-1)**((signs>>i)&1) for i in range(len(last['walsh_rows']))]
  need([s*d*s for s,d in zip(ss,last['converse_fourier_diagonal'])]==last['converse_fourier_diagonal'],'sign-gauge changed diagonal')
 print(json.dumps({'status':'PASS','convention':'lex edge bits1=up; canonical ordinary type minima; Phi orbit signs via edge inversion','cases':rows,'negative_diagonal_first_order':4,'star_edges':3,'sign_gauges_checked':16},indent=2,sort_keys=True))
if __name__=='__main__':main()
