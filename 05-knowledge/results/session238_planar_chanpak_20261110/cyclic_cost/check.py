#!/usr/bin/env python3
"""Exact cyclic prime-orbit floors; independent permutation controls through p=7."""
import itertools,json,sys
from collections import Counter
from pathlib import Path

def require(condition):
 if not condition:raise RuntimeError('verification gate failed')
def matrix(n,S):return [[(j-i)%n in S for j in range(n)] for i in range(n)]
def hp(a):
 n=len(a); d=[[0]*n for _ in range(1<<n)]
 for i in range(n):d[1<<i][i]=1
 for m in range(1,1<<n):
  for i in range(n):
   if d[m][i]:
    for j in range(n):
     if not m>>j&1 and a[i][j]:d[m|1<<j][j]+=d[m][i]
 return sum(d[-1])
def brute(a):return sum(all(a[s[i]][s[i+1]] for i in range(len(a)-1)) for s in itertools.permutations(range(len(a))))
def ordinal(blocks):
 ns=list(map(len,blocks)); n=sum(ns); a=[[False]*n for _ in range(n)]; k=0
 for b in blocks:
  for i,row in enumerate(b):
   for j,v in enumerate(row):a[k+i][k+j]=v
   for j in range(k+len(b),n):a[k+i][j]=True
  k+=len(b)
 return a
rows={}
for p in (3,5,7,11):
 vals=[]
 for signs in itertools.product((0,1),repeat=(p-1)//2):
  S=[p-i if s else i for i,s in enumerate(signs,1)]; a=matrix(p,S); h=hp(a)
  require(h%p==0)
  if p<=7:require(h==brute(a))
  vals.append({'S':S,'H':h})
 rows[str(p)]={'universe':2**((p-1)//2),'histogram':dict(sorted(Counter(v['H'] for v in vals).items())),'minimum':min(v['H'] for v in vals),'rows':vals}
expected7=177 if '--inject-failure' in sys.argv else 175
require(rows['3']['minimum']==3 and rows['5']['minimum']==15 and rows['7']['minimum']==expected7)
# Each minimizing seven-circulant is an affine relabeling of the carousel.
carousel=set(range(1,4)); forms={frozenset(u*i%7 for i in carousel) for u in range(1,7)}
require({frozenset(v['S']) for v in rows['7']['rows'] if v['H']==175}==forms)
q=matrix(7,{1,2,3}); paley=matrix(7,{1,2,4}); single=[[False]]
require(hp(ordinal([single,q,single]))==175)
require(hp(ordinal([q,q]))==175**2)
require(hp(paley)==189)
# A directed triangle between three blocks creates strict product growth.
a=ordinal([matrix(3,{1}),single,single]); a[4][0]=a[4][1]=a[4][2]=True
for i in range(3):a[i][4]=False
require(hp(a)>3)
out={'prime_orbit_census':rows,'controls':{'Q7':175,'P7':189,'two_Q7_blocks':175**2,'cross_triangle_H':hp(a)},'status':'PASS'}
Path(__file__).with_name('results.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({p:{k:v for k,v in row.items() if k!='rows'} for p,row in rows.items()},indent=2));print('PASS')
