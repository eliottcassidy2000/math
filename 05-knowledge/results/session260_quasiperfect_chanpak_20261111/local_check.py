#!/usr/bin/env python3
"""Focused exact controls for Session260's new all-edge theorem."""
from pathlib import Path
import json
from functools import lru_cache
from itertools import combinations
P=Path(__file__).resolve().parent
def graph(a,b):
 n=a+b+1; out=[0]*n
 for i,j in combinations(range(n),2):
  if j==n-1 and i<a:out[j]|=1<<i
  else:out[i]|=1<<j
 return out
def flip(out,i,j):
 z=out.copy();z[i]^=1<<j;z[j]^=1<<i;return z
def count_paths(out):
 n=len(out);full=(1<<n)-1
 @lru_cache(None)
 def dp(mask,last):
  if mask==full:return 1
  bits=out[last]&~mask;total=0
  while bits:
   bit=bits&-bits;bits-=bit;total+=dp(mask|bit,bit.bit_length()-1)
  return total
 return sum(dp(1<<i,i) for i in range(n))
def req(x,msg):
 if not x: raise RuntimeError(msg)
def f(t): return (1<<(t-1)) if t else 0
def g(l,t): return f(t) if not l else (1<<(l-1))*((1<<t)-1)
records=[];count=0
for n in range(2,6):
 b=n+1; T=graph(n,b); h=count_paths(T); changes=[]
 for i in range(len(T)):
  for j in range(i+1,len(T)):
   delta=count_paths(flip(T,i,j))-h; count+=1
   typ=''
   if j==n+b:
    if i<n:
     rank=i+1; expected=2*f(rank-1)-2*g(n-rank,b);typ='A-z'
    else:
     rank=i-n+1;expected=2*f(b-rank)-2*g(rank-1,n);typ='B-z'
    req(delta==expected,('singleton',n,i,j,delta,expected))
   elif (i<n)==(j<n):
    typ='internal';req((delta==0)==(j==i+1),('internal equality',n,i,j,delta));req(delta>=0,('internal monotone',n,i,j,delta))
   else:
    typ='crossing';boundary=(i==n-1 and j==n)
    req((delta<0)==boundary,('crossing sign',n,i,j,delta))
   changes.append({'edge':[i,j],'type':typ,'delta':delta})
 d=3*(1<<(n-1))-2
 req(min(-x['delta'] for x in changes if x['delta']<0)==d,('minimum',n))
 req(all(x['delta']!=-2 for x in changes),('neighbor',n))
 records.append({'a':n,'b':b,'H':h,'target':h-2,'minimum_decrease':d,'minimizing_edges':[x['edge'] for x in changes if x['delta']==-d],'changes':changes})
print(json.dumps({'status':'PASS','scope':'Every single edge of T(a,a+1),2<=a<=5; all-order theorem rests on proofs','edge_checks':count,'families':records},indent=2))
