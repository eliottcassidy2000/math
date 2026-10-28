#!/usr/bin/env python3
"""Exact hostile controls. Pure Python, no dependencies; run from any cwd."""
from itertools import combinations, permutations
import json
from pathlib import Path

def check(condition, message):
 if not condition:
  raise RuntimeError(message)

def tour(n,m):
 a=[[0]*n for _ in range(n)]
 for k,(i,j) in enumerate(combinations(range(n),2)):
  a[i][j]=(m>>k)&1;a[j][i]=1-a[i][j]
 return a

def path(a,s,t):
 q=[(s,)]; seen={s}
 for p in q:
  if p[-1]==t:return p
  for v,x in enumerate(a[p[-1]]):
   if x and v not in seen:seen.add(v);q.append(p+(v,))
 return None

def h(a):
 n=len(a);d=[[0]*n for _ in range(1<<n)]
 for i in range(n):d[1<<i][i]=1
 for m in range(1<<n):
  for j in range(n):
   if d[m][j]:
    for k in range(n):
     if a[j][k] and not(m>>k&1):d[m|1<<k][k]+=d[m][j]
 return sum(d[-1])

def aut(a):
 return sum(all(a[i][j]==a[p[i]][p[j]] for i in range(len(a)) for j in range(len(a))) for p in permutations(range(len(a))))

def root(n,i,j):return [int(k==i)-int(k==j) for k in range(n)]
def dot(a,b):return sum(x*y for x,y in zip(a,b))
def mm(a,b):return [[sum(x*y for x,y in zip(r,c)) for c in zip(*b)] for r in a]
def independence(n,edges):return sum(2**m.bit_count() for m in range(1<<n) if all(not(m>>i&1 and m>>j&1) for i,j in edges))

rows=[]
for n in range(2,6):
 strong=0;trees=0;hs=set()
 for m in range(1<<(n*(n-1)//2)):
  a=tour(n,m)
  paths=[path(a,s,t) for s in range(n) for t in range(n) if s!=t]
  if any(p is None for p in paths):continue
  strong+=1;hs.add(h(a))
  for p in paths:
   v=root(n,p[0],p[1])
   for u,w in zip(p[1:-1],p[2:]):
    r=root(n,u,w);check(dot(v,r)==-1, "Path merge must have beta -1")
    v=[x+y for x,y in zip(v,r)]
   check(v==root(n,p[0],p[-1]), "Path tree endpoint mismatch");trees+=1
 rows.append(dict(n=n,strong_labelled=strong,H_spectrum=sorted(hs),zero_cost_root_trees=trees))

special=[]
for n in (3,7):
 qr={k*k%n for k in range(1,n)}
 a=[[int(i!=j and (j-i)%n in qr) for j in range(n)] for i in range(n)]
 s=[[a[i][j]-a[j][i] for j in range(n)] for i in range(n)]
 check(all(sum(r)==0 for r in s), "Paley skew matrix not regular")
 check(mm(s,s)==[[1-n*int(i==j) for j in range(n)] for i in range(n)], "Paley skew-square identity failed")
 H=h(a);A=aut(a)
 special.append(dict(name=f'Paley{n}',H=H,Aut=A,tc=H//A,skew_square_on_sumzero=-n))
# Positive-debt A2 has same symmetric edge roots for both T3 classes.
check(h(tour(3,0))==1, "Transitive T3 control failed")
# Orthogonal helper certificate: 2alpha=(alpha+h)+(alpha-h).
a=(1,-1,0,0);b=(0,0,1,-1)
x=tuple(u+v for u,v in zip(a,b));y=tuple(u-v for u,v in zip(a,b))
check(dot(a,b)==dot(x,y)==0, "Perpendicular helper certificate failed")
out=dict(status='FINITE-EXACT, exhaustive labelled tournaments n<=5; symbolic proofs in ROOT_COST.md',
 census=rows,skew_controls=special,
 dynkin_independence_at2={f'A{n}':independence(n,[(i,i+1) for i in range(n-1)]) for n in range(1,7)},
 affine_A2_independence_at2=independence(3,[(0,1),(1,2),(0,2)]))
print(json.dumps(out,indent=2))
