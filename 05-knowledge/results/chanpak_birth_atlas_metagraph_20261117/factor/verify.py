"""Exact small physical census and independent min-plus word controls."""
from itertools import permutations
from math import inf,isqrt
from functools import cache
import json

def require(b):
 if not b:raise ArithmeticError('control failed')

def induced(rows,p):return tuple(sum(((rows[v]>>w)&1)<<j for j,w in enumerate(p)) for v in p)
def canon(rows):return min(induced(rows,p) for p in permutations(range(len(rows))))
def converse(rows):return tuple(sum(((rows[j]>>i)&1)<<j for j in range(len(rows))) for i in range(len(rows)))
def H(rows):
 n=len(rows);dp=[[0]*n for _ in range(1<<n)]
 for j in range(n):dp[1<<j][j]=1
 for s in range(1,1<<n):
  for v in range(n):
   if s>>v&1:
    r=s^(1<<v)
    if r:dp[s][v]=sum(dp[r][u] for u in range(n) if r>>u&1 and rows[u]>>v&1)
 return sum(dp[-1])
def strong(rows):
 for v in range(len(rows)):
  r=1<<v
  while True:
   q=r
   for j in range(len(rows)):
    if r>>j&1:q|=rows[j]
   if q==r:break
   r=q
  if r!=(1<<len(rows))-1:return False
 return True

def decoder(atoms):
 # atoms[d]=(minimum strong order, multiplicity, self-converse multiplicity)
 @cache
 def birth(h):
  if h==1:return 0,1
  best=inf;ways=0
  for d,(s,v,c) in atoms.items():
   if h%d:continue
   tail,w=birth(h//d);cost=s+tail
   if cost<best:best=cost;ways=v*w
   elif cost==best and cost<inf:ways+=v*w
  return best,ways
 def fixed(h):
  b,a=birth(h);ans=0
  if b==inf:return 0
  for d in range(1,isqrt(h)+1):
   if h%(d*d):continue
   g=h//(d*d);s,v,c=(0,1,1) if g==1 else atoms.get(g,(inf,0,0))
   bd,ad=birth(d)
   if 2*bd+s==b:ans+=ad*c
  return ans
 return birth,fixed

observed={};atoms={};labelled=0
for n in range(1,6):
 reps=set()
 for bits in range(1<<(n*(n-1)//2)):
  rows=[0]*n;i=0
  for u in range(n):
   for v in range(u+1,n):
    if bits>>i&1:rows[u]|=1<<v
    else:rows[v]|=1<<u
    i+=1
  reps.add(canon(rows));labelled+=1
 for rows in reps:
  h=H(rows);f=int(canon(converse(rows))==rows)
  if h>1:
   if h not in observed:observed[h]=[n,0,0]
   if observed[h][0]==n:observed[h][1]+=1;observed[h][2]+=f
   if strong(rows):
    if h not in atoms:atoms[h]=[n,0,0]
    if atoms[h][0]==n:atoms[h][1]+=1;atoms[h][2]+=f
birth,fixed=decoder(atoms)
for h,(n,a,f) in observed.items():require(birth(h)==(n,a) and fixed(h)==f)
# This artificial dictionary is only an algebraic test, not physical claims.
b,f=decoder({3:(3,1,1),5:(4,1,1),9:(6,1,1)})
require(b(9)==(6,2) and f(9)==2)
require(b(27)==(9,3) and f(27)==1)
require(b(15)==(7,2) and f(15)==0)
print(json.dumps({'labelled_tournaments':labelled,'physical_first_births':observed,'synthetic_27':{'order':9,'classes':3,'self_converse':1},'synthetic_squarefree_15':{'order':7,'classes':2,'self_converse':0}},indent=2))
