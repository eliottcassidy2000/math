#!/usr/bin/env python3
"""Independent audit; no producer imports. Permutations rather than subset DP."""
from itertools import product,permutations
from fractions import Fraction as F
from pathlib import Path
import json,sys,hashlib
ROOT=Path(__file__).resolve().parents[1]
def require(x,m):
 if not x:raise RuntimeError(m)
def words(a,n):
 for k in range(n+1):
  for w in product(a,repeat=k):yield ''.join(w)
def graph(w):
 arcs={(0,1)};r=0;n=2
 for c in w:
  arcs.update((n,v) if c in 'ZJ' and v==r else (v,n) for v in range(n))
  if c=='J':r=n
  n+=1
 return n,arcs,r
def h(n,arcs,omit=None,end=None):
 return sum(all((p[i],p[i+1]) in arcs for i in range(len(p)-1)) and (end is None or p[-1]==end) for p in permutations([v for v in range(n) if v!=omit]))
def move(s,c):
 D,E,U,d=s
 # Independent positive-coordinate realization inherited from Session199.
 a,f,e,p=D-E-d,U-D+E-d,E,d-E
 if c=='O':a,f,e,p=a+e,a+2*f+e+2*p,0,e+p
 elif c=='Z':a,f,e,p=2*a+f+e+p,f,e+p,0
 else:a,f,e,p=f+p,0,a+e+p,e
 return a+2*e+p,e,a+f+2*e+2*p,e+p
def state(w):
 s=(1,0,2,1)
 for c in w:s=move(s,c)
 return s
if '--fail' in sys.argv:require(False,'intentional audit failure')
profiles=0
for w in words('OZJ',4):
 n,a,r=graph(w);D=h(n,a);E=h(n,a,end=r);d=h(n,a,omit=r)
 nn,aa,rr=graph(w+'Z');U=h(nn,aa)-D+E
 require(state(w)==(D,E,U,d),('graph',w))
 profiles+=1
trajectories=0
for w in words('OZJ',7):
 s=(1,0,2,1);ds=[]
 for c in w:
  s=move(s,c);D,E,U,d=s
  require(0<=E<=d<=D<=U and U>=2*d,'cone')
  if c=='J':
   ds.append(d)
   if len(ds)>2:require(ds[-1]>=ds[-2]+ds[-3],'Fibonacci')
 require(len(ds)<=2*(s[3].bit_length()-1)+2,'log bound')
 trajectories+=1
shadow=0
for root in ('','J','ZJ','OJJ','JOZJ','ZZOJ','JJOJ'):
 s0=state(root);D0,E0,U0,d=s0;a=U0-d
 for w in words('OZ',8):
  s=state(root+w);n=len(w);A=sum(2**j for j,c in enumerate(w) if c=='Z')
  y=(F(D0,U0)+A)/2**n
  require(s[2]>=d+2**n*a and abs(F(s[0],s[2])-y)<=F(n*d,2**n*a),'shadow')
  shadow+=1
# Full rational row reduction independently checks shortest rank-one resets.

def matrix(w):
 cols=[]
 for j in range(4):
  s=tuple(int(i==j) for i in range(4))
  for c in w:s=move(s,c)
  cols.append(s)
 return [list(r) for r in zip(*cols)]
def rank(m):
 a=[[F(v) for v in row] for row in m];r=0
 for j in range(4):
  p=next((p for p in range(r,4) if a[p][j]),None)
  if p is None:continue
  a[r],a[p]=a[p],a[r];q=a[r][j];a[r]=[x/q for x in a[r]]
  for k in range(r+1,4):
   q=a[k][j];a[k]=[x-q*y for x,y in zip(a[k],a[r])]
  r+=1
 return r
resets=[w for w in words('OZJ',4) if rank(matrix(w))==1]
require(resets==['OJOJ','OJJZ','JZOJ','JZJZ'],'rank resets')
require(matrix('OJOJ')==[[3,0,6,0],[1,0,2,0],[3,0,6,0],[1,0,2,0]],'reset identity')
reset_grid=0
for m in range(6):
 for k in range(1,7):
  v=(1+2**k,1,1+2**k,1);lam=(1,0,2**(m+1),0)
  for prefix in ('OJ','JZ'):
   require(matrix(prefix+'Z'*m+'O'*k+'J')==[[x*y for y in lam] for x in v],'general reset')
   reset_grid+=1
n,arcs,r=graph('OJOJ')
require(n==6 and [sum(1<<v for u,v in arcs if u==i) for i in range(n)]==[54,60,56,17,32,8],'physical placement')
require(h(n,arcs)==15,'physical reset count')
seals={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for lane in ('descent','bridges') for p in sorted((ROOT/lane).glob('*')) if p.suffix in ('.md','.py')}
print(json.dumps(dict(status='PASS',permutation_profiles=profiles,ternary_trajectories=trajectories,closed_formula_shadow_tests=shadow,shortest_rank_one_words=resets,reset_grid=reset_grid,source_sha256=seals,scope='structural checks only; no coverage expansion'),indent=2,sort_keys=True))
