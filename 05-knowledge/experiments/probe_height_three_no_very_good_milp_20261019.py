#!/usr/bin/env python3
import sys,json,time
from itertools import combinations
import numpy as np
from scipy.optimize import milp,Bounds,LinearConstraint
from scipy.sparse import coo_matrix
n=int(sys.argv[1]);limit=float(sys.argv[2]) if len(sys.argv)>2 else 120
pairs=list(combinations(range(n),2));index={p:i for i,p in enumerate(pairs)}
lo=[0]*len(pairs);hi=[1]*len(pairs);rows=[];rl=[];ru=[]
def variable(upper=1):
 i=len(lo);lo.append(0);hi.append(upper);return i
rank=[variable(2) for _ in range(n)]
def constraint(d,l=-np.inf,u=np.inf):rows.append(d);rl.append(l);ru.append(u)
def rel(a,b):return (index[a,b],True) if a<b else False
def neg(t):return not t if isinstance(t,bool) else (t[0],not t[1])
def witness(lits):
 normalized={}
 for t in lits:
  if t is False:return None
  if t is True:continue
  k,s=t
  if k in normalized and normalized[k]!=s:return None
  normalized[k]=s
 if not normalized:return True
 w=variable()
 for k,s in normalized.items():constraint({w:1,k:-1 if s else 1},u=0 if s else 1)
 return w
for i,j in pairs:
 constraint({rank[j]:1,rank[i]:-1,index[i,j]:-3},l=-2)
for i,j,k in combinations(range(n),3):constraint({index[i,j]:1,index[j,k]:1,index[i,k]:-1},u=1)
for i in range(n):constraint({index[min(i,j),max(i,j)]:1 for j in range(n) if i!=j},l=n-8)
constraint({v:1 for v in index.values()},u=len(pairs)-1)
# Exclude a very good pair in each order orientation.
for dual in (False,True):
 def D(a,k):return rel(a,k) if dual else rel(k,a)
 def U(a,k):return rel(k,a) if dual else rel(a,k)
 for a,b in pairs:
  options={index[a,b]:1}
  # Unequal strict predecessor sets, or a nonchain directed successor difference.
  for k in range(n):
   for x,y in [(a,b),(b,a)]:
    w=witness([D(x,k),neg(D(y,k))])
    if w is not None:options[w]=1
  for x,y in [(a,b),(b,a)]:
   for k,l in pairs:
    w=witness([U(x,k),U(x,l),neg(U(y,k)),neg(U(y,l)),neg(rel(k,l))])
    if w is not None:options[w]=1
  constraint(options,l=1)
ii=[];jj=[];vv=[]
for i,row in enumerate(rows):
 for j,v in row.items():ii.append(i);jj.append(j);vv.append(v)
A=coo_matrix((vv,(ii,jj)),shape=(len(rows),len(lo))).tocsc()
start=time.time();print(json.dumps({'n':n,'variables':len(lo),'constraints':len(rows),'limit':limit}),flush=True)
r=milp(np.zeros(len(lo)),integrality=np.ones(len(lo)),bounds=Bounds(lo,hi),constraints=LinearConstraint(A,rl,ru),options={'time_limit':limit,'mip_rel_gap':0})
out={'n':n,'status':int(r.status),'message':r.message,'seconds':time.time()-start,'variables':len(lo),'constraints':len(rows)}
if r.x is not None:
 values=np.rint(r.x).astype(int);p=[sum(1<<i for i in range(j) if values[index[i,j]]) for j in range(n)];out['pred']=p;out['rank']=[int(values[v]) for v in rank]
open(sys.argv[3] if len(sys.argv)>3 else f'probe_{n}.json','w').write(json.dumps(out,indent=2));print(json.dumps(out),flush=True)
