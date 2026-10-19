from itertools import permutations
from math import comb
import sympy as s
import numpy as np
import json
from pathlib import Path
from probe_norm import rref_mod
PERMS=list(permutations(range(3)))
def comps(n):return [(a,b,n-a-b) for a in range(n+1) for b in range(n-a+1)]
def permute(t,pm):return tuple(t[i] for i in pm)+tuple(t[i+3] for i in pm)
def orbs(d,k,eps):
 seen=set();out=[]
 for a in comps(d+3):
  for b in comps(k):
   if any(bi==0 and ai<2 for ai,bi in zip(a,b)):continue
   t=a+b
   if t in seen:continue
   terms={};bad=False
   for pm in PERMS:
    v=permute(t,pm);inv=sum(pm[i]>pm[j] for i in range(3) for j in range(i+1,3));c=eps**inv
    if v in terms and terms[v]!=c:bad=True
    terms[v]=c
   seen.update(terms)
   if not bad:out.append(terms)
 return out

def a_monom(t,j):
 a,b,c,A,B,C=t;out={}
 if j==1:
  for u in range(b+1):
   for v in range(B+1):out[(u,a+b-u,c,v,A+B-v,C)]=comb(b,u)*comb(B,v)*(-1)**(a+b-u+A+B-v)
 else:
  for u in range(a+1):
   for v in range(A+1):out[(u+b,a-u,c,v+B,A-v,C)]=comb(a,u)*comb(A,v)*(-1)**(u+b+v+B)
 return out

def add(out,t,c):
 out[t]=out.get(t,0)+c
 if out[t]==0:del out[t]
def shift(t,j):return tuple(v+(i==j) for i,v in enumerate(t))
def cert_column(o):
 out={}
 for t,c in o.items():
  add(out,shift(t,0),c);add(out,shift(t,1),-c)
  for u,e in a_monom(t,1).items():add(out,shift(u,0),-c*e)
  for u,e in a_monom(t,2).items():add(out,shift(u,1),c*e)
 return out

def getmatrix(d,k,eps):
 o=orbs(d,k,eps);cols=[cert_column(v) for v in o];keys=sorted(set().union(*(v for v in cols))) if cols else [];idx={v:i for i,v in enumerate(keys)}
 mat=np.zeros((len(keys),len(cols)),dtype=np.int64)
 for j,v in enumerate(cols):
  for t,c in v.items():mat[idx[t],j]=c
 return mat,o
if __name__=='__main__':
 for k,maxd in [(2,12),(3,3)]:
  for d in range(maxd+1):
   for eps in (1,-1):
    a,o=getmatrix(d,k,eps);rank=len(rref_mod(a,100003)[1]);null=len(o)-rank
    print(d,k,eps,a.shape,'null mod100003',null,flush=True)
    if null:
     b=s.polys.matrices.DomainMatrix.from_Matrix(s.Matrix(a.tolist())).nullspace().to_Matrix()
     print('rational null',b.shape,flush=True)
     arr=[]
     for v in b.tolist():
      den=s.ilcm(*[s.denom(t) for t in v]);v=[int(t*den) for t in v];g=s.igcd(*v);v=[t//g for t in v];pol={}
      for j,oo in enumerate(o):
       for t,c in oo.items():
        if v[j]:pol[','.join(map(str,t))]=int(c*v[j])
      arr.append(pol)
     (Path(__file__).resolve().parent / f'rational_kernel_d{d}_k{k}_eps{eps}.json').write_text(json.dumps(arr,indent=2)+'\n')
