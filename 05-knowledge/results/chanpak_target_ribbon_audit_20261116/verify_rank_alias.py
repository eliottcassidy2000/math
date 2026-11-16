#!/usr/bin/env python3
"""Rebuild the literal span-three system; certify a suffix/rank alias over Q.
No floating-point or probabilistic rank claims are used in the final certificate.
The modular sieve only chooses candidate independent rows and columns.
"""
import sys,json
from pathlib import Path
from collections import deque
import numpy as np
from sympy import Matrix
from cap_check import transfer_system,need
states,_,raw,final=transfer_system(3)
d=len(states); prime=10007
h=np.array([int(i in final) for i in range(d)],dtype=np.int64)
A=[]
for mask in range(8):
 a=frozenset(i+1 for i in range(3) if mask>>i&1)
 mat=np.zeros((d,d),dtype=np.int64)
 for (i,j),c in raw[a].items():mat[i,j]=c
 A.append(mat)
class Selector:
 def __init__(self): self.rows=[];self.piv=[]
 def add(self,v):
  v=np.array(v,dtype=np.int64)%prime
  for p,b in zip(self.piv,self.rows): v=(v-v[p]*b)%prime
  nz=np.flatnonzero(v)
  if not len(nz):return False
  p=int(nz[0]);v=v*pow(int(v[p]),-1,prime)%prime
  self.piv.append(p);self.rows.append(v);return True

def closure(start,step):
 S=Selector();S.add(start);chosen=[np.array(start,dtype=np.int64)];words=[()];front=deque([0])
 while front:
  k=front.popleft()
  for a in range(8):
   v=step(chosen[k],A[a])
   if S.add(v):
    words.append(words[k]+(a,));chosen.append(v);front.append(len(chosen)-1)
 return chosen,words

# Full observable row span, and reachable column span AFTER the initial vertex.
obs,ow=closure(h,lambda v,B:v@B)
v0=A[0][:,0]
reach,rw=closure(v0,lambda v,B:B@v)
O=Matrix(np.array(obs));B=Matrix(np.array(reach).T)
print('sieved dimensions',O.shape,B.shape,flush=True)
def rank(X):return len(X.to_DM().rref()[1])
def piv(X):return X.to_DM().rref()[1]
need(rank(O)==len(obs) and rank(B)==len(reach))
# A modular sieve cannot prove closure; these identities do so over Q.
# Pivot-column right inverse for O, pivot-row left inverse for B.
op=piv(O); Ocols=O[:,list(op)]; Oinv=Ocols.inv(method="DM")
bp=piv(B.T); Brows=B[list(bp),:]; Binv=Brows.inv(method="DM")
for a in range(8):
 AA=Matrix(A[a]); OA=O*AA; AB=AA*B
 need(OA[:,list(op)]*Oinv*O == OA)
 need(B*Binv*AB[list(bp),:] == AB)
need(Ocols.det()!=0 and Brows.det()!=0)
rank_pair=rank(O*B)
need(rank_pair==61)
print('exact observable dimension:',len(obs),'reachable-after-first:',len(reach),'minimal:',rank_pair,flush=True)

def product(w):
 P=Matrix.eye(d)
 for a in w:P=Matrix(A[a])*P
 return P
u=(0,4,1,7,0,2);v=(0,7,0,2);tail=0
result={'alphabet':'bit mask of backward gaps 1,2,3','u':u,'v':v,'common_suffix':u[-3:],'appended_letter':tail,'minimal_dimension':rank_pair,'ranks':{}}
for name,w in [('u',u),('v',v),('u0',u+(tail,)),('v0',v+(tail,))]:
 K=O*product(w)*B
 r=rank(K);result['ranks'][name]=r
 print(name,w,'rank',r,flush=True)
 # Retain compact explicit nonsingular minor and exact rank factorization.
 cs=piv(K);C=K[:,list(cs)];rs=piv(C.T)
 D=C[list(rs),:]
 need(D.det()!=0)
 F=D.inv()*K[list(rs),:]
 need(C*F==K)
 result[name+'_minor']={'rows':list(rs),'cols':list(cs),'det':str(D.det())}
need(result['ranks']=={'u':3,'v':3,'u0':3,'v0':2})
need(json.dumps(result,indent=2)+'\n'==Path(__file__).with_name('RANK_ALIAS_CERTIFICATE.json').read_text())
print('PASS: same last three letters and same rank do not determine next rank.',flush=True)
