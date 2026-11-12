"""Optional compatibility-first candidate discovery; requires NumPy/SciPy.
Numerical output must pass a separate rational checker before any proof claim.
"""
from verify import block,mv,R,ROOT
from itertools import product
import numpy as np
from pathlib import Path
import argparse
class V:
    block=staticmethod(block)
    mv=staticmethod(mv)
    RAYS=R
v=V()
from scipy.optimize import linprog
import json
j=json.loads((ROOT/'metric151.json').read_text());words=j['words'][:56];Qs=[np.array([[a,b],[b,c]],dtype=float)/j['Q_den'] for a,b,c in j['Q_num'][:56]];mats=list(map(v.block,words));rho=.999
basis=[np.array([[1.,0],[0,0]]),np.array([[0.,1],[1,0]]),np.array([[0.,0],[0,1]])]
def Ks(U,V):
 q=V[0][1];t=V[1][1];ans=[]
 for r in v.RAYS:
  a,w,d=v.mv(V,v.mv(U,r));rr=w/a;s=d/a;ans.append(np.array([[t-q*rr,rr],[-q*s,s]]))
 return ans
def attempt(V):
 Ds=[]
 for U,Q in zip(mats,Qs):
  for K in Ks(U,V):Ds.append((rho*Q,[-K.T@E@K for E in basis]))
  for K in Ks(V,U):Ds.append((-K.T@Q@K,[rho*E for E in basis]))
 for K in Ks(V,V):Ds.append((np.zeros((2,2)),[rho*E-K.T@E@K for E in basis]))
 Ds.append((-.5*np.eye(2),basis))
 cons=[];bs=[]
 def add(D,z):
  A,B=D;cons.append([-z@E@z for E in B]);bs.append(z@A@z-1e-6)
 for D in Ds:
  for z in [np.array([1,0]),np.array([0,1]),np.array([1,1]),np.array([1,-1])]:add(D,z)
 for it in range(30):
  res=linprog([1,0,1],A_ub=cons,b_ub=bs,bounds=[(None,None)]*3,method='highs')
  if not res.success:return None
  bad=0
  for D in Ds:
   A,B=D;T=A+sum(x*E for x,E in zip(res.x,B));val,vec=np.linalg.eigh(T)
   if val[0]<-1e-8:add(D,vec[:,0]);bad+=1
  if not bad:return np.array([[res.x[0],res.x[1]],[res.x[1],res.x[2]]])
 return None
for n in range(8,10):
 acc=[]
 for mid in product('01',repeat=n-2):
  w='1'+''.join(mid)+'0'
  if sum(a!=b for a,b in zip(w,w[1:]))>3:continue
  M=v.block(w);Q=attempt(M)
  if Q is not None:words.append(w);mats.append(M);Qs.append(Q);acc.append(w)
 print(n,len(acc),len(words),acc,flush=True)
 json.dump({'status':'NUMERICAL CANDIDATE, NOT CERTIFIED','words':words,'Q':[q.tolist() for q in Qs],'rho2':rho},open(str(ROOT/'greedy_candidate.json'),'w'))
