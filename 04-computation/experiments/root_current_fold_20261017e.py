#!/usr/bin/env python3
"""Exact odd-source folding of native ROOT currents; no ROOT oracle.

Measures the two components together, never assumes their independence.
All integer clocks are authenticated by the inherited full edge audit.
"""
import argparse,json,sys,hashlib
from pathlib import Path
from itertools import product
import numpy as np
from collatz_compact_resets_20261015f import Graph
from collatz_clock_tail_20261017b import audit


def need(ok,msg):
 if not ok:raise RuntimeError(msg)


def chi(n):return np.where(n%3==0,0,np.where(n%3==1,1,-1))
def weight(n):return np.where(n%3==2,2,-1)
def digest(a):return hashlib.sha256(a.tobytes()).hexdigest()
def energy(a):return int(np.dot(a,a))


def fold(M,t):
 n=np.arange(M+1,dtype=np.int64);odd=(n&1)!=0
 jp=np.cumsum(np.where(odd,chi(n)*t,0));op=np.cumsum(np.where(odd,weight(n)*t,0))
 rows=np.arange(1,M//4+1,dtype=np.int64)
 inside=-3*jp[2*rows];rim=op[4*rows]-op[2*rows]
 weighted=np.cumsum(weight(n)*t)
 actual=weighted[4*rows]-weighted[rows]
 need(np.array_equal(actual,inside+rim),'exact odd-source folding on every prefix row')
 # Independently enumerate only the three source events in the increment.
 delta=weight(4*rows-3)*t[4*rows-3]+weight(4*rows-1)*t[4*rows-1]
 m=2*rows-1;delta+=(-3*chi(m)-weight(m))*t[m]
 need(np.array_equal(actual,np.cumsum(delta)),'three marked odd-source events')
 need(int(np.max(np.abs(delta)))<=4,'four-unit increment bound')
 return inside,rim,actual


def controls():
 count=0;M=20
 for bits in product((-1,0,1),repeat=6):
  t=np.zeros(M+1,dtype=np.int64)
  for n in range(1,M+1):
   m=n;a=0
   while m%2==0:m//=2;a+=1
   if m<=11:t[n]=(-1)**a*bits[(m-1)//2]
  fold(M,t);count+=1
 n=np.arange(129,dtype=np.int64);t=chi(n);t[0]=0
 _,_,actual=fold(128,t)
 need(np.array_equal(actual,-3*np.arange(1,33)),'ROOT-normalized hostile f=chi3 gives R(N)=-3N')
 return dict(exhaustive_profiles=count,hostile='f(n)=chi3(n), f(1)=1, f(2n)=-f(n), R(N)=-3N; violates odd Collatz rule')


def run(k,stages):
 need(4<=k<=18,'int64 energy envelope M=8X and k<=18')
 X=1<<k;g=Graph(8*X);out=[]
 for stage in range(1,max(stages)+1):
  cols,inv=np.unique(g.endpoint[1:],return_inverse=True);g.extend(cols);g.close(cols,inv)
  if stage not in stages:continue
  rec=audit(g)
  t=np.where(g.endpoint==1,np.where(g.clock%2==0,1,-1),0).astype(np.int64);t[0]=0
  need(np.array_equal(t[2::2],-t[1:4*X+1]),'native source profile is dyadically compatible')
  I,B,R=fold(g.M,t);sl=slice(X-1,2*X-1);I,B,R=I[sl],B[sl],R[sl]
  EI,EB,ER=energy(I),energy(B),energy(R);cross=int(np.dot(I,B))
  need(ER==EI+EB+2*cross,'exact covariance identity')
  # A marked binary residue decomposition; all pair moments retained.
  values=[]
  for r in (1,3,5):
   n=np.arange(g.M+1,dtype=np.int64);pref=np.cumsum(np.where(n%6==r,t,0));q=np.arange(X,2*X,dtype=np.int64)
   values.append(pref[4*q]-pref[2*q])
  G=[[int(np.dot(a,b)) for b in values] for a in values]
  need(EB==sum(weight(np.array([r]))[0]*weight(np.array([s]))[0]*G[i][j] for i,r in enumerate((1,3,5)) for j,s in enumerate((1,3,5))),'signed residue Gram reconstructs rim')
  out.append(dict(stage=stage,audit=rec,root_sources=int(np.count_nonzero(t)),
   interior_energy=EI,rim_energy=EB,cross=cross,root_energy=ER,
   normalized=dict(interior=EI/X**2,rim=EB/X**2,cross=cross/X**2,root=ER/X**2),
   signed_root_range=[int(R.min()),int(R.max())],residue_rim_gram=G,
   full_profile_sha256=digest(t),folded_rows_sha256=digest(R)))
  print(f'k={k} stage={stage}: E(R)/X²={ER/X**2:.6f}',file=sys.stderr,flush=True)
 return dict(k=k,X=X,rows=X,stages=out)

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--k',nargs='+',type=int,default=[10,14,18]);p.add_argument('--stages',nargs='+',type=int,default=[1,2,4,8,12,20,28]);p.add_argument('--output',required=True);a=p.parse_args()
 d=dict(status='FINITE-EXACT native source-current decomposition; uniform bound OPEN',controls=controls(),runs=[])
 for k in a.k:
  d['runs'].append(run(k,a.stages));Path(a.output).write_text(json.dumps(d,sort_keys=True,separators=(',',':'))+'\n')
