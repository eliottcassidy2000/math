#!/usr/bin/env python3
"""Two exact source folds, a sharp linear hostile, and native signed balance.

No unobserved trajectory receives a ROOT label. Native endpoints/clocks are
authenticated in full by the inherited observed-edge audit.
"""
import argparse,json,hashlib
from pathlib import Path
from fractions import Fraction
import numpy as np
from collatz_compact_resets_20261015f import Graph
from collatz_clock_tail_20261017b import audit
from root_current_fold_20261017e import fold,need


def chi(n):return 0 if n%3==0 else 1 if n%3==1 else -1
def w(n):return 2 if n%3==2 else -1
def seed(n):
 while n%8==5:n=(n-1)//4
 return n
def coefficient(N,s):
 out=0;m=s
 while m<=4*N:
  out+=-3*chi(m) if m<=2*N else w(m);m=4*m+1
 return out
def hostile(M,normalize=False):
 odd=np.arange(1,M+1,2,dtype=np.int64);seeds=odd.copy()
 while True:
  mask=seeds%8==5
  if not np.any(mask):break
  seeds[mask]=(seeds[mask]-1)//4
 F=np.where(seeds%3==2,1,-1)
 if normalize:F[seeds==1]=1
 t=np.zeros(M+1,dtype=np.int64);t[1::2]=F
 a=1;q=2
 while q<=M:
  t[q::2*q]=(-1)**a*F[:len(t[q::2*q])];q*=2;a+=1
 return t
def controls():
 count=0;observed={0:set(),1:set(),2:set()}
 for N in range(1,257):
  total=0;value=0
  for s in range(1,4*N+1,2):
   if s%8==5:continue
   c=coefficient(N,s);observed[s%3].add(c)
   need(c*w(s)>=0 and abs(c)<=3,'fixed column sign and sharp three-unit coefficient')
   total+=abs(c);value+=c*(1 if s%3==2 else -1);count+=1
  need(value==total,'one fixed profile attains every row norm')
  t=hostile(4*N);need(sum(w(n)*int(t[n]) for n in range(N+1,4*N+1))==total,'independent direct source sum')
 # Mean derived from the disjoint (a,j,seed mod24) progression partition.
 mean=sum(Fraction(((-1)**a)*(1 if r==2 else -1)*w((2**a)*(r+j)),8*2**a*4**j)
          for a in range(2) for j in range(3) for r in range(3))/((1-Fraction(1,4))*(1-Fraction(1,64)))
 need(mean==Fraction(16,21),'exact progression density')
 return dict(coefficient_checks=count,coefficient_values={str(k):sorted(v) for k,v in observed.items()},prefix_mean=str(mean),annular_linear_coefficient='16/7')


def run(k,stages):
 X=1<<k;M=8*X;g=Graph(M);bad=hostile(M);badroot=hostile(M,True)
 _,_,A=fold(M,bad);_,_,Ar=fold(M,badroot);sl=slice(X-1,2*X-1);A=A[sl];Ar=Ar[sl]
 need(np.all(A>=0) and np.max(np.abs(A-Ar))<=6,'normalizing ROOT changes only one folded column')
 rows=np.arange(X,2*X,dtype=np.int64);err=7*A-16*rows
 out=[]
 for stage in range(1,max(stages)+1):
  cols,inv=np.unique(g.endpoint[1:],return_inverse=True);g.extend(cols);g.close(cols,inv)
  if stage not in stages:continue
  receipt=audit(g);t=np.where(g.endpoint==1,np.where(g.clock%2==0,1,-1),0).astype(np.int64);t[0]=0
  src=np.arange(1,(M-1)//4+1,2,dtype=np.int64)
  need(np.array_equal(t[4*src+1],t[src]),'braid equality for every native pair in the universe')
  need(np.all(g.seen[2:]),'every interior nonROOT edge observed: defect is outside the source box')
  _,_,R=fold(M,t);_,_,cert=fold(M,bad*np.abs(t));R=R[sl];cert=cert[sl]
  need(np.all(cert>=np.abs(R)) and np.all(cert<=A),'signed certified mass and unsigned source envelope')
  need(np.all((cert+R)%2==0),'positive/negative signed masses are integers')
  out.append(dict(stage=stage,audit=receipt,root_energy=int(R@R),root_energy_over_X2=float(R@R)/X**2,
   certified_mass_range=[int(cert.min()),int(cert.max())],signed_current_range=[int(R.min()),int(R.max())],
   max_relative_imbalance=float(np.max(np.abs(R)/np.maximum(cert,1))),
   positive_mass_sum=int(((cert+R)//2).sum()),negative_mass_sum=int(((cert-R)//2).sum()),
   source_profile_sha256=hashlib.sha256(t.tobytes()).hexdigest(),current_sha256=hashlib.sha256(R.tobytes()).hexdigest()))
 return dict(k=k,X=X,hostile_linear_error_numerator_range=[int(err.min()),int(err.max())],
  hostile_current_over_N_range=[float(np.min(A/rows)),float(np.max(A/rows))],stages=out)


if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--k',nargs='+',type=int,default=[10,14,18]);p.add_argument('--stages',nargs='+',type=int,default=[1,4,12,28]);p.add_argument('--output',required=True);a=p.parse_args()
 d=dict(status='PROVED exact braid fold and asymptotic operator norm; FINITE-EXACT native balance; uniform ROOT bound OPEN',controls=controls(),runs=[])
 for k in a.k:
  d['runs'].append(run(k,a.stages));Path(a.output).write_text(json.dumps(d,sort_keys=True,separators=(',',':'))+'\n')
  print(json.dumps(d['runs'][-1],sort_keys=True),flush=True)
