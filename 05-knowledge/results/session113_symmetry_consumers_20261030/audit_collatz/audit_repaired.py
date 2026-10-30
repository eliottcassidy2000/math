#!/usr/bin/env python3
"""Independent NumPy integer checker of the five-mode repaired certificate.
All verification gates survive python -O. No producer/checker code imports.
"""
import argparse,hashlib,json,time
from pathlib import Path
import numpy as np

def require(test,detail):
 if not bool(test):raise ValueError(detail)
p=argparse.ArgumentParser();p.add_argument('E',type=int);p.add_argument('vector',type=Path);p.add_argument('--output',type=Path);p.add_argument('--chunk',type=int,default=100000);a=p.parse_args()
start=time.time();E=a.E;require(1<=E<=9,'supported E');q=5**E;D=8*q;N=5*q-5;cap=64;NEG=-(1<<60)
require(a.vector.stat().st_size==5*q*8,'vector byte length')
V=np.memmap(a.vector,dtype='<i8',mode='r',shape=(5*q,))
z=np.arange(q,dtype=np.int64);h=np.zeros(q,dtype=np.int64)
for e in range(1,E+1):h+=(z%5**e==0)
# Independent high-layer density = sum inverse periods at levels2..h.
rates=np.zeros(E+1,dtype=np.int64)
for e in range(2,E+1):rates[e]=rates[e-1]+D//(2*5**(e-1))
valid=np.concatenate((np.arange(4*q,dtype=np.int64),4*q+np.arange(0,q,5,dtype=np.int64)))
require(len(valid)==4*q+q//5,'valid cardinality')
require(np.all(V[valid]!=NEG),'missing valid-state potential')
require(np.all(np.abs(V[valid])<2**50),'safe exact signed64 arithmetic')
excluded=4*q+z[z%5!=0];require(np.all(V[excluded]==NEG),'unexpected excluded-state payload')
mn=int(np.min(V[valid]));mx=int(np.max(V[valid]));require(mx<=0,'nonpositive maximum')
least=None;edges=0
for letter in range(1,cap+1):
 inv=pow(2**letter,-1,q)
 for lo in range(0,len(valid),a.chunk):
  src=valid[lo:lo+a.chunk];mode=src//q;oldz=src%q;u=np.where(mode==4,4,mode+1)
  newu=(3*u*pow(2**letter,-1,5))%5
  newz=((3*oldz+pow(2,letter,q)-2)*inv)%q
  H=h[newz];pending=(letter>1)&(newu==4)&(H>0)
  newmode=np.where(pending,4,newu-1);dst=newmode*q+newz
  require(np.all((newmode!=4)|(newz%5==0)),('closure',letter,lo))
  require(np.all(V[dst]!=NEG),('destination finite',letter,lo))
  if letter>1:gross=np.where(newu==1,D*H,0)
  else:
   gross=np.where(mode==4,D*H,np.where(((newu==1)|(newu==4))&(H>0),np.where(newu==1,D,0)+rates[H],0))
  cost=gross-N*letter
  slack=V[dst]-V[src]-cost
  require(np.all(slack>=0),('edge violation',letter,int(src[np.argmin(slack)]),int(slack.min())))
  sl=int(slack.min());least=sl if least is None else min(least,sl);edges+=len(src)
 # Entry handled independently rather than selecting rows from edge arrays.
 if letter>1:
  for u in (1,4):
   nu=3*u*pow(2**letter,-1,5)%5;nz=((2**letter-2)*inv)%q;H=int(h[nz])
   mode=4 if nu==4 and H>0 else nu-1;dest=mode*q+nz
   c=(D*H if nu==1 else 0)-N*letter
   require(int(V[dest])>=c,('entry',u,letter))
require(edges==(4*q+q//5)*cap,'coverage')
et=mn-(D*E-65*N);tt=mn-(mx+D*E-65*N)
require(et>=0,'entry tail');require(tt>=0,'transition tail')
r=dict(E=E,q=q,valid_states=len(valid),excluded_pending_height_zero_states=len(excluded),edge_inequalities=edges,entry_inequalities=126,letters_checked=cap,minimum_V=mn,maximum_V=mx,minimum_edge_slack=least,entry_tail_margin=et,transition_tail_margin=tt,sha256=hashlib.sha256(a.vector.read_bytes()).hexdigest(),seconds=time.time()-start,scope='all finite positive words at this fixed E via repaired one-digit epoch envelope; not all E or convergence')
s=json.dumps(r,indent=2)+'\n';print(s,end='')
if a.output:a.output.write_text(s)
