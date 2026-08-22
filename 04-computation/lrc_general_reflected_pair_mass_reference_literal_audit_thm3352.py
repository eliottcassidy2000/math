#!/usr/bin/env python3
"""Deterministic hostile/random literal audit of the corrected mass engine.

Audit scope: L >= 168 with 14 | L, integral cells, distinct endpoint labels
1 <= e,f <= 14, and positive integral levels p,q.
"""
from fractions import Fraction as F
from hashlib import sha256
from importlib.util import spec_from_file_location,module_from_spec
from pathlib import Path
import os,random

ROOT=Path(os.environ.get('LRC_REPO',Path(__file__).resolve().parents[1])).resolve()
ENGINE=ROOT/'04-computation/lrc_general_reflected_pair_mass_reference_audit_thm3352.py'
EXPECTED_ENGINE='da941a4267147d5442be81ae81880742d2f6b901bfc1d20fb667822402a2950e'
EXPECTED_CHECKS=1044591
EXPECTED_SEMANTIC='b960b450d70fbfb5521b3d90eb4d8064fcb126fe7be4fd81f3c2beebbbffa4a0'

def load(path):
 s=spec_from_file_location('mass_engine_audit',path);m=module_from_spec(s);s.loader.exec_module(m);return m
if sha256(ENGINE.read_bytes().replace(b'\r\n',b'\n')).hexdigest()!=EXPECTED_ENGINE:
 raise RuntimeError(('engine hash',ENGINE))
M=load(ENGINE)

def clip(a,b):a=max(F(),a);b=min(F(1),b);return (a,b) if a<b else None
def arcs(L,e,q,j):
 z=q*L-e;r=e*j%L;rad=F(L,14*z);out=[]
 for k in range((-r)//L-1,(z-r)//L+2):
  c=F(r+k*L,z);x=clip(c-rad,c+rad)
  if x:out.append(x)
 return tuple(sorted(set(out)))
def inter(A,B):
 i=j=0;t=F()
 while i<len(A) and j<len(B):
  t+=max(F(),min(A[i][1],B[j][1])-max(A[i][0],B[j][0]))
  if A[i][1]<B[j][1]:i+=1
  else:j+=1
 return t
def literal(L,c,e,p,f,q):return inter(arcs(L,e,p,c),arcs(L,f,q,c))

HOSTILES=(
 (2520,1953,8,134,9,32),
 (1848,1552,7,133,1,495),
 (336,158,2,38,12,461),
 (2520,55,3,3,10,352),
 (336,144,8,284,7,1),
 (168,112,6,5,1,385),
 (168,33,3,1,8,1),
)

def main():
 sem=sha256();checks=0
 for row in HOSTILES:
  x=M.mass(*row);y=literal(*row)
  if x!=y:raise RuntimeError(('hostile',row,x,y))
  sem.update((repr((row,x))+'\n').encode());checks+=1
 random.seed(0)
 rulers=(168,336,840,1848,2520)
 # Exhaust the equal-level orientation boundary where ordering by p is wrong.
 for L in rulers:
  for cell in range(L):
   for e in range(1,15):
    for f in range(1,15):
     if e==f:continue
     row=(L,cell,e,1,f,1);x=M.mass(*row);y=literal(*row)
     if x!=y:raise RuntimeError(('equal',row,x,y,x-y))
     sem.update((repr((row,x))+'\n').encode());checks+=1
 for index in range(5000):
  L=random.choice(rulers);cell=random.randrange(L);e,f=random.sample(range(1,15),2)
  p=random.randrange(1,1000);q=random.randrange(1,1000);row=(L,cell,e,p,f,q)
  x=M.mass(*row);y=literal(*row)
  if x!=y:raise RuntimeError(('random',index,row,x,y,x-y))
  sem.update((repr((row,x))+'\n').encode());checks+=1
 if checks!=EXPECTED_CHECKS or sem.hexdigest()!=EXPECTED_SEMANTIC:
  raise RuntimeError(('frozen audit',checks,sem.hexdigest()))
 print('LRC GENERAL REFLECTED PAIR MASS AUDIT')
 print('status=VERIFIED-EXACT hostile + seeded literal interval audit')
 print('hostiles',len(HOSTILES),'equal_level_checks',sum(L*14*13 for L in rulers),'random_seed',0,'random_checks',5000,'checks',checks)
 print('semantic_sha256',sem.hexdigest())
 print('engine_sha256',sha256(ENGINE.read_bytes().replace(b'\r\n',b'\n')).hexdigest())
 print('audit_sha256',sha256(open(__file__,'rb').read().replace(b'\r\n',b'\n')).hexdigest())

if __name__=='__main__':main()
