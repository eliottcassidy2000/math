#!/usr/bin/env python3
import json,hashlib
from pathlib import Path
from functools import lru_cache
from inherited import profile,tournaments,require
from blocks import graph
from transfer import S,mv,ell,v

def induced(rows,vs):return [sum(1<<j for j,w in enumerate(vs) if rows[u]>>w&1) for u in vs]
def singleton_word(bits):
 n=len(bits); a=[0]*(n+1)
 for i,b in enumerate(bits,1):
  a[i]=sum(1<<j for j in range(i+1,n+1))|(not b)
  if b:a[0]|=1<<i
 return a
@lru_cache(None)
def bank(U):
 if U==3:return {(3,3,3,3)}
 if U<3 or U%3:return set()
 bits=[0,0]+list(map(int,bin(U//3-1)[2:]));return {profile(graph(bits,p),0) for p in range(len(bits)-1)}
def strong(a):
 n=len(a)
 for r in range(n):
  seen=1<<r
  while True:
   new=seen
   for i in range(n):
    if seen>>i&1:new|=a[i]
   if new==seen:break
   seen=new
  if seen!=(1<<n)-1:return False
 return True
stats={1:0,3:0};total=0
for n in range(1,7):
 for a in tournaments(n):
  x=profile(a,0);total+=1
  if x[3]==1:
   M=x[2]-1
   b=[] if M==0 else list(map(int,bin(M)[2:]))
   require(x==profile(singleton_word(b),0),'deletion-one exact representative');stats[1]+=1
  if x[3]==3:
   require(x in bank(x[2]),'deletion-three canonical membership');stats[3]+=1
# All canonical bank candidates are physical by construction. Check terminal E formula explicitly.
canon=0
for M in range(1,129):
 bits=[0,0]+list(map(int,bin(M)[2:]));L=len(bits)
 for p in range(L-1):
  x=profile(graph(bits,p),0);k=2*bits[p]+bits[p+1]
  E=(3-k) if p==L-2 else (0 if bits[-1] else 3)
  require(x[1:]==(E,3*(1+M),3),'canonical sidecar and E');canon+=1
strong_controls=[]
for r in range(3,6):
 B=(4**r+6*r-1)//3;vals=[]
 for p in range(1,2*r-3):
  a=graph([0,1]*r,p); a=induced(a,[0]+list(range(2,len(a)-1)))
  x=profile(a,0)
  require(strong(a),'strongness');require(x==(B+(-2)**p-2,3,2**(2*r-1)+1,3),'strong profile')
  vals.append(x[0])
 strong_controls.append({'r':r,'scalars':vals})
files=['MEMBERSHIP.md','STRONG.md','verify_membership.py']
print(json.dumps({'status':'PASS','labelled_tournaments_root0_orders1_to6':total,'membership_controls':stats,'canonical_candidates_M1_to128':canon,'strong_controls':strong_controls,'sha256':{f:hashlib.sha256(Path(__file__).with_name(f).read_bytes()).hexdigest() for f in files}},sort_keys=True,indent=2))
