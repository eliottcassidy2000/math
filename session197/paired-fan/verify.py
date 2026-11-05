#!/usr/bin/env python3
from itertools import combinations
from pathlib import Path
import json
from core import *
ROOT=Path(__file__).resolve().parent
C=[26,16,3,6,12];a=2;A=[90,16,3,6,12,31,62];B=[26,48,3,6,44,13]
require(plus2(C)[0]==A and plus4(C)[0]==B,'literal core gates')
expected=[(13,2,18,3,3,15),(41,6,54,9,3,99),(43,8,54,13,6,72)]
for T,p in zip([C,A,B],expected):require(profile(T,a)==literal_profile(T,a)==p,'independent core profiles');require(strong(T),'strong core')
S=[8,1,3,6];require(profile(S,1)==literal_profile(S,1)==(5,2,6,3,1,9),'unique-start source profile')
fan_results=[]
for m in range(1,7):
 P,b=fan(C,a,m);X,portsA=plus2(P);Y,portsB=plus4(P);PX,bX=fan(A,a,m);PY,bY=fan(B,a,m);N=len(P)
 mapX=list(range(5))+[N,N+1,5]+list(range(6,6+m))
 mapY=list(range(5))+[N,5]+list(range(6,6+m))
 require(relabel(PX,mapX)==X and relabel(PY,mapY)==Y,'graph-level gate/chain commutation')
 H=18*2**m+13
 require(profile(P,b)[:4]==(H,13,H,13),'equal-response parent')
 require(profile(X,b)[:4]==(3*H+2,41,3*H+2,41),'plus2 fan profile')
 require(profile(Y,b)[:4]==(3*H+4,43,3*H+4,43),'plus4 fan profile')
 require(strong(P) and strong(X) and strong(Y),'strong fan')
 require(portsA==[0,N] and portsB==[1,4],'fixed physical exposure')
 fan_results.append({'m':m,'parent_H':H,'plus2_H':3*H+2,'plus4_H':3*H+4,'parent_rows':P,'plus2_rows':X,'plus4_rows':Y,'parent_root':b,'plus2_ports':portsA,'plus4_ports':portsB})
# Universal prefix-triangle gate, including singleton endpoint correction.
checks=0
for n in range(1,5):
 edges=list(combinations(range(n),2))
 for bits in range(1<<len(edges)):
  T=[0]*n
  for k,(u,v) in enumerate(edges):
   if bits>>k&1:T[u]|=1<<v
   else:T[v]|=1<<u
  for r in range(n):
   x=n;T0=T[:]+[sum(1<<v for v in range(n) if v!=r)];T0[r]|=1<<x
   # Source u, then v that loses only to x; no port-index ambiguity.
   T1,ps=left(T0,[x,r],0,1);T1,ps=left(T1,ps,1,2)
   p=profile(T0,r);q=profile(T1,r);s=profile(T,r)[4]
   require(q[:4]==(3*p[0]+2*s,3*p[1]+2*(n==1),3*p[2],3*p[3]),'universal triangle gate')
   G,b=fan(T0,r,1);GG,bb=fan(T1,r,1)
   require(strong(G) and strong(GG),'general gate fan strongness')
   require(dp(GG)[0][-1]==3*dp(G)[0][-1]+2*s,'general fan digit gate')
   checks+=1
# Exact small prefix factory exclusion: no nonnegative endpoint count can give H43 at R54.
solutions=[(h,R) for h in range(1,10,2) for R in range(h,10) if h+2*R==9]
require(solutions==[(1,4),(3,3)],'prefix factory equation')
require(all(5*h+9*R<43 for h,R in solutions),'prefix factory H43 ceiling')
print(json.dumps({'core_profiles':[list(v) for v in expected],'source_profile':list(profile(S,1)),'positive_chain_controls':fan_results,'universal_triangle_gate_controls':checks,'prefix_factory_candidates':solutions},sort_keys=True,indent=2))
