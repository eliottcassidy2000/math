#!/usr/bin/env python3
"""Restricted all-order theorem controls and literal rerooting escape."""
import importlib.util,json
from pathlib import Path
from itertools import product
from math import comb
s=importlib.util.spec_from_file_location('v',Path(__file__).with_name('verify.py'));v=importlib.util.module_from_spec(s);s.loader.exec_module(v)
D=(18,8,3,5,14);Ts=[v.attach(D,c) for c in (9,3,24)];Ps=list(map(v.packet,Ts))
def strong(T):
 for i in range(len(T)):
  seen=1<<i
  while True:
   nxt=seen
   for j in range(len(T)):
    if seen>>j&1:nxt|=T[j]
   if nxt==seen:break
   seen=nxt
  if seen!=(1<<len(T))-1:return False
 return True
v.check(strong(D),'deletion strong')
def tiled(W):
 j=i=0
 while i<len(W):
  if W[i]==0:i+=1
  elif W[i:i+2]==(1,2):i+=2;j+=1
  else:return None
 return j
checks=0;rows=[]
for k in range(1,9):
 base=v.word([Ps[0]]*k);hist={}
 for W in product(range(3),repeat=k):
  delta=sum([0,-8,1][letter]*8**i for i,letter in enumerate(W));j=tiled(W)
  v.check((delta==0)==(j is not None),'exact receiver monomer-dimer language')
  if j is not None:
   P=v.word([Ps[t] for t in W]);v.check(all(P[z]==base[z] for z in P if z!='H'),'full shadow equality')
   v.check(P['H']==base['H'] if k==1 else P['H']==base['H']+2*j*13**(k-2),'exact count fiber')
   hist[j]=hist.get(j,0)+1
  checks+=1
 v.check(hist=={j:comb(k-j,j) for j in range(k//2+1)},'binomial multiplicity');rows.append(dict(k=k,multiplicities=hist))
def delete(T,r):
 ids=[i for i in range(len(T)) if i!=r];return tuple(sum(1<<j for j,b in enumerate(ids) if T[a]>>b&1) for a in ids)
menus=[]
for W in ((1,2,0),(0,1,2)):
 T=v.fuse(v.fuse(Ts[W[0]],Ts[W[1]]),Ts[W[2]]);v.check(strong(T),'physical strong word')
 H,end=v.hp_subsets(T);F=(1<<len(T))-1;menu=[H[F^(1<<i)] for i in range(len(T))]
 # Independent source-packet assembly of every nonroot deletion.
 algebra=[13**3]
 for i in range(3):
  for j in range(1,6):
   source=[Ps[t] for t in W];source[i]=v.packet(delete(Ts[W[i]],j));algebra.append(v.word(source)['H'])
 v.check(menu==algebra,'literal and component deletion menus')
 Vs=[sum(end(S,r)*H[F^S] for S in range(1,F+1) if S>>r&1) for r in range(len(T))]
 v.check(H[-1]==264251,'same scalar count')
 menus.append(dict(word=''.join('ABC'[i] for i in W),rows=list(T),deletion_menu=menu,ending_mass=Vs))
v.check(sorted(menus[0]['deletion_menu'])!=sorted(menus[1]['deletion_menu']),'unrooted nonisomorphism witness')
new=264251+2*menus[1]['ending_mass'][1];v.check(new==914613 and new%13==11,'physical reroot triangle escape')
print(json.dumps(dict(status='PASS',all_word_checks=checks,fiber_levels=rows,root_menus=menus,reroot_triangle_count=new,scope='Exact restricted monomer-dimer fiber; no universal target coverage'),sort_keys=True,indent=2))
