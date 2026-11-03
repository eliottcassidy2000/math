#!/usr/bin/env python3
"""Exhaustive three-vertex controls for operator-class restrictions."""
import itertools as it,json
from fractions import Fraction
checks=0
def need(ok,msg):
 global checks
 checks+=1
 if not ok:raise RuntimeError(msg)
def graph(n,a):
 return {(i,j):bool(a>>k&1) for k,(i,j) in enumerate(it.combinations(range(n),2))}
def h(n,a,cut=None):
 g=graph(n,a)
 if cut is not None:
  g.update({(i,n):not bool(cut>>i&1) for i in range(n)});n+=1
 return sum(all(g[min(u,v),max(u,v)]==(u<v) for u,v in zip(w,w[1:])) for w in it.permutations(range(n)))
H=[h(3,a) for a in range(8)];B=[h(3,a,1) for a in range(8)]
means=[]
for parity in (0,1):
 aa=[a for a in range(8) if a.bit_count()%2==parity]
 means.append({'parity':parity,'H':str(Fraction(sum(H[a] for a in aa),4)),'B0':str(Fraction(sum(B[a] for a in aa),4)),'condition_a01_1_H':str(Fraction(sum(H[a] for a in aa if a&1),2)),'adaptive_H':str(Fraction(sum((B if parity==0 else H)[a] for a in aa),4))})
need(means[0]['H']==means[1]['H']=='3/2','unconditioned equal H')
need(means[0]['B0']==means[1]['B0']=='7/2','fixed B0 equal H')
need(means[0]['condition_a01_1_H']=='2' and means[1]['condition_a01_1_H']=='1','conditioning is outside carrier')
need(means[0]['adaptive_H']=='7/2' and means[1]['adaptive_H']=='3/2','adaptive operation outside carrier')
for cut in range(8):
 vals=[h(3,a,cut) for a in range(8)]
 need(sum((1 if a.bit_count()%2==0 else -1)*vals[a] for a in range(8))==0,'every fixed one-birth cut parity invisible')
print(json.dumps({'checks':checks,'universe':'all 8 labelled tournaments on [3], lexicographic edge bits','H':H,'B0':B,'means':means,'adaptive_rule':'append B0 on even old-edge parity, a universal sink on odd parity','conditional_event':'a01=1; probability1/2 under both parity laws'},indent=2))
