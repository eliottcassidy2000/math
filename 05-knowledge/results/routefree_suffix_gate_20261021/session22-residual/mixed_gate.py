#!/usr/bin/env python3
"""Separate exact mixed-zonotope corollary; base receipts are unchanged."""
import argparse,hashlib,itertools,json
from pathlib import Path
import check as base
CHECKS=0
def ck(ok,msg):
 global CHECKS
 CHECKS+=1
 if not ok:raise ValueError(msg)
def vertices(count):
 t=tuple(sum(count[i]*base.VECTORS[i][k]for i in range(12))for k in (0,1))
 return sorted(set((2*z[0]-t[0],2*z[1]-t[1])for z in base.cuts(count)))
def mixed(old,new):
 return max(base.Q(n)+base.beta(o,n)for o in vertices(old)for n in vertices(new))
def counter(ids):return [ids.count(i)for i in range(12)]
def brute_vectors(ids):
 return [tuple(sum(s*base.VECTORS[i][k]for s,i in zip(signs,ids))for k in (0,1))for signs in itertools.product((-1,1),repeat=len(ids))]
def run():
 sets=[ids for n in range(3)for ids in itertools.combinations_with_replacement(range(12),n)]
 tested=0;strict=0;largest=0;witness=None
 for oi in sets:
  oc=counter(oi);ov=brute_vectors(oi);Do=base.zonotope(oc)
  for ni in sets:
   nc=counter(ni);nv=brute_vectors(ni);Dn=base.zonotope(nc);G=mixed(oc,nc)
   brute=max(base.Q(n)+base.beta(o,n)for o in ov for n in nv)
   ck(G==brute,'144-vertex exact maximum versus all independent signs')
   # Square form audits the stronger real-valued inequality, not just ceil.
   ck(G>=Dn and (G-Dn)**2<=4*Do*Dn,'Cauchy domination without floating point')
   old=Dn+base.ceilroot(4*Do*Dn)
   if G<old:strict+=1
   if old-G>largest:largest=old-G;witness=[list(oi),list(ni),G,old]
   if not oi:ck(G==Dn,'empty-old equality')
   if not ni:ck(G==0,'empty-new equality')
   tested+=1
 W={0:(1,-1),1:(-2,-1),2:(0,0)};H={0:(2,1),1:(0,-1),2:(0,0)}
 inputs=[(0,(1,0)),(0,(1,2)),(1,(-1,1)),(1,(0,1))]
 counts=[[0]*12,[0]*12]
 for k,v in inputs:counts[k][base.INDEX[v]]+=1
 G=mixed(*counts);bg=base.gate(W,H,*counts)
 bound=bg['delta']+G+bg['L_old']+bg['L_new'];maxactual=-10**9;maxquadratic=-10**9;assignments=0
 for choices in itertools.product(range(7),repeat=len(inputs)):
  O={u:(0,0)for u in W};N=O.copy()
  for (k,v),z in zip(inputs,choices):
   if not z:continue
   u=(z-1)//2;sg=1 if z%2 else -1;target=N if k else O
   target[u]=base.add(target[u],tuple(sg*x for x in v))
  quadratic=sum(base.Q(N[u])+base.beta(O[u],N[u])for u in W)
  delta=sum(base.Q(base.add(base.add(W[u],H[u]),base.add(O[u],N[u])))-base.Q(base.add(W[u],O[u]))for u in W)
  ck(quadratic<=G,'common-terminal sign grouping with absent sources')
  ck(delta<=bound,'strengthened all-future energy bound')
  maxactual=max(maxactual,delta);maxquadratic=max(maxquadratic,quadratic);assignments+=1
 a=counter((0,));b=counter((4,))
 ck(mixed(a,b)==2 and base.zonotope(b)+base.ceilroot(4*base.zonotope(a)*base.zonotope(b))==3,'strict one-old one-new angular witness')
 return {'status':'PROVED exact mixed-zonotope strengthening; base native receipts unchanged','multiset_pairs':tested,'strict_improvements_over_integer_Cauchy':strict,'largest_integer_improvement':largest,'largest_improvement_witness':witness,'destination_assignments':assignments,'mixed_test':{'G':G,'strengthened_upper':bound,'base_upper':bg['upper'],'actual_max':maxactual,'quadratic_actual_max':maxquadratic},'strict_simple_witness':{'old':[1,0],'new':[0,1],'G':2,'base_quadratic_upper':3},'checks':CHECKS}
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',type=Path);p.add_argument('--verify',type=Path);a=p.parse_args();result=run();blob=(json.dumps(result,sort_keys=True,separators=(',',':'))+'\n').encode()
 if a.output:a.output.write_bytes(blob)
 if a.verify:ck(a.verify.read_bytes()==blob,'immutable receipt replay')
 print(json.dumps({'receipt':result,'sha256':hashlib.sha256(blob).hexdigest(),'checks':CHECKS},sort_keys=True))
