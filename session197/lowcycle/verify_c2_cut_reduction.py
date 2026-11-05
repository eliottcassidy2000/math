#!/usr/bin/env python3
"""Finite controls for the partial C2 reduction; no census extrapolation."""
from pathlib import Path
import hashlib,json
from verify_c1_complete import graphs,cycles,strong,induced,components,path_counts,require
P=Path(__file__).resolve().parent
require(hashlib.sha256((P/'verify_c1_complete.py').read_bytes()).hexdigest()=='645ddfc551fd5ce2760f2d1f71318c08f57e1115e54edf9d28d7db3e6ad91d3b','dependency hash')

def bridge_list(r):
 ans=[]
 for u in range(len(r)):
  for v in range(len(r)):
   if r[u]>>v&1:
    s=r[:];s[u]^=1<<v
    if not strong(s):ans.append((u,v))
 return ans

def main():
 total=0;stats=[];cutchecks=0
 for n in range(1,7):
  two=0;minimum=None;hist={}
  for r in graphs(n):
   total+=1
   if not strong(r) or n<3:continue
   cuts=[]
   for v in range(n):
    keep=[u for u in range(n) if u!=v]
    if not strong(induced(r,keep)):cuts.append((v,keep))
   cs=cycles(r);c=len(cs)
   if not cuts:
    two+=1;minimum=c if minimum is None else min(minimum,c);hist[c]=hist.get(c,0)+1
   elif c<=3:
    root,keep=cuts[0];Ds=[tuple(keep[u] for u in D) for D in components(induced(r,keep))]
    _,s,_=path_counts(r,Ds[0]);_,_,e=path_counts(r,Ds[-1])
    aa=[u for u in Ds[0] if r[root]>>u&1];bb=[u for u in Ds[-1] if r[u]>>root&1]
    A=sum(s[u] for u in aa);B=sum(e[u] for u in bb)
    require(A==1 or B==1,'low-cycle endpoint factors')
    u,v=(root,aa[0]) if A==1 else (bb[0],root)
    z=r[:];z[u]^=1<<v;require(not strong(z),'cut endpoint edge not bridge');cutchecks+=1
  stats.append(dict(order=n,two_strong_labelled=two,minimum_cycles=minimum,cycle_histogram=hist))
 require(stats[4]['two_strong_labelled']==24 and stats[4]['minimum_cycles']==2,'regular5 bank')
 require(stats[5]['two_strong_labelled']==2640 and stats[5]['minimum_cycles']==4,'order6 bank')
 require(5**41<=2**96<5**42 and 5**65<=3**96<5**66,'Busch integer cutoffs')
 witness=[10,12,9,112,39,71,23]
 require(strong(witness) and len(cycles(witness))==9 and not bridge_list(witness),'C9 control')
 cuts=[v for v in range(7) if not strong(induced(witness,[u for u in range(7) if u!=v]))]
 require(cuts,'C9 lacks separating vertex')
 hostile=[14,20,56,50,33,3];hcs=cycles(hostile)
 require(all(strong(induced(hostile,[u for u in range(6) if u!=v])) for v in range(6)),'arc hostile not2strong')
 require(len(hcs)==5 and hostile[2]>>5&1,'arc hostile counts')
 require(all(not any(u==2 and v==5 for u,v in zip(w,w[1:]+w[:1])) for w in hcs),'arc hostile is Hamiltonian')
 result=dict(status='PASS',labelled_tournaments=total,maximum_labelled_order=6,low_cycle_cut_bridge_checks=cutchecks,two_strong_inventory=stats,bridgefree_C2_order_bound=73,bridgefree_C3_order_bound=97,C9_bridgefree_cut_example=dict(rows=witness,separating_vertices=cuts),non_Hamiltonian_arc_control=dict(rows=hostile,C=5,arc=[2,5]))
 data=json.dumps(result,sort_keys=True,indent=2)+'\n';(P/'c2_cut_reduction_certificate.json').write_text(data);print(data,end='')
if __name__=='__main__':main()
