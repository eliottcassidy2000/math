#!/usr/bin/env python3
"""Exact envelope vertices and arbitrary-core recursive tail transplant controls."""
from fractions import Fraction as F
from pathlib import Path
from math import gcd,lcm,ceil
from itertools import combinations,product
import importlib.util,json
ROOT=Path(__file__).resolve().parents[3]
p=ROOT/'05-knowledge/results/entry_cosets_resonant_echoes_20261030/lrc_audit.py'
spec=importlib.util.spec_from_file_location('windows',p);old=importlib.util.module_from_spec(spec);spec.loader.exec_module(old)
def need(b,msg):
 if not b:raise RuntimeError(str(msg))
def dist(x):return min(x%1,(-x)%1)
def envelope(V):
 V=tuple(sorted(set(V)));events={F(0),F(1)}
 denoms={2*v for v in V}|{a+b for a,b in combinations(V,2)}|{abs(a-b) for a,b in combinations(V,2)}
 for d in denoms:events.update(F(k,d) for k in range(d+1))
 T=sorted(events);Y=[min(dist(v*t) for v in V) for t in T]
 slopes=[(b-a)/(y-x) for x,y,a,b in zip(T,T[1:],Y,Y[1:])]
 corners=[i for i in range(1,len(T)-1) if slopes[i]!=slopes[i-1]]
 mu=max(Y);need(mu<=F(1,2),'distance radius')
 maxima=[i for i in corners if Y[i]==mu]
 low=max([F(0)]+[Y[i] for i in corners if Y[i]<mu]);Q=lcm(*(T[i].denominator for i in maxima));p=Q*mu
 need(p.denominator==1,'integral radius in common frame');p=int(p)
 sides=[(T[i],-1,int(slopes[i-1])) for i in maxima]+[(T[i],1,-int(slopes[i])) for i in maxima]
 need(all(s>0 for _,_,s in sides),'positive maximum owners');s=min(x[2] for x in sides)
 threshold=max(F(s,Q-1),low*s/((Q-1)*(mu-low)));m=threshold.numerator//threshold.denominator+1
 kmin=1 if mu<F(1,2) else ceil(F(Q*m+s,Q*s))
 return dict(V=V,mu=mu,low=low,Q=Q,p=p,s=s,m=m,kmin=kmin,sides=sides,event_count=len(T))
def check_core(core,Ks=(1,2,1),do_global=True):
 c=envelope(core);Q,p,s,m=(c[k] for k in ('Q','p','s','m'));V=list(c['V']);aa=m;bb=0;P=1;H=F(1)
 limit=F(p*m*(Q-1),Q*(m*(Q-1)+s));need(limit>c['low'],'full proof stays in chamber');records=[]
 for j in range(len(Ks)+1):
  V.append(Q*aa);rho=F(p*m,Q*m+s*H);x=p*H/(Q*m+s*H)
  pts=sorted({centre+sgn*x/Q for centre,sgn,owner in c['sides'] if owner==s})
  need(F(p*aa-s*bb,Q*aa+s)==rho,'native coefficient identity')
  need(all(min(dist(v*t) for v in V)==rho for t in pts),'actual all-speed witness')
  if do_global:
   need(old.safe(V,rho)==[(t,t) for t in pts],('complete max set',core,j))
   need(not old.safe(V,rho+F(1,8*max(V)**2)),('no larger maximum',core,j))
  records.append(dict(j=j,speeds=list(V),maximum=str(rho),points=list(map(str,pts))))
  if j<len(Ks):
   K=c['kmin']*Ks[j];aa,bb=K*(Q*aa+s),K*(Q*bb+p);P*=Q*K;H+=F(1,P)
 return c,records
def main():
 cases=[];total=0;native=0
 # Explicit finite universe: every primitive triple through12, including half-radius cores.
 for V in combinations(range(1,13),3):
  if gcd(*V)!=1:continue
  c,rows=check_core(V,(1,2),True);total+=1;native+=len(rows)
  if V in ((1,2,3),(1,3,5),(1,6,7),(2,3,7),(3,4,11)):cases.append(dict(core=V,Q=c['Q'],p=c['p'],s=c['s'],m=c['m'],kmin=c['kmin'],mu=str(c['mu']),chamber_floor=str(c['low']),rows=rows))
 c,rows=check_core((1,6),(1,2,7),True)
 need((c['mu'],c['low'],c['Q'],c['p'],c['s'],c['m'])==(F(3,7),F(2,5),7,3,1,3),'repaired local-germ hostile')
 reroute,rr=check_core((1,6,21,23),(2,1),True)
 need(reroute['mu']==F(9,22) and reroute['s']==21,'owner change reverses cheap side')
 # Independent pair-clock maxima on small first additions, not the same interval engine.
 pair_checks=0
 for V in ((1,2),(1,6),(2,3),(1,2,4),(2,3,7)):
  cc,rs=check_core(V,(),True);vv=rs[0]['speeds'];got,pts,_=old.pair.maximum(tuple(vv))
  need(str(got)==rs[0]['maximum'] and list(map(str,pts))==rs[0]['points'],'independent pair-clock maximum');pair_checks+=1
 half,rh=check_core((1,),(1,2,1),True)
 got,pts,_=old.pair.maximum((1,6,14,30))
 need(got==F(1,3)<F(12,31),'unguarded half-core folding hostile')
 out=dict(status='finite exact controls for all-core transplant theorem',primitive_triple_cores=total,complete_recursive_banks=native,
          pair_clock_checks=pair_checks,half_core_guard=half['kmin'],half_fold_hostile=dict(speeds=(1,6,14,30),predicted='12/31',actual=str(got)),examples=cases,germ_hostile=dict(core=(1,6),rows=rows),
          changed_owner=dict(Q=reroute['Q'],p=reroute['p'],s=reroute['s'],m=reroute['m'],floor=str(reroute['low']),rows=rr))
 Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n')
 print(json.dumps({k:v for k,v in out.items() if k not in ('examples','germ_hostile','changed_owner')},indent=2));print('reroute', {k:out['changed_owner'][k] for k in ('Q','p','s','m','floor')})
if __name__=='__main__':main()
