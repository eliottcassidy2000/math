"""Exact mixed-drift cancellation and paired periodic-return seed equations."""
from find_intersections import g,rev,append_coeff,HERE
from collections import defaultdict
import json

def drift_words(L):
 Q=(1<<L)-1;ds=defaultdict(list)
 for r in range(1,Q):
  a,b,c=append_coeff(r,L);N=(b*(r&1)+c)*Q-a*r;ds[N].append(r)
 return ds

def main():
 stats=[];pairs=[]
 for L in range(2,15):
  ds=drift_words(L);m=[]
  for D,ws in ds.items():
   if D%13==0 and 25*(D//13) in ds:
    for rb in ws:
     for rc in ds[25*(D//13)]:m.append((rb,rc))
  stats.append({'length':L,'ordered_word_pairs':len(m)})
  if L==14:pairs=m
 lookup=defaultdict(list)
 for D in range(1,32768):lookup[25*g(D)+244*D+186].append(D)
 seeds=[]
 for C in range(1,32768):
  H=13*g(C)+334*C+198
  for D in lookup.get(H,[]):seeds.append((C,D,H))
 Q=(1<<14)-1;hits=[]
 for pb,pc in [(0,0)]+[(D,C) for C,D,H in seeds]:
  for rb,rc in pairs:
   ab,bb,cb=append_coeff(rb,14);ac,bc,cc=append_coeff(rc,14)
   if bb*((pb&1)-(rb&1)) or bc*((pc&1)-(rc&1)):continue
   Kb=25*ab+244*Q;Kc=13*ac+334*Q
   if Kb*(Q*pb+rb)==Kc*(Q*pc+rc):hits.append({'new_prefix':pb,'old_prefix':pc,'new_loop':rb,'old_loop':rc})
 result={'word_search':stats,'chosen_pair':{'new_word':'00000000000101','old_word':'00001001110101','new_drift_numerator':18434,'old_drift_numerator':35450,'denominator':Q,'weighted_drift_numerator':460850},'synchronized_seed_box':32767,'synchronized_seeds':len(seeds),'exact_period14_returns':hits}
 (HERE/'period14.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
if __name__=='__main__':main()
