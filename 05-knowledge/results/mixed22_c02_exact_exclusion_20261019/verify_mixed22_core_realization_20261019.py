#!/usr/bin/env python3
"""Exact core-cone obstruction and actual exterior realization, stdlib only.
This proves a limitation of core-only screens, not a counterexample to balance.
"""
from fractions import Fraction
from functools import lru_cache
from itertools import combinations,product
from pathlib import Path
import hashlib,json
ROOT=Path(__file__).resolve().parent
CORE=(0,12,0,4,1,21)
ORDER=(2,0,3,4,1,5)
PAIRS=((0,1),(2,0),(0,3),(4,1),(1,5),(2,4),(3,4),(3,5))

def check(ok,msg):
 if not ok:raise RuntimeError(msg)

def closure(n,relations):
 r=set(relations)
 while True:
  nxt=r|{(a,d) for a,b in r for c,d in r if b==c}
  if nxt==r:break
  r=nxt
 check(all(a!=b for a,b in r),'cycle')
 return tuple(sum(1<<a for a,b in r if b==v) for v in range(n))

def extensions(p):
 out=[];n=len(p)
 def visit(w,mask):
  if len(w)==n:out.append(w);return
  for v,d in enumerate(p):
   if not(mask>>v&1) and not(d&~mask):visit(w+(v,),mask|1<<v)
 visit((),0);return out

def count(p):
 full=(1<<len(p))-1
 @lru_cache(None)
 def visit(mask):
  if mask==full:return 1
  return sum(visit(mask|1<<v) for v,d in enumerate(p) if not(mask>>v&1) and not(d&~mask))
 return visit(0)

def realize():
 q=(0,0,1,1,7);n=16
 rel={(i,j) for j,m in enumerate(CORE) for i in range(6) if m>>i&1}
 rel|={(6+i,6+j) for j,m in enumerate(q) for i in range(5) if m>>i&1}
 rel|={(6+i,j) for i in range(5) for j in (0,1,3,4,5)}
 rel|={(6+i,2) for i in (0,3)}
 dual={0:1,1:0,2:5,5:2,3:4,4:3,**{6+i:11+i for i in range(5)},**{11+i:6+i for i in range(5)}}
 rel|={(dual[b],dual[a]) for a,b in list(rel)}
 return closure(n,rel),q

def main():
 w=extensions(CORE);check(len(w)==19,'core extensions')
 columns=[]
 for i,j in product(range(3),repeat=2):
  prefix=(2,3)[:i];suffix=(4,5)[2-j:]
  selected=[s for s in w if s[:i]==prefix and (not j or s[-j:]==suffix)]
  numer=[sum(s.index(a)<s.index(b) for s in selected) for a,b in PAIRS]
  columns.append({'prefix':prefix,'suffix':suffix,'E':len(selected),'numerators':numer,'deficits':[3*v-2*len(selected) for v in numer]})
 by={(tuple(g['prefix']),tuple(g['suffix'])):g for g in columns}
 weights={((),()):1,((2,),()):5,((),(5,)):5}
 e=sum(weight*by[key]['E'] for key,weight in weights.items());nn=[sum(weight*by[key]['numerators'][i] for key,weight in weights.items()) for i in range(8)]
 check(e==119 and nn==[108,90,80,80,90,111,80,111],'three-generator obstruction')
 u=v=Fraction(15,11);weights={((),()):1,((2,),()):u,((),(5,)):v,((2,),(5,)):u*v}
 pe=sum(weight*by[key]['E'] for key,weight in weights.items());pn=[sum(weight*by[key]['numerators'][i] for key,weight in weights.items()) for i in range(8)]
 check(pe==Fraction(6949,121) and min(3*n-2*pe for n in pn)==Fraction(7,121),'rank-one obstruction')
 p,q=realize();check(count(q)==11 and count(q+(9,))==26,'left boundary insertion counts')
 words=extensions(p);total=len(words);check(total==count(p)==6949,'actual extension count')
 allpairs=list(combinations(range(16),2));counts={pair:0 for pair in allpairs};fibers={s:0 for s in w}
 for word in words:
  rank={a:i for i,a in enumerate(word)}
  for a,b in allpairs:counts[a,b]+=rank[a]<rank[b]
  fibers[tuple(a for a in word if a<6)]+=1
 for word,number in fibers.items():
  expected=(11+15*(word[0]==2))*(11+15*(word[-1]==5))
  check(number==expected,'fixed-core product fiber identity')
 rel={(a,b) for b,d in enumerate(p) for a in range(16) if d>>a&1}
 for (a,b),number in counts.items():
  if (b,a) in rel:expected=0
  elif (a,b) in rel:expected=total
  else:expected=count(closure(16,rel|{(a,b)}))
  check(number==expected,'independent subset-DP pair count')
 inc=[(a,b) for a,b in allpairs if not(p[b]>>a&1 or p[a]>>b&1)]
 degree=[sum(a==v or b==v for a,b in inc) for v in range(16)]
 corecounts=[counts[min(a,b),max(a,b)] if a<b else total-counts[b,a] for a,b in PAIRS]
 check(corecounts==[6273,5200,4635,4635,5200,6421,4645,6421],'actual core numerators')
 check(min(3*n-2*total for n in corecounts)==7 and max(degree)==5,'actual core and degree claims')
 anchored_lower=[total-counts[2,a] for a in range(6,11)]
 anchored_upper=[counts[5,b] for b in range(11,16)]
 check(anchored_lower==anchored_upper==[6949,6149,5949,6949,3949],'anchored contact count')
 balanced=[{'pair':[a,b],'numerator':counts[a,b],'denominator':total} for a,b in inc if total<=3*counts[a,b]<=2*total]
 check(len(balanced)==8,'balanced pair count')
 result={'status':'ALL_EXACT_CONTROLS_PASSED','scope':'A core-only method obstruction with actual five-thin realization containing balanced exterior pairs; no conjecture counterexample.',
 'core_pred':CORE,'core_majority_order':ORDER,'forward_core_pairs':PAIRS,'core_extensions':19,'all_nine_generator_columns':columns,
 'unrestricted_three_generator':{'coefficients':[1,5,5],'E':e,'numerators':nn,'minimum_probability':'80/119'},
 'rank_one':{'u':str(u),'v':str(v),'denominator':str(pe),'minimum_probability':'4635/6949','minimum_3p_minus_2':'7/6949'},
 'actual_poset':{'pred':p,'labels':['x','y','b1','b2','a1','a2']+[f'A{i}' for i in range(5)]+[f'B{i}' for i in range(5)],'left_pred':q,'left_b1_predecessor_mask':9,'e_left':11,'e_left_with_b1':26,'extensions':total,'incomparability_degrees':degree,'core_numerators':corecounts,'anchored_A_before_b1_numerators':anchored_lower,'anchored_a2_before_B_numerators':anchored_upper,'balanced_pairs':balanced},
 'controls':{'fixed_core_fibers':len(fibers),'independent_DP_pair_counts':len(allpairs),'generator_columns':len(columns)},
 'verifier_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
 target=ROOT/'mixed22_core_realization_receipt_20261019.json';target.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({'status':result['status'],'extensions':total,'balanced_pairs':len(balanced),'receipt':str(target)}))
if __name__=='__main__':main()
