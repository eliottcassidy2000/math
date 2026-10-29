#!/usr/bin/env python3
"""Independent Session94 exact checks; stdlib, no producer imports or asserts."""
from itertools import product
from fractions import Fraction
from math import comb, factorial
from collections import Counter
import json, random
checks=0
def check(x,why):
 global checks
 checks+=1
 if not x: raise RuntimeError(why)
def mul(p,q): return tuple(p[q[i]] for i in range(len(p)))
def inv(p):
 out=[0]*len(p)
 for i,j in enumerate(p):out[j]=i
 return tuple(out)
def power(p,n):
 out=tuple(range(len(p)))
 while n:
  if n&1:out=mul(p,out)
  p=mul(p,p);n//=2
 return out
def sqrt(p):
 out=list(range(len(p))); seen=set()
 for a in range(len(p)):
  if a in seen:continue
  cyc=[];b=a
  while b not in seen:seen.add(b);cyc.append(b);b=p[b]
  check(len(cyc)%2==1,'odd-root domain')
  for i,b in enumerate(cyc):out[b]=cyc[(i+(len(cyc)+1)//2)%len(cyc)]
 return tuple(out)
def neg(w):return tuple((-x)%3 for x in w)
def local(words,w):
 lookup={v:i for i,v in enumerate(words)};d=len(w)
 return tuple(lookup[v[:d]+((v[d]+1)%3,)+v[d+1:]] if v[:d]==w else i for i,v in enumerate(words))
def centralizer(k,coeff):
 words=list(product(range(3),repeat=k));lookup={v:i for i,v in enumerate(words)}
 return tuple(lookup[tuple((v[d]+coeff.get(v[:d],0))%3 for d in range(k))] for v in words)
def portrait_from_perm(k,p):
 words=list(product(range(3),repeat=k));out={}
 for d in range(k):
  for w in product(range(3),repeat=d):
   i=words.index(w+(0,)*(k-d));out[w]=words[p[i]][d]
 return out
def loop_checks():
 rows=[];rng=random.Random(94)
 for k in range(1,5):
  words=list(product(range(3),repeat=k));lookup={v:i for i,v in enumerate(words)};alpha=tuple(lookup[neg(v)] for v in words);idp=tuple(range(3**k))
  reps=[w for d in range(k) for w in product(range(3),repeat=d) if any(w) and next(x for x in w if x)!=2]
  check(len(reps)==(3**k-2*k-1)//4,'portrait dimension')
  bs={w:mul(local(words,w),inv(local(words,neg(w)))) for w in reps}
  for w,b in bs.items():check(mul(alpha,b)==mul(b,alpha),'paired centralizer')
  # Every canonical node-pair generator is a literal based triangle.
  for d in range(1,k):
   for a in range(d):
    w=(0,)*a+(1,)+(0,)*(d-a-1);x=local(words,(0,)*a);y=local(words,(0,)*d)
    t=mul(power(x,2),alpha);u=mul(power(y,2),alpha)
    h=mul(inv(sqrt(mul(u,alpha))),mul(sqrt(mul(u,t)),sqrt(mul(t,alpha))))
    check(h in (bs[w],inv(bs[w])),'canonical triangle')
  # Explicit centralizer conjugator, no rotations at or below target depth.
  for w in reps:
   d=len(w);a=next(i for i,v in enumerate(w) if v);canonical=(0,)*a+(1,)+(0,)*(d-a-1)
   coeff={}
   for j in range(a+1,d):
    pref=canonical[:j];coeff[pref]=w[j];coeff[neg(pref)]=(-w[j])%3
   h=centralizer(k,coeff)
   check(mul(alpha,h)==mul(h,alpha),'conjugator in centralizer')
   check(mul(h,mul(bs[canonical],inv(h)))==bs[w],'literal paired generator conjugacy')
  choices=list(product(range(3),repeat=len(reps))) if k<=3 else [tuple(rng.randrange(3) for _ in reps) for _ in range(200)]
  seen=set();cost=Counter()
  for vals in choices:
   coeff={}
   for w,v in zip(reps,vals):coeff[w]=v;coeff[neg(w)]=(-v)%3
   p=centralizer(k,coeff);check(mul(alpha,p)==mul(p,alpha),'portrait centralizer')
   # Descending-depth factors execute first; source-prefix rotations are applied last.
   out=idp
   for w in sorted(reps,key=lambda x:(-len(x),x)):
    out=mul(power(bs[w],coeff[w]),out)
   check(out==p,'descending normal form')
   seen.add(p);cost[sum(v!=0 for v in vals)]+=1
  check(len(seen)==len(choices),'faithful portrait')
  if k<=3:
   for j in range(len(reps)+1):check(cost[j]==comb(len(reps),j)*2**j,'binomial cost')
  rows.append({'k':k,'free_pairs':len(reps),'tested':len(choices),'max_steps':3*len(reps)})
 return rows

def moment_checks():
 # Enumerate all automorphisms A_2 (81) and C_3 (243), not recurrence outputs.
 words=list(product(range(3),repeat=2)); nodes=[()]+[(i,) for i in range(3)]
 A=[]
 for vals in product(range(3),repeat=4):A.append(centralizer(2,dict(zip(nodes,vals))))
 dist=Counter(sum(i==j for i,j in enumerate(p)) for p in A)
 check(dist=={0:62,3:12,6:6,9:1},'literal A2 distribution')
 # C3= C2 x A2, C2 fixed counts9,3,3; linked A2 gets doubled.
 y=Counter()
 for u,count in {9:1,3:2}.items():
  for v,m in dist.items():y[u+2*v]+=count*m
 for r in range(1,7):
  m=sum(Fraction(count*x**r,243) for x,count in y.items())
  check(m.denominator==1,'Burnside integral')
 check(sum(x*c for x,c in y.items())==7*243,'one ports')
 check(sum(x*x*c for x,c in y.items())==73*243,'two ports')
 # Hostile: same support does not determine weighted moments.
 support_mean=Fraction(sum(dist),len(dist))
 uniform_mean=sum(Fraction(x*c,81) for x,c in dist.items())
 check(support_mean!=uniform_mean and uniform_mean==1,'support multiplicity hostile')
 return {'A2_distribution':dict(sorted(dist.items())),'C3_distribution':dict(sorted(y.items())),'support_mean':str(support_mean),'group_mean':str(uniform_mean)}
def cycle_inventory_checks():
 words=list(product(range(3),repeat=2)); reps=[(1,)]
 group=[]
 for v in range(3):group.append(centralizer(2,{(1,):v,(2,):(-v)%3}))
 remaining=set(range(512)); inventory=[0]*10
 while remaining:
  mask=min(remaining)
  orbit={sum(((mask>>i)&1)<<g[i] for i in range(9)) for g in group}
  remaining-=orbit;inventory[mask.bit_count()]+=1
 check(inventory==[1,5,14,30,46,46,30,14,5,1],'subset orbit inventory')
 check(sum(inventory)==192,'binary orbit total')
 fixed_only=sum(Fraction(comb(sum(i==j for i,j in enumerate(g)),3),3) for g in group)
 check(fixed_only==Fraction(86,3) and fixed_only!=inventory[3],'fixed-points not subsets')
 # Literal Burnside inventory via invariant subsets must coincide.
 via=[Fraction(0) for _ in range(10)]
 for g in group:
  for mask in range(512):
   image=sum(((mask>>i)&1)<<g[i] for i in range(9))
   if mask==image:via[mask.bit_count()]+=Fraction(1,3)
 check(via==inventory,'literal weighted Burnside')
 groups=[[(tuple((i+a)%9 for i in range(9))) for a in range(9)],
         [tuple(3*((i//3+a)%3)+(i%3+b)%3 for i in range(9)) for a,b in product(range(3),repeat=2)]]
 counts=[]; laws=[]
 from itertools import combinations
 for G in groups:
  laws.append(Counter(sum(i==g[i] for i in range(9)) for g in G))
  triples=set(combinations(range(9),3));number=0
  while triples:
   S=min(triples);triples-={tuple(sorted(g[i] for i in S)) for g in G};number+=1
  counts.append(number)
 check(laws[0]==laws[1]=={9:1,0:8},'all fixedpoint moments hostile')
 check(counts==[10,12],'regular actions subset hostile')
 return {'subset_inventory':inventory,'false_fixed_only_triples':str(fixed_only),'same_all_moments_subset_counts':counts}

if __name__=='__main__':
 report={'loop':loop_checks(),'moments':moment_checks(),'cycles':cycle_inventory_checks()};report['checks']=checks
 print(json.dumps(report,sort_keys=True,indent=2))
