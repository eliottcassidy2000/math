#!/usr/bin/env python3
from itertools import combinations, product
from math import isqrt
import json,random
CHECKS=0

def check(v,msg):
 global CHECKS
 CHECKS+=1
 if not v: raise RuntimeError(msg)

def sq(x):return x>=0 and isqrt(x)**2==x

def realize(parents):
 labels=[1]
 for k in range(1,len(parents)):
  old=labels[:];H=max(old);lam=1
  while lam<32*(H+1):lam*=4
  labels=[4*lam*b for b in old]
  if not parents[k]:x=2
  else:
   vals=sorted(old[p] for p in parents[k])
   b=vals[0];c=vals[1] if len(vals)==2 else H+1
   x=(lam*(c-b)-1)**2-4*lam*b
  check(x>0 and x not in labels,'bad new label')
  labels.append(x)
  check(max(labels)<=65536*H**4,'step height')
  for i,j in combinations(range(k+1),2):
   check(sq(labels[i]+labels[j])==(i in parents[j]),('graph mismatch',parents,i,j))
 return max(labels).bit_length()
counts={}
for n in range(1,6):
 choices=[[()]+[(j,) for j in range(i)]+list(combinations(range(i),2)) for i in range(n)]
 c=0
 for parents in product(*choices):realize(parents);c+=1
 counts[n]=c
rng=random.Random(20261019)
bits=[]
for z in range(100):
 n=rng.choice([6,7])
 parents=[()]
 for i in range(1,n):parents.append(tuple(sorted(rng.sample(range(i),rng.randrange(min(i,2)+1)))))
 bits.append(realize(parents))
# Exhaust all possible common-neighbor candidates of scaled 1 and 2,
# using complete factor pairs for their difference S^2, then check 3.
common_candidates=0
for S in range(1,501):
 d=S*S
 for u in range(1,isqrt(d)+1):
  if d%u:continue
  v=d//u
  if v<=u or (v-u)%2:continue
  r=(v-u)//2;x=r*r-d
  if x<=0:continue
  check(sq(x+d) and sq(x+2*d),'factorization')
  check(not sq(x+3*d),'Fermat hostile')
  common_candidates+=1
print(json.dumps({'status':'PASS','always_active_checks':CHECKS,'exhaustive_ordered_two_parent_graphs':counts,'random_order_six_or_seven':100,'maximum_random_label_bit_length':max(bits),'triple_seed_scales_through':500,'all_common_neighbors_of_first_two_tested':common_candidates},indent=2,sort_keys=True))
