"""Exact maximum prefix-antichain size in the certified A/B/C schedule tree."""
from functools import lru_cache

def count_frontier(N):
 if N<1062:raise ValueError('cycle-count construction requires N >= 1062')
 budget=2*N+1
 @lru_cache(None)
 def count(R):
  if 85*25*R>budget:return 1
  return 2*count(25*R)+(count(49*R) if 85*49*R<=budget else 0)
 result=count(1)
 return {'N':N,'distinct_cycles_lower_bound':result,'product_states':count.cache_info().currsize}

def frontier(N):
 if N<1062:raise ValueError('cycle-count construction requires N >= 1062')
 todo=[('',1)]
 while todo:
  word,R=todo.pop()
  if 85*25*R<=2*N+1:
   todo.extend((word+a,25*R) for a in 'AB')
   if 85*49*R<=2*N+1:todo.append((word+'C',49*R))
  else:yield word,R

if __name__=='__main__':
 import argparse,json
 p=argparse.ArgumentParser();p.add_argument('N',type=int);p.add_argument('--list',action='store_true');a=p.parse_args();r=count_frontier(a.N)
 if a.list:r['prefix_frontier']=[{'word':w,'product':R} for w,R in frontier(a.N)]
 print(json.dumps(r,indent=2))
