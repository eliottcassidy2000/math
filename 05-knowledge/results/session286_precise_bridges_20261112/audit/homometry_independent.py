#!/usr/bin/env python3
"""Independent memoized first-vertex recursion for phase-breaking witnesses."""
from functools import lru_cache
import json

def must(x,m):
 if not x:raise RuntimeError(m)

def count(T):
 @lru_cache(None)
 def start(mask,v):
  rest=mask^(1<<v)
  if not rest:return 1
  total=0;ok=rest&T[v]
  while ok:
   bit=ok&-ok;ok-=bit;total+=start(rest,bit.bit_length()-1)
  return total
 return sum(start((1<<len(T))-1,v)for v in range(len(T)))

def run():
 n=13;D=[sum(1<<((i+d)%n)for d in range(1,7))for i in range(n)]
 A={0,1,3,9};B={1,2,5,7}
 corr=lambda C:[sum((i+d)%n in C for i in C)for d in range(n)]
 must(corr(A)==corr(B)==[4]+[1]*12,'cyclic homometry')
 values=[]
 for C,want in[(A,21508189),(B,21508189),(A|{8},23572029),(B|{8},23562137)]:
  T=[r|(0 if i in C else 1<<n)for i,r in enumerate(D)]+[sum(1<<i for i in C)]
  h=count(T);must(h==want,'literal path count');values.append(dict(cut=sorted(C),H=h))
 extraD=[sum(1<<((i+d)%13)for d in (12,2,10,9,8,7))for i in range(13)]
 extra=[]
 for mask,want in ((2135,24389195),(3180,24389151)):
  T=[r|(0 if mask>>i&1 else 1<<13)for i,r in enumerate(extraD)]+[mask]
  got=count(T);must(got==want,'gap44 witness');extra.append(got)
 return dict(gap44_pair=extra,status='PASS' ,core_H=count(D),witnesses=values,released_gap=values[2]['H']-values[3]['H'])
if __name__=='__main__':print(json.dumps(run(),indent=2))
