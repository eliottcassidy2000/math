#!/usr/bin/env python3
from math import isqrt
from verify import verify

def prime(n):return n>=2 and all(n%d for d in range(2,isqrt(n)+1))
rows={}
with open('witnesses.txt') as f:
 for line in f:
  h,c,cut=map(int,line.split())
  if h in rows or not (90003<=h<=100001) or h%2==0 or not (0<=c<2**55 and 0<=cut<2**11):raise ValueError(('bad row',h,c,cut))
  verify(h,c,cut,12)
  rows[h]=(c,cut)
for n in range(90003,100002,2):
 if n in rows:continue
 split=next(((d,n//d) for d in range(3,isqrt(n)+1,2) if n%d==0 and d not in (7,21) and n//d not in (7,21) and n//d<=90001),None)
 if split is None:raise ValueError(('uncovered',n))
print('PASS:',len(rows),'strong order-12 witnesses; all odd counts 90003 through 100001 attained')
