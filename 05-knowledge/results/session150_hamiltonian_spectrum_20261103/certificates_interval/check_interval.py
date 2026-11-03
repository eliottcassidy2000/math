#!/usr/bin/env python3
"""Independent exact verification of one order-12 strong witness per odd H."""
from verify import verify
expected=list(range(80407,80874,2))
rows=[]
with open('witnesses.txt') as f:
 for line in f:
  vals=list(map(int,line.split()))
  if len(vals)!=3:raise ValueError('wrong row width')
  h,c,cut=vals
  if not (0<=c<2**55 and 0<=cut<2**11):raise ValueError('out-of-range code')
  verify(h,c,cut,12)
  rows.append(h)
if rows!=expected:raise ValueError('missing, duplicate, out-of-order, or extra target')
print('PASS: 234 strong order-12 tournaments, H every odd integer 80407 through 80873')
