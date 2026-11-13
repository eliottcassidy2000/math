#!/usr/bin/env python3
import itertools,json,math
checks=0
def require(q):
 global checks;checks+=1
 if not q:raise RuntimeError('certificate failure')
def weights(bits):
 w=[0]*9
 def arc(a,b):
  if a==6:return bits[b]=='1'
  if b==6:return bits[a]=='0'
  return a<b
 for p in itertools.permutations(range(7)):
  optional=required=0;ok=True
  for g in range(8):
   legal=(g==0 or p[g-1]<3) and (g==7 or p[g]>=3)
   optional+=legal
   if 0<g<7 and not arc(p[g-1],p[g]):
    if not legal:ok=False;break
    required+=1
  if ok:
   for k in range(required,optional+1):w[k]+=math.comb(optional-required,k-required)
 return w
x=weights('101101');y=weights('011110')
require(x==[15,31,15,0,0,0,0,0,0]);require(y==[17,33,15,0,0,0,0,0,0])
print(json.dumps({'exterior_orders_per_scaffold':math.factorial(7),'X':x,'Y':y,'checks':checks},indent=2))
