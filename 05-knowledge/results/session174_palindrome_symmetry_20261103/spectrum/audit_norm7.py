#!/usr/bin/env python3
import itertools,json
from common import *
def valuation(n,p):
 require(n>0,'positive valuation argument');v=0
 while n%p==0:n//=p;v+=1
 return v
records=[]
for n in range(1,5):
 pairs=list(itertools.combinations(range(n),2))
 for code in range(1<<len(pairs)):
  a=[0]*n
  for k,(i,j) in enumerate(pairs):
   if code>>k&1:a[i]|=1<<j
   else:a[j]|=1<<i
  p=cover(a);x=evalF(p,1);y=evalF(p,2)-x;h=paths(wrap(a,1))[-1]
  require(h==y*y+2*x*x,'literal Hamilton norm formula')
  require(valuation(h,7)==2*min(valuation(x,7),valuation(y,7)),'exact norm valuation')
  records.append({'n':n,'code':code,'F1':x,'fixed':y,'H':h,'v7':valuation(h,7)})
require(5 not in {x*x%7 for x in range(7)},'minus two is nonsquare')
for a in range(49):
 for b in range(49):
  q=b*b+2*a*a
  require((q%7==0)==(a%7==b%7==0),'mod7 anisotropy')
  require(q%7!=0 or q%49==0,'no valuation one modulo49')
print(json.dumps({'literal_tournaments':len(records),'mod49_pairs':49**2,'records':records},sort_keys=True,indent=2))
