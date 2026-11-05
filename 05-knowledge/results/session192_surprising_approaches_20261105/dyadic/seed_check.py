#!/usr/bin/env python3
import json
from pathlib import Path

def dp(r):
 n=len(r);d=[[0]*n for _ in range(1<<n)]
 for i in range(n):d[1<<i][i]=1
 for s in range(1,1<<n):
  for v in range(n):
   if not (s>>v&1):continue
   p=s^(1<<v)
   if p:d[s][v]=sum(d[p][u] for u in range(n) if r[u]>>v&1)
 return sum(d[-1])

def graph(w,apices):
 r=[0]*(len(w)+apices)
 if apices==2:r[0]=2
 for j,c in enumerate(map(int,w),apices):
  for i in range(apices):
   if c>>i&1:r[i]|=1<<j
   else:r[j]|=1<<i
  for i in range(apices,j):r[i]|=1<<j
 return r
out=[]
for h,w,a in [(3,'10',1),(5,'100',1),(23,'11000',1),(63,'1001100',1),(67,'12301',2)]:
 r=graph(w,a);actual=dp(r)
 if actual!=h:raise RuntimeError((h,w,actual))
 out.append({'H':h,'word':w,'apices':a,'rows':r})
# Exhaust the exact conditional arithmetic descent, without claiming physical closure.
for n in range(1,100001,2):
 if n in (1,3,5,7,21,23,25,63,65,67):continue
 q=(n if n%3==0 else n-2 if n%3==2 else n-4)//3
 if not (1<=q<n and q%2==1 and q not in (7,21)):raise RuntimeError(('descent',n,q))
path=Path(__file__).with_name('seed_receipt.json');path.write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({'status':'PASS','seed_checks':len(out),'conditional_descent_odd_targets':50000},sort_keys=True))
