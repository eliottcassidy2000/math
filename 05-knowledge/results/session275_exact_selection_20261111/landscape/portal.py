from experiment import hcount
from itertools import combinations,product,permutations
from collections import deque
import json

def canon(a):
 n=len(a);blocks=[[i for i,x in enumerate(a) if x.bit_count()==d] for d in range(n)]
 best=None
 for pp in product(*(permutations(c) for c in blocks if c)):
  p=sum(pp,());r=tuple(sum(1<<j for j in range(n) if a[p[i]]>>p[j]&1) for i in range(n));code=sum(x<<(8*i) for i,x in enumerate(r))
  if best is None or code<best:best=code;out=r
 return out
seed=canon([16,5,17,7,74,31,175,63]);es=list(combinations(range(8),2));moves=[(i,) for i in range(28)]+list(combinations(range(28),2));queue=deque([seed]);seen={seed:None};counts={seed:63};tries=0;found=None
while queue and len(seen)<3000:
 a=queue.popleft();h=counts[a]
 if h<=61 and any(x.bit_count() in(0,7) for x in a):found=a;break
 for ids in moves:
  r=list(a)
  for z in ids:
   i,j=es[z];r[i]^=1<<j;r[j]^=1<<i
  hr=hcount(r);tries+=1
  if not(23<=hr<=63 and abs(hr-h)<=2):continue
  common=list(a)
  for z in ids:
   i,j=es[z];common[i]&=~(1<<j);common[j]&=~(1<<i)
  if not hcount(common):continue
  c=canon(r)
  if c not in seen:seen[c]=(a,ids,r);counts[c]=hr;queue.append(c)
chain=[]
if found:
 cur=found
 while cur!=seed:
  a,ids,raw=seen[cur];chain.append({'before':a,'before_H':counts[a],'flips':[es[z] for z in ids],'after_raw':raw,'after_canonical':cur,'after_H':counts[cur]});cur=a
 chain.reverse()
print(json.dumps({'seed':seed,'visited_classes':len(seen),'tested_moves':tries,'found':found,'chain':chain},indent=2))
