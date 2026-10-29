"""All labelled tournaments through5; count paths by permutation constraints."""
from itertools import combinations,permutations
import json
checks=0
for n in range(2,6):
 es=list(combinations(range(n),2));k=n-1
 for orient in range(1<<len(es)):
  E={(a,b) if orient>>i&1 else (b,a) for i,(a,b) in enumerate(es)}
  values=[0]*(1<<k)
  for p in permutations(range(n)):
   required={};valid=True
   for a,b in zip(p,p[1:]):
    if abs(a-b)==1:
     required[min(a,b)]=int((a,b) in E)
    elif (b,a) not in E:
     valid=False;break
   if valid:
    for mask in range(1<<k):
     if all((mask>>i&1)==bit for i,bit in required.items()):values[mask]+=1
  coeff=sum((-1)**(k-mask.bit_count())*h for mask,h in enumerate(values))
  sign=1
  for i in range(k):sign*=1 if (i,i+1) in E else -1
  if coeff!=(1+(-1)**k)*sign:raise RuntimeError('top coefficient')
  checks+=1
print(json.dumps({'all_labelled_tournaments_orders_2_to_5':checks,'status':'pass'}))
