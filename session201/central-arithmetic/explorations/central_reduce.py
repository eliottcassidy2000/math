from itertools import product
from collections import deque
import json

def component(w):
 seen={w};todo=[w]
 for x in todo:
  for i in range(len(x)-1):
   a=x[i:i+2]
   if a in ['JZ','OJ']:
    y=x[:i]+('OJ'if a=='JZ'else'JZ')+x[i+2:]
    if y not in seen:seen.add(y);todo.append(y)
 return seen
kept=[];seen=set();dominated=[]
for n in range(9):
 W=['']if n==0 else ['J']if n==1 else ['J'+''.join(w)+'J'for w in product('OZJ',repeat=n-2)]
 for w in W:
  if w in seen:continue
  cc=component(w);seen.update(cc)
  loose=[c for c in cc if c and(c[0]!='J'or c[-1]!='J')]
  if loose:
   x=min(loose);q=x.strip('OZ');dominated.append((w,x,q))
  else:kept.append(min(cc))
print('kept',len(kept),kept)
print('counts',{i:sum(len(w)==i for w in kept)for i in range(9)})
json.dump({'kept':kept,'dominated':dominated},open(__file__.replace('central_reduce.py','central_reduction.json'),'w'),indent=2)
