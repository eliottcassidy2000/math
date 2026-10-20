from itertools import permutations,product
from math import lcm

def order(p):
 s=set();o=1
 for i in range(len(p)):
  if i not in s:
   j=i;c=0
   while j not in s:s.add(j);j=p[j];c+=1
   o=lcm(o,c)
 return o
for w in range(3,10):
 patterns=list(product((0,1),repeat=w-2));covered={};odd={}
 for mid in permutations(range(w-1)):
  p=mid+(w-1,)
  o=order(p)
  if o==1:continue
  req={}
  for a,b in zip(p,p[1:]):
   d=b-a;v=int(d>0);d=abs(d)
   if (d==1 and not v) or (d in req and req[d]!=v):break
   req[d]=v
  else:
   for row in patterns:
    if all((1 if d==1 else row[d-2])==v for d,v in req.items()):
     if o%2==0:covered.setdefault(row,p)
     else:odd.setdefault(row,p)
 residual=[row for row in patterns if row not in covered]
 print(w,len(patterns),len(covered),'residual',len(residual),['1'+''.join(map(str,r)) for r in residual][:100],flush=True)
