from collections import deque
from math import gcd
def find(a,b,limit=None):
 g=gcd(a,b);a//=g;b//=g
 if (a-b)&1:return None,0
 p=((a-b)//2,b-a)
 seen={p:None};Q=deque([p])
 while Q:
  u,v=q=Q.popleft()
  midd=[m for m in (None,0,1) if (u==v if m is None else u+(a-b)*m==2*v)]
  if midd:
   L,R=[],[]
   while seen[q] is not None:
    q,x,y=seen[q];L.append(x);R.append(y)
   w=(1,*L[::-1],*((midd[0],) if midd[0] is not None else ()),*R,1)
   return sum(x<<i for i,x in enumerate(w)),len(seen)
  for x in (0,1):
   for y in (0,1):
    z=u+a*x-b*y
    if z&1:continue
    r=z//2,2*v-(a*y-b*x)
    if not -b<=r[1]<=a or r in seen:continue
    seen[r]=(q,x,y);Q.append(r)
  if limit and len(seen)>limit:return 'limit',len(seen)
 return None,len(seen)
if __name__=='__main__':
 for a,b in ((1847,1967),(3,5),(7,19),(13,31)):
  print(a,b,find(a,b,1000000))
