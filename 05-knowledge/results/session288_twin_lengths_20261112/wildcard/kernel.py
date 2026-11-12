from functools import lru_cache
from collections import defaultdict
from math import comb

def skeleton(rows,O,ports):
 fixed=[v for v in range(len(rows)) if v not in ports];pi={v:i for i,v in enumerate(ports)};fi={v:i for i,v in enumerate(fixed)};full=(1<<len(fixed))-1
 @lru_cache(None)
 def rec(mask,last,phase):
  out=defaultdict(int)
  if mask==full:out[(0,0,phase)]+=1
  for v in range(len(rows)):
   if v==last:continue
   if v in fi and mask>>fi[v]&1:continue
   shifts=(int(v==ports[0]),int(v==ports[1]));nm=mask|(1<<fi[v]) if v in fi else mask
   if rows[last]>>v&1:
    for (a,b,p),c in rec(nm,v,phase).items():out[a+shifts[0],b+shifts[1],p]+=c
   if phase==0 and last in O and v in O:
    for (a,b,p),c in rec(nm,v,1).items():out[a+shifts[0],b+shifts[1],p]+=c
  return dict(out)
 total=defaultdict(int)
 for v in range(len(rows)):
  mask=(1<<fi[v]) if v in fi else 0;sa=int(v==ports[0]);sb=int(v==ports[1])
  for (a,b,p),c in rec(mask,v,0).items():total[a+sa,b+sb,p]+=c
 d=defaultdict(int)
 for (a,b,p),c in total.items():
  if a and b:d[a,b]+=c*(1 if p==0 else -1)
 return {k:v for k,v in d.items() if v}
def exp(poly):
 out=defaultdict(int)
 for (r,s),c in poly.items():
  for a in range(1,r+1):
   for b in range(1,s+1):out[a,b]+=c*(-1)**(r+s-a-b)*comb(r,a)*comb(s,b)
 return {k:v for k,v in out.items() if v}
def evaluate(poly,a,b):return sum(c*x**a*y**b for (x,y),c in poly.items())

def observable_skeleton(rows,O,ports,kind):
 """kind H, U or W; ordered nonempty port runs, selected fixed endpoints."""
 fixed=[v for v in range(len(rows)) if v not in ports];fi={v:i for i,v in enumerate(fixed)};full=(1<<len(fixed))-1
 @lru_cache(None)
 def rec(mask,last,phase):
  out=defaultdict(int)
  if mask==full and ((kind=='H' and phase==0) or (kind in ('U','W') and phase==1)):out[(0,0)]+=1
  for v in range(len(rows)):
   if v==last or (v in fi and mask>>fi[v]&1):continue
   sa=int(v==ports[0]);sb=int(v==ports[1]);nm=mask|(1<<fi[v]) if v in fi else mask
   if rows[last]>>v&1:
    for (a,b),c in rec(nm,v,phase).items():out[a+sa,b+sb]+=c
   if kind!='H' and phase==0 and v in O and (kind=='U' or last in O):
    for (a,b),c in rec(nm,v,1).items():out[a+sa,b+sb]+=c
  return dict(out)
 out=defaultdict(int)
 for v in range(len(rows)):
  mask=1<<fi[v] if v in fi else 0;sa=int(v==ports[0]);sb=int(v==ports[1])
  for phase in ([0,1] if kind=='U' and v in O else [0]):
   for (a,b),c in rec(mask,v,phase).items():
    if a+sa and b+sb:out[a+sa,b+sb]+=c
 return dict(out)
