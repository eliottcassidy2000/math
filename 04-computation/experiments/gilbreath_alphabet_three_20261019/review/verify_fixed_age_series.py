from collections import deque
from fractions import Fraction
import json

def next_state(code,bit,width):
 out=0
 for j in range(width):
  bit=abs(bit-((code>>(2*j))&3));out|=bit<<(2*j)
 return out
rows=[];slope=Fraction(0)
for age,h,c in ((2,2,1),(4,4,3),(6,4,1),(8,8,7),(10,8,2),(12,8,3),(14,8,1)):
 width=age+1;start=(1<<(2*width))-1;known={start};queue=deque([start]);edges={};parents={}
 while queue:
  x=queue.popleft(); ys=tuple(next_state(x,b,width) for b in(0,1));edges[x]=ys
  for y in ys:
   parents.setdefault(y,set()).add(x)
   if y not in known:known.add(y);queue.append(y)
 acc={x for x in known if all(((x>>(2*j))&3)<=1 for j in range(age)) and ((x>>(2*age))&3)>=2}
 co=set(acc);queue=deque(acc)
 while queue:
  for x in parents.get(queue.popleft(),set()):
   if x not in co:co.add(x);queue.append(x)
 m=len(co);last=m+2*h
 # Backward counts from every state, independently of producer's forward propagation.
 vals={x:int(x in acc) for x in co}; observed=[]
 for L in range(last+1):
  actual=vals.get(start,0); expected=c*((L+h)//(2*h))
  if actual!=expected:raise RuntimeError((age,L,actual,expected))
  observed.append(actual)
  vals={x:sum(vals.get(y,0) for y in edges[x]) for x in co}
 slope+=Fraction(c,h*2**age)
 rows.append(dict(age=age,reachable=len(known),coaccessible=m,checked_through_length=last,coefficients=observed))
print(json.dumps(dict(status='PASS',method='Independent packed-state BFS and backward coefficient propagation, with exact generating-function degree certificate',rows=rows,slope=str(slope),ratio=str(1+slope/2)),indent=2))
