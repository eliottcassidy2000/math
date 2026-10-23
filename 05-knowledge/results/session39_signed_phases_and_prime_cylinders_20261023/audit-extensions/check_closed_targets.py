"""Independent backwards construction; no producer imports."""
from itertools import product
import json
cases=0
for k in range(4):
 for w in product(range(1,5),repeat=k):
  for s,y in product((-1,1),(1,5,17)):
   period=2*3**k;hits=[]
   # Construct backwards using exact integer divisibility at EACH step,
   # rather than injecting the producer's total numerator congruence.
   for b in range(1,period+1):
    z=2**b*y-s
    if z%3:continue
    z//=3
    for a in reversed(w):
     z=2**a*z-s
     if z%3:break
     z//=3
    else:hits.append((b,z))
   if len(hits)!=1:raise ValueError('phase uniqueness')
   b,n=hits[0]
   while n<=0:
    b+=period;z=y
    for a in reversed(w+(b,)):z=(2**a*z-s)//3
    n=z
   z=n;observed=[]
   for _ in range(k+1):
    v=3*z+s;a=0
    while v%2==0:v//=2;a+=1
    observed.append(a);z=v
   if tuple(observed)!=w+(b,) or z!=y:raise ValueError('forward replay')
   cases+=1
print(json.dumps({'status':'PASS','independent_backward_cases':cases,'max_prefix_length':3,'alphabet':[1,2,3,4],'signs':[-1,1],'targets':[1,5,17]},indent=2))
