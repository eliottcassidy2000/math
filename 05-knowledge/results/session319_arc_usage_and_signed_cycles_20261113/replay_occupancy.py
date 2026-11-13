from itertools import permutations, combinations
from math import factorial
from replay_layouts import A,B,t
def counts(R):
 out=[]
 for j in range(1,4):
  count=0; extras=tuple(range(7,6+j))
  for w in permutations(tuple(range(7))+extras):
   if [a for a in w if a==t or a>=7]!=[t]+list(extras):continue
   w=tuple(t if a>=7 else a for a in w)
   for i in range(len(w)+1):
    if all(R[a]>>b&1 for k,(a,b) in enumerate(zip(w,w[1:]),1) if k!=i):
     count+=1
  out.append(count)
 return out
a,b=counts(A),counts(B)
if [v-u for u,v in zip(a,b)]!=[0,12,4]:raise ValueError((a,b))
print('two-slot layout counts',a,b)
from replay_layouts import layouts
for r,want in ((1,6),(2,2)):
 for R,edge in ((A,(5,6)),(B,(6,1))):
  if sum(edge in set(zip(w,w[1:])) for w in layouts(R,r))!=want:raise ValueError((R,r,edge))
print('active flip weight verified in both slot layouts')
