from collections import deque
from itertools import product
import json

def front(v,a):
 w=[]
 for z in v:a=abs(a-z);w.append(a)
 return tuple(w)
def regex(s,a):
 if s<4:return (s+1)%4 if a==0 else (4 if s==0 else 7)
 if s==4:return 5 if a else 7
 if s==5:return 6 if a else 7
 if s==6:return 7 if a else 6
 return 7
start=((3,3,3),0);q=deque([start]);seen={start};trans=[]
while q:
 v,s=q.popleft()
 accept=max(v[:2])<=1 and v[2]>=2
 if accept!=(s in (5,6)):raise RuntimeError(('language',v,s))
 for a in(0,1):
  w=(front(v,a),regex(s,a));trans.append([list(v),s,a,list(w[0]),w[1]])
  if w not in seen:seen.add(w);q.append(w)
checks=[]
for L in range(0,17):
 c=0
 for w in product((0,1),repeat=L):
  v=(3,3,3)
  for a in w:v=front(v,a)
  c+=int(max(v[:2])<=1 and v[2]>=2)
 if c!=(L+2)//4:raise RuntimeError(('count',L,c))
 checks.append([L,c])
print(json.dumps(dict(status='PASS',claim='Exact DFA equivalence for all word lengths, plus bounded enumeration controls',product_states=len(seen),front_states=len({v for v,s in seen}),transitions=trans,finite_counts=checks),indent=2))
