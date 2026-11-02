"""Independent leaf: literal permutations, no import from any research engine."""
from itertools import permutations,combinations
import json,hashlib
A=(6,112,74,35,13,21,57)
B=(38,116,72,35,77,20,41)
PERMS=tuple(permutations(range(7)))

def check(v,msg):
 if not v:raise RuntimeError(msg)
def literal(a):
 return sum(all(a[p[i]]&(1<<p[i+1]) for i in range(6)) for p in PERMS)
def isomorphism(a,b):
 for p in PERMS:
  if all(bool(a[i]&(1<<j))==bool(b[p[i]]&(1<<p[j])) for i in range(7) for j in range(i+1,7)):return p
 return None
def strong(a):
 for r in range(7):
  seen=1<<r
  while True:
   new=seen
   for i in range(7):
    if seen&(1<<i):new|=a[i]
   if new==seen:break
   seen=new
  if seen!=127:return False
 return True
check(literal(A)==147 and literal(B)==145,'literal witness path counts')
check(strong(A) and strong(B),'strong positive witnesses')
zeros=[];image=set();cases=0
for r in range(7):
 others=[v for v in range(7) if v!=r]
 for bits in range(64):
  a=list(A)
  for i,v in enumerate(others):
   if bits&(1<<i):a[r]|=1<<v;a[v]&=~(1<<r)
   else:a[r]&=~(1<<v);a[v]|=1<<r
  a=tuple(a);h=literal(a);cases+=1
  if strong(a):
   image.add(h);check(h not in [145,149],'no adjacent strong replacement')
   if h==147:
    p=isomorphism(A,a);check(p is not None,'zero-response successor stays in same iso orbit');zeros.append((r,bits,a,p))
check(cases==448 and len(zeros)==11,'presentation and zero-response counts')
pairs=list(combinations(range(7),2));close=[]
for size in [1,2]:
 for edges in combinations(pairs,size):
  a=list(A)
  for i,j in edges:a[i]^=1<<j;a[j]^=1<<i
  h=literal(tuple(a));check(h!=145,'no one/two-arc target even without strong filter')
  close.append((edges,h))
row={'test':'independent literal closed-orbit leaf','base_H':147,'target_H':145,'one_vertex_presentations':cases,'zero_response_presentations':len(zeros),'literal_one_two_arc_comparisons':len(close),'strong_response_image':sorted(image),'zero_response_isomorphism_certificates':zeros}
print(json.dumps(row,sort_keys=True,separators=(',',':')))
print('PASS semantic_sha256='+hashlib.sha256(json.dumps(row,sort_keys=True,separators=(',',':')).encode()).hexdigest())
