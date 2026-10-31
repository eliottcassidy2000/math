"""Independent H=63 witness check, without constructor/profile imports."""
from itertools import permutations
import hashlib,json
# Seed A has one vertex, external x loses to it. Appended A vertices form an
# ordinal singleton chain; each bit is precisely x -> that appended vertex.
bits=[0,1,0,0,1,1,0,0]
n=len(bits)+1;x=n-1
out=[0]*n
for u in range(x):
 for v in range(u+1,x):out[u]|=1<<v
 if bits[u]:out[x]|=1<<u
 else:out[u]|=1<<x
count=0
for p in permutations(range(n)):
 if all(out[u]>>v&1 for u,v in zip(p,p[1:])):count+=1
D=[[0]*n for _ in range(1<<n)]
for u in range(n):D[1<<u][u]=1
for mask in range(1,1<<n):
 for u,value in enumerate(D[mask]):
  if value:
   todo=out[u]&~mask
   while todo:
    bit=todo&-todo;todo-=bit;D[mask|bit][bit.bit_length()-1]+=value
if count!=63 or sum(D[-1])!=63:raise RuntimeError('bad H63 witness')
# The displayed factorization branches reconstruct all finite sample targets;
# this is a control, not the universal proof.
branches=0
for N in range(1,10000,2):
 k=0;m=N
 while m%3==0:k+=1;m//=3
 if m==7:
  if k<2:
   if N not in (7,21):raise RuntimeError('wrong excluded case')
  elif 63*3**(k-2)!=N:raise RuntimeError('63 branch')
 elif m*3**k!=N:raise RuntimeError('coprime branch')
 branches+=1
print(json.dumps(dict(status='PASS',order=n,outgoing_masks=out,literal_permutation_H=count,independent_subset_dp_H=sum(D[-1]),factorization_controls=branches,adjacency_sha256=hashlib.sha256(json.dumps(out,separators=(',',':')).encode()).hexdigest()),sort_keys=True))
