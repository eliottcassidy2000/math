from math import isqrt
from itertools import permutations
from collections import defaultdict
import json
from pathlib import Path
checks=0
def need(x,msg):
 global checks
 checks+=1
 if not x:raise RuntimeError(msg)
def starts(s,N):
 c=[0]
 for a in s:c.append(a-c[-1])
 L=max(1-a if i%2==0 else a-N for i,a in enumerate(c))
 U=min(N-a if i%2==0 else a-1 for i,a in enumerate(c))
 E=set()
 for i in range(len(c)):
  for j in range(i):
   if i%2==j%2:
    if c[i]==c[j]:return []
   else:
    den=(-1)**i-(-1)**j;num=c[j]-c[i]
    if num%den==0:E.add(num//den)
 return [t for t in range(L,U+1) if t not in E]
for N in range(1,9):
 bank=defaultdict(list)
 for p in permutations(range(1,N+1)):bank[tuple(a+b for a,b in zip(p,p[1:]))].append(p)
 for s,ps in bank.items():
  need(len(ps)<=1+(N%2==0),'multiplicity')
  need(starts(s,N)==sorted(p[0] for p in ps),'enumerator')
  if len(ps)==2:
   p,q=sorted(ps);d=q[0]-p[0]
   need(N%(2*d)==0,'divisor')
   need(set(p[::2])=={v for v in range(1,N+1) if (v-1)//d%2==0},'blocks')
   need(tuple(v+d if i%2==0 else v-d for i,v in enumerate(p))==q,'twin')
cert=[]
for N in range(2,36,2):
 for d in range(1,N//2+1):
  if N%(2*d):continue
  adj={a:[b for b in range(1,N+1) if (a-1)//d%2!=(b-1)//d%2 and isqrt(a+b)**2==a+b] for a in range(1,N+1)}
  iso=[a for a in adj if not adj[a]];leaves=[a for a in adj if len(adj[a])==1]
  seen={1};todo=[1]
  for a in todo:
   for b in adj[a]:
    if b not in seen:seen.add(b);todo.append(b)
  need(bool(iso) or len(leaves)>2 or len(seen)<N,'lower threshold')
  cert.append(dict(N=N,d=d,isolated=iso,leaves=leaves,component_of_1=sorted(seen)))
w=json.loads((Path(__file__).parent/'twin_witness.json').read_text());N=w['N'];p=w['path'];q=w['twin']
need(sorted(p)==list(range(1,N+1))==sorted(q),'permutations')
s=[a+b for a,b in zip(p,p[1:])]
need(s==[a+b for a,b in zip(q,q[1:])],'word')
need(all(isqrt(a)**2==a for a in s),'square edges')
need(starts(s,N)==sorted([p[0],q[0]]),'exact twin fiber')
for m in range(1,101):
 for t in range(1,11):
  a=m*m*t*t;b=3*a
  need(a%m==b%m and isqrt(a+b)**2==a+b,'periodic hostile')
print(json.dumps({'status':'PASS','checks':checks,'permutation_orders':list(range(1,9)),'minimal_square_twin_order':36,'lower_certificates':cert,'square_clock':[isqrt(a) for a in s]},indent=2))
