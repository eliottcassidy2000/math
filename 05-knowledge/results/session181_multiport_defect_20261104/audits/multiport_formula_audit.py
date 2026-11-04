import itertools,json
from functools import lru_cache
from pathlib import Path

def profile(a,w):
 n=len(a);N=(1<<n)-1
 @lru_cache(None)
 def e(mask,v):
  rem=mask^(1<<v)
  return 1 if not rem else sum(e(rem,u) for u in range(n) if rem>>u&1 and a[u]>>v&1)
 @lru_cache(None)
 def h(mask):return sum(e(mask,v) for v in range(n) if mask>>v&1) if mask else 1
 return h(N),sum(e(mask,w)*h(N^mask) for mask in range(1,N+1) if mask>>w&1)

def suffix(a,W):
 n=len(a);F=J=0
 for p in itertools.permutations(range(n)):
  for t in range(1,n+1):
   if not W>>p[t-1]&1:continue
   if not all(a[p[i]]>>p[i+1]&1 for i in range(n-1) if i!=t-1):continue
   F+=1
   if t<n and W>>p[t]&1:J+=1
 return profile(a,0)[0],F,J
cases=0
for n in range(1,5):
 for code in range(1<<(n*(n-1)//2)):
  a=[0]*n
  for bit,(i,j) in enumerate(itertools.combinations(range(n),2)):
   if code>>bit&1:a[i]|=1<<j
   else:a[j]|=1<<i
  for W in range(1,1<<n):
   h,F,J=suffix(a,W)
   for s in range(1,4):
    root=n+s;b=[row|((1<<root) if W>>i&1 else 0) for i,row in enumerate(a)]
    b += [((1<<n)-1)|sum(1<<j for j in range(i+1,root)) for i in range(n,root)]
    b += [((1<<root)-1)^W]
    H,E=profile(b,root)
    if (H,E,E-H)!=(h-J+(1<<s)*F,h+(1<<s)*F,J):raise RuntimeError((a,W,s,H,E,h,F,J))
    cases+=1
out={'status':'PASS','independent_multiport_formula_cases':cases,'scope':'all suffix tournaments of orders1..4, every nonempty port set, prefix lengths1..3'}
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out))
