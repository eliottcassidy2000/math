def require(x):
 if not x: raise AssertionError("bridge check failed")
from itertools import combinations, permutations
from collections import Counter
from math import factorial
for n in range(1,5):
 pairs=list(combinations(range(n),2)); N=len(pairs); ix={e:i for i,e in enumerate(pairs)}
 def forward(s,a,b):
  k=ix[tuple(sorted((a,b)))]; return bool(s>>k&1)==(a<b)
 def paths(s): return [p for p in permutations(range(n)) if all(forward(s,a,b) for a,b in zip(p,p[1:]))]
 def partner(s,p):
  used={tuple(sorted((a,b))) for a,b in zip(p,p[1:])}; mask=sum(1<<i for i,e in enumerate(pairs) if e not in used); return s^mask
 def reflect(s):
  t=0
  for i,(a,b) in enumerate(pairs):
   if not forward(s,n-1-a,n-1-b): t|=1<<i
  return t
 trans=Counter(); total=0
 for s in range(1<<N):
  P=paths(s); total+=len(P)
  for p in P:
   t=partner(s,p); require(p in paths(t) and partner(t,p)==s)
   trans[s,t]+=1
   q=tuple(n-1-x for x in p[::-1]); require(q in paths(reflect(s)))
   require(partner(reflect(s),q)==reflect(t))
  require(len(P)==len(paths(reflect(s))))
 require(all(c==trans[t,s] for (s,t),c in trans.items()))
 require(total==factorial(n)*2**(N-n+1))
 print(f'n={n}: {1<<N} labelled tournaments, {total} marked pairs; transition balance and reflection verified')
