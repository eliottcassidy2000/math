from itertools import combinations, permutations
from collections import Counter
n=5;e=list(combinations(range(n),2));t=[frozenset(s) for s in combinations(range(n),3)]
yes=set()
for x in range(1<<len(e)):
 a={};Q=0
 for j,(i,k) in enumerate(e):a[i,k]=x>>j&1;a[k,i]=1-a[i,k]
 for j,s in enumerate(t):
  i,k,l=sorted(s)
  if a[i,k]==a[k,l]==a[l,i]:Q|=1<<j
 yes.add(Q)
good=[]
for q in range(1<<len(t)):
 if q in yes:continue
 b=[s for i,s in enumerate(t) if q>>i&1]
 if all(sum(s<=set(S) for s in b)<=2 for S in combinations(range(n),4)):good.append(b)
perm=[{s:t.index(frozenset(p[i] for i in s)) for s in t} for p in permutations(range(n))]
c=Counter(min(sum(1<<r[s] for s in b) for r in perm) for b in good)
for j,k in sorted(c.items()):print(k,[tuple(sorted(s)) for i,s in enumerate(t) if j>>i&1])
