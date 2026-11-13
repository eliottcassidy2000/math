from itertools import product
from collections import defaultdict

def lawful(scores):
    dp = {(0, 0): 0}
    for s in scores:
        new = {}
        for (k, prior), w in dp.items():
            for now in (0, 1):
                q = (k + now, now)
                val = w + now*s - prior*(1-now)
                new[q] = min(new.get(q, val), val)
        dp = new
    return all(w >= k*(k-1)//2 for (k, _), w in dp.items())

def laurent_support(n):
    m = n//2
    factors = []
    for i in range(m):
        for j in range(i+1, m):
            if j >= i+2:
                factors.append(tuple((k == j) - (k == i) for k in range(m)))
            factors.append(tuple(-int(k in (i, j)) for k in range(m)))
        for _ in range(int(i < m-1 or n%2) + int(n%2 and i < m-1)):
            factors.append(tuple(-int(k == i) for k in range(m)))
    supp = {tuple(n-1-i for i in range(m))}
    for step in factors:
        supp |= {tuple(a+b for a,b in zip(word,step)) for word in supp}
    return supp

from itertools import combinations_with_replacement
for m in range(1,10):
 V={}
 for A0 in combinations_with_replacement(range(m+1),m):
  A=A0[::-1];ss=[m-x for x in A]+[m]+[m+x for x in A[::-1]]
  ex=[sum(ss[:i])-i*(i-1)//2 for i in range(1,m+1)]
  if min(ex)>=0:V[A]=min(ex)>0
 bad=[];total=0;rev=0
 for A in V:
  for j in set(A)-{0}:
   L=list(A);L.remove(j); B=tuple(sorted(L+[j-1],reverse=True))
   if B not in V or not (V[A] or V[B]):continue
   total+=1
   v=[j]+L; S=tuple(m+x for x in v)+ (m,)+tuple(m-x for x in v[::-1])
   if not lawful(S):
    v=[1-j]+L;S=tuple(m+x for x in v)+(m,)+tuple(m-x for x in v[::-1])
    if not lawful(S):bad.append((A,j,B))
    else:rev+=1
 print(m,len(V),total,'rev',rev,'fail',len(bad),bad[:5],flush=True)
