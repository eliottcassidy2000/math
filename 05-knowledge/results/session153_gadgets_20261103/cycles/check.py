from itertools import combinations, permutations, product
from collections import Counter
for bits in product((0,1),repeat=3):
 S={i if bits[i-1] else 7-i for i in (1,2,3)}
 edge=lambda x,y: (y-x)%7 in S
 paths=sum(all(edge(q[i],q[i+1]) for i in range(6)) for q in permutations(range(7)))
 cyc={}
 for k in (3,5,7):
  cc=[]
  for subset in combinations(range(7),k):
   root=min(subset)
   for rest in permutations(x for x in subset if x!=root):
    q=(root,)+rest
    if all(edge(q[i],q[(i+1)%k]) for i in range(k)):cc.append(frozenset(q))
  cyc[k]=cc
 pairs=sum(not a&b for a,b in combinations(cyc[3],2))
 ocf=1+2*sum(map(len,cyc.values()))+4*pairs
 if paths!=ocf: raise AssertionError((paths,ocf))
 print(''.join(map(str,bits)),sorted(S),*[len(cyc[k]) for k in (3,5,7)],pairs,paths)
