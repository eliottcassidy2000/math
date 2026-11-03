from itertools import permutations
from collections import Counter

def proof(n):
    pairs = [(i,j) for i in range(n) for j in range(i+2,n)]
    F=Counter(); G=Counter()
    for mask in range(1 << len(pairs)):
        back={(j,i) for k,(i,j) in enumerate(pairs) if mask & 1<<k}
        def edge(i,j): return (i,j) in back or (i<j and (j,i) not in back)
        strong=all(any(edge(j,i) for i in range(k) for j in range(k,n)) for k in range(1,n))
        H=sum(all(edge(p[i],p[i+1]) for i in range(n-1)) for p in permutations(range(n)))
        F[H]+=1
        if strong: G[H]+=1
    print(n,'pairs',len(pairs),'all',sorted(F.items()),'strong',sorted(G.items()),'total',sum(G.values()))
for n in [3,4,5,6]:proof(n)
