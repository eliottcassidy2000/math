"""Check hepta-census geometric obstructions; pass the exact 9,5 census."""
from itertools import combinations
from pathlib import Path
import sys
P=list(combinations(range(9),2))
count=0
for line in Path(sys.argv[1]).read_text().splitlines():
    if not line[0].isdigit():continue
    q=int(line.split()[0]);h=[set() for _ in range(9)]
    for i,(u,v) in enumerate(P):
        if q>>(35-i)&1:h[u].add(v);h[v].add(u)
    # Diameter graph after adjoining the universal apex 9.
    g=[set(range(10))-{u}-(h[u] if u<9 else set()) for u in range(10)]
    forbidden=None
    for S in combinations(range(10),3):
        if not all(v in g[u] for u,v in combinations(S,2)):continue
        common=set.intersection(*(g[x] for x in S))
        E=[(u,v) for u,v in combinations(sorted(common),2) if v in g[u]]
        forbidden=next(((S,e,f) for e,f in combinations(E,2) if len(set(e+f))==4),None)
        if forbidden:break
    count+=1
    if forbidden:
        print(q,"triangle with double-edge common neighbors",forbidden)
        continue
    # Positive symbolic support of 2 Gram - J, including its diagonal.
    for a in range(9):h[a].add(a)
    V=tuple(combinations(range(9),6))
    def perms(A,B):
        signs=set()
        def dfs(i,M,sgn):
            if i==6:signs.add(sgn);return
            for j in range(6):
                if M>>j&1 and B[j] in h[A[i]]:
                    dfs(i+1,M^(1<<j),sgn^((M&((1<<j)-1)).bit_count()%2))
                    if len(signs)>1:return
        dfs(0,63,0);return signs
    cert=next(((A,B) for A in V for B in V if len(perms(A,B))==1),None)
    assert cert is not None, q
    print(q,"sign-nonsingular 6 by 6 minor",cert)
assert count==7
