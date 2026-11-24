#!/usr/bin/env python3
"""Direct graph spines, exact chord footprints, and block-contiguous Hamilton spines."""
from collections import Counter
from itertools import combinations, permutations
from math import comb

def graph(d):
    if d==2:
        return [{(i-1)%5,(i+1)%5} for i in range(5)]
    if d==3:
        a=list(map(set,combinations(range(5),2)))
        return [{j for j,y in enumerate(a) if not x&y} for x in a]
    V=[(z,i,x) for z in (0,1) for i in range(5) for x in range(5)]
    def edge(u,v):
        z,i,x=u;w,j,y=v
        if z==w:return i==j and (y-x)%5 in (1+z,4-z)
        if z:return edge(v,u)
        return (y-x-i*j)%5==0
    return [{j for j,w in enumerate(V) if edge(v,w)} for v in V]

def field_spine():
    P=[]; a=0
    for i in range(5):
        P.extend(5*i+(a+j)%5 for j in range(5))
        b=(a-1+i*i)%5
        P.extend(25+5*i+(b+2*j)%5 for j in range(5))
        a=(a-3-i)%5
    assert a==0
    return P

def tiles(G,p):
    n=len(G)
    assert len(set(p))==n and all(v in G[u] for u,v in zip(p,p[1:]))
    return {(i,j) for i in range(n) for j in range(i+2,n) if p[j] in G[p[i]]}

def check_packing(G,p,d):
    n=len(G); C=tiles(G,p)
    assert len(C)==n*d//2-n+1 and all(j-i>=4 for i,j in C)
    plus=set();wedge=set(); X=[set() for _ in G]
    for i,j in C:
        X[i].add(j);X[j].add(i)
        for a,b in [(i,j),(i-1,j),(i+1,j),(i,j-1),(i,j+1)]:
            if 0<=a<b<n:
                assert b-a>=3 and (a,b) not in plus
                plus.add((a,b))
    for R in X:
        for i,j in combinations(sorted(R),2):
            assert j-i>=3 and (i,j) not in plus|wedge
            wedge.add((i,j))
    assert len(plus)==5*len(C)-2*(d-1)
    assert len(wedge)==2*comb(d-1,2)+(n-2)*comb(d-2,2)
    assert plus|wedge=={(i,j) for i in range(n) for j in range(i+3,n)}
    reflect={(n-1-j,n-1-i) for i,j in C}
    print(f"d={d}: chord {len(C)}; plus cells {len(plus)}, wedge cells {len(wedge)}; reflection defect {len(C^reflect)}")
    # Tournament triangle receipts for the native complementary pair.
    out=[]
    for mode in (0,1):
        scores=[0]*n
        for i,j in combinations(range(n),2):
            rev=j>i+1 and (((i,j) in C)^bool(mode))
            scores[j if rev else i]+=1
        out.append(comb(n,3)-sum(comb(s,2) for s in scores))
    assert out[1]-out[0]==(d-1)**2
    print('directed triangles sparse, dense:',tuple(out))
    return C

for d,p in [(2,list(range(5))),(3,[0,7,3,4,9,1,6,2,8,5]),(7,field_spine())]:
    G=graph(d)
    assert all(len(G[i])==d for i in range(len(G)))
    assert all(len(G[i]&G[j])==(j not in G[i]) for i,j in combinations(range(len(G)),2))
    C=check_packing(G,p,d)
    if d!=3: assert p[0] in G[p[-1]]
    if d==7:assert len(C^{(49-j,49-i) for i,j in C})==200

def petersen_spines():
    G=graph(3); tab=Counter()
    def rec(P,B):
        if not B:
            assert P[-1] not in G[P[0]] # no graph Hamilton cycle in Petersen
            C=tiles(G,P)
            keys=[(i,j) for i in range(10) for j in range(i+2,10)]
            tab[sum(1<<k for k,a in enumerate(keys) if a in C)]+=1
            return
        for x in B&G[P[-1]]:
            rec(P+[x],B-{x})
    for x in range(10):rec([x],set(range(10))-{x})
    assert tab=={4304406596:120,4311892004:120}
    print('Petersen oriented graph spines and tile masks',dict(tab))

petersen_spines()
S=list(permutations(range(5)))
area=Counter(sum((u[i]-u[(i+1)%5])*v[i] for i in range(5))%5 for u in S for v in S)
eps=Counter({0:1})
for a in (1,)*5+(2,)*5:
    d=Counter()
    for t,y in eps.items():
        d[(t-a)%5]+=y;d[(t+a)%5]+=y
    eps=d
assert eps=={0:204,1:205,2:205,3:205,4:205}
ordered_U_first=5*sum(area[i]*eps[i] for i in range(5))
assert ordered_U_first%10==0
print('HS affine area distribution:',sorted(area.items()))
print('HS alternating whole-cohort spines:',2*len(S)**2*5*2**10,'oriented paths;',ordered_U_first//10,'unoriented cycles')
