#!/usr/bin/env python3
"""Exact new-construction probes; no inherited-claim census."""
from itertools import product, combinations
from collections import Counter
import json

def group(d):
    if d == 0: return [(0,)]
    small=group(d-1); m=3**(d-1); ans=[]
    for r in range(3):
        for gs in product(small,repeat=3):
            ans.append(tuple(((i//m+r)%3)*m+gs[i//m][i%m] for i in range(3*m)))
    return ans

def digits(i,d):
    return tuple((i//3**j)%3 for j in reversed(range(d)))
def arrow(i,j,d):
    a,b=digits(i,d),digits(j,d)
    return next((y-x)%3==1 for x,y in zip(a,b) if x!=y)
def polynomial_next(p):
    cube=Counter()
    for a,ca in p.items():
        for b,cb in p.items():
            for c,cc in p.items(): cube[a+b+c]+=ca*cb*cc
    for a,ca in p.items(): cube[3*a]-=ca
    assert all(v%3==0 for v in cube.values())
    out=Counter({a:v//3 for a,v in cube.items() if v})
    for a,ca in p.items(): out[3*a+1]+=ca
    return out

result={}
for d in (1,2):
    n=3**d; G=group(d); pairs=list(combinations(range(n),2))
    assert len(G)==3**((n-1)//2)
    assert all(all(arrow(i,j,d)==arrow(g[i],g[j],d) for i,j in pairs) for g in G)
    losses={}
    for l in range(d):
        i,j=0,3**(d-l-1)
        stab=sum(g[i]==i and g[j]==j for g in G)
        assert len(G)//stab==3**(2*d-l-1)
        losses[l]={'arc_orbit':len(G)//stab,'surviving_group':stab,'block_flip_cost':9**(d-l-1)}
    coloring_stabs=Counter()
    for bits in product(range(2),repeat=n):
        s=sum(all(bits[g[i]]==bits[i] for i in range(n)) for g in G)
        coloring_stabs[s]+=1
    assert coloring_stabs[1]==2*len(G)
    # Orbit counts from all labelings, independently of recursive polynomial.
    orbit_stabs={str(s):c*s//len(G) for s,c in sorted(coloring_stabs.items())}
    result[d]={'order':n,'group_order':len(G),'arc_scales':losses,'binary_coloring_orbits_by_stabilizer':orbit_stabs}

p=Counter({0:2})
for d in range(1,6):
    p=polynomial_next(p)
    assert p[0]==2
    assert sum(c*3**((3**d-1)//2-a) for a,c in p.items())==2**(3**d)
    result.setdefault(d,{})['binary_type_count']=sum(p.values())
    result[d]['rigid_type_count']=p[0]
    if d<=2: assert {str(3**a):c for a,c in sorted(p.items())}==result[d]['binary_coloring_orbits_by_stabilizer']

# Hierarchical rectangle code: disjoint coordinates partition all leaf pairs.
d=2; n=3**d; rectangles={}
for i,j in combinations(range(n),2):
    a,b=digits(i,d),digits(j,d)
    l=next(t for t in range(d) if a[t]!=b[t])
    key=(a[:l],min(a[l],b[l]),max(a[l],b[l]))
    rectangles.setdefault(key,set()).add((i,j))
assert len(rectangles)==3*(n-1)//2
assert sorted(map(len,rectangles.values()))==[1]*9+[9]*3
assert sum(map(len,rectangles.values()))==n*(n-1)//2
result['code_n9']={'dimension':len(rectangles),'coordinate_weights':sorted(map(len,rectangles.values())),'ambient_arc_dimension':n*(n-1)//2}
print(json.dumps(result,indent=2))
