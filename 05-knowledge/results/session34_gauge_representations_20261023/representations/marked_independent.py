#!/usr/bin/env python3
"""Direct labelled (tournament,path) census, independent of chord normalization."""
from itertools import combinations,permutations
from math import factorial
from collections import Counter
import json
rows=[]
for n in range(1,5):
    pairs=list(combinations(range(n),2));ix={e:k for k,e in enumerate(pairs)};ps=list(permutations(range(n)))
    def edge(mask,a,b):return bool((mask>>ix[tuple(sorted((a,b)))])&1)==(a<b)
    Y=[(mask,p) for mask in range(1<<len(pairs)) for p in ps if all(edge(mask,a,b) for a,b in zip(p,p[1:]))]
    m=(n-1)*(n-2)//2;q=(m+(n-1)//2)//2
    counts=[]
    for g in ps:
        f=sum(tuple(g[v] for v in p[::-1])==p and all(edge(mask,g[a],g[b])!=edge(mask,a,b) for a,b in pairs) for mask,p in Y)
        # g must be conjugate to reversal, hence involutive with n mod 2 fixed vertices.
        typematches=all(g[g[i]]==i for i in range(n)) and sum(g[i]==i for i in range(n))==n%2
        expected=(2**(n//2)*factorial(n//2)*2**q) if typematches else 0
        assert f==expected
        counts.append(f)
    assert len(Y)==factorial(n)*2**m
    rows.append({'n':n,'marked_objects':len(Y),'twisted_fixed_count_histogram':dict(Counter(counts))})
p=7;Q={i*i%p for i in range(1,p)}
paths=[v for v in permutations(range(p)) if all((b-a)%p in Q for a,b in zip(v,v[1:]))]
fixed=sum(tuple((-x)%p for x in v[::-1])==v for v in paths)
assert len(paths)==189 and fixed==9
anti=[tuple((a*x+b)%7 for x in range(7)) for a in range(1,7) if a not in Q for b in range(7)]
fixed_total=sum(tuple(g[x] for x in v[::-1])==v for g in anti for v in paths)
assert fixed_total==63
assert fixed_total//21==3
print(json.dumps({'marked_census':rows,'paley7':{'paths':len(paths),'fixed_by_negative_reverse':fixed,'automorphism_orbits':9,'converse_trace_on_aut_orbits':3,'converse_plus_multiplicity':6,'converse_minus_multiplicity':3}},indent=2))
