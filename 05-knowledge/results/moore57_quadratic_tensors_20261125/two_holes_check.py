"""Physical triple obstruction, its pair projections, and the label-free repair board."""
import json
import sys
from itertools import combinations, product
from collections import Counter
from pathlib import Path
import nested_core_check as v
from three_stars import SLOTS

C,S=v.orbit_chart('period_fourteen_site.txt',4,((10,45,2),(27,49,0)))
R=dict(json.loads(Path('rigid_51_contacts.json').read_text()))
b,x=1,1
centers=(56,0,28)
T=lambda h,y:tuple(C[h,c][y] for c in centers)
A=[set(range(56))-({C[b,c][x]}|{C[d,c][C[b,d][x]] for d in centers if d!=c}) for c in centers]
V=[a-{T(h,y)[i] for h,y in R.items()} for i,a in enumerate(A)]
assert V==[{3,16},{43,49},{24,45}]
assert set(SLOTS)-{b}-set(R)=={13,33}
slots=sorted(A[0])
used={h:{y for E in S.values() for k,y in E if k==h} for h in SLOTS}
assert used[13]==used[33]==set(range(2,56,4))
old={y:56*h+y for h,y in R.items()}


def options(y):
    for h in SLOTS:
        if h==b or y in used[h] or (y==7 and h!=2):
            continue
        t=T(h,y)
        if all(t[i] in A[i] for i in range(3)):
            yield [56*h+y,*t]


if '--board' in sys.argv:
    print(len(slots))
    for y in slots:
        D=list(options(y))
        print(old.get(y,-1),len(D))
        for a in D:print(*a)
else:
    all_trips={T(h,y):(h,y) for h in SLOTS for y in range(56)}
    assert len(all_trips)==3024
    assert not set(product(*V))&all_trips.keys()
    for i,j in combinations(range(3),2):
        proj={(t[i],t[j]) for t in all_trips}
        assert set(product(V[i],V[j]))<=proj
    assert [y for y in range(56) if sum(t in z for t,z in zip(T(13,y),V))>=2]==[16,29]
    assert not any(sum(t in z for t,z in zip(T(33,y),V))>=2 for y in range(56))
    sol=list(map(int,Path('relaxed_witness.txt').read_text().split()))
    assert len(sol)==53 and len(set(old.values())-set(sol))==7
    tt=[T(t//56,t%56) for t in sol]
    assert all({t[i] for t in tt}==A[i] for i in range(3))
    assert all(any(t==a[0] for a in options(x)) for x,t in zip(slots,sol))
    assert len({t//56 for t in sol})==46
    # Square-sum controls, independent of the chart addresses.
    from math import isqrt
    sq=lambda n:isqrt(n)**2==n
    assert not sq(13+33)
    assert [x for x in range(1,10001) if sq(x+13) and sq(x+33)]==[3]
    print("Eight absent triples and twelve present projected pairs; seven rewires in witness, seven repeated branch uses.")
