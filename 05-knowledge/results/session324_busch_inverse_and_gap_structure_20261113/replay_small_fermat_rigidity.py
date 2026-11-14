"""Exact circulant-orbit profiles for the small Fermat rigidity obstruction."""
from itertools import product
from replay_gap_structure import ends, req
for p,wants in ((5,{(15,120)}),(7,{(175,1988),(189,2016)})):
    got=set()
    for choice in product((0,1),repeat=(p-1)//2):
        offsets=[p-j if rev else j for j,rev in enumerate(choice,1)]
        T=tuple(sum(1<<((i+j)%p) for j in offsets) for i in range(p))
        H=[1]+[sum(e) for e in ends(T)[1:]]
        got.add((H[-1],sum(a*b for a,b in zip(H,H[::-1]))))
    req(got==wants)
print("All five- and seven-orbit profiles verified")
