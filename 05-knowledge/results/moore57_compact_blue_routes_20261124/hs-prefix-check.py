#!/usr/bin/env python3
"""Replay a native all-blue path out of the compact 50-vertex HS chart."""
from pathlib import Path
from string import digits,ascii_letters
from compact import field_host,forward
T,p=field_host(5);alpha=dict(zip(p,p[::-1]))
n=len(T); pairs=Path(__file__).with_name('hs-monotone-prefix.txt').read_text().split()
def debt(T):
    return sum(abs(i-s) for i,s in enumerate(sorted(u.bit_count() for u in T)))//2
assert debt(T)==106 and len(pairs)==91
for i,t in enumerate(pairs,1):
    first=[(digits+ascii_letters).index(c) for c in t]
    p=first+[alpha[j] for j in first[::-1]]
    assert set(p)==set(range(n)) and len(p)==n
    assert all(forward(T,u,v) for u,v in zip(p,p[1:]))
    for u,v in zip(p,p[1:]):
        T[u]^=1<<v;T[v]^=1<<u
    assert debt(T)==106-i
# In literal native walk, use p or p[::-1] alternately as its orientation
# alternates between T and its converse. Each p is reflection-symmetric.
assert debt(T)==15
print("compact HS: 91 blue moves, score debt 106 -> 15; transitive not reached")
