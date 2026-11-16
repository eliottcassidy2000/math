#!/usr/bin/env python3
"""conjugacy_triples.py: (1) the Collatz conjugacy h = 2x + 1: T(x) = 3x+1 (x odd), x/2 (x even) becomes h -> 3h (h = 3 mod 4)
and h -> (h+1)/2 (h = 1 mod 4); orbit of x = 3 in count coordinates; (2) the sum graphs of the triples (5,12,13), (7,24,25),
(10,24,26): components of the graph on 1..m-1 with sums s = a^2, m = b^2, m + s = c^2 (FENCE_AND_SQUARE orbit decomposition:
one path on L-1 vertices, a second path on L if d even, floor((d-1)/2) cycles of length 2L, d = gcd(s, m), L = m/d)."""
def T(x): return 3 * x + 1 if x % 2 else x // 2
x = 3; orbit = [x]
while x != 1: x = T(x); orbit.append(x)
print("Collatz orbit of 3:", orbit)
print("count coordinates h = 2x + 1:", [2 * v + 1 for v in orbit], "(7 and 21 are the first two; all later ones are realized counts)")
def Th(h): return 3 * h if h % 4 == 3 else (h + 1) // 2
h = 7; horb = [h]
while h != 3: h = Th(h); horb.append(h)
print("conjugated map on h from 7:", horb)
from math import gcd
for (a, b, c) in ((3, 4, 5), (5, 12, 13), (7, 24, 25), (10, 24, 26), (24, 45, 51)):
    s, m = a * a, b * b; assert s + m == c * c
    d = gcd(s, m); L = m // d
    V = range(1, m); adj = {v: set() for v in V}
    for u in V:
        for t in (s, m, m + s):
            w = t - u
            if 1 <= w < m and w != u: adj[u].add(w)
    seen = set(); comps = []
    for v in V:
        if v in seen: continue
        st = [v]; seen.add(v); comp = []
        while st:
            y = st.pop(); comp.append(y)
            for z in adj[y]:
                if z not in seen: seen.add(z); st.append(z)
        comps.append(comp)
    kinds = sorted(('path' if min(len(adj[y]) for y in comp) <= 1 else 'cycle', len(comp)) for comp in comps)
    ends = sorted(y for comp in comps for y in comp if len(adj[y]) == 1)
    print(f"triple {(a, b, c)}: s={s} m={m} d={d} L={L}: pieces {kinds}; predicted path(L-1)={L-1}" + (f", path(L)={L}" if d % 2 == 0 else "") + f", {(d-1)//2} cycles of {2*L}; degree-1 vertices {ends[:8]}")
