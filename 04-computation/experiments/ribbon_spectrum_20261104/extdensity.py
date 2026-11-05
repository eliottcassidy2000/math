#!/usr/bin/env python3
"""extdensity.py -- two-dimensional density of extensions in the 3-letter span-2 family:
for words u of length n, the pairs (H(u), H(ua)) per letter a; for each h in the solid interval,
the set of attainable H(ua) values: is it an interval of odd numbers?"""
import sys, itertools
from collections import defaultdict
from ribbon import transfer_system, letters, word_str
w = 2
order, states, mats, final = transfer_system(w)
F = frozenset
alph = [F([1]), F([2]), F([1, 2])]
finalset = set(final)
def apply(M, vd):
    nv = {}
    for (j, i), mult in M.items():
        if i in vd:
            nv[j] = nv.get(j, 0) + mult * vd[i]
    return nv
def Hof(vd):
    return sum(c for j, c in vd.items() if j in finalset)
n = int(sys.argv[1]) if len(sys.argv) > 1 else 10
cur = {((0, 1),)}
for k in range(n):
    nxt = set()
    for vec in cur:
        vd = dict(vec)
        for a in alph:
            nxt.add(tuple(sorted(apply(mats[a], vd).items())))
    cur = nxt
print(f"n={n}: {len(cur)} distinct state vectors")
ext = {a: defaultdict(set) for a in alph}
allext = defaultdict(set)
hs = set()
for vec in cur:
    vd = dict(vec); h = Hof(vd); hs.add(h)
    for a in alph:
        h2 = Hof(apply(mats[a], vd))
        ext[a][h].add(h2); allext[h].add(h2)
hs = sorted(hs)
print("H range at length n:", hs[0], hs[-1], "count", len(hs))
# for each h: extension set size, min, max, ratio range, and gaps
rows = []
for h in hs:
    s = sorted(allext[h])
    gaps = sum(1 for x in range(s[0], s[-1] + 1, 2) if x not in allext[h])
    rows.append((h, len(s), s[0], s[-1], gaps))
for r in rows[:: max(1, len(rows)//25)]:
    h, k, lo, hi, gaps = r
    print(f"  h={h:5d}: {k:3d} extension values in [{lo},{hi}] (ratios {lo/h:.3f}..{hi/h:.3f}), gaps {gaps}")
# global: union of all extension values vs the next length's solid interval
U = set()
for h in hs: U |= allext[h]
U = sorted(U)
miss = [x for x in range(U[0], U[-1] + 1, 2) if x not in set(U)]
print(f"union of extensions: [{U[0]},{U[-1]}], {len(U)} values, misses inside {len(miss)}: {miss[:30]}")
# per letter residual ranges relative to h
for a in alph:
    lo = min(min(v) - 2 * h for h, v in ext[a].items()); hi = max(max(v) - 2 * h for h, v in ext[a].items())
    print(f"  letter {sorted(a)}: H(ua) - 2H(u) in [{lo},{hi}]")
