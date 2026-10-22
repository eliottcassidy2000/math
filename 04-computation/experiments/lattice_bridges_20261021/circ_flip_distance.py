# Circulant tournaments Cay(Z_n, S): H against the flip distance to the rotational tournament R_n = {1..m}.
# Distance d(S) = min over units a of #{k in [1,m] : n-k in aS}  (number of difference classes to reverse).
# Data: the complete multiplier-class tables of circ_counts.c (tournament_pseudoprimes_20261021/circ_*.out).
# mac-mini-2026-10-21-fanopsp.
import re, glob
from math import gcd
from collections import defaultdict
rows = defaultdict(list)
for f in sorted(glob.glob('../tournament_pseudoprimes_20261021/circ_*.out')):
    for line in open(f):
        m_ = re.match(r'n=(\d+) S=\{([\d,]+)\} mult=(\d+) H=(\d+)', line)
        if m_:
            n = int(m_.group(1)); S = frozenset(int(x) for x in m_.group(2).split(',')); H = int(m_.group(4))
            rows[n].append((S, H))
def dist(n, S):
    m = (n - 1) // 2; best = m
    for a in range(1, n):
        if gcd(a, n) != 1: continue
        aS = {(a * s) % n for s in S}
        best = min(best, sum(1 for k in range(1, m + 1) if (n - k) in aS))
    return best
for n in sorted(rows):
    data = sorted(((H, S) for S, H in rows[n]), key=lambda t: -t[0])
    ds = [(H, dist(n, S), S) for H, S in data]
    Hmax, dmax, Smax = ds[0]
    bydist = defaultdict(list)
    for H, d, S in ds: bydist[d].append(H)
    mono = all(max(bydist[d]) > max(bydist[d + 1]) for d in sorted(bydist) if d + 1 in bydist)
    print(f"n={n:2d}: classes={len(ds):3d}  H-max class S={sorted(Smax)} at flip distance {dmax} from R_n;  "
          f"max H by distance: " + ", ".join(f"d={d}:{max(bydist[d])}" for d in sorted(bydist)) + f"  monotone-decreasing={mono}")
