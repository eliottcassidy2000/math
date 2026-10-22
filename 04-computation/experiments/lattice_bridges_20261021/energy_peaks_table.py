# For each odd n with a complete circulant class table (n <= 25): H-rank of the energy maximum (R_n) and of the energy
# minimum, and whether each is a strict local maximum of H on the flip graph.  mac-mini-2026-10-21-fanopsp.
import re, glob
from math import gcd
from collections import defaultdict
rows = defaultdict(dict)
for f in sorted(glob.glob('../tournament_pseudoprimes_20261021/circ_*.out')):
    for line in open(f):
        mm = re.match(r'n=(\d+) S=\{([\d,]+)\} mult=(\d+) H=(\d+)', line)
        if mm: rows[int(mm.group(1))][tuple(sorted(int(x) for x in mm.group(2).split(',')))] = int(mm.group(4))
for line in open('circ25_H.txt'):
    mm = re.match(r'S=\{([\d,]+)\} H=(\d+)', line)
    if mm: rows[25][tuple(int(x) for x in mm.group(1).split(','))] = int(mm.group(2))
def canon(n, S): return min(tuple(sorted((a * s) % n for s in S)) for a in range(1, n) if gcd(a, n) == 1)
for n in sorted(rows):
    if n < 7: continue
    m = (n - 1) // 2
    H = {canon(n, S): h for S, h in rows[n].items()}
    def E(S):
        Ss = set(S); return sum(1 for a in S for b in S for c in S if (a + b - c) % n in Ss)
    def flip(S, k):
        S2 = set(S)
        if k in S2: S2.discard(k); S2.add(n - k)
        else: S2.discard(n - k); S2.add(k)
        return canon(n, S2)
    def is_loc(S): return all(H[S] > H[flip(S, k)] for k in range(1, m + 1))
    order = sorted(H, key=lambda S: -H[S]); rank = {S: i + 1 for i, S in enumerate(order)}
    Ev = {S: E(S) for S in H}; emin = min(Ev.values()); emax = max(Ev.values())
    mins = [S for S in H if Ev[S] == emin]; maxs = [S for S in H if Ev[S] == emax]
    locs = [S for S in H if is_loc(S)]
    desc_min = "; ".join(f"H-rank {rank[S]}, local max {is_loc(S)}" for S in mins)
    print(f"n={n:2d} ({len(H):3d} classes, {len(locs)} strict local maxima): E-max H-rank {rank[maxs[0]]} (local max {is_loc(maxs[0])}) | "
          f"E-min x{len(mins)}: {desc_min}")
