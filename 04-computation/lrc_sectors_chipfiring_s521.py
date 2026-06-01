#!/usr/bin/env python3
"""
lrc_sectors_chipfiring_s521.py   claudebox-2026-06-01-S521

"Sectors as nodes" encoding of LRC (reflection:
07-reflections/lrc-sectors-as-nodes-token-chipfiring-s521.md).

Nodes = n fixed sectors S_k=[k/n,(k+1)/n). Runners are m tokens with occupancy
o_k=#runners in S_k. strict-LRC <=> some t has o_0=o_{n-1}=0 (clear observer's two
forbidden sectors). A boundary crossing v_i t=k/n moves one token to the adjacent
sector: the dynamics is m tokens CHIP-FIRING on the n-cycle. At t=a/n the
occupancy equals the residue distribution {v_i a mod n} (THM-369 picture).
Captures STRICT loneliness; tight extremizers are lonely only at the boundary.
"""
from fractions import Fraction as F
from math import gcd, comb
import random

def fr(x): return x % 1
def sector(p, n): return int(fr(p) * n)
def cells_w(sp, n):
    W = set([F(0)])
    for v in sp:
        for k in range(n*v+1): W.add(F(k, n*v) % 1)
    W = sorted(w for w in W if 0 <= w < 1); W2 = W + [F(1)]
    return W + [(a+b)/2 for a, b in zip(W, W2[1:])]
def occ(sp, t, n):
    o = [0]*n
    for v in sp: o[sector(F(v)*t, n)] += 1
    return tuple(o)
def adjacent_move(a, b, n):
    d = [b[k]-a[k] for k in range(n)]
    if sum(abs(x) for x in d) != 2: return False
    plus = [k for k in range(n) if d[k] == 1]; minus = [k for k in range(n) if d[k] == -1]
    return len(plus) == 1 and len(minus) == 1 and (plus[0]-minus[0]) % n in (1, n-1)

def main():
    print("Sector-occupancy encoding (strict loneliness = clear sectors {0,n-1}):\n")
    print(f"{'speeds':18} {'realizable occ':14} {'#compositions':13} {'strict-lonely?':14}")
    for sp in [(1,2,3,4),(1,2,4,7),(1,3,4,5,9),(2,3,5,7,11)]:
        n = len(sp)+1; m = len(sp)
        R = set(occ(list(sp), t, n) for t in cells_w(list(sp), n))
        lonely = any(o[0] == 0 and o[n-1] == 0 for o in R)
        print(f"{str(sp):18} {len(R):14} {comb(m+n-1,n-1):13} {str(lonely):14}")
    print("\nResidue identity at t=a/n: occupancy sector of runner i = v_i*a mod n.")
    for sp in [(1,2,3,4),(1,2,4,7)]:
        n = len(sp)+1
        hits = [a for a in range(1, n) if all((v*a) % n not in (0, n-1) for v in sp)]
        print(f"  {sp} n={n}: a with strict-lonely at a/n: {hits}  (empty for tight sets)")
    print("\nChip-firing dynamics: fraction of occupancy-changes that are single adjacent-token hops:")
    random.seed(0)
    for sp in [(1,2,4,7),(2,3,5,7,11)]:
        n = len(sp)+1
        seq = []
        for t in sorted(set(cells_w(list(sp), n))):
            o = occ(list(sp), t, n)
            if not seq or o != seq[-1]: seq.append(o)
        adj = sum(1 for a, b in zip(seq, seq[1:]) if adjacent_move(a, b, n))
        print(f"  {sp}: {adj}/{len(seq)-1} single adjacent-token hops (rest = simultaneous resonant crossings)")

if __name__ == "__main__":
    main()
