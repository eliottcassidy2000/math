#!/usr/bin/env python3
"""
lrc_encoding_taxonomy_s521.py   claudebox-2026-06-01-S521

Survey of LRC-to-tournament encodings by FAITHFULNESS x RESTRICTION.
(Reflection: 07-reflections/lrc-encoding-taxonomy-faithfulness-restriction-s521.md)

For speeds v (observer 0), n=m+1, threshold 1/n.  Each encoding maps (v,t) to an
iso-class; the exhibitable set is over t in [0,1) (CLOSED cells incl. walls, since
tight families are lonely only at the boundary t=1/n).  Encodings:

  E0  round tournament (runner half-turn)  -- UNFAITHFUL (observer-blind) = A000016(m)
  E5  two-gap {LL,LS,SL,SS}                -- minimal faithful; LRC <=> LL exhibitable
  EB  gap-type necklace (tight/loose gaps) -- #tight<=n-1; LRC <=> both obs-gaps loose
  E1  finest faithful (round+tight+observer rooted)
  E12 residue-orbit at n  -- UNFAITHFUL (single-modulus): q=n loneliness only

Plus two structural checks: gap-necklace #tight<=n-1; residue-orbit unfaithfulness.
"""
from fractions import Fraction as F
from math import gcd
from itertools import permutations
import random

def fr(x): return x % 1
def dist(x):
    x = x % 1; return min(x, 1 - x)

def walls(sp, n):
    W = set([F(0)])
    for v in sp:
        for k in range(2*v+1): W.add(F(k, 2*v))
        for k in range(v+1):
            W.add(fr(F(k, v) + F(1, n*v))); W.add(fr(F(k, v) - F(1, n*v)))
    for i in range(len(sp)):
        for j in range(len(sp)):
            d = abs(sp[i]-sp[j])
            if d:
                for k in range(2*d+1): W.add(F(k, 2*d))
    return sorted(w for w in W if 0 <= w < 1)

def times(sp, n):
    W = walls(sp, n); W2 = W + [F(1)]
    return W + [(a+b)/2 for a, b in zip(W, W2[1:])]

def canon_round(sp, t):
    m = len(sp); adj = [[1 if (a != b and 0 < fr((sp[a]-sp[b])*t) < F(1,2)) else 0
                         for b in range(m)] for a in range(m)]
    best = None
    for p in permutations(range(m)):
        f = tuple(adj[p[i]][p[j]] for i in range(m) for j in range(m) if i != j)
        if best is None or f < best: best = f
    return best

def gap_necklace(sp, t, n):
    pos = sorted([F(0)] + [fr(F(v)*t) for v in sp]); k = len(pos)
    oi = pos.index(F(0))
    g = [(pos[(i+1) % k]-pos[i]) % 1 for i in range(k)]
    tl = tuple('T' if g[(oi+i) % k] < F(1, n) else 'L' for i in range(k))  # from observer
    return tl

def two_gap(sp, t, n):
    pos = [fr(F(v)*t) for v in sp]
    cw = min((p for p in pos if p > 0), default=F(1))
    ccw = min((1-p for p in pos if p > 0), default=F(1))
    return ('L' if cw >= F(1, n) else 'S') + ('L' if ccw >= F(1, n) else 'S')

def survey(sp):
    n = len(sp)+1
    R0, R5, RB = set(), set(), set()
    for t in times(list(sp), n):
        R0.add(canon_round(list(sp), t))
        R5.add(two_gap(list(sp), t, n))
        RB.add(gap_necklace(list(sp), t, n))
    return len(R0), len(R5), len(RB), max(s.count('T') for s in RB)

def lonely_some_t(sp, n):
    return any(all(dist(F(v)*t) >= F(1, n) for v in sp) for t in times(list(sp), n))
def lonely_at_qn(sp, n):
    return any(all((v*a) % n != 0 for v in sp) for a in range(1, n))

def main():
    print("Encoding survey (closed cells): E0 round (obs-blind), E5 two-gap, EB gap-necklace\n")
    print(f"{'speeds':22} {'E0':>4} {'E5':>4} {'EB':>4} {'max#tight':>9} {'n-1':>4}")
    for sp in [(1,2,3,4),(1,2,4,7),(1,3,4,5,9),(2,3,5,7,11),(1,5,6,11,16)]:
        n = len(sp)+1; r0, r5, rb, mt = survey(sp)
        print(f"{str(sp):22} {r0:>4} {r5:>4} {rb:>4} {mt:>9} {n-1:>4}")
    print("\nStructural law: max #tight gaps <= n-1 (all-tight impossible: n gaps <1/n can't sum to 1).")
    print("\nResidue-orbit-at-n UNFAITHFULNESS (lonely overall but NOT at q=n):")
    random.seed(0)
    for n in (6, 8, 9):
        m = n-1; unf = tested = 0
        for _ in range(2000):
            sp = sorted(random.sample(range(1, 3*n), m))
            g = 0
            for v in sp: g = gcd(g, v)
            if g != 1: continue
            tested += 1
            if lonely_some_t(sp, n) and not lonely_at_qn(sp, n): unf += 1
        print(f"  n={n}: {unf}/{tested} sets lonely only away from q=n -> residue-orbit misses them")

if __name__ == "__main__":
    main()
