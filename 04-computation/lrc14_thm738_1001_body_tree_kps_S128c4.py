#!/usr/bin/env python3
"""
lrc14_thm738_1001_body_tree_kps_S128c4.py
=========================================
kind-pasteur-2026-07-13-S128 (cont.4).  THM-738: the near-AP THREE-SLOT closure.

THEOREM: every 13-speed family with AT LEAST 10 speeds in {1..14} satisfies LRC(14).
Equivalently: for every 10-element body E subset {1..14} (all C(14,10)=1001) and all c<a<b not in E,
{E,c,a,b} is lonely.

Per body (THM-735's Bonferroni tree, Q-general):
  Q(E) = {q in 2..14 : no multiple of q in E}    [general bodies may miss SMALL q too]
  LEG J3: all c,a,b >= V1(E)  [j=3 Bonferroni vs G_E: 3/V1 < 4 m_E/((99/70) r_E)]
  LEG J2: c < V1, both a,b >= V2(c)  [j=2 vs exact G_{E,c}]
  LEG J1: a < V2(c), b > v0(c,a)  [THM-732(iii) tail vs exact G_{E,c,a}]
  BOTTOM: b <= v0(c,a): the family is covering iff b covers Qb = {q in Q(E): q ndiv c, q ndiv a};
          enumerate b as multiples of lcm(Qb) (all of Qb must divide b) -- or ALL b when Qb is empty;
          exact-Q sweep each; non-covering b -> THM-366 (t=1/q).
Any L=0 among swept (necessarily covering) triples would BLOCK the theorem -> flagged loudly.
REGRESSION: body {1..10} must reproduce THM-735(iii): V1=154, 143 J2 bodies, 7537 J1 pairs, 27 sweeps.
"""
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations
from math import gcd
import time

ONE = F(1)
SQRT2_UP = F(99, 70)

@lru_cache(maxsize=16384)
def bad_pieces(w):
    r = F(1, 14 * w)
    out = []
    for k in range(w):
        c = F(k, w)
        lo, hi = c - r, c + r
        if lo < 0:
            out.append((F(0), hi)); out.append((lo + 1, ONE))
        elif hi > 1:
            out.append((lo, ONE)); out.append((F(0), hi - 1))
        else:
            out.append((lo, hi))
    return tuple(sorted(out))

def good_norm(speeds):
    pieces = []
    for w in speeds:
        pieces += bad_pieces(w)
    pieces = sorted(pieces)
    comps = []
    for lo, hi in pieces:
        if comps and lo <= comps[-1][1]:
            if hi > comps[-1][1]:
                comps[-1][1] = hi
        else:
            comps.append([lo, hi])
    out = []
    for i in range(len(comps)):
        a = comps[i][1]
        j = (i + 1) % len(comps)
        b = comps[j][0] + (ONE if j == 0 else 0)
        if a < b:
            if b <= 1:
                out.append((a, b))
            else:
                out.append((a, ONE))
                if b - 1 > 0:
                    out.append((F(0), b - 1))
    out.sort()
    return out, len(out), sum(b - a for a, b in out)

def subtract(G, w):
    B = bad_pieces(w)
    out = []
    j = 0
    for a, b in G:
        cur = a
        while j > 0 and B[j - 1][1] > cur:
            j -= 1
        k = j
        while k < len(B) and B[k][0] < b:
            lo, hi = B[k]
            if hi <= cur:
                k += 1; continue
            if lo > cur:
                out.append((cur, min(lo, b)))
            cur = max(cur, hi)
            if cur >= b:
                break
            k += 1
        if cur < b:
            out.append((cur, b))
        j = k
    return len(out), sum(b - a for a, b in out), out

def lcm(xs):
    L = 1
    for x in xs:
        L = L * x // gcd(L, x)
    return L

print("=" * 102)
print("THM-738: all 1001 ten-element bodies E in {1..14}; {E,c,a,b} for all c<a<b not in E  (j=3 tree)")
print("=" * 102)
t0 = time.time()
bodies = list(combinations(range(1, 15), 10))
tot_J1 = tot_sweeps = tot_366 = 0
covering_fails = []
tights = {}
maxV1 = (0, None)
worst_body = (0.0, None)
regression = None

for bi, E in enumerate(bodies):
    tb = time.time()
    Eset = set(E)
    QE = [q for q in range(2, 15) if not any(w % q == 0 for w in E)]
    GE, rE, mE = good_norm(E)
    assert mE > 0
    V1 = 1
    while F(3, V1) >= 4 * mE / (SQRT2_UP * rE):
        V1 += 1
    if V1 > maxV1[0]:
        maxV1 = (V1, E)
    nJ2 = nJ1 = nsw = n366 = 0
    for c in range(1, V1):
        if c in Eset:
            continue
        r_c, m_c, G_c = subtract(GE, c)
        assert m_c > 0
        nJ2 += 1
        V2 = c + 1
        while F(2, V2) >= 5 * m_c / (SQRT2_UP * r_c):
            V2 += 1
        for a in range(c + 1, V2):
            if a in Eset:
                continue
            r_ca, m_ca, G_ca = subtract(G_c, a)
            assert m_ca > 0
            nJ1 += 1
            v0 = SQRT2_UP * F(r_ca) / (6 * m_ca)
            bmax = v0.numerator // v0.denominator + 1
            while F(bmax) > v0:
                bmax -= 1
            if bmax <= a:
                continue
            Qb = [q for q in QE if c % q != 0 and a % q != 0]
            if Qb:
                L0 = lcm(Qb)
                bs = range(((a // L0) + 1) * L0, bmax + 1, L0)
                n366 += (bmax - a) - len([b for b in bs])   # rough count of THM-366-dispatched b (display only)
            else:
                bs = range(a + 1, bmax + 1)
            for b in bs:
                if b in Eset:
                    continue
                nsw += 1
                _, L, _ = subtract(G_ca, b)
                if L == 0:
                    fam = tuple(sorted(E + (c, a, b)))
                    mq = [q for q in range(2, 15) if not any(w % q == 0 for w in fam)]
                    tights[fam] = mq
                    if not mq:
                        covering_fails.append(fam)
    tot_J1 += nJ1; tot_sweeps += nsw; tot_366 += n366
    dt = time.time() - tb
    if dt > worst_body[0]:
        worst_body = (dt, E)
    if E == tuple(range(1, 11)):
        regression = (V1, nJ2, nJ1, nsw)
    if bi % 77 == 0 or dt > 20:
        print("  body %4d/1001 %-38s Q=%-12s r=%2d m=%.4f V1=%4d J1-pairs=%6d sweeps=%5d  [%.1fs tot %.0fs]"
              % (bi + 1, "{" + ",".join(map(str, E)) + "}", QE, rE, float(mE), V1, nJ1, nsw, dt, time.time() - t0))

print("\n" + "=" * 102)
print("GLOBAL (%.0f s): 1001 bodies; J1 pairs total %d ; bottom exact sweeps %d ; max V1=%d at %s"
      % (time.time() - t0, tot_J1, tot_sweeps, maxV1[0], maxV1[1]))
print("worst body time %.1fs at %s" % worst_body)
print("REGRESSION body {1..10}: V1,J2,J1,sweeps = %s   (THM-735(iii) expects (154, 143, 7537, 27))" % (regression,))
assert regression == (154, 143, 7537, 27), "regression mismatch vs THM-735(iii)!"
print("\nTIGHTS among swept covering triples: %s"
      % (", ".join("%s missing q=%s" % (f, q) for f, q in tights.items()) if tights else "NONE"))
print("COVERING L=0 (would block): %s" % (covering_fails if covering_fails else "NONE"))
if not covering_fails:
    print("\nTHM-738 ESTABLISHED: every 13-speed family with >=10 speeds in {1..14} satisfies LRC(14).")
else:
    print("\nBLOCKED -- investigate.")
