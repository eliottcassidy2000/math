#!/usr/bin/env python3
"""
lrc14_thm735_three_slot_closure_kps_S128c3.py
=============================================
kind-pasteur-2026-07-13-S128 (cont.3).  THM-735(iii): the THREE-SLOT closure.

THEOREM: for ALL integers 10 < c < a < b, the 13-speed family {1..10, c, a, b} satisfies LRC(14).

Partition (each leg exact; sqrt2 replaced by the safe rational 99/70 > sqrt2):
  LEG J3 (c >= V1):        Bonferroni j=3 vs body {1..10}: Sum 1/v <= 3/V1 < 4 m0/(sqrt2 r0). One inequality.
  LEG J2 (c < V1, a >= V2(c)): Bonferroni j=2 vs body {1..10,c} (exact r_c, m_c): 1/a+1/b <= 2/V2 < 5 m_c/(sqrt2 r_c).
  LEG J1 (a < V2(c), b > v0(c,a)): THM-732(iii) tail vs body {1..10,c,a} (exact r_ca, m_ca).
  BOTTOM (b <= v0(c,a)): covering triple -> exact rational sweep of L; non-covering -> THM-366 (t=1/q).
Every eps in J3/J2 is measured against the FIXED body good set -- no carving, clustering-immune
(klein-S289's isolation wall never applies). MISTAKE-141: all thresholds exact.
"""
from fractions import Fraction as F
from functools import lru_cache
from math import isqrt
import time

ONE = F(1)
SQRT2_UP = F(99, 70)   # > sqrt2

@lru_cache(maxsize=8192)
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

def good_intervals_norm(speeds):
    """Good set as sorted linear pieces in [0,1] (wrap split), plus piece count and measure."""
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
    prev_hi = None
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
    """G (sorted disjoint pieces) minus D_w -> (count, measure, pieces)."""
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

def covering_1114(c, a, b):
    return all(c % q == 0 or a % q == 0 or b % q == 0 for q in (11, 12, 13, 14))

print("=" * 100)
print("THM-735(iii): {1..10,c,a,b} for ALL 10<c<a<b   (legs J3 / J2 / J1 / bottom, all exact)")
print("=" * 100)
t0 = time.time()
E0 = list(range(1, 11))
G0, r0, m0 = good_intervals_norm(E0)
# V1: minimal integer with 3/V1 < 4 m0/(SQRT2_UP r0)  (then Sum 1/v <= 3/c <= 3/V1 < threshold)
V1 = 1
while F(3, V1) >= 4 * m0 / (SQRT2_UP * r0):
    V1 += 1
print("body {1..10}: r0=%d m0=%s~%.5f   =>  LEG J3: all c>=V1=%d closed (j=3 Bonferroni, one inequality)"
      % (r0, m0, float(m0), V1))

tot_J2bodies = tot_pairs = tot_bottom = tot_sweeps = 0
tights = {}
covering_fails = []
maxV2 = (0, None)
maxv0 = (0, None)

for c in range(11, V1):
    r_c, m_c, G_c = subtract(G0, c)
    assert m_c > 0
    # V2(c): minimal integer with 2/V2 < 5 m_c/(SQRT2_UP r_c)
    V2 = c + 1
    while F(2, V2) >= 5 * m_c / (SQRT2_UP * r_c):
        V2 += 1
    tot_J2bodies += 1
    if V2 > maxV2[0]:
        maxV2 = (V2, c)
    for a in range(c + 1, V2):
        r_ca, m_ca, G_ca = subtract(G_c, a)
        assert m_ca > 0
        tot_pairs += 1
        # v0(c,a): tail fires for b > v0 ; v0 = SQRT2_UP*r_ca/(6*m_ca) rational
        v0 = SQRT2_UP * F(r_ca) / (6 * m_ca)
        bmax = v0.numerator // v0.denominator + 1
        while F(bmax) > v0:
            bmax -= 1
        if bmax > maxv0[0]:
            maxv0 = (bmax, (c, a))
        for b in range(a + 1, bmax + 1):
            tot_bottom += 1
            if not covering_1114(c, a, b):
                continue   # THM-366: some q in 11..14 has no multiple -> t=1/q witness
            tot_sweeps += 1
            _, L, _ = subtract(G_ca, b)
            if L == 0:
                fam = tuple(E0 + [c, a, b])
                mq = [q for q in range(2, 15) if not any(w % q == 0 for w in fam)]
                tights[fam] = mq
                if not mq:
                    covering_fails.append(fam)

print("LEG J2: %d bodies {1..10,c} (11<=c<%d), max V2=%d at c=%s" % (tot_J2bodies, V1, maxV2[0], maxV2[1]))
print("LEG J1: %d bodies {1..10,c,a}, max bottom-b=%d at (c,a)=%s" % (tot_pairs, maxv0[0], maxv0[1]))
print("BOTTOM: %d (c,a,b) triples below tail thresholds; %d covering (exact-swept); %d non-covering (THM-366)"
      % (tot_bottom, tot_sweeps, tot_bottom - tot_sweeps))
print("\nTIGHTS among swept covering triples (L=0): %s"
      % (", ".join("%s missing q=%s" % (f, q) for f, q in tights.items()) if tights else "NONE"))
print("COVERING L=0 (would block): %s" % (covering_fails if covering_fails else "NONE"))
print("\nelapsed %.1f s" % (time.time() - t0))
if not covering_fails:
    print("\nTHM-735(iii) ESTABLISHED: every {1..10,c,a,b} (10<c<a<b) satisfies LRC(14).")
    print("  J3: one inequality (c>=%d). J2: %d exact bodies. J1: %d exact bodies. Bottom: %d exact sweeps."
          % (V1, tot_J2bodies, tot_pairs, tot_sweeps))
else:
    print("\nBLOCKED -- investigate.")
