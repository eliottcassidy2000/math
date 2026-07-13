#!/usr/bin/env python3
"""
lrc14_thm734_body_sweep_kps_S128c2.py
=====================================
kind-pasteur-2026-07-13-S128 (cont.2).  THM-734: the near-AP two-slot closure.

For EVERY 11-element body E subset {1..14} (all C(14,11)=364) and every pair a<b of positive
integers not in E, the 13-speed family {E,a,b} satisfies LRC(14).

Per body (THM-733's method, body-parametric):
  LEG 1 (uniform, a >= A0(E)): peel lemmas P1 (r_a <= a m + (15/7) r) and P2
        (m_a >= (6/7) m - 8r/(49a)) + THM-732(iii) tail. f(a) = K M(a)(a+1) - R(a) with
        K = 2121/500 < sqrt(18) is EXACTLY increasing (slope m(6K/7-1) > 0 plus an increasing
        term), so f(A0) > 0 certifies all a >= A0: b >= a+1 > v0(base) => certificate => L > 0.
  LEG 2 (per-a, a < A0, a not in E): exact (r_a, m_a) of E u {a}; every b with
        b^2 > v0sq = r_a^2/(18 m_a^2) closes by the tail.
  LEG 3 (box): remaining pairs a < b <= sqrt(v0sq), b not in E: exact rational sweep of
        L({E,a,b}). L > 0 => proved. L = 0 => must be NON-covering (some q in 2..14 with no
        multiple) => THM-366 (t=1/q); a COVERING L=0 pair would block the theorem (flag).

m_E > 0 and m_a > 0 are automatic (LRC(<=13) settled: <=12 speeds have max-min >= 1/13 > 1/14,
giving an open neighborhood of good times) — asserted anyway.

Output: per-body line (A0, box size, tights); global tight census (deduped families) + verdict.
All arithmetic exact (fractions.Fraction); floats display-only.
"""
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations
from math import isqrt
import time

ONE = F(1)
K = F(2121, 500)          # K^2 = 4498641/250000 < 18
assert K * K < 18

@lru_cache(maxsize=4096)
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
    return tuple(out)

def good_intervals(speeds):
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
    n = len(comps)
    goods = []
    for i in range(n):
        a = comps[i][1]
        j = (i + 1) % n
        b = comps[j][0] + (ONE if j == 0 else 0)
        if a < b:
            goods.append((a, b))
    meas = sum(b - a for a, b in goods)
    return goods, len(goods), meas

def normalize(goods):
    out = []
    for a, b in goods:
        if b <= 1:
            out.append((a, b))
        else:
            out.append((a, ONE)); out.append((F(0), b - 1))
    out.sort()
    return out

def subtract_arcs(G, w):
    """G (sorted disjoint pieces in [0,1]) minus D_w -> (piece_count, measure, pieces)."""
    B = bad_pieces(w)   # already sorted by construction? ensure:
    B = sorted(B)
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
    m = sum(b - a for a, b in out)
    return len(out), m, out

def find_A0(r, m):
    """Minimal a with f(a) > 0; f exactly increasing (proved)."""
    lo, hi = 1, 2
    def f(a):
        return K * (F(6, 7) * m - F(8 * r, 49 * a)) * (a + 1) - (a * m + F(15, 7) * r)
    while f(hi) <= 0:
        hi *= 2
        assert hi < 10**9
    while lo + 1 < hi:
        mid = (lo + hi) // 2
        if f(mid) > 0:
            hi = mid
        else:
            lo = mid
    assert f(hi) > 0 and (hi == 1 or f(hi - 1) <= 0)
    return hi

def missing_q(S):
    return [q for q in range(2, 15) if not any(w % q == 0 for w in S)]

print("=" * 100)
print("THM-734 BODY-BY-BODY SWEEP: all C(14,11)=364 bodies E in {1..14}; close {E,a,b} for all a<b not in E")
print("=" * 100)

t0 = time.time()
bodies = list(combinations(range(1, 15), 11))
tights = {}          # frozenset(family) -> missing q list
covering_fails = []
tot_box = 0
maxA0 = (0, None)
minm = (F(10), None)
worst_box = (0, None)

for bi, E in enumerate(bodies):
    goodsE, rE, mE = good_intervals(E)
    assert mE > 0, "body %s has m=0" % (E,)
    if mE < minm[0]:
        minm = (mE, E)
    GE = normalize(goodsE)
    A0 = find_A0(rE, mE)
    if A0 > maxA0[0]:
        maxA0 = (A0, E)
    box_ct = 0
    tight_here = 0
    for a in range(1, A0):
        if a in E:
            continue
        r_a, m_a, base = subtract_arcs(GE, a)
        assert m_a > 0, "base %s+{%d} has m=0" % (E, a)
        v0sq = F(r_a * r_a, 18) / (m_a * m_a)
        bmax = isqrt(v0sq.numerator // v0sq.denominator) + 2
        while F(bmax * bmax) > v0sq:
            bmax -= 1
        for b in range(a + 1, bmax + 1):
            if b in E:
                continue
            box_ct += 1
            _, L, _ = subtract_arcs(base, b)
            if L == 0:
                fam = frozenset(E) | {a, b}
                mq = missing_q(sorted(fam))
                tights[fam] = mq
                tight_here += 1
                if not mq:
                    covering_fails.append((E, a, b))
    tot_box += box_ct
    if box_ct > worst_box[0]:
        worst_box = (box_ct, E)
    if bi % 28 == 0 or tight_here:
        print("  body %3d/364 %-46s r=%2d m=%.5f A0=%4d box=%5d tights=%d  [%.0fs]"
              % (bi + 1, "{" + ",".join(map(str, E)) + "}", rE, float(mE), A0, box_ct,
                 tight_here, time.time() - t0))

print("\n" + "=" * 100)
print("GLOBAL SUMMARY  (%.1f s total)" % (time.time() - t0))
print("  bodies: %d ; total box pairs swept exactly: %d ; worst box %d at %s"
      % (len(bodies), tot_box, worst_box[0], worst_box[1]))
print("  max A0 = %d at %s ; min body measure m = %s ~ %.5f at %s"
      % (maxA0[0], maxA0[1], minm[0], float(minm[0]), minm[1]))
print("\nTIGHT-FAMILY CENSUS (all L=0 families found in any box), deduped: %d" % len(tights))
for fam, mq in sorted(tights.items(), key=lambda kv: sorted(kv[0])):
    print("   %-40s missing q=%s -> THM-366 t=1/%d"
          % ("{" + ",".join(map(str, sorted(fam))) + "}", mq, mq[0] if mq else 0))
print("\nCOVERING pairs with L=0 (would block THM-734): %s" % (covering_fails if covering_fails else "NONE"))
if not covering_fails:
    print("\nTHM-734 ESTABLISHED: every 13-speed family with >=11 speeds in {1..14} satisfies LRC(14).")
    print("  (L>0 except at the census tights above, each non-covering, dispatched by THM-366.)")
else:
    print("\nTHM-734 BLOCKED -- INVESTIGATE the covering L=0 pairs above.")
