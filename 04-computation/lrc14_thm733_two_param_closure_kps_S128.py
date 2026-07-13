#!/usr/bin/env python3
"""
lrc14_thm733_two_param_closure_kps_S128.py
==========================================
kind-pasteur-2026-07-13-S128.  THM-733: for EVERY 11 < a < b, the 13-speed family
{1,...,11, a, b} satisfies LRC(14).

Proof structure (everything exact rational arithmetic; THM-731 chain + THM-732 tail):

  UNIFORM LEG (a >= A0). Peel lemmas (elementary, proved in THM-733 file):
     P1: r_a := r({1..11,a}) <= a*m + (15/7)*r
     P2: m_a := m({1..11,a}) >= (6/7)*m - 8r/(49a)          (r, m = data of {1..11})
  THM-732(iii) tail: the certificate for {1..11,a,b} peeling b holds when b^2 > r_a^2/(18 m_a^2).
  Sufficient: R(a) <= K * M(a) * (a+1) with K = 2121/500 < sqrt(18), R(a)=a*m+(15/7)r,
  M(a)=(6/7)m-8r/(49a).  f(a) := K*M(a)*(a+1) - R(a) is EXACTLY increasing
  (f = [K*(6/7)m*(a+1) - a*m] + [-(8Kr/49)*(1+1/a)] - (15/7)r: first bracket has slope
  m*(6K/7 - 1) > 0, second is increasing since 1/a decreases). So ONE exact check f(A0) > 0
  certifies ALL a >= A0, giving b >= a+1 > v0(base_a) => L > 0.

  PER-a LEG (11 < a < A0). Exact (r_a, m_a); every b with b^2 > v0sq(a) := r_a^2/(18 m_a^2)
  closes by the tail. Remaining finite box {(a,b): a < b <= sqrt(v0sq(a))}: exact sweep of
  L({1..11,a,b}). L > 0 => proved. L = 0 => must be non-covering (THM-366 handles LRC), expected
  ONLY at the tight AP (a,b)=(12,13); any covering L=0 pair would be flagged LOUDLY.
"""
from fractions import Fraction as F
from math import isqrt
import sys, time

sys.setrecursionlimit(10000)
ONE = F(1)

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
    return out

def good_intervals(speeds):
    pieces = []
    for w in speeds:
        pieces += bad_pieces(w)
    pieces.sort()
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

def subtract_arcs(goods, w):
    """Good intervals minus D_w: returns (r_new, m_new). goods may contain a wrap interval (b>1)."""
    # normalize goods into [0,1] pieces
    G = []
    for a, b in goods:
        if b <= 1:
            G.append((a, b))
        else:
            G.append((a, ONE)); G.append((F(0), b - 1))
    G.sort()
    B = sorted(bad_pieces(w))
    out = []
    j = 0
    for a, b in G:
        cur = a
        # advance bad pointer
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
    # note: wrap-rejoin across 0/1 affects r by at most 1 (two pieces are one circular interval);
    # for the tail bound a SMALLER r is what we need, so we conservatively count pieces (an OVERcount
    # of r) — safe for v0; measure is exact.
    m_new = sum(b - a for a, b in out)
    return len(out), m_new, out

# ------------------------------------------------------------------
print("=" * 100)
print("THM-733: {1..11, a, b} satisfies LRC(14) for ALL 11 < a < b")
print("=" * 100)
t0 = time.time()

BODY = list(range(1, 12))
goods0, r0, m0 = good_intervals(BODY)
print("body {1..11}: r = %d intervals, m = %s ~ %.6f" % (r0, m0, float(m0)))

# ---- uniform leg: exact A0 ----
K = F(2121, 500)               # K < sqrt(18) since K^2 = 4498641/250000 < 18 = 4500000/250000
assert K * K < 18
def f(a):
    R = a * m0 + F(15, 7) * r0
    M = F(6, 7) * m0 - F(8 * r0, 49 * a)
    return K * M * (a + 1) - R

A0 = None
for a in range(12, 100000):
    if f(a) > 0:
        A0 = a
        break
assert A0 is not None
assert f(A0) > 0 and f(A0 + 1) > 0 and f(10 * A0) > 0   # spot confirms of proved monotonicity
print("UNIFORM LEG: K=2121/500 < sqrt18 exact; f increasing (proved); minimal A0 with f(A0)>0: A0 = %d" % A0)
print("  => for all a >= %d and all b > a: certificate holds (peel b), L > 0.   [f(A0)=%s>0]"
      % (A0, str(f(A0))[:40]))

# ---- per-a leg ----
print("\nPER-a LEG: 11 < a < %d.  v0sq(a) = r_a^2/(18 m_a^2) exact; box = {b: a < b <= sqrt(v0sq)}"
      % A0)
box = []
maxv0 = 0
for a in range(12, A0):
    r_a, m_a, _ = subtract_arcs(goods0, a)
    assert m_a > 0, "m_a=0 at a=%d" % a
    v0sq = F(r_a * r_a, 18) / (m_a * m_a)
    # integer bmax with b^2 <= v0sq
    bmax = isqrt(v0sq.numerator // v0sq.denominator)
    while F((bmax + 1) * (bmax + 1)) <= v0sq:
        bmax += 1
    while F(bmax * bmax) > v0sq:
        bmax -= 1
    maxv0 = max(maxv0, bmax)
    for b in range(a + 1, bmax + 1):
        box.append((a, b))
print("  box size = %d pairs ; max needed b = %d ; (per-a exact thresholds computed for %d values of a)"
      % (len(box), maxv0, A0 - 12))

# ---- exact sweep of the box ----
fails = []
zeroL = []
for (a, b) in box:
    _, _, L = good_intervals(BODY + [a, b])
    if L == 0:
        # covering check: does every q in 2..14 divide some element?
        S = BODY + [a, b]
        missing = [q for q in range(2, 15) if not any(w % q == 0 for w in S)]
        zeroL.append((a, b, missing))
        if not missing:
            fails.append((a, b))
print("\nBOX SWEEP: %d exact rational sweeps done." % len(box))
print("  pairs with L = 0: %s" % (", ".join("(%d,%d) missing q=%s -> THM-366" % (a, b, mq)
                                             for a, b, mq in zeroL) if zeroL else "NONE"))
print("  covering pairs with L = 0 (would block the theorem): %s" % (fails if fails else "NONE"))

ok = not fails
print("\n" + "=" * 100)
if ok:
    print("THM-733 ESTABLISHED: for every 11 < a < b, {1,...,11,a,b} satisfies LRC(14).")
    print("  - a >= %d: uniform peel (P1/P2 + THM-732 tail), L > 0" % A0)
    print("  - a < %d, b > sqrt(v0sq(a)): exact per-a tail, L > 0" % A0)
    print("  - remaining %d pairs: exact rational sweep, L > 0 except %s (non-covering, THM-366)"
          % (len(box), [(a, b) for a, b, _ in zeroL] if zeroL else "none"))
else:
    print("THM-733 BLOCKED at covering L=0 pairs: %s -- INVESTIGATE" % fails)
print("elapsed %.1f s" % (time.time() - t0))
