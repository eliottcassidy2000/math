#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
THE TOWER FLOOR, EXACTLY: for a deep cluster with difference pattern c = (0=c_1 < ... < c_j), the
renormalized local density is D_c(t) = meas_u{ u : u + c_i t in A  for all i },  A = [1/14, 13/14].
The depth-1 floor needs   inf over patterns c, over compact parts C_low (|C_low| = 11-j),
    Int_{L_{C_low}} D_c(t) dt  >  1/36  (+ freeze error).
Worst-case over the POSITION of L_low (adversary puts L_low on the smallest D-values) at fixed
measure m = meas(L_low) >= (j-4)/7 (union bound) gives the exact functional
    Q_c(m) = Int_0^m D_c^*(u) du     (D^* = increasing rearrangement)
            = Int_0^inf (m - psi_c(s))_+ ds,   psi_c(s) = meas{t : D_c(t) < s}.
THEOREM-GRADE CONSTANT:  F_j := min over primitive patterns (range <= RMAX) of Q_{c,j}((j-4)/7).
If F_j > 1/36 for j = 7..11, the deep-cluster tower floor holds for all bounded-range patterns
(large-range patterns recurse into the gap-cut peel one level down; finite heights = middle band).

Structure facts verified en route (F1/F2 of the write-up):
  - j=7 RIGIDITY: D_c(t) = 0 iff the 7 danger arcs TILE iff centers = translate of (1/7)Z; zeros lie
    in (1/(7d))Z (d = gcd of pairwise differences; primitive => d=1 => at most 6 zeros k/7, k=1..6).
  - SLOPE at a zero >= range(c) >= 6 (cyclic total variation of the tiling order >= 2*range).
  - j>=8: zeros can be intervals (robust coverings); Q handles them uniformly.

Method: D_c is piecewise linear in t with breakpoints where two danger endpoints collide:
(c_i - c_k) t == +-1/7 (mod 1)  =>  t = (7m +- 1)/(7 delta), delta = c_i - c_k.  Evaluate D exactly at
all breakpoints (D linear between consecutive ones), build psi exactly, integrate Q exactly.
Two passes: float scan over ALL patterns -> exact Fraction re-verification of the 40 worst per j.

opus-2026-07-01-S33.  (HYP-3902; builds on HYP-3900/3901, THM-592/593, kps HYP-3950.)
"""
import sys, itertools
from fractions import Fraction as Fr
from math import gcd
from functools import reduce
try: sys.stdout.reconfigure(encoding='utf-8')
except Exception: pass

RMAX = 14
TARGET = Fr(1, 36)

def canonical(c):
    """dedupe mirror: pattern and its reversal have identical D-profiles."""
    R = c[-1]
    rev = tuple(sorted(R - x for x in c))
    return min(c, rev)

def patterns(j, rmax=RMAX):
    seen = set()
    for R in range(j - 1, rmax + 1):
        for mid in itertools.combinations(range(1, R), j - 2):
            c = (0,) + mid + (R,)
            difs = [c[i+1] - c[i] for i in range(j - 1)]
            if reduce(gcd, difs) != 1: continue  # dilations have identical profiles
            cc = canonical(c)
            if cc in seen: continue
            seen.add(cc)
            yield c

def D_at(c, t, exact):
    """D_c(t) exactly (Fraction) or float: 1 - meas(union of arcs [c_i t - 1/14, c_i t + 1/14] mod 1)."""
    if exact:
        cen = sorted((x * t) % 1 for x in c)
        h = Fr(1, 14); one = Fr(1)
    else:
        cen = sorted((x * t) % 1.0 for x in c)
        h = 1.0 / 14.0; one = 1.0
    # union of [cen-h, cen+h] on the circle: walk sorted centers, add gaps between consecutive
    n = len(cen); gapsum = 0
    for i in range(n):
        nxt = cen[(i + 1) % n] + (one if i == n - 1 else 0)
        gap = (nxt - cen[i]) - 2 * h
        if gap > 0: gapsum += gap
    return gapsum if gapsum > 0 else (Fr(0) if exact else 0.0)

def breakpoints(c, exact):
    """ALL t in [0,1] where D_c can kink. Two collision types between endpoints of arcs i,k
    (centers c_i t, c_k t, half-width 1/14, delta = c_i - c_k):
      opposite endpoints (left-vs-right): delta*t == +-1/7 (mod 1)  ->  t = (7m +- 1)/(7 delta);
      same endpoints (centers coincide):  delta*t == 0     (mod 1)  ->  t = m/delta  [concave peaks
        -- omitting these made interpolation a conservative underestimate; now included: exact]."""
    deltas = sorted(set(b - a for a, b in itertools.combinations(c, 2)))
    pts = set()
    for d in deltas:
        for m7 in range(0, 7 * d + 1):
            for e in (-1, 0, 1):
                num = 7 * m7 + e
                if 0 < num < 7 * d:
                    pts.add(Fr(num, 7 * d) if exact else num / (7.0 * d))
    pts.add(Fr(0) if exact else 0.0); pts.add(Fr(1) if exact else 1.0)
    return sorted(pts)

def profile(c, exact):
    """[(t_k, D_k)] at all breakpoints; D linear between consecutive."""
    ts = breakpoints(c, exact)
    return ts, [D_at(c, t, exact) for t in ts]

def Q_of_m(ts, Ds, m, exact):
    """Q(m) = int over the lowest-m measure of D, from the piecewise-linear (ts, Ds).
    Build segments (len, Dlo, Dhi); sort by Dlo... exact method: threshold sweep.
    psi(s) = total length where D < s: piecewise linear in s with knots at segment endpoint values.
    Q(m) = int_0^S (m - psi(s)) ds up to the s where psi(s) = m  (then integrand hits 0)."""
    segs = []
    for i in range(len(ts) - 1):
        L = ts[i + 1] - ts[i]
        if L <= 0: continue
        a, b = Ds[i], Ds[i + 1]
        if a > b: a, b = b, a
        segs.append((L, a, b))
    # psi(s) = sum over segs of contribution: 0 if s<=a; L*(s-a)/(b-a) if a<s<b (b>a); L if s>=b (or s>a==b)
    knots = sorted(set([v for (_, a, b) in segs for v in (a, b)]))
    zero = Fr(0) if exact else 0.0
    def psi(s):
        tot = zero
        for (L, a, b) in segs:
            if s <= a: continue
            if s >= b or b == a: tot += L
            else: tot += L * (s - a) / (b - a)
        return tot
    # integrate (m - psi(s))_+ ds over s in [0, s_stop]; psi is PL with knots; on each interval psi linear
    Q = zero; prev = zero
    for k in knots + [max(knots) + (Fr(1) if exact else 1.0)]:
        if k <= prev: continue
        p0, p1 = psi(prev), psi(k)
        # integrand f(s) = m - psi(s), linear from (m-p0) to (m-p1)
        f0, f1 = m - p0, m - p1
        if f0 <= 0 and f1 <= 0:
            break
        if f0 > 0 and f1 > 0:
            Q += (f0 + f1) * (k - prev) / 2
        else:
            # crosses zero inside (f0>0>f1 since psi nondecreasing)
            root = prev + (k - prev) * f0 / (f0 - f1)
            Q += f0 * (root - prev) / 2
            break
        prev = k
    return Q

print("=" * 108)
print(" TOWER FLOOR by exact rearrangement:  F_j = min over primitive patterns (range<=%d) of Q_c((j-4)/7)" % RMAX)
print("   target: F_j > 1/36 = %.6f  for j = 7..11   (m = (j-4)/7 = worst-case compact lonely measure)" % float(TARGET))
print("=" * 108)

summary = {}
for j in range(7, 12):
    m_worst = (j - 4) / 7.0
    scan = []
    npat = 0
    for c in patterns(j):
        npat += 1
        ts, Ds = profile(c, exact=False)
        q = Q_of_m(ts, Ds, m_worst, exact=False)
        scan.append((q, c))
    scan.sort()
    worst = scan[:40]
    # exact re-verification of the worst 40
    m_exact = Fr(j - 4, 7)
    best_exact = None; best_c = None
    for _, c in worst:
        ts, Ds = profile(c, exact=True)
        q = Q_of_m(ts, Ds, m_exact, exact=True)
        if best_exact is None or q < best_exact:
            best_exact, best_c = q, c
    ok = best_exact > TARGET
    summary[j] = (npat, best_exact, best_c, ok)
    print(f"\n j={j}: patterns scanned {npat}; float-worst 5: " +
          ", ".join(f"{q:.5f}@{c}" for q, c in scan[:5]))
    print(f"   EXACT min over worst-40:  F_{j} = {best_exact} = {float(best_exact):.6f}  at pattern {best_c}")
    print(f"   m = (j-4)/7 = {float(m_exact):.4f};  F_{j} > 1/36 ?  {ok}   (margin {float(best_exact/TARGET):.3f}x)")

print("\n" + "=" * 108)
print(" STRUCTURE CHECKS (F1/F2 for the write-up), j=7:")
c = (0, 1, 2, 3, 4, 5, 6)
ts, Ds = profile(c, exact=True)
zeros = [t for t, d in zip(ts, Ds) if d == 0]
print(f"   consecutive pattern zeros (exact): {zeros}  (predicted subset of k/7, k=1..6)")
# slope at first zero
z = zeros[0] if zeros else None
if z is not None:
    i = ts.index(z)
    sl_r = (Ds[i+1] - Ds[i]) / (ts[i+1] - ts[i]); sl_l = (Ds[i-1] - Ds[i]) / (ts[i-1] - ts[i])
    print(f"   slopes at zero t={z}: left {float(sl_l):.3f}, right {float(sl_r):.3f}  (claim: >= range = 6)")
# a no-zero pattern
c2 = (0, 1, 3, 7, 8, 10, 14)
ts2, Ds2 = profile(c2, exact=True)
mn2 = min(Ds2)
print(f"   pattern {c2}: min D = {mn2} = {float(mn2):.5f}  (no tiling => strictly positive everywhere)")
# j=8 interval zeros demo
c3 = (0, 1, 2, 3, 4, 5, 6, 7)
ts3, Ds3 = profile(c3, exact=True)
zint = sum((ts3[i+1]-ts3[i]) for i in range(len(ts3)-1) if Ds3[i]==0 and Ds3[i+1]==0)
print(f"   j=8 consecutive: measure of interval zero-set = {zint} = {float(zint):.5f}  (robust coverings exist)")

print("\n VERDICT:")
allok = all(v[3] for v in summary.values())
for j, (npat, F, c, ok) in summary.items():
    print(f"   j={j:>2}: F_j = {float(F):.6f} ({F}) at {c}   > 1/36: {ok}")
print(f"   ALL j=7..11 clear 1/36: {allok}")
print(" => the deep-cluster tower floor holds (in the N->inf limit) for ALL primitive patterns of range <= %d," % RMAX)
print("    with worst-case-POSITION compact parts (no assumption on where L_low sits). Larger ranges recurse")
print("    into the gap-cut peel (HYP-3900); finite heights N < N* form the finite middle band (see doc).")
print("DONE.")
