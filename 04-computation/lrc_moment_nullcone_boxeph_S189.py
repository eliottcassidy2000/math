#!/usr/bin/env python3
"""
lrc_moment_nullcone_boxeph_S189.py  (HYP-8615, THM-1820)

LRC AS A MOMENT-NULLCONE PROBLEM.

(B1) THE BRIDGE IDENTITY (Weyl/Poisson): for integer speeds v and 1-periodic
     f_j: int_0^1 prod_j f_j(v_j t) dt = sum_{k: k.v=0} prod_j fhat_j(k_j).
     The RELATION LATTICE {k.v = 0} = the charge lattice of GMC; the
     product pairing = the moment functional. Verified numerically here
     (truncated Fourier vs exact interval arithmetic).
     One-sided match: Q-independent speeds -> only k = 0 -> |G_delta| =
     (1-2 delta)^n > 0: loneliness free = charge-imbalance triviality.

(B2) delta* vs relation-lattice richness: delta*(v) = max_t min_j ||v_j t||
     computed exactly-ish (fine sweep + refinement) for integer speed sets
     n = 3, 4; tabulated against N_R(v) = #{k != 0 : k.v = 0, |k|_inf <= 3}
     (small-relation count). Reframe prediction: tight families (delta*
     near 1/(n+1)) are relation-RICH; Kronecker-free direction: richer
     lattice = more conspiracy. Report the correlation honestly.

(B3) THE delta-LADDER (the S183 stacking analog in LRC): sweep delta and
     count connected components of the good set G_delta (exact interval
     union on the circle; integer speeds => rational endpoints). The
     component-death events are the FOLDS; the tight family v = (1..n)
     should show STACKED deaths (symmetry: many components die at the same
     delta) vs spread deaths for generic speeds. Measured.

boxeph-2026-07-20-S189. Pure python.
"""
from fractions import Fraction
import math
import itertools

# ---------- (B1) bridge identity ----------
print("=" * 78)
print("(B1) bridge identity: int prod f(v_j t) dt = sum_{k.v=0} prod fhat(k_j)")
v = (1, 2, 3)
delta = 0.11
# LHS: measure of {t : all ||v_j t|| <= delta} by exact interval intersection
def all_close_measure(v, delta):
    # intervals for runner j: ||v_j t|| <= delta <=> t in union_m [(m-delta)/v_j, (m+delta)/v_j]
    # intersect over j on [0,1): do it by fine exact sweep over breakpoints
    pts = set([Fraction(0), Fraction(1)])
    dl = Fraction(delta).limit_denominator(10 ** 6)
    for vj in v:
        for m in range(0, vj + 1):
            pts.add(Fraction(m - dl, vj) % 1)
            pts.add(Fraction(m + dl, vj) % 1)
    pts = sorted(pts)
    tot = Fraction(0)
    for i in range(len(pts) - 1):
        mid = (pts[i] + pts[i + 1]) / 2
        ok = True
        for vj in v:
            x = (vj * mid) % 1
            if min(x, 1 - x) > dl:
                ok = False
                break
        if ok:
            tot += pts[i + 1] - pts[i]
    return float(tot)

lhs = all_close_measure(v, delta)
# RHS: sum over relation lattice k.v = 0, |k| <= K, of prod sinc weights
def fhat(kj, delta):
    if kj == 0:
        return 2 * delta
    return math.sin(2 * math.pi * kj * delta) / (math.pi * kj)

K = 40
rhs = 0.0
for k1 in range(-K, K + 1):
    for k2 in range(-K, K + 1):
        # k3 determined: k1*1 + k2*2 + k3*3 = 0
        num = -(k1 * v[0] + k2 * v[1])
        if num % v[2]:
            continue
        k3 = num // v[2]
        if abs(k3) > K:
            continue
        rhs += fhat(k1, delta) * fhat(k2, delta) * fhat(k3, delta)
print("  v=%s delta=%.2f: LHS (exact intervals) = %.6f ; RHS (lattice sum, K=%d) = %.6f"
      % (v, delta, lhs, K, rhs))
print("  Q-independent analog: only k=0 term = (2 delta)^n = %.6f (the 'one-sided' free case)"
      % ((2 * delta) ** 3))

# ---------- shared: exact good set for integer speeds ----------
def good_set_components(v, dl):
    # G = {t in [0,1) : ||v_j t|| >= dl for all j}: complement of union of
    # open intervals; exact endpoints as Fractions; return component count
    # and total measure
    bad = []
    for vj in v:
        for m in range(0, vj + 1):
            lo = Fraction(m - dl, vj)
            hi = Fraction(m + dl, vj)
            bad.append((lo, hi))
    # clip to [0,1) with wrap: normalize pieces
    pieces = []
    for lo, hi in bad:
        lo_ = lo % 1
        length = hi - lo
        if length >= 1:
            return 0, Fraction(0)
        hi_ = lo_ + length
        if hi_ <= 1:
            pieces.append((lo_, hi_))
        else:
            pieces.append((lo_, Fraction(1)))
            pieces.append((Fraction(0), hi_ - 1))
    pieces.sort()
    merged = []
    for lo, hi in pieces:
        if merged and lo <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], hi))
        else:
            merged.append((lo, hi))
    # wrap-merge first/last
    if len(merged) >= 2 and merged[0][0] == 0 and merged[-1][1] == 1:
        pass  # components computed on the circle below
    # good components on the circle:
    total_bad = sum(hi - lo for lo, hi in merged)
    if total_bad >= 1:
        return 0, Fraction(0)
    ncomp = len(merged)
    # on the circle, gaps between consecutive merged bad intervals = good comps
    if ncomp == 0:
        return 1, Fraction(1)
    # check wrap adjacency
    wraps = (merged[0][0] == 0 and merged[-1][1] == 1)
    goodm = Fraction(1) - total_bad
    comps = ncomp - (1 if wraps else 0)
    if comps == 0:
        comps = 1 if goodm > 0 else 0
    return comps, goodm

# ---------- (B2) delta* vs relation richness ----------
print()
print("=" * 78)
print("(B2) delta* vs relation-lattice richness (n = 3, 4 batteries)")

def delta_star(v):
    lo, hi = Fraction(0), Fraction(1, 2)
    for _ in range(40):
        mid = (lo + hi) / 2
        comps, gm = good_set_components(v, mid)
        if gm > 0:
            lo = mid
        else:
            hi = mid
    return float(lo)

def relation_count(v, box=3):
    n = len(v)
    cnt = 0
    for k in itertools.product(range(-box, box + 1), repeat=n):
        if any(k) and sum(ki * vi for ki, vi in zip(k, v)) == 0:
            cnt += 1
    return cnt // 2  # up to sign

print("  n=3 (threshold 1/4 = 0.25):")
rows = []
for v3 in itertools.combinations(range(1, 10), 3):
    ds = delta_star(v3)
    nr = relation_count(v3)
    rows.append((ds, nr, v3))
rows.sort()
for ds, nr, v3 in rows[:6]:
    print("    v=%-12s delta* = %.4f  N_R = %d" % (str(v3), ds, nr))
print("    ... (loosest:)")
for ds, nr, v3 in rows[-3:]:
    print("    v=%-12s delta* = %.4f  N_R = %d" % (str(v3), ds, nr))
tight = [r for r in rows if r[0] < 0.26]
loose = [r for r in rows if r[0] > 0.30]
avg_t = sum(r[1] for r in tight) / max(1, len(tight))
avg_l = sum(r[1] for r in loose) / max(1, len(loose))
print("    mean N_R: tight(delta*<0.26) = %.2f (%d sets) vs loose(>0.30) = %.2f (%d sets)"
      % (avg_t, len(tight), avg_l, len(loose)))

# ---------- (B3) the delta-ladder: stacked vs spread component deaths ----------
print()
print("=" * 78)
print("(B3) delta-ladder: component-death events (folds), tight vs generic")
for v_ in ((1, 2, 3, 4), (1, 2, 5, 7), (2, 3, 7, 11)):
    # sweep delta upward, record component count changes
    events = []
    prev = None
    steps = 400
    for i in range(1, steps):
        dl = Fraction(i, 2 * steps)  # delta in (0, 1/2)
        comps, gm = good_set_components(v_, dl)
        if prev is not None and comps != prev:
            events.append((float(dl), prev, comps))
        prev = comps
        if comps == 0:
            break
    drops = [(d, a - b) for d, a, b in events if a > b]
    big = [e for e in drops if e[1] >= 2]
    print("  v=%s: delta* = %.4f ; death events: %d ; multi-deaths (STACKED, >=2 at once): %d %s"
          % (v_, delta_star(v_), len(drops), len(big),
             ["%d@%.3f" % (m, d) for d, m in big[:5]]))
print()
print("  (prediction: the tight consecutive family stacks its deaths — symmetry;")
print("   generic speeds die one at a time — the LRC face of S183's stacked germs)")
print("\nDONE.")
