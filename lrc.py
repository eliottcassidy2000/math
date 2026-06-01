"""
Lonely Runner Conjecture: moment method / anti-concentration investigation.
Exact rational arithmetic throughout (fractions.Fraction).

Setup:
  distinct positive integer speeds v_1..v_m, observer 0, n = m+1, threshold 1/n.
  B_i = { t in [0,1) : ||v_i t|| < 1/n }, ||x|| = dist to nearest integer.
  N(t) = #{ i : t in B_i }.  Observer lonely at t iff N(t)=0.
  LRC: exists t with N(t)=0  <=>  the B_i do NOT cover [0,1).

B_i is a union of v_i open intervals each of width 2/(n v_i), total measure 2/n.
  ||v_i t|| < 1/n  <=>  for some integer k, |v_i t - k| < 1/n
                  <=>  t in ( (k - 1/n)/v_i , (k + 1/n)/v_i ),  k=0..v_i-1 within [0,1)
  (k=0 gives an arc straddling 0, wrapping around 1).

All endpoints are rationals with denominator n*v_i. So we work on the circle R/Z
with exact Fractions. The "lonely set" is a finite union of open rational intervals.
"""

from fractions import Fraction as F
from math import gcd, exp
from itertools import combinations


def arc_endpoints(v, n):
    """Return list of (lo, hi) open intervals (mod 1) for B_v, as Fractions in [0,1).
    Each interval ((k - 1/n)/v, (k + 1/n)/v) reduced mod 1 and split at 0 if it wraps."""
    half = F(1, n)  # half-width in 'k - units': ||vt||<1/n
    raw = []
    for k in range(v):
        lo = F(k, v) - half / v
        hi = F(k, v) + half / v
        raw.append((lo, hi))
    # reduce mod 1, splitting wraps
    out = []
    for lo, hi in raw:
        lo_m = lo % 1
        hi_m = hi % 1
        if lo % 1 <= hi % 1 and hi - lo < 1:
            # check wrap: if lo<0 then lo_m large, hi_m small
            if lo_m <= hi_m:
                out.append((lo_m, hi_m))
            else:
                out.append((lo_m, F(1)))
                out.append((F(0), hi_m))
        else:
            out.append((lo_m, F(1)))
            out.append((F(0), hi_m))
    return out


def measure_union_intervals(intervals):
    """Exact measure of union of (lo,hi) open intervals on [0,1)."""
    iv = sorted(intervals)
    total = F(0)
    cur_lo = None
    cur_hi = None
    for lo, hi in iv:
        if cur_lo is None:
            cur_lo, cur_hi = lo, hi
        elif lo <= cur_hi:
            if hi > cur_hi:
                cur_hi = hi
        else:
            total += cur_hi - cur_lo
            cur_lo, cur_hi = lo, hi
    if cur_lo is not None:
        total += cur_hi - cur_lo
    return total


def measure_B(v, n):
    return measure_union_intervals(arc_endpoints(v, n))


def intersect_two(ints1, ints2):
    """Intersection of two unions of intervals -> list of intervals."""
    out = []
    for a, b in ints1:
        for c, d in ints2:
            lo = max(a, c)
            hi = min(b, d)
            if lo < hi:
                out.append((lo, hi))
    return out


def measure_intersection(speeds, n):
    """Exact measure of intersection of B_{v} over speeds list."""
    cur = arc_endpoints(speeds[0], n)
    for v in speeds[1:]:
        cur = intersect_two(cur, arc_endpoints(v, n))
        if not cur:
            return F(0)
    return measure_union_intervals(cur)


def lonely_measure(speeds, n):
    """Exact measure of { t : N(t)=0 } = measure of complement of union B_i.
    Use a sweep over all endpoints (cell decomposition)."""
    all_arcs = [arc_endpoints(v, n) for v in speeds]
    # collect breakpoints
    pts = set([F(0), F(1)])
    for arcs in all_arcs:
        for lo, hi in arcs:
            pts.add(lo)
            pts.add(hi)
    pts = sorted(pts)
    lonely = F(0)
    for a, b in zip(pts, pts[1:]):
        mid = (a + b) / 2
        covered = False
        for arcs in all_arcs:
            for lo, hi in arcs:
                if lo < mid < hi:
                    covered = True
                    break
            if covered:
                break
        if not covered:
            lonely += b - a
    return lonely


if __name__ == "__main__":
    # sanity: measure_B == 2/n
    for n, v in [(5, 3), (5, 7), (7, 11)]:
        m = n - 1
        assert measure_B(v, n) == F(2, n), (v, n, measure_B(v, n))
    print("measure(B_v) == 2/n verified.")
