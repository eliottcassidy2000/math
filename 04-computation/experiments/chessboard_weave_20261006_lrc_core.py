#!/usr/bin/env python3
"""Exact LRC primitives for the chessboard-weave session (2026-10-06).

All arithmetic is exact (Python ints / fractions.Fraction).  Conventions:
  v = tuple of n distinct positive integer speeds, N = n + 1, threshold c = 1/N,
  ||x|| = distance from x to the nearest integer,
  f_v(t) = min_i ||t v_i||   (1-periodic, f_v(1-t) = f_v(t)).

  delta(v)            loneliness  max_t f_v(t)
  safe_components(v)  the closed set Safe(v) = {t in [0,1] : f_v(t) >= 1/N}
                      as a sorted list of closed intervals (lo, hi) of integer
                      numerators over the common denominator D = N*lcm(v);
                      degenerate intervals lo == hi are isolated safe POINTS
  tau(v)              first lonely time  min Safe(v) cap (0,1)
  good_period(v)      least q >= 1 such that some k/q (0<k<q) lies in Safe(v)

delta algorithm (proved in the A-script output): a local maximum of the
min of tent functions is either the peak of one tent or a crossing of an
increasing piece of one tent with a decreasing piece of another; both are of
the form t = m/(v_i + v_j) with i <= j.  Crossings of two increasing (or two
decreasing) pieces, t = m/|v_i - v_j|, therefore need not be added; the
brute-force checker below uses the full breakpoint lattice (including the
difference crossings) as an independent test.

Run as a script: self-tests (brute force cross-validation) -> .out file.
"""
from fractions import Fraction
from math import gcd
from functools import reduce
from itertools import combinations


def lcm(a, b):
    return a // gcd(a, b) * b


def lcm_list(v):
    return reduce(lcm, v, 1)


def gcd_list(v):
    return reduce(gcd, v, 0)


def primitive_sets(n, B):
    """All n-subsets of {1..B} (sorted tuples) with gcd 1."""
    for c in combinations(range(1, B + 1), n):
        if gcd_list(c) == 1:
            yield c


def dist_num(num, den):
    """den * ||num/den|| as an integer."""
    r = num % den
    return den - r if den - r < r else r


# ---------------------------------------------------------------- delta
def delta(v):
    """Exact loneliness of the speed set v (Fraction).  Candidates m/(v_i+v_j), i<=j."""
    vs = sorted(set(v))
    bestR, bestS = 0, 1
    k = len(vs)
    for a in range(k):
        for b in range(a, k):
            s = vs[a] + vs[b]
            for m in range(1, s // 2 + 1):
                R = s
                ok = True
                for w in vs:
                    r = (m * w) % s
                    if s - r < r:
                        r = s - r
                    if r < R:
                        R = r
                        if R * bestS <= bestR * s:
                            ok = False
                            break
                if ok and R * bestS > bestR * s:
                    bestR, bestS = R, s
    return Fraction(bestR, bestS)


def delta_argmax(v):
    """All t in (0,1/2] attaining delta (searched on the candidate set; with the
    proof that every local max is a candidate, this is the full argmax in (0,1/2])."""
    d = delta(v)
    vs = sorted(set(v))
    pts = set()
    for a in range(len(vs)):
        for b in range(a, len(vs)):
            s = vs[a] + vs[b]
            for m in range(1, s // 2 + 1):
                R = min(dist_num(m * w, s) for w in vs)
                if Fraction(R, s) == d:
                    pts.add(Fraction(m, s))
    return d, sorted(pts)


def delta_brute(v, Lmax=200000):
    """Independent check.  Every breakpoint of f lies in the union U of
    {k/(2v_i)} (tent peaks and zeros), {m/(v_i+v_j)} and {m/|v_i-v_j|} (crossings),
    and f is linear between consecutive breakpoints, so max f = max over U.
    If L = lcm(all 2v_i, v_i+v_j, |v_i-v_j|) <= Lmax we instead scan the whole
    lattice (1/L)Z (a superset of U); otherwise we scan U itself."""
    vs = sorted(set(v))
    dens = set()
    for i, a in enumerate(vs):
        dens.add(2 * a)
        for b in vs[i + 1:]:
            dens.add(a + b)
            dens.add(b - a)
    L = lcm_list(dens)
    if L <= Lmax:
        best = 0
        for k in range(0, L // 2 + 1):
            R = min(dist_num(k * w, L) for w in vs)
            if R > best:
                best = R
        return Fraction(best, L), "lattice"
    best = Fraction(0)
    for s in dens:
        for m in range(0, s // 2 + 1):
            R = min(dist_num(m * w, s) for w in vs)
            if Fraction(R, s) > best:
                best = Fraction(R, s)
    return best, "union"


# ---------------------------------------------------------------- safe set
def _intersect(A, B):
    out = []
    i = j = 0
    la, lb = len(A), len(B)
    while i < la and j < lb:
        a0, a1 = A[i]
        b0, b1 = B[j]
        lo = a0 if a0 > b0 else b0
        hi = a1 if a1 < b1 else b1
        if lo <= hi:
            out.append((lo, hi))
        if a1 < b1:
            i += 1
        else:
            j += 1
    return out


def safe_components(v, N=None):
    """Components of Safe(v) = {t in [0,1]: ||t w|| >= 1/N for all w in v}.
    Returns (comps, D): comps = sorted list of (lo, hi) integer numerators over D."""
    vs = sorted(set(v))
    if N is None:
        N = len(vs) + 1
    L = lcm_list(vs)
    D = N * L
    # process large speeds first (more, shorter intervals prune fastest)?  Order is
    # irrelevant for correctness; start from the smallest speed (fewest intervals).
    comps = None
    for w in vs:
        u = L // w
        ivs = [((m * N + 1) * u, (m * N + N - 1) * u) for m in range(w)]
        comps = ivs if comps is None else _intersect(comps, ivs)
        if not comps:
            break
    return comps, D


def tau(v, N=None):
    comps, D = safe_components(v, N)
    if not comps:
        return None
    return Fraction(comps[0][0], D)


def simplest_in(xn, xd, yn, yd):
    """Stern-Brocot simplest rational p/q in the CLOSED interval [xn/xd, yn/yd],
    0 <= x <= y.  It has the least denominator (and the least numerator) of all
    rationals in the interval."""
    if xn == 0:
        return (0, 1)
    fl = xn // xd
    if fl * xd == xn:
        return (fl, 1)
    if (fl + 1) * yd <= yn:
        return (fl + 1, 1)
    # fl < x <= y < fl + 1 ; x = fl + 1/z with z in [1/(y-fl), 1/(x-fl)]
    zp, zq = simplest_in(yd, yn - fl * yd, xd, xn - fl * xd)
    return (fl * zp + zq, zp)


def good_period_from(comps, D):
    best = None
    bestt = None
    for lo, hi in comps:
        p, q = simplest_in(lo, D, hi, D)
        if best is None or q < best:
            best, bestt = q, Fraction(p, q)
    return best, bestt


def good_period(v, N=None):
    comps, D = safe_components(v, N)
    if not comps:
        return None, None
    return good_period_from(comps, D)


# ---------------------------------------------------------------- brute forces
def tau_brute(v):
    vs = sorted(set(v))
    N = len(vs) + 1
    D = N * lcm_list(vs)
    for k in range(1, D):
        if all(N * dist_num(k * w, D) >= D for w in vs):
            return Fraction(k, D)
    return None


def tau_chase(v):
    """Third method: the fixed-point chase t <- max_i next_i(t)."""
    vs = sorted(set(v))
    N = len(vs) + 1
    t = Fraction(0)
    c = Fraction(1, N)
    while t < 1:
        nt = t
        for w in vs:
            x = t * w
            fl = x.numerator // x.denominator
            fr = x - fl
            if fr < c:
                cand = (fl + c) / w
            elif fr > 1 - c:
                cand = (fl + 1 + c) / w
            else:
                cand = t
            if cand > nt:
                nt = cand
        if nt == t:
            return t
        t = nt
    return None


def good_period_brute(v, qmax=10000):
    vs = sorted(set(v))
    N = len(vs) + 1
    for q in range(2, qmax):
        for k in range(1, q):
            if all(N * dist_num(k * w, q) >= q for w in vs):
                return q, Fraction(k, q)
    return None, None


def check(cond, msg=""):
    """Explicit runtime check (not an assert, so it also runs under python -O)."""
    if not cond:
        raise RuntimeError("CHECK FAILED: " + str(msg))


# ---------------------------------------------------------------- self-test
def _selftest():
    import random
    print("chessboard_weave_20261006_lrc_core self-test (exact arithmetic)")
    # 1. simplest_in vs brute force on all intervals with endpoints a/b, b <= 30
    cnt = 0
    fr = sorted({Fraction(a, b) for b in range(1, 31) for a in range(1, b)})
    rnd = random.Random(20261006)
    for _ in range(20000):
        x, y = sorted(rnd.sample(fr, 2))
        if rnd.random() < 0.1:
            y = x
        p, q = simplest_in(x.numerator, x.denominator, y.numerator, y.denominator)
        # brute force least denominator
        qq = 1
        while True:
            # least k with k/qq >= x
            k = -((-x.numerator * qq) // x.denominator)
            if Fraction(k, qq) <= y:
                break
            qq += 1
        check(q == qq and x <= Fraction(p, q) <= y and Fraction(p, q) == Fraction(k, qq), (x, y, p, q, qq))
        cnt += 1
    print(f"[1] simplest_in == brute least-denominator on {cnt} random closed intervals: OK")

    # 2. delta vs brute (full breakpoint lattice incl. difference crossings)
    cnt = 0
    modes = {"lattice": 0, "union": 0}
    for n in range(1, 6):
        Bn = {1: 25, 2: 25, 3: 16, 4: 12, 5: 11}[n]
        for v in combinations(range(1, Bn + 1), n):
            d1 = delta(v)
            d2, mode = delta_brute(v)
            check(d1 == d2, (v, d1, d2))
            modes[mode] += 1
            cnt += 1
    print(f"[2] delta(candidates m/(v_i+v_j), i<=j) == delta_brute on {cnt} sets "
          f"(all n-subsets of [1,B], n=1..5, B=25,25,16,12,11); brute mode counts {modes}: OK")

    # 3. tau: interval sweep == brute grid == chase
    cnt = 0
    for n in range(1, 6):
        Bn = {1: 30, 2: 30, 3: 18, 4: 11, 5: 10}[n]
        for v in combinations(range(1, Bn + 1), n):
            t1 = tau(v)
            t2 = tau_brute(v)
            t3 = tau_chase(v)
            check(t1 == t2 == t3, (v, t1, t2, t3))
            cnt += 1
    print(f"[3] tau(sweep) == tau(brute grid 1/((n+1)lcm)) == tau(chase) on {cnt} sets "
          f"(all n-subsets of [1,B], n=1..5, B=30,30,18,11,10): OK")

    # 4. good period: interval/Stern-Brocot == brute force over k/q
    cnt = 0
    for n in range(1, 6):
        Bn = {1: 30, 2: 40, 3: 24, 4: 16, 5: 13}[n]
        for v in combinations(range(1, Bn + 1), n):
            q1, t1 = good_period(v)
            q2, t2 = good_period_brute(v)
            check(q1 == q2, (v, q1, q2))
            cnt += 1
    print(f"[4] good_period(Stern-Brocot on components) == brute k/q search on {cnt} sets "
          f"(all n-subsets of [1,B], n=1..5, B=30,40,24,16,13): OK")
    print("self-test PASSED")


if __name__ == "__main__":
    _selftest()
