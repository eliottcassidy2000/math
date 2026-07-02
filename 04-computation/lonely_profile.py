#!/usr/bin/env python3
"""
lonely_profile -- exact piecewise-linear survival profile of the LRC lonely measure.

ENGINEERING DELIVERABLE (mac-mini-2026-07-01-S93, THM-592/THM-593 companion library).

For a finite speed set S ⊂ Z_{>0}, the clearance process f_S(t) = min_{v∈S} ||v t||
on the circle has survival function m_S(r) = |{t : f_S(t) >= r}| (the "lonely measure").
THM-592: m_S is piecewise linear in r with breakpoints on {d/(v+w)} ∪ {d/(w-v)} and
slope -Σ_components(1/v_L + 1/v_R) (the co-area formula).  This module computes the
ENTIRE exact profile in one sweep, plus derived quantities.

API:
    profile(S, rmax=None) -> Profile
        .knots          list[(r, m)] exact Fractions, the full piecewise-linear graph
        .measure(r)     exact m_S(r)
        .slope(r)       exact left-slope at r (co-area density, negated)
        .M()            exact covering-min M(S) (root of the profile / max of f_S)
        .convex_on(a,b) True iff no concave kink in (a,b)
        .concave_kinks((a,b)) list of (r, defect) overtaking events with slope defects
        .tangent_reach(r0)  root of the tangent at r0 (ladder certificate radius)
        .critical_slope(q)  c with m(r) = c*(1-q*r) on the last cell below 1/q

Use cases: exact M(S) certificates; tangent-ladder covering-floor certificates
(THM-592(v)); tight-set slope floors (THM-593); frequency-planning "danger profiles"
(jitter margin of a frequency set = the profile of its beat spectrum).

Everything is exact rational arithmetic (fractions.Fraction); no floats internally.
"""
from fractions import Fraction as F
from math import gcd
from bisect import bisect_right

__all__ = ["profile", "Profile"]


def _danger_intervals(S, r):
    iv = []
    for v in S:
        half = F(r) / v
        for k in range(v):
            c = F(k, v)
            a, b = c - half, c + half
            if a < 0:
                iv.append((a + 1, F(1))); iv.append((F(0), b))
            elif b > 1:
                iv.append((a, F(1))); iv.append((F(0), b - 1))
            else:
                iv.append((a, b))
    return sorted(iv)


def _measure(S, r):
    tot = F(0); ca = cb = None
    for a, b in _danger_intervals(S, r):
        if cb is None:
            ca, cb = a, b
        elif a <= cb:
            cb = max(cb, b)
        else:
            tot += cb - ca; ca, cb = a, b
    if cb is not None:
        tot += cb - ca
    return 1 - tot


def _breakpoints(S, rmax):
    pts = set()
    Sl = sorted(set(S))
    for i, v in enumerate(Sl):
        for j in range(i, len(Sl)):
            w = Sl[j]
            for den in ({v + w} | ({w - v} if w > v else set())):
                if den <= 0:
                    continue
                d = 1
                while F(d, den) <= rmax:
                    pts.add(F(d, den)); d += 1
    return sorted(pts)


class Profile:
    def __init__(self, S, rmax):
        self.S = sorted(set(S))
        self.rmax = F(rmax)
        bps = [p for p in _breakpoints(self.S, self.rmax) if p < self.rmax]
        xs = [F(0)] + bps + [self.rmax]
        # exact measure at each breakpoint; slopes per cell from two interior samples
        knots = [(F(0), F(1))]
        self._cells = []  # (a, b, slope, m_at_a)
        for i in range(len(xs) - 1):
            a, b = xs[i], xs[i + 1]
            if b <= a:
                continue
            r1 = a + (b - a) / 4
            r2 = a + 3 * (b - a) / 4
            m1, m2 = _measure(self.S, r1), _measure(self.S, r2)
            sl = (m2 - m1) / (r2 - r1)
            m_a = m1 - sl * (r1 - a)
            self._cells.append((a, b, sl, m_a))
            knots.append((b, m_a + sl * (b - a)))
        self.knots = knots

    # -- queries ------------------------------------------------------------
    def measure(self, r):
        r = F(r)
        for a, b, sl, m_a in self._cells:
            if a <= r <= b:
                return m_a + sl * (r - a)
        raise ValueError("r outside profiled range")

    def slope(self, r):
        r = F(r)
        for a, b, sl, m_a in self._cells:
            if a < r <= b:
                return sl
        return self._cells[0][2] if self._cells else F(0)

    def M(self):
        """Exact covering-min: the radius where the profile hits 0 (if <= rmax)."""
        for a, b, sl, m_a in self._cells:
            m_b = m_a + sl * (b - a)
            if m_a > 0 and m_b <= 0:
                return a + m_a / (-sl)
            if m_a == 0:
                return a
        return None  # M(S) > rmax

    def concave_kinks(self, lo=None, hi=None):
        lo = F(lo) if lo is not None else F(0)
        hi = F(hi) if hi is not None else self.rmax
        out = []
        for i in range(1, len(self._cells)):
            r = self._cells[i][0]
            if lo < r < hi:
                d = self._cells[i][2] - self._cells[i - 1][2]
                if d < 0:
                    out.append((r, -d))
        return out

    def convex_on(self, lo=None, hi=None):
        return not self.concave_kinks(lo, hi)

    def tangent_reach(self, r0):
        """Root of the tangent at r0: certified positivity reach under convexity."""
        r0 = F(r0)
        m, sl = self.measure(r0), self.slope(r0 + F(1, 10**12))
        if sl >= 0:
            return None
        return r0 + m / (-sl)

    def defect(self, r0, r1):
        """K(r0,r1): total concave slope defect (THM-592(v) ladder correction)."""
        return sum(d for _, d in self.concave_kinks(r0, r1))

    def critical_slope(self, q):
        """c with m(r) = c*(1 - q r) on the last linear cell below 1/q (tight sets)."""
        rq = F(1, q)
        cells = [c for c in self._cells if c[1] <= rq and c[0] < rq]
        a, b, sl, m_a = cells[-1]
        r1 = a + (b - a) / 3
        return self.measure(r1) / (1 - q * r1)


def profile(S, rmax=None):
    if rmax is None:
        rmax = F(1, 2)
    return Profile(S, rmax)


def convexity_criterion(S, rho):
    """THM-592(iv): sufficient condition for convexity of m_S on (0, rho]."""
    Sl = sorted(set(S))
    bad = [(v, w, F(gcd(v, w), w - v)) for i, v in enumerate(Sl)
           for w in Sl[i + 1:] if F(gcd(v, w), w - v) < rho]
    return (not bad), bad


if __name__ == "__main__":
    # self-test against THM-592/593 verified values
    AP = list(range(1, 14))
    p = profile(AP, F(1, 14))
    assert p.critical_slope(14) == F(1666, 6435), p.critical_slope(14)
    assert p.convex_on(), "AP must be convex below 1/14"
    ok, bad = convexity_criterion(AP, F(1, 14))
    assert ok
    DW = list(range(1, 13)) + [182]
    pd = profile(DW, F(1, 12))  # window must contain M(DW)=14/183 > 1/14
    assert pd.M() == F(14, 183), pd.M()
    kk = pd.concave_kinks()
    assert kk and kk[0][0] == F(1, 181), kk[:1]
    assert pd.tangent_reach(F(1, 16)) >= F(1, 14)
    p8 = profile([1, 4, 5, 6, 7, 11, 13], F(1, 8))
    assert p8.critical_slope(8) == F(328, 1001)
    assert p8.M() == F(1, 8)
    print("lonely_profile self-test: ALL PASS")
    print(f"  AP{{1..13}} slope 1666/6435, convex; DW M=14/183, first kink 1/181, "
          f"tangent@1/16 reaches {float(pd.tangent_reach(F(1,16))):.5f} >= 1/14")
