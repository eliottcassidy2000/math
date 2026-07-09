#!/usr/bin/env python3
"""
HYP-5722: the mu-LEVEL robust threshold transfer (the piece THM-670 disclaims).
boxeph-2026-07-09-S2.  Pure Python, exact rationals.

THM-670 (monad-explorer-S5) transfers FIRST MOMENTS: E[W_{th+d}] >= E[W_th] - 6d,
and explicitly disclaims any transfer for the LEVEL-SET measure mu_th.
LEM-014's (H1) needs mu_{1/7+delta}(E) >= bar_k.  This script verifies:

  LEMMA (mu-level transfer).  With W = W_{1/7}, Wd = W_{1/7+delta}:
    (i)   W - 6*delta <= Wd <= W pointwise            [THM-670's core]
    (ii)  m1' >= m1 - 6d;  m2' <= m2;  m3' >= m3 - 18*d*m2
          (a^3 - (a-b)^3 <= 3*a^2*b for 0<=b<=a, + Wd >= (W-6d)_+)
    (iii) D3(m1,m2,m3) = m1/M + (m1 - m2/M)^2/(m2 - m3/M), M = 6/7, is
          INCREASING in m1, DECREASING in m2, INCREASING in m3 on the cone
          {m1 >= m2/M, m2 > m3/M} (which holds since 0 <= W <= M):
            dD3/dm1 = 1/M + 2(m1-m2/M)/(m2-m3/M)          >= 0
            dD3/dm2 = -2(m1-m2/M)/(M(m2-m3/M)) - (m1-m2/M)^2/(m2-m3/M)^2 <= 0
            dD3/dm3 = (m1-m2/M)^2/(M(m2-m3/M)^2)          >= 0
    =>  mu_{1/7+delta}(E) = P(Wd > 0) >= D3(m1-6d, m2, m3-18d*m2) =: D3_d
        and (PZ)  mu_{1/7+delta}(E) >= (m1-6d)^2/m2 =: PZ_d.

So the CLOSED theta=1/7 D3/PZ legs (THM-660/661) yield (H1)'s robust floor
with an explicit O(delta) perturbation -- no re-proof of the legs.

Table per extremizer shape: exact m1,m2,m3; D3_0 (cross-check vs published);
true mu_{1/7} and mu_{1/7+d}; D3_d, PZ_d at the LEM-014 ladder d = 3s/V for
V/s in {400,200,100,50}; and d_max (bisection) where D3_d hits bar_k.
"""
from fractions import Fraction as F
import itertools

ONE7 = F(1, 7)
M = F(6, 7)

# published union-bound bars (THM-663 table), k=8..13
BARS = {8: 0.6750, 9: 0.5622, 10: 0.4521, 11: 0.3312, 12: 0.1993,
        13: float(F(14249, 252252))}


def breakpoints(E):
    pts = {F(0), F(1)}
    for e in E:
        if e > 0:
            for m_ in range(e + 1):
                pts.add(F(m_, e))
    for a, b in itertools.combinations(E, 2):
        d = abs(a - b)
        if d > 0:
            for m_ in range(d + 1):
                pts.add(F(m_, d))
    return sorted(pts)


def cell_gaps(E, lo, hi):
    mid = (lo + hi) / 2
    ph = []
    for e in E:
        c = (e * mid).numerator // (e * mid).denominator
        ph.append((e * mid - c, F(-c), F(e)))
    ph.sort()
    k = len(ph)
    gaps = []
    for i in range(k - 1):
        _, c1, s1 = ph[i]
        _, c2, s2 = ph[i + 1]
        gaps.append((c2 - c1, s2 - s1))
    _, cb, sb = ph[0]
    _, ct, st = ph[-1]
    gaps.append((1 + cb - ct, sb - st))
    return gaps


class Engine:
    def __init__(self, E):
        self.E = sorted(E)
        bps = breakpoints(self.E)
        self.cells = [(bps[i], bps[i + 1]) for i in range(len(bps) - 1)
                      if bps[i] < bps[i + 1]]
        self._gaps = [cell_gaps(self.E, lo, hi) for lo, hi in self.cells]

    def moments(self, theta):
        """Exact (m1, m2, m3) of W_theta = Sum (g_i - theta)_+."""
        m1 = m2 = m3 = F(0)
        for (lo, hi), gaps in zip(self.cells, self._gaps):
            # subdivide the cell at each gap's crossing of theta
            cuts = {lo, hi}
            for c, s in gaps:
                if s != 0:
                    r = (theta - c) / s
                    if lo < r < hi:
                        cuts.add(r)
            pts = sorted(cuts)
            for u, v in zip(pts, pts[1:]):
                mid = (u + v) / 2
                a = F(0); b = F(0)
                for c, s in gaps:
                    if c + s * mid > theta:
                        a += c - theta
                        b += s
                # integrate (a + b x)^p over (u, v)
                if b == 0:
                    if a > 0:
                        m1 += a * (v - u)
                        m2 += a * a * (v - u)
                        m3 += a ** 3 * (v - u)
                else:
                    for p, acc in ((1, 'm1'), (2, 'm2'), (3, 'm3')):
                        val = ((a + b * v) ** (p + 1) - (a + b * u) ** (p + 1)) \
                              / (b * (p + 1))
                        if p == 1: m1 += val
                        elif p == 2: m2 += val
                        else: m3 += val
        return m1, m2, m3

    def mu_gt(self, theta):
        tot = F(0)
        for (lo, hi), gaps in zip(self.cells, self._gaps):
            ivs = []
            for c, s in gaps:
                if s == 0:
                    if c > theta:
                        ivs.append((lo, hi))
                else:
                    r = (theta - c) / s
                    if s > 0:
                        aa, bb = max(lo, r), hi
                    else:
                        aa, bb = lo, min(hi, r)
                    if aa < bb:
                        ivs.append((aa, bb))
            ivs.sort()
            cur = None
            for aa, bb in ivs:
                if cur is None:
                    cur = [aa, bb]
                elif aa <= cur[1]:
                    cur[1] = max(cur[1], bb)
                else:
                    tot += cur[1] - cur[0]
                    cur = [aa, bb]
            if cur:
                tot += cur[1] - cur[0]
        return tot


def D3(m1, m2, m3):
    if m2 - m3 / M <= 0:
        return None
    return m1 / M + (m1 - m2 / M) ** 2 / (m2 - m3 / M)


def D3_robust(m1, m2, m3, d):
    # VALID DOMAIN: the monotone-transfer argument needs the perturbed point
    # on the cone with nonnegative bracket: m1 - 6d >= m2/M. Outside, the
    # closed form is NOT a bound (the square term grows spuriously).
    if m1 - 6 * d < m2 / M:
        return None
    return D3(m1 - 6 * d, m2, m3 - 18 * d * m2)


SHAPES = [
    ("block8  {0..7}",  8,  list(range(8))),
    ("block9  {0..8}",  9,  list(range(9))),
    ("block10 {0..9}",  10, list(range(10))),
    ("block11 {0..10}", 11, list(range(11))),
    ("A11 tail-min",    11, [0, 3, 6, 8, 9, 12, 15, 18, 21, 24, 27]),
    ("block12 {0..11}", 12, list(range(12))),
    ("block13 AP{0..12}", 13, list(range(13))),
    ("7-struct M128 k13", 13, [0, 7, 14, 21, 26, 29, 37, 44, 51, 58, 67, 75, 82]),
]

LADDER = [F(3, 400), F(3, 200), F(3, 100), F(3, 50)]  # d = 3s/V at V/s=400..50 (s-normalized: d per unit slope ratio)


def main():
    print("HYP-5722: mu-level robust transfer table (exact rationals)")
    print("THM-670 gives E[W] transfer; this supplies the disclaimed mu-level "
          "transfer via D3/PZ monotonicity.\n")
    for name, k, E in SHAPES:
        s = max(E)
        eng = Engine(E)
        m1, m2, m3 = eng.moments(ONE7)
        mu0 = eng.mu_gt(ONE7)
        d30 = D3(m1, m2, m3)
        bar = BARS[k]
        print(f"=== {name}  (k={k}, s={s}, bar_k={bar:.4f}) ===")
        print(f"  m1={float(m1):.6f} m2={float(m2):.6f} m3={float(m3):.6f}")
        print(f"  D3_0={float(d30):.6f}  true mu_1/7={float(mu0):.6f}  "
              f"(D3_0 {'>=':s} bar: {float(d30) >= bar})")
        # sanity: pointwise lemma checks via true mu at lifted threshold
        for d in LADDER:
            Vs = float(3 / d)          # V/s corresponding to delta = 3s/V, s-normalized delta = 3/(V/s)
            d3d = D3_robust(m1, m2, m3, d)
            pzd = (m1 - 6 * d) ** 2 / m2 if m1 > 6 * d else F(0)
            mud = eng.mu_gt(ONE7 + d)
            m1d, m2d, m3d = eng.moments(ONE7 + d)
            # verify lemma (ii) exactly
            ok2 = (m1d >= m1 - 6 * d) and (m2d <= m2) and (m3d >= m3 - 18 * d * m2)
            d3s = f"{float(d3d):.4f}" if d3d is not None else " n/a "
            flag = "CLEARS" if (d3d is not None and float(d3d) >= bar) else "below"
            print(f"  d={float(d):.4f} (V/s={Vs:4.0f}): D3_d={d3s} [{flag} bar]"
                  f"  PZ_d={float(pzd):.4f}  TRUE mu_d={float(mud):.4f}"
                  f"  (ii)holds={ok2}")
        # d_max via bisection on D3_d = bar (rational bisection, 40 rounds)
        lo_, hi_ = F(0), (m1 - m2 / M) / 6   # bisect inside the valid domain only
        if d30 is not None and float(d30) > bar and hi_ > 0:
            for _ in range(40):
                mid = (lo_ + hi_) / 2
                v = D3_robust(m1, m2, m3, mid)
                if v is not None and float(v) >= bar:
                    lo_ = mid
                else:
                    hi_ = mid
            print(f"  d_max(D3_d=bar) ~ {float(lo_):.5f}  "
                  f"=> min V/s ~ {float(3/lo_):.1f}" if lo_ > 0 else
                  "  d_max ~ 0")
        print()


if __name__ == '__main__':
    main()
