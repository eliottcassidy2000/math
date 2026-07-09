#!/usr/bin/env python3
"""
LEM-014 verification: the P-separated composed realization (wide regime).
boxeph-2026-07-09-S1.  Pure Python, exact rationals (fractions.Fraction).

For a covering 13-set S = P u L (P = small part <= 13, L = cluster),
Vmax = max L, co-offsets E = {Vmax - u} (0 in E), s = spread(E) = max(E):

  1. Build the exact cell decomposition of x in [0,1): cells are open
     intervals on which every phase frac(e_i x) is linear and their sorted
     order is constant.  Breakpoints: m/e_i and m/|e_i - e_j|.
  2. Compute exactly:  mu_theta(E) = meas{maxgap > theta}  for theta = 1/7
     and theta = 1/7 + delta;  E[W] = int_0^1 W,  W = Sum (g_i - 1/7)_+.
  3. Build G_P^eps = {x : ||p x|| >= 1/14 + eps  for all p in P} exactly.
  4. Robust set R = {x in G_P^eps : maxgap > 1/7 + delta}; exact measure;
     pick x* = midpoint of a longest interval of R.
  5. Composed construction: j = round(Vmax x*), teeth t_i = frac(e_i j/Vmax),
     phi* = midpoint of the largest circular gap of the teeth,
     tau = (j + phi*)/Vmax.  Verify EXACTLY: ||v tau|| >= 1/14 for all v in S.
  6. Report chain margins (gap survival, drift budget, P-erosion budget).

delta = 3 s / Vmax, eps = 20 / Vmax  (the LEM-014 constants).
"""
from fractions import Fraction as F
import sys, itertools

ONE7 = F(1, 7)
ONE14 = F(1, 14)


def frac(x):
    return x - (x.numerator // x.denominator)


def norm_dist(x):
    """||x|| = distance to nearest integer, exact."""
    f = frac(x)
    return min(f, 1 - f)


# ---------------------------------------------------------------- cells
def breakpoints(E):
    """All x in [0,1] where some phase hits 0 or two phases collide."""
    pts = {F(0), F(1)}
    for e in E:
        if e > 0:
            for m in range(e + 1):
                pts.add(F(m, e))
    for a, b in itertools.combinations(E, 2):
        d = abs(a - b)
        if d > 0:
            for m in range(d + 1):
                pts.add(F(m, d))
    return sorted(pts)


def cell_gaps(E, lo, hi):
    """On the open cell (lo,hi): return list of gaps as (const, slope) with
    g(x) = const + slope*x, for the circular gaps of {frac(e_i x)}."""
    mid = (lo + hi) / 2
    ph = []  # (value_at_mid, const, slope)
    for e in E:
        c = (e * mid).numerator // (e * mid).denominator  # floor(e*mid)
        ph.append((e * mid - c, F(-c), F(e)))
    ph.sort()
    k = len(ph)
    gaps = []
    for i in range(k - 1):
        _, c1, s1 = ph[i]
        _, c2, s2 = ph[i + 1]
        gaps.append((c2 - c1, s2 - s1))
    # wraparound gap: 1 - top + bottom
    _, cb, sb = ph[0]
    _, ct, st = ph[-1]
    gaps.append((1 + cb - ct, sb - st))
    return gaps


def solve_gt(c, s, theta, lo, hi):
    """{x in (lo,hi) : c + s x > theta} as interval or None (exact)."""
    if s == 0:
        return (lo, hi) if c > theta else None
    r = (theta - c) / s
    if s > 0:
        a, b = max(lo, r), hi
    else:
        a, b = lo, min(hi, r)
    return (a, b) if a < b else None


def merge(intervals):
    """Merge sorted-or-not interval list into disjoint sorted list."""
    if not intervals:
        return []
    iv = sorted(intervals)
    out = [list(iv[0])]
    for a, b in iv[1:]:
        if a <= out[-1][1]:
            out[-1][1] = max(out[-1][1], b)
        else:
            out.append([a, b])
    return [(a, b) for a, b in out]


def measure(intervals):
    return sum(b - a for a, b in intervals)


def intersect(A, B):
    """Intersect two disjoint sorted interval lists."""
    out, i, j = [], 0, 0
    while i < len(A) and j < len(B):
        a = max(A[i][0], B[j][0])
        b = min(A[i][1], B[j][1])
        if a < b:
            out.append((a, b))
        if A[i][1] < B[j][1]:
            i += 1
        else:
            j += 1
    return out


class ClusterEngine:
    def __init__(self, E):
        self.E = sorted(E)
        assert 0 in E, "cluster must contain the observer co-offset 0"
        bps = breakpoints(self.E)
        self.cells = [(bps[i], bps[i + 1]) for i in range(len(bps) - 1)
                      if bps[i] < bps[i + 1]]
        self._gaps = [cell_gaps(self.E, lo, hi) for lo, hi in self.cells]

    def mu_gt(self, theta):
        """meas{ x : maxgap > theta } and the interval list, exact."""
        ivs = []
        for (lo, hi), gaps in zip(self.cells, self._gaps):
            cell_ivs = []
            for c, s in gaps:
                r = solve_gt(c, s, theta, lo, hi)
                if r:
                    cell_ivs.append(r)
            ivs.extend(merge(cell_ivs))
        ivs = merge(ivs)
        return measure(ivs), ivs

    def EW(self):
        """int_0^1 Sum_i (g_i - 1/7)_+ dx, exact."""
        tot = F(0)
        for (lo, hi), gaps in zip(self.cells, self._gaps):
            for c, s in gaps:
                r = solve_gt(c, s, ONE7, lo, hi)
                if not r:
                    continue
                a, b = r
                # integrate (c + s x - 1/7) over (a,b): linear, positive there
                tot += (c - ONE7) * (b - a) + s * (b * b - a * a) / 2
        return tot


def gp_eps(P, eps):
    """G_P^eps = {x in [0,1] : ||p x|| >= 1/14 + eps for all p} exact."""
    cur = [(F(0), F(1))]
    thr = ONE14 + eps
    if thr >= F(1, 2):
        return []
    for p in P:
        good = []
        for m in range(p):
            a, b = F(m + thr, p) if False else ((m + thr) / p, (m + 1 - thr) / p)
            good.append((F(a), F(b)))
        cur = intersect(cur, merge(good))
        if not cur:
            return []
    return cur


def compose_tau(E, Vmax, xstar):
    """The LEM-014 construction: j = round(Vmax x*), phi* = gap midpoint,
    tau = (j + phi*)/Vmax.  Returns (tau, j, phi*, grid_gap)."""
    j = (Vmax * xstar + F(1, 2)).numerator // (Vmax * xstar + F(1, 2)).denominator
    teeth = sorted(frac(F(e * j, Vmax)) for e in E)
    k = len(teeth)
    # largest circular gap
    best_len, best_mid = F(-1), None
    for i in range(k):
        a = teeth[i]
        b = teeth[(i + 1) % k] + (1 if i == k - 1 else 0)
        if b - a > best_len:
            best_len, best_mid = b - a, frac((a + b) / 2)
    phi = best_mid
    tau = (j + phi) / Vmax
    return tau, j, phi, best_len


def verify(S, tau):
    """Exact min clearance over all speeds."""
    return min(norm_dist(v * tau) for v in S)


def is_covering(S):
    return all(any(v % q == 0 for v in S) for q in range(2, 15))


def run_instance(name, P, E, Vmax, scan_note=""):
    P, E = sorted(P), sorted(E)
    S = sorted(P + [Vmax - e for e in E])
    k, s = len(E), max(E)
    assert len(S) == 13, f"{name}: need 13 speeds, got {len(S)}"
    cov = is_covering(S)
    delta = F(3 * s, Vmax)
    eps = F(20, Vmax)
    eng = ClusterEngine(E)
    mu0, _ = eng.mu_gt(ONE7)
    mud, ivd = eng.mu_gt(ONE7 + delta)
    ew = eng.EW()
    gp = gp_eps(P, eps)
    R = intersect(ivd, gp) if P else ivd
    mR = measure(R)
    print(f"\n=== {name} ===  {scan_note}")
    print(f"  S = {S}")
    print(f"  covering: {cov} | P={P} k={k} s={s} Vmax={Vmax} "
          f"Vmax/s={Vmax/s:.1f}")
    print(f"  delta=3s/V={float(delta):.6f}  eps=20/V={float(eps):.6f}")
    print(f"  mu_1/7(E)          = {mu0} = {float(mu0):.6f}")
    print(f"  mu_1/7+delta(E)    = {mud} = {float(mud):.6f}   "
          f"shell={float(mu0-mud):.6f}")
    print(f"  E[W]               = {ew} = {float(ew):.6f}")
    print(f"  one-line bound (7/6)(E[W]-6d) = {float(F(7,6)*(ew-6*delta)):.6f}"
          f"   (vs actual mu_robust {float(mud):.6f})")
    print(f"  meas(G_P^eps)      = {float(measure(gp)):.6f}" if P else
          "  P empty: G_P = [0,1]")
    print(f"  meas(R robust set) = {mR} = {float(mR):.6f}")
    if mR == 0:
        print("  ROBUST SET EMPTY -> construction not attempted "
              "(outside wide regime or floor fails)")
        return None
    # x* = midpoint of longest interval of R
    a, b = max(R, key=lambda iv: iv[1] - iv[0])
    xstar = (a + b) / 2
    tau, j, phi, grid_gap = compose_tau(E, Vmax, xstar)
    clr = verify(S, tau)
    ok = clr >= ONE14
    print(f"  x* = {xstar} (interval len {float(b-a):.6g})")
    print(f"  j = {j}   grid teeth maxgap = {float(grid_gap):.6f} "
          f"(need > 1/7 + drift: {float(ONE7 + 2*F(s,Vmax)):.6f})")
    print(f"  phi* = {float(phi):.6f}   tau = {tau}")
    print(f"  min clearance = {clr} = {float(clr):.8f}  "
          f"{'>= 1/14 OK' if ok else '< 1/14 FAIL'}   "
          f"margin = {float(clr - ONE14):+.8f}")
    return ok, clr


def find_covering_vmax(P, E, start):
    """Smallest Vmax >= start making S = P u {Vmax - e} covering."""
    V = start
    while True:
        S = list(P) + [V - e for e in E]
        if len(set(S)) == 13 and min(S) > 0 and is_covering(S):
            return V
        V += 1


def main():
    print("LEM-014 composed realization - exact verification "
          "(boxeph-2026-07-09-S1)")

    # ---- 1. tight-AP cluster, embedded WIDE (k=13, P empty) --------------
    E_ap = list(range(13))
    V = find_covering_vmax([], E_ap, 372 * 12 + 3540)   # ~ 8000
    run_instance("AP{0..12} wide", [], E_ap, V)

    # ---- 2. the 7-structured hard cluster (MISTAKE-128), WIDE ----------
    E_7s = [0, 7, 14, 21, 26, 29, 37, 44, 51, 58, 67, 75, 82]
    V = find_covering_vmax([], E_7s, 372 * 82 + 3540)   # ~ 34000
    run_instance("7-structured (M128) wide", [], E_7s, V)

    # ---- 3. knife-edge shape (MISTAKE-130), WIDE + covering ------------
    E_ke = [0, 7, 10, 14, 18, 20, 21, 26, 28, 35, 36, 37, 42]
    V = find_covering_vmax([], E_ke, 372 * 42 + 3540)   # ~ 19000
    run_instance("knife-edge (M130) wide", [], E_ke, V)

    # ---- 4. P nonempty: k=10, P = {11,12,13} ---------------------------
    E_10 = [0, 1, 3, 7, 12, 20, 30, 43, 65, 82]     # scattered, s=82
    V = find_covering_vmax([11, 12, 13], E_10, 372 * 82 + 3540)
    run_instance("k=10 P={11,12,13} wide", [11, 12, 13], E_10, V)

    # ---- 5. P nonempty: k=8, P = {8, 11, 12, 13, 9} --------------------
    E_8 = [0, 2, 5, 11, 17, 25, 36, 42]
    V = find_covering_vmax([7, 9, 11, 12, 13], E_8, 372 * 42 + 3540)
    run_instance("k=8 P={7,9,11,12,13} wide", [7, 9, 11, 12, 13], E_8, V)

    # ---- 6. empirical wide/compressed boundary on the AP shape ---------
    print("\n=== Vmax descent scan: AP{0..12}, P empty -- where does the "
          "composed construction first fail? ===")
    eng = ClusterEngine(E_ap)
    lastfail, firstok_streak = None, None
    for V in [8000, 4000, 2000, 1000, 500, 250, 125, 64, 48, 32, 24, 16, 14]:
        if not is_covering([V - e for e in E_ap]):
            V0 = find_covering_vmax([], E_ap, V)
        else:
            V0 = V
        delta, eps = F(3 * 12, V0), F(20, V0)
        mud, ivd = eng.mu_gt(ONE7 + delta)
        if mud == 0:
            print(f"  Vmax={V0:6d} (V/s={V0/12:6.1f}): robust set EMPTY")
            continue
        a, b = max(ivd, key=lambda iv: iv[1] - iv[0])
        xstar = (a + b) / 2
        tau, j, phi, gg = compose_tau(E_ap, V0, xstar)
        clr = verify([V0 - e for e in E_ap], tau)
        print(f"  Vmax={V0:6d} (V/s={V0/12:6.1f}): mu_rob={float(mud):.4f}  "
              f"min clr={float(clr):.6f}  "
              f"{'OK' if clr >= ONE14 else 'FAIL'} "
              f"(margin {float(clr-ONE14):+.6f})")

    # ---- 7. compressed sanity: the real 7-structured @91 (in scope) ----
    # covering-derived cluster at Vmax=91 (mac-mini/klein), all 13 co-offsets
    E91 = [0, 7, 14, 21, 26, 29, 37, 44, 51, 58, 67, 75, 82]
    S91 = [91 - e for e in E91]
    if is_covering(S91):
        run_instance("7-structured @91 (compressed)", [], E91, 91,
                     "expected outside wide regime")
    else:
        print("\n7-structured @91: NOT covering as bare cluster "
              "(needs its own P) -- skipped")


if __name__ == '__main__':
    sys.setrecursionlimit(10000)
    main()
