#!/usr/bin/env python3
r"""
lrc14_cubic_moment_gate_monad_S11.py   (monad-explorer-2026-07-07-S11, HYP-5117)

THE CUBIC MOMENT GATE -- build and stress-test mac-mini-S49's named-but-never-run tool:
the joint (Sum g^2, Sum g^3) moment bound on the Bad event, as a PER-SHAPE EXACT
certificate machine for the k=8..13 tail floors (A').

FRAME.  The load-bearing (A') is mu_{1/7}(E) >= T_k for every single-scale k-cluster E
(bars: klein ledger 0.675 at k=8 ... 0.056 at k=13; kps T_8 = 0.6185).  Equivalently
P(Bad_E) <= 1 - T_k where Bad_E = {x : all k gaps of {frac(e x)} <= 1/7}.

  * S49's 2-moment Cantelli on Y = Sum g^2 - 1/8: iid floor 0.667, census worst 0.4754
    -- SHORT.  Named next tool: joint (Sum g^2, Sum g^3).  Never run.
  * monad-S10 degree gap: the LRC deficit's quadratic part vanishes identically; leading
    degree 3.  PREDICTION: the second moment (a quadratic instrument) must fail and the
    third moment is the first that can grip.
  * HYP-4987: any proof must consume weight->=3 relation structure.  The moment route
    does (E[Y] already consumes k-body hole correlations; E[Y^2] 4-tuples; E[Y^3] 6-tuples).
  * HYP-3789 (forgotten factoid): the covering-min is a truncated moment problem with
    atomic extremals (Curto-Fialkow).  This is the SAME mathematics aimed at the tail side.

MACHINE.
  1. Exact order-cell engine: for integer offsets E (|E| = k), the circular gap vector
     g(x) is piecewise affine in x with breakpoints m/d, d in {|e_i - e_j|} u {e_i}.
     Per cell every gap is A x + B (A, B integers).  Exact rational:
       - moments  E[Y^m], Y = Sum g^2 - 1/8   (polynomial integration, Fractions)
       - mixed    E[Y^a Z^b], Z = Sum g^3 - 1/64 (centered at uniform value 8*(1/8)^3)
       - P(Bad)   (per-cell interval intersection, exact; validates vs known mu values)
       - exact per-shape range [Ymin, Ymax] of Y (support for the moment LP)
  2. Bounds on P(Y <= ythr), ythr = 1/56 (Bad => Y <= 1/7 - 1/8 = 1/56):
       - 2-moment one-sided Cantelli (reproduces S49)
       - sharp 3-moment bound: dual LP over cubics q with q >= 1 on [0, ythr],
         q >= 0 on [ythr, Ymax]; P <= sum c_j m_j.  Grid LP + continuous verification.
       - joint (Y, Z) linear-combination bound: for lam >= 0, Bad => Y + lam*Zc <= b(lam);
         apply the 1D machinery to W = Y + lam*Zc (moments expand from mixed moments).
  3. Adversarial battery (MISTAKE-119 discipline: structured families + JUMP descent,
     not local nudges): AP_8, near-APs, parity interlace, two-block, single-far,
     geometric, Sidon, random, + annealing directly on the certificate value.

HONESTY.  Per-shape certificates only; the uniform-over-shapes statement is NOT proved
here.  The output to watch: does the 3-moment (or joint) bound clear the bar at the
shapes where the 2-moment bound fails -- and if not, WHERE exactly does it lose
(event relaxation vs moment truncation)?  Loss decomposition reported per shape.
"""

from fractions import Fraction
from itertools import combinations
from math import comb, factorial, sqrt
import random

F = Fraction
ONE7 = F(1, 7)

# ----------------------------------------------------------------------------------
# 1. EXACT ORDER-CELL ENGINE
# ----------------------------------------------------------------------------------

def breakpoints(E):
    """All x in [0,1] where the circular order/wrapping of {frac(e x)} can change."""
    ds = set()
    for a, b in combinations(E, 2):
        ds.add(abs(a - b))
    for e in E:
        if e != 0:
            ds.add(abs(e))
    bps = {F(0), F(1)}
    for d in ds:
        for m in range(1, d + 1):
            bps.add(F(m, d))
    return sorted(bps)


def cell_gaps(E, x0, x1):
    """Gap vector as affine forms (A, B) meaning g = A x + B, valid on (x0, x1).
    Uses the midpoint to fix floors and the circular order (constant on the cell)."""
    xm = (x0 + x1) / 2
    k = len(E)
    # floor offsets constant on the cell
    fl = [ (e * xm).__floor__() for e in E ]
    # phase_i(x) = e_i x - fl_i ; sort by value at midpoint
    idx = sorted(range(k), key=lambda i: E[i] * xm - fl[i])
    gaps = []
    for r in range(k):
        i, j = idx[r], idx[(r + 1) % k]
        A = E[j] - E[i]
        B = F(fl[i] - fl[j])
        if r == k - 1:  # wrap gap
            B += 1
        gaps.append((F(A), B))
    return gaps


def poly_mul(p, q):
    r = [F(0)] * (len(p) + len(q) - 1)
    for i, a in enumerate(p):
        if a:
            for j, b in enumerate(q):
                if b:
                    r[i + j] += a * b
    return r


def poly_add(p, q):
    n = max(len(p), len(q))
    return [ (p[i] if i < len(p) else F(0)) + (q[i] if i < len(q) else F(0))
             for i in range(n) ]


def poly_int(p, x0, x1):
    """Integrate polynomial (coeff list, ascending) from x0 to x1, exactly."""
    tot = F(0)
    xp0, xp1 = x0, x1
    for i, c in enumerate(p):
        if c:
            tot += c * (xp1 - xp0) / (i + 1)
        xp0 *= x0
        xp1 *= x1
    return tot


def poly_eval(p, x):
    v = F(0)
    for c in reversed(p):
        v = v * x + c
    return v


class ShapeData:
    """Exact per-shape data: moments of (Y, Z), P(Bad), Y-range."""
    __slots__ = ("E", "k", "mom", "pbad", "ymin", "ymax", "p_y_le")

    def __init__(self, E, ythr=F(1, 56), max_mixed=3):
        E = sorted(set(int(e) for e in E))
        self.E = E
        self.k = len(E)
        k = self.k
        bps = breakpoints(E)
        # Y = sum g^2 - 1/8 ; Zc = sum g^3 - 1/64  (k=8 uniform values; general: 1/k, 1/k^2)
        y_shift = F(1, k)
        z_shift = F(1, k * k)
        mom = {}          # (a,b) -> E[Y^a Zc^b] for a + b <= max_mixed
        pbad = F(0)
        p_y_le = F(0)     # exact P(Y <= ythr)  (quadratic root isolation per cell)
        ymin, ymax = None, None
        for ci in range(len(bps) - 1):
            x0, x1 = bps[ci], bps[ci + 1]
            if x0 == x1:
                continue
            gaps = cell_gaps(E, x0, x1)
            # Y(x) = sum (A x + B)^2 - 1/k  as poly
            Yp = [F(-1, 1) * y_shift]
            Zp = [F(-1, 1) * z_shift]
            for (A, B) in gaps:
                Yp = poly_add(Yp, [B * B, 2 * A * B, A * A])
                Zp = poly_add(Zp, [B**3, 3 * A * B * B, 3 * A * A * B, A**3])
            # accumulate mixed moments
            Ypow = {0: [F(1)]}
            for a in range(1, max_mixed + 1):
                Ypow[a] = poly_mul(Ypow[a - 1], Yp)
            Zpow = {0: [F(1)]}
            for b in range(1, max_mixed + 1):
                Zpow[b] = poly_mul(Zpow[b - 1], Zp)
            for a in range(max_mixed + 1):
                for b in range(max_mixed + 1 - a):
                    if a == b == 0:
                        continue
                    mom[(a, b)] = mom.get((a, b), F(0)) + poly_int(
                        poly_mul(Ypow[a], Zpow[b]), x0, x1)
            # P(Bad): all gaps <= 1/7 -> interval intersection
            lo, hi = x0, x1
            for (A, B) in gaps:
                # A x + B <= 1/7
                if A == 0:
                    if B > ONE7:
                        lo, hi = x1, x0  # empty
                        break
                elif A > 0:
                    hi = min(hi, (ONE7 - B) / A)
                else:
                    lo = max(lo, (ONE7 - B) / A)
            if hi > lo:
                pbad += hi - lo
            # exact Y-range on the cell (quadratic: check endpoints + vertex)
            c0 = poly_eval(Yp, x0); c1 = poly_eval(Yp, x1)
            cand = [c0, c1]
            if len(Yp) >= 3 and Yp[2] != 0:
                xv = -Yp[1] / (2 * Yp[2])
                if x0 < xv < x1:
                    cand.append(poly_eval(Yp, xv))
            mn, mx = min(cand), max(cand)
            ymin = mn if ymin is None else min(ymin, mn)
            ymax = mx if ymax is None else max(ymax, mx)
            # exact-ish P(Y <= ythr): quadratic sublevel per cell
            p_y_le += _quad_sublevel_measure(Yp, ythr, x0, x1)
        self.mom = mom
        self.pbad = pbad
        self.ymin, self.ymax = ymin, ymax
        self.p_y_le = p_y_le

    # convenience
    def m(self, a, b=0):
        if a == 0 and b == 0:
            return F(1)
        return self.mom[(a, b)]


def _quad_sublevel_measure(Yp, thr, x0, x1):
    """Measure of {x in [x0,x1] : Y(x) <= thr} for quadratic Y. Exact rational when
    roots rational; else uses Fraction-approximated sqrt (denominator 10^12) --
    reported quantity is diagnostic, not certificate."""
    a = Yp[2] if len(Yp) >= 3 else F(0)
    b = Yp[1] if len(Yp) >= 2 else F(0)
    c = Yp[0] - thr
    if a == 0:
        if b == 0:
            return (x1 - x0) if c <= 0 else F(0)
        r = -c / b
        if b > 0:
            lo, hi = x0, min(x1, r)
        else:
            lo, hi = max(x0, r), x1
        return max(F(0), hi - lo)
    disc = b * b - 4 * a * c
    if disc <= 0:
        return (x1 - x0) if a < 0 else F(0)
    sq = _frac_sqrt(disc)
    r1 = (-b - sq) / (2 * a)
    r2 = (-b + sq) / (2 * a)
    r1, r2 = min(r1, r2), max(r1, r2)
    if a > 0:   # <= thr between roots
        lo, hi = max(x0, r1), min(x1, r2)
        return max(F(0), hi - lo)
    else:       # <= thr outside roots
        return max(F(0), min(x1, r1) - x0) + max(F(0), x1 - max(x0, r2))


def _frac_sqrt(fr, prec=10**12):
    from math import isqrt
    num, den = fr.numerator, fr.denominator
    s = isqrt(num * den * prec * prec)
    return F(s, den * prec)


# ----------------------------------------------------------------------------------
# 1b. EXCESS-MASS FUNCTIONALS  V_theta(x) = sum_i (g_i(x) - theta)_+
#     V_theta = 0 EXACTLY on Bad (all gaps <= 1/7 <= theta) -- no event relaxation.
#     Exact moments and cross-moments via cell subdivision at gap-crossing kinks.
# ----------------------------------------------------------------------------------

def excess_moments(E, thetas, max_pow=3):
    """Exact E[V_a], E[V_a V_b], E[V_1^3] (first theta) for V_t = sum_i (g_i - t)_+.
    Returns (a, M, m3) with a[i] = E[V_i] (Fraction), M[i][j] = E[V_i V_j], m3 = E[V_0^3],
    plus vmax[i] = exact max of V_i over x."""
    E = sorted(set(int(e) for e in E))
    bps = breakpoints(E)
    nT = len(thetas)
    a = [F(0)] * nT
    M = [[F(0)] * nT for _ in range(nT)]
    m3 = F(0)
    vmax = [F(0)] * nT
    for ci in range(len(bps) - 1):
        x0, x1 = bps[ci], bps[ci + 1]
        if x0 == x1:
            continue
        gaps = cell_gaps(E, x0, x1)
        # subdivision points: g_r(x) = theta crossings inside the cell
        cuts = {x0, x1}
        for (A, B) in gaps:
            if A != 0:
                for t in thetas:
                    xc = (t - B) / A
                    if x0 < xc < x1:
                        cuts.add(xc)
        cuts = sorted(cuts)
        for si in range(len(cuts) - 1):
            u0, u1 = cuts[si], cuts[si + 1]
            um = (u0 + u1) / 2
            # V_t affine on (u0,u1): V = P_t x + Q_t
            V = []
            for t in thetas:
                P, Q = F(0), F(0)
                for (A, B) in gaps:
                    if A * um + B > t:
                        P += A
                        Q += B - t
                V.append((P, Q))
            for i in range(nT):
                Pi, Qi = V[i]
                a[i] += poly_int([Qi, Pi], u0, u1)
                vmax[i] = max(vmax[i], Pi * u0 + Qi, Pi * u1 + Qi)
                for j in range(i, nT):
                    Pj, Qj = V[j]
                    M[i][j] += poly_int(poly_mul([Qi, Pi], [Qj, Pj]), u0, u1)
            P0, Q0 = V[0]
            m3 += poly_int(poly_mul(poly_mul([Q0, P0], [Q0, P0]), [Q0, P0]), u0, u1)
    for i in range(nT):
        for j in range(i):
            M[i][j] = M[j][i]
    return a, M, m3, vmax


def rayleigh_floor(a, M):
    """Exact rational max_c (c.a)^2 / (c^T M c) = a^T M^{-1} a (if M invertible):
    mu >= this value.  Gaussian elimination over Fractions; returns (floor, coeffs)."""
    n = len(a)
    # solve M c = a
    A = [row[:] + [a[i]] for i, row in enumerate(M)]
    # forward elimination with partial pivot (exact)
    piv = []
    r = 0
    for c in range(n):
        pr = None
        for rr in range(r, n):
            if A[rr][c] != 0:
                pr = rr
                break
        if pr is None:
            continue
        A[r], A[pr] = A[pr], A[r]
        pv = A[r][c]
        A[r] = [v / pv for v in A[r]]
        for rr in range(n):
            if rr != r and A[rr][c] != 0:
                f = A[rr][c]
                A[rr] = [v - f * w for v, w in zip(A[rr], A[r])]
        piv.append(c)
        r += 1
        if r == n:
            break
    c = [F(0)] * n
    for r_i, c_i in enumerate(piv):
        c[c_i] = A[r_i][n]
    val = sum(ci * ai for ci, ai in zip(c, a))
    return val, c   # val = a^T M^{-1} a  (Fraction)


def atom_zero_bound_3mom(m1, m2, m3, vmax):
    """Sharp max of P(V = 0) given exact (m1,m2,m3) of V supported on [0, vmax]:
    q = 1 - p must admit a measure on (0, vmax] with moments (q, m1, m2, m3).
    Conditions: q*m2 >= m1^2   and   (q*vmax - m1)(m2*vmax - m3) >= (m1*vmax - m2)^2.
    Returns exact Fraction upper bound on P(V=0)."""
    if m2 == 0:
        return F(1)
    q1 = m1 * m1 / m2
    den = m2 * vmax - m3
    if den <= 0:
        q2 = F(0)
    else:
        q2 = ((m1 * vmax - m2) ** 2 / den + m1) / vmax
    q = max(q1, q2)
    return 1 - min(q, F(1))


# ----------------------------------------------------------------------------------
# 2. MOMENT BOUNDS
# ----------------------------------------------------------------------------------

def cantelli2(m1, m2, thr):
    """One-sided Cantelli: P(Y <= thr) <= Var/(Var + (E[Y]-thr)^2) when E[Y] > thr."""
    var = m2 - m1 * m1
    if m1 <= thr:
        return 1.0
    v = float(var); d = float(m1 - thr)
    return v / (v + d * d)


def dual_cubic_bound(moms, thr, ymax, ngrid=4000):
    """Sharp-ish 3-moment bound on P(Y <= thr), support [0, ymax]:
    min c.m  s.t.  q(y) = c0 + c1 y + c2 y^2 + c3 y^3 >= 1 on [0, thr], >= 0 on [thr, ymax].
    Grid LP (scipy linprog) + continuous a-posteriori verification with margin report.
    moms = (1, m1, m2, m3) floats. Returns (bound, verified_ok, worst_violation)."""
    import numpy as np
    from scipy.optimize import linprog
    thr = float(thr); ymax = float(ymax)
    ys1 = np.linspace(0.0, thr, max(60, ngrid // 10))
    ys2 = np.linspace(thr, ymax, ngrid)
    # constraints: -q(y) <= -1 on ys1 ; -q(y) <= 0 on ys2
    def vand(ys):
        return np.stack([np.ones_like(ys), ys, ys ** 2, ys ** 3], axis=1)
    A_ub = np.vstack([-vand(ys1), -vand(ys2)])
    b_ub = np.concatenate([-np.ones(len(ys1)), np.zeros(len(ys2))])
    cvec = np.array([float(m) for m in moms])
    res = linprog(cvec, A_ub=A_ub, b_ub=b_ub, bounds=[(None, None)] * 4,
                  method="highs")
    if not res.success:
        return (1.0, False, None)
    c = res.x
    # continuous verification: min of q-1 on [0,thr], min of q on [thr,ymax]
    qd = np.polynomial.polynomial.polyder(c)
    crit = [r.real for r in np.polynomial.polynomial.polyroots(qd)
            if abs(r.imag) < 1e-12]
    def qv(y):
        return ((c[3] * y + c[2]) * y + c[1]) * y + c[0]
    v1 = min([qv(0), qv(thr)] + [qv(r) for r in crit if 0 < r < thr]) - 1.0
    v2 = min([qv(thr), qv(ymax)] + [qv(r) for r in crit if thr < r < ymax])
    worst = min(v1, v2)
    bound = float(np.dot(cvec, c))
    # tiny violations from grid LP -> inflate bound by exact compensation:
    # q + |worst| is feasible, costs + |worst| * m0
    if worst < 0:
        bound += -worst
    return (min(bound, 1.0), worst > -1e-9, worst)


def joint_bound(sd, lam, bar_thr=None, ngrid=4000):
    """Linear-combination joint bound: W = Y + lam * Zc.
    Bad => all g <= 1/7 => Y <= 1/56 and Zc <= sup(sum g^3 - 1/64 | g <= 1/7, sum g = 1).
    For k=8: max sum g^3 under g_i <= 1/7 is 7 * (1/7)^3 = 1/49; Zc <= 1/49 - 1/64.
    So Bad => W <= 1/56 + lam * (1/49 - 1/64) =: wthr  (lam >= 0).
    Moments of W from mixed moments. 3-moment dual-cubic on W."""
    k = sd.k
    y_bad = F(1, 7) - F(1, k)                      # max Y on Bad
    z_bad = F(1, 49) - F(1, k * k) if k == 8 else F(k, 7**3) - F(1, k * k)
    lamF = F(lam).limit_denominator(10**6)
    wthr = y_bad + lamF * z_bad
    # W moments up to 3
    def mW(n):
        tot = F(0)
        for a in range(n + 1):
            tot += comb(n, a) * lamF ** (n - a) * sd.m(a, n - a)
        return tot
    m1, m2, m3 = mW(1), mW(2), mW(3)
    # support: W >= Wmin; shift to [0, Wmax - Wmin]
    # crude exact Wmax/Wmin: use per-shape Y-range and |Zc| <= max(|z|) bound via ranges:
    # do it numerically -- sample cells? For the dual LP support use conservative
    # [wmin_est, wmax_est] from coarse sampling; conservativeness checked by margin.
    wmin, wmax = _w_range_numeric(sd, float(lamF))
    sh = F(wmin).limit_denominator(10**9)
    moms = (1.0, float(m1 - sh), float(m2 - 2 * sh * m1 + sh * sh),
            float(m3 - 3 * sh * m2 + 3 * sh * sh * m1 - sh ** 3))
    if m1 <= wthr:
        return 1.0
    b, ok, worst = dual_cubic_bound(moms, float(wthr - sh), float(F(wmax) - sh) * 1.0001 + 1e-9,
                                    ngrid=ngrid)
    return b if ok or worst > -1e-7 else 1.0


def _w_range_numeric(sd, lam, ngrid=20011):
    """Numeric range of W = Y + lam Zc over x (fine grid; support padding applied by caller)."""
    E = sd.E; k = sd.k
    wmin, wmax = None, None
    for r in range(ngrid):
        x = (r + 0.5) / ngrid
        ph = sorted((e * x) % 1.0 for e in E)
        gs = [ph[i + 1] - ph[i] for i in range(k - 1)] + [ph[0] + 1 - ph[-1]]
        y = sum(g * g for g in gs) - 1.0 / k
        z = sum(g ** 3 for g in gs) - 1.0 / k ** 2
        w = y + lam * z
        wmin = w if wmin is None else min(wmin, w)
        wmax = w if wmax is None else max(wmax, w)
    return wmin - 1e-6, wmax + 1e-6


# ----------------------------------------------------------------------------------
# 3. BASELINES AND BATTERY
# ----------------------------------------------------------------------------------

def dirichlet_moments(k, max_order=3):
    """Exact iid baseline: E[(sum g^2)^m] for spacings ~ Dirichlet(1^k), then centered.
    E[prod g_i^{a_i}] = prod(a_i!) (k-1)! / (k-1+sum a)!."""
    # (sum g^2)^m = sum over multisets; expand via compositions of m into parts,
    # each part j meaning exponent 2j on a distinct coordinate.
    from itertools import product as iproduct
    def esym(m):
        # sum over ways to write the multinomial: iterate partitions of m into parts
        total = F(0)
        for part in _partitions(m):
            # part = list of part sizes (each part -> one coordinate with exponent 2*size)
            r = len(part)
            if r > k:
                continue
            # multinomial coefficient m! / prod(part!) / prod(mult!)
            coef = factorial(m)
            for p in part:
                coef //= factorial(p)
            from collections import Counter
            for c in Counter(part).values():
                coef //= factorial(c)
            # choose ordered distinct coordinates: k(k-1)...(k-r+1)
            nch = 1
            for t in range(r):
                nch *= (k - t)
            aa = [2 * p for p in part]
            emom = F(factorial(k - 1))
            num = 1
            for a in aa:
                num *= factorial(a)
            emom = F(num * factorial(k - 1), factorial(k - 1 + sum(aa)))
            total += coef * nch * emom
        return total
    raw = [F(1)] + [esym(m) for m in range(1, max_order + 1)]
    # center: Y = S - 1/k
    sh = F(1, k)
    cen = []
    for m in range(max_order + 1):
        v = F(0)
        for j in range(m + 1):
            v += comb(m, j) * (-sh) ** (m - j) * raw[j]
        cen.append(v)
    return cen  # [1, E[Y], E[Y^2], E[Y^3]]


def _partitions(m, mx=None):
    if mx is None:
        mx = m
    if m == 0:
        yield []
        return
    for p in range(min(m, mx), 0, -1):
        for rest in _partitions(m - p, p):
            yield [p] + rest


def certify(E, bars=(0.675, 0.6185), verbose=True, ngrid=4000):
    """Full per-shape report. Returns dict of floats."""
    sd = ShapeData(E)
    k = sd.k
    thr = F(1, 7) - F(1, k)
    m1, m2, m3 = sd.m(1), sd.m(2), sd.m(3)
    pb = float(sd.pbad)
    mu = 1 - pb
    c2 = cantelli2(float(m1), float(m2), float(thr))
    b3, ok3, worst3 = dual_cubic_bound((1.0, float(m1), float(m2), float(m3)),
                                       float(thr), float(sd.ymax) + 1e-12, ngrid=ngrid)
    # joint bound, optimize lam over a small grid
    best_joint, best_lam = b3, 0.0
    for lam in (0.5, 1.0, 2.0, 4.0, 8.0, 16.0, -0.5, -1.0):
        # NOTE: lam < 0 needs Zc's Bad-max replaced by Bad-min; skip (only lam >= 0 valid)
        if lam < 0:
            continue
        bj = joint_bound(sd, lam, ngrid=ngrid // 2)
        if bj < best_joint:
            best_joint, best_lam = bj, lam
    res = dict(E=sd.E, k=k, mu=mu, pbad=pb, p_y_le=float(sd.p_y_le),
               m1=float(m1), var=float(m2 - m1 * m1), m3=float(m3),
               ymin=float(sd.ymin), ymax=float(sd.ymax),
               cant2=1 - c2, cube3=1 - b3, joint=1 - best_joint, lam=best_lam)
    if verbose:
        line = (f"E={sd.E}  mu={mu:.4f} [P(Bad)={pb:.5f}] P(Y<=thr)={res['p_y_le']:.5f} | "
                f"EY={res['m1']:.5f} sd={sqrt(res['var']):.5f} | floors: "
                f"C2={res['cant2']:.4f} M3={res['cube3']:.4f} J={res['joint']:.4f}"
                f"(lam={best_lam}) | bars {bars}")
        print(line)
    return res


# ----------------------------------------------------------------------------------
if __name__ == "__main__":
    import sys
    print("=" * 100)
    print("PART 0 -- VALIDATION")
    print("=" * 100)
    # AP_8: known mu(AP_8) = 691/735 (fleet ledger) => P(Bad) = 44/735
    sd = ShapeData(range(8))
    print(f"AP_8 exact P(Bad) = {sd.pbad} = {float(sd.pbad):.6f}  "
          f"(expect 44/735 = {44/735:.6f})  mu = {1-float(sd.pbad):.6f}")
    print(f"AP_8 exact moments: E[Y]={sd.m(1)} = {float(sd.m(1)):.6f}, "
          f"E[Y^2]={sd.m(2)}, E[Y^3]={sd.m(3)}")
    print(f"AP_8 Y-range: [{float(sd.ymin):.6f}, {float(sd.ymax):.6f}]")
    dm = dirichlet_moments(8)
    print(f"Dirichlet(1^8) baseline: E[Y]={dm[1]} = {float(dm[1]):.6f} (expect 7/72={7/72:.6f}), "
          f"E[Y^2]={dm[2]} = {float(dm[2]):.6f} (expect 133/10560={133/10560:.6f}), "
          f"E[Y^3]={dm[3]} = {float(dm[3]):.6f}")
    thr = 1/7 - 1/8
    c2 = cantelli2(float(dm[1]), float(dm[2]), thr)
    print(f"iid 2-moment Cantelli floor: mu >= {1-c2:.4f} (S49 reported 0.6672)")
    b3, ok, worst = dual_cubic_bound(tuple(float(v) for v in dm), thr, 7/8)
    print(f"iid 3-moment dual-cubic floor: mu >= {1-b3:.4f} (verified={ok}, margin={worst:.2e})")
    sys.stdout.flush()

    print()
    print("=" * 100)
    print("PART 1 -- STRUCTURED BATTERY (k=8; bars: klein 0.675 / kps T_8 0.6185)")
    print("     floors reported: C2 = 2-moment Cantelli, M3 = 3-moment dual cubic, J = joint (Y,Z)")
    print("=" * 100)
    battery = [
        ("AP_8 (conjectured mu-min)", list(range(8))),
        ("near-AP one-swap 7->8",     [0,1,2,3,4,5,6,8]),
        ("near-AP one-swap 7->13",    [0,1,2,3,4,5,6,13]),
        ("near-AP two-defect",        [0,1,2,3,5,6,7,9]),
        ("parity interlace (M119)",   [0,2,4,6,8,10,11,12]),
        ("parity interlace b",        [0,2,4,6,7,8,10,12]),
        ("two-block 4+4 far",         [0,1,2,3,40,41,42,43]),
        ("two-block 4+4 near",        [0,1,2,3,10,11,12,13]),
        ("single-far",                [0,1,2,3,4,5,6,40]),
        ("geometric",                 [0,1,2,4,8,16,32,64]),
        ("Sidon (perfect diff)",      [0,1,3,7,12,20,30,44]),
        ("GW-like 12->24 analog",     [0,1,2,3,4,5,7,14]),
        ("dilAP d=3 + defect",        [0,3,6,9,12,15,18,22]),
        ("dense-random",              [0,2,3,7,9,13,16,17]),
    ]
    rows = []
    for name, E in battery:
        print(f"[{name}]")
        rows.append((name, certify(E)))
        sys.stdout.flush()

    print()
    print("=" * 100)
    print("PART 2 -- LOSS DECOMPOSITION at the binding shapes")
    print("  P(Bad) <= P(Y<=thr) <= moment bound.  loss1 = P(Y<=thr) - P(Bad) (event relaxation),")
    print("  loss2 = bound - P(Y<=thr) (moment truncation).")
    print("=" * 100)
    for name, r in rows:
        l1 = r['p_y_le'] - r['pbad']
        l2 = (1 - r['joint']) - r['p_y_le']
        print(f"{name:28s} P(Bad)={r['pbad']:.4f} P(Y<=thr)={r['p_y_le']:.4f} "
              f"bestbound={1-r['joint']:.4f}  loss1={l1:.4f} loss2={l2:.4f}")
    sys.stdout.flush()

    print()
    print("=" * 100)
    print("PART 3 -- ADVERSARIAL DESCENT on the joint floor (JUMP moves, MISTAKE-119 discipline)")
    print("=" * 100)
    rng = random.Random(20260707)
    def floor_of(E):
        try:
            r = certify(E, verbose=False, ngrid=1500)
            return r['joint'], r
        except Exception as ex:
            return None, None
    cur = list(range(8)); curf, curr = floor_of(cur)
    best = (curf, cur, curr)
    seen = {tuple(cur): curf}
    for step in range(120):
        cand = sorted(set(cur))
        i = rng.randrange(8)
        newv = rng.randint(0, 45)
        cand2 = sorted(set(cand[:i] + [newv] + cand[i+1:]))
        if len(cand2) != 8:
            continue
        # translate to min 0, reduce by gcd (mu invariances)
        m0 = cand2[0]; cand2 = [c - m0 for c in cand2]
        from math import gcd
        g = 0
        for c in cand2[1:]:
            g = gcd(g, c)
        if g > 1:
            cand2 = [c // g for c in cand2]
        t = tuple(cand2)
        if t in seen:
            f = seen[t]
        else:
            f, r = floor_of(cand2)
            seen[t] = f
            if f is not None and f < best[0]:
                best = (f, cand2, r)
                print(f"  step {step}: new worst joint floor {f:.4f} at {cand2}")
        # hill-descend: accept if lower floor (we hunt the WORST shape for the method)
        if f is not None and (f < curf or rng.random() < 0.25):
            cur, curf = cand2, f
    print(f"DESCENT WORST: joint floor {best[0]:.4f} at {best[1]}")
    if best[2] is not None:
        r = best[2]
        print(f"  detail: mu={r['mu']:.4f} P(Bad)={r['pbad']:.4f} EY={r['m1']:.5f} "
              f"Var={r['var']:.6f} M3floor={r['cube3']:.4f}")
    sys.stdout.flush()

    print()
    print("=" * 100)
    print("PART 4 -- THE EXCESS-MASS ROUTE: V_theta = sum (g - theta)_+ vanishes EXACTLY on Bad")
    print("  PZ = E[V]^2/E[V^2] (exact rational); AT3 = 3-moment atom-at-zero bound (exact);")
    print("  RAY = Rayleigh over multi-width family {V_theta_j} = a^T M^{-1} a (exact rational).")
    print("=" * 100)
    T0 = F(1, 7)
    widths = [F(1,7), F(31,210), F(8,49), F(6,35), F(19,105), F(1,5), F(23,105), F(1,4), F(2,7)]
    for name, E in battery:
        aV, MV, m3V, vmaxV = excess_moments(E, widths)
        pz = float(aV[0] * aV[0] / MV[0][0]) if MV[0][0] else 0.0
        at3 = 1 - float(atom_zero_bound_3mom(aV[0], MV[0][0], m3V, vmaxV[0]))
        ray, coef = rayleigh_floor(aV, MV)
        sd0 = ShapeData(E)
        print(f"{name:28s} mu={1-float(sd0.pbad):.4f} | PZ={pz:.4f} AT3={at3:.4f} "
              f"RAY={float(ray):.4f}  (bars 0.675/0.6185)")
        sys.stdout.flush()
