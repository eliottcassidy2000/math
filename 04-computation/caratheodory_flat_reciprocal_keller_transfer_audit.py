#!/usr/bin/env python3
"""Exact/finite audit of a flat reciprocal umbilic construction and JC2 gates.

Source integrity
----------------
The formula audited here is the public announcement

    https://x.com/__alpoge__/status/2089971359921156203

not arXiv:2608.19068.  That arXiv identifier is an unrelated paper by
Brendle--Hung.  The announcement says that a PDF writeup circulated privately;
no public copy was available to this audit.

Status lanes
------------
PROVED/SYMBOLIC:
  * the seed trace-free Hessian is nowhere zero, with its exact sharp lower
    bound;
  * the conductor-order-two, degree-at-most-two-in-u nodal cell over Q has no
    Keller survivor (indeed the forced residual is -87/32);
  * algebraic identities underlying the polynomial spinor-degree gate.

FINITE-NUMERICAL:
  * k=2 spherical trace-free-Hessian winding on six specified circles, at
    2^17 equally spaced samples per circle.  Sampled nonvanishing is not a
    proof of nonvanishing between mesh points or of a unique global umbilic.

The script deliberately uses ``require`` rather than ``assert`` so that the
ordinary and ``python -O`` executions have identical validity checks.
"""

from __future__ import annotations

import math

import numpy as np
import sympy as sp


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def coeff(expr: sp.Expr, x: sp.Symbol, xp: int, u: sp.Symbol, up: int) -> sp.Expr:
    return sp.expand(expr).coeff(x, xp).coeff(u, up)


def exact_seed_audit() -> None:
    x, y = sp.symbols("x y", real=True)
    f = (
        -sp.cos(2 * x) / 4
        + sp.Rational(3, 10) * sp.cos(2 * y)
        - sp.cos(4 * y) / 32
        + sp.sin(x) * sp.sin(y)
    )
    d = sp.trigsimp(sp.diff(f, x, 2) - sp.diff(f, y, 2), method="fu")
    c = sp.trigsimp(2 * sp.diff(f, x, y), method="fu")
    d_expected = sp.cos(2 * x) + sp.Rational(6, 5) * sp.cos(2 * y) - sp.cos(4 * y) / 2
    c_expected = 2 * sp.cos(x) * sp.cos(y)
    require(sp.trigsimp(d - d_expected) == 0, "seed D identity failed")
    require(sp.trigsimp(c - c_expected) == 0, "seed C identity failed")

    # Put a=cos(x)^2 and b=cos(y)^2.  Then |4 f_{bar w bar w}|^2=M.
    a, b = sp.symbols("a b", real=True)
    d_ab = 2 * a - 4 * b**2 + sp.Rational(32, 5) * b - sp.Rational(27, 10)
    m_ab = sp.expand(d_ab**2 + 4 * a * b)
    a_star = sp.factor(-sp.diff(m_ab, a).subs(a, 0) / 8)
    require(
        sp.expand(a_star - (40 * b**2 - 74 * b + 27) / 20) == 0,
        "a minimizer failed",
    )

    # The constrained a-minimizer is a=1 on b in [0,1/10], a=a_star on
    # [1/10,1/2], and a=0 on [1/2,1].  These breakpoints follow by solving
    # a_star=1 and a_star=0.  The following identities certify each piece.
    require(
        sp.expand(
            sp.factor(a_star - 1)
            - sp.Rational(1, 20) * (10 * b - 1) * (4 * b - 7)
        )
        == 0,
        "a_star=1 breakpoints failed",
    )
    require(
        sp.expand(
            sp.factor(a_star)
            - sp.Rational(1, 20) * (2 * b - 1) * (20 * b - 27)
        )
        == 0,
        "a_star=0 breakpoints failed",
    )

    target = sp.Rational(49, 2500)
    m_at_zero = sp.factor(m_ab.subs(a, 0))
    require(
        sp.expand(
            sp.factor(m_at_zero - target)
            - sp.Rational(4, 625)
            * (5 * b - 4) ** 2
            * (100 * b**2 - 160 * b + 71)
        )
        == 0,
        "large-b lower-bound factorization failed",
    )
    require(
        sp.expand(
            sp.factor(100 * b**2 - 160 * b + 71)
            - (100 * (b - sp.Rational(4, 5)) ** 2 + 7)
        )
        == 0,
        "positive quadratic certificate failed",
    )

    m_at_star = sp.factor(m_ab.subs(a, a_star))
    require(
        sp.expand(m_at_star - b * (5 * b - 3) * (8 * b - 9) / 5) == 0,
        "middle piece failed",
    )
    require(
        sp.expand(
            sp.factor(sp.diff(m_at_star, b))
            - sp.Rational(3, 5) * (4 * b - 1) * (10 * b - 9)
        )
        == 0,
        "middle derivative factorization failed",
    )
    middle_values = tuple(
        sp.factor(m_at_star.subs(b, q))
        for q in (sp.Rational(1, 10), sp.Rational(1, 4), sp.Rational(1, 2))
    )
    require(min(middle_values) == sp.Rational(1, 4), "middle minimum failed")

    # On the small-b piece a=1.  If b>=49/10000, the term 4ab alone is
    # target-sized.  Below that cutoff D is increasing and remains far from 0.
    cutoff = sp.Rational(49, 10000)
    d_at_one = sp.factor(d_ab.subs(a, 1))
    require(
        sp.expand(sp.diff(d_at_one, b) - (-8 * b + sp.Rational(32, 5))) == 0,
        "D derivative failed",
    )
    d_cutoff = sp.factor(d_at_one.subs(b, cutoff))
    require(d_cutoff < 0 and d_cutoff**2 > target, "small-b rational bound failed")

    # Equality is attained: a=0, b=4/5 is realizable by cos(x)=0 and
    # cos(y)^2=4/5.  Since 4 f_{bar w bar w}=D+iC, divide by four.
    require(m_at_zero.subs(b, sp.Rational(4, 5)) == target, "sharp point failed")

    print("[PROVED/SYMBOLIC] seed spinor")
    print("D=f_xx-f_yy=cos(2x)+(6/5)cos(2y)-(1/2)cos(4y)")
    print("C=2f_xy=2cos(x)cos(y)")
    print("min_R2 (D^2+C^2)=49/2500; hence min |f_barwbarw|=7/200")
    print("seed phase has winding zero on every closed plane loop")


def exact_hessian_spinor_identity() -> None:
    a, b, c, mu = sp.symbols("a b c mu")
    q_plus = a - c + 2 * sp.I * b
    q_minus = a - c - 2 * sp.I * b
    identity = sp.expand(q_plus * q_minus - ((a + c) ** 2 - 4 * (a * c - b**2)))
    require(identity == 0, "Hessian spinor factorization failed")
    specialized = sp.expand(q_plus * q_minus - ((a + c) ** 2 - 4 * mu))
    require(sp.expand(specialized.subs(a * c, mu + b**2)) == 0, "constant-Hessian specialization failed")
    print("[PROVED/SYMBOLIC] polynomial Hessian spinor gate")
    print("q_+ q_-=(Delta H)^2-4 det(Hess H)")
    print("circle Laurent support for deg(q)<=D lies in [-D,D], so |wind(q)|<=D")
    print("therefore a degree-d Hessian carrier has |line-index|<=(d-2)/2")
    print("matching index 1+k/2 requires d>=k+4 (necessary, not sufficient)")


def exact_quadratic_i2_nodal_no_go() -> None:
    """All-Q no-go in the quadratic-u conductor-square THM-3586 cell."""

    u, x = sp.symbols("u X")
    h0, h1, h2, k0, k1, k2 = sp.symbols("h0 h1 h2 k0 k1 k2")
    z0, z1, z2, w0, w1, w2 = sp.symbols("z0 z1 z2 w0 w1 w2")
    h = h0 + h1 * u + h2 * u**2
    k = k0 + k1 * u + k2 * u**2
    z = z0 + z1 * u + z2 * u**2
    w = w0 + w1 * u + w2 * u**2

    f_aux = (
        -(3 * u**2 + 1) * (h**2 + k**2)
        - sp.Rational(3, 2) * u * k
        + sp.diff(k, u) / 4
        - sp.Rational(3, 16)
    )
    g_aux = (
        -2 * (3 * u**2 + 1) * h * k
        - sp.Rational(3, 2) * u * h
        + sp.diff(h, u) / 4
        + sp.Rational(3, 8)
    )
    p = -f_aux + 2 * u * z
    q = -sp.Rational(3, 2) * u * f_aux + (3 * u**2 - 1) * z
    r = -g_aux + 2 * u * w
    ss = -sp.Rational(3, 2) * u * g_aux + (3 * u**2 - 1) * w

    aa = u**2 - 1 + x * 2 * u * h + x**2 * p
    bb = u**3 - u + x * (3 * u**2 - 1) * h + x**2 * q
    cc = sp.Rational(1, 2) + 2 * u * k + x * r
    dd = sp.Rational(3, 4) * u + (3 * u**2 - 1) * k + x * ss

    def lop(expr: sp.Expr) -> sp.Expr:
        return (3 * x + 2) * expr + 2 * x * (x + 1) * sp.diff(expr, x)

    even = sp.expand(
        sp.diff(aa, u) * lop(dd)
        - lop(cc) * sp.diff(bb, u)
        + 2
        * x
        * (x + 1)
        * (sp.diff(cc, u) * sp.diff(bb, x) - sp.diff(aa, x) * sp.diff(dd, u))
    )
    odd = sp.expand(
        2 * (sp.diff(aa, u) * sp.diff(bb, x) - sp.diff(aa, x) * sp.diff(bb, u))
        + x * (sp.diff(cc, u) * lop(dd) - lop(cc) * sp.diff(dd, u))
    )
    require(coeff(even - 1, x, 0, u, 0) == 0, "boundary E equation failed")
    require(sp.expand(even - 1).coeff(x, 1) == 0, "first E equation failed")
    require(sp.expand(odd).coeff(x, 0) == 0, "boundary O equation failed")
    require(sp.expand(odd).coeff(x, 1) == 0, "first O equation failed")

    # Positive and hostile controls for the solved first layer.  The generic
    # p,q,r,s above make the X^0 and X^1 equations vanish identically.  If p
    # is perturbed by 1, the X^1 odd equation changes by -4(3u^2-1).
    aa_hostile = aa + x**2
    odd_hostile = sp.expand(
        2
        * (
            sp.diff(aa_hostile, u) * sp.diff(bb, x)
            - sp.diff(aa_hostile, x) * sp.diff(bb, u)
        )
        + x * (sp.diff(cc, u) * lop(dd) - lop(cc) * sp.diff(dd, u))
    )
    require(
        sp.expand(odd_hostile - odd).coeff(x, 1) == -4 * (3 * u**2 - 1),
        "p->p+1 hostile was not rejected by the first O equation",
    )

    # Compute and record the conductor period at the initial quadratic stage
    # through the closed THM-3586 coefficient pairing.  The Keller PDE will
    # prove more; the period is an independent necessary check and is not
    # inserted into the substitutions used by the elimination below.
    def beta(n: int) -> sp.Rational:
        return sp.Rational(
            (-1) ** n * 2 ** (2 * n + 1) * math.factorial(n) ** 2,
            math.factorial(2 * n + 1),
        )

    period_closed = sp.Integer(0)
    for i in range(3):
        for j in range(2):
            kij = -sp.Rational(2 * i, 2 * (i + j) + 3) * beta(i + j)
            period_closed += kij * (
                sp.expand(aa).coeff(x, i) * sp.expand(dd).coeff(x, j)
                - sp.expand(bb).coeff(x, i) * sp.expand(cc).coeff(x, j)
            )
    period_closed = sp.expand(period_closed)
    period_defect_full = sp.expand(sp.diff(period_closed, u) - 2)
    period_top_expected = sp.Rational(512, 105) * (
        -2 * h2**2 * w2
        + 4 * h2 * k2 * z2
        + 3 * k2**3
        - 2 * k2**2 * w2
    )
    require(
        sp.expand(period_defect_full).coeff(u, 7) == period_top_expected,
        "recorded quadratic period coefficient failed",
    )

    # These three definite leading forms are the new quadratic-u squeeze.
    # Their implications use the order on Q (or any ordered field), and are
    # not asserted over C.
    gate1 = sp.factor(coeff(odd, x, 3, u, 12))
    gate1_expected = -54 * (h2**4 + 6 * h2**2 * k2**2 + k2**4)
    require(sp.expand(gate1 - gate1_expected) == 0, "quadratic h,k gate failed")
    sub = {h2: 0, k2: 0}
    gate2 = sp.factor(coeff(odd.subs(sub), x, 3, u, 8))
    gate2_expected = -54 * (h1**4 + 6 * h1**2 * k1**2 + k1**4)
    require(sp.expand(gate2 - gate2_expected) == 0, "linear h,k gate failed")
    sub.update({h1: 0, k1: 0})
    gate3 = sp.factor(coeff(odd.subs(sub), x, 3, u, 6))
    require(sp.expand(gate3 - (-24 * (z2**2 + w2**2))) == 0, "quadratic gauge gate failed")
    sub.update({z2: 0, w2: 0})

    # The surviving equations are exactly the previously audited affine
    # ladder, now reached without assuming the quadratic gauges vanish.
    gate4 = sp.factor(coeff(odd.subs(sub), x, 4, u, 4))
    require(
        sp.expand(gate4 - (-30 * (3 * h0 * k0 + w1) ** 2)) == 0,
        "first collapsed affine gate failed",
    )
    sub[w1] = -3 * h0 * k0
    gate5 = sp.factor(coeff(odd.subs(sub), x, 3, u, 4))
    require(
        sp.expand(gate5 - (-6 * (3 * h0**2 + 3 * k0**2 + 2 * z1) ** 2)) == 0,
        "second collapsed affine gate failed",
    )
    sub[z1] = -sp.Rational(3, 2) * (h0**2 + k0**2)
    gate6 = sp.factor(coeff(odd.subs(sub), x, 4, u, 2))
    require(
        sp.expand(gate6 - (-sp.Rational(15, 8) * (3 * h0 + 4 * w0) ** 2)) == 0,
        "third collapsed affine gate failed",
    )
    sub[w0] = -sp.Rational(3, 4) * h0
    gate7 = sp.factor(coeff(odd.subs(sub), x, 3, u, 2))
    require(
        sp.expand(gate7 - (-sp.Rational(3, 2) * (3 * k0 + 4 * z0) ** 2)) == 0,
        "fourth collapsed affine gate failed",
    )
    sub[z0] = -sp.Rational(3, 4) * k0

    residual = sp.factor(coeff((even - 1).subs(sub), x, 2, u, 0))
    require(residual == -sp.Rational(87, 32), "terminal residual changed")

    collapsed = tuple(sp.factor(expr.subs(sub)) for expr in (p, q, r, ss))
    collapsed_expected = (
        (16 * h0**2 + 16 * k0**2 + 3) / 16,
        sp.Rational(3, 32) * ((32 * h0**2 + 32 * k0**2 + 3) * u + 8 * k0),
        (16 * h0 * k0 - 3) / 8,
        sp.Rational(3, 16) * (32 * h0 * k0 * u + 4 * h0 - 3 * u),
    )
    require(all(sp.expand(v - e) == 0 for v, e in zip(collapsed, collapsed_expected)), "collapsed cell failed")

    # Cross-check the closed coefficient period against an independent direct
    # pullback and termwise t-moment calculation on the collapsed cell.
    t = sp.symbols("t", real=True)
    xx = t**2 - 1
    yy = t * xx
    p0, q0, r0, s0 = collapsed
    a_pull = u**2 - 1 + xx * 2 * u * h0 + xx**2 * p0 + yy * (sp.Rational(1, 2) + 2 * u * k0 + xx * r0)
    b_pull = (
        u**3
        - u
        + xx * (3 * u**2 - 1) * h0
        + xx**2 * q0
        + yy * (sp.Rational(3, 4) * u + (3 * u**2 - 1) * k0 + xx * s0)
    )
    period_direct = sp.Integer(0)
    direct_integrand = sp.Poly(sp.expand(a_pull * sp.diff(b_pull, t)), t)
    for (power,), coefficient in direct_integrand.terms():
        if power % 2 == 0:
            period_direct += coefficient * sp.Rational(2, power + 1)
    period_direct = sp.expand(period_direct)
    period_from_pairing = sp.expand(period_closed.subs(sub))
    require(
        sp.expand(period_direct - period_from_pairing) == 0,
        "closed conductor pairing disagrees with direct t-moments",
    )
    period_defect = sp.factor(sp.diff(period_direct, u) - 2)
    period_expected = sp.Rational(2, 35) * (
        4 * h0**2 + 4 * h0 * k0 + 64 * k0**3 * u + 28 * k0**2 - 35
    )
    require(sp.expand(period_defect - period_expected) == 0, "period defect failed")

    print("[PROVED/SYMBOLIC over Q] quadratic-u conductor-square nodal cell")
    print("universe: deg_u(h,k,z,w)<=2; conductor order<=2; no coefficient-height cap")
    print("ordered-field gates force h2=k2=0, then h1=k1=0, then z2=w2=0")
    print("[X^3 u^12]O=-54(h2^4+6h2^2k2^2+k2^4)")
    print("[X^3 u^8]O=-54(h1^4+6h1^2k1^2+k1^4) after h2=k2=0")
    print("[X^3 u^6]O=-24(z2^2+w2^2) after h1=k1=h2=k2=0")
    print("successive square gates force w1=-3h0k0, z1=-3(h0^2+k0^2)/2,")
    print("w0=-3h0/4, z0=-3k0/4")
    print("terminal residual [X^2 u^0](E-1)=-87/32")
    print("controls: generic first-layer solve PASS; p->p+1 hostile REJECTED")
    print("period cross-check: closed coefficient pairing == direct exact t-moments")
    print("period defect=2(4h0^2+4h0k0+64k0^3u+28k0^2-35)/35")
    print("period on the collapsed cell forces k0=0, h0^2=35/4: no Q-point")
    print("scope: definite quartics/sums of squares use Q/order; no complex no-go claimed")


def build_k2_spinors():
    x, y = sp.symbols("x y", real=True)
    r2 = x**2 + y**2
    u = 100 * x / r2
    v = 100 * y / r2
    seed = (
        -sp.cos(2 * u) / 4
        + sp.Rational(3, 10) * sp.cos(2 * v)
        - sp.cos(4 * v) / 32
        + sp.sin(u) * sp.sin(v)
    )
    envelope = r2 / (1 + r2) * sp.exp(-(r2 ** sp.Rational(-1, 8)) * sp.exp(-r2))
    g = envelope * seed  # the +10^10 shift has zero trace-free Hessian
    gx = sp.diff(g, x)
    gy = sp.diff(g, y)
    euclidean_q4 = sp.diff(g, x, 2) - sp.diff(g, y, 2) - 2 * sp.I * sp.diff(g, x, y)
    spherical_q4 = euclidean_q4 + 4 * (x - sp.I * y) / (1 + r2) * (gx - sp.I * gy)
    # Q=g_zz+2*bar(z)/(1+|z|^2)g_z, and these expressions are 4Q.
    return sp.lambdify((x, y), (spherical_q4, euclidean_q4), modules="numpy", cse=True)


def winding_stats(values: np.ndarray, zero_tol: float = 1e-12):
    values = np.asarray(values, dtype=np.complex128)
    require(np.all(np.isfinite(values)), "non-finite spinor sample")
    minimum = float(np.min(np.abs(values)))
    require(minimum > zero_tol, f"sampled zero/near-zero detected: {minimum}")
    steps = np.angle(np.roll(values, -1) / values)
    winding_float = float(np.sum(steps) / (2 * np.pi))
    winding = int(round(winding_float))
    require(abs(winding_float - winding) < 2e-9, f"nonintegral sampled winding: {winding_float}")
    return winding, minimum, float(np.max(np.abs(steps)))


def finite_k2_circle_audit() -> None:
    n = 1 << 17
    theta = np.arange(n, dtype=np.float64) * (2 * np.pi / n)

    positive = np.exp(-4j * theta)
    require(winding_stats(positive)[0] == -4, "positive winding control failed")
    hostile_rejected = False
    try:
        winding_stats(np.cos(theta).astype(np.complex128), zero_tol=1e-10)
    except RuntimeError:
        hostile_rejected = True
    require(hostile_rejected, "zero-crossing hostile control was accepted")

    fn = build_k2_spinors()
    print("[FINITE-NUMERICAL] k=2 circle audit, N=131072 equally spaced samples")
    print("tensor: 4Q=(g_xx-g_yy)-2i g_xy+4 bar(z)(g_x-i g_y)/(1+|z|^2)")
    print("radius spherical_wind min|4Q_sph| max_phase_step euclidean_wind min|4Q_euc|")
    for radius in (0.3, 0.5, 1.0, 2.0, 5.0, 10.0):
        xx = radius * np.cos(theta)
        yy = radius * np.sin(theta)
        spherical, euclidean = fn(xx, yy)
        sw, smin, sstep = winding_stats(spherical)
        ew, emin, _ = winding_stats(euclidean)
        require(sw == -4, f"unexpected spherical winding at radius {radius}")
        print(f"{radius:4.1f} {sw:3d} {smin:.12e} {sstep:.12e} {ew:3d} {emin:.12e}")
    print("positive control exp(-4i theta): PASS; hostile cos(theta) zero-crossing: REJECTED")
    print("line-index convention is -wind(Q)/2, hence sampled spherical line index=2")
    print("scope: mesh nonvanishing only; unique global umbilic and convexity are NOT proved here")


def main() -> None:
    print("source=https://x.com/__alpoge__/status/2089971359921156203")
    print("source_correction=arXiv:2608.19068 is unrelated (Brendle--Hung, S2xS2 curvature)")
    exact_seed_audit()
    exact_hessian_spinor_identity()
    exact_quadratic_i2_nodal_no_go()
    finite_k2_circle_audit()


if __name__ == "__main__":
    main()
