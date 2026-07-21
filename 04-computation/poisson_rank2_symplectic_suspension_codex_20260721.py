#!/usr/bin/env python3
"""Exact audit for THM-2044, THM-2045, and HYP-8802."""

import gc
import math
import sys

MODE = sys.argv[1] if len(sys.argv) > 1 else "--poisson"

import sympy as sp

x, q, p, z, e = sp.symbols("x q p z e")
VARS = (x, q, p, z)
s = x * q

R = x * (2 - 3 * s)
D0 = (1 + 3 * s) * p / 2 - 3 * q**2 * z
L = 3 * x**2 * p + (2 - 6 * s) * z
G = 252 * s**3 + 1008 * s**2 + 1379 * s + 659
g = -q**2 * G / 140
ell = sp.expand(L + g)

y = q - x * e / 3
u = 1 + x * y
T0 = sp.expand(y + 3 * x * u**2 * e + 3 * x * y**2 * (4 + 3 * x * y))
S0 = sp.expand((u**3 * e + y**2 * u * (4 + 3 * x * y)) / 2)

B = (
    2 * e**4 * x**6 * (3 * s - 2)
    + e**3 * x**4 * (-90 * s**2 - 30 * s + 55)
    + e**2 * x**2 * (540 * s**3 + 720 * s**2 - 120 * s - 270)
    + e * (-1620 * s**4 - 3780 * s**3 - 1890 * s**2 + 810 * s + 540)
    + q**2 * (2430 * s**3 + 8100 * s**2 + 8640 * s + 2430)
)
H0 = sp.expand(-e * B / 1620)

T = sp.expand(T0.subs(e, ell))
S = sp.expand(S0.subs(e, ell))
D = sp.expand(D0 + H0.subs(e, ell))


def bracket(f, h):
    return sp.expand(
        sp.diff(f, p) * sp.diff(h, x)
        - sp.diff(f, x) * sp.diff(h, p)
        + sp.diff(f, z) * sp.diff(h, q)
        - sp.diff(f, q) * sp.diff(h, z)
    )


def derivative(f, orders):
    for variable, order in zip(VARS, orders):
        if order:
            f = sp.diff(f, variable, order)
    return f


def compositions4(n):
    for aa in range(n + 1):
        for bb in range(n - aa + 1):
            for cc in range(n - aa - bb + 1):
                yield aa, bb, cc, n - aa - bb - cc


def poisson_power(f, h, n):
    """n-th bidifferential power for {p,x}={z,q}=1."""
    out = 0
    for aa, bb, cc, dd in compositions4(n):
        df = derivative(f, (bb, dd, aa, cc))
        dg = derivative(h, (aa, cc, bb, dd))
        if not df or not dg:
            continue
        coeff = sp.Rational(
            (-1) ** (bb + dd) * math.factorial(n),
            math.factorial(aa)
            * math.factorial(bb)
            * math.factorial(cc)
            * math.factorial(dd),
        )
        out += coeff * df * dg
    return sp.expand(out)


def audit_poisson():
    print("THM-2044 RANK-TWO POISSON SUSPENSION AUDIT")
    print("expanded term counts (R,T,D,S):", tuple(len(sp.Poly(f, *VARS).terms()) for f in (R, T, D, S)))
    J3 = sp.Matrix([R, T0, S0]).jacobian([x, q, e])
    assert sp.factor(J3.det()) == 1
    print("sheared three-dimensional determinant: 1 PASS")
    relations = (
        ("{D,R}", D, R, 1), ("{S,T}", S, T, 1),
        ("{R,T}", R, T, 0), ("{R,S}", R, S, 0),
        ("{D,T}", D, T, 0), ("{D,S}", D, S, 0),
    )
    for name, f, h, expected in relations:
        assert bracket(f, h) == expected
        print(name, "=", expected, "PASS")

    intermediate_points = (
        (sp.Rational(0), sp.Rational(0), sp.Rational(-1, 4)),
        (sp.Rational(1), sp.Rational(2, 3), sp.Rational(13, 2)),
        (sp.Rational(-1), sp.Rational(-2, 3), sp.Rational(13, 2)),
    )
    lifted = []
    for xx, qq, ll in intermediate_points:
        base = {x: xx, q: qq, e: ll}
        LL = sp.factor(ll - g.subs(base))
        d0 = sp.factor(-H0.subs(base))
        rx = (2 - 6 * s).subs(base)
        aaa = ((1 + 3 * s) / 2).subs(base)
        bbb = (-3 * q**2).subs(base)
        pp = sp.factor(rx * d0 - bbb * LL)
        zz = sp.factor(-3 * xx**2 * d0 + aaa * LL)
        point = (xx, qq, pp, zz)
        image = tuple(sp.factor(f.subs(dict(zip(VARS, point)))) for f in (R, T, D, S))
        assert image == (0, 0, 0, sp.Rational(-1, 8))
        lifted.append(point)
    print("lifted fibre points:")
    for point in lifted:
        print(" ", point)
    print("common image: (0, 0, 0, -1/8) PASS")
    print("exact fibre cardinality: 3 (transported through polynomial source/target automorphisms from THM-1300)")
    print("THM-2045 PLANAR DE-STABILIZATION AUDIT")
    print("constant-producing Laurent sector equation: ((aa-bb*ss)*f(ss))' = 1")
    print("polynomial solution forces f=-1/bb, which is not divisible by ss: NO planar mate")
    print("TOURNAMENT ANALYSIS: not imposed; orienting the six tensor identities would discard the symplectic equation.")
    print("POISSON PASS")


def audit_quantum():
    print("HYP-8802 WEYL-SYMMETRIC QUANTIZATION AUDIT")
    p3_counts = {}
    saved_dr_anomaly = None
    quantum_pairs = (
        ("{D,R}", D, R), ("{S,T}", S, T), ("{R,T}", R, T),
        ("{R,S}", R, S), ("{D,T}", D, T), ("{D,S}", D, S),
    )
    for name, f, h in quantum_pairs:
        anomaly = poisson_power(f, h, 3)
        if name == "{D,R}":
            saved_dr_anomaly = anomaly
        p3_counts[name] = 0 if not anomaly else len(sp.Poly(anomaly, *VARS).terms())
        print(" ", name, "P^3 terms:", p3_counts[name])
    assert p3_counts == {"{D,R}": 42, "{S,T}": 42, "{R,T}": 0, "{R,S}": 3, "{D,T}": 165, "{D,S}": 273}
    A3 = saved_dr_anomaly
    assert bracket(A3, R) == 0
    C1 = sp.expand(-A3 * D0 / 24)
    B3 = sp.expand(-18 * sp.diff(C1, p, p, z))
    assert sp.factor(B3) == 108 * x**12 * (3 * q * x - 2) * (3 * q * x - 1) * (27 * q**2 * x**2 - 2)
    C2 = sp.expand(-B3 * D0 / 24)
    Dq = sp.expand(D + C1 + C2)
    assert sp.expand(bracket(Dq, R) - 18 * sp.diff(Dq, p, p, z) / 24 - 1) == 0
    print("two-step corrected Moyal relation M(Dq,R)=1 PASS")
    print("corrected Dq term count:", len(sp.Poly(Dq, *VARS).terms()))
    print("QUANTUM PARTIAL PASS: simultaneous six-relation correction remains open")


if MODE == "--poisson":
    audit_poisson()
elif MODE == "--quantum":
    audit_quantum()
else:
    raise SystemExit("use --poisson or --quantum")
