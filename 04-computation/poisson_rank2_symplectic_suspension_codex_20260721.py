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
    aa, bb = sp.symbols("aa bb", nonzero=True)
    weighted_checks = 0
    for rr in range(1, 5):
        for ss in range(1, 5):
            gg = math.gcd(rr, ss)
            rr0, ss0 = rr // gg, ss // gg
            vv = x**rr0 * q**ss0
            RR = x * (aa - bb * x**rr * q**ss)
            for nn in range(4):
                sector_monomial = q * vv**nn
                lhs = sp.expand(sp.diff(RR, x) * sp.diff(sector_monomial, q) - sp.diff(RR, q) * sp.diff(sector_monomial, x))
                rhs = sp.expand(aa * (1 + ss0 * nn) * vv**nn - bb * (rr + 1 + ss0 * nn) * vv**(nn + gg))
                assert sp.expand(lhs - rhs) == 0
                weighted_checks += 1
    print("weighted-sector monomial identities:", weighted_checks, "PASS")
    print("top coefficient -bb*(r+1+s0*N)*c_N is nonzero: NO planar mate")
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

    ASR = poisson_power(S, R, 3)
    assert bracket(ASR, R) == 0
    CS = sp.expand(-ASR * D0 / 24)
    assert poisson_power(CS, R, 3) == 0
    Sq = sp.expand(S + CS)
    assert sp.expand(bracket(Sq, R) + poisson_power(Sq, R, 3) / 24) == 0
    print("one-step corrected Moyal relation M(Sq,R)=0 PASS")
    print("corrected Sq term count:", len(sp.Poly(Sq, *VARS).terms()))
    print("EXACT R-COLUMN: M(Dq,R)=1, M(T,R)=M(Sq,R)=0")
    print("QUANTUM R-COLUMN PASS: run --quantum-t for the coupled T-column audit")


def audit_quantum_t():
    print("HYP-8802 T-COLUMN CENTRALIZER REPAIR AUDIT")
    ASR = poisson_power(S, R, 3)
    CS = sp.expand(-ASR * D0 / 24)
    Sq = sp.expand(S + CS)
    assert sp.expand(bracket(Sq, R) + poisson_power(Sq, R, 3) / 24) == 0
    st_layers = [poisson_power(Sq, T, odd_order) for odd_order in (1, 3, 5)]
    EST = sp.expand(st_layers[0] + st_layers[1] / 24 + st_layers[2] / 1920 - 1)
    est_terms = sp.Poly(EST, *VARS).terms()
    est_momentum_degree = max(monomial[2] + monomial[3] for monomial, _ in est_terms)
    assert sp.expand(bracket(EST, R) + poisson_power(EST, R, 3) / 24) == 0
    print("initial Weyl-symbol residual terms / momentum degree:", len(est_terms), "/", est_momentum_degree)
    print("initial residual star-commutes with R: PASS")

    dd = sp.symbols("dd")
    sigma = sp.symbols("sigma")
    pp_sub = (2 - 6 * s) * dd + 3 * q**2 * (e - g)
    zz_sub = -3 * x**2 * dd + (1 + 3 * s) * (e - g) / 2
    est_adapted = sp.expand(EST.subs({p: pp_sub, z: zz_sub}, simultaneous=True))
    residual_e2 = sp.Poly(est_adapted, e).coeff_monomial(e**2)
    print("adapted residual terms / ell-degree:", len(sp.Poly(est_adapted, x, q, e).terms()), "/", sp.Poly(est_adapted, e).degree())
    print("ell^2 residual coefficient:", sp.factor(residual_e2))

    FF = -(108 * s**3 - 72 * s**2 - 93 * s + 58) / 16
    ff = sp.expand(x**6 * FF)
    C_ST_1 = sp.expand(ff * ell + bracket(ff, ell) / 2)
    assert sp.expand(bracket(C_ST_1, R) + poisson_power(C_ST_1, R, 3) / 24) == 0
    CST1_T = sp.expand(bracket(C_ST_1, T) + poisson_power(C_ST_1, T, 3) / 24)
    EST1 = sp.expand(EST + CST1_T)
    est1_adapted = sp.expand(EST1.subs({p: pp_sub, z: zz_sub}, simultaneous=True))
    assert sp.Poly(est1_adapted, e).degree() <= 1
    print("first star-central correction kills ell^2: PASS")
    post_e1 = sp.factor(sp.Poly(est1_adapted, e).coeff_monomial(e))
    reduced_e1 = sp.expand(sp.cancel(post_e1 / (x**8 * (3 * s - 2))).subs(q, sigma / x))
    reduced_poly = sp.Poly(reduced_e1, x)
    assert set(monomial[0] for monomial, _ in reduced_poly.terms()) <= {0, 3}
    AA = reduced_poly.coeff_monomial(1)
    BB = reduced_poly.coeff_monomial(x**3)

    def solve_weight_ode(rhs, weight):
        degree = sp.degree(rhs, sigma)
        unknowns = sp.symbols("k0:%d" % (degree + 1))
        trial = sum(unknowns[ii] * sigma**ii for ii in range(degree + 1))
        equation = sp.Poly(sp.expand((2 - 3 * sigma) * sp.diff(trial, sigma) + 3 * weight * trial - rhs), sigma)
        solution = sp.solve(equation.all_coeffs(), unknowns, dict=True)
        return None if not solution else sp.expand(trial.subs(solution[0]))

    K4 = solve_weight_ode(AA * sp.Rational(3, 2), 4)
    K7 = solve_weight_ode(BB * sp.Rational(3, 2), 7)
    assert K4 is not None and K7 is not None
    kk = sp.expand(x**4 * K4.subs(sigma, s) + x**7 * K7.subs(sigma, s))
    assert bracket(kk, R) == 0 and poisson_power(kk, T, 3) == 0
    EST2 = sp.expand(EST1 + bracket(kk, T))
    est2_adapted = sp.expand(EST2.subs({p: pp_sub, z: zz_sub}, simultaneous=True))
    assert sp.Poly(est2_adapted, e).degree() <= 0
    print("second central correction kills ell: PASS")
    print("constant residual:", sp.factor(est2_adapted))

    def integrate_W(rhs_xq):
        rhs_poly = sp.Poly(sp.expand(rhs_xq.subs(q, sigma / x)), x)
        answer = 0
        for (x_degree,), coefficient in rhs_poly.terms():
            weight = x_degree - 1
            if weight < 0:
                return None
            piece = solve_weight_ode(coefficient, weight)
            if piece is None:
                return None
            answer += x**weight * piece.subs(sigma, s)
        return sp.expand(answer)

    kernel_modes = (6, 12, 18)
    kernel_adapted_effects = []
    for mm in kernel_modes:
        f_kernel = sp.expand(x * R**(mm - 1))
        c_kernel = sp.expand(f_kernel * ell + bracket(f_kernel, ell) / 2)
        effect = sp.expand(bracket(c_kernel, T) + poisson_power(c_kernel, T, 3) / 24)
        effect_adapted = sp.expand(effect.subs({p: pp_sub, z: zz_sub}, simultaneous=True))
        effect_e = sp.Poly(effect_adapted, e).coeff_monomial(e)
        k_kernel = integrate_W(sp.cancel(sp.Rational(3, 2) * effect_e / (x**3 * (3 * s - 2))))
        assert k_kernel is not None
        kernel_effect = sp.expand(effect + bracket(k_kernel, T))
        kernel_adapted = sp.expand(kernel_effect.subs({p: pp_sub, z: zz_sub}, simultaneous=True))
        assert sp.Poly(kernel_adapted, e).degree() <= 0
        kernel_adapted_effects.append(kernel_adapted)
        kernel_weight_profile = sorted({monomial[0] for monomial, _ in sp.Poly(sp.expand(kernel_adapted.subs(q, sigma / x)), x, sigma).terms()})
        print("  homogeneous mode", mm, "constant weight profile:", kernel_weight_profile)
    closure_coeffs = sp.symbols("c6 c12 c18")
    closure = sp.Poly(sp.expand(est2_adapted + sum(cc * effect for cc, effect in zip(closure_coeffs, kernel_adapted_effects))).subs(q, sigma / x), x, sigma)
    solution = sp.solve(closure.coeffs(), closure_coeffs, dict=True)
    print("weight-6/12/18 homogeneous freedom closes constant residual:", bool(solution))
    print("T-COLUMN AUDIT COMPLETE")


if MODE == "--poisson":
    audit_poisson()
elif MODE == "--quantum":
    audit_quantum()
elif MODE == "--quantum-t":
    audit_quantum_t()
else:
    raise SystemExit("use --poisson, --quantum, or --quantum-t")
