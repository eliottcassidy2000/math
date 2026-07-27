#!/usr/bin/env python3
"""Independent hostile referee for THM-2468 (degree-22 B-W plane closure).

Fresh implementation, written without consulting the primary companion
script.  It re-derives everything from the THM-2411 flux system:

1. Specializes the two normalized fluxes N_1, N_2 of THM-2411 (12)--(16)
   to the B-W plane C=D=E=0, substitutes p=B/y^2, v=u/y^2, zeta=Z/y^3,
   W=lambda*B^3, and checks the quotient fluxes (5)-(6) of THM-2468.
2. Eliminates zeta by the closed linear-vs-quadratic resultant formula
   (cross-checked against sympy's resultant), extracts the content
   2^4*11^6, and compares the primitive part with the 24-term quintic
   pencil R_lambda typed from THM-2468 (9)-(10), including the exact
   boundary-square pencil identity (11).
3. Redoes the factor-shape audit with its own Newton polygon machinery:
   monotone-chain hull, primitive edge decomposition, and exhaustive
   balanced Minkowski-summand enumeration, proving every proper factor
   polygon is k*Delta and hence every candidate factor is a monic-in-p
   line or conic.  The line ideal and conic type I ideal must be unit;
   conic types II and III must contain a power of h, with h=0 excluded
   because the top cubic has constant term 231.
4. Rechecks the fixed section R(0,v)=-567*L_5(v): lambda-free, degree
   exactly five, squarefree, transverse derivative, nonzero v^5 leading
   coefficient (no extra point at infinity on p=0).
5. Rechecks the branch-parity and Riemann-Hurwitz arithmetic for the
   square lift Y^2=1/p, and the y=0 boundary of the first flux.

All truth-bearing checks raise explicit exceptions and stay live under
``python -O``.
"""

from __future__ import annotations

from itertools import product as iproduct
from math import gcd

import sympy as sp


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError("REFEREE FAILURE: " + message)


# ----------------------------------------------------------------------
# Own Newton polygon machinery (monotone chain + Minkowski enumeration).
# ----------------------------------------------------------------------

Point = tuple[int, int]


def turn(o: Point, a: Point, b: Point) -> int:
    return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])


def hull_ccw(points: set[Point]) -> list[Point]:
    """Counterclockwise convex hull, first vertex not repeated."""

    pts = sorted(points)
    if len(pts) <= 2:
        return pts
    lo: list[Point] = []
    for q in pts:
        while len(lo) >= 2 and turn(lo[-2], lo[-1], q) <= 0:
            lo.pop()
        lo.append(q)
    hi: list[Point] = []
    for q in reversed(pts):
        while len(hi) >= 2 and turn(hi[-2], hi[-1], q) <= 0:
            hi.pop()
        hi.append(q)
    return lo[:-1] + hi[:-1]


def edge_data(hull: list[Point]) -> tuple[list[Point], list[int]]:
    dirs: list[Point] = []
    lens: list[int] = []
    for i, a in enumerate(hull):
        b = hull[(i + 1) % len(hull)]
        dx, dy = b[0] - a[0], b[1] - a[1]
        g = gcd(abs(dx), abs(dy))
        dirs.append((dx // g, dy // g))
        lens.append(g)
    return dirs, lens


def balanced_subtuples(dirs: list[Point], lens: list[int]) -> set[tuple[int, ...]]:
    """All Minkowski-summand edge-length allocations (closed sub-polygons)."""

    out = set()
    for cand in iproduct(*(range(n + 1) for n in lens)):
        sx = sum(k * d[0] for k, d in zip(cand, dirs))
        sy = sum(k * d[1] for k, d in zip(cand, dirs))
        if sx == 0 and sy == 0:
            out.add(cand)
    return out


# ----------------------------------------------------------------------
# Ideal helpers.
# ----------------------------------------------------------------------


def unit_basis(gb: sp.GroebnerBasis) -> bool:
    return len(gb.polys) == 1 and gb.polys[0].as_expr() == 1


def reduces_to_zero(gb: sp.GroebnerBasis, expr: sp.Expr) -> bool:
    return sp.expand(gb.reduce(expr)[1]) == 0


def same_ideal(gens_a: list[sp.Expr], gens_b: list[sp.Expr], variables) -> bool:
    ga = sp.groebner(gens_a, *variables, order="grevlex")
    gb = sp.groebner(gens_b, *variables, order="grevlex")
    return all(reduces_to_zero(ga, g) for g in gens_b) and all(
        reduces_to_zero(gb, g.as_expr()) for g in ga.polys
    )


def cleared_numerators(equations, ring_vars) -> list[sp.Expr]:
    out = []
    for eq in equations:
        num = sp.together(sp.expand(eq)).as_numer_denom()[0]
        num = sp.expand(num)
        if num == 0:
            continue
        out.append(sp.Poly(num, *ring_vars).primitive()[1].as_expr())
    return out


def main() -> None:
    y, u, z = sp.symbols("y u Z")
    B, C, D, E, W = sp.symbols("B C D E W")
    p, v, zeta, lam = sp.symbols("p v zeta lambda")
    a, b, c, d, e, h, t = sp.symbols("a b c d e h t")

    # ------------------------------------------------------------------
    # (1) THM-2411 flux system, eqs (12)-(16), transcribed in full.
    # ------------------------------------------------------------------
    cal_A = 616 * B - 1089 * u + 63 * y**2
    cal_K = (
        -745360 * B * u * y
        + 6160 * B * y**3
        + 2342560 * C * u
        - 58080 * C * y**2
        + 511104 * D * y
        - 3748096 * E
        + 922383 * u**2 * y
        - 25410 * u * y**3
        + 63 * y**5
    )
    N1 = 1331 * cal_A * z + 4 * cal_K
    N2 = (
        15944049 * z**2
        + 65591680 * B * z * y
        - 206145280 * C * z
        - 162339408 * z * u * y
        + 2236080 * z * y**3
        + 1443016960 * B * u**2
        - 71554560 * B * u * y**2
        + 98560 * B * y**4
        + 449771520 * C * u * y
        - 1239040 * C * y**3
        - 1978994688 * D * u
        + 16355328 * D * y**2
        - 239878144 * E * y
        - 1319329792 * W
        - 1190488992 * u**3
        + 147581280 * u**2 * y**2
        - 1219680 * u * y**4
        + 672 * y**6
    )

    # Weighted homogeneity audit: wt(y,u,Z,B,C,D,E,W)=(1,2,3,2,3,4,5,6).
    weights = {y: 1, u: 2, z: 3, B: 2, C: 3, D: 4, E: 5, W: 6}
    gens = (y, u, z, B, C, D, E, W)
    wt_vec = [weights[g] for g in gens]
    for name, poly, wt in (("N1", N1, 5), ("N2", N2, 6)):
        for mono, coeff in sp.Poly(poly, *gens).terms():
            total = sum(m * w for m, w in zip(mono, wt_vec))
            require(
                total == wt,
                f"{name} term {mono} has weight {total}, expected {wt}",
            )

    # B-W plane: C=D=E=0, W=lambda*B^3, quotient chart.
    plane = {C: 0, D: 0, E: 0}
    N1_bw = N1.subs(plane)
    N2_bw = N2.subs(plane)
    chart = {B: p * y**2, u: v * y**2, z: zeta * y**3, W: lam * p**3 * y**6}
    f1 = sp.expand(N1_bw.subs(chart) / y**5)
    f2 = sp.expand(N2_bw.subs(chart) / y**6)
    require(not f1.has(y) and not f2.has(y), "quotient fluxes failed to descend")

    f1_claim = (
        -2981440 * p * v
        + 819896 * p * zeta
        + 24640 * p
        + 3689532 * v**2
        - 1449459 * v * zeta
        - 101640 * v
        + 83853 * zeta
        + 252
    )
    f2_claim = (
        -1319329792 * lam * p**3
        + 1443016960 * p * v**2
        - 71554560 * p * v
        + 65591680 * p * zeta
        + 98560 * p
        - 1190488992 * v**3
        + 147581280 * v**2
        - 162339408 * v * zeta
        - 1219680 * v
        + 15944049 * zeta**2
        + 2236080 * zeta
        + 672
    )
    require(sp.expand(f1 - f1_claim) == 0, "THM-2468 eq (5) f_1 mismatch")
    require(sp.expand(f2 - f2_claim) == 0, "THM-2468 eq (6) f_2 mismatch")

    wall = 616 * p - 1089 * v + 63
    require(
        sp.expand(sp.Poly(f1, zeta).coeff_monomial(zeta) - 1331 * wall) == 0,
        "zeta coefficient of f_1 is not 1331*(616p-1089v+63)",
    )

    # ------------------------------------------------------------------
    # (2) Elimination: closed resultant of (linear, quadratic) in zeta.
    # ------------------------------------------------------------------
    lin = sp.Poly(f1, zeta)
    qua = sp.Poly(f2, zeta)
    require(lin.degree() == 1 and qua.degree() == 2, "zeta degrees changed")
    a1, b1 = lin.all_coeffs()
    A2, B2, C2 = qua.all_coeffs()
    direct_res = sp.expand(A2 * b1**2 - a1 * B2 * b1 + a1**2 * C2)
    require(
        sp.expand(direct_res - sp.resultant(f1, f2, zeta)) == 0,
        "closed resultant disagrees with sympy resultant",
    )

    # THM-2468 (9)-(10), typed independently from the theorem statement.
    L5 = (
        155624547606 * v**5
        + 3215383215 * v**4
        - 1700698560 * v**3
        + 58124770 * v**2
        - 855470 * v
        + 2583
    )
    R = (
        -31289225347072 * lam * p**5
        + (110629761048576 * lam * v - 6400068820992 * lam) * p**4
        + (
            -97788806641152 * lam * v**2
            + 11314407379968 * lam * v
            - 327276246528 * lam
            + 34222590223360 * v**2
            + 3959638538240 * v
            - 44411530240
        )
        * p**3
        + (
            -149234938081152 * v**3
            - 9500102156160 * v**2
            + 695599766400 * v
            - 6017413248
        )
        * p**2
        + (
            206782580709936 * v**4
            + 6246495741024 * v**3
            - 1509756494400 * v**2
            + 34466937120 * v
            - 193496688
        )
        * p
        - 567 * L5
    )
    R = sp.expand(R)
    require(
        sp.expand(direct_res - 2**4 * 11**6 * R) == 0,
        "Res_zeta(f1,f2) != 2^4*11^6*R_lambda: eq (8)/(9) mismatch",
    )
    require(len(sp.Poly(R, p, v, lam).terms()) == 24, "R is not the 24-term quintic")
    require(sp.Poly(R, p).degree() == 5 and sp.Poly(R, p, v).total_degree() == 5,
            "R degree changed")

    # Boundary-square pencil identity (11).
    R0 = R.subs(lam, 0)
    require(
        sp.expand(R - R0 + 82458112 * lam * p**3 * wall**2) == 0,
        "pencil identity R_lambda = R_0 - 82458112*lambda*p^3*A^2 fails",
    )
    # Structural cross-check of the pencil constant: lambda enters the
    # resultant only through C2, with square prefactor a1^2=(1331*A)^2,
    # and 1331^2*1319329792/(2^4*11^6) = 82458112.
    require(
        sp.Rational(1331**2 * 1319329792, 2**4 * 11**6) == 82458112,
        "pencil constant 82458112 arithmetic fails",
    )

    # ------------------------------------------------------------------
    # (3) Uniform absolute irreducibility for lambda != 0.
    # ------------------------------------------------------------------
    # Top homogeneous quintic and its claimed factorization (12).
    top5 = sum(
        coeff * p ** m[0] * v ** m[1]
        for m, coeff in sp.Poly(R, p, v).terms()
        if m[0] + m[1] == 5
    )
    top_claim = -38974342 * (56 * p - 99 * v) ** 2 * (
        256 * lam * p**3 - 280 * p * v**2 + 231 * v**3
    )
    require(sp.expand(top5 - top_claim) == 0, "top-form factorization (12) fails")
    require(
        sp.expand(top5.subs(v, 0) + 31289225347072 * lam * p**5) == 0,
        "top form at v=0 (13) fails",
    )

    # Vertex coefficients of the Newton triangle: p^5 needs lambda != 0,
    # while v^5 and the constant vertex are lambda-free and nonzero.
    Rpv = sp.Poly(R, p, v)
    c50 = Rpv.coeff_monomial(p**5)
    c05 = Rpv.coeff_monomial(v**5)
    c00 = Rpv.coeff_monomial(1)
    require(sp.expand(c50 + 31289225347072 * lam) == 0, "p^5 vertex changed")
    require(c05 == -88239118492602, "v^5 vertex changed")
    require(c00 == -567 * 2583, "constant vertex changed")

    # Own Newton polygon audit at two generic specializations.
    for probe in (1, 17):
        support = {
            m for m, coeff in sp.Poly(R.subs(lam, probe), p, v).terms() if coeff != 0
        }
        hull = hull_ccw(support)
        require(
            set(hull) == {(0, 0), (5, 0), (0, 5)},
            f"Newton hull at lambda={probe} is not the full 5-simplex",
        )
        dirs, lens = edge_data(hull)
        summands = balanced_subtuples(dirs, lens)
        # Every Minkowski summand of 5*Delta must be k*Delta: all three
        # edge lengths equal.  This forces every proper factor to be a
        # monic-in-p line (k=1) or conic (k=2) up to constants.
        require(
            summands == {(k, k, k) for k in range(6)},
            "Minkowski summands of the 5-simplex are not k*Delta",
        )

    # Line factor audit: polygon 1*Delta => normalized p - a*v - b.
    line_rem = sp.expand(R.subs(p, a * v + b))
    line_eqs = [coeff for _, coeff in sp.Poly(line_rem, v).terms()]
    gb_line = sp.groebner(line_eqs, a, b, lam, order="grevlex")
    require(unit_basis(gb_line), "line factor ideal is not unit")

    # Conic factor audit: polygon 2*Delta => normalized
    # F = p^2 + (a*v+b)*p + c*v^2 + d*v + e.
    conic = p**2 + (a * v + b) * p + c * v**2 + d * v + e
    conic_rem = sp.rem(sp.Poly(R, p), sp.Poly(conic, p)).as_expr()
    conic_eqs = [coeff for _, coeff in sp.Poly(conic_rem, p, v).terms()]

    # The conic top p^2 + a*p*v + c*v^2 must divide the top quintic,
    # i.e. select two of the five linear factors of (12).  Three types.

    # Type I: wall squared, (p - 99v/56)^2.
    sub_I = {a: -sp.Rational(99, 28), c: sp.Rational(9801, 3136)}
    require(
        sp.expand(
            (p - sp.Rational(99, 56) * v) ** 2
            - (p**2 + a * p * v + c * v**2).subs(sub_I)
        )
        == 0,
        "type-I top expansion fails",
    )
    gb_I = sp.groebner([eq.subs(sub_I) for eq in conic_eqs], b, d, e, lam,
                       order="grevlex")
    require(unit_basis(gb_I), "conic type-I ideal is not unit")

    # Parametrize the cubic roots: p = h*v with 256*lam*h^3 - 280*h + 231 = 0,
    # equivalently lam = (280h-231)/(256h^3), legal iff h != 0.  h = 0 is
    # never a root because the cubic has constant term 231 != 0, and
    # lam != 0 iff 280h - 231 != 0.
    lam_h = (280 * h - 231) / (256 * h**3)
    cubic_at_zero = (256 * lam * h**3 - 280 * h + 231).subs(h, 0)
    require(cubic_at_zero == 231, "cubic constant term changed; h=0 exclusion fails")
    # Deflation check: dividing the cubic by (X - h) after lam -> lam_h
    # leaves remainder 0 and quotient with root-sum -h and root-product
    # -231*h^2/(280h-231) (Vieta for the complementary pair).
    X = sp.symbols("X")
    cubic_X = sp.together(256 * lam_h * X**3 - 280 * X + 231)
    quo_X, rem_X = sp.div(
        sp.numer(cubic_X), sp.expand((X - h) * sp.denom(cubic_X)), X
    )
    require(sp.simplify(rem_X) == 0, "h is not a root of the deflated cubic")
    pair_sum = sp.cancel(-sp.Poly(quo_X, X).all_coeffs()[1] / sp.Poly(quo_X, X).all_coeffs()[0])
    pair_prod = sp.cancel(sp.Poly(quo_X, X).all_coeffs()[2] / sp.Poly(quo_X, X).all_coeffs()[0])
    require(sp.cancel(pair_sum + h) == 0, "complementary root sum is not -h")
    require(
        sp.cancel(pair_prod + 231 * h**2 / (280 * h - 231)) == 0,
        "complementary root product is not -231h^2/(280h-231)",
    )

    # Type II: one wall line and one cubic-root line:
    # top = (p - 99v/56)(p - h v) => a = -(99/56 + h), c = 99h/56.
    sub_II = {lam: lam_h, a: -(sp.Rational(99, 56) + h), c: sp.Rational(99, 56) * h}
    eqs_II = cleared_numerators(
        [eq.subs(sub_II) for eq in conic_eqs], (b, d, e, h)
    )
    gb_II = sp.groebner(eqs_II, b, d, e, h, order="grevlex")
    require(not unit_basis(gb_II), "type-II ideal unexpectedly unit (theorem says h=0)")
    require(reduces_to_zero(gb_II, h), "type-II ideal does not force h = 0")
    claim_II = [
        7744 * b**2 - 1584 * b + 81,
        88 * b * e - 9 * e,
        e**2,
        d,
        h,
    ]
    match_II = same_ideal(eqs_II, claim_II, (b, d, e, h))
    require(match_II, "type-II ideal differs from THM-2468 eq (18)")

    # Type III: two cubic-root lines, omitting root h:
    # top = p^2 + h*p*v - 231h^2/(280h-231) v^2 (Vieta, verified above).
    sub_III = {lam: lam_h, a: h, c: -231 * h**2 / (280 * h - 231)}
    eqs_III = cleared_numerators(
        [eq.subs(sub_III) for eq in conic_eqs], (b, d, e, h)
    )
    gb_III = sp.groebner(eqs_III, b, d, e, h, order="grevlex")
    require(not unit_basis(gb_III), "type-III ideal unexpectedly unit")
    h_power = None
    for k in range(1, 9):
        if reduces_to_zero(gb_III, h**k):
            h_power = k
            break
    require(h_power is not None, "type-III ideal does not force h = 0")
    claim_III = [
        b**2 - e,
        b * d,
        d**2,
        b * e,
        d * e,
        e**2,
        2 * b * h - d,
        d * h,
        e * h,
        h**2 - 11 * d,
    ]
    match_III = same_ideal(eqs_III, claim_III, (b, d, e, h))
    require(match_III, "type-III ideal differs from THM-2468 eq (20)")

    # ------------------------------------------------------------------
    # (4) Fixed section R(0,v) = -567*L_5(v), squarefree, transverse.
    # ------------------------------------------------------------------
    section = sp.expand(R.subs(p, 0))
    require(not section.has(lam), "p=0 section moves with lambda")
    require(sp.expand(section + 567 * L5) == 0, "p=0 section is not -567*L_5")
    L5p = sp.Poly(L5, v)
    require(L5p.degree() == 5, "L_5 degree changed")
    require(L5p.LC() == 155624547606, "L_5 leading coefficient changed")
    require(
        sp.gcd(L5p, sp.Poly(sp.diff(L5, v), v)).degree() == 0,
        "L_5 is not squarefree",
    )
    disc = sp.discriminant(L5, v)
    require(disc != 0, "L_5 discriminant vanishes")
    require(
        sp.expand(sp.diff(R, v).subs(p, 0) + 567 * sp.diff(L5, v)) == 0,
        "transverse derivative identity at p=0 fails",
    )
    # No extra p=0 point at infinity: v^5 coefficient of the section is
    # nonzero, so the projective line p=0 meets the curve only at the
    # five affine zeros of L_5, each smooth with p a local uniformizer.
    require(
        sp.Poly(section, v).LC() == -88239118492602,
        "p=0 section top coefficient changed",
    )

    # ------------------------------------------------------------------
    # (5) Square-lift parity and genus floor; y=0 boundary.
    # ------------------------------------------------------------------
    # 1/p has odd valuation (exactly 1) at the five simple section
    # places; the number of odd valuations of a principal divisor pulled
    # into a quadratic extension is even, hence at least 6 branch places.
    visible_odd_places = 5
    branch_floor = visible_odd_places + (visible_odd_places % 2)
    require(branch_floor == 6, "branch parity floor arithmetic fails")
    # Riemann-Hurwitz with g(C) >= 0:
    # 2g(D)-2 = 2(2g(C)-2) + deg Ram >= -4 + 6 = 2  =>  g(D) >= 2.
    two_gD_minus_2_floor = 2 * (2 * 0 - 2) + branch_floor
    require(two_gD_minus_2_floor == 2, "Riemann-Hurwitz floor arithmetic fails")
    genus_floor = (two_gD_minus_2_floor + 2) // 2
    require(genus_floor == 2, "genus floor is not 2")
    # Physical square coordinate: p*Y^2 = 1 under Y=y/b0, b0^2=B, p=B/y^2.
    b0, Y = sp.symbols("b0 Y", nonzero=True)
    require(
        sp.cancel((p * Y**2).subs({p: b0**2 / y**2, Y: y / b0})) == 1,
        "physical square-lift identity fails",
    )
    # y=0 boundary of the first flux on the plane, THM-2468 (29).
    require(
        sp.expand(N1_bw.subs(y, 0) - 1331 * (616 * B - 1089 * u) * z) == 0,
        "y=0 boundary (29) fails",
    )

    # ------------------------------------------------------------------
    # Transcript.
    # ------------------------------------------------------------------
    print("THM-2468 independent hostile referee (opus, 2026-07-26)")
    print("flux_source=THM-2411_eqs_12-16_transcribed_with_C_D_E_W")
    print("weighted_homogeneity=N1_wt5_N2_wt6_checked")
    print("quotient_fluxes=MATCH_THM-2468_eqs_5_6")
    print("elimination=closed_linear_quadratic_resultant+sympy_cross_check")
    print("resultant_content=2^4*11^6")
    print("R_lambda=MATCH_24_term_quintic_eq_9_10")
    print("pencil_identity=R0_minus_82458112_lambda_p3_A2_EXACT")
    print("pencil_constant_check=1331^2*1319329792/(2^4*11^6)=82458112")
    print("top_form_factorization=eq_12_EXACT")
    print("newton_hull=full_5_simplex_at_lambda_1_and_17")
    print("minkowski_summands=k_Delta_only_(k=0..5)")
    print("line_ideal=UNIT")
    print("conic_type_I_ideal=UNIT")
    print(f"conic_type_II_ideal=CONTAINS_h_matches_eq_18={match_II}")
    print(
        "conic_type_III_ideal="
        f"CONTAINS_h^{h_power}_matches_eq_20={match_III}"
    )
    print("h_zero_excluded=cubic_constant_term_231")
    print("p0_section=-567*L5_lambda_free")
    print("L5_squarefree=YES_discriminant_nonzero")
    print("p0_smooth_simple_places=5")
    print("p0_point_at_infinity=NONE_v5_coeff_-88239118492602")
    print("branch_floor=6")
    print("genus_floor=2")
    print("square_lift=physical_Y=y/sqrt(B)")
    print("y0_boundary=1331*(616B-1089u)*Z")
    print("ALL INDEPENDENT CHECKS PASSED")


if __name__ == "__main__":
    main()
