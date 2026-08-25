#!/usr/bin/env python3
"""Exact certificate for THM-4120's target and DE-response ledgers.

The script checks polynomial identities and finite arithmetic.  The cited
geometric input is Shioda--Tate (or, equivalently here, the classified
II*+I1+I1 rational elliptic surface): a rational elliptic surface whose
section/fibre and II* component lattices already have rank ten and absolute
discriminant one has Mordell--Weil group {O}.

Universe
========

* k is algebraically closed of characteristic zero, a,gamma,theta are nonzero,
  K=k(q), and gamma=-a^3/2;
* the hypothetical pair lies in THM-4103's smooth nonresonant theta-only
  reduced (2,3) cell;
* the seven source-infinity points have indices AB:1, BC:2+2, CD:3,
  DE:7, EF:3+3, with only BC non-rational over K;
* THM-4053's generic degree is an Eisenstein norm;
* THM-3992's DE support depths hold.

Nothing here asserts entry of an arbitrary planar Keller counterexample into
that cell, and the script does not prove Shioda--Tate.
"""

from __future__ import annotations

import sympy as sp


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def polynomial_valuation(poly: sp.Expr, variable: sp.Symbol) -> int:
    terms = sp.Poly(sp.expand(poly), variable).terms()
    return min(monomial[0] for monomial, coefficient in terms if coefficient)


def is_eisenstein_norm(value: int) -> bool:
    """The local prime-exponent criterion for x^2-xy+y^2."""
    remaining = value
    prime = 2
    while prime * prime <= remaining:
        exponent = 0
        while remaining % prime == 0:
            remaining //= prime
            exponent += 1
        if prime % 3 == 2 and exponent % 2:
            return False
        prime += 1
    return not (remaining > 1 and remaining % 3 == 2)


def layer_solutions(weight: int, depth: int, limit: int = 8) -> list[tuple[int, int]]:
    """Solutions 7m-6n=weight, m<=n+depth, in increasing total degree."""
    solutions = []
    for m in range(0, 4 * weight + 80):
        for n in range(0, 4 * weight + 80):
            if 7 * m - 6 * n == weight and m <= n + depth:
                solutions.append((m, n))
    return sorted(solutions, key=lambda pair: (sum(pair), pair))[:limit]


def main() -> None:
    a, gamma, theta, phi, q, u, r = sp.symbols(
        "a gamma theta phi q u r", nonzero=True
    )
    x, t, w = sp.symbols("x t w")

    # Target invariants and its three singular fibres.
    a4 = -sp.Rational(3, 4) * a**2
    a6 = q - sp.Rational(1, 4) * a**3
    c4 = sp.expand(-48 * a4)
    c6 = sp.expand(-864 * a6)
    discriminant = sp.factor(-16 * (4 * a4**3 + 27 * a6**2))
    require(c4 == 36 * a**2, "finite c4 identity changed")
    require(c6 == -864 * q + 216 * a**3, "finite c6 identity changed")
    require(
        sp.expand(discriminant + 432 * q * (q - a**3 / 2)) == 0,
        "target discriminant identity changed",
    )

    # q=1/u and U=u^-2 X,V=u^-3 Y give the minimal integral infinity model.
    a4_inf = sp.expand(a4 * u**4)
    a6_inf = sp.expand(u**5 - sp.Rational(1, 4) * a**3 * u**6)
    c4_inf = sp.expand(-48 * a4_inf)
    c6_inf = sp.expand(-864 * a6_inf)
    delta_inf = sp.factor(-16 * (4 * a4_inf**3 + 27 * a6_inf**2))
    infinity_triple = (
        polynomial_valuation(c4_inf, u),
        polynomial_valuation(c6_inf, u),
        polynomial_valuation(delta_inf, u),
    )
    require(infinity_triple == (4, 5, 10), "infinity is no longer II*")

    # Exact Shioda--Tate/discriminant bookkeeping.  The section/fibre block is
    # odd unimodular with Gram [[-1,1],[1,0]], not the even hyperbolic plane.
    section_fibre_gram = sp.Matrix([[-1, 1], [1, 0]])
    require(section_fibre_gram.det() == -1, "section/fibre block changed")
    rho_rational_surface = 10
    e8_rank = 8
    trivial_rank = section_fibre_gram.rows + e8_rank
    trivial_abs_discriminant = abs(int(section_fibre_gram.det())) * 1
    ns_abs_discriminant = 1
    mw_rank = rho_rational_surface - trivial_rank
    torsion_order_squared = trivial_abs_discriminant // ns_abs_discriminant
    require((mw_rank, torsion_order_squared) == (0, 1), "MW ledger changed")

    # Boundary residue fields and the target-infinity response.
    rational_punctures = (
        ("AB", 1),
        ("CD", 3),
        ("DE", 7),
        ("EF+", 3),
        ("EF-", 3),
    )
    rational_weight = sum(index for _, index in rational_punctures)
    require(rational_weight == 17, "K-rational boundary weight changed")
    pre_norm_degrees = (rational_weight, rational_weight + 2 + 2)
    norm_degrees = tuple(n for n in pre_norm_degrees if is_eisenstein_norm(n))
    require(pre_norm_degrees == (17, 21), "pre-norm degree pair changed")
    require(norm_degrees == (21,), "degree 17 was not excluded")
    pullback_O = (1, 2, 2, 3, 7, 3, 3)
    require(sum(pullback_O) == 21, "pullback divisor degree changed")
    pole_degrees = (2 * sum(pullback_O), 3 * sum(pullback_O))
    require(pole_degrees == (42, 63), "coordinate pole degrees changed")

    # Sharp quadratic-field hostile: the target need not remain point-free
    # after the BC base change.
    q_quadratic = a**3 / 2 + r**2
    hostile_residual = sp.expand(
        r**2
        - ((a / 2) ** 3 + a4 * (a / 2) + q_quadratic - a**3 / 4)
    )
    require(hostile_residual == 0, "quadratic target point hostile disappeared")

    # DE edge expansion.  With z=1/s and w=s^6 t=x^6 t^7,
    # z^8 H=theta+phi*z+(epsilon+kappa)*z^2+O(z^3).
    epsilon, kappa, z = sp.symbols("epsilon kappa z")
    T = -gamma / theta
    w_series = sp.series(
        -gamma / (theta + phi * z + (epsilon + kappa) * z**2), z, 0, 3
    ).removeO()
    c1 = sp.expand(w_series).coeff(z, 1)
    c2 = sp.expand(w_series).coeff(z, 2)
    require(sp.simplify(c1 - gamma * phi / theta**2) == 0, "DE c1 changed")
    require(
        sp.simplify(c2 - gamma * (epsilon + kappa) / theta**2
                    + gamma * phi**2 / theta**3) == 0,
        "DE c2 changed",
    )
    A5 = sp.symbols("A5", nonzero=True)
    epsilon_forced = sp.Rational(2752, 135) * gamma / A5**3
    kappa_forced = -sp.Rational(5696, 45) * gamma / A5**3
    sum_forced = sp.factor(epsilon_forced + kappa_forced)
    require(
        sum_forced == -sp.Rational(14336, 135) * gamma / A5**3,
        "phi=0 second re-entry coefficient changed",
    )

    # Direct top-layer bracket audit with generic degree-two/three f,g.
    f0, f1, f2, g0, g1, g2 = sp.symbols("f0 f1 f2 g0 g1 g2")
    f = f0 + f1 * w + f2 * w**2
    g = g0 + g1 * w + g2 * w**2
    f_xt = f.subs(w, x**6 * t**7)
    g_xt = g.subs(w, x**6 * t**7)
    A14 = x**2 * f_xt
    C21 = x**3 * g_xt
    bracket = sp.expand(sp.diff(A14, x) * sp.diff(C21, t)
                        - sp.diff(A14, t) * sp.diff(C21, x))
    expected_bracket = sp.expand(
        7 * x**10 * t**6
        * (2 * f * sp.diff(g, w) - 3 * sp.diff(f, w) * g).subs(
            w, x**6 * t**7
        )
    )
    require(sp.expand(bracket - expected_bracket) == 0, "DE bracket changed")

    # A second, non-baseline layer checks the general fixed-weight identity
    # rather than inferring it from the square/cube chamber.
    aa, bb, cc, dd = 9, 8, 10, 8
    DA, DC = 7 * aa - 6 * bb, 7 * cc - 6 * dd
    A_general = x**aa * t**bb * f_xt
    C_general = x**cc * t**dd * g_xt
    bracket_general = sp.expand(
        sp.diff(A_general, x) * sp.diff(C_general, t)
        - sp.diff(A_general, t) * sp.diff(C_general, x)
    )
    coefficient_general = (
        (aa * dd - bb * cc) * f * g
        + DA * w * f * sp.diff(g, w)
        - DC * w * sp.diff(f, w) * g
    )
    expected_general = sp.expand(
        x ** (aa + cc - 1) * t ** (bb + dd - 1)
        * coefficient_general.subs(w, x**6 * t**7)
    )
    require(
        sp.expand(bracket_general - expected_general) == 0,
        "general fixed-layer bracket changed",
    )
    require((DA, DC) == (15, 22), "general layer weights changed")

    # Fixed-layer normal forms and the first possible higher layers.
    # The first entries are THM-3992's two explicit exceptions to the depth
    # inequalities; every later entry is governed by those inequalities.
    A14_solutions = [(2, 0)] + layer_solutions(14, 1, 3)
    C21_solutions = [(3, 0)] + layer_solutions(21, 2, 3)
    A15_solutions = layer_solutions(15, 1, 3)
    C22_solutions = layer_solutions(22, 2, 3)
    require(A14_solutions == [(2, 0), (8, 7), (14, 14), (20, 21)],
            "A weight-14 layer changed")
    require(C21_solutions == [(3, 0), (9, 7), (15, 14), (21, 21)],
            "C weight-21 layer changed")
    require(A15_solutions[:2] == [(9, 8), (15, 15)],
            "first higher A layer changed")
    require(C22_solutions[:2] == [(10, 8), (16, 15)],
            "first higher C layer changed")

    # In the (D_A,D_C)=(14,21) cell the UFD equation forces
    # f=gamma^2 h^2,g=gamma^3 h^3,h(0)=1.  The unique simultaneous target-O
    # response is h(T)=-7.  A linear h realizes the minimal support floors.
    h = 1 - 8 * w / T
    require(sp.simplify(h.subs(w, 0) - 1) == 0, "h(0) normalization changed")
    require(sp.simplify(h.subs(w, T) + 7) == 0, "h(T) response changed")
    target_A_coefficient = sp.simplify(T**-2 * gamma**2 * h.subs(w, T) ** 2)
    target_C_coefficient = sp.simplify(T**-3 * gamma**3 * h.subs(w, T) ** 3)
    require(target_A_coefficient == 49 * theta**2, "target A coefficient changed")
    require(target_C_coefficient == 343 * theta**3, "target C coefficient changed")
    baseline_degree_floor = (2 + 2 * (6 + 7), 3 + 3 * (6 + 7))
    higher_degree_floor = (2 * 15, 2 * 22 - 13)
    combined_degree_floor = (
        min(baseline_degree_floor[0], higher_degree_floor[0]),
        min(baseline_degree_floor[1], higher_degree_floor[1]),
    )
    require(baseline_degree_floor == (28, 42), "baseline floor changed")
    require(higher_degree_floor == (30, 31), "higher-layer floor changed")
    require(combined_degree_floor == (28, 31), "combined floor changed")

    # Re-entry hostiles: a weight-15 layer with one (w-T) factor when phi!=0,
    # or a weight-16 layer when phi=0, can re-enter at effective pole 14.
    # Hence the theorem must retain the principal-kernel sidecar and must not
    # claim that the raw weight-14 coefficient alone is the target response.
    reentry_effective_weights = (15 - 1, 16 - 2)
    require(reentry_effective_weights == (14, 14), "re-entry hostile changed")

    print("THM4120 JC23 EXTREMAL TARGET AND DE RESPONSE CERTIFICATE")
    print("scope=THM4103_smooth_nonresonant_theta_only;JC2=OPEN")
    print(f"target_c4={c4}")
    print(f"target_discriminant={discriminant}")
    print(f"infinity_valuations={infinity_triple};fibres=I1,I1,IIstar")
    print("NS=rho10;section_fibre_det=-1;E8_rank8;MW={O}")
    print(f"K_rational_punctures={rational_punctures};weight={rational_weight}")
    print("BC=one_quadratic_closed_point;joint_weight=4")
    print(f"degree_before_norm={pre_norm_degrees};after_norm={norm_degrees}")
    print(f"pullback_O_indices={pullback_O};pole_degrees_A_C={pole_degrees}")
    print("quadratic_hostile=(U,V)=(a/2,r),q=a^3/2+r^2")
    print(f"DE_T={T};c1={c1};c2={c2}")
    print(f"phi0_epsilon_plus_kappa={sum_forced};reentry_costs=1,2")
    print("DE_layer_kernel=(w-T);w=x^6*t^7")
    print("DE_top_bracket=7*x^10*t^6*(2*f*g'-3*f'*g)")
    print(f"A14_layers={A14_solutions};C21_layers={C21_solutions}")
    print(f"A15_layers={A15_solutions};C22_layers={C22_solutions}")
    print("baseline_response=f=gamma^2*h^2,g=gamma^3*h^3,h(0)=1,h(T)=-7")
    print(f"baseline_degree_floor={baseline_degree_floor}")
    print(f"higher_degree_floor={higher_degree_floor}")
    print(f"combined_degree_floor={combined_degree_floor}")
    print(f"cross_layer_reentry_hostile_effective_weights={reentry_effective_weights}")
    print("RESULT=PASS")


if __name__ == "__main__":
    main()
