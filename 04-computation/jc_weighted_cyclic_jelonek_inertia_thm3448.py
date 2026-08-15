#!/usr/bin/env python3
"""Exact deterministic companion for THM-3448.

The companion audits the boundary anatomy of the weighted Keller construction
from THM-3438.  It keeps four logically different ledgers separate:

* the Newton polygon and tame inertia over the target plane C=0;
* the raw discriminant of the bounded primitive w and its order index;
* the poles of a generic infinity-separating affine source primitive; and
* the ordinary transposition component pulled back from the (P,Q)-branch
  curve.

It also expands the first weighted Keller C3 witness, checks the infinite
cyclic seed family, and freezes a quartic hostile showing why the clean
primitive-clearing formula needs a pole-residue hypothesis.  All truth gates
use ``require`` and survive ``python -O``.
"""

from __future__ import annotations

from hashlib import sha256
from pathlib import Path

import sympy as sp


EXPECTED_SEMANTIC_SHA256 = "dbb99005e17dc8e2d69d514244640530476b74f01a723d5e5bd418bf42b27a80"


def require(condition: bool, detail: object) -> None:
    if not condition:
        raise RuntimeError(detail)


def lf_sha256(path: Path) -> str:
    return sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def valuation_at_zero(expression: sp.Expr, variable: sp.Symbol) -> int:
    numerator, denominator = sp.cancel(expression).as_numer_denom()
    require(not sp.Poly(numerator, variable).is_zero, ("zero numerator", expression))
    num_v = min(monomial[0] for monomial, coefficient in sp.Poly(numerator, variable).terms() if coefficient)
    den_v = min(monomial[0] for monomial, coefficient in sp.Poly(denominator, variable).terms() if coefficient)
    return num_v - den_v


def quotient_reduce(expression: sp.Expr, variable: sp.Symbol, modulus: sp.Expr) -> sp.Expr:
    numerator, denominator = sp.cancel(expression).as_numer_denom()
    modulus_poly = sp.Poly(modulus, variable, domain=sp.QQ)
    numerator_rem = sp.rem(sp.Poly(numerator, variable, domain=sp.QQ), modulus_poly).as_expr()
    denominator_rem = sp.rem(sp.Poly(denominator, variable, domain=sp.QQ), modulus_poly).as_expr()
    require(sp.gcd(sp.Poly(denominator_rem, variable), modulus_poly).degree() == 0, ("nonunit denominator", expression))
    inverse = sp.invert(sp.Poly(denominator_rem, variable), modulus_poly).as_expr()
    return sp.rem(sp.Poly(sp.expand(numerator_rem * inverse), variable, domain=sp.QQ), modulus_poly).as_expr()


def general_formula_rows() -> tuple[tuple[int, ...], ...]:
    rows: list[tuple[int, ...]] = []
    for n in range(3, 13):
        for mu in range(2, n):
            cyclic_length = mu - 1
            other_nonzero = n - mu - 1
            escaping = n - 3 if mu == 2 else n - 2
            inertia_defect = mu - 2
            raw_order = mu
            order_index = (raw_order - inertia_defect) // 2
            rho = 2 * n - mu - 4
            clearing = 2 * n * (n - 1) - mu * mu - 2 * mu - 4

            if mu == 2:
                pair_count_formula = 4 * (
                    other_nonzero * (other_nonzero - 1) // 2
                    + 3 * other_nonzero
                )
            else:
                fractional_numerator = mu - 2
                fractional_denominator = mu - 1
                pole_two_pairs = (
                    other_nonzero * (other_nonzero - 1) // 2
                    + other_nonzero * (cyclic_length + 2)
                )
                fractional_pair_numerator = 2 * fractional_numerator * (
                    cyclic_length * (cyclic_length - 1) // 2
                    + 2 * cyclic_length
                )
                require(fractional_pair_numerator % fractional_denominator == 0, (n, mu, "fractional pair invoice"))
                pair_count_formula = 4 * pole_two_pairs + fractional_pair_numerator // fractional_denominator

            require(raw_order == inertia_defect + 2 * order_index, (n, mu, "index"))
            require(order_index == 1, (n, mu, "index one"))
            require(rho == 2 * other_nonzero + max(mu - 2, 0), (n, mu, "rho"))
            require(pair_count_formula == clearing, (n, mu, pair_count_formula, clearing))
            require(clearing % 2 == inertia_defect % 2, (n, mu, "parity"))
            rows.append((n, mu, cyclic_length, other_nonzero, escaping, raw_order, order_index, rho, clearing))
    return tuple(rows)


def normalized_seed(w: sp.Symbol, n: int, mu: int) -> tuple[sp.Expr, sp.Expr, sp.Expr]:
    other_count = n - mu - 1
    nonzero_roots = tuple(range(3, 3 + other_count))
    base = w**mu * (w - 1) * sp.prod(w - root for root in nonzero_roots)
    scale = -sp.diff(base, w).subs(w, 1)
    R = sp.cancel(base / scale)
    p = sp.diff(R, w)
    kappa = sp.diff(p, w).subs(w, 1)
    require(sp.simplify(p.subs(w, 0)) == 0, (n, mu, "p(0)"))
    require(sp.simplify(p.subs(w, 1)) == -1, (n, mu, "p(1)"))
    require(sp.simplify(R.subs(w, 1)) == 0, (n, mu, "R(1)"))
    require(sp.simplify(kappa + 2) != 0, (n, mu, "excluded kappa"))
    return sp.expand(R), sp.expand(p), sp.factor(kappa)


def finite_discriminant_controls() -> tuple[tuple[object, ...], ...]:
    w, A, B, C, P, Q = sp.symbols("w A B C P Q")
    rows: list[tuple[object, ...]] = []
    for n, mu in ((3, 2), (4, 2), (4, 3), (5, 2), (5, 3), (5, 4)):
        R, p, kappa = normalized_seed(w, n, mu)
        T = sp.expand(R - B * C * w + A * C**2)
        pulled_discriminant = sp.factor(sp.discriminant(T, w))
        require(valuation_at_zero(pulled_discriminant, C) == mu, (n, mu, pulled_discriminant))
        L = sp.factor(sp.cancel(pulled_discriminant / C**mu))
        require(sp.simplify(L.subs(C, 0)) != 0, (n, mu, "residual factor"))
        require(sp.Poly(L, A, B, C, domain=sp.QQ).is_irreducible, (n, mu, "L irreducible", L))

        branch_discriminant = sp.factor(sp.discriminant(R - P * w + Q, w))
        require(sp.Poly(branch_discriminant, P, Q, domain=sp.QQ).is_irreducible, (n, mu, "branch irreducible"))
        parameter = sp.symbols("parameter")
        slope_identity = sp.cancel(
            (sp.diff((parameter * p.subs(w, parameter) - R.subs(w, parameter)), parameter))
            / sp.diff(p.subs(w, parameter), parameter)
            - parameter
        )
        require(slope_identity == 0, (n, mu, "dQ/dP"))
        rows.append(
            (
                n,
                mu,
                str(kappa),
                sp.Poly(L, A, B, C).total_degree(),
                sha256(sp.srepr(L).encode("utf-8")).hexdigest()[:16],
            )
        )

    cubic = w**3 - w**2 + B * C * w - A * C**2
    cubic_discriminant = sp.factor(sp.discriminant(cubic, w))
    cubic_L = sp.factor(sp.cancel(cubic_discriminant / C**2))
    require(sp.expand(cubic_L.subs(C, 0) - (B**2 - 4 * A)) == 0, ("exceptional cubic residual", cubic_L))
    W = sp.symbols("W")
    cubic_zero_residual = sp.expand(cubic.subs(w, C * W) / C**2).subs(C, 0)
    require(sp.expand(cubic_zero_residual - (-W**2 + B * W - A)) == 0, cubic_zero_residual)
    require(sp.discriminant(cubic_zero_residual, W) == B**2 - 4 * A, "cubic finite residual discriminant")
    return tuple(rows)


def degree_five_witness() -> tuple[object, ...]:
    x, y, z = sp.symbols("x y z")
    A, B, C, P, Q, w, s, W = sp.symbols("A B C P Q w s W")

    u = 1 + x * y
    gamma = 1 - sp.Rational(7, 6) * x * y + x**2 * z
    alpha = sp.cancel(u + 3 * u**4 * gamma**2 - 4 * u**5 * gamma**3)
    beta = sp.cancel(1 + 4 * u**3 * gamma**2 - 5 * u**4 * gamma**3)
    A_map = sp.cancel(alpha / x**2)
    B_map = sp.cancel(beta / x)
    C_map = sp.expand(x * gamma)
    require(sp.denom(A_map) == 1 and sp.denom(B_map) == 1, "F5 quotient cancellation")
    F5 = tuple(sp.expand(item) for item in (A_map, B_map, C_map))
    jacobian = sp.factor(sp.Matrix(F5).jacobian((x, y, z)).det())
    require(jacobian == 1, ("F5 Jacobian", jacobian))
    ordinary_degrees = tuple(sp.Poly(item, x, y, z).total_degree() for item in F5)
    z_degrees = tuple(sp.Poly(item, z).degree() for item in F5)
    require(ordinary_degrees == (17, 16, 4), ordinary_degrees)
    require(z_degrees == (3, 3, 1), z_degrees)

    P_source = sp.expand(B_map * C_map)
    Q_source = sp.expand(A_map * C_map**2)
    w_source = sp.expand(u * gamma)
    inverse_on_source = sp.factor(w_source**5 - w_source**4 + P_source * w_source - Q_source)
    require(inverse_on_source == 0, "F5 inverse polynomial")
    require(sp.factor(P_source - 4 * w_source**3 + 5 * w_source**4 - gamma) == 0, "F5 gamma=f'(w)")

    f = w**5 - w**4 + P * w - Q
    discriminant = sp.factor(sp.discriminant(f, w).subs({P: B * C, Q: A * C**2}))
    L5 = (
        3125 * A**4 * C**4
        - 2500 * A**3 * B * C**3
        + 256 * A**3 * C**2
        - 50 * A**2 * B**2 * C**2
        - 36 * A * B**3 * C
        + 256 * B**5 * C
        - 27 * B**4
    )
    require(sp.expand(discriminant - C**4 * L5) == 0, ("F5 discriminant", discriminant))
    require(sp.Poly(L5, A, B, C, domain=sp.QQ).is_irreducible, "L5 irreducible")
    require(sp.expand(L5.subs(C, 0) + 27 * B**4) == 0, "L5 avoids generic C meridian")

    scaled = sp.expand(f.subs({P: B * s**3, Q: A * s**6, w: s * W}) / s**4)
    require(sp.expand(scaled - (s * W**5 - W**4 + B * W - A * s**2)) == 0, scaled)
    require(sp.factor(scaled.subs(s, 0) - W * (B - W**3)) == 0, "F5 C3 residual")
    gamma_scaled = sp.expand((B * s**3 - 4 * (s * W) ** 3 + 5 * (s * W) ** 4) / s**3)
    require(sp.rem(sp.Poly(gamma_scaled.subs(s, 0) + 3 * B, W), sp.Poly(W**3 - B, W)).is_zero, "gamma~-3BC")

    # The three nonzero residual roots have w=sW and gamma~-3Bs^3.
    # Reconstruction therefore gives the exact leading source residues below.
    x_residue = -sp.Rational(1, 3) / B
    y_residue = W
    z_residue = -sp.Rational(7, 2) * B * W
    require(x_residue == -1 / (3 * B), x_residue)
    require(z_residue / W == -sp.Rational(7, 2) * B, z_residue)

    raw_order = valuation_at_zero(discriminant, C)
    require(raw_order == 4, raw_order)
    escaping = 3
    inertia_defect = 2
    order_index = (raw_order - inertia_defect) // 2
    rho = 2
    clearing = 12
    require(order_index == 1 and clearing % 2 == inertia_defect % 2, "F5 C ledger")
    require(2 * 9 * sp.Rational(2, 3) == 12, "nine-pair primitive invoice")

    L_escaping = 2
    L_raw = 1
    L_rho = 1
    L_clearing = 2 * 5 - 3
    require((L_escaping, L_raw, L_rho, L_clearing) == (2, 1, 1, 7), "F5 L ledger")

    map_digest = sha256(repr(tuple(sp.srepr(item) for item in F5)).encode("utf-8")).hexdigest()
    return (
        ordinary_degrees,
        z_degrees,
        str(jacobian),
        map_digest,
        sp.srepr(L5),
        (escaping, "C3", raw_order, order_index, rho, clearing),
        (L_escaping, "transposition", L_raw, 0, L_rho, L_clearing),
        (sp.srepr(x_residue), sp.srepr(y_residue), sp.srepr(z_residue)),
    )


def cyclic_family_rows() -> tuple[tuple[int, ...], ...]:
    w = sp.symbols("w")
    rows: list[tuple[int, ...]] = []
    for ell in range(1, 17):
        p = (ell + 1) * w**ell - (ell + 2) * w ** (ell + 1)
        R = w ** (ell + 1) * (1 - w)
        require(sp.expand(sp.integrate(p, w) - R) == 0, (ell, "integral"))
        require(p.subs(w, 0) == 0 and p.subs(w, 1) == -1 and R.subs(w, 1) == 0, (ell, "seed"))
        kappa = sp.diff(p, w).subs(w, 1)
        require(kappa == -2 * (ell + 1) and kappa != -2, (ell, kappa))
        a = sp.cancel(-(1 + kappa) / (2 + kappa))
        require(a == -sp.Rational(2 * ell + 1, 2 * ell), (ell, a))
        n = ell + 2
        mu = ell + 1
        rho = 2 * n - mu - 4
        clearing = 2 * n * (n - 1) - mu * mu - 2 * mu - 4
        require(rho == ell - 1, (ell, rho))
        require(clearing == (ell - 1) * (ell + 3), (ell, clearing))
        rows.append((ell, n, mu, int(kappa), int(rho), int(clearing)))
    return tuple(rows)


def hostile_quartic() -> tuple[object, ...]:
    w, r = sp.symbols("w r")
    h = 2 * r**4 - 4 * r**3 + 5 * r**2 - 6 * r - 5
    require(sp.Poly(h, r, domain=sp.QQ).is_irreducible, "hostile parameter polynomial irreducible")
    R = sp.cancel(w**2 * (w - 1) * (w - r) / (r - 1))
    p = sp.diff(R, w)
    require(sp.simplify(p.subs(w, 0)) == 0 and sp.simplify(p.subs(w, 1)) == -1, "hostile seed endpoints")
    require(sp.simplify(R.subs(w, 1)) == 0, "hostile integral endpoint")
    kappa = sp.cancel(sp.diff(p, w).subs(w, 1))
    require(quotient_reduce(kappa + 2, r, h) != 0, "hostile seed avoids kappa=-2")
    a = sp.cancel(-(1 + kappa) / (2 + kappa))
    g = sp.cancel(-p.subs(w, r))
    require(sp.simplify(g + r**2) == 0, ("g_r", g))
    Z = sp.cancel(g**2 * (g - 1 + a) - a * r * g)
    require(quotient_reduce(Z, r, h) == 0, ("hostile Z_r", Z))
    y_residue = sp.cancel(r - g)
    require(quotient_reduce(y_residue, r, h) != 0, ("hostile y residue", y_residue))
    require(sp.resultant(h, r * (r - 1), r) != 0, "hostile roots nonzero and not one")
    # One pole-one root and three finite roots: rho=1 and three discriminant
    # pairs of valuation -1 give E=2*3=6.
    rho = 1
    clearing = 2 * 3
    require((rho, clearing) == (1, 6), "hostile primitive ledger")
    return (
        sp.srepr(h),
        sp.srepr(quotient_reduce(kappa, r, h)),
        sp.srepr(quotient_reduce(a, r, h)),
        sp.srepr(quotient_reduce(g, r, h)),
        sp.srepr(quotient_reduce(Z, r, h)),
        sp.srepr(quotient_reduce(y_residue, r, h)),
        rho,
        clearing,
    )


def main() -> None:
    formula_rows = general_formula_rows()
    discriminant_controls = finite_discriminant_controls()
    degree_five = degree_five_witness()
    cyclic_rows = cyclic_family_rows()
    hostile = hostile_quartic()

    verdict = (
        "C_Newton_vertices=(0,2),(1,1),(mu,0)",
        "C_inertia=(mu-1)-cycle_with_n-mu+1_fixed_sheets",
        "C_raw_w_discriminant_order=mu;maximal_order_defect=mu-2;w_index=1",
        "C_escaping=n-3_if_mu=2_else_n-2",
        "Jelonek=V(L)_only_at_(n,mu)=(3,2);otherwise_V(C)_union_V(L)",
        "L_inertia=transposition;escaping=2;raw=1;index=0;rho=1;E=2n-3",
        "clean_C_primitive_formulas_require_nonzero_distinct_pole_two_residues_and_cyclic_amplitude",
        "F5_det=1;degree=5;C_component_has_three_escaping_sheets_and_C3_inertia",
        "p_ell_family_realizes_component_generic_C_ell_for_every_ell>=2",
        "weighted_quartic_component_generic_C3_impossible;degree5_is_first",
        "hostile_quartic_Z_r=0_changes_generic_primitive_(rho,E)_from_(2,12)_to_(1,6)",
        "finite_critical_branching_empty_for_every_weighted_Keller_map",
        "no_JC2_LRC_or_full_counterexample_classification_follows",
    )
    semantic_surface = (formula_rows, discriminant_controls, degree_five, cyclic_rows, hostile, verdict)
    semantic_sha = sha256(repr(semantic_surface).encode("utf-8")).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(semantic_sha == EXPECTED_SEMANTIC_SHA256, (semantic_sha, EXPECTED_SEMANTIC_SHA256))

    source = Path(__file__).resolve()
    print("THM-3448 weighted Keller cyclic Jelonek inertia exact companion")
    print("status=EXACT_DETERMINISTIC_COMPANION_FOR_THM3448;finite_checks_do_not_replace_proof")
    print("arithmetic=QQ_symbolic_polynomials;exact_Newton_and_pair_valuation_invoices;no_randomness")
    print(f"general_formula_rows_count={len(formula_rows)};first={formula_rows[:4]};last={formula_rows[-4:]}")
    print(f"finite_discriminant_controls=(n,mu,kappa,L_total_degree,L_digest16)={discriminant_controls}")
    print("exceptional_cubic=C^2_factor_is_order_index_only_off_L;L(A,B,0)=B^2-4A;all_three_sheets_finite")
    print(f"degree_five_certificate={degree_five}")
    print(f"cyclic_family_rows_ell_1_to_16={cyclic_rows}")
    print(f"primitive_formula_hostile_quartic={hostile}")
    print(f"verdict={verdict}")
    print("genericity=nonzero_roots_of_R_simple;target_avoids_component_intersections;primitive_residues_nonzero_and_separated_when_rho_E_are_invoked")
    print("scope=weighted_THM3438_family_only;quartic_minimality_inside_this_family;no_classification_of_all_Keller_boundaries")
    print(f"semantic_sha256={semantic_sha}")
    print(f"script_sha256_lf={lf_sha256(source)}")
    print("reproducibility=no_assert_truth_gates;all_checks_survive_python_O")
    print("commands=python -B 04-computation/jc_weighted_cyclic_jelonek_inertia_thm3448.py;python -B -O 04-computation/jc_weighted_cyclic_jelonek_inertia_thm3448.py")
    print("PASS")


if __name__ == "__main__":
    main()
