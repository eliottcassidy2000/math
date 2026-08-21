#!/usr/bin/env python3
"""Exact controls for the consecutive-factor/integrability Keller gates.

The all-degree UFD, curl, and coefficient-span arguments are proved in the
companion reflection.  This script checks their load-bearing identities and
the sharp tame/hostile examples over exact rational polynomial rings.  Every
gate raises explicitly, so ``python`` and ``python -O`` have identical logic.
"""

from __future__ import annotations

import sympy as sp
from flint import fmpq_mpoly_ctx


x, y = sp.symbols("x y")


def require(label: str, condition: object, checks: list[str]) -> None:
    """Record one exact check without relying on ``assert``."""

    if condition is not True and condition != sp.S.true:
        raise RuntimeError(f"FAIL: {label}: {condition!r}")
    checks.append(label)


def matrix_is_zero(matrix: sp.MatrixBase) -> bool:
    return all(sp.expand(entry) == 0 for entry in matrix)


def vector_substitute(
    vector: sp.MatrixBase, substitutions: dict[sp.Symbol, sp.Expr]
) -> sp.Matrix:
    return sp.Matrix(
        [sp.expand(entry.subs(substitutions, simultaneous=True)) for entry in vector]
    )


def jacobian_table(P: sp.Expr, Q: sp.Expr) -> tuple[sp.Expr, ...]:
    return tuple(sp.expand(entry) for entry in (sp.diff(P, x), sp.diff(P, y),
                                                sp.diff(Q, x), sp.diff(Q, y)))


def distinct_factor_count(poly: sp.Expr) -> int:
    return len(sp.factor_list(sp.Poly(poly, x, y))[1])


def multiplicity_factor_count(poly: sp.Expr) -> int:
    return sum(exponent for _, exponent in sp.factor_list(sp.Poly(poly, x, y))[1])


def coefficient_vectors(entries: tuple[sp.Expr, ...]) -> dict[tuple[int, int], sp.Matrix]:
    polys = [sp.Poly(entry, x, y) for entry in entries]
    support = set().union(*(poly.monoms() for poly in polys))
    support.discard((0, 0))
    return {
        monomial: sp.Matrix([poly.coeff_monomial(monomial) for poly in polys])
        for monomial in sorted(support)
    }


def check_keller(
    label: str, P: sp.Expr, Q: sp.Expr, kappa: sp.Expr, checks: list[str]
) -> tuple[sp.Expr, ...]:
    a, b, c, d = jacobian_table(P, Q)
    require(f"{label}: determinant", sp.expand(a * d - b * c - kappa) == 0, checks)
    require(f"{label}: first curl", sp.expand(sp.diff(a, y) - sp.diff(b, x)) == 0, checks)
    require(f"{label}: second curl", sp.expand(sp.diff(c, y) - sp.diff(d, x)) == 0, checks)
    return a, b, c, d


def main() -> None:
    checks: list[str] = []

    # Sharp tame control for the multiplicity floor Omega=4.  It also has
    # only three distinct factors, showing that the nonautomorphic omega>=4
    # theorem genuinely uses nonautomorphy rather than the Keller equation.
    ray = x + y
    sharp_Q = y + ray**2 / 2
    sharp_P = ray + sharp_Q
    sharp_entries = check_keller("sharp tame factor floor", sharp_P, sharp_Q, 1, checks)
    sharp_a, sharp_b, sharp_c, sharp_d = sharp_entries
    sharp_T = sp.expand(sharp_b * sharp_c)
    require("sharp T factorization", sp.factor(sharp_T) == ray * (ray + 2), checks)
    require("sharp T+1 factorization", sp.factor(sharp_a * sharp_d) == (ray + 1) ** 2,
            checks)
    require("sharp Omega total is four",
            multiplicity_factor_count(sharp_T)
            + multiplicity_factor_count(sharp_T + 1) == 4, checks)
    require("sharp omega total is three",
            distinct_factor_count(sharp_T) + distinct_factor_count(sharp_T + 1) == 3,
            checks)
    U, V = sp.symbols("U V")
    inverse_ray = U - V
    sharp_inverse = sp.Matrix(
        [inverse_ray - (V - inverse_ray**2 / 2), V - inverse_ray**2 / 2]
    )
    require(
        "sharp inverse: forward after inverse",
        matrix_is_zero(
            vector_substitute(
                sp.Matrix([sharp_P, sharp_Q]),
                {x: sharp_inverse[0], y: sharp_inverse[1]},
            )
            - sp.Matrix([U, V])
        ),
        checks,
    )
    require(
        "sharp inverse: inverse after forward",
        matrix_is_zero(
            vector_substitute(sharp_inverse, {U: sharp_P, V: sharp_Q})
            - sp.Matrix([x, y])
        ),
        checks,
    )

    # Consecutive factorization alone is not integrability.  This SL_2
    # matrix over C[xy] has determinant one and exactly the sharp factor
    # pattern, but both row curls fail by x-y.
    z = x * y
    hostile = sp.Matrix([[z + 1, z + 2], [z, z + 1]])
    require("toric hostile determinant", sp.expand(hostile.det()) == 1, checks)
    require(
        "toric hostile first curl defect",
        sp.expand(sp.diff(hostile[0, 0], y) - sp.diff(hostile[0, 1], x)) == x - y,
        checks,
    )
    require(
        "toric hostile second curl defect",
        sp.expand(sp.diff(hostile[1, 0], y) - sp.diff(hostile[1, 1], x)) == x - y,
        checks,
    )

    # In the zero-level branch of the matching-ray theorem, curl transport
    # fixes the two top ratios.  Their determinant compatibility misses by
    # the strictly positive integer p*K+q*L+1.
    p, q, K, L = sp.symbols("p q K L", positive=True, integer=True)
    toric_gap = sp.expand((p * K + 1) * (q * L + 1) - p * q * K * L)
    require("symbolic toric leading gap", toric_gap == p * K + q * L + 1, checks)
    require(
        "finite toric leading-gap hostile sweep",
        all(
            int(toric_gap.subs({p: pp, q: qq, K: kk, L: ll})) > 0
            for pp in range(1, 8)
            for qq in range(1, 8)
            for kk in range(1, 8)
            for ll in range(1, 8)
        ),
        checks,
    )

    # Multiplicity-safe identities for the only unbalanced omega=4 cell.
    # Here R=S*E and Q0=R+2.  If Q0 belongs to C, the two curls reduce to
    # Q0*E_y=S*E_x.  If it belongs to B, they reduce to Q0*S_x=E*S_y.
    S = sp.Function("S")(x, y)
    E = sp.Function("E")(x, y)
    R = S * E
    Q0 = R + 2
    B_left = S
    C_left = E * Q0
    left_remainder = sp.expand(
        sp.diff(C_left, y)
        - sp.diff(R, x)
        - E * (sp.diff(R, y) - sp.diff(B_left, x))
        - (Q0 * sp.diff(E, y) - S * sp.diff(E, x))
    )
    require("multiplicity-safe left allocation identity", left_remainder == 0, checks)
    B_right = S * Q0
    C_right = E
    right_remainder = sp.expand(
        sp.diff(B_right, x)
        - sp.diff(R, y)
        - S * (sp.diff(R, x) - sp.diff(C_right, y))
        - (Q0 * sp.diff(S, x) - E * sp.diff(S, y))
    )
    require("multiplicity-safe right allocation identity", right_remainder == 0, checks)

    # Universal three-element Gauss/cluster chart for the balanced (2,2)
    # factor cell.  Its determinant is automatic; the two displayed curl
    # equations are the entire remaining global-integrability burden.
    U0 = sp.Function("U")(x, y)
    V0 = sp.Function("V")(x, y)
    W0 = sp.Function("W")(x, y)
    E_plus_V = sp.Matrix([[1, V0], [0, 1]])
    E_minus_U = sp.Matrix([[1, 0], [U0, 1]])
    E_plus_W = sp.Matrix([[1, W0], [0, 1]])
    gauss_matrix = sp.simplify(E_plus_V * E_minus_U * E_plus_W)
    gauss_expected = sp.Matrix(
        [[1 + U0 * V0, V0 + W0 * (1 + U0 * V0)], [U0, 1 + U0 * W0]]
    )
    require("three-element Gauss identity", gauss_matrix == gauss_expected, checks)
    require("three-element Gauss determinant", sp.expand(gauss_matrix.det()) == 1, checks)
    require(
        "three-element Gauss lower curl",
        sp.expand(sp.diff(gauss_matrix[1, 0], y) - sp.diff(gauss_matrix[1, 1], x))
        == sp.expand(sp.diff(U0, y) - sp.diff(U0 * W0, x)),
        checks,
    )
    require(
        "three-element Gauss upper curl",
        sp.expand(sp.diff(gauss_matrix[0, 0], y) - sp.diff(gauss_matrix[0, 1], x))
        == sp.expand(
            sp.diff(1 + U0 * V0, y)
            - sp.diff(V0 + W0 * (1 + U0 * V0), x)
        ),
        checks,
    )

    # Positive boundary-only control from the normalized nodal frame.  It
    # factors in the opposite elementary orientation, and branch holonomy is
    # literally the final right shear.  This does not assert a global lift.
    node_u, node_k = sp.symbols("node_u node_k")
    E_plus = lambda value: sp.Matrix([[1, value], [0, 1]])
    E_minus = lambda value: sp.Matrix([[1, 0], [value, 1]])
    nodal_frame = sp.Matrix(
        [
            [2 * node_u, 1 + 2 * node_k * node_u],
            [3 * node_u**2 - 1,
             3 * node_u / 2 + node_k * (3 * node_u**2 - 1)],
        ]
    )
    nodal_factorization = (
        E_minus(3 * node_u / 2 - 1)
        * E_plus(1)
        * E_minus(2 * node_u - 1)
        * E_plus(node_k)
    )
    require("nodal boundary elementary factorization",
            matrix_is_zero(nodal_factorization - nodal_frame), checks)

    # A two-shear gradient-coprime integrable seed.  Both consecutive fibres
    # split, but a rank-one coefficient observer exposes a constant
    # derivative.  Gradient-coprime is not the unit-gradient-ideal/no-critical
    # condition: this example has critical points for T.
    shear_u = x + y**2 / 2
    shear_v = y + shear_u**2 / 2
    two_Q = shear_v
    two_P = shear_u + shear_v
    two_entries = check_keller("two-shear gradient-coprime seed", two_P, two_Q, 1,
                               checks)
    two_a, two_b, two_c, two_d = two_entries
    two_T = sp.expand(two_b * two_c)
    require(
        "two-shear T factorization",
        sp.expand(two_T - shear_u * (1 + y + shear_u * y)) == 0,
        checks,
    )
    require(
        "two-shear T+1 factorization",
        sp.expand(two_a * two_d - (1 + shear_u) * (1 + shear_u * y)) == 0,
        checks,
    )
    require("two-shear T has two distinct factors",
            distinct_factor_count(two_T) == 2, checks)
    require("two-shear T+1 has two distinct factors",
            distinct_factor_count(two_T + 1) == 2, checks)
    require("two-shear gradient gcd is one in SymPy",
            sp.Poly(sp.gcd(sp.diff(two_T, x), sp.diff(two_T, y)), x, y).total_degree() == 0,
            checks)
    require(
        "two-shear T has critical points",
        all(
            sp.diff(two_T, variable).subs({x: x0, y: y0}) == 0
            for x0, y0 in ((-sp.Rational(3, 2), 1), (-sp.Rational(1, 2), -1))
            for variable in (x, y)
        ),
        checks,
    )
    two_vectors = coefficient_vectors(two_entries)
    two_matrix = sp.Matrix.hstack(*two_vectors.values())
    require("two-shear coefficient span is three", two_matrix.rank() == 3, checks)
    observer_W = sp.Matrix([[-1, 0], [1, 0]])
    observer_flat = sp.Matrix(list(observer_W))
    require("rank-one observer orientation", observer_W.rank() == 1, checks)
    require(
        "rank-one observer kills nonconstant coefficients",
        all((observer_flat.T * vector)[0] == 0 for vector in two_vectors.values()),
        checks,
    )
    require(
        "rank-one observer reads c-a=-1",
        sp.expand(-two_a + two_c) == -1,
        checks,
    )
    two_inverse_u = U - V
    two_inverse_y = V - two_inverse_u**2 / 2
    two_inverse = sp.Matrix([two_inverse_u - two_inverse_y**2 / 2, two_inverse_y])
    require(
        "two-shear inverse: forward after inverse",
        matrix_is_zero(
            vector_substitute(
                sp.Matrix([two_P, two_Q]), {x: two_inverse[0], y: two_inverse[1]}
            )
            - sp.Matrix([U, V])
        ),
        checks,
    )

    # A three-shear tame hostile has gradient-coprime T and reaches full
    # coefficient span four.  Thus neither condition is a counterexample
    # certificate.  Four displayed coefficient columns give a transparent
    # exact rank witness.  Its T still has a critical point at the origin.
    three_Q = shear_v
    three_P = shear_u + shear_v**2 / 2
    three_entries = check_keller("three-shear gradient-coprime rank-four seed", three_P,
                                 three_Q, 1, checks)
    three_a, three_b, three_c, three_d = three_entries
    three_T = sp.expand(three_b * three_c)
    require(
        "three-shear T factorization",
        sp.expand(three_T - shear_u * (y + shear_v * (1 + shear_u * y))) == 0,
        checks,
    )
    require(
        "three-shear T+1 factorization",
        sp.expand(three_a * three_d - (1 + shear_u * shear_v) * (1 + shear_u * y))
        == 0,
        checks,
    )
    require("three-shear T has two distinct factors",
            distinct_factor_count(three_T) == 2, checks)
    require("three-shear T+1 has two distinct factors",
            distinct_factor_count(three_T + 1) == 2, checks)
    sympy_gradient_gcd = sp.Poly(
        sp.gcd(sp.diff(three_T, x), sp.diff(three_T, y)), x, y
    )
    require("three-shear gradient gcd is one in SymPy",
            sympy_gradient_gcd.total_degree() == 0,
            checks)
    require(
        "three-shear T has a critical point",
        sp.diff(three_T, x).subs({x: 0, y: 0}) == 0
        and sp.diff(three_T, y).subs({x: 0, y: 0}) == 0,
        checks,
    )
    three_vectors = coefficient_vectors(three_entries)
    three_matrix = sp.Matrix.hstack(*three_vectors.values())
    require("three-shear coefficient span is four", three_matrix.rank() == 4, checks)
    pivot_monomials = ((0, 1), (0, 2), (0, 3), (0, 6))
    expected_pivots = (
        sp.Matrix([0, 2, 0, 0]),
        sp.Matrix([0, 0, sp.Rational(1, 2), 0]),
        sp.Matrix([sp.Rational(1, 2), 0, 0, sp.Rational(1, 2)]),
        sp.Matrix([sp.Rational(1, 16), 0, 0, 0]),
    )
    require(
        "three-shear displayed pivot columns",
        all(three_vectors[monomial] == expected
            for monomial, expected in zip(pivot_monomials, expected_pivots)),
        checks,
    )
    require(
        "three-shear displayed pivot determinant",
        sp.Matrix.hstack(*(three_vectors[m] for m in pivot_monomials)).det()
        == -sp.Rational(1, 32),
        checks,
    )
    three_inverse_u = U - V**2 / 2
    three_inverse_y = V - three_inverse_u**2 / 2
    three_inverse = sp.Matrix([three_inverse_u - three_inverse_y**2 / 2,
                               three_inverse_y])
    require(
        "three-shear inverse: forward after inverse",
        matrix_is_zero(
            vector_substitute(
                sp.Matrix([three_P, three_Q]),
                {x: three_inverse[0], y: three_inverse[1]},
            )
            - sp.Matrix([U, V])
        ),
        checks,
    )

    # Independent python-FLINT computation of the gradient gcd.
    ctx = fmpq_mpoly_ctx.get(["x", "y"])
    fx, fy = ctx.gens()
    flint_u = fx + fy**2 / 2
    flint_v = fy + flint_u**2 / 2
    flint_T = flint_u * (fy + flint_v * (1 + flint_u * fy))
    flint_gcd = flint_T.derivative(0).gcd(flint_T.derivative(1))
    require("three-shear gradient gcd is one in FLINT", flint_gcd.is_constant(), checks)

    # The rank-three trace/Hamiltonian router and its exact matching carrier.
    lam = sp.symbols("lambda")
    H = x**4 + 2 * x**2 * y**2 - 3 * x * y**3 + 5 * y**4 + x**2 * y
    Hxx = sp.diff(H, x, 2)
    Hxy = sp.diff(H, x, y)
    Hyy = sp.diff(H, y, 2)
    routed_F = sp.Matrix([lam * x + sp.diff(H, y), lam * y - sp.diff(H, x)])
    routed_DF = routed_F.jacobian([x, y])
    hessian_mu = sp.expand(Hxx * Hyy - Hxy**2)
    require("rank-three router has constant trace", sp.trace(routed_DF) == 2 * lam, checks)
    require(
        "rank-three determinant router",
        sp.expand(routed_DF.det() - (lam**2 + hessian_mu)) == 0,
        checks,
    )
    routed_T = sp.expand(routed_DF[0, 1] * routed_DF[1, 0])
    require(
        "rank-three matching carrier T",
        sp.expand(routed_T + Hxy**2 + hessian_mu) == 0,
        checks,
    )
    require(
        "rank-three matching carrier T+kappa",
        sp.expand(
            routed_T + (lam**2 + hessian_mu) - (lam**2 - Hxy**2)
        )
        == 0,
        checks,
    )

    print("consecutive-factor/integrability exact scout")
    print("universe=C[x,y], kappa=1 controls; exact Q and Q(lambda) arithmetic")
    print("factor_floor_control=Omega_total:4,omega_total:3,tame")
    print("toric_gap=(p*K+1)(q*L+1)-p*q*K*L=p*K+q*L+1")
    print("toric_gap_sweep=1<=p,q,K,L<=7")
    print("balanced_generator=three-element_Gauss_chart;two_curl_PDEs_remain")
    print("nodal_boundary=three_nontrivial_shears+holonomy_right_shear")
    print("two_shear=balanced(2,2),gradient-coprime,critical,span:3,observer:c-a=-1")
    print("three_shear=balanced(2,2),gradient-coprime,critical,span:4,pivot_det:-1/32")
    print("gradient_gcd_crosscheck=SymPy+python-FLINT")
    print("rank3_matching=T:-(R^2+mu),T+kappa:lambda^2-R^2")
    print(f"checks={len(checks)}")
    print("PASS")


if __name__ == "__main__":
    main()
