#!/usr/bin/env python3
"""Independent symbolic audit of the Russell debt/contact correspondence.

This companion imports no production computation.  SymPy expands the formal
collision equations directly, while the primary high-order scout uses a
hand-written truncated-series algebra over a four-dimensional Fraction ring.

Universe: Q=Q1+x^2(x^2-1)^2(a+b*x+c*x^2), the constant pencil Q+s, three
labelled sections near -1,0,1, a common (C,E) target, and one compensating x
coefficient.  We successively restrict to the old zero-D2 plane, zero-D4
parabola, and irreducible exceptional quartic.  No global collision locus,
polynomial termination, target lift, Keller pair, or JC(2) claim is tested.
"""

from __future__ import annotations

import hashlib
import sys

import sympy as sp


if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(newline="\n")


CHECKS = 0


def require(condition: bool, label: str) -> None:
    global CHECKS
    if not condition:
        raise RuntimeError(label)
    CHECKS += 1


x, s = sp.symbols("x s")
a, b, c = sp.symbols("a b c")
p, r = sp.symbols("p r")
alpha = sp.symbols("alpha")

P = x**2 * (x**2 - 1) ** 2
Q1 = (
    x**5
    + sp.Rational(9, 2) * x**4
    - 2 * x**3
    - sp.Rational(27, 4) * x**2
    + x
    - sp.Rational(3, 4)
)
Q_ABC = sp.expand(Q1 + P * (a + b * x + c * x**2))

F4 = (
    72783360 * alpha**4
    - 77822208 * alpha**3
    - 28419741 * alpha**2
    + 7849770 * alpha
    - 1276420
)
F4_POLY = sp.Poly(F4, alpha, domain=sp.QQ)


def identity_reduce(value: sp.Expr) -> sp.Expr:
    return sp.cancel(value)


def alpha_reduce(value: sp.Expr) -> sp.Expr:
    value = sp.cancel(value)
    numerator, denominator = sp.fraction(value)
    require(not denominator.has(alpha), "alpha reduction has scalar denominator")
    # During state specialization the coefficients can still be polynomials
    # in the formal parameter s; regard those as coefficient-ring elements.
    coefficient_domain = sp.QQ.poly_ring(s)
    alpha_polynomial = sp.Poly(sp.expand(numerator), alpha, domain=coefficient_domain)
    alpha_modulus = sp.Poly(F4, alpha, domain=coefficient_domain)
    remainder = alpha_polynomial.rem(alpha_modulus)
    return sp.cancel(remainder.as_expr() / denominator)


J_EXPECTED = sp.Matrix(
    (
        (3, 0, 0, -1, 0, -2),
        (-9, 0, 0, 0, -1, 2),
        (0, 3, 0, -1, 0, 0),
        (0, 4, 0, 0, -1, 0),
        (0, 0, 3, -1, 0, -2),
        (0, 0, 9, 0, -1, -2),
    )
)


def compiler_coordinates(q: sp.Expr) -> tuple[sp.Expr, sp.Expr]:
    D = 1 + x**2 * q
    C = x * D * (D + 2)
    E = q * (D + 3)
    return sp.expand(C), sp.expand(E)


def derive_implicit_jacobian(Q: sp.Expr) -> sp.Matrix:
    C, E = compiler_coordinates(Q)
    tangent = [sp.Matrix((sp.diff(C, x).subs(x, z), sp.diff(E, x).subs(x, z))) for z in (-1, 0, 1)]
    u = sp.symbols("u")
    Cu, Eu = compiler_coordinates(Q + u * x)
    normal = [
        sp.Matrix((sp.diff(Cu, u).subs({x: z, u: 0}), sp.diff(Eu, u).subs({x: z, u: 0})))
        for z in (-1, 0, 1)
    ]
    matrix = sp.zeros(6, 6)
    for index in range(3):
        matrix[2 * index : 2 * index + 2, index] = tangent[index]
        matrix[2 * index : 2 * index + 2, 3:5] = -sp.eye(2)
        matrix[2 * index : 2 * index + 2, 5] = normal[index]
    return matrix.applyfunc(sp.expand)


State = tuple[list[sp.Expr], sp.Expr, sp.Expr, sp.Expr]


def collision_residuals(
    Q: sp.Expr,
    state: State,
    order: int,
    reducer,
) -> sp.Matrix:
    endpoints, common_c, common_e, compensator = state
    rows: list[sp.Expr] = []
    for endpoint in endpoints:
        q = sp.expand(Q.subs(x, endpoint) + s + compensator * endpoint)
        D = 1 + endpoint**2 * q
        C = endpoint * D * (D + 2) - common_c
        E = q * (D + 3) - common_e
        rows.extend(
            (
                reducer(sp.Poly(C, s).coeff_monomial(s**order)),
                reducer(sp.Poly(E, s).coeff_monomial(s**order)),
            )
        )
    return sp.Matrix(rows)


def advance(Q: sp.Expr, state: State, order: int, reducer) -> tuple[State, sp.Expr]:
    endpoints, common_c, common_e, compensator = state
    residual = collision_residuals(Q, state, order, reducer)
    solution = [reducer(value) for value in (-J_EXPECTED.inv() * residual)]
    repaired = J_EXPECTED * sp.Matrix(solution) + residual
    require(
        all(reducer(value) == 0 for value in repaired),
        f"order {order} exact repair",
    )
    next_state: State = (
        [sp.expand(value + solution[index] * s**order) for index, value in enumerate(endpoints)],
        sp.expand(common_c + solution[3] * s**order),
        sp.expand(common_e + solution[4] * s**order),
        sp.expand(compensator + solution[5] * s**order),
    )
    return next_state, solution[5]


def specialize_state(state: State, substitutions: dict, reducer) -> State:
    endpoints, common_c, common_e, compensator = state
    return (
        [reducer(value.subs(substitutions)) for value in endpoints],
        reducer(common_c.subs(substitutions)),
        reducer(common_e.subs(substitutions)),
        reducer(compensator.subs(substitutions)),
    )


def first_constant_period_coefficient(
    Q: sp.Expr,
    state: State,
    order: int,
    reducer,
) -> sp.Expr:
    """First live coefficient of the normalized tangent-relation sum.

    If t_-,t_0,t_+ are the moving (C,E) tangents, their canonical raw
    relation is (det(t_0,t_+),det(t_+,t_-),det(t_-,t_0)).  At the base point
    it is (15,-54,39), so division by 54 normalizes the middle weight to -1.
    When all lower sums vanish, the requested coefficient of the normalized
    constant-density period is just the raw-sum coefficient divided by 54.
    """
    endpoints, _, _, compensator = state

    def truncate(value: sp.Expr) -> sp.Expr:
        polynomial = sp.Poly(value, s)
        return sp.expand(
            sum(polynomial.coeff_monomial(s**degree) * s**degree for degree in range(order + 1))
        )

    def evaluate_truncated(polynomial: sp.Expr, value: sp.Expr) -> sp.Expr:
        coefficients = sp.Poly(polynomial, x).all_coeffs()
        answer = sp.Integer(0)
        for coefficient in coefficients:
            answer = truncate(answer * value + coefficient)
        return answer

    Q_derivative = sp.diff(Q, x)
    tangents: list[sp.Matrix] = []
    for endpoint in endpoints:
        endpoint = truncate(endpoint)
        q = truncate(evaluate_truncated(Q, endpoint) + s + compensator * endpoint)
        qx = truncate(evaluate_truncated(Q_derivative, endpoint) + compensator)
        endpoint_squared = truncate(endpoint**2)
        D = truncate(1 + endpoint_squared * q)
        Dx = truncate(2 * endpoint * q + endpoint_squared * qx)
        Cx = truncate(D * (D + 2) + 2 * endpoint * (D + 1) * Dx)
        Ex = truncate(qx * (D + 3) + q * Dx)
        tangents.append(sp.Matrix((Cx, Ex)))
    def tangent_determinant(left: sp.Matrix, right: sp.Matrix) -> sp.Expr:
        return truncate(left[0] * right[1] - left[1] * right[0])

    raw = (
        tangent_determinant(tangents[1], tangents[2]),
        tangent_determinant(tangents[2], tangents[0]),
        tangent_determinant(tangents[0], tangents[1]),
    )
    require(tuple(reducer(value.subs(s, 0)) for value in raw) == (15, -54, 39), "base tangent relation")
    raw_sum = sum(raw)
    for lower in range(order):
        require(
            reducer(sp.Poly(raw_sum, s).coeff_monomial(s**lower)) == 0,
            f"constant period vanishes below order {order}",
        )
    return reducer(sp.Poly(raw_sum, s).coeff_monomial(s**order) / 54)


def main() -> None:
    J = derive_implicit_jacobian(Q_ABC)
    require(J == J_EXPECTED, "compiler-derived implicit Jacobian")
    require(J.det() == -288, "implicit Jacobian determinant")
    require(J[:, :5].rank() == 5, "uncompensated collision rank")
    require(F4_POLY.is_irreducible, "exceptional quartic irreducible")

    base_state: State = ([sp.Integer(-1), sp.Integer(0), sp.Integer(1)], sp.Integer(0), sp.Integer(-3), sp.Integer(0))
    state, contact_1 = advance(Q_ABC, base_state, 1, identity_reduce)
    curvature_minus = sp.diff(Q_ABC, x, 2).subs(x, -1)
    curvature_plus = sp.diff(Q_ABC, x, 2).subs(x, 1)
    debt_2 = sp.factor(-sp.Rational(2, 81) * (5 * curvature_minus + 13 * curvature_plus + 243))
    require(contact_1 == 0, "universal constant first-order tangency")
    period_1 = first_constant_period_coefficient(Q_ABC, state, 1, identity_reduce)
    require(sp.factor(period_1 + debt_2) == 0, "first constant period=-D2")
    state, contact_2 = advance(Q_ABC, state, 2, identity_reduce)
    require(sp.factor(contact_2 - sp.Rational(9, 8) * debt_2) == 0, "contact2=(9/8)D2")
    require(sp.factor(2 * contact_2 + sp.Rational(9, 4) * period_1) == 0, "first adjoint contact identity")

    zero_d2 = {
        a: -sp.Rational(259, 36) + p + 4 * r,
        b: -9 * r,
        c: -p,
    }
    Q_PR = sp.expand(Q_ABC.subs(zero_d2))
    state = specialize_state(state, zero_d2, identity_reduce)
    require(sp.cancel(state[3].coeff(s, 2)) == 0, "zero-D2 plane kills contact2")
    parabola_equation = 729 * p - 42120 * r**2 + 15192 * r + 5717
    debt_4 = sp.Rational(64, 6561) * parabola_equation
    period_2 = first_constant_period_coefficient(Q_PR, state, 2, identity_reduce)
    require(sp.factor(period_2 + debt_4) == 0, "second constant period=-D4")
    state, contact_3 = advance(Q_PR, state, 3, identity_reduce)
    require(sp.factor(contact_3 - sp.Rational(3, 4) * debt_4) == 0, "contact3=(3/4)D4")
    require(sp.factor(3 * contact_3 + sp.Rational(9, 4) * period_2) == 0, "second adjoint contact identity")

    p_of_r = sp.Rational(520, 9) * r**2 - sp.Rational(1688, 81) * r - sp.Rational(5717, 729)
    Q_R = sp.expand(Q_PR.subs(p, p_of_r))
    state = specialize_state(state, {p: p_of_r}, identity_reduce)
    require(sp.cancel(state[3].coeff(s, 3)) == 0, "zero-D4 parabola kills contact3")
    quartic_r = F4.subs(alpha, r)
    debt_6 = -sp.Rational(256, 1594323) * quartic_r
    delta_6 = -debt_6
    period_3 = first_constant_period_coefficient(Q_R, state, 3, identity_reduce)
    require(sp.factor(period_3 + delta_6) == 0, "third constant period=-delta6")
    state, contact_4 = advance(Q_R, state, 4, identity_reduce)
    require(sp.factor(contact_4 - sp.Rational(9, 16) * delta_6) == 0, "contact4=(9/16)delta6")
    require(sp.factor(4 * contact_4 + sp.Rational(9, 4) * period_3) == 0, "third adjoint contact identity")

    Q_ALPHA = sp.expand(Q_R.subs(r, alpha))
    state = specialize_state(state, {r: alpha}, alpha_reduce)
    require(alpha_reduce(state[3].coeff(s, 4)) == 0, "zero-D6 quartic kills contact4")
    period_4 = first_constant_period_coefficient(Q_ALPHA, state, 4, alpha_reduce)
    state, contact_5 = advance(Q_ALPHA, state, 5, alpha_reduce)
    expected_contact_5 = (
        sp.Rational(259188338368, 129140163)
        - sp.Rational(46584993664, 4782969) * alpha
        + sp.Rational(23019960448, 531441) * alpha**2
        + sp.Rational(9180348416, 177147) * alpha**3
    )
    require(alpha_reduce(contact_5 - expected_contact_5) == 0, "independent exceptional contact5")
    kappa_8 = (
        -sp.Rational(5183766767360, 3**19)
        + sp.Rational(931699873280, 3**16) * alpha
        - sp.Rational(460399208960, 3**14) * alpha**2
        - sp.Rational(183606968320, 3**13) * alpha**3
    )
    delta_8 = -kappa_8
    require(alpha_reduce(contact_5 - sp.Rational(9, 20) * delta_8) == 0, "contact5=(9/20)delta8")
    require(contact_5 != 0, "exceptional contact5 nonzero")
    require(alpha_reduce(period_4 + delta_8) == 0, "fourth constant period=-delta8")
    require(alpha_reduce(5 * contact_5 + sp.Rational(9, 4) * period_4) == 0, "fourth adjoint contact identity")

    contact_payload = "|".join(
        sp.sstr(item)
        for item in (contact_1, contact_2, contact_3, contact_4, contact_5, period_1, period_2, period_3, period_4)
    )
    print("constant_normal_debt_contact_independent_audit")
    print("engine=SymPy_direct_series;imports_primary=no;formal_unknowns=endpoints,common_CE,x_compensator")
    print("implicit_J=compiler_derived;det=-288;uncompensated_rank=5")
    print(f"general_contact_s1={sp.factor(contact_1)}")
    print(f"general_contact_s2={sp.factor(contact_2)}=(9/8)*D2")
    print(f"zero_D2_contact_s3={sp.factor(contact_3)}=(3/4)*D4")
    print(f"zero_D2_D4_contact_s4={sp.factor(contact_4)}=(9/16)*delta6")
    print(f"exceptional_contact_s5={sp.factor(contact_5)}=(9/20)*delta8=-9/20*kappa8")
    print("constant_period_coefficients=-delta2,-delta4,-delta6,-delta8")
    print("adjoint_first_exit=n*contact_s^n=-(9/4)*constant_period_s^(n-1);n=2,3,4,5")
    print(f"contact_payload_sha256={hashlib.sha256(contact_payload.encode('ascii')).hexdigest()}")
    print("nested_strata=D2_zero_plane,D4_zero_parabola,D6_exceptional_quartic,D8_exit")
    print("scope=formal_local_labelled_triple_and_inherited_retained_debts_only")
    print("NO_CLAIM=all_order_identity_or_global_discriminant_or_target_lift_or_Keller_pair_or_JC2")
    print(f"checks={CHECKS};result=PASS")


if __name__ == "__main__":
    main()
