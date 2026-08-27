#!/usr/bin/env python3
"""Exact certificate for THM-4259's concrete W=0 glue dictionary.

The characteristic-zero part verifies the visible projection of H_lambda and
the pullback-differential formulas used to locate its hidden projection.  The
finite part does not infer a characteristic-zero vanishing from one modular
vanishing.  THM-4241 first places the hidden projection in L with degree 12;
the coercive bound and exact enumeration here produce the complete 24-element
shell.  Reduction modulo 397 is injective on that finite set, and an
invertible differential system selects one member.

The script freezes the identity

    h = H_lambda + v - omega^2*g,
    2h = v + omega^2*f + g.

It does not evaluate the remaining map-ratio incidences or prove JC(2).

The algebraic branch is the one selected by the simple triangular point over
F_397 checked below.  Hensel lifting supplies a characteristic-zero algebraic
point with that reduction; an embedding of its number field into C supplies
the complex maps.  This makes the later finite-candidate reduction a map on
the chosen characteristic-zero coefficient field, rather than an unrelated
finite-field test point.
"""

from __future__ import annotations

from itertools import product

import sympy as sp


E = tuple[int, int]  # m+n*omega, omega^2+omega+1=0

ZERO: E = (0, 0)
ONE: E = (1, 0)
OMEGA: E = (0, 1)
OMEGA2: E = (-1, -1)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def e_add(left: E, right: E) -> E:
    return left[0] + right[0], left[1] + right[1]


def e_neg(value: E) -> E:
    return -value[0], -value[1]


def e_sub(left: E, right: E) -> E:
    return left[0] - right[0], left[1] - right[1]


def e_mul(left: E, right: E) -> E:
    m, n = left
    r, s = right
    return m * r - n * s, m * s + n * r - n * s


def e_conjugate(value: E) -> E:
    return value[0] - value[1], -value[1]


def e_trace(value: E) -> int:
    return 2 * value[0] - value[1]


def e_norm(value: E) -> int:
    return value[0] * value[0] - value[0] * value[1] + value[1] * value[1]


def hidden_degree(left: E, right: E) -> int:
    cross = e_trace(e_mul(e_mul(left, e_conjugate(right)), (-4, -2)))
    return 6 * e_norm(left) + 6 * e_norm(right) + cross


def reduce_numerator(groebner: sp.GroebnerBasis, expression: sp.Expr) -> sp.Expr:
    numerator = sp.Poly(sp.together(expression).as_numer_denom()[0], *groebner.gens)
    return sp.factor(groebner.reduce(numerator.as_expr())[1])


def verify_visible_projection() -> tuple[str, str]:
    """Verify H_lambda+H_lambda o iota=-v in characteristic zero."""

    alpha, lam, r = sp.symbols("alpha lambda r")
    relations = sp.groebner(
        [
            r**2 - 12 * r + 24,
            lam**4 * (r - 9) - 1,
            alpha * lam * (8 - r) + 4,
        ],
        alpha, lam, r,
        order="lex",
        domain=sp.QQ,
    )

    require(reduce_numerator(relations, alpha**3 - r * lam) == 0,
            "normalized alpha lost alpha^3=r*lambda")
    require(reduce_numerator(relations, alpha**2 - r * lam**2 + 6 * lam**2) == 0,
            "visible Y constant no longer cancels")

    # Adding H(y) and H(-y) on Y^2=X^3+1 has slope
    # -epsilon*alpha^2/(2x).  Its X coefficient and residual Y constant are:
    x_coefficient = 2 * alpha * lam - alpha**4 / 4
    y_constant_twice = alpha**2 - alpha**3 * lam + 6 * lam**2
    require(reduce_numerator(relations, x_coefficient + 1) == 0,
            "visible projection X is not -x^-2")
    require(reduce_numerator(relations, y_constant_twice) == 0,
            "visible projection Y retained a constant term")
    return "X=-x^-2", "Y=-epsilon*y^2/x^3"


def verify_hidden_differentials() -> tuple[str, str, str]:
    """Verify the f,g and anti-invariant H_lambda differential formulas."""

    a, y = sp.symbols("a y")
    quartic = a**4 - 2 * a**3 - 2 * a + 1
    quartic_basis = sp.groebner([quartic], a, y, order="lex", domain=sp.QQ)

    left_x = 1 - a**2
    right_x = 1 + a**2
    left_y = 1 + a**3
    right_y = 1 - a**3
    a_zero = -a**3 + 2 * a**2 + 3
    b_zero = a**3 - 2 * a**2 - 1

    differential_remainder = quartic_basis.reduce(sp.expand(
        3 * right_x + 2 * left_x * y**2 - right_x * y**4
        - (a_zero + b_zero * y**2) * (left_y + right_y * y**2)
    ))[1]
    require(sp.factor(differential_remainder) == 0,
            "hidden f differential quotient changed")

    # The quartic root gives t=a+a^-1 with t^2-2t-2=0 and
    # r=4+2t, hence r^2-12r+24=0.
    r_from_a = 4 + 2 * (a + 1 / a)
    r_relation = sp.together(r_from_a**2 - 12 * r_from_a + 24)
    require(sp.rem(
        sp.Poly(r_relation.as_numer_denom()[0], a), sp.Poly(quartic, a)
    ).as_expr() == 0, "hidden quartic no longer supplies the H-lambda r")

    lam, r, c = sp.symbols("lambda r c")
    q_plus = y**2 - c * y / 2 + 3 * lam**2
    q_minus = y**2 + c * y / 2 + 3 * lam**2
    f_plus = (3 + y**4 - 4 * lam * y**3) / q_plus
    f_minus = (3 + y**4 + 4 * lam * y**3) / q_minus
    h_relations = sp.groebner(
        [c - r * lam, r**2 - 12 * r + 24, lam**4 * (r - 9) - 1],
        c, r, lam, y,
        order="lex",
        domain=sp.QQ,
    )
    target_even = 2 * (y**2 + lam**2 * (r - 9))
    require(reduce_numerator(h_relations, f_plus + f_minus - target_even) == 0,
            "anti-invariant H-lambda differential changed")

    return (
        "eta_f=(2/(3s))*x^-5*(A0+B0*y^2)dy",
        "eta_g=kappa*(2/(3s))*x^-5*(A0-B0*y^2)dy",
        "eta_ell=(2alpha/(3epsilon))*x^-5*(y^2+lambda^2(r-9))dy",
    )


def enumerate_hidden_degree_twelve() -> set[tuple[E, E]]:
    # The Hermitian off-diagonal has norm 12, so Cauchy--Schwarz gives
    #
    #   q(Af+Bg) >= (6-2*sqrt(3))*(N(A)+N(B)).
    #
    # If q=12 and the integer S=N(A)+N(B) were at least 5, then the
    # right side would be at least 30-10*sqrt(3)>12 (square 18>10*sqrt(3)).
    # Hence S<=4.  Also N(m+n*omega)=(m-n/2)^2+3n^2/4, and by symmetry
    # N(m+n*omega)=(n-m/2)^2+3m^2/4, so norm at most four forces
    # |m|,|n|<=2.  Thus the following is a proved complete universe, not a
    # boundary heuristic.
    eisenstein_ball = {
        (m, n)
        for m, n in product(range(-2, 3), repeat=2)
        if e_norm((m, n)) <= 4
    }
    vectors = {
        (left, right)
        for left, right in product(eisenstein_ball, repeat=2)
        if e_norm(left) + e_norm(right) <= 4
        and hidden_degree(left, right) == 12
    }
    require(len(vectors) == 24, "degree-12 hidden shell changed")
    require(max(e_norm(left) + e_norm(right) for left, right in vectors) <= 4,
            "degree-12 vector escaped the proved norm bound")
    return vectors


def verify_good_prime_and_lift() -> dict[str, object]:
    q = 397
    require(bool(sp.isprime(q)), "the reduction modulus is not prime")
    zeta, a, scale = 157, 161, 27
    r, lam, alpha, epsilon = 363, 30, 92, 334
    inv = lambda value: pow(value % q, -1, q)

    require((zeta**4 - zeta**2 + 1) % q == 0, "zeta_12 equation failed")
    require(pow(zeta, 12, q) == 1
            and all(pow(zeta, divisor, q) != 1 for divisor in (1, 2, 3, 4, 6)),
            "zeta lost exact order twelve")
    omega = pow(zeta, 4, q)
    kappa = pow(zeta, 5, q)
    require((omega * omega + omega + 1) % q == 0 and omega != 1,
            "Eisenstein embedding changed")
    require((kappa * kappa + omega) % q == 0, "kappa^2=-omega changed")

    quartic = a**4 - 2 * a**3 - 2 * a + 1
    scale_denominator = 2 * a**3 + 3 * a**2 - 1
    require(quartic % q == 0, "hidden quartic root changed")
    require(scale_denominator % q != 0
            and (scale**6 * scale_denominator - 4) % q == 0,
            "hidden scale relation changed")
    require((r - (4 + 2 * (a + inv(a)))) % q == 0,
            "r no longer matches the hidden root")
    require((r**2 - 12 * r + 24) % q == 0, "r quadratic changed")
    require((lam**4 * (r - 9) - 1) % q == 0, "lambda relation changed")
    require((alpha * lam * (8 - r) + 4) % q == 0,
            "alpha normalization changed")
    require((alpha**3 - r * lam) % q == 0, "alpha cube changed")
    require(epsilon == pow(zeta, 3, q) and (epsilon**2 + 1) % q == 0,
            "epsilon=i embedding changed")

    simple_derivatives = (
        4 * zeta**3 - 2 * zeta,
        4 * a**3 - 6 * a**2 - 2,
        6 * scale**5 * scale_denominator,
        4 * lam**3 * (r - 9),
        lam * (8 - r),
    )
    require(all(value % q != 0 for value in simple_derivatives),
            "good-prime parameter system lost a simple coordinate")

    a_zero = (-a**3 + 2 * a**2 + 3) % q
    b_zero = (a**3 - 2 * a**2 - 1) % q
    determinant = (-2 * kappa * a_zero * b_zero) % q
    require(determinant != 0, "differential coefficient system became singular")

    right_constant = (
        scale * alpha * inv(epsilon) * lam**2 * (r - 9) * inv(a_zero)
    ) % q
    right_quadratic = (scale * alpha * inv(epsilon) * inv(b_zero)) % q
    coefficient_f = ((right_constant + right_quadratic) * inv(2)) % q
    coefficient_g = (
        (right_constant - right_quadratic) * inv(2) * inv(kappa)
    ) % q
    require((coefficient_f, coefficient_g) == (362, 328),
            "modular hidden coefficients changed")

    vectors = enumerate_hidden_degree_twelve()
    reduce_e = lambda value: (value[0] + value[1] * omega) % q
    reductions = {
        (reduce_e(left), reduce_e(right)): (left, right)
        for left, right in vectors
    }
    require(len(reductions) == 24, "degree-12 reductions collided")
    require((coefficient_f, coefficient_g) in reductions,
            "modular coefficients missed the exact degree-12 shell")
    lift = reductions[(coefficient_f, coefficient_g)]
    expected_lift = (OMEGA2, e_sub(OMEGA2, OMEGA))
    require(lift == expected_lift, "hidden coefficient lift changed")
    require(hidden_degree(*lift) == 12, "lift lost hidden degree twelve")

    return {
        "q": q,
        "parameters": (zeta, a, scale, r, lam, alpha, epsilon),
        "omega": omega,
        "kappa": kappa,
        "system_determinant": determinant,
        "residue_coefficients": (coefficient_f, coefficient_g),
        "shell_size": len(vectors),
        "distinct_reductions": len(reductions),
        "lift": lift,
    }


def verify_normalized_dictionary() -> tuple[str, str]:
    # 2H=-v+omega^2*f+(omega^2-omega)*g and
    # h=H+v-omega^2*g imply 2h=v+omega^2*f+g.
    hidden_g = e_sub(OMEGA2, OMEGA)
    doubled_h_g = e_add(hidden_g, e_neg(e_add(OMEGA2, OMEGA2)))
    require(doubled_h_g == ONE, "normalized h relation lost its g coefficient")

    # The relation also forces Th=omega^2*h-omega*f from
    # Tv=omega^2*v, Tf=g, Tg=-omega*f. Compare after doubling.
    t_doubled_h = (OMEGA2, e_neg(OMEGA), OMEGA2)  # coefficients of (v,f,g)
    rhs_doubled_h = (
        OMEGA2,
        e_sub(e_mul(OMEGA2, OMEGA2), e_add(OMEGA, OMEGA)),
        OMEGA2,
    )
    require(t_doubled_h == rhs_doubled_h, "normalized T action changed")
    return "h=H_lambda+v-omega^2*g", "T*h=omega^2*h-omega*f"


def main() -> None:
    visible = verify_visible_projection()
    differentials = verify_hidden_differentials()
    finite = verify_good_prime_and_lift()
    dictionary, t_action = verify_normalized_dictionary()

    print("THM-4259 W=0 explicit H_lambda glue dictionary exact certificate")
    print("symbolic_visible_projection=H_lambda+H_lambda_o_iota=-v PASS "
          f"coordinates={visible}")
    print(f"symbolic_hidden_differentials=PASS formulas={differentials}")
    print(f"good_prime={finite['q']} parameters(zeta12,a,s,r,lambda,alpha,epsilon)={finite['parameters']} omega={finite['omega']} kappa={finite['kappa']}")
    print(f"differential_system_determinant={finite['system_determinant']} nonzero=PASS residue_coefficients={finite['residue_coefficients']}")
    print(f"hidden_degree12_shell={finite['shell_size']} distinct_mod397={finite['distinct_reductions']} unique_lift={finite['lift']}")
    print("hidden_projection=H_lambda-H_lambda_o_iota=omega^2*f+(omega^2-omega)*g")
    print("doubled_map=2*H_lambda=-v+omega^2*f+(omega^2-omega)*g")
    print(f"normalized_dictionary={dictionary} relation=2h=v+omega^2*f+g")
    print(f"normalized_T_action={t_action} PASS")
    print("modular_logic=FINITE_CANDIDATE_IDENTIFICATION_NOT_ZERO_LIFT")
    print("normalization_blocker=REMOVED mixed_1512_incidences=UNEVALUATED")
    print("verdict=CONCRETE_GLUE_GENERATOR_PROVED JC2_OPEN W0_OPEN")


if __name__ == "__main__":
    main()
