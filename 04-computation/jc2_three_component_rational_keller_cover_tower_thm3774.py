#!/usr/bin/env python3
"""Exact companion for THM-3774's rational Keller cover tower."""

from __future__ import annotations

import ast
import hashlib
from pathlib import Path

import sympy as sp


CHECKS = 0


def gate(condition: bool, message: str) -> None:
    global CHECKS
    CHECKS += 1
    if condition is not True and condition != sp.S.true:
        raise RuntimeError(message)


def jac(first: sp.Expr, second: sp.Expr, left: sp.Symbol, right: sp.Symbol) -> sp.Expr:
    return sp.expand(
        sp.diff(first, left) * sp.diff(second, right)
        - sp.diff(first, right) * sp.diff(second, left)
    )


x, y = sp.symbols("x y")
X, A = sp.symbols("X A", nonzero=True)
T, Ptar, L, r = sp.symbols("T P L r", nonzero=True)
m, q = sp.symbols("m q", integer=True, positive=True)

# Work first in the etale chart (x,A), where A=1+xy and
# J_(x,y)(F,G)=x J_(x,A)(F,G).  Keeping m and q symbolic proves the
# all-degree and off-exponent identities rather than interpolating them.
Bq = 1 + X**q * A**m
Uq = X * A * Bq
Pq = ((m + 1) * Bq - 1) / (m * X)
chart_response = sp.factor(
    X * (sp.diff(Uq, X) * sp.diff(Pq, A) - sp.diff(Uq, A) * sp.diff(Pq, X))
)
response_error = sp.factor(chart_response - 1)
expected_error = sp.factor(
    -(m + 1) * (q - 2 * m - 1) * X**q * A**m * Bq / m
)
gate(sp.factor(response_error - expected_error) == 0, "symbolic exponent gate")
gate(
    sp.factor(
        response_error / Uq
        + (m + 1) * (q - 2 * m - 1) * X ** (q - 1) * A ** (m - 1) / m
    )
    == 0,
    "symbolic U-multiple gate",
)

q_exact = 2 * m + 1
B_exact = Bq.subs(q, q_exact)
U_exact = Uq.subs(q, q_exact)
P_exact = Pq.subs(q, q_exact)
gate(sp.factor(chart_response.subs(q, q_exact) - 1) == 0, "uniform response")
gate(
    sp.factor(sp.diff(U_exact, X) - A * (1 + 2 * (m + 1) * (B_exact - 1))) == 0,
    "U_X chart derivative",
)
gate(
    sp.factor(sp.diff(U_exact, A) - X * (1 + (m + 1) * (B_exact - 1))) == 0,
    "U_A chart derivative",
)
gate(
    sp.factor(sp.diff(P_exact, X) - (-1 + 2 * (m + 1) * (B_exact - 1)) / X**2)
    == 0,
    "P_X chart derivative",
)
gate(
    sp.factor(sp.diff(P_exact, A) - (m + 1) * (B_exact - 1) / (X * A))
    == 0,
    "P_A chart derivative",
)
s_chart = X * A
t_chart = X**2 * A
U_birational = s_chart + t_chart ** (m + 1)
gate(
    sp.factor(
        X
        * (
            sp.diff(t_chart, X) * sp.diff(U_birational, A)
            - sp.diff(t_chart, A) * sp.diff(U_birational, X)
        )
        - X * t_chart
    )
    == 0,
    "rational constant-field derivation",
)

direct_controls = []
resultant_controls = []
discriminant_controls = []
irreducibility_controls = []
hostile_controls = []
packet_rows = []
Axy = 1 + x * y
t_source = x**2 * Axy
s_source = x * Axy

for degree_parameter in range(1, 8):
    exponent = 2 * degree_parameter + 1
    Bxy = sp.expand(1 + x**exponent * Axy**degree_parameter)
    Uxy = sp.expand(x * Axy * Bxy)
    Pxy = sp.cancel(
        ((degree_parameter + 1) * Bxy - 1) / (degree_parameter * x)
    )
    Wxy = sp.expand(Uxy * Pxy)

    gate(sp.cancel(jac(Uxy, Pxy, x, y) - 1) == 0,
         f"direct source response m={degree_parameter}")
    gate(sp.expand(jac(Uxy, Wxy, x, y) - Uxy) == 0,
         f"direct log-canonical response m={degree_parameter}")
    gate(
        sp.expand(Uxy - s_source - t_source ** (degree_parameter + 1)) == 0,
        f"generic-fibre s,t identity m={degree_parameter}",
    )
    gate(
        sp.cancel(Pxy - Uxy / t_source - t_source**degree_parameter / degree_parameter)
        == 0,
        f"generic-fibre P identity m={degree_parameter}",
    )
    gate(
        sp.cancel(
            t_source ** (degree_parameter + 1)
            - degree_parameter * Pxy * t_source
            + degree_parameter * Uxy
        )
        == 0,
        f"minimal relation on source m={degree_parameter}",
    )

    # The affine zero fibre has the three displayed addresses.  W=UP records
    # their first vertical coefficients without choosing local parameters.
    gate(sp.expand(Wxy.subs(x, 0)) == 1, f"axis spectrum m={degree_parameter}")
    gate(sp.cancel(Wxy.subs(y, -1 / x)) == 0,
         f"A-component spectrum m={degree_parameter}")
    gate(sp.rem(sp.Poly(Wxy, y), sp.Poly(Bxy, y)).as_expr() == 0,
         f"B-component spectrum m={degree_parameter}")
    gate(sp.expand(sp.diff(Uxy, x).subs(x, 0)) == 1,
         f"safe axis m={degree_parameter}")

    # Finite irreducibility controls accompany the all-m primitive-binomial
    # proof.  The specialization at P=0 is the Eisenstein polynomial used in
    # the theorem's generic irreducibility argument.
    b_factorization = sp.factor_list(Bxy)[1]
    gate(len(b_factorization) == 1 and b_factorization[0][1] == 1,
         f"B irreducibility control m={degree_parameter}")
    eisenstein_specialization = T ** (degree_parameter + 1) + degree_parameter * L
    gate(sp.factor(eisenstein_specialization) == eisenstein_specialization,
         f"Eisenstein specialization control m={degree_parameter}")
    irreducibility_controls.append(str(degree_parameter))

    minimal_polynomial = (
        T ** (degree_parameter + 1)
        - degree_parameter * Ptar * T
        + degree_parameter * L
    )
    discriminant = sp.factor(sp.discriminant(minimal_polynomial, T))
    expected_discriminant = sp.factor(
        (-1) ** (degree_parameter * (degree_parameter + 1) // 2)
        * degree_parameter**degree_parameter
        * (
            (degree_parameter + 1) ** (degree_parameter + 1) * L**degree_parameter
            - degree_parameter ** (degree_parameter + 1)
            * Ptar ** (degree_parameter + 1)
        )
    )
    gate(sp.expand(discriminant - expected_discriminant) == 0,
         f"discriminant formula m={degree_parameter}")
    branch_factorization = sp.factor_list(
        (degree_parameter + 1) ** (degree_parameter + 1) * L**degree_parameter
        - degree_parameter ** (degree_parameter + 1) * Ptar ** (degree_parameter + 1)
    )[1]
    gate(len(branch_factorization) == 1 and branch_factorization[0][1] == 1,
         f"branch irreducibility control m={degree_parameter}")
    discriminant_controls.append(f"m{degree_parameter}:{sp.sstr(discriminant)}")

    # q=2m and q=2m+2 are the nearest hostile exponents.
    for hostile_q in (2 * degree_parameter, 2 * degree_parameter + 2):
        hostile_b = 1 + x**hostile_q * Axy**degree_parameter
        hostile_u = sp.expand(x * Axy * hostile_b)
        hostile_p = sp.cancel(
            ((degree_parameter + 1) * hostile_b - 1) / (degree_parameter * x)
        )
        hostile_error = sp.factor(jac(hostile_u, hostile_p, x, y) - 1)
        expected_hostile_error = sp.factor(
            -(degree_parameter + 1)
            * (hostile_q - 2 * degree_parameter - 1)
            * x**hostile_q
            * Axy**degree_parameter
            * hostile_b
            / degree_parameter
        )
        gate(sp.expand(hostile_error - expected_hostile_error) == 0,
             f"hostile exponent m={degree_parameter},q={hostile_q}")
        gate(hostile_error != 0,
             f"hostile exponent nonzero m={degree_parameter},q={hostile_q}")
        hostile_controls.append(f"m{degree_parameter}q{hostile_q}")

    if degree_parameter <= 5:
        gradient_resultant = sp.factor(
            sp.resultant(sp.diff(Uxy, x), sp.diff(Uxy, y), y)
        )
        expected_resultant = (
            (degree_parameter + 1) ** degree_parameter
            * x ** (3 * degree_parameter**2 + 4 * degree_parameter + 2)
        )
        gate(gradient_resultant == expected_resultant,
             f"gradient resultant m={degree_parameter}")
        resultant_controls.append(f"m{degree_parameter}:{sp.sstr(gradient_resultant)}")

    # Killing the constant term in A equalizes only by creating a critical
    # origin; deleting the nonlinear B term recovers the two-component seed.
    A_degenerate = x * y
    B_degenerate = 1 + x**exponent * A_degenerate**degree_parameter
    U_degenerate = sp.expand(x * A_degenerate * B_degenerate)
    gate(
        tuple(
            sp.diff(U_degenerate, variable).subs({x: 0, y: 0})
            for variable in (x, y)
        )
        == (0, 0),
        f"critical degeneration m={degree_parameter}",
    )
    two_component_u = x * Axy
    two_component_p = 1 / x
    gate(sp.cancel(jac(two_component_u, two_component_p, x, y) - 1) == 0,
         f"two-component boundary m={degree_parameter}")

    direct_controls.append(str(degree_parameter))
    packet_rows.extend(
        (
            f"m={degree_parameter};U={sp.sstr(Uxy)}",
            f"m={degree_parameter};W={sp.sstr(Wxy)}",
            f"m={degree_parameter};f={sp.sstr(minimal_polynomial)}",
            f"m={degree_parameter};disc={sp.sstr(discriminant)}",
        )
    )

# The degree-three member is the first non-Galois cover.  Its splitting-field
# group is S3 because the cubic is irreducible and its discriminant has a
# squarefree cuspidal factor.
f2 = T**3 - 2 * Ptar * T + 2 * L
disc2 = sp.factor(sp.discriminant(f2, T))
gate(sp.expand(disc2 + 4 * (27 * L**2 - 8 * Ptar**3)) == 0,
     "m=2 cuspidal discriminant")
gate(sp.expand(sp.factor(27 * L**2 - 8 * Ptar**3) - (27 * L**2 - 8 * Ptar**3)) == 0,
     "m=2 cusp irreducibility control")
gate(
    tuple(
        sp.diff(27 * L**2 - 8 * Ptar**3, variable).subs({L: 0, Ptar: 0})
        for variable in (L, Ptar)
    )
    == (0, 0),
    "m=2 cusp singular origin",
)

# The normalization of the branch binomial exhibits one double root.  The
# theorem turns this exact calculation into the transposition-inertia proof
# of non-Galoisness for every m>=2.
f_symbolic = T ** (m + 1) - m * Ptar * T + m * L
branch_substitution = {Ptar: (m + 1) * r**m / m, L: r ** (m + 1), T: r}
gate(sp.factor(f_symbolic.subs(branch_substitution)) == 0,
     "branch double-root equation")
gate(sp.factor(sp.diff(f_symbolic, T).subs(branch_substitution)) == 0,
     "branch derivative equation")
gate(
    sp.factor(sp.diff(f_symbolic, T, 2).subs(branch_substitution)
              - m * (m + 1) * r ** (m - 1))
    == 0,
    "branch ordinary-double-root equation",
)

semantic_rows = (
    "family:A=1+xy;B=1+x^(2m+1)A^m;U=xAB;P=((m+1)B-1)/(mx);W=UP",
    "response:J(U,P)=1;J(U,W)=U;exponent_unique:q=2m+1",
    "zero_fibre:three_smooth_disjoint_components;x,A,B;spectrum=(1,0,0)",
    "torsor:P+k(U);vertical_equalizer_fails_at_one_axis_simple_coefficient",
    "birational_chart:s=xA;t=x^2A;U=s+t^(m+1);P=U/t+t^m/m",
    "minimal:t^(m+1)-mPt+mU;degree=m+1;irreducible_by_P0_Eisenstein",
    "disc:sign*m^m*((m+1)^(m+1)U^m-m^(m+1)P^(m+1))",
    "m2:cubic_degree3;S3_closure;A2_cusp=27U^2-8P^3;not_polynomial",
)
semantic = hashlib.sha256("\n".join(semantic_rows).encode()).hexdigest()
packet = hashlib.sha256("\n".join(packet_rows).encode()).hexdigest()

source = Path(__file__).read_text(encoding="utf-8")
gate(
    not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source))),
    "inactive Python assert",
)

print("theorem=THM-3774-three-component-rational-keller-cover-tower")
print("field=algebraically_closed_characteristic_zero;m_integer_at_least_1")
print("family=A=1+x*y;B=1+x^(2m+1)*A^m;U=x*A*B;P=((m+1)*B-1)/(m*x);W=U*P")
print("uniform_response=J(U,P)=1;log_canonical=J(U,W)=U")
print("exponent_gate=J(U_q,P_q)-1=-(m+1)*(q-2m-1)*x^q*A^m*B_q/m")
print("zero_fibre=three_smooth_pairwise_disjoint_irreducibles_(x,A,B)")
print("component_spectrum_W=(1,0,0);principal_coefficients_P_in_1/U=(1,0,0)")
print("complete_rational_torsor=P+k(U);vertical_equalizer=FAIL_at_axis_simple_row")
print("generic_fibre=k(U)(t);s=x*A;t=x^2*A;U=s+t^(m+1);P=U/t+t^m/m")
print("minimal_polynomial=t^(m+1)-m*P*t+m*U;function_field_degree=m+1")
print("discriminant=(-1)^(m*(m+1)/2)*m^m*((m+1)^(m+1)*U^m-m^(m+1)*P^(m+1))")
print("monodromy=m>=2_non_Galois_by_generic_transposition_inertia")
print("m2=t^3-2*P*t+2*U;disc=-4*(27*U^2-8*P^3);closure_group=S3;cusp=A2")
print("blowdown=W=U*P_clears_axis_pole_but_retains_spectrum_(1,0,0)")
print("hostiles=q!=2m+1_nonconstant_response;A=xy_critical_origin;B=1_two_components")
print("direct_controls_m=" + ",".join(direct_controls))
print("irreducibility_controls_m=" + ",".join(irreducibility_controls))
print("resultant_controls=" + ",".join(resultant_controls))
print("discriminant_controls=" + ",".join(discriminant_controls))
print("hostile_exponent_controls=" + ",".join(hostile_controls))
print(f"packet_sha256={packet}")
print(f"semantic_sha256={semantic}")
print(f"CHECKS={CHECKS}")
print("NO_CLAIM=polynomial_Keller_pair_or_planar_Jacobian_counterexample")
print("RESULT=PASS")
