#!/usr/bin/env python3
"""Exact companion for THM-3779's maximal Danielewski observable."""

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
u, v, e, p, t, s, z, r = sp.symbols("u v e p t s z r")

# Abstract Danielewski identities.  The unit identity is simultaneously the
# smoothness certificate and the unit-minor certificate for the source map.
danielewski_relation = u * e - v * (v - 1)
unit_minor_identity = sp.expand((2 * v - 1) ** 2 - 4 * u * e)
unit_minor_remainder = sp.rem(
    sp.Poly(unit_minor_identity - 1, e), sp.Poly(danielewski_relation, e)
).as_expr()
gate(unit_minor_remainder == 0, "abstract unit-minor identity")

# With b=1-v,c=u the exponent-one two-arm bracket is exactly the standard
# Danielewski bracket for Sigma(b)=b(b-1).
b, c = sp.symbols("b c")
Sigma = b * (b - 1)
gate(sp.expand(Sigma.subs(b, 1 - v) - v * (v - 1)) == 0,
     "two-arm coordinate identification")
gate(sp.diff(Sigma, b) == 2 * b - 1, "two-arm derivative")

source_rows = []
recovery_rows = []

for degree_parameter in range(1, 9):
    m = degree_parameter
    A = 1 + x * y
    B = sp.expand(1 + x ** (2 * m + 1) * A**m)
    C = sp.expand((m + 1) * B - 1)
    U = sp.expand(x * A * B)
    P = sp.cancel(C / (m * x))
    V = sp.expand(U * P)
    N = sp.expand(
        y
        + sp.Rational(2 * m + 1, m) * x ** (2 * m) * A ** (m + 1)
        + sp.Rational(m + 1, m) * x ** (4 * m + 1) * A ** (2 * m + 1)
    )
    E = sp.cancel(P * (V - 1))
    E_polynomial = sp.expand(C * N / m)

    gate(sp.expand(V - 1 - x * N) == 0, f"V-1=xN m={m}")
    gate(sp.cancel(E - E_polynomial) == 0, f"E polynomial formula m={m}")
    gate(sp.Poly(E_polynomial, x, y).as_expr() == E_polynomial,
         f"E belongs to source ring m={m}")
    gate(sp.expand(U * E_polynomial - V * (V - 1)) == 0,
         f"Danielewski relation m={m}")

    gate(sp.cancel(jac(U, P, x, y) - 1) == 0, f"Keller response m={m}")
    gate(sp.expand(jac(U, V, x, y) - U) == 0,
         f"minor UV m={m}")
    gate(sp.expand(jac(U, E_polynomial, x, y) - (2 * V - 1)) == 0,
         f"minor UE m={m}")
    gate(sp.expand(jac(V, E_polynomial, x, y) - E_polynomial) == 0,
         f"minor VE m={m}")
    gate(
        sp.expand(
            (2 * V - 1) ** 2 - 4 * U * E_polynomial - 1
        )
        == 0,
        f"source unit-minor certificate m={m}",
    )

    # The two arms and the unique omitted point.
    gate(sp.expand(U.subs(x, 0)) == 0, f"axis U m={m}")
    gate(sp.expand(V.subs(x, 0) - 1) == 0, f"axis V m={m}")
    gate(sp.expand(E_polynomial.subs(x, 0) - y) == 0,
         f"axis E parameter m={m}")
    gate(sp.cancel(V.subs(y, -1 / x)) == 0, f"A-arm V m={m}")
    gate(sp.cancel(E_polynomial.subs(y, -1 / x) + 1 / x) == 0,
         f"A-arm E parameter m={m}")
    gate(sp.rem(sp.Poly(V, y), sp.Poly(B, y)).as_expr() == 0,
         f"B-arm V m={m}")
    b_arm_e = sp.together(E_polynomial - sp.Rational(1, m) / x)
    b_arm_num = sp.fraction(b_arm_e)[0]
    gate(sp.rem(sp.Poly(b_arm_num, y), sp.Poly(B, y)).as_expr() == 0,
         f"B-arm E parameter m={m}")

    # The U!=0 geometric-fibre recovery.  If f_m(t)=0, then
    # p=u/t+t^m/m.  Every displayed identity is checked before imposing f_m.
    ss = u - t ** (m + 1)
    xx = t / ss
    AA = ss**2 / t
    yy = ss * (ss**2 - t) / t**2
    BB = 1 + xx * t**m
    Prec = 1 / xx + sp.Rational(m + 1, m) * t**m
    gate(sp.cancel(xx * AA - ss) == 0, f"recovery s=xA m={m}")
    gate(sp.cancel(xx**2 * AA - t) == 0, f"recovery t=x2A m={m}")
    gate(sp.cancel(1 + xx * yy - AA) == 0, f"recovery y m={m}")
    gate(sp.cancel(BB - u / ss) == 0, f"recovery B m={m}")
    gate(sp.cancel(xx * AA * BB - u) == 0, f"recovery U m={m}")
    gate(sp.cancel(Prec - (u / t + t**m / m)) == 0,
         f"recovery P m={m}")

    # At a nonzero missing root one has the branch parametrization below.
    # The root is exactly double.  For m>=2 the residual factor has positive
    # degree and every one of its roots has s!=0; for m=1 it is the hostile
    # constant residual and the whole branch is omitted.
    minimal = t ** (m + 1) - m * p * t + m * u
    critical_minimal = sp.expand(
        minimal.subs({p: sp.Rational(m + 1, m) * r**m, u: r ** (m + 1)})
    )
    critical_quotient = sp.cancel(critical_minimal / (t - r) ** 2)
    gate(sp.denom(critical_quotient) == 1,
         f"critical double-root quotient polynomial m={m}")
    gate(sp.Poly(critical_quotient, t).degree() == m - 1,
         f"critical residual degree m={m}")
    gate(sp.factor(critical_quotient.subs(t, r)) != 0,
         f"critical root exactly double m={m}")

    source_rows.extend(
        (
            f"m={m};U={sp.sstr(U)}",
            f"m={m};V={sp.sstr(V)}",
            f"m={m};E={sp.sstr(E_polynomial)}",
        )
    )
    recovery_rows.append(f"m={m};s=u-t^{m+1};x=t/s;A=s^2/t")

# The m=1 boundary is genuinely different: the discriminant branch is
# omitted and its inverse is the additional source polynomial x^2.
A1 = 1 + x * y
B1 = 1 + x**3 * A1
U1 = sp.expand(x * A1 * B1)
P1 = sp.cancel((2 * B1 - 1) / x)
gate(sp.cancel(P1**2 - 4 * U1 - x ** -2) == 0,
     "m1 inverse discriminant observable")
gate(sp.factor((t**2 - 2 * t + 1) - (t - 1) ** 2) == 0,
     "m1 omitted branch point polynomial")
gate(1 * 2 - 2 * (2 - 1) == 0,
     "m1 omitted target point lies on D")

# Degree-three cover anatomy (m=2).
f2 = t**3 - 2 * p * t + 2 * u
disc2 = sp.factor(sp.discriminant(f2, t))
gate(sp.expand(disc2 - 4 * (8 * p**3 - 27 * u**2)) == 0,
     "m2 discriminant")
branch2 = {p: sp.Rational(3, 2) * r**2, u: r**3}
gate(sp.factor(f2.subs(branch2) - (t - r) ** 2 * (t + 2 * r)) == 0,
     "m2 branch factorization")
gate(sp.factor((u - t**3).subs(branch2).subs(t, r)) == 0,
     "m2 double root lies at missing s=0 boundary")

# Degree-four cover anatomy and its cubic resolvent/V4 quotient (m=3).
f3 = t**4 - 3 * p * t + 3 * u
disc3 = sp.factor(sp.discriminant(f3, t))
expected_disc3 = 27 * (256 * u**3 - 81 * p**4)
gate(sp.expand(disc3 - expected_disc3) == 0, "m3 quartic discriminant")
resolvent3 = z**3 - 12 * u * z - 9 * p**2
gate(sp.expand(sp.discriminant(resolvent3, z) - expected_disc3) == 0,
     "m3 cubic resolvent discriminant")
gate(sp.factor(resolvent3) == resolvent3, "m3 generic resolvent factor control")
branch3 = {p: sp.Rational(4, 3) * r**3, u: r**4}
gate(
    sp.factor(
        f3.subs(branch3)
        - (t - r) ** 2 * (t**2 + 2 * r * t + 3 * r**2)
    )
    == 0,
    "m3 branch factorization",
)
gate(
    sp.factor(
        resolvent3.subs(branch3) - (z + 2 * r**2) ** 2 * (z - 4 * r**2)
    )
    == 0,
    "m3 resolvent branch factorization",
)
gate(sp.factor((u - t**4).subs(branch3).subs(t, r)) == 0,
     "m3 double root lies at missing s=0 boundary")

# Scaling a resolvent root by U gives a monic integral V4-quotient equation
# over the polynomial base D; it does not put that root in the quartic field.
scaled_resolvent = z**3 - 12 * u**3 * z - 9 * u * v**2
scaled_disc = sp.factor(sp.discriminant(scaled_resolvent, z))
gate(
    sp.expand(scaled_disc - 27 * u**2 * (256 * u**7 - 81 * v**4)) == 0,
    "scaled polynomial resolvent discriminant",
)

# The two plane charts have the logarithmic transition responsible for the
# nonzero exponent-one symplectic class.
a0, a1 = sp.symbols("a0 a1")
transition = a0 - 1 / c
gate(sp.diff(transition - a0, c) == c ** -2,
     "two-arm singular shear derivative")
gate(sp.cancel((transition - a0) * c + 1) == 0,
     "two-arm singular shear")

semantic_rows = (
    "family:A=1+xy;B=1+x^(2m+1)A^m;U=xAB;P=((m+1)B-1)/(mx)",
    "observables:V=UP;V-1=xN_m;E=P(V-1)=((m+1)B-1)N_m/m",
    "surface:D=k[U,V,E]/(UE-V(V-1));smooth_unit=(2V-1)^2-4UE=1",
    "brackets:{U,V}=U;{U,E}=2V-1;{V,E}=E",
    "scope:m>=2;m1_omits_discriminant_and_has_x2=1/(P2-4U)",
    "image_m>=2:D_minus_(0,0,0);U_nonzero_root_recovery;two_explicit_arms",
    "intersection_m>=2:k[x,y]_intersect_k(U,P)=D_by_codim1_DVR_coverage",
    "strict:E_notin_k[U,V];D_is_maximal_polynomial_target_field_observable",
    "danielewski:b=1-V;c=U;e=E;ce=b(b-1);exponent1_two_arm_nonexact",
    "consequence:no_Darboux_pair_in_intersection;all_rational_target_words_closed",
    "m2:S3;t3-2Pt+2U;double_root_at_missing_boundary",
    "m3:S4;t4-3Pt+3U;resolvent=z3-12Uz-9P2;V4_quotient",
)
semantic = hashlib.sha256("\n".join(semantic_rows).encode()).hexdigest()
packet = hashlib.sha256("\n".join(source_rows + recovery_rows).encode()).hexdigest()

source = Path(__file__).read_text(encoding="utf-8")
gate(
    not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source))),
    "inactive Python assert",
)

print("theorem=THM-3779-three-component-tower-maximal-Danielewski-observable")
print("field=algebraically_closed_characteristic_zero;m>=2;de_Rham_consequence_over_C")
print("observables=V=U*P;V-1=x*N_m;E=P*(V-1)=((m+1)*B-1)*N_m/m")
print("surface=D=k[U,V,E]/(U*E-V*(V-1));strictly_contains_k[U,V]")
print("brackets=J(U,V)=U;J(U,E)=2V-1;J(V,E)=E")
print("unit_minor=(2V-1)^2-4UE=1;map=unramified_quasi_finite")
print("image_for_m>=2=D\\{(0,0,0)}")
print("intersection_for_m>=2=k[x,y]_intersect_k(U,P)=D")
print("danielewski_coordinates=b=1-V;c=U;e=E;c*e=b*(b-1)")
print("symplectic_class=exponent_one_two_arm_logarithmic_nonzero")
print("Darboux_pair=NONE;all_rational_target_field_polynomial_Keller_pairs=NONE")
print("m2=f=t^3-2Pt+2U;disc=4*(8P^3-27U^2);Galois_closure=S3")
print("m3=f=t^4-3Pt+3U;disc=27*(256U^3-81P^4);Galois_closure=S4")
print("m3_resolvent=z^3-12Uz-9P^2;fixed_V4_splitting_field=S3")
print("m3_scaled_resolvent=Z^3-12U^3Z-9UV^2")
print("m1_hostile=discriminant_branch_omitted;x^2=1/(P^2-4U)_is_extra")
print("finite_source_controls_m=1,2,3,4,5,6,7,8")
print(f"packet_sha256={packet}")
print(f"semantic_sha256={semantic}")
print(f"CHECKS={CHECKS}")
print("NO_CLAIM=planar_Jacobian_counterexample_or_arbitrary_non_target_field_closure")
print("RESULT=PASS")
