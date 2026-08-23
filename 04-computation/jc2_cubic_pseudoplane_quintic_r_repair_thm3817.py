#!/usr/bin/env python3
"""Exact companion for THM-3817's uniform quintic r-repair obstruction."""

from __future__ import annotations

import ast
import hashlib
import json
from pathlib import Path

import sympy as sp


CHECKS = 0


def gate(condition: bool, message: str) -> None:
    global CHECKS
    CHECKS += 1
    if condition is not True and condition != sp.S.true:
        raise RuntimeError(message)


def zero(expression: sp.Expr, message: str) -> None:
    gate(sp.cancel(expression) == 0, message)


r, z, e, u = sp.symbols("r z e u")
b0, b1, b2, b3, b4, b5 = sp.symbols("b0 b1 b2 b3 b4 b5")
variables = (r, z, e)
surface = r**2 * e - z**3 + r
poisson = sp.Matrix(
    [
        [0, 3 * r**2, 9 * z**2],
        [-3 * r**2, 0, 3 + 6 * r * e],
        [-9 * z**2, -3 - 6 * r * e, 0],
    ]
)


def bracket(left: sp.Expr, right: sp.Expr) -> sp.Expr:
    dl = sp.Matrix([sp.diff(left, q) for q in variables])
    dr = sp.Matrix([sp.diff(right, q) for q in variables])
    return sp.expand((dl.T * poisson * dr)[0])


g = b0 + b1 * e + b2 * e**2 + b3 * e**3 + b4 * e**4 + b5 * e**5
gp = sp.diff(g, e)
A = e**2 - z / 3 + r * g
K_source = 1 + 2 * r * e
critical = [sp.factor(bracket(A, q)) for q in variables]
zero(
    critical[0] - (r**2 - 9 * z**2 * (2 * e + r * gp)),
    "quintic r Hamiltonian",
)
zero(
    critical[1] - (3 * g * r**2 - 3 * K_source * (2 * e + r * gp)),
    "quintic z Hamiltonian",
)
zero(critical[2] - (9 * g * z**2 - K_source), "quintic e Hamiltonian")
zero(bracket(A, surface), "surface Casimir")


# Universal two-equation compression and residual resultant.
G, D = sp.symbols("G D")
K = 1 + 2 * u
P = sp.expand(G * u**2 - K * (2 * e**3 + u * e * D))
Q = sp.expand(e**2 * K**3 - 729 * G**3 * u**2 * (1 + u) ** 2)
resultant_generic = sp.factor(sp.resultant(P, Q, u))
H_generic = sp.factor(resultant_generic / (G**3 * e**4))
zero(resultant_generic - G**3 * e**4 * H_generic, "universal resultant factor")

H = sp.Poly(sp.expand(H_generic.subs({G: g, D: gp})), e)
gate(H.degree() == 27, "quintic residual degree twenty-seven")
gate(H.LC() == 34012224 * b5**5, "quintic residual leading coefficient")
gate(
    H.as_expr().subs({b0: 2, b1: -3, b2: 5, b3: -7, b4: 11, b5: -13})
    != 0,
    "positive residual control",
)


# Boundary-only support, including repeated roots, forces H|e*g*H'.  Compute
# the quotient over the parameter fraction field and record its exact global
# denominator before any specialization.
quotient, remainder = sp.div(
    sp.Poly(sp.expand(e * g * sp.diff(H.as_expr(), e)), e, domain="EX"),
    sp.Poly(H.as_expr(), e, domain="EX"),
)
q_numerator, q_denominator = sp.fraction(sp.together(quotient.as_expr()))
zero(q_denominator - 512 * b5**4, "quintic quotient denominator")
gate(quotient.degree() == 5, "quintic logarithmic quotient degree")
gate(quotient.LC() == 27 * b5, "quintic logarithmic quotient lead")
R = sp.Poly(
    sp.expand(q_denominator * e * g * sp.diff(H.as_expr(), e)
              - q_numerator * H.as_expr()),
    e,
)
zero(R.as_expr() - q_denominator * remainder.as_expr(), "cleared remainder")
gate(R.degree() == 26, "cleared quintic remainder degree")


# The quintic leading coefficient has grading weight -2.  Over the
# algebraically closed ground field choose s with b5=s^2.  This finite cover
# is surjective on b5!=0 and restores integral weight-zero coordinates.
C0, C1, C2, C3, C4, T, s = sp.symbols("C0 C1 C2 C3 C4 T s")
parameter_pullback = {
    b0: C0 / s**3,
    b1: C1 / s**2,
    b2: C2 / s,
    b3: C3,
    b4: C4 * s,
    b5: s**2,
}


def invariant_primitive(expression: sp.Expr) -> tuple[int, sp.Expr]:
    pulled = sp.expand(sp.cancel(expression.subs(parameter_pullback)))
    terms = sp.Add.make_args(pulled)
    powers = [int(term.as_powers_dict().get(s, 0)) for term in terms]
    minimum = min(powers)
    normalized = sp.expand(pulled / s**minimum)
    invariant = sp.Integer(0)
    for term in sp.Add.make_args(normalized):
        power = int(term.as_powers_dict().get(s, 0))
        gate(power % 7 == 0, "quintic pullback has one mod-seven class")
        invariant += term / s**power * T ** (power // 7)
    polynomial = sp.Poly(sp.expand(invariant), C0, C1, C2, C3, C4, T)
    content, primitive = polynomial.primitive()
    gate(content != 0, "nonzero quintic invariant primitive content")
    return minimum, primitive.as_expr()


top_degrees = tuple(range(26, 14, -1))
top_packet: list[sp.Expr] = []
laurent_powers: list[int] = []
for degree in top_degrees:
    power, primitive = invariant_primitive(R.coeff_monomial(e**degree))
    laurent_powers.append(power)
    top_packet.append(primitive)
gate(laurent_powers == list(range(14, 2, -1)), "quintic Laurent powers")

# Verify directly that primitive normalization loses no equation.
for degree, power, primitive in zip(top_degrees, laurent_powers, top_packet):
    raw_pullback = sp.expand(
        sp.cancel(R.coeff_monomial(e**degree).subs(parameter_pullback) / s**power)
    )
    primitive_lift = sp.expand(primitive.subs(T, s**7))
    ratio = sp.cancel(raw_pullback / primitive_lift)
    gate(
        not ratio.has(C0, C1, C2, C3, C4, s),
        "quintic primitive differs by scalar only",
    )
    gate(ratio != 0, "quintic primitive scalar is nonzero")


# The first five equations form the fixed exact ideal I.  Reduce the next
# seven rows against one Groebner basis rather than repeatedly recomputing
# enlarged ideals.
I_basis = sp.groebner(
    top_packet[:5], C0, C1, C2, C3, C4, T, order="grevlex"
)
gate(len(I_basis.polys) == 149, "quintic five-row Groebner basis size")
normal = {
    degree: sp.factor(I_basis.reduce(primitive)[1])
    for degree, primitive in zip(top_degrees[5:], top_packet[5:])
}

expected_21 = -512 * T * (210 * C3 - 83 * C4**2 + 678)
expected_20 = -128 * T * (
    6060 * C2 + 5322 * C3 * C4 - 2567 * C4**3 + 25238 * C4
)
expected_19 = -128 * T * (
    9920 * C1 + 16317 * C2 * C4 + 4548 * C3**2
    + 436 * C3 * C4**2 + 16256 * C3 - 2201 * C4**4
    + 30018 * C4**2 - 9680
)
expected_18 = -sp.Rational(256, 153) * T * (
    4088250 * C0 + 1340840 * C1 * C4 + 693900 * C2 * C3
    + 1400082 * C2 * C4**2 + 666198 * C2 + 732067 * C3**2 * C4
    - 642798 * C3 * C4**3 + 2909169 * C3 * C4
    + 1283409 * C4**3 - 1975219 * C4
)
zero(normal[21] - expected_21, "quintic row21 normal form")
zero(normal[20] - expected_20, "quintic row20 normal form")
zero(normal[19] - expected_19, "quintic row19 normal form")
zero(normal[18] - expected_18, "quintic row18 normal form")

# Since T=s^7!=0, the four rows solve triangularly without division by a
# profile parameter.
C3_solution = (83 * C4**2 - 678) / 210
C2_solution = 26 * C4 * (156 * C4**2 - 2711) / 53025
C1_solution = (
    13008795 * C4**4 - 322587262 * C4**2 + 2738675802
) / 1841028000
C0_solution = C4 * (
    148537264577 * C4**4 - 4086004213596 * C4**2 - 24178143226221
) / 282246852037500
triangular_substitution = {
    C3: C3_solution,
    C2: C2_solution,
    C1: C1_solution,
    C0: C0_solution,
}
zero(normal[21].subs(C3, C3_solution), "solve C3")
zero(normal[20].subs({C3: C3_solution, C2: C2_solution}), "solve C2")
zero(
    normal[19].subs({C3: C3_solution, C2: C2_solution, C1: C1_solution}),
    "solve C1",
)
zero(normal[18].subs(triangular_substitution), "solve C0")


def branch_primitive(degree: int) -> sp.Expr:
    numerator = sp.together(normal[degree].subs(triangular_substitution)).as_numer_denom()[0]
    return sp.Poly(numerator, C4, T).primitive()[1].as_expr()


P17 = branch_primitive(17)
P16 = branch_primitive(16)
P15 = branch_primitive(15)
x = sp.symbols("x")
S = (
    133257064516197275229295 * x**3
    - 3208348664856734112057306 * x**2
    + 349275415597792003907999514 * x
    - 100800589473950881879245888
)
zero(P16 - C4 * T * S.subs(x, C4**2), "quintic row16 branch factor")
P17_at_zero = 49208702469888916532372952734335038991872656250000 * T
zero(P17.subs(C4, 0) - P17_at_zero, "quintic C4 zero contradiction")

# On C4!=0 and S(C4^2)=0, eliminate T between the two linear equations
# P17=P15=0.  Its exact resultant factors into three profiles coprime to S.
Q4 = 148537264577 * x**2 - 4086004213596 * x - 24178143226221
Q8 = (
    22172090448333650297319377595200 * x**4
    - 782476464930749584609729836216645 * x**3
    + 43455557119423145069424296532772242 * x**2
    - 196699263766307709411772113438437964 * x
    + 637608548591074605888432844835874192
)
Q10 = (
    709557651508549958433446867 * x**5
    - 38269012170154817355067751208 * x**4
    + 413899168413250905311964034530 * x**3
    + 3370561945779184898686058483676 * x**2
    - 51355055679882251411436522510600 * x
    + 80149644817765957262829141790560
)
branch_resultant = sp.factor(sp.resultant(P17, P15, T))
expected_branch_resultant = (
    -3363550544979296875 * C4
    * Q4.subs(x, C4**2) * Q8.subs(x, C4**2) * Q10.subs(x, C4**2)
)
zero(branch_resultant - expected_branch_resultant, "quintic T-resultant factorization")

gcd_resultants: list[int] = []
expected_mod_103 = (20, 37, 37)
for polynomial, residue, label in zip((Q4, Q8, Q10), expected_mod_103, ("Q4", "Q8", "Q10")):
    gate(sp.gcd(sp.Poly(S, x), sp.Poly(polynomial, x)).degree() == 0,
         f"gcd(S,{label}) is one")
    certificate = int(sp.resultant(S, polynomial, x))
    gate(certificate % 103 == residue, f"nonzero resultant mod 103 for {label}")
    gcd_resultants.append(certificate)


# Exact reconstruction from any residual root eta with eta*g(eta)!=0.
r_rec = u / e
z_rec = 9 * G * u * (1 + u) / (e * K)
zero(z_rec**2 - K / (9 * G) + Q / (9 * G * e**2 * K**2),
     "Q reconstructs z square")
zero(
    surface.subs({r: r_rec, z: z_rec})
    - u * (1 + u) * Q / (e**3 * K**3),
    "Q reconstructs surface",
)
A_e_symbol = 2 * e + r_rec * D
C_e_rec = sp.factor(9 * G * z_rec**2 - K)
C_z_rec = sp.factor(3 * G * r_rec**2 - 3 * K * A_e_symbol)
zero(C_e_rec + Q / (e**2 * K**2), "Q kills e Hamiltonian")
zero(C_z_rec - 3 * P / e**2, "P kills z Hamiltonian")
C_r_rec = sp.factor(r_rec**2 - 9 * z_rec**2 * A_e_symbol)
zero(
    C_r_rec - P / (G * e**2)
    - Q * (u**2 - P / G) / (e**4 * K**3),
    "P and Q kill r Hamiltonian",
)
zero(Q.subs(u, 0) - e**2, "exclude u=0")
zero(Q.subs(u, -1) + e**2, "exclude u=-1")
zero(Q.subs(u, -sp.Rational(1, 2)) + 729 * G**3 / 16, "exclude K=0")
zero(sp.LC(sp.Poly(Q, u)) + 729 * G**3, "Q excludes infinity")


# Sparse, repeated-root, and dense hostile controls.
for profile in (
    (-1, 5, -10, 10, -5, 1),   # (e-1)^5
    (0, 0, 0, 0, 0, 1),        # e^5
    (0, 0, 1, -2, 1, 1),       # repeated lower roots plus quintic lead
    (2, -3, 5, -7, 11, -13),   # dense hostile control
):
    values = dict(zip((b0, b1, b2, b3, b4, b5), profile))
    evaluated = [sp.expand(R.coeff_monomial(e**degree).subs(values))
                 for degree in top_degrees]
    gate(any(value != 0 for value in evaluated), f"hostile quintic {profile}")

packet_blob = "\n".join(sp.srepr(sp.expand(item)) for item in top_packet).encode()
basis_blob = "\n".join(
    sp.srepr(sp.expand(poly.as_expr())) for poly in I_basis.polys
).encode()
normal_blob = "\n".join(
    f"{degree}:{sp.srepr(sp.expand(normal[degree]))}" for degree in range(21, 14, -1)
).encode()
gcd_blob = "\n".join(str(value) for value in gcd_resultants).encode()
semantic = {
    "carrier": "A=e^2-z/3+r*sum(b_i e^i,0<=i<=5); c=1",
    "new_case": "b5!=0; degree<=4 inherited from THM-3813",
    "resultant": "Res_u(P,Q)=e^4*g^3*H; deg(H)=27; LC=34012224*b5^5",
    "criterion": "boundary-only roots, with repetitions, imply H divides e*g*Hprime",
    "cover": "b5=s^2; C_i=b_i*s^(3-i),0<=i<=4;T=s^7!=0",
    "top_ideal": "F26..F22; 149-element exact Groebner basis",
    "triangular": "rows21..18 solve C3,C2,C1,C0 over C4",
    "branch": "row16=C4*T*S(C4^2); C4=0 dies in row17; S branch dies by Res_T(row17,row15)",
    "gcd_mod_103": expected_mod_103,
    "recovery": "r=u/e; z=9g*u*(u+1)/(e*(1+2u))",
    "open": "degree>=6 g; mixed z2/r corrections; different arm profiles",
}
semantic_blob = json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode()

source = Path(__file__).read_text(encoding="utf-8")
gate(
    not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source))),
    "inactive Python assert",
)

print("theorem=THM-3817-quintic-r-repairs-of-nodal-carriers-have-critical-points")
print("surface=r^2e-z^3+r;c=1;field=algebraically_closed_characteristic_zero")
print("carrier=A=e^2-z/3+r*sum(b_i*e^i,0<=i<=5)")
print("new_case=b5!=0;degree_le_4_inherited=THM-3813")
print("resultant=Res_u(P,Q)=e^4*g^3*H;deg_H=27;LC_H=34012224*b5^5")
print("quotient=degree_5;denominator=512*b5^4")
print("boundary_criterion=all_H_roots_in_V(e*g)_implies_H_divides_e*g*Hprime")
print("cover=b5:s^2;C_i:b_i*s^(3-i);T:s^7_nonzero")
print("top_ideal=F26_to_F22;groebner_basis_size=149")
print("triangular_rows=21_to_18;solve=C3,C2,C1,C0_over_C4")
print("row16=C4*T*S(C4^2);degree_S=3")
print("C4_zero=row17_nonzero_scalar_times_T")
print("S_branch=Res_T(row17,row15)=scalar*C4*Q4*Q8*Q10")
print("gcd_resultant_mod103=20,37,37")
print("recovery=r=u/e;z=9g*u*(u+1)/(e*(1+2u))")
print("open=degree_ge_6_g;mixed_z2_r;different_arm_profile;Darboux_pair")
print(f"invariant_packet_sha256={hashlib.sha256(packet_blob).hexdigest()}")
print(f"groebner_basis_sha256={hashlib.sha256(basis_blob).hexdigest()}")
print(f"normal_packet_sha256={hashlib.sha256(normal_blob).hexdigest()}")
print(f"gcd_certificate_sha256={hashlib.sha256(gcd_blob).hexdigest()}")
print(f"semantic_sha256={hashlib.sha256(semantic_blob).hexdigest()}")
print(f"CHECKS={CHECKS}")
print("RESULT=PASS")
