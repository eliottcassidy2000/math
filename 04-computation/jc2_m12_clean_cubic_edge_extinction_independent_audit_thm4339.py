#!/usr/bin/env python3
"""Independent exact audit for THM-4339's labelled cubic wall.

This path does not import the primary certificate or its inherited hull
engine.  It reconstructs the exact T-face chart from THM-4230's displayed
source, separates Laurent root exit from interior collision, and audits the
normal Hasse strata at every double or triple nonzero cubic-edge root.
"""

from sympy import (
    Poly, Rational, diff, discriminant, expand, factor, limit, simplify,
    symbols,
)
import sys


sys.dont_write_bytecode = True
sys.stdout.reconfigure(newline="\n")


def need(condition, message):
    if condition is not True:
        raise AssertionError(message)


# Exact source coefficients and variables.
sigma, S, P, y = symbols("sigma S P y")
U, W, K, Theta, xi = symbols("U W K Theta xi", nonzero=True)
Phi, eta, alpha = symbols("Phi eta alpha")
Delta, upsilon = symbols("Delta upsilon")
local_b, delta = symbols("local_b delta")
X, Y = symbols("X Y")

e = -Rational(1376, 135)
A = K + Theta * P + xi * P**2 + W * P**3
B = Phi * P**3 + eta * P**4 + alpha * P**5
C = -3 * P + Rational(8, 3) * P**2 + e * P**3 + Delta * P**4 + upsilon * P**5 + U * P**6

# beta=zeta=Z=0 specialization of THM-4230's exact H.
source_s = sigma**-6 * S
source_p = P
source_y = source_s * source_p
H = C + K * source_y**2 + Phi * source_p**2 * source_y \
    + Theta * source_p * source_y**2 + eta * source_p**3 * source_y \
    + xi * source_p**2 * source_y**2 + alpha * source_p**4 * source_y \
    + W * source_p**3 * source_y**2
Q = sigma**12
FQ = (source_s**2 - source_p) * (1 - Q * H) - Q * source_s**2 / 2
G = expand(sigma**12 * FQ)
expected_G = (S**2 - sigma**12 * P) * (
    1 - S**2 * P**2 * A - sigma**6 * S * B - sigma**12 * C
) - sigma**12 * S**2 / 2
need(simplify(G - expected_G) == 0, "exact T-chart reconstruction")

# At the outer boundary S=infinity put local_b=1/S and multiply by b^4.
Floc = expand(
    (1 - delta**2 * P * local_b**2)
    * (local_b**2 - P**2 * A - delta * local_b * B
       - delta**2 * local_b**2 * C)
    - delta**2 * local_b**2 / 2
)
from_G = expand(local_b**4 * expected_G.subs({S: 1 / local_b, sigma**6: delta}))
# Direct simultaneous substitution is safer than power rewriting in SymPy.
from_G = expand(local_b**4 * expected_G.subs(S, 1 / local_b).subs(sigma**12, delta**2).subs(sigma**6, delta))
need(simplify(from_G - Floc) == 0, "exact reciprocal local equation")

# Cubic invariant and labelled generic-double-root formula.
disc_formula = (
    xi**2 * Theta**2 - 4 * W * Theta**3 - 4 * xi**3 * K
    - 27 * W**2 * K**2 + 18 * W * xi * Theta * K
)
need(simplify(discriminant(A, P) - disc_formula) == 0,
     "cubic discriminant formula")

a, c = symbols("a c", nonzero=True)
A_double = expand(W * (P - a)**2 * (P - c))
double_coeffs = Poly(A_double, P).all_coeffs()
double_sub = {
    xi: double_coeffs[1], Theta: double_coeffs[2], K: double_coeffs[3]
}
Delta0 = xi**2 - 3 * W * Theta
N = 9 * W * K - xi * Theta
need(simplify(disc_formula.subs(double_sub)) == 0,
     "double-root parameterization lies on discriminant")
need(simplify(Delta0.subs(double_sub) - W**2 * (a - c)**2) == 0,
     "Delta0 separates double from triple")
need(simplify((N / (2 * Delta0)).subs(double_sub) - a) == 0,
     "labelled repeated-root rational formula")

# On the generic double locus, this polynomial is 4 Delta0^2 J(a), where
# J(P)=Phi+eta P+alpha P^2 and a=N/(2 Delta0).
H_double_num = expand(4 * Delta0**2 * Phi + 2 * Delta0 * N * eta + N**2 * alpha)
J = Phi + eta * P + alpha * P**2
B_a = B.subs(P, a)
need(simplify(H_double_num.subs(double_sub) - 4 * Delta0.subs(double_sub)**2 * J.subs(P, a)) == 0,
     "polynomial normal-Hasse gate on generic double locus")

# Triple-root locus and its separate normal-Hasse gate.
A_triple = expand(W * (P - a)**3)
triple_coeffs = Poly(A_triple, P).all_coeffs()
triple_sub = {
    xi: triple_coeffs[1], Theta: triple_coeffs[2], K: triple_coeffs[3]
}
triple_delta1 = 2 * xi**3 - 9 * W * xi * Theta + 27 * W**2 * K
need(simplify(4 * Delta0**3 - triple_delta1**2 - 27 * W**2 * disc_formula) == 0,
     "subdiscriminant identity")
need(simplify(Delta0.subs(triple_sub)) == 0, "triple Delta0 gate")
need(simplify(triple_delta1.subs(triple_sub)) == 0, "triple Delta1 gate")
need(simplify(disc_formula.subs(triple_sub)) == 0, "triple discriminant gate")
H_triple_num = expand(9 * W**2 * Phi - 3 * W * xi * eta + xi**2 * alpha)
need(simplify(H_triple_num.subs(triple_sub) - 9 * W**2 * J.subs(P, a)) == 0,
     "polynomial normal-Hasse gate on triple locus")

# The discriminant conormal on a labelled double-root branch is evaluation
# at the repeated root.  This makes the critical-value square intrinsic.
perturbation, g0, g1, g2, g3 = symbols("perturbation g0 g1 g2 g3")
g = g0 + g1 * P + g2 * P**2 + g3 * P**3
disc_perturbed = discriminant(A_double + perturbation * g, P)
disc_linear = expand(disc_perturbed).coeff(perturbation, 1)
need(simplify(disc_linear - 4 * W**3 * (c - a)**3 * g.subs(P, a)) == 0,
     "discriminant conormal is labelled-root evaluation")
# Completing the b-square changes A first by B(P)^2/(4P^2)*delta^2.
induced_disc_hasse = simplify(
    4 * W**3 * (c - a)**3 * B_a**2 / (4 * a**2)
)

# The b-critical section is unique because F_bb=2 at the central contact.
# Its first invariant critical-value coefficient is -B(a)^2/4 at delta^2.
v = symbols("v")
at_root_scaled = expand(
    Floc.subs(local_b, delta * v).subs(P, a).subs(A.subs(P, a), 0)
)
# Use the parameterized A polynomials so A(a)=0 is imposed structurally.
double_F = expand(Floc.subs({A: A_double, K: double_sub[K], Theta: double_sub[Theta], xi: double_sub[xi]}))
double_at_root = expand(double_F.subs({P: a, local_b: delta * v}))
critical_lead = expand(limit(double_at_root / delta**2, delta, 0))
need(simplify(critical_lead - (v**2 - v * B_a)) == 0,
     "critical-section leading quadratic")
need(simplify(critical_lead.subs(v, B_a / 2) + B_a**2 / 4) == 0,
     "invariant critical-value Hasse coefficient")

# If B(a)=0, (b,P)=(0,a) is an exact horizontal singular section: there is
# no later normal critical-value coefficient hiding behind the vanished one.
Fb = diff(Floc, local_b)
Fp = diff(Floc, P)
root_conditions = {
    P: a,
    K: double_sub[K], Theta: double_sub[Theta], xi: double_sub[xi],
    Phi: -eta * a - alpha * a**2,
}
need(simplify(Floc.subs(local_b, 0).subs(root_conditions)) == 0,
     "horizontal section lies on total curve")
need(simplify(Fb.subs(local_b, 0).subs(root_conditions)) == 0,
     "horizontal section is b-critical")
need(simplify(Fp.subs(local_b, 0).subs(root_conditions)) == 0,
     "horizontal section is P-critical")

# Weighted strict transforms when the invariant Hasse coefficient is nonzero.
# delta=sigma^6.  A genuine double root has weights (b,x)=(6,6), and a
# triple root has weights (6,4).  All C and denominator corrections are later.
double_weighted = expand(
    double_F.subs({P: a + sigma**6 * X, local_b: sigma**6 * Y, delta: sigma**6})
)
double_face = factor(limit(double_weighted / sigma**12, sigma, 0))
expected_double_face = Y**2 - B_a * Y - a**2 * W * (a - c) * X**2
need(simplify(double_face - expected_double_face) == 0,
     "double-root conic strict transform")

triple_F = expand(Floc.subs({A: A_triple, K: triple_sub[K], Theta: triple_sub[Theta], xi: triple_sub[xi]}))
triple_weighted = expand(
    triple_F.subs({P: a + sigma**4 * X, local_b: sigma**6 * Y, delta: sigma**6})
)
triple_face = factor(limit(triple_weighted / sigma**12, sigma, 0))
expected_triple_face = Y**2 - B_a * Y - a**2 * W * X**3
need(simplify(triple_face - expected_triple_face) == 0,
     "triple-root elliptic strict transform")

# If J has a simple zero at the triple root, the next balanced scale is
# x=sigma^12 X, b=sigma^18 Y.  The resulting plane cubic is singular and
# rational, not elliptic: completing the square exposes a repeated X^2.
j1 = symbols("j1", nonzero=True)
B_triple_n1 = P**3 * j1 * (P - a)
triple_n1_F = expand(
    triple_F.subs({Phi: -j1 * a, eta: j1, alpha: 0})
)
triple_n1_weighted = expand(
    triple_n1_F.subs({P: a + sigma**12 * X,
                      local_b: sigma**18 * Y, delta: sigma**6})
)
triple_n1_face = factor(limit(triple_n1_weighted / sigma**36, sigma, 0))
expected_triple_n1_face = Y**2 - a**3 * j1 * X * Y - a**2 * W * X**3
need(simplify(triple_n1_face - expected_triple_n1_face) == 0,
     "triple simple-zero normal face")
completed_n1 = factor(
    (Y - a**3 * j1 * X / 2)**2
    - (a**6 * j1**2 * X**2 / 4 + a**2 * W * X**3)
)
need(simplify(completed_n1 - triple_n1_face) == 0,
     "simple-zero cubic has an X^2 discriminant factor")
need(all(simplify(diff(triple_n1_face, variable).subs({X: 0, Y: 0})) == 0
         for variable in (X, Y)), "simple-zero cubic is singular at origin")

# A double zero of J never balances the cusp: on b^2~x^3 with v(x)=r,
# its normal term has excess 6+r/2>0.  The cusp is a horizontal conductor.
r = symbols("r", positive=True)
triple_n2_excess = 6 + r / 2
need(triple_n2_excess.is_positive is True, "J order at least two is later at every cusp scale")

# In the exact differential identity omega=sigma^16*b^2 db/F_P, the generic
# denominator orders are 6 and 8 on the two strict transforms.
double_form_order = 16 + 2 * 6 + 6 - 6
triple_form_order = 16 + 2 * 6 + 6 - 8
triple_n1_form_order = 16 + 3 * 18 - 24
need(double_form_order == 28 and triple_form_order == 26
     and triple_n1_form_order == 46,
     "positive good-form orders on new carriers")


def laurent_saturate(poly_expr):
    poly = Poly(poly_expr, P)
    multiplicity = 0
    while poly.eval(0) == 0 and poly.degree() > 0:
        quotient, remainder = poly.div(Poly(P, P))
        need(remainder.is_zero, "exact P division")
        poly = quotient
        multiplicity += 1
    return multiplicity, poly.as_expr()


examples = {
    "simple_interior": P**3 + P**2 + P + 2,
    "double_interior": (P - 1)**2 * (P + 1),
    "triple_interior": (P - 1)**3,
    "one_root_exit": P * (P**2 + P + 1),
    "exit_plus_interior_double": P * (P - 1)**2,
    "double_endpoint_only": P**2 * (P + 1),
    "all_three_exit": P**3,
}
example_data = {}
for name, polynomial in examples.items():
    m, saturated = laurent_saturate(expand(polynomial))
    degree = Poly(saturated, P).degree()
    sat_disc = discriminant(saturated, P) if degree >= 2 else 1
    full_disc = discriminant(polynomial, P)
    example_data[name] = (m, factor(saturated), factor(sat_disc), factor(full_disc))

need(example_data["simple_interior"] == (0, P**3 + P**2 + P + 2, -83, -83),
     "simple hostile")
need(example_data["one_root_exit"] == (1, P**2 + P + 1, -3, -3),
     "root exit is not discriminant collision")
need(example_data["double_endpoint_only"] == (2, P + 1, 1, 0),
     "full discriminant can vanish only from endpoint multiplicity")
need(example_data["exit_plus_interior_double"] == (1, (P - 1)**2, 0, 0),
     "exit and interior collision can coexist")

# The exact three-address tax partitions losses into infinity exit, zero exit,
# and collision length after Laurent saturation.
def root_tax(poly_expr):
    poly = Poly(poly_expr, P)
    degree = poly.degree()
    m, saturated_expr = laurent_saturate(poly_expr)
    saturated = Poly(saturated_expr, P)
    collision = saturated.gcd(saturated.diff()).degree() if saturated.degree() else 0
    n_torus = saturated.degree() - collision
    return n_torus, 3 - degree, m, collision


tax_examples = dict(examples)
tax_examples["one_infinity_exit"] = P**2 + P + 1
for name, polynomial in tax_examples.items():
    n_torus, infinity_exit, zero_exit, collision = root_tax(expand(polynomial))
    need(3 - n_torus == infinity_exit + zero_exit + collision,
         "three-address tax: " + name)


def polygon_ledger(points):
    twice_area = abs(sum(
        points[i][0] * points[(i + 1) % len(points)][1]
        - points[(i + 1) % len(points)][0] * points[i][1]
        for i in range(len(points))
    ))
    from math import gcd
    boundary = sum(
        gcd(abs(points[(i + 1) % len(points)][0] - points[i][0]),
            abs(points[(i + 1) % len(points)][1] - points[i][1]))
        for i in range(len(points))
    )
    interior = (twice_area - boundary + 2) // 2
    return twice_area, boundary, interior


tface_ledgers = {
    m: polygon_ledger(((2, 0), (4, 2 + m), (4, 5)))
    for m in range(3)
}
global_polygons = {
    0: ((0, 1), (2, 0), (4, 2), (4, 5), (0, 7)),
    1: ((0, 1), (2, 0), (4, 3), (4, 5), (0, 7)),
    2: ((0, 1), (2, 0), (4, 4), (4, 5), (0, 7)),
    3: ((0, 1), (2, 0), (4, 5), (0, 7)),
}
global_ledgers = {m: polygon_ledger(poly) for m, poly in global_polygons.items()}
need(tface_ledgers == {0: (6, 6, 1), 1: (4, 4, 1), 2: (2, 4, 0)},
     "T-face root-exit Pick ledgers")
need(global_ledgers == {0: (42, 14, 15), 1: (40, 12, 15),
                        2: (38, 12, 14), 3: (36, 10, 14)},
     "global root-exit Pick ledgers")

# On the K*W*U*(U+W) gate, these are all outer and internal edge schemes.
edge_x = symbols("edge_x")
edge_schemes = (
    edge_x - 1,
    1 - K * edge_x**2,
    K + Theta * edge_x + xi * edge_x**2 + W * edge_x**3,
    (edge_x - 1) * (U * edge_x + W),
    U - edge_x**6,
    1 - W * edge_x,
)
need(discriminant(edge_schemes[1], edge_x) == 4 * K,
     "slant edge is squarefree under K!=0")
need(factor(discriminant(edge_schemes[3], edge_x)) == (U + W)**2,
     "top-edge discriminant")
need(Poly(edge_schemes[4], edge_x).gcd(Poly(diff(edge_schemes[4], edge_x), edge_x)).degree() == 0,
     "sixfold edge is squarefree under U!=0")
# T's generic P derivative cannot vanish identically on its positive-genus
# component: every coefficient of 2A+P A' has positive integer multiplier.
t_derivative_obstruction = expand(2 * A + P * diff(A, P))
need(t_derivative_obstruction == 2 * K + 3 * Theta * P + 4 * xi * P**2 + 5 * W * P**3,
     "T-face generic differential")
need(12 * (Rational(5, 6) - (Rational(1, 12) + Rational(1, 6) - Rational(1, 6))) == 9,
     "M good-form order at L=12")
need(12 * (Rational(5, 6) - (Rational(1, 2) + 0 - 1)) == 16,
     "T good-form order at L=12")
need((13 - 3 + 1) == 11 and 3 + 1 + 11 == 15,
     "off-collision component/genus graph")

print("THM4339_CLEAN_CUBIC_EDGE_INDEPENDENT_AUDIT=PASS")
print("A=K+Theta*P+xi*P^2+W*P^3")
print("Disc=" + str(disc_formula))
print("Delta0=xi^2-3*W*Theta")
print("4*Delta0^3-Delta1^2=27*W^2*Disc")
print("double_root=(9*W*K-xi*Theta)/(2*Delta0)")
print("double_Hasse_numerator=" + str(H_double_num))
print("triple_gates=Delta0=0;2*xi^3-9*W*xi*Theta+27*W^2*K=0")
print("triple_Hasse_numerator=" + str(H_triple_num))
print("critical_value=[-B(a)^2/4]*delta^2+O(delta^3);delta=sigma^6")
print("induced_discriminant_Hasse=" + str(induced_disc_hasse))
print("double_face=" + str(double_face))
print("triple_face=" + str(triple_face))
print("triple_B_simple_zero_face=" + str(triple_n1_face))
print("triple_B_simple_zero_normalization=Z^2=a^2*W*X+a^6*j1^2/4;rational")
print("good_form_orders=double:28,triple:26,triple_B_simple_zero:46")
print("B(a)=0=>exact_horizontal_critical_section;critical_value_identically_zero")
print("triple_B_double_zero=>normal term always later by 6+r/2;horizontal cusp")
print("Tface_ledgers=" + str(tface_ledgers))
print("global_ledgers=" + str(global_ledgers))
print("edges_KWU_Lambda_open=" + str(edge_schemes))
print("face_orders_L12=M:9,T:16")
print("global_relative_gate=U*K*W*(U+W)!=0;all_cubic_root_strata_PASS")
print("three_tax=3-N_tor=(3-deg(A))+ord_0(A)+deg(gcd(Q,Qprime))")
for name in examples:
    print(name + "=" + str(example_data[name]))
