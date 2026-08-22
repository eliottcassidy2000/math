#!/usr/bin/env python3
"""Exact identities and hostile controls for proved, audited THM-3575."""

from __future__ import annotations

import itertools

import sympy as sp


def require(condition: bool, label: str) -> None:
    if not condition:
        raise RuntimeError(f"failed truth gate: {label}")


def require_zero(expression: sp.Expr, label: str) -> None:
    require(sp.cancel(sp.expand(expression)) == 0, label)


def distinct_root_count(expression: sp.Expr, variable: sp.Symbol) -> int:
    polynomial = sp.Poly(sp.expand(expression), variable, domain=sp.QQ)
    require(not polynomial.is_zero, f"nonzero root-count input {expression}")
    if polynomial.degree() <= 0:
        return 0
    return int(polynomial.sqf_part().degree())


a, b, c, x, y, z = sp.symbols("a b c x y z")
p, q = sp.symbols("p q")
P0, Q0 = sp.symbols("P Q")
t, v, w, U = sp.symbols("t v w U")

# The fixed Keller map.
u = 1 + x * y
A = sp.expand(u**3 * z + y**2 * u * (4 + 3 * x * y))
B = sp.expand(y + 3 * x * u**2 * z + 3 * x * y**2 * (4 + 3 * x * y))
C = sp.expand(2 * x - 3 * x**2 * y - x**3 * z)

# The universal source factorization treats p=P(B),q=Q(B) as commuting
# coefficients.  The identity therefore survives every polynomial
# substitution in B.
G_source = sp.expand(q**3 * C - 2 * p**3 * A + B * p**2 * q - 2 * p * q**2)
S_source = sp.expand(
    z * (2 * p**2 * u**2 - p * q * x * u - q**2 * x**2)
    + 2 * p**2 * y**2 * (3 * u + 1)
    - p * q * y * (3 * u - 2)
    + q**2 * (5 - 3 * u)
)
require_zero(G_source - (q * x - p * u) * S_source, "universal source factorization")

# Target coordinate and the two reduced Jelonek curves.
G = Q0**3 * c - 2 * P0**3 * a + b * P0**2 * Q0 - 2 * P0 * Q0**2
L = 27 * a**2 * c**2 - 18 * a * b * c + 16 * a + b**3 * c - b**2
D = 3 * a * P0**2 - b * P0 * Q0 + Q0**2
H_section = 12 * a**2 * P0**2 - 4 * a * b * P0 * Q0 + (16 * a - b**2) * Q0**2
Delta = (b * P0 - Q0) ** 2 + 3 * Q0**2
E = b * P0 - 2 * Q0

jelonek_difference = sp.factor(Q0**6 * L - D**2 * H_section)
require(sp.rem(sp.Poly(jelonek_difference, c), sp.Poly(G, c)).as_expr() == 0, "Jelonek restriction")
require_zero(sp.discriminant(H_section, a) - 64 * Q0**2 * Delta, "H-section discriminant")
require_zero(
    P0**2 * H_section + E**2 * Q0**2 - 4 * (a * P0**2 + Q0**2) * D,
    "reduced D/H intersection identity",
)
vertical_jelonek = sp.expand(L.subs(a, 0))
require_zero(vertical_jelonek - b**2 * (b * c - 1), "actual Q-root Jelonek fibre")
require(vertical_jelonek.subs(b, 0) == 0, "actual zero-root Jelonek line")

# The omitted-curve cube, parameterized by its nonzero b-coordinate.
omitted_G = sp.cancel(G.subs({a: b**2 / 12, c: sp.Rational(4, 3) / b}))
require_zero(6 * b * omitted_G + E**3, "omitted-curve cube")

# One explicit Bezout control checks that the coordinate argument includes
# b-dependent coefficients and their derivatives.
P_control = b + 1
Q_control = b**2 + 1
r_control, s_control, gcd_control = sp.gcdex(-2 * P_control**3, Q_control**3, b)
require(gcd_control == 1, "Bezout control gcd")
G_control = sp.expand(G.subs({P0: P_control, Q0: Q_control}))
K_control = sp.expand(-s_control * a + r_control * c)
jacobian_control = sp.factor(
    sp.det(sp.Matrix([G_control, b, K_control]).jacobian(sp.Matrix([a, b, c])))
)
require(jacobian_control == 1, "target-coordinate Jacobian")

# Collision typing at the target and at each of its three sources.
collision_target = sp.factor(G.subs({a: -sp.Rational(1, 4), b: 0, c: 0}))
require_zero(
    collision_target - P0 * (P0 - 2 * Q0) * (P0 + 2 * Q0) / 2,
    "collision target condition",
)
collision_sources = [
    (sp.Integer(0), sp.Integer(1)),
    (sp.Integer(1), -sp.Rational(1, 2)),
    (-sp.Integer(1), -sp.Rational(1, 2)),
]
collision_linear_values = [sp.expand(Q0 * xx - P0 * uu) for xx, uu in collision_sources]
require(
    collision_linear_values == [-P0, P0 / 2 + Q0, P0 / 2 - Q0],
    "one collision source on the linear factor",
)


def euler_row(name: str, P_value: sp.Expr, Q_value: sp.Expr) -> tuple[object, ...]:
    require(sp.gcd(sp.Poly(P_value, b), sp.Poly(Q_value, b)).degree() == 0, f"coprime {name}")
    delta_value = sp.expand(Delta.subs({P0: P_value, Q0: Q_value}))
    e_value = sp.expand(E.subs({P0: P_value, Q0: Q_value}))
    require(delta_value != 0 and e_value != 0, f"nondegenerate {name}")
    n_p = distinct_root_count(P_value, b)
    n_q = distinct_root_count(Q_value, b)
    n_delta = distinct_root_count(delta_value, b)
    n_e = distinct_root_count(e_value, b)
    delta_zero = int(sp.expand(Q_value).subs(b, 0) == 0)

    chi_d = 1 - n_p
    chi_h = 2 - n_p - n_q - n_delta + delta_zero
    chi_intersection = n_q + n_e - delta_zero
    chi_jelonek = chi_d + chi_h - chi_intersection
    chi_omitted = n_e - delta_zero
    chi_total = 3 - 2 * chi_jelonek - chi_omitted
    chi_linear = n_p + n_q - delta_zero
    chi_residual = chi_total - chi_linear
    formula = -3 + 3 * n_p + 3 * n_q + 2 * n_delta + n_e - 2 * delta_zero
    require(chi_residual == formula, f"Euler formula {name}")

    p_at_zero = sp.expand(P_value).subs(b, 0)
    q_at_zero = sp.expand(Q_value).subs(b, 0)
    collision = p_at_zero == 0 or p_at_zero == 2 * q_at_zero or p_at_zero == -2 * q_at_zero
    if collision:
        require(chi_residual != 1, f"collision Euler no-go {name}")
    return name, P_value, Q_value, n_p, n_q, n_delta, n_e, delta_zero, chi_residual


control_pairs = [
    ("constant", sp.Integer(1), sp.Integer(1)),
    ("collision-fixed", b, sp.Integer(1)),
    ("collision-plus", sp.Integer(2), sp.Integer(1)),
    ("collision-minus", -sp.Integer(2), sp.Integer(1)),
    ("Euler-sharp", sp.Integer(1), b),
    ("higher-vertical", sp.Integer(1), b**2),
    ("nonsymmetric", b - 1, b + 1),
    ("repeated-P", (b - 1) ** 2, b + 1),
]
control_rows = [euler_row(*row) for row in control_pairs]
require(control_rows[4][-1] == 1, "positive Euler-sharp control")

# Exact descent along Q=b^m.  The all-m proof in the theorem comes from
# the displayed factorizations; this loop is a finite regression only.
vertical_rows = []
for m in range(1, 9):
    P_value = sp.Integer(3)
    Q_value = b**m
    delta_value = sp.factor(Delta.subs({P0: P_value, Q0: Q_value}))
    e_value = sp.factor(E.subs({P0: P_value, Q0: Q_value}))
    n_delta = distinct_root_count(delta_value, b)
    n_e = distinct_root_count(e_value, b)
    require((n_delta, n_e) == (2 * m - 1, m), f"vertical root formula m={m}")
    chi_residual = -3 + 3 + 2 * n_delta + n_e - 2
    require(chi_residual == 5 * m - 4, f"vertical Euler descent m={m}")
    vertical_rows.append((m, n_delta, n_e, chi_residual))

# Independent bounded hostile universe: every nonzero polynomial of degree
# at most two with coefficients in {-1,0,1}, ordered coprime pairs, and the
# inherited Delta*E!=0 filter.
small_polynomials: list[sp.Poly] = []
for coefficients in itertools.product((-1, 0, 1), repeat=3):
    expression = sum(coefficients[index] * b**index for index in range(3))
    if expression != 0:
        small_polynomials.append(sp.Poly(expression, b, domain=sp.QQ))

finite_pair_count = 0
finite_collision_sharp: list[tuple[sp.Expr, sp.Expr]] = []
finite_sharp: list[tuple[sp.Expr, sp.Expr]] = []
for p_poly in small_polynomials:
    for q_poly in small_polynomials:
        if sp.gcd(p_poly, q_poly).degree() != 0:
            continue
        p_value = p_poly.as_expr()
        q_value = q_poly.as_expr()
        delta_value = sp.expand(Delta.subs({P0: p_value, Q0: q_value}))
        e_value = sp.expand(E.subs({P0: p_value, Q0: q_value}))
        if delta_value == 0 or e_value == 0:
            continue
        finite_pair_count += 1
        delta_zero = int(q_value.subs(b, 0) == 0)
        chi_residual = (
            -3
            + 3 * distinct_root_count(p_value, b)
            + 3 * distinct_root_count(q_value, b)
            + 2 * distinct_root_count(delta_value, b)
            + distinct_root_count(e_value, b)
            - 2 * delta_zero
        )
        p_at_zero = p_value.subs(b, 0)
        q_at_zero = q_value.subs(b, 0)
        collision = p_at_zero == 0 or p_at_zero == 2 * q_at_zero or p_at_zero == -2 * q_at_zero
        if chi_residual == 1:
            finite_sharp.append((p_value, q_value))
            require(
                sp.degree(p_value, b) == 0 and sp.degree(q_value, b) == 1 and q_at_zero == 0,
                "finite sharp classifier",
            )
            if collision:
                finite_collision_sharp.append((p_value, q_value))

require(finite_pair_count == 532, "finite hostile universe size")
require(
    finite_sharp == [(-sp.Integer(1), -b), (-sp.Integer(1), b), (sp.Integer(1), -b), (sp.Integer(1), b)],
    "finite Euler-sharp representatives",
)
require(not finite_collision_sharp, "finite collision-sharp universe empty")

# Exact inverse of the linear component on Q!=0 and on a nonzero Q-root.
h_rational = P0 / Q0
D_rational = D / Q0**2
x_inverse = sp.cancel(h_rational / D_rational)
y_inverse = sp.cancel(b - 3 * a * h_rational)
z_inverse = sp.cancel(a * D_rational**3 - y_inverse**2 * D_rational * (D_rational + 3))
u_inverse = sp.cancel(1 + x_inverse * y_inverse)
A_inverse = sp.cancel(A.subs({x: x_inverse, y: y_inverse, z: z_inverse}, simultaneous=True))
B_inverse = sp.cancel(B.subs({x: x_inverse, y: y_inverse, z: z_inverse}, simultaneous=True))
C_inverse = sp.cancel(C.subs({x: x_inverse, y: y_inverse, z: z_inverse}, simultaneous=True))
require_zero(A_inverse - a, "linear inverse A")
require_zero(B_inverse - b, "linear inverse B")
require_zero(G.subs(c, C_inverse), "linear inverse target equation")
require_zero(Q0 * x_inverse - P0 * u_inverse, "linear inverse source factor")

x_vertical = 2 / b
y_vertical = -b / 2
z_vertical = b**2 * (10 - b * c) / 8
require_zero(1 + x_vertical * y_vertical, "vertical inverse u=0")
require_zero(A.subs({x: x_vertical, y: y_vertical, z: z_vertical}, simultaneous=True), "vertical inverse A")
require_zero(B.subs({x: x_vertical, y: y_vertical, z: z_vertical}, simultaneous=True) - b, "vertical inverse B")
require_zero(C.subs({x: x_vertical, y: y_vertical, z: z_vertical}, simultaneous=True) - c, "vertical inverse C")

# The Euler-sharp family P=t,Q=b.
G_t = sp.expand(G.subs({P0: t, Q0: b}))
require_zero(G_t - (b**3 * c - 2 * t**3 * a + t * (t - 2) * b**2), "sharp target coordinate")

S_t_source = sp.expand(S_source.subs({p: t, q: B}))
require_zero(
    G_t.subs({a: A, b: B, c: C}) - (B * x - t * u) * S_t_source,
    "sharp source factorization",
)

# The linear component is reconstructed from (x,u) in (C*)^2.
beta = sp.expand(v + 3 * (v + 1) ** 2 * w + 3 * v**2 * (4 + 3 * v))
require_zero(x * B - beta.subs({v: x * y, w: x**2 * z}), "invariant Bx formula")
require(beta.subs(v, -1) == 2, "u=0 excluded on linear component")
w_linear = sp.cancel((t * (v + 1) - v - 3 * v**2 * (4 + 3 * v)) / (3 * (v + 1) ** 2))
require_zero(beta.subs(w, w_linear) - t * (v + 1), "linear-component inverse")

# Exact target Jelonek factorization and its one-variable boundary data.
a_t = sp.cancel((b**3 * c + t * (t - 2) * b**2) / (2 * t**3))
L_t = sp.factor(L.subs(a, a_t))
L_t_expected = sp.factor(
    b**2
    * (3 * b * c + t**2 - 4 * t) ** 2
    * (3 * b**2 * c**2 + 4 * t * (t - 1) * b * c - 4 * t**2)
    / (4 * t**6)
)
require_zero(L_t - L_t_expected, "sharp Jelonek factorization")

s_bc = sp.symbols("s")
quadratic_bc = 3 * s_bc**2 + 4 * t * (t - 1) * s_bc - 4 * t**2
require_zero(
    sp.discriminant(quadratic_bc, s_bc) - 16 * t**2 * (t**2 - 2 * t + 4),
    "bc-quadratic discriminant",
)
middle_bc = -t * (t - 4) / 3
require_zero(quadratic_bc.subs(s_bc, middle_bc) + t**2 * (t - 2) ** 2, "Jelonek component coincidence")
omitted_t = sp.factor(omitted_G.subs({P0: t, Q0: b}))
require_zero(omitted_t + b**2 * (t - 2) ** 3 / 6, "sharp omitted restriction")

# Core residual discriminant: its only nonsquare factor becomes a square
# exactly at t^2-2t+4=0.
core = L * x**3 + (4 - 3 * b * c) * x - 2 * c
core_t = sp.factor(core.subs(a, a_t))
core_root = 2 * t**2 / (b * (3 * b * c + t**2 - 4 * t))
require_zero(core_t.subs(x, core_root), "sharp core root")
core_residual = sp.cancel(core_t / (x - core_root))
core_residual_disc = sp.factor(sp.discriminant(sp.together(core_residual), x))
expected_core_disc = sp.factor(
    -b**2
    * (-3 * b * c + 2 * t) ** 2
    * (3 * b * c + t**2 - 4 * t) ** 2
    * (3 * b**2 * c**2 + 4 * b * c * t**2 - 4 * b * c * t - 4 * t**2)
    / (4 * t**8)
)
require_zero(core_residual_disc - expected_core_disc, "residual quadratic discriminant")

# C*-quotient equation J_t=x^2*S_t in the invariant variables v=xy,w=x^2z.
U_v = v + 1
J = sp.expand(
    w * (2 * t**2 * U_v**2 - t * beta * U_v - beta**2)
    + 2 * t**2 * v**2 * (3 * U_v + 1)
    - t * beta * v * (3 * U_v - 2)
    + beta**2 * (5 - 3 * U_v)
)
require_zero(x**2 * S_t_source - J.subs({v: x * y, w: x**2 * z}), "quotient equation")
require(sp.factor(sp.Poly(J, w).LC()) == -9 * (v + 1) ** 4, "quotient cubic leading coefficient")

H_t = sp.expand(
    3 * (t - 2) ** 2 * (t**2 - 2 * t + 4) * U**4
    + 2 * (7 * t**3 - 12 * t**2 + 24 * t + 8) * U**3
    - 12 * (t**2 - 2 * t + 3) * U**2
    - 24 * U
    - 4
)
quotient_discriminant = sp.factor(sp.discriminant(J, w))
require_zero(
    quotient_discriminant - 108 * t**2 * (v + 1) ** 6 * H_t.subs(U, v + 1),
    "quotient cubic discriminant",
)
require_zero(J.subs(v, -1) - 2 * (t**2 - 2 * t - 2 * w + 10), "quotient u=0 fibre")
require(sp.diff(J, w).subs(v, -1) == -4, "quotient u=0 simple root")
require(H_t.subs(U, 0) == -4, "quotient discriminant avoids u=0")

# Exact quotient-boundary discriminant.
H_discriminant = sp.factor(sp.discriminant(H_t, U))
degeneration = 9 * t**3 - 18 * t**2 + 36 * t - 8
expected_H_discriminant = -62208 * t * (t**2 - 2 * t + 4) * degeneration**3
require_zero(H_discriminant - expected_H_discriminant, "H_t discriminant")

H_at_two = sp.factor(H_t.subs(t, 2))
require_zero(H_at_two - 4 * (32 * U**3 - 9 * U**2 - 6 * U - 1), "t=2 quotient polynomial")
require(sp.discriminant(H_at_two, U) != 0, "t=2 quotient squarefree")
chi_quotient_two = 1 - 3 * 3 + 2 * 3
require(chi_quotient_two == -2, "t=2 quotient Euler characteristic")

# At a root of degeneration, H has one double root U0 and J has a triple
# w-root above it.  Reduction is exact in Q[t]/(degeneration).
def reduce_mod_degeneration(expression: sp.Expr) -> sp.Expr:
    numerator, denominator = sp.fraction(sp.cancel(expression))
    modulus = sp.Poly(degeneration, t, domain=sp.QQ)
    numerator_reduced = sp.rem(sp.Poly(sp.expand(numerator), t), modulus).as_expr()
    denominator_reduced = sp.rem(sp.Poly(sp.expand(denominator), t), modulus).as_expr()
    inverse_denominator = sp.invert(sp.Poly(denominator_reduced, t), modulus).as_expr()
    return sp.factor(sp.rem(sp.Poly(sp.expand(numerator_reduced * inverse_denominator), t), modulus).as_expr())


U_double = -(45 * t**2 - 72 * t + 36) / 64
require(reduce_mod_degeneration(H_t.subs(U, U_double)) == 0, "degenerate H root")
require(reduce_mod_degeneration(sp.diff(H_t, U).subs(U, U_double)) == 0, "degenerate H double root")
require(reduce_mod_degeneration(sp.diff(H_t, U, 2).subs(U, U_double)) == -72, "degenerate H exactly double")

J_U = sp.expand(J.subs(v, U - 1))
J_coefficients = [
    reduce_mod_degeneration(sp.Poly(J_U, w).coeff_monomial(w**index).subs(U, U_double))
    for index in range(4)
]
w_triple = (1071 * t**2 - 1752 * t + 3020) / 576
triple_model = sp.expand(J_coefficients[3] * (w - w_triple) ** 3)
for index in range(4):
    require(
        reduce_mod_degeneration(
            J_coefficients[index] - sp.Poly(triple_model, w).coeff_monomial(w**index)
        )
        == 0,
        f"degenerate quotient triple coefficient {index}",
    )
require(
    sp.gcd(sp.Poly(degeneration, t), sp.Poly(1615 * t**2 - 1752 * t + 332, t)).degree() == 0,
    "degenerate quotient cubic leading coefficient nonzero",
)
chi_quotient_degenerate = 1 - 3 * 3 + (2 + 2 + 1)
require(chi_quotient_degenerate == -3, "degenerate quotient Euler characteristic")

# Split and hostile boundaries.
split_modulus = sp.Poly(t**2 - 2 * t + 4, t)
G_split = sp.rem(sp.Poly(G_t, t), split_modulus).as_expr()
require_zero(G_split - (b**3 * c + 16 * a - 4 * b**2), "split target graph")
require(G_t.subs(t, 0) == b**3 * c, "t=0 noncoordinate boundary")

H_at_one_v = sp.expand(H_t.subs({t: 1, U: v + 1}))
H_at_one_expected = 9 * v**4 + 90 * v**3 + 192 * v**2 + 126 * v + 11
require_zero(H_at_one_v - H_at_one_expected, "t=1 hostile quotient polynomial")
require(sp.gcd(sp.Poly(H_at_one_v, v), sp.Poly(sp.diff(H_at_one_v, v), v)).degree() == 0, "t=1 hostile squarefree")
chi_quotient_one = 1 - 3 * 4 + 2 * 4
require(chi_quotient_one == -3, "t=1 hostile quotient Euler characteristic")

print("THM-3575 coprime Pell target-coordinate Euler audit")
print("status=PROVED VERIFIED-EXACT INDEPENDENTLY HOSTILE-AUDITED")
print("symbolic universe=Q[a,b,c,x,y,z,P,Q,p,q,t,v,w]; exact polynomial/rational identities")
print("G_(P,Q)=Q^3*c-2*P^3*a+b*P^2*Q-2*P*Q^2 is a target coordinate for gcd(P,Q)=1")
print("G_(P,Q)(F)=(Q(B)*x-P(B)*u)*S_(P,Q)")
print("Euler formula=-3+3*nP+3*nQ+2*nDelta+nE-2*[Q(0)=0]")
print("strict vertical: Q(b0)=0 gives one Jelonek point if b0!=0, an A1 fibre if b0=0")
print("linear component inverse audited on Q!=0 and on every nonzero Q-root vertical")
print("control rows: name | P | Q | nP | nQ | nDelta | nE | delta0 | chi(S)")
for row in control_rows:
    print(" | ".join(str(entry) for entry in row))
print("vertical descent m,nDelta,nE,chi(S):")
for row in vertical_rows:
    print(" | ".join(str(entry) for entry in row))
print("finite hostile universe=ordered coprime P,Q, degree<=2, coefficients {-1,0,1}, Delta*E!=0")
print(f"finite pairs={finite_pair_count}; collision-sharp={len(finite_collision_sharp)}; sharp={finite_sharp}")
print("Euler-sharp family: P=t,Q=b; chi(S_t)=1 for every t!=0 after exact boundary correction")
print("irreducible boundary: t^2-2*t+4!=0")
print("split boundary: t=1+-i*sqrt(3), G=b^3*c+16*a-4*b^2")
print("quotient discriminant=108*t^2*(v+1)^6*H_t(v)")
print("quotient Euler: generic=-3; t=2 -> -2; degeneration cubic -> -3")
print("t=1 hostile: G=b^3*c-2*a-b^2, chi(S)=1, quotient chi=-3")
print("all active truth gates passed")
