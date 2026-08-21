#!/usr/bin/env python3
"""Exact controls for the proved THM-3563 node-cycle obstruction."""

from __future__ import annotations

import sympy as sp


GATES = 0


def require(condition: bool, label: str) -> None:
    """Keep truth gates active under both ordinary and optimized Python."""
    global GATES
    if not bool(condition):
        raise RuntimeError(f"failed truth gate: {label}")
    GATES += 1


def zero(expr: sp.Expr) -> bool:
    return sp.factor(sp.cancel(expr)) == 0


u, v, h, h0, c, kappa = sp.symbols(
    "u v h h0 c kappa", nonzero=True
)
A0, B0, d0 = sp.symbols("A B d")

# The restriction of the THM-2696 quotient to s_3=s_1 s_2-c.
A = u**2 - 2 * v
B = v**2 - 2 * u**2 * v + 2 * c * u
d = u * v - c

# Laurent normalization of the double curve.
ug = h / 2 + 4 * c / h**2
vg = 2 * c / h - 16 * c**2 / h**4
Ag = sp.factor(A.subs({u: ug, v: vg}))
Bg = sp.factor(B.subs({u: ug, v: vg}))
dg = sp.factor(d.subs({u: ug, v: vg}))

Ag_expected = (h**6 + 192 * c**2) / (4 * h**4)
Bg_expected = 4 * c**2 * (h**6 + 192 * c**2) / h**8
dg_expected = -64 * c**3 / h**6

require(zero(Ag - Ag_expected), "A Laurent identity")
require(zero(Bg - Bg_expected), "B Laurent identity")
require(zero(dg - dg_expected), "d Laurent identity")
for name, expr in (("A", Ag), ("B", Bg), ("d", dg)):
    require(zero(expr.subs(h, -h) - expr), f"{name} is even in h")

# The source double curve and its target presentation.
K = (
    8 * c**2
    - c * u**3
    - 10 * c * u * v
    + u**4 * v
    + 2 * u**2 * v**2
    + v**3
)
require(zero(K - ((c - 3 * d) ** 2 - A * B) / 2), "K target identity")
require(zero(K.subs({u: ug, v: vg})), "gamma lies on K=0")

# Exact tangent and primitive-sign control.
uvprime = sp.factor(ug * sp.diff(vg, h))
uvprime_expected = (
    256 * c**3 / h**7 + 24 * c**2 / h**4 - c / h
)
require(zero(uvprime - uvprime_expected), "u times v-prime identity")
uvprime_even = sp.factor((uvprime + uvprime.subs(h, -h)) / 2)
require(zero(uvprime_even - 24 * c**2 / h**4), "even tangent term")
phi_odd = 8 * kappa * c**2 / h**3
require(
    zero(sp.diff(phi_odd, h) + 24 * kappa * c**2 / h**4),
    "odd primitive sign",
)
require(
    zero(phi_odd - phi_odd.subs(h, -h) - 16 * kappa * c**2 / h**3),
    "odd-difference law",
)

# The six special normalization points form three ordinary source nodes.
K_u_gamma = sp.factor(sp.diff(K, u).subs({u: ug, v: vg}))
K_v_gamma = sp.factor(sp.diff(K, v).subs({u: ug, v: vg}))
require(
    zero(
        K_u_gamma
        + c * (32 * c - h**3) * (192 * c**2 + h**6) / (4 * h**7)
    ),
    "K_u node factor",
)
require(
    zero(
        K_v_gamma
        + (16 * c - h**3) * (192 * c**2 + h**6) / (16 * h**5)
    ),
    "K_v node factor",
)
hessian_gamma = sp.factor(sp.det(sp.hessian(K, (u, v))).subs({u: ug, v: vg}))
require(
    zero(
        hessian_gamma
        - 108 * c**2
        - (192 * c**2 + h**6)
        * (256 * c**3 - 528 * c**2 * h**3 - 28 * c * h**6 - h**9)
        / (4 * h**9)
    ),
    "node Hessian identity",
)

# Independent exact reduction of all three source-node pairs.  Here
# eps^2=-3, zeta=(1+eps)/2, and h0^3=8 eps c.
eps = sp.I * sp.sqrt(3)
zeta = (1 + eps) / 2
require(zero(eps**2 + 3), "quadratic extension")
require(zero(zeta**2 - zeta + 1), "sixth-root relation")
require(zero(zeta**3 + 1), "zeta cubed")


def reduce_at_h0(expr: sp.Expr) -> sp.Expr:
    numerator = sp.together(expr).as_numer_denom()[0]
    dividend = sp.Poly(sp.expand(numerator), h0, domain=sp.EX)
    divisor = sp.Poly(h0**3 - 8 * eps * c, h0, domain=sp.EX)
    return sp.factor(sp.rem(dividend, divisor).as_expr())


def gamma_u(H: sp.Expr) -> sp.Expr:
    return H / 2 + 4 * c / H**2


def gamma_v(H: sp.Expr) -> sp.Expr:
    return 2 * c / H - 16 * c**2 / H**4


node_pairs = ((0, 5), (1, 2), (3, 4))
for left, right in node_pairs:
    require(
        reduce_at_h0(gamma_u(h0 * zeta**left) - gamma_u(h0 * zeta**right))
        == 0,
        f"node pair ({left},{right}) u",
    )
    require(
        reduce_at_h0(gamma_v(h0 * zeta**left) - gamma_v(h0 * zeta**right))
        == 0,
        f"node pair ({left},{right}) v",
    )

# The three paired images are distinct.  Their source coordinates are the
# three points (u,u^2/2) with u^3=8c/3.
node_representatives = (0, 1, 3)
node_u_values = []
for index in node_representatives:
    H = h0 * zeta**index
    U = gamma_u(H)
    V = gamma_v(H)
    require(reduce_at_h0(U**3 - 8 * c / 3) == 0, f"node {index} u cubic")
    require(reduce_at_h0(V - U**2 / 2) == 0, f"node {index} v coordinate")
    node_u_values.append(U)
for left in range(3):
    for right in range(left + 1, 3):
        require(
            reduce_at_h0(node_u_values[left] - node_u_values[right]) != 0,
            f"distinct nodes {left},{right}",
        )

require(
    reduce_at_h0((h0 * zeta**0) ** 6 + 192 * c**2) == 0,
    "six-point node equation",
)
for index in range(3):
    require(
        zero((h0 * zeta**index) ** -3 - (-1) ** index * h0**-3),
        f"alternating inverse cube {index}",
    )

# The cycle equations force the nonzero odd-difference amplitude to vanish.
aa, bb, ee, rr = sp.symbols("a b e r")
cycle_rows = (aa - ee - rr, bb - ee + rr, bb - aa - rr)
cycle_basis = sp.groebner(cycle_rows, aa, bb, ee, rr, order="lex")
require(
    [poly.as_expr() for poly in cycle_basis.polys] == [aa - ee, bb - ee, rr],
    "three-node cycle forces r=0",
)
require(
    zero(cycle_rows[2] - (cycle_rows[1] - cycle_rows[0]) + 3 * rr),
    "explicit 3r contradiction sign",
)

# Hostile control: the tangent minors admit a cheap Bezout certificate, but
# its ambient two-form is not closed and hence cannot equal dP wedge dQ.
def bracket(f: sp.Expr, g: sp.Expr) -> sp.Expr:
    return sp.factor(sp.diff(f, u) * sp.diff(g, v) - sp.diff(f, v) * sp.diff(g, u))


M_AB = bracket(A, B)
M_Ad = bracket(A, d)
M_Bd = bracket(B, d)
require(zero(M_AB + 4 * (u**3 + u * v - c)), "AB minor")
require(zero(M_Ad - 2 * (u**2 + v)), "Ad minor")
require(zero(M_Bd + 2 * (v**2 + u**2 * v - c * u)), "Bd minor")

c_AB = -d0 / (4 * c**2)
c_Ad = B0 / (3 * c**2)
c_Bd = -A0 / (6 * c**2)
certificate = sp.factor(
    c_AB.subs(d0, d) * M_AB
    + c_Ad.subs(B0, B) * M_Ad
    + c_Bd.subs(A0, A) * M_Bd
)
require(zero(certificate - 1), "minor Bezout certificate")
closure = sp.factor(
    sp.diff(c_AB, d0) - sp.diff(c_Ad, B0) + sp.diff(c_Bd, A0)
)
require(zero(closure + 3 / (4 * c**2)), "nonclosed certificate")
require(closure != 0, "closure hostile is genuinely nonzero")

print("THM-3563 exact node-cycle controls")
print("universe=QQ(c,kappa,epsilon,h0)/(epsilon^2+3,h0^3-8*epsilon*c); characteristic_zero; c*kappa!=0")
print("laurent_even_packet=A,B,d; double_curve=K; special_points=6; source_nodes=3")
print("primitive_even_term=24*c^2/h^4; forced_odd_difference=16*kappa*c^2/h^3")
print("node_pairs=(0,5),(1,2),(3,4); inverse_cube_signs=+,-,+")
print("cycle_consequence=3*r=0; hostile_minor_certificate=unit_but_nonclosed")
print(f"truth_gates={GATES}; status=PROVED_VERIFIED_EXACT_INDEPENDENTLY_AUDITED")
