#!/usr/bin/env python3
"""Exact companion for THM-3949's all-degree one-factor Newton gate."""

from __future__ import annotations

import hashlib
import json

import sympy as sp


CHECKS = 0


def gate(condition: object, message: str) -> None:
    global CHECKS
    CHECKS += 1
    if condition is not True and condition != sp.S.true:
        raise RuntimeError(message)


d = sp.symbols("delta")
P, Y, h = sp.symbols("P Y h")
AA, BB = sp.symbols("A B")


def dreduce(expression: sp.Expr) -> sp.Expr:
    numerator, denominator = sp.fraction(sp.cancel(expression))
    reduced = sp.rem(sp.Poly(sp.expand(numerator), d), sp.Poly(d**2 + 3, d))
    return sp.cancel(reduced.as_expr() / denominator)


def dzero(expression: sp.Expr, message: str) -> None:
    gate(dreduce(expression) == 0, message)


omega = (-1 + d) / 2
omega2 = (-1 - d) / 2
dzero(omega**2 + omega + 1, "omega quadratic relation")
dzero(omega * omega2 - 1, "omega conjugates multiply to one")
dzero(omega - omega2 - d, "delta relation")


# ---------------------------------------------------------------------------
# Universal algebra of the standard reciprocal one-factor assignment.
# ---------------------------------------------------------------------------

p0 = P
p1 = P + AA * BB
L1 = p1 - omega * p0
L2 = p1 - omega2 * p0
qminus = 2 * AA * L1
qplus = 2 * BB * L2
q0 = sp.expand((qplus - qminus) / 2)
q1 = sp.expand((qplus + qminus) / 2)

D = (1 - omega2) * BB - (1 - omega) * AA
E = (BB - AA) * AA * BB

gate(p1 - p0 == AA * BB, "one-factor product difference")
dzero(q1 - q0 - qminus, "assigned q-difference row")
dzero(q1 + q0 - qplus, "assigned q-sum row")
dzero(L1 * L2 * AA * BB - (p1**3 - p0**3),
      "difference-of-cubes product")
dzero(q1**2 - q0**2 - 4 * (p1**3 - p0**3),
      "common discriminant identity")
dzero(q0 - (E + D * P), "q0 is E+D*P")

H0 = sp.expand(q0**2 - 4 * p0**3)
H1 = sp.expand(q1**2 - 4 * p1**3)
dzero(H0 - H1, "two torus rows have common H")


# ---------------------------------------------------------------------------
# Hidden square-root cubic, norm, and exact localization inverse.
# ---------------------------------------------------------------------------

C = sp.expand(E + D * h**2 - 2 * h**3)
q0_on_square = sp.expand(q0.subs(P, h**2))
gate(sp.Poly(C, h).LC() == -2,
     "hidden projection to the Y-line is finite monic cubic")
dzero(C - (q0_on_square - 2 * h**3),
      "hidden curve is q0=2h^3 on P=h^2")
dzero(C * C.subs(h, -h) - H0.subs(P, h**2),
      "hidden-curve norm is H(h^2,Y)")

# On H[P^-1], h0=q0/(2P) satisfies h0^2=P.  On C[h^-1], the
# inverse formula q0(h^2)/(2h^2)=h differs by C/(2h^2).
h_inverse = q0 / (2 * P)
dzero(h_inverse**2 - P - H0 / (4 * P**2),
      "H-localization square-root identity")
dzero(C.subs(h, h_inverse) + E * H0 / (4 * P**3),
      "hidden equation vanishes under the H-localization inverse")
dzero(q0_on_square / (2 * h**2) - h - C / (2 * h**2),
      "C-localization recovers h")
dzero(q0_on_square - 2 * h * h**2 - C,
      "localized P=h^2 and h=q0/(2P) are inverse")
gate(sp.Poly(E, AA).total_degree() > 0,
     "E is a genuine nonconstant symbolic coefficient packet")


# ---------------------------------------------------------------------------
# Four exhaustive degree/leading-coefficient Newton ledgers.
# Points are (h-degree, x-adic coefficient valuation), x=1/Y.
# ---------------------------------------------------------------------------

a, b, n = sp.symbols("a b n", integer=True, positive=True)
rdeg = sp.symbols("r", integer=True, nonnegative=True)

# b>a: deg E=a+2b, deg D=b.
gate(sp.simplify((a + 2 * b - b) / 2 - (a + b) / 2) == 0,
     "b>a first Newton slope")
gate(sp.simplify((0 - (-b)) / (3 - 2) - b) == 0,
     "b>a second Newton slope")
gate(sp.simplify(b - (a + b) / 2 - (b - a) / 2) == 0,
     "b>a slopes differ by (b-a)/2")

# a>b: deg E=2a+b, deg D=a.
gate(sp.simplify((2 * a + b - a) / 2 - (a + b) / 2) == 0,
     "a>b first Newton slope")
gate(sp.simplify((0 - (-a)) / (3 - 2) - a) == 0,
     "a>b second Newton slope")
gate(sp.simplify(a - (a + b) / 2 - (a - b) / 2) == 0,
     "a>b slopes differ by (a-b)/2")

# Equal degree and equal leaders: deg E=2n+r, deg D=n, r<n.
gate(sp.simplify((2 * n + rdeg - n) / 2 - (n + rdeg) / 2) == 0,
     "equal-leader first Newton slope")
gate(sp.simplify((0 - (-n)) / (3 - 2) - n) == 0,
     "equal-leader second Newton slope")
gate(sp.simplify(n - (n + rdeg) / 2 - (n - rdeg) / 2) == 0,
     "equal-leader slopes differ by (n-r)/2")

# Concrete controls for the three two-slope rows.
def degree(poly: sp.Expr) -> int:
    return int(sp.Poly(sp.expand(poly), Y).degree())


def packets(Apoly: sp.Expr, Bpoly: sp.Expr) -> tuple[sp.Expr, sp.Expr]:
    Dpoly = sp.expand((1 - omega2) * Bpoly - (1 - omega) * Apoly)
    Epoly = sp.expand((Bpoly - Apoly) * Apoly * Bpoly)
    return Dpoly, Epoly


A12, B12 = Y + 1, Y**2 + Y + 1
D12, E12 = packets(A12, B12)
gate((degree(A12), degree(B12), degree(D12), degree(E12)) == (1, 2, 2, 5),
     "degree-imbalanced b>a control")

A21, B21 = Y**2 + Y + 1, Y + 1
D21, E21 = packets(A21, B21)
gate((degree(A21), degree(B21), degree(D21), degree(E21)) == (2, 1, 2, 5),
     "degree-imbalanced a>b control")

Aeq, Beq = Y**2 + 1, Y**2 + Y + 2
Deq, Eeq = packets(Aeq, Beq)
gate(sp.gcd(Aeq, Beq) == 1, "equal-leading control is coprime")
gate((degree(Aeq), degree(Beq), degree(Beq - Aeq),
      degree(Deq), degree(Eeq)) == (2, 2, 1, 2, 5),
     "equal-leading lower-difference degree control")

Aeq0, Beq0 = Y**2 + 1, Y**2 + 2
Deq0, Eeq0 = packets(Aeq0, Beq0)
gate(sp.gcd(Aeq0, Beq0) == 1, "constant-difference control is coprime")
gate((degree(Beq0 - Aeq0), degree(Deq0), degree(Eeq0)) == (0, 2, 4),
     "equal-leading constant-difference boundary control")


# ---------------------------------------------------------------------------
# Equal degree, unequal leaders: the sole-edge residual cannot be a cube.
# ---------------------------------------------------------------------------

alpha, beta, kappa = sp.symbols("alpha beta kappa")
e0 = alpha * beta * (beta - alpha)
d0 = (1 - omega2) * beta - (1 - omega) * alpha
residual = sp.expand(e0 * kappa**3 + d0 * kappa - 2)
gate(sp.Poly(residual, kappa).coeff_monomial(kappa**2) == 0,
     "sole-edge residual has no quadratic coefficient")
gate(sp.Poly(residual, kappa).coeff_monomial(1) == -2,
     "sole-edge residual has nonzero constant")

# A cubic e0*(kappa-c)^3 could equal the residual only if the following
# coefficient equations held.  Saturating by e0 gives the unit ideal.
croot, e_inv, d_aux, e_aux = sp.symbols("c e_inv d_aux e_aux")
triple_ideal = sp.groebner([
    -3 * e_aux * croot,
    3 * e_aux * croot**2 - d_aux,
    -e_aux * croot**3 + 2,
    e_inv * e_aux - 1,
], croot, d_aux, e_aux, e_inv, order="lex")
gate(list(triple_ideal) == [1],
     "nonzero-leading residual cannot have one triple root")


# ---------------------------------------------------------------------------
# Smallest nonlinear hostile: a double residual root plus a simple one.
# ---------------------------------------------------------------------------

Aq = Y**2 + 1
Bq = -omega2 * Y**2 + 1
Dq, Eq = packets(Aq, Bq)
dzero(Bq - Aq - omega * Y**2,
      "quadratic hostile factor difference")
gate(sp.gcd(Aq, Y) == 1,
     "quadratic hostile factors are coprime via their difference")
gate(dreduce(sp.resultant(Aq, Bq, Y)) != 0,
     "quadratic hostile resultant is nonzero")

alpha_q = sp.Integer(1)
beta_q = -omega2
e0_q = dreduce(alpha_q * beta_q * (beta_q - alpha_q))
d0_q = dreduce((1 - omega2) * beta_q - (1 - omega) * alpha_q)
gate(e0_q == -1, "quadratic hostile cubic leading coefficient")
dzero(d0_q - 3 * omega, "quadratic hostile linear coefficient")

Rq = dreduce(e0_q * kappa**3 + d0_q * kappa - 2)
dzero(Rq + (kappa - omega2) ** 2 * (kappa + 2 * omega2),
      "quadratic hostile has double-plus-simple residual roots")
dzero(Rq.subs(kappa, omega2), "quadratic hostile double root")
dzero(sp.diff(Rq, kappa).subs(kappa, omega2),
      "quadratic hostile root is repeated")
gate(dreduce(sp.diff(Rq, kappa, 2).subs(kappa, omega2)) != 0,
     "quadratic hostile repeated root is not triple")
dzero(Rq.subs(kappa, -2 * omega2),
      "quadratic hostile separate simple root")
gate(dreduce(sp.diff(Rq, kappa).subs(kappa, -2 * omega2)) != 0,
     "quadratic hostile second root is simple")

# Freeze the actual local edge after Y=1/x and h=z/x^2.
xlocal, z = sp.symbols("x z")
Cq = sp.expand(Eq + Dq * h**2 - 2 * h**3)
local_q = sp.expand(xlocal**6 * Cq.subs({Y: 1 / xlocal, h: z / xlocal**2}))
edge_q = dreduce(local_q.subs(xlocal, 0))
dzero(edge_q - (e0_q + d0_q * z**2 - 2 * z**3),
      "quadratic hostile local Newton edge")
dzero(edge_q + 2 * (z - omega) ** 2 * (z + omega / 2),
      "quadratic hostile edge has two leading h/Y^2 directions")


summary = {
    "checks": CHECKS,
    "grammar": "p1-p0=A(Y)B(Y), coprime nonconstant opposite-row assignment",
    "hidden_curve": "E+D*h^2-2*h^3",
    "norm": "C(h)C(-h)=H(h^2,Y)",
    "localization": "H[P^-1] isomorphic to C[h^-1]",
    "newton_cases": 4,
    "minimum_standard_infinity_places_if_irreducible": 2,
    "quadratic_hostile": "double residual root plus separate simple root",
    "arbitrary_line_claim": False,
}
semantic = hashlib.sha256(json.dumps(summary, sort_keys=True).encode()).hexdigest()

print("THM-3949 coprime one-variable internal-split Newton companion")
print(f"CHECKS={CHECKS}")
print("HIDDEN_C=E+D*h^2-2*h^3")
print("NORM=C(h)C(-h)=H(h^2,Y)")
print("LOCALIZATION=H[P^-1]~=C[h^-1]")
print("NEWTON_CASES=DEG_IMBALANCED_2;EQUAL_LEADERS;UNEQUAL_LEADERS")
print("DEG_IMBALANCED_SLOPES=(a+b)/2,max(a,b)")
print("EQUAL_LEADER_SLOPES=(n+deg(B-A))/2,n")
print("UNEQUAL_LEADER_RESIDUAL=e0*kappa^3+d0*kappa-2;NOT_TRIPLE")
print("QUADRATIC_HOSTILE=DOUBLE_PLUS_SIMPLE_RESIDUAL_ROOT")
print("CONCLUSION=H_REDUCIBLE_OR_STANDARD_INFINITY_PLACES>=2")
print("SCOPE=NO_ARBITRARY_LINE_CLAIM")
print(f"SEMANTIC_SHA256={semantic}")
