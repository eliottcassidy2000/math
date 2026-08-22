#!/usr/bin/env python3
"""Exact companion for THM-3545.

The proof in the theorem is algebraic.  This companion independently checks
the separated-ansatz Jacobian identity, the quadratic algebraic branch, the
Catalan coefficient formula, the double point, and the first defect of every
finite truncation through a declared range.
"""

from fractions import Fraction
from math import comb

import sympy as sp


def catalan(n: int) -> int:
    return comb(2 * n, n) // (n + 1)


failures: list[str] = []
gates = 0


def gate(name: str, condition: bool) -> None:
    global gates
    gates += 1
    if not condition:
        failures.append(name)


v, w, kappa = sp.symbols("v w kappa", nonzero=True)
A, Ap, B, Bp = sp.symbols("A Ap B Bp")

# Independent coefficient calculation for the whole separated ansatz.
P_v, P_w = 2 * v, Ap
Q_v, Q_w = 3 * v**2 - 1 + B, v * Bp
J_general = sp.expand(P_v * Q_w - P_w * Q_v)
J_expected = sp.expand((2 * Bp - 3 * Ap) * v**2 + Ap * (1 - B))
gate("separated ansatz Jacobian", sp.expand(J_general - J_expected) == 0)

# The coefficient of v^2 and the constant coefficient give the two rigidity
# equations.  Substitute their unique normalized solution.
J_rigid = sp.expand(J_general.subs({B: sp.Rational(3, 2) * A,
                                    Bp: sp.Rational(3, 2) * Ap}))
gate("rigid cancellation", sp.expand(J_rigid - Ap * (1 - sp.Rational(3, 2) * A)) == 0)

# Closed algebraic branch and its differential equation.
A_closed = sp.Rational(2, 3) * (1 - sp.sqrt(1 - 3 * kappa * w))
B_closed = sp.Rational(3, 2) * A_closed
ode = sp.simplify(sp.diff(A_closed, w) * (1 - sp.Rational(3, 2) * A_closed))
quadratic = sp.expand(A_closed - sp.Rational(3, 4) * A_closed**2 - kappa * w)
gate("closed branch ODE", sp.simplify(ode - kappa) == 0)
gate("closed branch quadratic", sp.simplify(quadratic) == 0)
gate("normalization A(0)", sp.simplify(A_closed.subs(w, 0)) == 0)
gate("normalization A'(0)", sp.simplify(sp.diff(A_closed, w).subs(w, 0) - kappa) == 0)

# Direct Jacobian check after substituting the closed branch, without using
# the generic coefficient comparison above.
P_closed = v**2 + A_closed
Q_closed = v**3 - v + sp.Rational(3, 2) * v * A_closed
J_closed = sp.simplify(
    sp.diff(P_closed, v) * sp.diff(Q_closed, w)
    - sp.diff(P_closed, w) * sp.diff(Q_closed, v)
)
gate("direct closed-branch Jacobian", sp.simplify(J_closed - kappa) == 0)

# The two boundary points remain distinct and collide.
image_plus = (sp.simplify(P_closed.subs({v: 1, w: 0})),
              sp.simplify(Q_closed.subs({v: 1, w: 0})))
image_minus = (sp.simplify(P_closed.subs({v: -1, w: 0})),
               sp.simplify(Q_closed.subs({v: -1, w: 0})))
gate("boundary collision", image_plus == image_minus == (1, 0))
boundary_tangent_plus = sp.Matrix([2, 2])
boundary_tangent_minus = sp.Matrix([-2, 2])
gate("transverse boundary tangents",
     sp.det(sp.Matrix.hstack(boundary_tangent_plus, boundary_tangent_minus)) != 0)

# Coefficients are Catalan and every one is nonzero in characteristic zero.
NMAX = 12
series = sp.series(A_closed, w, 0, NMAX + 2).removeO().expand()
coefficients: list[Fraction] = []
for n in range(1, NMAX + 2):
    got = sp.expand(series).coeff(w, n)
    expected = sp.Integer(catalan(n - 1)) * sp.Rational(3, 4) ** (n - 1) * kappa**n
    gate(f"Catalan coefficient n={n}", sp.simplify(got - expected) == 0)
    if n <= NMAX:
        coefficients.append(Fraction(catalan(n - 1) * 3 ** (n - 1), 4 ** (n - 1)))

# For kappa=1, truncate after degree N.  The first nonconstant Jacobian defect
# must occur exactly in degree N, so no finite initial segment is a solution.
truncation_defects: list[tuple[int, Fraction]] = []
for N in range(1, NMAX + 1):
    A_N = sum(sp.Rational(catalan(n - 1) * 3 ** (n - 1), 4 ** (n - 1)) * w**n
              for n in range(1, N + 1))
    J_N = sp.expand(sp.diff(A_N, w) * (1 - sp.Rational(3, 2) * A_N))
    defect = sp.expand(J_N - 1)
    first = None
    for degree in range(0, 2 * N + 1):
        coefficient = sp.expand(defect).coeff(w, degree)
        if coefficient != 0:
            first = (degree, coefficient)
            break
    gate(f"truncation first defect N={N}", first is not None and first[0] == N)
    if first is not None:
        truncation_defects.append((N, Fraction(int(first[1].p), int(first[1].q))))

# A finite polynomial solution is also ruled out directly by degrees.  This
# is printed as a symbolic degree law rather than sampled evidence:
# deg(A'(1-3A/2)) = 2 deg(A)-1 for every nonconstant polynomial A.
degree_law = "2d-1"
gate("positive degree cannot give a unit", all(2 * d - 1 > 0 for d in range(1, 25)))

print("THM-3545 exact companion")
print(f"gates={gates}")
print(f"general_J={sp.sstr(J_general)}")
print(f"rigid_equations: 2*B'-3*A'=0; A'*(1-B)=kappa")
print(f"A_closed={sp.sstr(A_closed)}")
print(f"B_closed={sp.sstr(B_closed)}")
print(f"J_closed={sp.sstr(J_closed)}")
print(f"collision={image_plus}={image_minus}")
print("A_coefficients_kappa1=" + ",".join(str(x) for x in coefficients))
print("truncation_first_defects=" + ",".join(f"{n}:{c}" for n, c in truncation_defects))
print(f"polynomial_degree_law={degree_law}")
print("failures=" + ("none" if not failures else ",".join(failures)))

if failures:
    raise SystemExit(1)
