#!/usr/bin/env python3
"""Exact b=2 centralizer/cusp Cech audit for the S6 manuscript.

The integer calculations below use only the displayed monodromy, filling
vectors, and toric cusp shear.  The final H^1 vanishing is the standard
cohomological consequence of the manuscript's displayed exact sequence

    0 -> O -> V -> O(-1) -> 0.

Thus the completed-orbit conclusion is conditional on the manuscript's
analytic pieces existing as stated and on V being the translation sheaf for
their overlaps.  Nothing here verifies compactness, topology, or a complex
structure on S6.
"""

from sympy import Matrix, eye, symbols


def require(label, condition):
    if not condition:
        raise RuntimeError(f"FAIL  {label}")
    print(f"PASS  {label}")


T1 = Matrix(
    [[1, 0, -6, 2], [0, -1, 1, 1], [0, -1, 0, 1], [0, 0, 0, 1]]
)
T2 = Matrix(
    [[1, 6, 0, -3], [0, 0, -1, 1], [0, 1, 0, 0], [0, 0, 0, 1]]
)
A1 = T1.inv().T
A2 = T2.inv().T
A0 = (A1 * A2).inv()
I4 = eye(4)

gamma_hat = Matrix([1, 0, 0, 0])
u_hat = Matrix([0, 1, 0, 0])
w_hat = Matrix([0, 0, 1, 0])
delta_hat = Matrix([0, 0, 0, 1])

E41 = Matrix.zeros(4)
E41[3, 0] = 1
k = symbols("k", integer=True)
Q_even = I4 + 2 * k * E41

require("A1*A2*A0=I", A1 * A2 * A0 == I4)
require(
    "Q_(2k) centralizes both elliptic monodromies",
    Q_even * A1 == A1 * Q_even and Q_even * A2 == A2 * Q_even,
)
require("Q_(2k) is unimodular", Q_even.det() == 1)
require(
    "Q_(2k) fixes Lambda_tor pointwise",
    Q_even * w_hat == w_hat and Q_even * delta_hat == delta_hat,
)
require(
    "Q_(2k) is identity on Lambda/Lambda_tor",
    Q_even * gamma_hat - gamma_hat == 2 * k * delta_hat
    and Q_even * u_hat == u_hat,
)

# beta -> beta+2k changes C(t) by 2k E21.  On an integral quotient
# coordinate lambda_bar=(a,b), the change is (0,2ka), whose componentwise
# exponential is one.  Hence c_lambda(t), Psi_lambda, N0, E0, and G0 are
# literally unchanged.
E21 = Matrix.zeros(2)
E21[1, 0] = 1
a, b = symbols("a b", integer=True)
cusp_shift = 2 * k * E21 * Matrix([a, b])
require("cusp C-shift is the integral vector (0,2ka)", cusp_shift == Matrix([0, 2 * k * a]))
print("RESULT integer c0 shifts leave every exp(2*pi*i*C(t)*lambda_bar) unchanged")
print("RESULT the marked cusp quotient N0 and overlap G0 are literally unchanged")

# The source uses v1=epsilon and v2=-epsilon'.  Rewriting the c0+2k
# filling through the common punctured family replaces vj by Q_(2k)vj.
v1 = Matrix([1, 2, -4, 0])
v2 = Matrix([-1, -3, 3, 0])
d1 = Q_even * v1 - v1
d2 = Q_even * v2 - v2
require("order-three affine numerator shifts by 2k delta", d1 == 2 * k * delta_hat)
require("order-four affine numerator shifts by -2k delta", d2 == -2 * k * delta_hat)

# Translation by r conjugates x -> A*x+v/m to the target action exactly
# when d/m=(I-A)r+lambda with lambda integral.  These choices also expose
# the boundary period cocycles kappa=(I-A)r-d/m.
r1 = k * Matrix([1, -6, 0, 0]) / 12
r2 = k * Matrix([1, 0, 0, 0]) / 6
lambda1 = k * u_hat
lambda2 = -k * w_hat
require(
    "explicit order-three translation conjugator",
    d1 / 3 == (I4 - A1) * r1 + lambda1,
)
require(
    "explicit order-four translation conjugator",
    d2 / 4 == (I4 - A2) * r2 + lambda2,
)
kappa1 = (I4 - A1) * r1 - d1 / 3
kappa2 = (I4 - A2) * r2 - d2 / 4
require("first elliptic overlap cocycle is -k u_hat", kappa1 == -k * u_hat)
require("second elliptic overlap cocycle is k w_hat", kappa2 == k * w_hat)

# For g1*g2*g0=1 and left cocycles, the relation is
# kappa1+A1*kappa2+A1*A2*kappa0=0.
kappa0 = -A0 * (kappa1 + A1 * kappa2)
require("triangle relation forces cusp cocycle k w_hat", kappa0 == k * w_hat)
require(
    "three boundary cocycles satisfy the exact triangle relation",
    kappa1 + A1 * kappa2 + A1 * A2 * kappa0 == Matrix.zeros(4, 1),
)

# In quotient coordinates (gamma_bar,u_bar) and torus coordinates
# (w_hat,delta_hat), B0 sends -k*u_bar to -k*w_hat.  Since s o g0=s-1,
# Phi_(-k*u_bar) has boundary lift -k*s*w_hat and cocycle +k*w_hat.
B0 = Matrix([[0, 1], [-1, 0]])
require("A0 fixes the cusp vanishing cycle w_hat", A0 * w_hat == w_hat)
require(
    "toric shear Phi_(-k u_bar) has exponent -k w_hat",
    B0 * Matrix([0, -k]) == Matrix([-k, 0]),
)
print("RESULT Phi_(-k*u_bar) absorbs the forced cusp cocycle k*w_hat")
print("RESULT the discrete completed-overlap Cech obstruction vanishes for every b=2k")

# The local order-four invariant covector has psi2(delta)=2.  With ell2 odd,
# its residue is 2*b*ell2 mod 4: zero exactly for even b.
require("even shifts pass the order-four residue", (2 * (2 * k)) % 4 == 0)
require("odd shifts fail the order-four residue by 2", (2 * (2 * k + 1)) % 4 == 2)
print("RESULT odd centralizer shifts fail the marked order-four filling residue by 2 mod 4")

# Standard P1 line-bundle cohomology: h^1(O(d))=max(-d-1,0).
def h1_line(d):
    return max(-d - 1, 0)


require("H^1(P1,O)=0", h1_line(0) == 0)
require("H^1(P1,O(-1))=0", h1_line(-1) == 0)
print("CONDITIONAL INPUT manuscript exact sequence 0 -> O -> V -> O(-1) -> 0")
print("CONDITIONAL RESULT H^1(P1,V)=0 kills the connected additive overlap class")
print("CONDITIONAL RESULT the centralizer-generated completed orbit is exactly 2Z")
print("NONCONSEQUENCE no compactness, topology, homology, diffeomorphism, or S6 claim is verified")
print("ALL EXACT CHECKS PASSED")
