"""Exact referee for repaired THM-3025 (MISTAKE-422).

Typed notation:
    J0 = Jac(P,Q) in k*              (the full Keller constant),
    j  = Jac(H,Q_(m-1))              (the possibly-zero subleading form).

THM-3016 proves j*Jac(W,H)=0.  This replay does not divide by j in its zero
branch.  It checks the common-form degree obstruction by exact linear
nullities and checks the relevant coprimality identities.
"""

from math import gcd

import sympy as sp


x, y = sp.symbols("x y")


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def jac(f, g):
    return sp.expand(sp.diff(f, x) * sp.diff(g, y) - sp.diff(f, y) * sp.diff(g, x))


def homogeneous_nullity(H, degree):
    coeffs = sp.symbols(f"c0:{degree + 1}")
    F = sum(coeffs[i] * x**i * y ** (degree - i) for i in range(degree + 1))
    polynomial = sp.Poly(jac(F, H), x, y)
    equations = polynomial.coeffs()
    matrix, _ = sp.linear_eq_to_matrix(equations, coeffs)
    return len(coeffs) - matrix.rank()


print("STEP 1. Keep the two Jacobians typed.")
print("   J0=Jac(P,Q) is the nonzero Keller constant.")
print("   j=Jac(H,Q_(m-1)) is a subleading form and may be zero.")
print("   THM-3016 gives j*Jac(W,H)=0, so a j=0/j!=0 split is mandatory.")

print("\nSTEP 2. Exact coprime-degree gate.")
bad_n = []
bad_m = []
for g in range(1, 60):
    for a in range(1, 60):
        if gcd(g, g * a - 1) != 1:
            bad_n.append((g, a))
    for b in range(1, 60):
        if gcd(g, g * b - 1) != 1:
            bad_m.append((g, b))
require(not bad_n and not bad_m, "coprime-degree gate failed")
print("   gcd(g,ga-1)=gcd(g,gb-1)=1 for 1<=g,a,b<60: PASS")

print("\nSTEP 3. K>=2 hostile nullities in degree ga-1 (the j!=0 branch).")
hostiles = [
    ((x + y) * (x - y), 2, 2, "(x+y)(x-y), g=2, a=2"),
    ((x + y) * (x - y), 2, 3, "(x+y)(x-y), g=2, a=3"),
    ((x + y) ** 2 * (x - y), 3, 2, "(x+y)^2(x-y), g=3, a=2"),
    (x * y * (x + y), 3, 3, "xy(x+y), g=3, a=3"),
    ((x + y) * (x - y) * (x + 2 * y), 3, 4, "three roots, g=3, a=4"),
]
for H, g, a, label in hostiles:
    degree = g * a - 1
    nullity = homogeneous_nullity(H, degree)
    require(nullity == 0, f"unexpected K>=2 kernel for {label}")
    print(f"   {label:38s} degree={degree:2d}, nullity={nullity}: PASS")

print("\nSTEP 4. The j=0 branch closes both subleading forms.")
H = x * y
q_nullity = homogeneous_nullity(H, 3)  # g=2, m=4, degree m-1
p_nullity = homogeneous_nullity(H, 5)  # g=2, n=6, degree n-1
require(q_nullity == 0 and p_nullity == 0, "j=0 branch hostile failed")
print("   H=xy: Jac(H,Q_3)=0 has nullity 0: Q_3=0")
print("   H=xy: Jac(P_5,H)=0 has nullity 0: P_5=0")
print("   hence W=0 without dividing by j: PASS")

print("\nSTEP 5. K=1 sharp controls.")
controls = [((x + y) ** 2, 3, 1), ((x + y) ** 3, 5, 1)]
for H, degree, expected in controls:
    nullity = homogeneous_nullity(H, degree)
    require(nullity == expected, "K=1 sharp control failed")
    print(f"   H={sp.factor(H)}, degree={degree}: nullity={nullity} (nonzero): PASS")

print("\nSTEP 6. Symbolic product identity controls.")
j, w0, w1 = sp.symbols("j w0 w1")
W = w0 * x**3 + w1 * x**2 * y
H = x * y
product = sp.expand(j * jac(W, H))
require(product.subs(j, 0) == 0, "zero branch is not tautological")
require(product.subs({j: 1, w0: 1, w1: 0}) != 0, "nonzero branch hostile failed")
print("   j=0 annihilates the product for arbitrary W: PASS")
print("   j!=0 does not annihilate a hostile W automatically: PASS")

print("\nALL CHECKS PASSED")
