#!/usr/bin/env python3
"""Exact controls for THM-3116.

This checks the three-variable radial/projective Jacobian and Gamma shifts,
the E_{D,3}/E_{D,4} coefficient bookkeeping, the finite multinomial radial
formula, the affine-simplex exponential formula (including all confluent
cases), a derivative-aligned quadratic family, and two hostile controls:

* q=u+v has the algebraic period 1, so a blanket two-dimensional
  transcendence extension is false.  Its parameter integral is an
  E-function, but the exceptional value at z=1 is already explained by an
  exact functional relation with exp(z); this is a hostile control against
  using E-function status without auditing the relation module;
* q=a(u+v), a=1+W_1(-e^{-1}), has zero simplex exponential integral, so
  positivity of the representing measure alone does not prevent complex
  cancellation.  Hermite--Lindemann forces this nonzero a to be
  transcendental.

Reproduce:
  python3 04-computation/fc3_flat_top_simplex_asymptotic_thm3116.py
  python3 -O 04-computation/fc3_flat_top_simplex_asymptotic_thm3116.py
"""

from __future__ import annotations

from math import factorial

import mpmath as mp
import sympy as sp


u, v, r, zz = sp.symbols("u v r zz")
w = 1 - u - v


def check(condition: bool, message: str) -> None:
    """Optimization-stable exact check (unlike Python's assert)."""
    if not condition:
        raise RuntimeError(message)


def simplex_integral(poly: sp.Expr) -> sp.Expr:
    """Exact coordinate-area integral on {(u,v)>=0:u+v<=1}."""
    return sp.simplify(sp.integrate(sp.integrate(poly, (v, 0, 1 - u)), (u, 0, 1)))


def factorial_functional(poly: sp.Expr, xs: tuple[sp.Symbol, ...]) -> sp.Expr:
    p = sp.Poly(sp.expand(poly), *xs)
    total = sp.Integer(0)
    for powers, coeff in p.terms():
        weight = sp.Integer(1)
        for exponent in powers:
            weight *= factorial(exponent)
        total += coeff * weight
    return sp.expand(total)


print("THM-3116 FC(3) FLAT-TOP SIMPLEX AUDIT")

# C1. Jacobian and the monomial/Dirichlet normalization.
x_expr, y_expr, z_expr = r * u, r * v, r * w
jac = sp.factor(sp.det(sp.Matrix([
    [sp.diff(x_expr, t) for t in (r, u, v)],
    [sp.diff(y_expr, t) for t in (r, u, v)],
    [sp.diff(z_expr, t) for t in (r, u, v)],
])))
check(jac == r**2, "three-variable radial Jacobian")

dirichlet_checks = 0
for aa in range(5):
    for bb in range(5):
        for cc in range(5):
            angular = simplex_integral(u**aa * v**bb * w**cc)
            expected = sp.Rational(
                factorial(aa) * factorial(bb) * factorial(cc),
                factorial(aa + bb + cc + 2),
            )
            check(angular == expected, f"Dirichlet integral {(aa, bb, cc)}")
            radial = factorial(aa + bb + cc + 2)
            check(
                sp.expand(radial * angular) == factorial(aa) * factorial(bb) * factorial(cc),
                f"radial times angular monomial {(aa, bb, cc)}",
            )
            dirichlet_checks += 1

print(f"C1 jacobian={jac}; exact Dirichlet/Gamma monomials={dirichlet_checks}")

# C2. Mittag--Leffler coefficient shifts.  Multiplication by the radial
# Gamma integral must return L(f^n)/Gamma(Dn+3), and one extra divisor gives
# Gamma(Dn+4).
mittag_checks = 0
for D in range(1, 7):
    for n in range(0, 9):
        radial_gamma = sp.gamma(D * n + 3)
        check(
            sp.simplify(radial_gamma / sp.gamma(D * n + 3)) == 1,
            f"E_(D,3) shift {(D, n)}",
        )
        check(
            sp.simplify(
                radial_gamma / sp.gamma(D * n + 4) - sp.Rational(1, D * n + 3)
            ) == 0,
            f"E_(D,4) shift {(D, n)}",
        )
        mittag_checks += 1
print(f"C2 Gamma(Dn+3), E_(D,3), E_(D,4) coefficient checks={mittag_checks}")

# C3. Exact finite radial multinomial formula against a direct factorial
# expansion.  This is deliberately nonradial and contains every lower layer.
x, y, z = sp.symbols("x y z")
D = 2
f = (x + y + z) ** 2 + (2 * x - y + 3 * z) + 5
A = [
    sp.Integer(1),
    sp.expand(2 * u - v + 3 * w),
    sp.Integer(5),
]
radial_formula_checks = 0
for n in range(0, 7):
    direct = factorial_functional(f**n, (x, y, z))
    angular_sum = sp.Integer(0)
    for k1 in range(n + 1):
        for k2 in range(n - k1 + 1):
            k0 = n - k1 - k2
            J = k1 + 2 * k2
            multinomial = sp.Rational(factorial(n), factorial(k0) * factorial(k1) * factorial(k2))
            angular_sum += (
                multinomial
                * sp.gamma(D * n + 3 - J)
                * simplex_integral(A[1] ** k1 * A[2] ** k2)
            )
    check(sp.expand(angular_sum) == direct, f"finite radial formula n={n}")
    radial_formula_checks += 1
print(f"C3 direct factorial == finite radial-layer formula, n=0..6: {radial_formula_checks}")

# C4. Affine exponential integral formulas, derived by exact integration.
a, b, c = sp.symbols("a b c", nonzero=True)
generic = 1 / (a * b) + sp.exp(a) / (a * (a - b)) + sp.exp(b) / (b * (b - a))
one_axis = (sp.exp(a) - 1 - a) / a**2
diagonal = ((a - 1) * sp.exp(a) + 1) / a**2

# Check generic numerical algebraic points against direct symbolic integration.
affine_checks = 0
for av, bv in [(1, 2), (2, -1), (sp.Rational(1, 2), 3), (-2, -3)]:
    direct = simplex_integral(sp.exp(av * u + bv * v))
    # SymPy leaves some exponential integrals in a different elementary form.
    check(
        sp.simplify(direct - generic.subs({a: av, b: bv})) == 0,
        f"generic affine integral {(av, bv)}",
    )
    affine_checks += 1
check(sp.simplify(simplex_integral(sp.exp(a * u)) - one_axis) == 0, "one-axis confluence")
check(
    sp.simplify(simplex_integral(sp.exp(a * (u + v))) - diagonal) == 0,
    "diagonal confluence",
)
print(f"C4 affine divided-difference formulas exact; generic controls={affine_checks}")

# C4b.  For any barycentric coordinate lambda_i, s=1-lambda_i has
# pushforward density s ds.  The phase c+b*s^2 is therefore aligned with the
# derivative of its exponent and reduces to a boundary difference.  Check all
# three barycentric choices directly on the triangle.
c, b = sp.symbols("c b", nonzero=True)
s, t = sp.symbols("s t", real=True)
quadratic_boundary = sp.exp(c) * (sp.exp(b) - 1) / (2 * b)
charts = ((t, s - t), (1 - s, t), (t, 1 - s))
for i, ((up, vp), si) in enumerate(zip(charts, (u + v, 1 - u, 1 - v))):
    check(
        sp.simplify(si.subs({u: up, v: vp}) - s) == 0,
        f"derivative-aligned phase coordinate i={i}",
    )
    check(
        abs(sp.det(sp.Matrix([[sp.diff(up, s), sp.diff(up, t)],
                              [sp.diff(vp, s), sp.diff(vp, t)]]))) == 1,
        f"derivative-aligned chart Jacobian i={i}",
    )
quadratic_antiderivative = sp.exp(c + b * s**2) / (2 * b)
check(
    sp.simplify(sp.diff(quadratic_antiderivative, s) - s * sp.exp(c + b * s**2))
    == 0,
    "derivative-aligned quadratic antiderivative",
)
check(
    sp.simplify(quadratic_antiderivative.subs(s, 1)
                - quadratic_antiderivative.subs(s, 0)
                - quadratic_boundary) == 0,
    "derivative-aligned quadratic boundary value",
)
print("C4b derivative-aligned quadratics: 3 barycentric boundary formulas exact")

# C5. The forced homogeneous level is Area=1/2 at Q=0.  The
# two-dimensional transcendence extension is nevertheless false even for an
# affine algebraic phase: q=u+v gives exactly 1 (not merely a decimal).
forced_level_at_zero = simplex_integral(sp.Integer(1))
check(forced_level_at_zero == sp.Rational(1, 2), "simplex forced level")
algebraic_period = sp.simplify(sp.limit(diagonal, a, 1))
check(algebraic_period == 1, "algebraic nontranscendence hostile")
print(
    "C5 affine forced level: Q=0 gives 1/2; hostile nonconstant q=u+v "
    f"gives {algebraic_period} exactly"
)

# C5b.  The whole parameter integral for q=u+v is an E-function,
#
#   F(z)=sum_{m>=0} z^m/(m! (m+2))=((z-1)e^z+1)/z^2.
#
# Its algebraic specialization F(1)=1 is not a transcendence contradiction:
# the value relation is the specialization of a functional relation involving
# e^z.  Equivalently F satisfies a second-order polynomial-coefficient ODE.
# This isolates the extra datum needed by any E-function attack: the full
# functional-relation module, not holonomicity alone.
Fzz = ((zz - 1) * sp.exp(zz) + 1) / zz**2
check(
    sp.simplify(zz**2 * Fzz - (zz - 1) * sp.exp(zz) - 1) == 0,
    "affine E-function relation",
)
check(
    sp.simplify(zz * sp.diff(Fzz, zz, 2) + (3 - zz) * sp.diff(Fzz, zz) - 2 * Fzz)
    == 0,
    "affine E-function differential equation",
)
for m in range(9):
    check(
        sp.expand(Fzz.series(zz, 0, 10).removeO()).coeff(zz, m)
        == sp.Rational(1, factorial(m) * (m + 2)),
        f"affine E-function coefficient m={m}",
    )
check(sp.limit(Fzz, zz, 1) == 1, "affine E-function value at one")
print("C5b E-function hostile: F(1)=1 lifts from z^2 F-(z-1)e^z-1=0; ODE exact")

# C6. Positivity of the measure does not prevent cancellation for a complex
# phase.  The Lambert-W equation makes the numerator vanish exactly.
mp.mp.dps = 80
alpha = 1 + mp.lambertw(-mp.e ** -1, 1)
lambert_residual = (alpha - 1) * mp.e**alpha + 1
zero_integral = lambert_residual / alpha**2
check(abs(lambert_residual) < mp.mpf("1e-75"), "Lambert numerator hostile")
check(abs(zero_integral) < mp.mpf("1e-75"), "Lambert simplex integral hostile")
print("C6 Lambert-W hostile a=1+W_1(-e^-1)")
print(f"   a={mp.nstr(alpha, 32)}")
print(f"   numerator residual={mp.nstr(abs(lambert_residual), 8)}")
print(f"   simplex integral residual={mp.nstr(abs(zero_integral), 8)}")

# C7. A radial flat-top family gives an independently computable coefficient
# limit.  For f=(x+y+z)^D+A1*(x+y+z)^(D-1), all projective layers are
# constant, and the exact finite sum converges to Area*exp(A1/D).
mp.mp.dps = 60
Dnum = 3
A1num = mp.mpf("2")
target = mp.e ** (A1num / Dnum) / 2


def radial_normalized_moment(n: int) -> mp.mpf:
    total = mp.mpf("0")
    for k in range(n + 1):
        total += (
            mp.binomial(n, k)
            * A1num**k
            * mp.gamma(Dnum * n + 3 - k)
            / mp.gamma(Dnum * n + 3)
        )
    return total / 2


rows = []
for n in (10, 40, 160, 640):
    value = radial_normalized_moment(n)
    rows.append((n, value, abs(value - target)))
print(f"C7 flat coefficient target=Area*exp(A1/D)={mp.nstr(target, 24)}")
for n, value, error in rows:
    print(f"   n={n:3d} value={mp.nstr(value, 24)} error={mp.nstr(error, 8)}")
check(rows[-1][2] < rows[-2][2] < rows[-3][2] < rows[-4][2], "flat coefficient convergence")

print("ALL THM-3116 CONTROLS PASSED")
