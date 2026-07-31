from __future__ import annotations

from math import factorial

import sympy as sp


s = sp.symbols("s")


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def factorial_functional(poly: sp.Expr) -> sp.Expr:
    expanded = sp.Poly(sp.expand(poly), s)
    return sp.simplify(
        sum(
            coefficient * factorial(degree[0])
            for degree, coefficient in expanded.terms()
        )
    )


def monic_laguerre(D: int) -> sp.Expr:
    """Return (-1)^D D! times the standard Laguerre polynomial L_D."""
    return sp.expand(
        sum(
            (-1) ** (D + k)
            * factorial(D)
            * sp.binomial(D, k)
            * s**k
            / factorial(k)
            for k in range(D + 1)
        )
    )


laguerre_checks = []
for D in range(1, 13):
    ell = monic_laguerre(D)
    require(sp.Poly(ell, s).LC() == 1, f"Laguerre monicity failed at D={D}")
    require(
        tuple(factorial_functional(s**r * ell) for r in range(D)) == (0,) * D,
        f"Laguerre orthogonality failed at D={D}",
    )
    norm = factorial_functional(ell**2)
    require(norm == factorial(D) ** 2, f"Laguerre norm failed at D={D}")

    for degree in range(2 * D):
        _, remainder = sp.div(s**degree, ell, domain=sp.QQ)
        require(
            sp.degree(remainder, s) < D,
            f"quotient degree failed at D={D}, k={degree}",
        )
        require(
            factorial_functional(s**degree) == factorial_functional(remainder),
            f"quotient readout failed at D={D}, k={degree}",
        )

    laguerre_checks.append((D, norm))


def inverse_hankel_entry(d: int, i: int, j: int) -> sp.Rational:
    return sp.Rational(
        (-1) ** (i + j)
        * sum(
            sp.binomial(k, i) * sp.binomial(k, j)
            for k in range(max(i, j), d + 1)
        ),
        factorial(i) * factorial(j),
    )


selector_checks = []
for d in range(0, 11):
    hankel = sp.Matrix(
        [[factorial(i + j) for j in range(d + 1)] for i in range(d + 1)]
    )
    inverse = sp.Matrix(
        [
            [inverse_hankel_entry(d, i, j) for j in range(d + 1)]
            for i in range(d + 1)
        ]
    )
    require(hankel * inverse == sp.eye(d + 1), f"inverse formula failed at d={d}")
    for j in range(d + 1):
        phi = sum(inverse[i, j] * s**i for i in range(d + 1))
        require(
            tuple(factorial_functional(s**r * phi) for r in range(d + 1))
            == tuple(1 if r == j else 0 for r in range(d + 1)),
            f"dual selector failed at d={d}, j={j}",
        )
        expected_norm = sp.Rational(
            sum(sp.binomial(k, j) ** 2 for k in range(j, d + 1)),
            factorial(j) ** 2,
        )
        require(
            factorial_functional(phi**2) == expected_norm,
            f"selector norm failed at d={d}, j={j}",
        )
    selector_checks.append(d)

# The MISTAKE-211 two-height shell is scalar-null but coefficient-nonzero.
a = b = sp.Integer(1)
c = -sp.Rational(factorial(6), factorial(18))
G4 = 4 * a * b**3 * s**6 + 4 * a**3 * c * s**18
require(factorial_functional(G4) == 0, "MISTAKE-211 scalar hostile failed")
require(
    sp.Poly(G4, s).coeff_monomial(s**6) != 0
    and sp.Poly(G4, s).coeff_monomial(s**18) != 0,
    "MISTAKE-211 private coefficients unexpectedly vanished",
)

H18 = sp.Matrix([[factorial(i + j) for j in range(19)] for i in range(19)])
H18_inverse = sp.Matrix(
    [[inverse_hankel_entry(18, i, j) for j in range(19)] for i in range(19)]
)
require(H18 * H18_inverse == sp.eye(19), "degree-18 inverse formula failed")
recovered = []
for j in (6, 18):
    phi = sum(H18_inverse[i, j] * s**i for i in range(19))
    recovered.append(factorial_functional(G4 * phi))
expected_recovered = (
    sp.Poly(G4, s).coeff_monomial(s**6),
    sp.Poly(G4, s).coeff_monomial(s**18),
)
require(tuple(recovered) == expected_recovered, "private selector recovery failed")

print("THM-2815 OPTIMAL FINITE LAGUERRE CARRIER")
print(f"laguerre_D_range=1..{laguerre_checks[-1][0]}")
print(f"laguerre_last_first_failure={laguerre_checks[-1][1]}")
print("carrier_exact_degree=2D-1; minimum_dimension=D; unique_relation=ell_D")
print(f"selector_d_range=0..{selector_checks[-1]}")
print(f"mistake211_scalar={factorial_functional(G4)}")
print(f"mistake211_private_coefficients={expected_recovered}")
print(f"mistake211_selector_recovery={tuple(recovered)}")
print(
    "scope=finite-height optimal carrier and extra multiplier selectors; "
    "scalar moment nullity alone does not provide selector access"
)
print("ALL EXACT CHECKS PASSED")
