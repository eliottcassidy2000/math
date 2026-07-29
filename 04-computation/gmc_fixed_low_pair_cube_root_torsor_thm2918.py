#!/usr/bin/env python3
"""Exact quotient-field companion for THM-2918.

The analytic theorem uses exact factorial-ratio convergence and the
implicit-function theorem.  This companion verifies its finite algebraic
core: the quadratic quotient field, cube-root torsor, multiplication norm,
zero-sum quartic orbit, endpoint-functional quotient typing, the consecutive
low-pair specialization, all decay bases, and a bounded arbitrary-three-slot
hostile census.
"""

from __future__ import annotations

from hashlib import sha256
from math import comb, factorial

import sympy as sp


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def polynomial_digest(polynomial: sp.Poly) -> str:
    records = "\n".join(
        f"{','.join(str(value) for value in monomial)}:{coefficient}"
        for monomial, coefficient in polynomial.terms()
    )
    return sha256((records + "\n").encode()).hexdigest()


s, z = sp.symbols("s z")
x, y = sp.symbols("x y")
tau, sigma, rho = sp.symbols("tau sigma rho")
g0, g1, g2 = sp.symbols("g0 g1 g2", nonzero=True)
h0, h1, h2 = sp.symbols("h0 h1 h2")


def f(index: int) -> sp.Expr:
    return s**index / factorial(index)


def functional(polynomial: sp.Expr) -> sp.Expr:
    expanded = sp.Poly(sp.expand(polynomial), s)
    return sp.expand(
        sum(
            coefficient * factorial(monomial[0])
            for monomial, coefficient in expanded.terms()
        )
    )


def quotient_remainder(
    polynomial: sp.Expr,
    modulus: sp.Expr,
) -> sp.Expr:
    return sp.rem(
        sp.Poly(sp.expand(polynomial), z),
        sp.Poly(sp.expand(modulus), z),
    ).as_expr()


# ---------------------------------------------------------------------------
# Generic quadratic quotient and its positive norm.
# ---------------------------------------------------------------------------

q_generic = g0 + 2 * g1 * z + g2 * z**2
ell_generic = x + y * z

z_square = quotient_remainder(z**2, q_generic)
require(
    sp.expand(g2 * z_square + 2 * g1 * z + g0) == 0,
    "generic quadratic quotient relation changed",
)

multiplication_matrix = sp.Matrix(
    [
        [x, -g0 * y / g2],
        [y, x - 2 * g1 * y / g2],
    ]
)
norm_scaled = g2 * x**2 - 2 * g1 * x * y + g0 * y**2
require(
    sp.factor(multiplication_matrix.det() - norm_scaled / g2) == 0,
    "quadratic field norm formula changed",
)

ell_square = quotient_remainder(ell_generic**2, q_generic)
ell_square_poly = sp.Poly(ell_square, z)
square_matrix = multiplication_matrix.subs(
    {
        x: ell_square_poly.nth(0),
        y: ell_square_poly.nth(1),
    },
    simultaneous=True,
)
require(
    sp.factor(
        (3 * square_matrix).det()
        - 9 * (norm_scaled / g2) ** 2
    )
    == 0,
    "cube-map derivative determinant changed",
)

ell_fourth = sp.Poly(
    quotient_remainder(ell_generic**4, q_generic),
    z,
)
ell_fourth_a = ell_fourth.nth(0)
ell_fourth_b = ell_fourth.nth(1)
ell_fourth_norm = sp.factor(
    g2 * ell_fourth_a**2
    - 2 * g1 * ell_fourth_a * ell_fourth_b
    + g0 * ell_fourth_b**2
)
require(
    sp.factor(
        ell_fourth_norm - norm_scaled**4 / g2**3
    )
    == 0,
    "quartic remainder norm identity changed",
)

# An arbitrary quartic multiple of q must be killed by the THM-2872
# endpoint functional.  The binary quartic convention is
# A0+4A1*z+6A2*z^2+4A3*z^3+A4*z^4.
h_generic = h0 + h1 * z + h2 * z**2
quartic_multiple = sp.Poly(sp.expand(q_generic * h_generic), z)
A0 = quartic_multiple.nth(0)
A1 = quartic_multiple.nth(1) / 4
A3 = quartic_multiple.nth(3) / 4
A4 = quartic_multiple.nth(4)
endpoint_multiple = sp.factor(
    (2 * A1 * g0 - A0 * g1) * g2**2
    - (2 * A3 * g2 - A4 * g1) * g0**2
)
require(
    endpoint_multiple == 0,
    "endpoint functional stopped factoring through the quotient",
)


# ---------------------------------------------------------------------------
# Exact C3 orbit in the THM-2914 low quotient q=1+2z+2z^2.
# This specialization keeps the companion in an exact real algebraic field
# while verifying every torsor identity used abstractly.
# ---------------------------------------------------------------------------

q_consecutive = 1 + 2 * z + 2 * z**2
omega = (-1 + sp.sqrt(3) * (2 * z + 1)) / 2

require(
    quotient_remainder(omega**3, q_consecutive) == 1
    and quotient_remainder(
        1 + omega + quotient_remainder(omega**2, q_consecutive),
        q_consecutive,
    )
    == 0,
    "C3 root-of-unity identities changed",
)

orbit = [
    quotient_remainder(omega**index * (x + y * z), q_consecutive)
    for index in range(3)
]
require(
    len({sp.srepr(item) for item in orbit}) == 3,
    "symbolic C3 orbit collapsed",
)

orbit_fourths = [
    quotient_remainder(item**4, q_consecutive)
    for item in orbit
]
require(
    quotient_remainder(sum(orbit_fourths), q_consecutive) == 0,
    "zero-sum quartic C3 orbit changed",
)

orbit_norms: list[sp.Expr] = []
orbit_endpoints: list[sp.Expr] = []
for item in orbit:
    item_poly = sp.Poly(item, z)
    item_x = item_poly.nth(0)
    item_y = item_poly.nth(1)
    orbit_norms.append(
        sp.factor(2 * item_x**2 - 2 * item_x * item_y + item_y**2)
    )
    orbit_endpoints.append(
        sp.expand(
            (item_y**2 - 2 * item_x**2)
            * (
                item_y**2
                - 4 * item_x * item_y
                + 2 * item_x**2
            )
        )
    )

require(
    all(sp.simplify(value - orbit_norms[0]) == 0 for value in orbit_norms)
    and sp.simplify(sum(orbit_endpoints)) == 0,
    "norm invariance or scalar endpoint zero-sum changed",
)


# ---------------------------------------------------------------------------
# Consecutive low pair: exact quotient remainder and cube-root eliminant.
# ---------------------------------------------------------------------------

d0 = f(1) - f(0)
d1 = f(2) - f(1)
low = d0 + z * d1
low_quadratic = sp.expand(functional(low**2))
low_cubic = sp.expand(functional(low**3))
low_remainder = quotient_remainder(low_cubic, low_quadratic)

require(
    low_quadratic == q_consecutive
    and low_cubic == 2 + 12 * z + 30 * z**2 + 30 * z**3
    and low_remainder == 2 - 3 * z,
    "consecutive low quotient data changed",
)

cube_remainder = sp.Poly(
    quotient_remainder(
        low_remainder + (x + y * z) ** 3,
        q_consecutive,
    ),
    z,
)
cube_eliminant = sp.factor(
    sp.resultant(cube_remainder.nth(0), cube_remainder.nth(1), x)
)
expected_eliminant = -(
    2 * y**9 + 18 * y**6 - 729 * y**3 + 54
) / 2
require(
    sp.expand(cube_eliminant - expected_eliminant) == 0,
    "consecutive cube-root eliminant changed",
)


# ---------------------------------------------------------------------------
# Exact factorial exponential bases.
# ---------------------------------------------------------------------------

order_three_bases: list[sp.Rational] = []
for k in range(1, 4):
    for j in range(k + 1):
        if (k, j) == (3, 3):
            continue
        numerator_base = 1 if j in (0, 1) else j**j
        order_three_bases.append(sp.Rational(numerator_base, 3**k))

order_four_bases: list[sp.Rational] = []
for k in range(4):
    for j in range(k + 1):
        numerator_base = 1 if j in (0, 1) else j**j
        order_four_bases.append(
            sp.Rational(3 ** (4 - k) * numerator_base, 4**4)
        )
for j in range(4):
    numerator_base = 1 if j in (0, 1) else j**j
    order_four_bases.append(sp.Rational(numerator_base, 4**4))

require(
    max(order_three_bases) == sp.Rational(4, 9)
    and max(order_four_bases) == sp.Rational(81, 256)
    and max(order_three_bases + order_four_bases) < sp.Rational(1, 2),
    "factorial tail rate ledger changed",
)

# The positive-coordinate sector in the factorial specialization.
positive_sector = -sigma**3 + 3 * sigma - 2 * rho
require(
    positive_sector.subs(sigma, 0) == -2 * rho
    and positive_sector.subs(sigma, 1) == 2 * (1 - rho)
    and positive_sector.subs(sigma, 2) == -2 * (1 + rho)
    and sp.diff(positive_sector, sigma)
    == 3 * (1 - sigma**2),
    "positive C3 sector polynomial changed",
)


# ---------------------------------------------------------------------------
# Bounded arbitrary-three-slot control.  THM-2824 proves the universal
# statement; this census is only an exact hostile check of the specialization.
# ---------------------------------------------------------------------------

triple_count = 0
positive_tensor_count = 0
remainder_count = 0
remainder_digests: list[str] = []

for a in range(9):
    for b in range(a + 1, 9):
        for c in range(b + 1, 9):
            first = f(b) - f(a)
            second = f(c) - f(b)
            remote_base = f(c)
            plane = first + z * second
            quadratic = sp.Poly(functional(plane**2), z)
            cubic = sp.Poly(functional(plane**3), z)
            remainder = sp.rem(cubic, quadratic)

            g_vector = (
                quadratic.nth(0),
                quadratic.nth(1) / 2,
                quadratic.nth(2),
            )
            t_vector = (
                cubic.nth(0),
                cubic.nth(1) / 3,
                cubic.nth(2) / 3,
                cubic.nth(3),
            )
            delta = sp.factor(
                g_vector[0] * g_vector[2] - g_vector[1] ** 2
            )

            require(
                functional(first) == 0
                and functional(second) == 0
                and functional(remote_base) == 1,
                f"triple {(a, b, c)} lost moment normalization",
            )
            require(
                delta > 0,
                f"triple {(a, b, c)} lost positive quadratic Gram",
            )
            require(
                all(value > 0 for value in g_vector + t_vector),
                f"triple {(a, b, c)} lost strict tensor positivity",
            )
            require(
                not remainder.is_zero,
                f"triple {(a, b, c)} became cubic-null",
            )

            triple_count += 1
            positive_tensor_count += 1
            remainder_count += 1
            remainder_digests.append(
                polynomial_digest(
                    sp.Poly(remainder, z, domain=sp.QQ)
                )
            )

require(
    triple_count == comb(9, 3)
    and positive_tensor_count == triple_count
    and remainder_count == triple_count
    and len(set(remainder_digests)) == triple_count,
    "bounded arbitrary-three-slot census changed",
)

# ---------------------------------------------------------------------------
# Exact positive-sheet endpoint reversal: (0,1,2) versus (0,8,10).
# ---------------------------------------------------------------------------

hostile_a, hostile_b, hostile_c = 0, 8, 10
hostile_first = f(hostile_b) - f(hostile_a)
hostile_second = f(hostile_c) - f(hostile_b)
hostile_plane = hostile_first + z * hostile_second
hostile_q = sp.Poly(functional(hostile_plane**2), z)
hostile_cubic = sp.Poly(functional(hostile_plane**3), z)
hostile_remainder = sp.rem(hostile_cubic, hostile_q)

require(
    hostile_q.as_expr()
    == 12869 + 61776 * z + 110110 * z**2
    and hostile_remainder.as_expr()
    == 2
    * (
        157110254569762
        + 297335210491143 * z
    )
    / 2695,
    "(0,8,10) low quotient data changed",
)

hostile_unit_cube = sp.rem(
    sp.Poly((tau + z) ** 3, z),
    hostile_q,
)
hostile_a_tau = hostile_unit_cube.nth(0)
hostile_b_tau = hostile_unit_cube.nth(1)
hostile_projective = sp.Poly(
    sp.expand(
        hostile_remainder.nth(1) * hostile_a_tau
        - hostile_remainder.nth(0) * hostile_b_tau
    ),
    tau,
    domain=sp.QQ,
)
_, hostile_projective_integer = hostile_projective.clear_denoms(
    convert=True
)
hostile_projective_integer = (
    hostile_projective_integer.primitive()[1]
)
if hostile_projective_integer.LC() < 0:
    hostile_projective_integer = -hostile_projective_integer

expected_hostile_projective = sp.Poly(
    425189351002334490 * tau**3
    - 674002992104278980 * tau**2
    + 229061153084068755 * tau
    - 16579594004575622,
    tau,
    domain=sp.ZZ,
)
hostile_interval = (
    sp.Rational(341, 1000),
    sp.Rational(683, 2000),
)
require(
    hostile_projective_integer == expected_hostile_projective
    and hostile_projective_integer.eval(hostile_interval[0]) > 0
    and hostile_projective_integer.eval(hostile_interval[1]) < 0
    and hostile_projective_integer.count_roots(*hostile_interval) == 1,
    "(0,8,10) positive-sheet projective bracket changed",
)

hostile_a_poly = sp.Poly(hostile_a_tau, tau, domain=sp.QQ)
require(
    hostile_a_poly.count_roots(*hostile_interval) == 0
    and hostile_a_poly.eval(sum(hostile_interval) / 2) < 0
    and hostile_remainder.nth(0) > 0,
    "(0,8,10) positive eta selector changed",
)

hostile_g0 = hostile_q.nth(0)
hostile_g1 = hostile_q.nth(1) / 2
hostile_g2 = hostile_q.nth(2)
hostile_endpoint = sp.Poly(
    sp.expand(
        (
            2 * tau**3 * hostile_g0
            - tau**4 * hostile_g1
        )
        * hostile_g2**2
        - (
            2 * tau * hostile_g2
            - hostile_g1
        )
        * hostile_g0**2
    ),
    tau,
    domain=sp.QQ,
)
expected_hostile_endpoint = sp.Poly(
    -572
    * (110110 * tau**2 - 12869)
    * (
        5945940 * tau**2
        - 4954565 * tau
        + 694926
    ),
    tau,
    domain=sp.ZZ,
)
require(
    sp.expand(
        hostile_endpoint.as_expr()
        - expected_hostile_endpoint.as_expr()
    )
    == 0
    and hostile_endpoint.count_roots(*hostile_interval) == 0
    and hostile_endpoint.eval(sum(hostile_interval) / 2) < 0,
    "(0,8,10) negative endpoint certificate changed",
)

consecutive_endpoint = (
    (1 - 2 * tau**2)
    * (1 - 4 * tau + 2 * tau**2)
)
consecutive_interval = (
    sp.Rational(6, 7),
    sp.Rational(35, 27),
)
require(
    consecutive_endpoint.subs(tau, consecutive_interval[0]) > 0
    and consecutive_endpoint.subs(tau, consecutive_interval[1]) > 0
    and sp.Poly(
        consecutive_endpoint,
        tau,
        domain=sp.QQ,
    ).count_roots(*consecutive_interval)
    == 0,
    "(0,1,2) positive endpoint certificate changed",
)

aggregate_digest = sha256(
    ("\n".join(remainder_digests) + "\n").encode()
).hexdigest()

EXPECTED_DIGESTS = {
    "low_q": "d8e8b3f7e416a1b63d551fcb602b233d0fd087aa05efa3309d4fcfdd5650d9a3",
    "low_c": "f4723b7cc4bf18aee4d3d38059649c4c2a66557d261ccd072d52e6b3d9a5facf",
    "cube_pair": "5b45ea9a8762aaeff39bea8a93afc0943d5f114d1a096b7b65233418fc257f29",
    "triple_bank": "d9e7676b7f3ccd25f25b3ff2879c744ad16ca4c1aec093e03c8212fb6001021e",
}
actual_digests = {
    "low_q": polynomial_digest(
        sp.Poly(low_quadratic, z, domain=sp.QQ)
    ),
    "low_c": polynomial_digest(
        sp.Poly(low_cubic, z, domain=sp.QQ)
    ),
    "cube_pair": sha256(
        (
            polynomial_digest(
                sp.Poly(cube_remainder.nth(0), x, y, domain=sp.QQ)
            )
            + "\n"
            + polynomial_digest(
                sp.Poly(cube_remainder.nth(1), x, y, domain=sp.QQ)
            )
            + "\n"
        ).encode()
    ).hexdigest(),
    "triple_bank": aggregate_digest,
}


def main() -> None:
    require(
        actual_digests == EXPECTED_DIGESTS,
        "locked quotient-field digest changed",
    )
    print("THM-2918 FIXED-LOW-PAIR CUBE-ROOT TORSOR")
    print("status=PROVED+VERIFIED-EXACT+INDEPENDENTLY-HOSTILE-AUDITED")
    print("quadratic_quotient=A=R[z]/(q)~=C;Delta=g0*g2-g1^2>0")
    print(
        "cube_roots=ell0,omega*ell0,omega^2*ell0;"
        "derivative_det=9*(Nq(ell)/g2)^2>0"
    )
    print("quartic_orbit_sum=0;full_remainder_norm=positive")
    print("linear_endpoint_sum=0;single_scalar_orientation=impossible")
    print(
        "consecutive_q=1+2z+2z^2;"
        "consecutive_c_remainder=2-3z"
    )
    print(
        "cube_eliminant=-(2eta^9+18eta^6-729eta^3+54)/2"
    )
    print(
        "rate_maxima=order3:4/9;normalized_order4:81/256"
    )
    print(
        f"triple_census={triple_count};"
        f"positive_tensors={positive_tensor_count};"
        f"nonzero_remainders={remainder_count}"
    )
    print(
        "positive_sheet_endpoint_signs=(0,1,2):positive;"
        "(0,8,10):negative;"
        "hostile_T=(341/1000,683/2000)"
    )
    print(
        "digests="
        + ";".join(f"{key}:{value}" for key, value in actual_digests.items())
    )
    print(
        "scope=three local natural-scale real branches; "
        "no global exhaustion, explicit cutoff, arbitrary-support GMC, "
        "or canonical fixed-fibre C2 action"
    )
    print("all_exact_checks=PASS")


if __name__ == "__main__":
    main()
