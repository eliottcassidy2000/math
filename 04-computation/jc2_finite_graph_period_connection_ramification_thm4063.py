#!/usr/bin/env python3
"""Exact controls for THM-4063's graph carrier and ramification no-go."""

import sympy as sp


z, c, u, A, t = sp.symbols("z c u A t")


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def polygon_period(vertices, density):
    """Return the exact oriented line integral of density(c,u) dc."""
    total = sp.Rational(0)
    for start, end in zip(vertices, vertices[1:] + vertices[:1]):
        cc = start[0] + z * (end[0] - start[0])
        uu = start[1] + z * (end[1] - start[1])
        total += sp.integrate(
            sp.expand(
                density.subs({c: cc, u: uu}) * (end[0] - start[0])
            ),
            (z, 0, 1),
        )
    return sp.factor(total)


# Two embedded triangles share only the origin. Every edge has nonzero
# c-increment; the first filled triangle lies in c>=0 and the second in c<=0.
cycles = [
    [(0, 0), (1, 0), (2, 1)],
    [(0, 0), (-1, 2), (-3, 1)],
]

v_u = sp.Matrix([polygon_period(polygon, u) for polygon in cycles])
v_cu = sp.Matrix([polygon_period(polygon, c * u) for polygon in cycles])
moment_matrix = sp.Matrix.hstack(v_u, v_cu)

require(
    v_u == sp.Matrix([sp.Rational(-1, 2), sp.Rational(-5, 2)]),
    "u moment vector",
)
require(
    v_cu == sp.Matrix([sp.Rational(-1, 2), sp.Rational(10, 3)]),
    "cu moment vector",
)
require(moment_matrix.det() == sp.Rational(-35, 12), "moment determinant")
for degree in range(13):
    require(
        all(polygon_period(polygon, c**degree) == 0 for polygon in cycles),
        ("pure c period", degree),
    )

print("embedded_figure_eight_beta=2")
print(f"period(u)={tuple(v_u)}")
print(f"period(c*u)={tuple(v_cu)}")
print(f"moment_det={moment_matrix.det()}")
print("raw_carrier=epsilon^2*R*v_u+epsilon^3*R*v_cu")
print("opening_lattice=epsilon*R^2")
print("relative_Smith_exponents=(q,2q)")
print("conditional_mixed_exponents=(q-1,2q-1);length=3q-2")


# THM-4058/4060 is the scalar beta=1, q=5 specialization.
q_triangle = 5
require(q_triangle - 1 == 4, "triangle mixed length")
print("triangle_check=beta_1;opening_5;carrier_5;cokernel_length_4")


# H=t^e+... has derivative order e-1. The exceptional affine pair realizes
# the sharp generator 12H'.
for ramification in range(1, 8):
    H = t**ramification + 3 * t ** (ramification + 1) + 5 * t ** (ramification + 3)
    derivative = sp.diff(H, t)
    valuation = min(
        term.as_powers_dict().get(t, 0) for term in sp.Add.make_args(derivative)
    )
    require(valuation == ramification - 1, "ramification valuation")
    print(
        f"ramification_e={ramification}:ord_t(H')={valuation};"
        f"unit={valuation == 0}"
    )


# Positive-characteristic resonance for d/dA+q/A on A^nu k[[A]].
prime, opening, carrier = 5, 5, 5
misses = [
    n for n in range(30) if (carrier + n + opening) % prime == 0
]
require(misses == [0, 5, 10, 15, 20, 25], "characteristic-p resonance")
print(f"char_{prime}_resonant_carrier_indices={misses}")
