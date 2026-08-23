#!/usr/bin/env python3
"""Exact additive edge-chart cocycle for the THM-3466 primitive.

The control is the identity Keller map g=x+iy on a rational, centred,
positively oriented triangle.  It tests only the single closure moment needed
for g^2 d(conjugate g), not the full HFC null orbit.
"""

from fractions import Fraction
import hashlib


ZERO = (Fraction(0), Fraction(0))
ONE_THIRD = Fraction(1, 3)


def add(left, right):
    return (left[0] + right[0], left[1] + right[1])


def neg(value):
    return (-value[0], -value[1])


def sub(left, right):
    return add(left, neg(right))


def mul(left, right):
    return (
        left[0] * right[0] - left[1] * right[1],
        left[0] * right[1] + left[1] * right[0],
    )


def scale(value, scalar):
    return (scalar * value[0], scalar * value[1])


def conjugate(value):
    return (value[0], -value[1])


def edge_integral(start, end):
    """Integral z(t)^2 d(conjugate z(t)) on one affine edge."""
    direction = sub(end, start)
    integrated_square = add(
        add(mul(start, start), mul(start, direction)),
        scale(mul(direction, direction), ONE_THIRD),
    )
    return mul(conjugate(direction), integrated_square)


def digest(value):
    return hashlib.sha256(repr(value).encode("ascii")).hexdigest()


def main():
    vertices = (
        (Fraction(-1), Fraction(-1)),
        (Fraction(2), Fraction(-1)),
        (Fraction(-1), Fraction(2)),
    )
    edges = tuple(
        (vertices[index], vertices[(index + 1) % len(vertices)])
        for index in range(len(vertices))
    )
    increments = tuple(edge_integral(*edge) for edge in edges)
    holonomy = ZERO
    for increment in increments:
        holonomy = add(holonomy, increment)
    if holonomy != ZERO:
        raise RuntimeError(("primitive did not close", holonomy))

    # Offsets turn edgewise primitives, each normalized to zero at its own
    # initial vertex, into the primitive based at vertex zero.
    offsets = [ZERO]
    for increment in increments[:-1]:
        offsets.append(add(offsets[-1], increment))
    closure = add(offsets[-1], increments[-1])
    if closure != ZERO:
        raise RuntimeError(("offset clutch did not close", closure))

    # The local transition at the common endpoint of edge e and e+1 is the
    # translation needed to identify their independently normalized values.
    transition_constants = tuple(
        increments[index] for index in range(len(increments))
    )

    print("vertices=" + repr(vertices))
    print("edge_primitive_increments=" + repr(increments))
    print("transition_constants=" + repr(transition_constants))
    print("global_basepoint_offsets=" + repr(tuple(offsets)))
    print("cocycle_holonomy=" + repr(holonomy))
    print("record_digest=" + digest((
        vertices, increments, tuple(offsets), holonomy
    )))
    print("typed_scope=identity_constant_J_map;centred_rational_triangle;"
          "single_g_squared_boundary_current_closure;not_full_HFC")
    print("all_exact_checks=PASS")


if __name__ == "__main__":
    main()
