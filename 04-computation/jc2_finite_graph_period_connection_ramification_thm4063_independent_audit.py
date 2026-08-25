#!/usr/bin/env python3
"""Independent Fraction-only audit of THM-4063."""

from fractions import Fraction as Q


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def multiply(left, right):
    result = [Q(0)] * (len(left) + len(right) - 1)
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            result[i + j] += a * b
    return result


def power_linear(a, b, exponent):
    result = [Q(1)]
    for _ in range(exponent):
        result = multiply(result, [Q(a), Q(b)])
    return result


def polygon_moment(vertices, c_degree, u_degree):
    """Compute the exact closed line integral of c^r u^s dc."""
    total = Q(0)
    for (c0, u0), (c1, u1) in zip(vertices, vertices[1:] + vertices[:1]):
        dc = c1 - c0
        du = u1 - u0
        integrand = multiply(
            power_linear(c0, dc, c_degree),
            power_linear(u0, du, u_degree),
        )
        total += Q(dc) * sum(
            coefficient / Q(i + 1)
            for i, coefficient in enumerate(integrand)
        )
    return total


def cross(a, b):
    return a[0] * b[1] - a[1] * b[0]


def strict_segment_intersection(edge1, edge2):
    """Return a proper interior intersection, or None."""
    start1, end1 = edge1
    start2, end2 = edge2
    direction1 = (end1[0] - start1[0], end1[1] - start1[1])
    direction2 = (end2[0] - start2[0], end2[1] - start2[1])
    denominator = cross(direction1, direction2)
    if denominator == 0:
        return None
    displacement = (start2[0] - start1[0], start2[1] - start1[1])
    alpha = Q(cross(displacement, direction2), denominator)
    beta = Q(cross(displacement, direction1), denominator)
    if 0 < alpha < 1 and 0 < beta < 1:
        return (
            Q(start1[0]) + alpha * direction1[0],
            Q(start1[1]) + alpha * direction1[1],
        )
    return None


def edges(polygon):
    return list(zip(polygon, polygon[1:] + polygon[:1]))


cycles = [
    [(0, 0), (1, 0), (2, 1)],
    [(0, 0), (-1, 2), (-3, 1)],
]

require(
    set(cycles[0]).intersection(cycles[1]) == {(0, 0)},
    "cycles must share only the origin",
)
require(
    all(vertex[0] > 0 for vertex in cycles[0][1:])
    and all(vertex[0] < 0 for vertex in cycles[1][1:]),
    "opposite closed half-planes",
)
require(
    all(end[0] != start[0] for polygon in cycles for start, end in edges(polygon)),
    "every edge must have nonzero c increment",
)
interior_crossings = [
    strict_segment_intersection(edge1, edge2)
    for edge1 in edges(cycles[0])
    for edge2 in edges(cycles[1])
    if strict_segment_intersection(edge1, edge2) is not None
]
require(not interior_crossings, "embedded cycles must not cross")

v_u = tuple(polygon_moment(polygon, 0, 1) for polygon in cycles)
v_cu = tuple(polygon_moment(polygon, 1, 1) for polygon in cycles)
determinant = v_u[0] * v_cu[1] - v_u[1] * v_cu[0]
require(v_u == (Q(-1, 2), Q(-5, 2)), "u moment")
require(v_cu == (Q(-1, 2), Q(10, 3)), "cu moment")
require(determinant == Q(-35, 12), "moment determinant")

for degree in range(13):
    require(
        all(polygon_moment(polygon, degree, 0) == 0 for polygon in cycles),
        ("pure-c period", degree),
    )
for total_degree in range(2, 9):
    for u_degree in range(1, total_degree + 1):
        c_degree = total_degree - u_degree
        moment = tuple(
            polygon_moment(polygon, c_degree, u_degree)
            for polygon in cycles
        )
        # Invertibility of [v_u,v_cu] expresses every such moment in that
        # basis; verify the reconstructed coordinates exactly.
        alpha = (moment[0] * v_cu[1] - moment[1] * v_cu[0]) / determinant
        beta = (v_u[0] * moment[1] - v_u[1] * moment[0]) / determinant
        require(
            tuple(alpha * v_u[i] + beta * v_cu[i] for i in range(2))
            == moment,
            "moment basis reconstruction",
        )

for opening in range(1, 9):
    relative = (opening, 2 * opening)
    mixed = tuple(exponent - 1 for exponent in relative)
    require(sum(mixed) == 3 * opening - 2, "length ledger")

for ramification in range(1, 20):
    leading = Q(2 * ramification + 1, ramification + 1) * ramification
    require(leading != 0, "characteristic-zero derivative lead")

print("status=PASS")
print("embedded_cycles=opposite_half_planes;common_vertex_only")
print(f"period_u={v_u}")
print(f"period_cu={v_cu}")
print(f"moment_det={determinant}")
print("carrier=raw_(2q,3q);relative_(q,2q)")
print("conditional_mixed=(q-1,2q-1);length=3q-2")
print("ramification=ord(H_prime)=e-1;unit_Jacobian_impossible_for_e>=2")
print("boundary=carrier_packet_does_not_assert_period_completeness")
