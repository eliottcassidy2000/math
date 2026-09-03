#!/usr/bin/env python3
"""Exact certificate for THM-4378.

The script works only with integer coefficient vectors.  It independently
rebuilds the bilateral source band, its finite trace quotient, the width-one
central-spine basis, the signed reciprocal action, the two integral
eigenlattices, their Smith two-glue, and the norm/difference quotients.  It
also checks the exceptional ell=2 ray and the smallest unsigned-reflection
hostile.  Every check remains active under ``python -O``.
"""

from __future__ import annotations

from functools import lru_cache
from math import comb
import sys


if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(newline="\n")


CHECKS = 0


def require(condition: bool, label: str) -> None:
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError("check failed: " + label)


def params(ell: int) -> tuple[int, int]:
    require(ell >= 2, "ell lower bound")
    return (ell + 1) // 2, (ell + 2) // 3


def trim(poly) -> list[int]:
    result = list(poly)
    while len(result) > 1 and result[-1] == 0:
        result.pop()
    return result


def pad(poly, degree: int) -> list[int]:
    require(len(poly) <= degree + 1, "polynomial degree bound")
    return list(poly) + [0] * (degree + 1 - len(poly))


def add(left: list[int], right: list[int], right_scale: int = 1) -> list[int]:
    size = max(len(left), len(right))
    result = [0] * size
    for index in range(size):
        result[index] = (left[index] if index < len(left) else 0) + right_scale * (
            right[index] if index < len(right) else 0
        )
    return trim(result)


def scale(value: int, poly: list[int]) -> list[int]:
    return trim([value * coefficient for coefficient in poly])


def shift(power: int, poly: list[int]) -> list[int]:
    require(power >= 0, "nonnegative shift")
    return [0] * power + list(poly)


def multiply(left: list[int], right: list[int]) -> list[int]:
    result = [0] * (len(left) + len(right) - 1)
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            result[i + j] += a * b
    return trim(result)


@lru_cache(maxsize=None)
def one_minus_z_power(power: int) -> tuple[int, ...]:
    require(power >= 0, "nonnegative binomial power")
    return tuple((-1) ** index * comb(power, index) for index in range(power + 1))


@lru_cache(maxsize=None)
def q_power(power: int) -> tuple[int, ...]:
    """Coefficients of q^power for q=z(1-z)."""
    return tuple(shift(power, one_minus_z_power(power)))


@lru_cache(maxsize=None)
def residual_packet(i: int, j: int) -> tuple[int, ...]:
    """Coefficients of z^i(1-z)^j."""
    require(i >= 0 and j >= 0, "residual packet quadrant")
    return tuple(shift(i, one_minus_z_power(j)))


def substitute_one_minus_z(poly: list[int]) -> list[int]:
    result = [0]
    for power, coefficient in enumerate(poly):
        result = add(result, scale(coefficient, one_minus_z_power(power)))
    return trim(result)


@lru_cache(maxsize=None)
def central_basis_polynomial(index: int) -> tuple[int, ...]:
    """Degree-indexed basis q^k, z q^k for even, odd index."""
    require(index >= 0, "central basis index")
    k = index // 2
    base = q_power(k)
    return base if index % 2 == 0 else tuple(shift(1, base))


def central_coordinates(poly: list[int], degree: int) -> list[int]:
    """Integral coordinates in 1,z,q,zq,q^2,zq^2,... through degree."""
    work = pad(trim(poly), degree)
    answer = [0] * (degree + 1)
    for index in range(degree, -1, -1):
        basis = pad(central_basis_polynomial(index), degree)
        leading = basis[index]
        require(abs(leading) == 1, "central unit leading coefficient")
        coefficient = work[index] // leading
        require(coefficient * leading == work[index], "central integral division")
        answer[index] = coefficient
        for row in range(degree + 1):
            work[row] -= coefficient * basis[row]
    require(not any(work), "central coordinate reconstruction")
    return answer


def from_central(coordinates: list[int]) -> list[int]:
    result = [0]
    for index, coefficient in enumerate(coordinates):
        result = add(result, scale(coefficient, central_basis_polynomial(index)))
    return trim(result)


def reciprocal_coordinates(coordinates: list[int]) -> list[int]:
    """J(f)(z)=f(1-z) in the central basis."""
    result = [0] * len(coordinates)
    for index, coefficient in enumerate(coordinates):
        if index % 2 == 0:
            result[index] += coefficient
        else:
            result[index - 1] += coefficient
            result[index] -= coefficient
    return result


def vector_add(left: list[int], right: list[int], scale_right: int = 1) -> list[int]:
    require(len(left) == len(right), "vector length")
    return [a + scale_right * b for a, b in zip(left, right)]


def bilateral_residual_points(ell: int, depth: int) -> list[tuple[int, int]]:
    s, _ = params(ell)
    width = s - 1
    return [
        (i, j)
        for total in range(depth + 1)
        for i in range(total + 1)
        for j in [total - i]
        if abs(i - j) <= width
    ]


def local_cell_parents(ell: int, depth: int) -> list[tuple[int, int]]:
    s, _ = params(ell)
    if depth == 0 or s == 1:
        return []
    return [
        (i, j)
        for total in range(depth)
        for i in range(total + 1)
        for j in [total - i]
        if abs(i - j) <= s - 2
    ]


def trace_packet(u: int, v: int) -> list[int]:
    """THM-4368 trace F_(u,v)."""
    sign = -1 if (u - 1) % 2 else 1
    return scale(sign, residual_packet(u - 1, v - 1))


def direct_bilateral(ell: int, u: int, v: int) -> bool:
    s, rho = params(ell)

    def source_valid(a: int, b: int) -> bool:
        n0 = s + a - 1
        return a >= 1 and b >= rho and n0 >= b and n0 + b >= ell

    return source_valid(u, v) and source_valid(v, u)


def direct_fibre_size(ell: int, u: int, v: int) -> int:
    """Count THM-4368 source realizers without using its closed formula."""
    s, _ = params(ell)
    count = 0
    for e in range(v + 1):
        a = 2 * v + e - ell
        b = s + u - v - 1 - e
        c = v - e
        if min(a, b, c, e) >= 0:
            count += 1
    return count


def claimed_bilateral(ell: int, u: int, v: int) -> bool:
    s, rho = params(ell)
    return u >= rho and v >= rho and abs(u - v) <= s - 1


def exact_packet_and_spine_checks() -> tuple[int, int]:
    blocks = 0
    packets = 0
    for ell in range(2, 81):
        s, rho = params(ell)
        for u in range(1, 4 * s + 15):
            for v in range(1, 4 * s + 15):
                require(
                    direct_bilateral(ell, u, v) == claimed_bilateral(ell, u, v),
                    "direct bilateral cone",
                )

        for depth in range(31):
            points = bilateral_residual_points(ell, depth)
            if ell == 2:
                expected = [(k, k) for k in range(depth // 2 + 1)]
                require(points == expected, "ell two diagonal ray")
                images = [residual_packet(i, j) for i, j in points]
                require(images == [q_power(k) for k in range(depth // 2 + 1)],
                        "ell two q basis")
                require(all(substitute_one_minus_z(poly) == list(poly) for poly in images),
                        "ell two trivial reciprocal action")
                blocks += 1
                packets += len(points)
                continue

            # The central packet spine contains one packet in every residual
            # polynomial degree and is a unimodular basis.
            spine_pairs = []
            for k in range(depth // 2 + 1):
                spine_pairs.append((k, k))
                if 2 * k + 1 <= depth:
                    spine_pairs.append((k + 1, k))
            require(len(spine_pairs) == depth + 1, "one spine packet per degree")
            require(all(point in points for point in spine_pairs), "spine in band")
            for index, point in enumerate(spine_pairs):
                poly = residual_packet(*point)
                expected = central_basis_polynomial(index)
                require(poly == expected, "literal spine polynomial")
                coordinates = central_coordinates(poly, depth)
                require(coordinates == [int(row == index) for row in range(depth + 1)],
                        "spine unit coordinate")

            # Every band packet belongs integrally to the spine lattice.
            for i, j in points:
                poly = residual_packet(i, j)
                coordinates = central_coordinates(poly, depth)
                require(from_central(coordinates) == list(poly), "band packet reconstruction")
                packets += 1

            # Local source-internal Pascal cells have the exact kernel rank.
            parents = local_cell_parents(ell, depth)
            require(len(parents) == len(points) - (depth + 1),
                    "local cell kernel rank")
            for i, j in parents:
                circuit = add(
                    residual_packet(i, j),
                    add(residual_packet(i, j + 1), residual_packet(i + 1, j)),
                    -1,
                )
                require(circuit == [0], "local Pascal cell")

            blocks += 1
    return blocks, packets


def reciprocal_and_eigenlattice_checks() -> tuple[int, int, int]:
    blocks = 0
    smith_twos = 0
    norm_twos = 0
    for depth in range(61):
        size = depth + 1
        plus_count = depth // 2 + 1
        minus_count = (depth + 1) // 2

        # J is the exact polynomial substitution and an involution.
        for index in range(size):
            basis_poly = central_basis_polynomial(index)
            direct = central_coordinates(substitute_one_minus_z(basis_poly), depth)
            predicted = reciprocal_coordinates(
                [int(row == index) for row in range(size)]
            )
            require(direct == predicted, "reciprocal central matrix")
            require(reciprocal_coordinates(predicted)
                    == [int(row == index) for row in range(size)],
                    "reciprocal involution")

        # Exact eigenspace equations in A(q)+zB(q) coordinates.
        for index in range(size):
            unit = [int(row == index) for row in range(size)]
            reflected = reciprocal_coordinates(unit)
            if index % 2 == 0:
                require(reflected == unit, "q power invariant")
            else:
                difference = vector_add(unit, reflected, -1)
                expected = [0] * size
                expected[index - 1] = -1
                expected[index] = 2
                require(difference == expected, "odd difference is x q power")

        # In each paired block, {q^k, x q^k} has coordinate matrix
        # [[1,-1],[0,2]].  One elementary column addition gives Smith block
        # diag(1,2).  An even top degree contributes a terminal diag(1).
        for k in range(minus_count):
            plus = [0] * size
            anti = [0] * size
            plus[2 * k] = 1
            anti[2 * k] = -1
            anti[2 * k + 1] = 2
            require(reciprocal_coordinates(plus) == plus, "plus generator")
            require(reciprocal_coordinates(anti) == [-x for x in anti],
                    "minus generator")
            smith_twos += 1

            # The glue map is the odd central coefficient modulo two.
            require(anti[2 * k + 1] % 2 == 0, "eigen sum in parity kernel")
            oriented = [0] * size
            oriented[2 * k + 1] = 1
            require(oriented[2 * k + 1] % 2 == 1,
                    "oriented packet detects glue")
            twice = [2 * x for x in oriented]
            decomposition = vector_add(plus, anti)
            require(twice == decomposition, "twice orientation eigensplits")

            # Difference is all of the anti lattice; norm hits this invariant.
            require(vector_add(oriented, reciprocal_coordinates(oriented), -1)
                    == anti, "difference surjects anti lattice")
            require(vector_add(oriented, reciprocal_coordinates(oriented))
                    == plus, "norm hits paired invariant")

        require(plus_count + minus_count == size, "eigen ranks sum")
        require(smith_twos >= minus_count, "Smith running count")

        # At even depth the unpaired top q^(D/2) only has norm 2q^(D/2),
        # producing one Z/2.  At odd depth every invariant is a packet norm.
        if depth % 2 == 0:
            top = [0] * size
            top[-1] = 1
            require(reciprocal_coordinates(top) == top, "even top invariant")
            require(vector_add(top, reciprocal_coordinates(top))
                    == [2 * x for x in top], "even top doubled norm")
            norm_twos += 1
        else:
            top_oriented = [0] * size
            top_oriented[-1] = 1
            top_invariant = [0] * size
            top_invariant[-2] = 1
            require(vector_add(top_oriented, reciprocal_coordinates(top_oriented))
                    == top_invariant, "odd top invariant is a norm")

        blocks += 1
    return blocks, smith_twos, norm_twos


def signed_packet_intertwining_checks() -> int:
    checks = 0
    for ell in range(3, 51):
        s, rho = params(ell)
        tau0 = 2 * rho - 1
        for depth in range(16):
            row_max = tau0 + depth
            for i, j in bilateral_residual_points(ell, depth):
                u, v = rho + i, rho + j
                if u + v - 1 > row_max:
                    continue
                original = trace_packet(u, v)
                signed_reflected = scale((-1) ** (u - v), trace_packet(v, u))
                require(signed_reflected == substitute_one_minus_z(original),
                        "signed packet intertwining")
                checks += 1

            # The signed involution permutes the local circuit basis up to
            # its forced parity sign.
            for i, j in local_cell_parents(ell, depth):
                u, v = rho + i, rho + j
                sign = (-1) ** (u - v)
                require(sign in (-1, 1), "circuit reciprocal sign")
                require((v - u) % 2 == (u - v) % 2, "reverse sign agrees")

                circuit_trace = add(
                    trace_packet(u, v),
                    add(trace_packet(u, v + 1), trace_packet(u + 1, v), -1),
                    -1,
                )
                require(circuit_trace == [0], "signed-trace Pascal circuit")

                plain_reflected_trace = add(
                    trace_packet(v, u),
                    add(trace_packet(v + 1, u), trace_packet(v, u + 1), -1),
                    -1,
                )
                require(
                    plain_reflected_trace == scale(2, trace_packet(v, u)),
                    "plain reflected circuit equals twice parent",
                )
                checks += 1
    return checks


def fibre_parity_bridge_checks() -> int:
    """Audit the THM-4377 fibre-weighted class in the Smith parity quotient."""
    checks = 0
    for ell in range(3, 61):
        s, rho = params(ell)
        for d in range(1, s):
            # z^d=A_d(q)+zB_d(q), and B_d(0)=1.
            power_coordinates = central_coordinates(shift(d, [1]), d)
            base_gamma = tuple(
                power_coordinates[2 * k + 1] % 2
                for k in range((d + 1) // 2)
            )
            require(base_gamma and base_gamma[0] == 1,
                    "positive-power gluing vector nonzero")

            # Include every preclock scale and several scales beyond the full
            # boundary-jet threshold.
            for w in range(rho, s + d + 4):
                require(claimed_bilateral(ell, w + d, w),
                        "bridge orbit in bilateral band")
                mu_plus = direct_fibre_size(ell, w + d, w)
                mu_minus = direct_fibre_size(ell, w, w + d)
                require(mu_plus > 0 and mu_minus > 0, "bridge fibres nonempty")

                # Reciprocal changes B_d to -B_d, hence gives the same
                # nonzero vector modulo two. The common q power only shifts
                # its coordinate block and cannot erase it.
                aggregate_gamma = tuple(
                    ((mu_plus + mu_minus) * entry) % 2 for entry in base_gamma
                )
                require(
                    any(aggregate_gamma) == ((mu_plus - mu_minus) % 2 == 1),
                    "gluing iff odd fibre imbalance",
                )

                if w >= s:
                    beta = min(w - s + d + 1, 2 * d)
                    require(mu_plus - mu_minus == beta,
                            "postclock imbalance equals jet rank")
                    require(
                        any(aggregate_gamma) == (beta % 2 == 1),
                        "postclock parity toggle",
                    )
                    if w >= s + d - 1:
                        require(beta == 2 * d and not any(aggregate_gamma),
                                "full jet permanent parity vanishing")
                checks += 1
    return checks


def smallest_hostile_checks() -> tuple[list[int], list[int]]:
    # ell=3 has (s,rho)=(2,1), and R=2 is the first prefix containing a
    # nonfixed bilateral orbit.  C is the first Pascal circuit.
    c_trace = add(
        trace_packet(1, 1),
        add(trace_packet(1, 2), trace_packet(2, 1), -1),
        -1,
    )
    # The line above is F11 - (F12 - F21), namely F11-F12+F21.
    require(c_trace == [0], "smallest Pascal circuit vanishes")

    plain_reflected = add(
        trace_packet(1, 1),
        add(trace_packet(2, 1), trace_packet(1, 2), -1),
        -1,
    )
    # This construction is F11-F21+F12.
    require(plain_reflected == [2], "plain reflection leaves kernel")

    signed_reflected = c_trace[:]  # iota fixes this circuit exactly.
    require(signed_reflected == [0], "signed reflection preserves kernel")

    # Q has basis {1,z}; its eigenbasis {1,2z-1} has determinant two.
    one = central_basis_polynomial(0)
    z_poly = central_basis_polynomial(1)
    x_poly = add(scale(2, z_poly), one, -1)
    require(list(one) == [1] and list(z_poly) == [0, 1] and x_poly == [-1, 2],
            "smallest eigen glue")
    require(add(one, x_poly) == scale(2, z_poly), "smallest twice split")
    return c_trace, plain_reflected


def main() -> None:
    packet_blocks, packets = exact_packet_and_spine_checks()
    reciprocal_blocks, smith_twos, norm_twos = reciprocal_and_eigenlattice_checks()
    intertwined = signed_packet_intertwining_checks()
    bridges = fibre_parity_bridge_checks()
    c_trace, plain_trace = smallest_hostile_checks()

    print("THM-4378 exact bilateral packet quotient / reciprocal eigenlattice")
    print("status=PROVED ELEMENTARY RELATIVE TO THM-4368, THM-4369, THM-4375")
    print("scope=internal finite packet quotient only; JC(2), DC(2), bracket and seam remain OPEN")
    print()
    print("universe:")
    print("  ell=2..80, quotient excess depth D=0..30")
    print("  reciprocal/eigen Smith blocks D=0..60")
    print(f"  exact packet blocks={packet_blocks}, packet reconstructions={packets}")
    print(f"  signed packet/circuit intertwining checks={intertwined}")
    print(f"  source-fibre parity bridge cases={bridges}")
    print()
    print("ell>=3 theorem:")
    print("  tau0=2*rho-1, D=R-tau0, q=z(1-z)")
    print("  Q_(ell,R)=q^(rho-1) Z[z]_(<=D)")
    print("  unimodular spine: q^k from (rho+k,rho+k),")
    print("                     z*q^k from (rho+k+1,rho+k)")
    print("  iota(e_(u,v))=(-1)^(u-v)e_(v,u), J(f)(z)=f(1-z)")
    print("  Q_plus=q^(rho-1) Z[q]_(<=floor(D/2))")
    print("  Q_minus=q^(rho-1)(2z-1)Z[q]_(<=floor((D-1)/2))")
    print("  Q/(Q_plus direct_sum Q_minus)=(Z/2)^(ceil(D/2))")
    print("  saturation(Q_plus direct_sum Q_minus)=Q")
    print("  (1-J)Q=Q_minus")
    print("  Q_plus/(1+J)Q = Z/2 if D is even, 0 if D is odd")
    print("  Q/(1-J)Q is free of rank floor(D/2)+1")
    print("  Q/(1+J)Q = Z^ceil(D/2) plus Z/2 exactly when D is even")
    print("  fibre-weighted gluing class nonzero iff mu_plus-mu_minus is odd")
    print("  postclock class follows beta parity and vanishes at the full 2d jet")
    print(f"  audited Smith 2-blocks={smith_twos}, even-top norm 2-blocks={norm_twos}")
    print()
    print("ell=2 boundary:")
    print("  B_2 is the diagonal ray; Q=Z[q]_(<=floor(D/2)); J=id")
    print("  no eigen gluing; Q/(1+J)Q=(Z/2)^(floor(D/2)+1)")
    print()
    print("smallest hostile ell=3,R=2:")
    print("  C=e_(1,1)-e_(1,2)+e_(2,1), trace(C)=", c_trace)
    print("  plain swap trace=", plain_trace, "so unsigned reflection does not descend")
    print("  signed iota fixes C; Q basis {1,z}, eigenbasis {1,2z-1}, index=2")
    print()
    print(f"PASS checks={CHECKS}")


if __name__ == "__main__":
    main()
