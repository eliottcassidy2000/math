#!/usr/bin/env python3
"""Exact multi-terminal carry-cocycle stopping probe at profile (1,2,5).

The fixed speeds are the hostile local row used in THM-2246 and in the
horizon-two full-simplex probe.  This script refines that single terminal
phase to:

* the original 112-obligation terminal;
* two adjacent membership chambers separated only by a carry wall; and
* the carry wall itself.

It computes the sheet-slope quotient, the relative carry transition, and a
three-terminal overlap loop.  The row is only a local hostile control and is
not asserted to extend to a global scalar cover.
"""

from collections import Counter
from fractions import Fraction
from hashlib import sha256
from math import comb


P = 13
N = P**2
TAU = Fraction(1, 14)
H = 1
C1 = 2 * P
C2 = 2 * P**2
C3 = 2 * P**5
QS = tuple(1 + N * 700000 * i for i in range(1, 6))

Z0 = Fraction(325007, 700000)
Z_WALL = Fraction(157, 338)  # (C3/N) * Z_WALL = 2041


def require(test, message):
    if not test:
        raise RuntimeError(message)


def floor_fraction(value):
    return value.numerator // value.denominator


def norm(value):
    residue = value % 1
    return min(residue, 1 - residue)


def point(z, n):
    return (z + n) / N


def obligation(z, n):
    x = point(z, n)
    return (
        norm(H * x) > Fraction(1, 7)
        and norm(C1 * x) >= TAU
        and all(norm(q * x) >= TAU for q in QS)
    )


def owners(z, n):
    x = point(z, n)
    return tuple(c for c in (C2, C3) if norm(c * x) < TAU)


def nearest_lattice_distance(value, offsets):
    """Distance from value to Z+offsets, with an exact closest witness."""
    base = floor_fraction(value)
    candidates = []
    for offset in offsets:
        for integer in range(base - 2, base + 3):
            candidates.append((abs(value - (integer + offset)), integer, offset))
    return min(candidates)


def membership_boundary_radius(z):
    """Nearest terminal displacement changing any retained membership bit."""
    tests = [
        ("guard", H, (Fraction(1, 7), Fraction(6, 7))),
        ("peeled", C1, (TAU, 1 - TAU)),
        *((f"unit_{i}", q, (TAU, 1 - TAU)) for i, q in enumerate(QS, 1)),
        ("owner_2", C2, (TAU, 1 - TAU)),
        ("owner_3", C3, (TAU, 1 - TAU)),
    ]
    candidates = []
    for label, coefficient, offsets in tests:
        for n in range(N):
            value = coefficient * point(z, n)
            distance, integer, offset = nearest_lattice_distance(value, offsets)
            terminal_distance = N * distance / coefficient
            require(terminal_distance > 0, f"membership wall at {label},{n}")
            candidates.append(
                (terminal_distance, label, n, integer, offset)
            )
    return min(candidates)


def other_carry_boundary_radius(z):
    """Nearest raw floor-carry wall other than the chosen simultaneous C3 wall."""
    coefficients = [
        ("guard", H),
        ("peeled", C1),
        ("owner_2", C2),
        ("owner_3", C3),
        *((f"unit_{i}", q) for i, q in enumerate(QS, 1)),
    ]
    candidates = []
    chosen_wall_count = 0
    for label, coefficient in coefficients:
        for n in range(N):
            value = coefficient * point(z, n)
            base = floor_fraction(value)
            for integer in (base, base + 1):
                distance = abs(value - integer)
                if distance == 0:
                    require(label == "owner_3", "unexpected coincident carry wall")
                    chosen_wall_count += 1
                    continue
                terminal_distance = N * distance / coefficient
                candidates.append((terminal_distance, label, n, integer))
    require(chosen_wall_count == N, "chosen simultaneous carry wall drift")
    return min(candidates)


def slope_quotient(z, c):
    require(c % N == 0, "terminal owner lacks depth two")
    slope = c // N
    return slope, floor_fraction(slope * z)


def carry_potential(z, c, n):
    return floor_fraction(c * point(z, n))


def relative_potential(z):
    slope2, gauge2 = slope_quotient(z, C2)
    slope3, gauge3 = slope_quotient(z, C3)
    return slope2 * gauge3 - slope3 * gauge2


def serialize_sets(sets):
    return "\n".join(
        f"{name}:" + ",".join(map(str, values))
        for name, values in sets
    )


def main():
    require((C1, C2, C3) == (26, 338, 742586), "profile speeds drift")
    require(
        all(q > 0 and q % P for q in QS) and H % 2 and H % P,
        "scalar scope drift",
    )
    require((C1 // P, C2 // N, C3 // N) == (2, 2, 4394), "depth drift")
    require((C3 // N) * Z_WALL == 2041, "chosen carry wall drift")

    nearest = membership_boundary_radius(Z_WALL)
    radius, boundary_label, boundary_n, boundary_integer, boundary_offset = nearest
    require(
        nearest
        == (
            Fraction(155, 233248167061),
            "unit_5",
            72,
            253625740,
            TAU,
        ),
        "nearest membership boundary drift",
    )
    epsilon = radius / 3
    nearest_other_carry = other_carry_boundary_radius(Z_WALL)
    require(
        nearest_other_carry
        == (
            Fraction(157, 199927000338),
            "unit_5",
            60,
            211625740,
        ),
        "nearest other carry wall drift",
    )
    require(epsilon < nearest_other_carry[0], "another carry wall crossed")
    z_minus = Z_WALL - epsilon
    z_plus = Z_WALL + epsilon

    terminals = (
        ("original", Z0),
        ("minus", z_minus),
        ("wall", Z_WALL),
        ("plus", z_plus),
    )
    obligation_sets = {
        name: tuple(n for n in range(N) if obligation(z, n))
        for name, z in terminals
    }
    require(len(obligation_sets["original"]) == 112, "original occupancy drift")
    require(
        obligation_sets["minus"]
        == obligation_sets["wall"]
        == obligation_sets["plus"],
        "adjacent membership chamber drift",
    )
    common = obligation_sets["minus"]
    require(len(common) == 68, "common stalk size drift")
    require(set(common) <= set(obligation_sets["original"]), "overlap drift")

    owner_words = {}
    for name, z in terminals:
        owner_words[name] = Counter(owners(z, n) for n in obligation_sets[name])
    require(owner_words["original"] == Counter({(C2,): 112}), "original owner")
    for name in ("minus", "wall", "plus"):
        require(owner_words[name] == Counter({(C2, C3): 68}), f"{name} owner")

    # The exact sheet-potential identity is what is meant by quotienting the
    # individual linear slopes.  The remainder is a terminal 0-cochain.
    gauges = {}
    for name, z in terminals:
        gauge = []
        for c in (C2, C3):
            slope, residual = slope_quotient(z, c)
            for n in range(N):
                require(
                    carry_potential(z, c, n) == slope * n + residual,
                    f"sheet quotient drift at {name},{c},{n}",
                )
            gauge.append(residual)
        gauges[name] = tuple(gauge)

    require(gauges["original"] == (0, 2040), "original gauge drift")
    require(gauges["minus"] == (0, 2040), "minus gauge drift")
    require(gauges["wall"] == (0, 2041), "wall gauge drift")
    require(gauges["plus"] == (0, 2041), "plus gauge drift")

    # The adjacent open chambers minus/plus are separated only by the C3 carry
    # wall.  No retained membership bit changes in the radius certified above.
    adjacent_delta = tuple(
        right - left for left, right in zip(gauges["minus"], gauges["plus"])
    )
    require(adjacent_delta == (0, 1), "adjacent carry transition drift")
    require(
        floor_fraction((C3 // N) * z_minus) == 2040
        and floor_fraction((C3 // N) * z_plus) == 2041,
        "carry wall crossing drift",
    )
    require(
        epsilon < Fraction(1, C3 // N),
        "more than one C3 carry wall crossed",
    )

    relative = {name: relative_potential(z) for name, z in terminals}
    require(relative == {
        "original": 4080,
        "minus": 4080,
        "wall": 4082,
        "plus": 4082,
    }, "relative potential drift")
    require(
        all(relative_potential(z + 1) == relative_potential(z) for _, z in terminals),
        "relative potential is not periodic",
    )

    # Three-terminal overlap loop: original -> minus -> plus -> original.
    # Every edge has the same 68-obligation common stalk.  The adjacent edge
    # has nonzero relative transition, but the loop sum is zero because the
    # transition is the coboundary of the displayed relative potential.
    loop = ("original", "minus", "plus", "original")
    loop_deltas = tuple(
        relative[right] - relative[left]
        for left, right in zip(loop, loop[1:])
    )
    require(loop_deltas == (0, 2, -2), "loop transition drift")
    require(sum(loop_deltas) == 0, "nonzero loop holonomy")

    # C2 is a common legal section over all three terminals and all common
    # obligations.  Hence the common-stalk coverability complex is a simplex.
    require(
        all(C2 in owners(z, n) for _, z in terminals for n in common),
        "persistent common owner lost",
    )
    switches = 0
    common_faces = 2 ** len(common)
    serialization = serialize_sets(
        tuple((name, obligation_sets[name]) for name, _ in terminals)
    )

    print("scope=FINITE-EXACT local hostile control; not a global cover")
    print("probe=profile_(1,2,5)_multi_terminal_overlap_cocycle")
    print(f"H={H} c1={C1} c2={C2} c3={C3}")
    print(f"unit_speeds={QS}")
    print(f"original_terminal={Z0} carry_wall={Z_WALL}")
    print(
        "nearest_membership_boundary="
        f"{radius} label={boundary_label} n={boundary_n} "
        f"target={boundary_integer}+{boundary_offset}"
    )
    print(
        "nearest_other_carry_boundary="
        f"{nearest_other_carry[0]} label={nearest_other_carry[1]} "
        f"n={nearest_other_carry[2]} target={nearest_other_carry[3]}"
    )
    print(f"epsilon={epsilon} minus={z_minus} plus={z_plus}")
    for name, z in terminals:
        print(
            f"terminal={name} z={z} obligations={len(obligation_sets[name])} "
            f"owners={dict(owner_words[name])} gauge={gauges[name]} "
            f"relative={relative[name]}"
        )
    print(
        f"adjacent_chambers=minus,plus common_obligations={len(common)} "
        f"membership_symmetric_difference=0 carry_delta={adjacent_delta} "
        "relative_jump=2"
    )
    print(
        f"overlap_loop={loop} relative_edge_values={loop_deltas} "
        f"holonomy={sum(loop_deltas)}"
    )
    print(
        f"persistent_section={C2} owner_switches={switches} "
        f"common_complex=full_{len(common)-1}_simplex faces={common_faces} "
        "minimal_nonfaces=()"
    )
    print(
        f"common_pairs={comb(len(common), 2)} "
        f"common_triples={comb(len(common), 3)}"
    )
    print(f"diagram_sha256={sha256(serialization.encode()).hexdigest()}")
    print(
        "general_identity=for owners 169a,169b, "
        "F(z)=a*floor(bz)-b*floor(az) is 1-periodic and dF is exact"
    )
    print(
        "verdict=nonzero adjacent relative carry does not force owner switching; "
        "the three-terminal holonomy vanishes"
    )


if __name__ == "__main__":
    main()
