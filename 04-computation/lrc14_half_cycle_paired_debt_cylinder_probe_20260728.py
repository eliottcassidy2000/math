#!/usr/bin/env python3
"""Exact paired-shift probe on the THM-2698 half-cycle cylinders.

This dependency-free scout tests equation (18) of
``lrc14-guard-source-debt-cone-is-target-null-20260728.md`` on the two
canonical THM-2698 event packets.  It inserts the deepest translated danger
probe and the first lawful blocker--graft safe pair on one open ancestry
interval.  The pivot ``k_a`` is not selected by the retained half-cycle
labels, so every canonically admissible ordinary-unit candidate is tested.

All support and wall calculations use exact ``Fraction`` arithmetic.  The
nonvanishing character test is exact reduction modulo
Phi_13(z)=1+z+...+z^12; no floating-point Fourier evaluation is used.
"""

from fractions import Fraction


P = 13
R = P**6

GUARD = 1
UNITS = (14, 27, 40, 53, 66)
SOURCE = P
A = P**3
C = 2 * P**5
K_B = 14
K_A_CANDIDATES = tuple(unit for unit in UNITS if unit != K_B)

EVENTS = (
    (Fraction(55_232_507, 24 * R), 7),
    (Fraction(58_313_459, 24 * R), 3),
)
PACKET_RADIUS = Fraction(1, 301_082_946_198_216)
THREE_STATE_RADIUS = Fraction(1, 50_883_017_907_498_504)
EXPECTED_DEEP_RADIUS = Fraction(71, 810_903_912)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def floor_fraction(value):
    return value.numerator // value.denominator


def fractional_part(value):
    return value - floor_fraction(value)


def centered(value):
    """The half-open representative in [-1/2,1/2)."""
    residue = fractional_part(value)
    return residue if residue < Fraction(1, 2) else residue - 1


def danger_phase(value):
    """The strict d_1 tooth: distance to Z is less than 1/14."""
    return abs(centered(value)) < Fraction(1, 14)


def wall_distance(value):
    """Distance in phase space to the nearest +/-1/14 wall modulo Z."""
    base = floor_fraction(value)
    walls = (
        Fraction(integer) + sign * Fraction(1, 14)
        for integer in range(base - 2, base + 3)
        for sign in (-1, 1)
    )
    return min(abs(value - wall) for wall in walls)


def x_radius_to_wall(coefficient, phase):
    require(coefficient != 0, "zero factor coefficient")
    return wall_distance(phase) / abs(coefficient)


def cyclotomic_reduction(coefficients):
    """Reduce a degree-at-most-12 vector modulo Phi_13."""
    require(len(coefficients) == P, "wrong cyclotomic coefficient length")
    top = coefficients[-1]
    return tuple(value - top for value in coefficients[:-1])


def exact_character_reductions(support):
    """Return reductions of sum_(s in support) z^(h*s), h=1,...,12."""
    reductions = []
    for frequency in range(1, P):
        coefficients = [0] * P
        for shift in support:
            coefficients[(frequency * shift) % P] += 1
        reduced = cyclotomic_reduction(coefficients)
        require(any(reduced),
                "a claimed primitive 13-character vanished")
        reductions.append(reduced)
    return tuple(reductions)


def paired_support_and_radius(x, root, k_a):
    deep_phase = C * x - Fraction(root, P)
    require(danger_phase(deep_phase), "selected deep root is not dangerous")
    deep_radius = x_radius_to_wall(C, deep_phase)

    support = []
    all_cell_radii = []
    for shift in range(P):
        a_phase = A * x - Fraction(shift, P)
        graft_phase = k_a * x + Fraction(shift, P)
        cell_radius = min(
            deep_radius,
            x_radius_to_wall(A, a_phase),
            x_radius_to_wall(k_a, graft_phase),
        )
        all_cell_radii.append(cell_radius)
        if not danger_phase(a_phase) and not danger_phase(graft_phase):
            support.append(shift)

    return tuple(support), min(all_cell_radii), deep_phase, deep_radius


def main():
    require((GUARD, UNITS, SOURCE, A, C) == (
        1, (14, 27, 40, 53, 66), 13, 2_197, 742_586,
    ), "canonical typed row changed")
    require(K_A_CANDIDATES == (27, 40, 53, 66),
            "unresolved k_a candidate set changed")
    require(PACKET_RADIUS == 169 * THREE_STATE_RADIUS,
            "THM-2698 packet/three-state radius relation changed")

    expected_supports = (
        {
            27: (0, 3, 4, 5, 8, 9, 10, 11, 12),
            40: (0, 1, 2, 3, 4, 5, 8, 9, 10, 11),
            53: (0, 1, 2, 3, 4, 5, 8, 11, 12),
            66: (0, 1, 2, 3, 4, 5, 8, 9, 10, 11, 12),
        },
        {
            27: (1, 2, 3, 4, 7, 8, 9, 10, 11),
            40: (1, 2, 3, 4, 5, 6, 7, 8, 9, 10),
            53: (1, 2, 3, 6, 7, 8, 9, 10, 11),
            66: (1, 2, 3, 4, 5, 6, 7, 8, 9),
        },
    )

    rows = []
    global_radius = None
    for event_index, (x, root) in enumerate(EVENTS):
        require(0 < x - PACKET_RADIUS < x + PACKET_RADIUS < 1,
                "open packet interval wrapped around the circle")
        require(all(not danger_phase(unit * x)
                    for unit in K_A_CANDIDATES),
                "an unresolved k_a candidate acquired the old-head danger anchor")

        event_rows = []
        for k_a in K_A_CANDIDATES:
            support, radius, deep_phase, deep_radius = (
                paired_support_and_radius(x, root, k_a)
            )
            require(centered(deep_phase) == Fraction(-1, 156),
                    "selected deep-root phase changed")
            require(deep_radius == EXPECTED_DEEP_RADIUS,
                    "selected deep-root wall radius changed")
            require(radius == EXPECTED_DEEP_RADIUS,
                    "closest three-factor wall radius changed")
            require(PACKET_RADIUS < radius,
                    "canonical packet cylinder crosses a paired-factor wall")
            require(support == expected_supports[event_index][k_a],
                    "paired-shift support changed")
            require(0 < len(support) < P,
                    "paired-shift support became constant")

            reductions = exact_character_reductions(support)
            require(len(reductions) == P - 1,
                    "primitive character census changed")

            # On the open interval |x-x_i|<PACKET_RADIUS all three factors
            # are constant for every shift.  Thus the exact integral is
            # 2*PACKET_RADIUS on support and zero off support.
            table = tuple(
                2 * PACKET_RADIUS if shift in support else Fraction(0)
                for shift in range(P)
            )
            require(sum(value > 0 for value in table) == len(support),
                    "positive table support mismatch")

            event_rows.append((
                k_a,
                support,
                len(support),
                0 in support,
                len(reductions),
            ))
            global_radius = radius if global_radius is None else min(
                global_radius, radius
            )
        rows.append(tuple(event_rows))

    require(global_radius == EXPECTED_DEEP_RADIUS,
            "global paired-factor radius changed")
    require(global_radius / PACKET_RADIUS == 26_361_803,
            "packet-to-paired-wall safety ratio changed")
    require(all(row[3] for row in rows[0])
            and all(not row[3] for row in rows[1]),
            "zero-shift occupancy pattern changed")

    print("LRC14 THM2698 PAIRED DEBT-CYLINDER EXACT PROBE")
    print("status=FINITE-EXACT positive moving-endpoint partial result; "
          "endpoint-current typing remains open")
    print(f"typed_row=(H={GUARD},units={UNITS},j={SOURCE},a={A},c={C}); "
          f"k_b={K_B}; unresolved_k_a={K_A_CANDIDATES}")
    print(f"events=(x,r)={EVENTS}")
    print(f"open_packet_radius={PACKET_RADIUS}; "
          f"open_packet_mass={2 * PACKET_RADIUS}; "
          f"three_state_radius={THREE_STATE_RADIUS}")
    print(f"closest_inserted_factor_wall_radius={global_radius}; "
          f"wall_radius_over_packet_radius={global_radius/PACKET_RADIUS}")
    print("table_law=K_(i,k_a)(s)=2*open_packet_radius on support, 0 off support")
    print(f"event_rows=(k_a,support,size,zero_shift_present,"
          f"nonzero_primitive_character_count)={tuple(rows)}")
    print("character_certificate=for every row and every h in F_13^*, "
          "sum_(s in support) zeta_13^(h*s) is nonzero by exact reduction "
          "mod Phi_13")
    print("zero_plane_boundary=all four k_a candidates are safe at both old "
          "heads; event0 contains s=0 and event1 omits it only because the "
          "a-factor changes, not because k_a is dangerous")
    print("first_typing_obstruction=THM2698 retains no atomic left relation "
          "index u, so eta.u is undefined and the surviving character types "
          "only the moving residue -eta.v, not eta.(u-v)")
    print("secondary_typing_obstructions=k_a is not selected among four "
          "ordinary roles; no canonical ancestry map identifies the THM2701 "
          "BABA debt SCC with the THM2698 B1 half-cycle")
    print("SCOPE=no completed THM2334 target current, semantic endpoint, "
          "scalar-row exclusion, or LRC14 conclusion")


if __name__ == "__main__":
    main()
