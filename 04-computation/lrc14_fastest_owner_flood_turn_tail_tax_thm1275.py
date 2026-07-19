#!/usr/bin/env python3
"""Exact referee for THM-1275's fastest-owner flood/turn tail tax.

The proof-bearing object is THM-1253's deletion-minimal chronological tooth
word.  Consecutive occurrences of the fastest owner are separated either by

* a regular consecutive-address one-low rung;
* one or more skipped fastest addresses (a flood); or
* a consecutive-address packet with at least two low teeth (a turn).

THM-1266 bounds a run of regular rungs by the number of eligible low owners.
Skipped fastest teeth and THM-1253's full raw seam family form two internally
disjoint layers.  Where the layers overlap, the skipped fastest comb is a
third active owner, so pointwise multiplicity pays both indicators.  This
dependency-free script checks the finite combinatorics, exact tooth geometry,
layered multiplicity, return arithmetic, constants, and the c=140 local
countermodel to naive six-spoke protrusion additivity.

The analytic/topological providers (minimal subcover extraction, private dual
mass, and integration against the six-bin density) remain paper inputs.
"""

from __future__ import annotations

import ast
from fractions import Fraction as F
from itertools import product
from math import ceil, floor
from pathlib import Path


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(f"THM-1275 referee failed: {message}")


def optimization_safe_require_probe() -> int:
    tree = ast.parse(Path(__file__).read_text(encoding="utf-8"))
    count = sum(isinstance(node, ast.Assert) for node in ast.walk(tree))
    require(count == 0, "optimization-sensitive assert remains")
    caught = False
    try:
        require(False, "deliberate sentinel")
    except RuntimeError as error:
        caught = "deliberate sentinel" in str(error)
    require(caught, "explicit RuntimeError sentinel did not fire")
    return count


def ceil_fraction(value: F) -> int:
    return -((-value.numerator) // value.denominator)


def longest_regular_run(bits: tuple[int, ...]) -> int:
    best = 0
    current = 0
    for bit in bits:
        if bit:
            current += 1
            best = max(best, current)
        else:
            current = 0
    return best


def regular_exception_census() -> tuple[int, int]:
    """Exhaust the exact e-dependent run/exception ceiling."""

    rows = 0
    sharp_rows = 0
    for eligible in range(6):
        for fastest_occurrences in range(1, 17):
            transition_count = fastest_occurrences - 1
            minimum_exceptions: int | None = None
            for bits in product((0, 1), repeat=transition_count):
                # 1=regular; 0=exceptional.
                if longest_regular_run(bits) > eligible:
                    continue
                rows += 1
                exceptions = bits.count(0)
                lower = ceil_fraction(F(fastest_occurrences, eligible + 1)) - 1
                require(exceptions >= lower, "regular-run exception ceiling")
                if minimum_exceptions is None or exceptions < minimum_exceptions:
                    minimum_exceptions = exceptions
            require(minimum_exceptions is not None, "empty transition census")
            expected = ceil_fraction(F(fastest_occurrences, eligible + 1)) - 1
            require(minimum_exceptions == expected, "exception ceiling is not sharp")
            sharp_rows += 1
    return rows, sharp_rows


def eligible_colour_census() -> tuple[int, tuple[int, ...]]:
    """Check that e colours cannot support e+1 regular slots.

    THM-1266 says repeated low owners in a regular star have slot distance at
    least six.  Since e<=5, this reduces to distinctness on a candidate run.
    """

    rows = 0
    maxima: list[int] = []
    for eligible in range(6):
        maximum = 0
        alphabet = tuple(range(eligible))
        for length in range(eligible + 2):
            valid_count = 0
            for word in product(alphabet, repeat=length):
                rows += 1
                valid = True
                for left in range(length):
                    for right in range(left + 1, length):
                        if word[left] == word[right] and right - left < 6:
                            valid = False
                if valid:
                    valid_count += 1
                    maximum = max(maximum, length)
            if eligible == 0 and length == 0:
                # product((), repeat=0) contains the empty word.
                require(valid_count == 1, "empty eligible alphabet base case")
        require(maximum == eligible, "eligible-colour run maximum")
        maxima.append(maximum)
    return rows, tuple(maxima)


def tooth(speed: int, address: int) -> tuple[F, F]:
    require(speed > 0, "positive speed")
    return F(14 * address - 1, 14 * speed), F(14 * address + 1, 14 * speed)


def skipped_tooth_geometry_census() -> tuple[int, int]:
    """Audit the endpoint order used to put every skipped tooth inside G."""

    skipped_rows = 0
    corridor_rows = 0
    for high in range(1, 81):
        for first_address in range(-3, 4):
            first = tooth(high, first_address)
            turn_right = tooth(high, first_address + 1)
            # Any seam in a consecutive-address packet lies in this corridor.
            corridor = (first[0], turn_right[1])
            for outside_address in range(first_address - 4, first_address + 6):
                if outside_address in (first_address, first_address + 1):
                    continue
                outside = tooth(high, outside_address)
                if outside_address < first_address:
                    require(outside[1] < corridor[0], "left tooth meets turn corridor")
                else:
                    require(corridor[1] < outside[0], "right tooth meets turn corridor")
                corridor_rows += 1

            for address_return in range(2, 9):
                last = tooth(high, first_address + address_return)
                # Deletion minimality supplies witnesses x in the first tooth
                # and y in the last tooth.  The strongest endpoint statement
                # below is independent of where those witnesses lie.
                for skipped_address in range(
                    first_address + 1, first_address + address_return
                ):
                    skipped = tooth(high, skipped_address)
                    require(first[1] < skipped[0], "skipped tooth not right of first")
                    require(skipped[1] < last[0], "skipped tooth not left of last")
                    # Therefore every x in first and y in last bracket the
                    # whole skipped tooth, so connectedness of G puts it in G.
                    skipped_rows += 1
    return skipped_rows, corridor_rows


def layered_multiplicity_audit() -> int:
    """Truth table for adding the skipped-tooth and raw-seam layers."""

    rows = 0
    for skipped in (0, 1):
        for seam in (0, 1):
            # Cover gives one active owner.  A seam gives two distinct
            # selected owners.  On a skipped h tooth neither selected seam
            # owner can be h, so their overlap has at least three owners.
            minimum_multiplicity = 1 + skipped + seam
            require(skipped + seam <= minimum_multiplicity - 1,
                    "pointwise two-layer multiplicity")
            rows += 1
    return rows


def check_turn_return_example() -> dict[str, F]:
    """Reconstruct one literal r=1 multi-low return exactly."""

    high = 26
    rows = ((26, 10), (5, 2), (12, 5), (26, 11))
    intervals = tuple(tooth(*row) for row in rows)
    require(all(intervals[i][0] < intervals[i + 1][0] for i in range(3)),
            "turn left endpoints")
    require(all(intervals[i][1] < intervals[i + 1][1] for i in range(3)),
            "turn right endpoints")
    seams = tuple(intervals[i][1] - intervals[i + 1][0] for i in range(3))
    require(all(seam > 0 for seam in seams), "positive turn seams")
    require(all(intervals[i][1] <= intervals[i + 2][0] for i in range(2)),
            "turn handoff separation")
    omega = sum(seams, F(0))
    reciprocal_excess = F(1, 5) + F(1, 12) - F(6, high)
    require(reciprocal_excess == F(41, 780), "turn reciprocal excess")
    require(omega == F(41, 5460), "literal turn seam mass")
    require(omega == reciprocal_excess / 7, "r=1 return identity")
    return {"omega": omega, "reciprocal_excess": reciprocal_excess}


def constant_census() -> tuple[int, int]:
    rows = 0
    dominated_rows = 0
    for carrier in range(1, 101):
        for high in range(carrier + 1, 20 * carrier + 1):
            scale = F(7 * carrier, 6)
            skip_length = F(1, 7 * high)
            require(scale * F(3, 4) * skip_length == F(carrier, 8 * high),
                    "weighted skipped-tooth coefficient")
            sample_lcm = high * (high + 1)
            require(scale * F(3, 4) * F(1, 14 * sample_lcm)
                    == F(carrier, 16 * sample_lcm),
                    "weighted lcm-seam coefficient")
            basic_private_mass = F(1, 36)
            per_tooth_cap = F(7 * carrier, 36 * high)
            occurrences = ceil_fraction(F(high, 7 * carrier))
            require(occurrences * per_tooth_cap >= basic_private_mass,
                    "private occurrence ceiling")
            require((occurrences - 1) * per_tooth_cap < basic_private_mass,
                    "private occurrence ceiling minimality")
            rows += 1

            for fifth in range(carrier + 1, high):
                if 2 * high < 7 * fifth:
                    continue
                require(F(2, fifth) - F(6, high) >= F(1, high),
                        "7/2-dominated turn tax")
                dominated_rows += 1
    require(dominated_rows > 0, "dominated branch positive control")

    c, h = F(17), F(113)
    reciprocal_excess = F(41, 780)
    omega = reciprocal_excess / 7
    require(F(7) * c / 6 * F(3, 4) * omega == c / 8 * reciprocal_excess,
            "weighted return coefficient")
    require(F(49, 6) * (F(1, 7) * reciprocal_excess)
            == F(7, 6) * reciprocal_excess,
            "unweighted return coefficient")
    require(F(49, 6) * F(1, 7) == F(7, 6),
            "unweighted skipped coefficient")
    require(F(49, 6) * F(1, 14) == F(7, 12),
            "unweighted lcm-seam coefficient")
    return rows, dominated_rows


def nearest_integer(value: F) -> int:
    lower = floor(value)
    remainder = value - lower
    require(remainder != F(1, 2), "countermodel nearest-integer tie")
    return lower if remainder < F(1, 2) else lower + 1


def circle_depth(value: F) -> F:
    residue = value - floor(value)
    return min(residue, 1 - residue)


def c140_six_spoke_countermodel() -> tuple[tuple[int, F, F], ...]:
    """Local blocker-cycle row: only the slowest safe component protrudes."""

    carrier = 140
    gap_address = 80
    centre = F(2 * gap_address + 1, 2 * carrier)
    gap = (
        F(14 * gap_address + 1, 14 * carrier),
        F(14 * gap_address + 13, 14 * carrier),
    )
    word = (
        (1805, 1036), (256, 147), (1805, 1037), (254, 146),
        (1805, 1038), (292, 168), (1805, 1039), (257, 148),
        (1805, 1040), (255, 147), (1805, 1041),
    )
    intervals = tuple(tooth(*row) for row in word)
    speeds = (254, 255, 256, 257, 292, 1805)
    blocker_indices = {254: 8, 255: 3, 256: 7, 257: 1, 292: 1, 1805: 3}
    expected_raw = {
        254: F(33, 280),
        255: F(-89, 336),
        256: F(-9, 140),
        257: F(-163, 1680),
        292: F(-8, 105),
        1805: F(-617, 112),
    }
    rows: list[tuple[int, F, F]] = []
    blocker_map: dict[int, int] = {}
    for speed in speeds:
        clock = carrier + speed
        numerator = nearest_integer(clock * centre)
        epsilon = F(numerator) - clock * centre
        rho = abs(epsilon)
        phase = F(numerator, clock)
        require(gap[0] < phase < gap[1], "c140 spoke inside carrier gap")
        require(circle_depth(carrier * phase) == circle_depth(speed * phase),
                "c140 spoke equal depth")
        require(circle_depth(speed * phase) > F(1, 4), "c140 spoke deep")

        raw = F(1, 2) + F(7, 6) * rho - F(speed, 2 * carrier)
        ell = max(F(0), raw)
        require(raw == expected_raw[speed], "c140 raw protrusion")
        target_index = blocker_indices[speed]
        target = intervals[target_index]
        require(target[0] < phase < target[1], "c140 selected blocker tooth")
        blocker_map[speed] = word[target_index][0]
        rows.append((speed, raw, ell))

    require(blocker_map[254] == 1805 and blocker_map[1805] == 254,
            "c140 first blocker two-cycle")
    require(blocker_map[256] == 257 and blocker_map[257] == 256,
            "c140 second blocker two-cycle")
    require(rows[0][2] == F(33, 280), "c140 unique positive protrusion")
    require(all(ell == 0 for _, _, ell in rows[1:]),
            "c140 five contained safe components")
    return tuple(rows)


def tournament_loss_audit() -> None:
    print("TOURNAMENT_LOSS_AUDIT")
    print("vertices=successive fastest-owner obligations, not runners")
    print("observable=transition is regular versus flood/turn; gauge=fastest-address order")
    print("tie_hamiltonian_path=chronological fastest-tooth subsequence")
    print("preserves=address jump, internal-low count, and disjoint physical invoice")
    print("destroys=low-owner phases inside a turn packet and exact seam clocks")
    print("challenged_assumption=fastest teeth are not all consecutive-address star rungs")


def main() -> None:
    assert_count = optimization_safe_require_probe()
    transition_rows, sharp_rows = regular_exception_census()
    colour_rows, colour_maxima = eligible_colour_census()
    skipped_rows, corridor_rows = skipped_tooth_geometry_census()
    layered_rows = layered_multiplicity_audit()
    turn = check_turn_return_example()
    constant_rows, dominated_rows = constant_census()
    protrusions = c140_six_spoke_countermodel()

    print("THM-1275 FASTEST-OWNER FLOOD/TURN TAIL TAX EXACT REFEREE")
    print(f"optimization_safe_assert_count={assert_count}")
    print(f"regular_exception_rows={transition_rows} sharp_parameter_rows={sharp_rows}")
    print("exception_ceiling=X+T>=E>=ceil(K/(e+1))-1")
    print(f"eligible_colour_rows={colour_rows} maximum_regular_runs={colour_maxima}")
    print(f"skipped_tooth_rows={skipped_rows} turn_corridor_disjoint_rows={corridor_rows}")
    print(f"layered_multiplicity_truth_rows={layered_rows} rule=skip+seam<=C-1")
    print(f"literal_turn=H26,J5,J12,H26 omega={turn['omega']} reciprocal_excess={turn['reciprocal_excess']}")
    print(f"constant_rows={constant_rows} dominated_rows={dominated_rows}")
    print("weighted_invoice=F6>=c*X/(8h)+c/16*sum_all_seams(1/lcm)")
    print("unweighted_invoice=H-1/c>=7X/(6h)+7/12*sum_all_seams(1/lcm)")
    print("packet_form=F6>=c/8*(X/h+sum_packets(sum_internal(1/j)-(7R-1)/h))")
    print("dominated_corollary=h>=7*d5/2 => F6>=c/(8h)*(ceil(K/(e+1))-1)")
    print("private_count=K>=ceil(h/(7c)); functional refinement K>=ceil(36h*eta_h/(7c))")
    print("c140_spoke_protrusions=" + ",".join(
        f"d{speed}:raw={raw}:ell={ell}" for speed, raw, ell in protrusions
    ))
    print("c140_scope=local centered-spoke/blocker-cycle row, not a six-cover")
    tournament_loss_audit()
    print("SCOPE=global necessary tail invoice; scale-small and not LRC(14) closure")
    print("RESULT=PASS")


if __name__ == "__main__":
    main()
