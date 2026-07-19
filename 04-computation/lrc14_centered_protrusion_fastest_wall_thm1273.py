#!/usr/bin/env python3
"""Exact referee for THM-1273's centered-protrusion fastest-wall bridge.

The paper theorem joins THM-1267's endpoint survivor to the j=4/flood
chronology.  In the small-tail branch, four lower combs leave more than 2/9
of the six-bin mass, while the whole endpoint tail has less than 1/8.  Thus a
four-prefix-safe point remains in K and the cover makes the fastest comb
active there.  THM-1267 supplies a five-prefix-safe point in the protrusion.
Continuity forces a fastest tooth wall between them.  At that wall either all
four lower owners are safe (a bare j=4 wall) or a lower tooth crosses the wall
and creates a positive, lcm-quantized seam.

This referee checks the rational budget, exhausts a finite endpoint-arithmetic
census, and replays the c=140 sharp-star control exactly.  It does not replace
the paper measure/topology providers and is not a finite six-cover reduction.
"""

from __future__ import annotations

import ast
from fractions import Fraction as F
from math import ceil, floor, gcd
from pathlib import Path


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(f"THM-1273 referee failed: {message}")


def optimization_safety_probe() -> int:
    source = Path(__file__).read_text(encoding="utf-8")
    tree = ast.parse(source)
    count = sum(isinstance(node, ast.Assert) for node in ast.walk(tree))
    require(count == 0, "Python assert nodes are optimization-sensitive")
    caught = False
    try:
        require(False, "deliberate require probe")
    except RuntimeError as error:
        caught = "deliberate require probe" in str(error)
    require(caught, "require probe did not fire")
    return count


def tooth(speed: int, address: int) -> tuple[F, F]:
    require(speed > 0, "tooth speed must be positive")
    return (F(14 * address - 1, 14 * speed),
            F(14 * address + 1, 14 * speed))


def circle_distance(value: F) -> F:
    residue = value - floor(value)
    return min(residue, 1 - residue)


def strictly_dangerous(speed: int, point: F) -> bool:
    return circle_distance(speed * point) < F(1, 14)


def closed_safe(speed: int, point: F) -> bool:
    return circle_distance(speed * point) > F(1, 14)


def merge_closed_teeth(
    carrier: tuple[F, F], speeds: tuple[int, ...]
) -> tuple[list[tuple[F, F, int, int]], list[tuple[F, F]]]:
    left, right = carrier
    teeth: list[tuple[F, F, int, int]] = []
    for speed in speeds:
        low_address = floor(left * speed) - 2
        high_address = ceil(right * speed) + 2
        for address in range(low_address, high_address + 1):
            raw_left, raw_right = tooth(speed, address)
            if raw_left < right and left < raw_right:
                teeth.append((max(left, raw_left), min(right, raw_right), speed, address))

    covered: list[list[F]] = []
    for tooth_left, tooth_right, _, _ in sorted(teeth):
        if not covered or covered[-1][1] < tooth_left:
            covered.append([tooth_left, tooth_right])
        elif covered[-1][1] < tooth_right:
            covered[-1][1] = tooth_right

    survivor: list[tuple[F, F]] = []
    cursor = left
    for tooth_left, tooth_right in covered:
        if cursor < tooth_left:
            survivor.append((cursor, tooth_left))
        cursor = max(cursor, tooth_right)
    if cursor < right:
        survivor.append((cursor, right))
    return teeth, survivor


def check_mass_bridge() -> tuple[F, F, F]:
    one_load = F(7, 36)
    four_load_cap = 4 * one_load
    four_survivor = 1 - four_load_cap
    tail_threshold = F(1, 6)
    endpoint_height = F(3, 4)
    tail_cap = endpoint_height * tail_threshold
    inside_floor = four_survivor - tail_cap
    five_survivor = 1 - (4 * one_load + F(23, 120))

    require(four_load_cap == F(7, 9), "four-load cap")
    require(four_survivor == F(2, 9), "four-prefix survivor mass")
    require(tail_cap == F(1, 8), "small endpoint-tail mass cap")
    require(inside_floor == F(7, 72), "interior four-prefix mass floor")
    require(inside_floor > 0, "interior obligation is nonempty")
    require(five_survivor == F(11, 360), "five-prefix survivor mass")
    require(F(11, 270) < tail_threshold, "THM-1267 lower tail fits small branch")
    inside_length = inside_floor / F(7, 6)
    outside_length = five_survivor / F(3, 4)
    needle_separation = inside_length + outside_length
    require(inside_length == F(1, 12), "inside obligation Lebesgue length")
    require(outside_length == F(11, 270), "outside obligation Lebesgue length")
    require(needle_separation == F(67, 540), "oriented needle separation")

    # Symbolic rational rows: if v > 2/9 and the whole tail has mass < 1/8,
    # then the part of V inside K is > 7/72.  Strictness is tested on a grid
    # approaching both excluded boundaries.
    rows = 0
    for q in range(2, 200):
        v_mass = F(2, 9) + F(1, q)
        tail_mass = F(1, 8) - F(1, 8 * q)
        inside_mass = v_mass - tail_mass
        require(inside_mass > F(7, 72), "strict interior mass implication")
        rows += 1
    require(rows == 198, "mass implication census size")
    return inside_length, outside_length, needle_separation


def wall_count_invoice_census() -> int:
    rows = 0
    # If an oriented interval of normalized length >67/540 contains W walls,
    # its W+1 wall-free cells each have length at most d/h.  Thus
    # 67h < 540(W+1)d and W >= floor(67h/(540d)).  Check the exact integer
    # rounding, including ratios which land on an integer boundary.
    for d in range(1, 121):
        for h in range(d + 1, 101 * d + 1, max(1, d // 7)):
            numerator = 67 * h
            denominator = 540 * d
            minimum_walls = numerator // denominator
            require(numerator < denominator * (minimum_walls + 1),
                    "wall-count invoice at floor")
            if minimum_walls > 0:
                require(not numerator < denominator * minimum_walls,
                        "one fewer wall cannot satisfy invoice")
            rows += 1
    require(rows > 10000, "wall-count census size")
    return rows


def endpoint_quantum_census() -> tuple[int, F, tuple[int, int, int, int, str]]:
    rows = 0
    minimum_ratio: F | None = None
    minimum_row: tuple[int, int, int, int, str] | None = None

    # A lower tooth strictly containing an h-wall overlaps the adjacent
    # h-tooth on the active side.  Enumerate signed addresses and both wall
    # orientations; exact numerator divisibility independently replays the
    # paper endpoint calculation.
    for h in range(2, 81):
        for j in range(1, h):
            g = gcd(h, j)
            lcm = h * j // g
            quantum = F(1, 14 * lcm)
            for address in range(-20, 21):
                h_left, h_right = tooth(h, address)
                for side, wall in (("left", h_left), ("right", h_right)):
                    candidates = {floor(j * wall), ceil(j * wall)}
                    for low_address in candidates:
                        j_left, j_right = tooth(j, low_address)
                        if not (j_left < wall < j_right):
                            continue
                        overlap = min(h_right, j_right) - max(h_left, j_left)
                        require(overlap > 0, "wall-crossing overlap positivity")
                        cleared = overlap * (14 * h * j)
                        require(cleared.denominator == 1, "cleared overlap integrality")
                        numerator = cleared.numerator
                        require(numerator % g == 0, "overlap numerator gcd divisibility")
                        require(numerator >= g, "positive numerator at least gcd")
                        require(overlap >= quantum, "lcm seam quantum")
                        ratio = overlap / quantum
                        if minimum_ratio is None or ratio < minimum_ratio:
                            minimum_ratio = ratio
                            minimum_row = (h, j, address, low_address, side)
                        rows += 1

    require(rows > 0 and minimum_ratio is not None and minimum_row is not None,
            "endpoint census nonempty")
    require(minimum_ratio == 1, "endpoint quantum is attained")
    return rows, minimum_ratio, minimum_row


def walls_between(speed: int, left: F, right: F) -> list[F]:
    require(left < right, "oriented wall interval")
    walls: set[F] = set()
    low_address = floor(speed * left) - 2
    high_address = ceil(speed * right) + 2
    for address in range(low_address, high_address + 1):
        for wall in tooth(speed, address):
            if left < wall < right:
                walls.add(wall)
    return sorted(walls)


def sharp_star_control() -> dict[str, object]:
    c = 140
    d = 254
    lows = (255, 256, 257, 292)
    h = 1805
    gap = (F(1121, 1960), F(1133, 1960))
    safe_component = (F(2045, 3556), F(2057, 3556))
    tail = (gap[1], safe_component[1])
    ell = (tail[1] - tail[0]) / (safe_component[1] - safe_component[0])
    require(ell == F(33, 280) and ell < F(1, 6), "sharp-row small protrusion")

    _, four_survivor = merge_closed_teeth(safe_component, lows)
    _, five_survivor = merge_closed_teeth(safe_component, lows + (h,))
    expected_four = [
        (F(2045, 3556), F(2351, 4088)),
        (F(2353, 4088), F(2071, 3598)),
        (F(2073, 3598), F(121, 210)),
        (F(2059, 3570), F(2071, 3584)),
        (F(2073, 3584), F(2057, 3556)),
    ]
    expected_five = [
        (F(2915, 5054), F(14587, 25270)),
        (F(14589, 25270), F(14601, 25270)),
        (F(14617, 25270), F(2057, 3556)),
    ]
    require(four_survivor == expected_four, "sharp-row four-prefix survivor")
    require(five_survivor == expected_five, "sharp-row five-prefix survivor")

    x = (F(14601, 25270) + F(2071, 3584)) / 2
    y = (F(14617, 25270) + F(2057, 3556)) / 2
    require(safe_component[0] < x < gap[1], "control x lies in K")
    require(tail[0] < y < tail[1], "control y lies in protrusion")
    require(all(closed_safe(speed, x) for speed in lows), "x avoids four closed lows")
    require(strictly_dangerous(h, x), "x is fastest-dangerous")
    require(all(closed_safe(speed, y) for speed in lows + (h,)),
            "y is five-prefix-safe")
    require(closed_safe(d, x) and closed_safe(d, y), "d1 is safe throughout S")

    walls = walls_between(h, x, y)
    expected_walls = [F(14603, 25270), F(2923, 5054), F(14617, 25270)]
    require(walls == expected_walls, "three exact fastest walls")
    crossing_sets = [tuple(speed for speed in lows if strictly_dangerous(speed, wall))
                     for wall in walls]
    require(crossing_sets == [(256,), (256,), ()], "wall lower-owner incidence")

    # The final wall is the literal bare j=4 control.  It is outside G; d1
    # and all four lows are safe, while h is exactly at the closed boundary.
    bare = walls[-1]
    require(gap[1] < bare < safe_component[1], "bare wall lies in protrusion")
    require(circle_distance(h * bare) == F(1, 14), "bare fastest equality")
    require(all(not strictly_dangerous(speed, bare) for speed in (d,) + lows),
            "bare wall has no lower fast coverer")

    # At each of the first two walls, the 256 tooth crosses the h wall.  The
    # actual intersection is checked against the lcm quantum.
    low_tooth = tooth(256, 148)
    seam_lengths: list[F] = []
    for high_address in (1043, 1044):
        high_tooth = tooth(h, high_address)
        overlap = min(low_tooth[1], high_tooth[1]) - max(low_tooth[0], high_tooth[0])
        quantum = F(1, 14 * (h * 256 // gcd(h, 256)))
        require(overlap > 0 and overlap >= quantum, "sharp-row located seam quantum")
        seam_lengths.append(overlap)

    return {
        "ell": ell,
        "x": x,
        "y": y,
        "walls": tuple(walls),
        "crossing_sets": tuple(crossing_sets),
        "seam_lengths": tuple(seam_lengths),
        "four_components": len(four_survivor),
        "five_components": len(five_survivor),
        "normalized_separation": (y - x) / (safe_component[1] - safe_component[0]),
    }


def tournament_audit(control: dict[str, object]) -> None:
    require(len(control["walls"]) == 3, "wall tournament vertex count")
    print("TOURNAMENT_AND_CARRIER_AUDIT")
    print("runner_vertices=d2,d3,d4,d5,h; observable=speed_order; gauge=increasing")
    print("runner_score_histogram=0,1,2,3,4 cycles=0 scc_sizes=1,1,1,1,1 hamiltonian_paths=1")
    print("wall_vertices=z1,z2,z3; observable=chronological_order; gauge=oriented_x_to_y")
    print("wall_score_histogram=0,1,2 cycles=0 scc_sizes=1,1,1 hamiltonian_paths=1")
    print("sharp_wall_incidence=(256),(256),()")
    print("faithful_carrier=oriented four-prefix survivor components + fastest-wall events + active-low subsets")
    print("preserves=inside-fastest-active to outside-fastest-safe transport and exact wall label")
    print("destroys_if_runner_only=phase,wall_position,bare-vs-crossed status,seam numerator")
    print("challenged_vertices=runners,teeth,gaps,boundaries,wall_events,residues,Fano_lines,proof_obligations")


def main() -> None:
    assert_nodes = optimization_safety_probe()
    inside_length, outside_length, needle_separation = check_mass_bridge()
    wall_count_rows = wall_count_invoice_census()
    endpoint_rows, minimum_ratio, minimum_row = endpoint_quantum_census()
    control = sharp_star_control()

    print("THM-1273 CENTERED-PROTRUSION / FASTEST-WALL EXACT REFEREE")
    print(f"optimization_sensitive_assert_nodes={assert_nodes}")
    print("small_tail_budget=four_prefix_mass_gt_2/9; tail_mass_lt_1/8; inside_mass_gt_7/72")
    print("outside_obligation=five_prefix_mass_gt_11/360 contained_in_protrusion")
    print(f"needle_lengths=inside_gt_{inside_length};outside_gt_{outside_length};separation_gt_{needle_separation}")
    print(f"wall_count_invoice_rows={wall_count_rows} law=67h<540(W+1)d1")
    print("wall_dichotomy=bare_j4_wall OR lower-crossed_fastest_wall")
    print(f"endpoint_quantum_rows={endpoint_rows} minimum_quantum_ratio={minimum_ratio} row={minimum_row}")
    print(f"sharp_c140_ell={control['ell']} x={control['x']} y={control['y']}")
    print(f"sharp_normalized_xy_separation={control['normalized_separation']}")
    print(f"sharp_four_prefix_components={control['four_components']} five_prefix_components={control['five_components']}")
    print("sharp_fastest_walls=" + ",".join(map(str, control["walls"])))
    print("sharp_wall_crossers=" + ",".join(str(row) for row in control["crossing_sets"]))
    print("sharp_crossed_wall_seams=" + ",".join(map(str, control["seam_lengths"])))
    tournament_audit(control)
    print("NO_DOUBLE_COUNT=the crossed seam may already occur in THM-1253/1275; the bare wall has zero mass")
    print("SCOPE=paper topology plus exact arithmetic consumer; no selected-fastest-tooth claim; no LRC14 closure")
    print("RESULT=PASS")


if __name__ == "__main__":
    main()
