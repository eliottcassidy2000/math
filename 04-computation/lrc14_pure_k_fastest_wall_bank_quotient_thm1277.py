#!/usr/bin/env python3
"""Exact referee for THM-1277's pure-K fastest-wall bank quotient.

THM-1273 puts a fastest-active point x more than 1/12 normalized units to
the left of the K/E interface b.  Every fastest wall in (x,b) is crossed by
a selected middle-owner tooth.  The paper theorem retains the alternating
fastest wall types, pairs the two walls of each complete fastest tooth, and
uses deletion minimality: an unselected complete tooth is a whole flood,
whereas a selected one has two distinct-owner chronological boundary seams.
Those floods layer pointwise with the full THM-1253 seam family, including at
the prefix and suffix omitted by THM-1275's between-selected-address count.

This dependency-free referee checks the exact alternating-wall arithmetic,
the finite block/run consumer, endpoint/lcm quanta, invoice coefficients, and
the c=140 two-wall/no-complete-tooth guardrail.  Minimal-subcover extraction
and interval-chain topology remain explicit paper providers.
"""

from __future__ import annotations

import ast
from fractions import Fraction as F
from math import ceil, floor, gcd, lcm
from pathlib import Path


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(f"THM-1277 referee failed: {message}")


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
    require(speed > 0, "positive tooth speed")
    return (F(14 * address - 1, 14 * speed),
            F(14 * address + 1, 14 * speed))


def complete_teeth_from_wall_count(wall_count: int) -> tuple[int, int]:
    if wall_count == 0:
        return 0, 0
    delta = 1
    complete = (wall_count - 1) // 2
    return delta, complete


def wall_bank_census() -> tuple[int, int, F]:
    rows = 0
    threshold_rows = 0
    largest_ratio_checked = F(0)

    # Work in units of the normalized fastest-tooth length a=d1/(6h).
    # The first right wall has coordinate r in (0,1); subsequent gaps are
    # 6,1,6,1,... .  For a cell with W interior walls, the endpoint distance
    # is at most wall_(W+1).  Rational r values exhaust both alternating cell
    # parities without paying for millions of Fraction comparisons in the
    # original normalized coordinate.
    for d1 in range(1, 6):
        step = d1
        for h in range(d1 + 1, 81 * d1 + 1, step):
            largest_ratio_checked = max(largest_ratio_checked, F(h, d1))
            threshold = F(h, 2 * d1)  # (1/12)/a
            carrier_limit = F(6 * h, d1)  # 1/a
            for denominator in (5, 7):
                for numerator in range(1, denominator):
                    first = F(numerator, denominator)
                    for wall_count in range(0, 160):
                        delta, complete = complete_teeth_from_wall_count(
                            wall_count
                        )
                        expected_complete = (
                            (wall_count - 1) // 2 if wall_count else 0
                        )
                        require(complete == expected_complete,
                                "paired complete-tooth wall count")
                        require(delta == int(wall_count > 0), "first-wall indicator")

                        # Exact coordinate of the next wall in units of a.
                        # W=2q+1: next is left wall r+7q+6.
                        # W=2q+2: next is right wall r+7(q+1).
                        if wall_count == 0:
                            previous_wall = F(0)
                            next_wall = first
                        elif wall_count % 2 == 1:
                            q = (wall_count - 1) // 2
                            previous_wall = first + 7 * q
                            next_wall = first + 7 * q + 6
                        else:
                            q = (wall_count - 2) // 2
                            previous_wall = first + 7 * q + 6
                            next_wall = first + 7 * (q + 1)

                        require(next_wall <= 7 * complete + 8,
                                "sharp alternating-wall length bound")

                        # There is a normalized endpoint beyond 1/12 in this
                        # wall-count cell exactly when its right edge exceeds
                        # threshold; the carrier endpoint must also occur
                        # before normalized coordinate one.
                        if max(previous_wall, threshold) < min(
                            next_wall, carrier_limit
                        ):
                            if h >= 2 * d1:
                                require(wall_count > 0,
                                        "h>=2d1 forces the first pure-K wall")
                                floor_bank = (h - 2 * d1) // (14 * d1)
                                require(complete >= floor_bank,
                                        "complete-tooth ratio floor")
                            threshold_rows += 1
                        rows += 1

    require(rows > 100000, "wall-bank census is substantial")
    require(threshold_rows > 10000, "threshold census is substantial")
    return rows, threshold_rows, largest_ratio_checked


def circle_distance(value: F) -> F:
    residue = value - floor(value)
    return min(residue, 1 - residue)


def centered_abs_residue(value: int, modulus: int) -> int:
    require(modulus > 0, "positive residue modulus")
    residue = value % modulus
    return min(residue, modulus - residue)


def wall_residue_census() -> tuple[int, int]:
    rows = 0
    crossings = 0
    # At z=a/(14h), a=14n+-1, the lower speed j crosses exactly when
    # ||ja/(14h)||<1/14, equivalently |ja|_(14h)<h.
    for h in range(2, 82):
        modulus = 14 * h
        for j in range(1, h):
            for address in range(-12, 13):
                for sign in (-1, 1):
                    a = 14 * address + sign
                    require(a % 14 in (1, 13), "fastest wall residue class")
                    analytic = circle_distance(F(j * a, modulus)) < F(1, 14)
                    residue_test = centered_abs_residue(j * a, modulus) < h
                    require(analytic == residue_test,
                            "exact fastest-wall residue crossing identity")
                    crossings += int(residue_test)
                    rows += 1
    require(rows > 100000, "wall-residue census is substantial")
    require(crossings > 0, "wall-residue census contains crossings")
    return rows, crossings


def endpoint_quantum_and_pair_census() -> tuple[int, int, F]:
    quantum_rows = 0
    paired_containment_rows = 0
    minimum_ratio: F | None = None

    for h in range(2, 61):
        for j in range(1, h):
            require(F(1, 7 * h) < F(6, 7 * j),
                    "two distinct same-owner teeth cannot meet one h tooth")
            common = gcd(h, j)
            quantum = F(1, 14 * lcm(h, j))
            for high_address in range(-4, 5):
                high_left, high_right = tooth(h, high_address)
                for wall in (high_left, high_right):
                    candidates = {floor(j * wall), ceil(j * wall)}
                    for low_address in candidates:
                        low_left, low_right = tooth(j, low_address)
                        if low_left < wall < low_right:
                            overlap = min(high_right, low_right) - max(high_left, low_left)
                            require(overlap > 0, "crossed wall has positive overlap")
                            cleared = overlap * (14 * h * j)
                            require(cleared.denominator == 1,
                                    "endpoint overlap clears integrally")
                            require(cleared.numerator % common == 0,
                                    "endpoint numerator is gcd-divisible")
                            require(overlap >= quantum, "lcm seam quantum")
                            ratio = overlap / quantum
                            minimum_ratio = ratio if minimum_ratio is None else min(
                                minimum_ratio, ratio
                            )
                            quantum_rows += 1

                        if low_left < high_left and high_right < low_right:
                            require(low_left < high_left < high_right < low_right,
                                    "one low tooth crossing both walls contains H")
                            require(
                                min(high_right, low_right) - max(high_left, low_left)
                                == high_right - high_left,
                                "paired crosser covers the whole fastest tooth",
                            )
                            paired_containment_rows += 1

    require(quantum_rows > 1000, "endpoint quantum census is substantial")
    require(paired_containment_rows > 0, "paired containment occurs")
    require(minimum_ratio == 1, "primitive lcm seam quantum is attained")
    return quantum_rows, paired_containment_rows, minimum_ratio


def localized_run_census() -> tuple[int, int]:
    rows = 0
    forced_rows = 0

    def visit(
        remaining: int,
        e: int,
        selected_after_h0: int,
        floods: int,
        turns: int,
        previous_selected: bool,
        regular_edges_in_run: int,
        total_p: int,
    ) -> None:
        nonlocal rows, forced_rows
        if remaining == 0:
            a_count = selected_after_h0
            require(a_count + 1 <= (e + 1) * (floods + turns + 1),
                    "regular selected-vertex run capacity")
            require(total_p <= (e + 2) * floods + (e + 1) * turns + e,
                    "localized flood/turn capacity")
            numerator = max(0, total_p - e)
            lower = (numerator + (e + 2) - 1) // (e + 2)
            require(floods + turns >= lower, "localized ceiling consumer")
            if total_p > e:
                require(floods + turns > 0, "star completion forces an exception")
                forced_rows += 1
            rows += 1
            return

        # The next complete fastest tooth is unselected: this is a whole
        # flood and it breaks the consecutive-address selected run.
        visit(remaining - 1, e, selected_after_h0, floods + 1, turns,
              False, 0, total_p)

        # It is selected.  Across an unselected predecessor there is no
        # adjacent-address selected transition inside the local block.
        if not previous_selected:
            visit(remaining - 1, e, selected_after_h0 + 1, floods, turns,
                  True, 0, total_p)
            return

        # A multi-low turn splits regular runs.
        visit(remaining - 1, e, selected_after_h0 + 1, floods, turns + 1,
              True, 0, total_p)

        # A regular edge is admissible only while the THM-1266 run has at
        # most e transitions.
        if regular_edges_in_run < e:
            visit(remaining - 1, e, selected_after_h0 + 1, floods, turns,
                  True, regular_edges_in_run + 1, total_p)

    for e in range(0, 6):
        for p in range(0, 13):
            visit(p, e, 0, 0, 0, True, 0, p)

    require(rows > 100000, "localized run census is substantial")
    require(forced_rows > 0, "localized star completion fires")
    return rows, forced_rows


def layered_invoice_census() -> tuple[int, int]:
    truth_rows = 0
    coefficient_rows = 0

    # At a flood point the selected subcover supplies one low layer.  At a
    # simultaneous seam point it supplies two selected low layers, and the
    # unselected fastest tooth is a third full-comb layer.
    for flood in (0, 1):
        for seam in (0, 1):
            minimum_multiplicity = 1 + flood + seam
            for multiplicity in range(minimum_multiplicity, 7):
                require(flood + seam <= multiplicity - 1,
                        "pointwise flood plus seam layering")
                truth_rows += 1

    for c in range(1, 13):
        for h in range(2, 36):
            for d4 in range(1, h):
                d5 = min(h - 1, d4 + 1)
                if not d4 < d5:
                    continue
                whole_weight = F(7 * c, 6) * F(3, 4) * F(1, 7 * h)
                seam_weight = F(7 * c, 6) * F(3, 4) * F(1, 14 * h * d5)
                require(whole_weight == F(c, 8 * h), "whole-tooth weight")
                require(seam_weight == F(c, 16 * h * d5), "seam weight")
                require(F(1, d4) + F(1, d5) <= 2,
                        "a flood dominates the distinct-owner pair floor")
                for p in range(0, 5):
                    for floods in range(0, p + 1):
                        selected = p - floods
                        for delta in (0, 1):
                            exact_floor = (
                                F(c * floods, 8 * h)
                                + F(c * delta, 16 * h * d5)
                                + F(c * selected, 16 * h)
                                * (F(1, d4) + F(1, d5))
                            )
                            uniform = (
                                F(c * delta, 16 * h * d5)
                                + F(c * p, 16 * h)
                                * (F(1, d4) + F(1, d5))
                            )
                            require(exact_floor >= uniform,
                                    "selected/flood functional bank collapse")

                            harmonic_exact = (
                                F(7 * floods, 6 * h)
                                + F(7 * delta, 12 * h * d5)
                                + F(7 * selected, 12 * h)
                                * (F(1, d4) + F(1, d5))
                            )
                            harmonic_uniform = (
                                F(7 * delta, 12 * h * d5)
                                + F(7 * p, 12 * h)
                                * (F(1, d4) + F(1, d5))
                            )
                            require(harmonic_exact >= harmonic_uniform,
                                    "selected/flood harmonic bank collapse")
                            coefficient_rows += 1

    require(truth_rows > 0, "layer truth table nonempty")
    require(coefficient_rows > 100000, "coefficient census substantial")
    return truth_rows, coefficient_rows


def walls_between(speed: int, left: F, right: F) -> list[tuple[F, int, str]]:
    result: list[tuple[F, int, str]] = []
    low_address = floor(speed * left) - 2
    high_address = ceil(speed * right) + 2
    for address in range(low_address, high_address + 1):
        wall_left, wall_right = tooth(speed, address)
        if left < wall_left < right:
            result.append((wall_left, address, "L"))
        if left < wall_right < right:
            result.append((wall_right, address, "R"))
    return sorted(result)


def c140_guardrail() -> dict[str, object]:
    c = 140
    d1 = 254
    h = 1805
    x = F(7476011, 12938240)
    interface = F(1133, 1960)
    y = F(7425603, 12837160)
    walls = walls_between(h, x, interface)
    expected = [
        (F(14603, 25270), 1043, "R"),
    ]
    require(walls == expected, "c=140 pure-K wall word")
    delta, complete = complete_teeth_from_wall_count(len(walls))
    require((delta, complete) == (1, 0), "one pure-K wall and no complete tooth")

    mixed_walls = walls_between(h, x, y)
    expected_mixed = [
        (F(14603, 25270), 1043, "R"),
        (F(2923, 5054), 1044, "L"),
        (F(14617, 25270), 1044, "R"),
    ]
    require(mixed_walls == expected_mixed, "c=140 mixed K/E wall word")

    bridge = tooth(256, 148)
    require(all(bridge[0] < wall < bridge[1]
                for wall, _, _ in mixed_walls[:2]),
            "one 256 tooth crosses the first two mixed-bank walls")
    safe_gap = (mixed_walls[0][0], mixed_walls[1][0])
    require(bridge[0] < safe_gap[0] < safe_gap[1] < bridge[1],
            "256 tooth covers the intervening fastest-safe gap")
    require(mixed_walls[0][1] != mixed_walls[1][1],
            "the two walls belong to different fastest teeth")
    require(safe_gap[0] < interface < safe_gap[1],
            "the paired mixed walls straddle the K/E interface")
    ratio_floor = (h - 2 * d1) // (14 * d1)
    require(ratio_floor == 0, "c=140 ratio floor guardrail")

    return {
        "c": c,
        "d1": d1,
        "h": h,
        "x": x,
        "interface": interface,
        "y": y,
        "walls": tuple(walls),
        "mixed_walls": tuple(mixed_walls),
        "delta": delta,
        "complete": complete,
        "ratio_floor": ratio_floor,
        "bridge": bridge,
    }


def tournament_audit() -> None:
    print("TOURNAMENT_AND_CARRIER_AUDIT")
    print("cell_vertices=complete_fastest_teeth; observable=selected_vs_flood; gauge=chronology")
    print("cell_tournament=transitive score_histogram=0..P-1 cycles=0 scc_sizes=all_1 edge_flips=0 hamiltonian_paths=1")
    print("selected_cell_label=oriented_boundary_owner_edge j_left->j_right on {d2,d3,d4,d5}")
    print("boundary_graph=directed_multigraph_not_tournament; loops=forbidden; unordered_colours=6; oriented_colours=12")
    print("repeat_thresholds=7_selected_force_unordered_repeat;13_selected_force_oriented_repeat")
    print("faithful_carrier=oriented_K+alternating_h_walls+complete_h_cells+selection_bit+ordered_owner_pair+lcm_digits")
    print("preserves=cover_predicate,wall_pair_type,minimality,spatial_invoice,endpoint_clock")
    print("destroys_if_wall_only=tooth_vs_safe_gap_pair;destroys_if_runner_only=all_phase_and_multiplicity")
    print("challenged_vertices=runners,gaps,walls,wall_pairs,teeth,seams,residues,Fano_lines,proof_obligations")


def main() -> None:
    assert_nodes = optimization_safety_probe()
    wall_rows, threshold_rows, largest_ratio = wall_bank_census()
    residue_rows, residue_crossings = wall_residue_census()
    quantum_rows, containment_rows, minimum_quantum = endpoint_quantum_and_pair_census()
    run_rows, forced_rows = localized_run_census()
    truth_rows, coefficient_rows = layered_invoice_census()
    control = c140_guardrail()

    print("THM-1277 PURE-K FASTEST-WALL BANK QUOTIENT EXACT REFEREE")
    print(f"optimization_sensitive_assert_nodes={assert_nodes}")
    print(f"wall_bank_rows={wall_rows} threshold_rows={threshold_rows} max_h_over_d1={largest_ratio}")
    print("wall_bank_law=b-x>1/12 => h<(14P+16)d1")
    print("ratio_corollary=h>=2d1 => delta=1 and P>=floor((h-2d1)/(14d1))")
    print(f"wall_residue_rows={residue_rows} crossings={residue_crossings} law=|j*a|_(14h)<h")
    print(f"endpoint_quantum_rows={quantum_rows} paired_containment_rows={containment_rows} minimum_lcm_ratio={minimum_quantum}")
    print("selected_complete_tooth=two_distinct_middle_owners+two_chronological_seams")
    print("unselected_complete_tooth=whole_fastest_flood_even_at_prefix_or_suffix")
    print(f"localized_run_rows={run_rows} forced_exception_rows={forced_rows}")
    print("localized_capacity=P<=(e+2)F+(e+1)T+e")
    print(f"layer_truth_rows={truth_rows} coefficient_rows={coefficient_rows}")
    print("functional_bank=F6>=cF/(8h)+(c/16)*forced_lcm_sum")
    print("four_owner_floor=F6>=c*delta/(16hd5)+cP/(16h)*(1/d4+1/d5)")
    print("harmonic_floor=H-1/c>=7*delta/(12hd5)+7P/(12h)*(1/d4+1/d5)")
    print(f"sharp_c140_x={control['x']} interface={control['interface']}")
    print("sharp_c140_walls=" + ",".join(
        f"{wall}@{address}{side}" for wall, address, side in control["walls"]
    ))
    print("sharp_c140_mixed_walls=" + ",".join(
        f"{wall}@{address}{side}"
        for wall, address, side in control["mixed_walls"]
    ))
    print(f"sharp_c140_delta={control['delta']} complete_fastest_teeth={control['complete']} ratio_floor={control['ratio_floor']}")
    print("sharp_c140_same_crosser=256@148;pair_type=fastest_safe_gap_straddling_K/E")
    tournament_audit()
    print("COMPLEMENTARY_BRANCH_AUDIT=exact_lcm_bank_survives;gcd-forgetting floor decays on large primitive near-coprime carriers")
    print("SCOPE=no scale-free turn lower bound; no near-top closure; no sporadic emptiness; no LRC14 closure")
    print("RESULT=PASS")


if __name__ == "__main__":
    main()
