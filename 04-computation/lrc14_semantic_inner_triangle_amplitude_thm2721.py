#!/usr/bin/env python3
"""Exact amplitude/reanchoring audit for the THM-2712 inner triangle.

The three addresses 0,13,26 carry equal translated copies of the frozen
THM-2680 following atom.  This script distinguishes that positive raw
parallel-arm amplitude from chronological reanchoring into the fixed current
atom.  Every calculation uses exact rational/integer interval arithmetic.
"""

from bisect import bisect_right
from fractions import Fraction
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "04-computation"))

import lrc14_semantic_following_congruence_lock_thm2712 as semantic


m = semantic.m
old = semantic.old


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def strict_contains(value, starts, intervals):
    index = bisect_right(starts, value) - 1
    return index >= 0 and intervals[index][0] < value < intervals[index][1]


def interval_intersection(left, right):
    """Return the strict intersection of two sorted ordinary interval unions."""
    out = []
    i = 0
    j = 0
    while i < len(left) and j < len(right):
        a, b = left[i]
        c, d = right[j]
        if max(a, c) < min(b, d):
            out.append((max(a, c), min(b, d)))
        if b < d:
            i += 1
        else:
            j += 1
    return tuple(out)


def weighted_piece(value, pieces):
    hits = tuple(piece for piece in pieces if piece[0] < value < piece[1])
    require(len(hits) == 1, "strict weighted point lost uniqueness")
    return hits[0]


def main():
    p = 13
    R = p**6
    S = p**5
    T = m.T
    x = Fraction(649039434905733, 1304692766858936)
    z = Fraction(46873542509301, 100360982066072)
    I = (
        Fraction(960117507257, 1930018885886),
        Fraction(324519717452867, 652346383429468),
    )
    radius = Fraction(1, 1304692766858936)
    require(I == (x - radius, x + radius), "inherited cylinder changed")

    module, prefixes, _, _, rails, present, starts = m.core.build_carrier_data()
    rows = m.shard((0, 1))[6][0]
    current = old.d.build_atom(
        module, prefixes, present, starts, rails[2], 3, 2, 0, 0
    )
    following = old.d.build_atom(
        module, prefixes, present, starts, rails[0], 1, 6, 1, 1
    )
    require(
        (current["j"], current["h"], following["j"], following["h"])
        == (5, 2, 2, 6),
        "frozen current/following labels changed",
    )

    # Rebuild the 304 semantic addresses rather than reading the prior output.
    rail_support = old.merge_support(rails[0][3])
    present_support = tuple(present[1, (-6) % p])
    following_support = old.merge_support(following["pieces"])
    denominator = (z * T).denominator
    modulus = T * denominator
    point = (z * T).numerator
    step13 = (91 * T // R) * denominator
    require(91 * T % R == 0, "inner translation left the exact grid")
    scaled_rail = tuple(
        (left * denominator, right * denominator)
        for left, right in rail_support
    )
    scaled_present = tuple(
        (left * denominator, right * denominator)
        for left, right in present_support
    )
    scaled_following = tuple(
        (left * denominator, right * denominator)
        for left, right in following_support
    )
    scaled_rail_starts = tuple(left for left, _ in scaled_rail)
    scaled_present_starts = tuple(left for left, _ in scaled_present)
    scaled_following_starts = tuple(left for left, _ in scaled_following)
    semantic_nodes = []
    for j in range(S):
        n = 13 * j
        root = (6 + n) % p
        carry = (2 + 7 * n) % p
        hit = (
            strict_contains(point, scaled_rail_starts, scaled_rail)
            and strict_contains(point, scaled_present_starts, scaled_present)
            and strict_contains(
                point, scaled_following_starts, scaled_following
            )
            and root
            and m.is_unit(rows[0][0][carry][1][6], root, 26)
        )
        if hit:
            semantic_nodes.append(n)
        point = (point + step13) % modulus
    semantic_nodes = tuple(semantic_nodes)
    require(len(semantic_nodes) == 304, "semantic address count changed")

    triangle = (0, 13, 26)
    require(all(n in semantic_nodes for n in triangle), "inner triangle vanished")
    q_centres = tuple(
        semantic.frac(z + Fraction(7 * n, R)) * T for n in triangle
    )
    following_pieces = tuple(
        weighted_piece(value, following["pieces"]) for value in q_centres
    )
    shift = 91 * T // R
    expected_following_pieces = (
        (139104562225860, 139104590874480, 27582102210),
        (139110177355380, 139110206004000, 27582102210),
        (139115792484900, 139115821133520, 27582102210),
    )
    require(following_pieces == expected_following_pieces,
            "triangle following pieces changed")
    require(
        all(
            following_pieces[i + 1][0] - following_pieces[i][0] == shift
            and following_pieces[i + 1][1] - following_pieces[i][1] == shift
            and following_pieces[i + 1][2] == following_pieces[i][2]
            for i in range(2)
        ),
        "triangle pieces stopped being equal translates",
    )

    current_piece = weighted_piece(x * T, current["pieces"])
    require(
        current_piece
        == (148163522522400, 148163551171020, 27580222516),
        "current weighted piece changed",
    )
    current_weight = current_piece[2]
    following_weight = following_pieces[0][2]
    raw_density = current_weight * following_weight
    raw_integral = raw_density * (I[1] - I[0])
    require(
        raw_density == 760720516410855360360
        and raw_integral == Fraction(2089891528601250990, 1792160394037),
        "equal raw overlap amplitude changed",
    )

    # In Z[omega]/(omega^2+omega+1), both nontrivial transforms of
    # (A,A,A) vanish.  Pairs encode a+b*omega.
    omega_powers = ((1, 0), (0, 1), (-1, -1))
    omega2_powers = ((1, 0), (-1, -1), (0, 1))
    dft1 = tuple(
        raw_integral * sum(power[i] for power in omega_powers)
        for i in (0, 1)
    )
    dft2 = tuple(
        raw_integral * sum(power[i] for power in omega2_powers)
        for i in (0, 1)
    )
    require(dft1 == dft2 == (0, 0), "equal-arm C3 transform changed")

    # The three comparison gains are a coboundary and telescope in Z.
    gains = (91, 91, -182)
    require(
        tuple(7 * (b - a) for a, b in zip(
            (0, 13, 26), (13, 26, 0)
        )) == gains
        and sum(gains) == 0,
        "inner gain telescope changed",
    )

    # Chronological reanchoring fails before any coefficient cancellation:
    # the complete following base support is disjoint from current rail 2.
    current_rail_support = old.merge_support(rails[2][3])
    support_intersection = interval_intersection(
        following_support, current_rail_support
    )
    require(not support_intersection, "following support met current rail 2")
    current_base_support = old.merge_support(current["pieces"])
    require(
        not interval_intersection(following_support, current_base_support),
        "following support met the fixed current base support",
    )

    # The failure is not repairable by merely changing rail or future clock
    # while retaining current depth h=2.  The following atom occupies the
    # epsilon=1 half of root 6, whereas the two possible h=2 halves belong to
    # root 10; their exact 182-grid intervals are disjoint.
    following_half = (14 * 6 - 13, 14 * 6)
    h2_halves = ((14 * 10 - 13, 14 * 10), (14 * 10, 14 * 10 + 13))
    require(
        following_half == (71, 84)
        and h2_halves == ((127, 140), (140, 153))
        and all(max(following_half[0], left)
                >= min(following_half[1], right)
                for left, right in h2_halves),
        "private half/root disjointness changed",
    )

    endpoint_current_midpoint_hits = []
    endpoint_current_whole_hits = []
    for n in semantic_nodes:
        q = semantic.frac(z + Fraction(7 * n, R))
        base = q * T
        delayed = semantic.frac(R * q) * T
        if any(left < base < right
               for left, right in current_base_support):
            endpoint_current_midpoint_hits.append(n)
        if semantic.atom_whole_I(
            current,
            base,
            13 * T * radius,
            delayed,
            13 * R * T * radius,
            prefixes,
        ):
            endpoint_current_whole_hits.append(n)
    require(
        not endpoint_current_midpoint_hits and not endpoint_current_whole_hits,
        "a semantic endpoint unexpectedly reanchored into the fixed current",
    )

    # If the frozen following atom is promoted to the next current, then a
    # chronological successor in the same sharp-graph grammar must have
    # predecessor j equal to its h=6.  The graph equation has one solution,
    # whose private root-one half is disjoint from every semantic endpoint.
    second_generation_labels = []
    for h in range(1, 12):
        for epsilon in (0, 1):
            for kappa in (0, 1):
                root = (-h - 1) % p
                j = old.d.INV2 * (root - epsilon - kappa) % p
                if j == following["h"]:
                    second_generation_labels.append(
                        (h, epsilon, kappa, root)
                    )
    require(
        tuple(second_generation_labels) == ((11, 1, 1, 1),),
        "second-generation label solution changed",
    )
    successor_half = (1, 14)
    semantic_deep_phases = tuple(
        semantic.frac(
            module.C3 * semantic.frac(z + Fraction(7 * n, R))
        ) * 182
        for n in semantic_nodes
    )
    require(
        all(Fraction(71) < phase < Fraction(84)
            and not Fraction(successor_half[0]) < phase
            < Fraction(successor_half[1])
            for phase in semantic_deep_phases),
        "a semantic endpoint entered the unique second-generation root half",
    )
    second_generation_placements = len(rails) * 7
    require(second_generation_placements == 1134,
            "rail/future placement count changed")

    # Sharp changed-source control.  For n=13j,
    # q_n(x)=D(x+7j/R).  Rebuilding the source on that translated cylinder
    # retains the fixed current for exactly four of the 304 semantic arms.
    shifted_source_hits = []
    for n in semantic_nodes:
        require(n % 13 == 0, "semantic address left the locked congruence")
        j = n // 13
        shifted_x = semantic.frac(x + Fraction(7 * j, R))
        q = semantic.frac(z + Fraction(7 * n, R))
        require(semantic.frac(13 * shifted_x) == q,
                "changed-source factorization failed")
        if semantic.atom_whole_I(
            current,
            shifted_x * T,
            T * radius,
            semantic.frac(R * shifted_x) * T,
            R * T * radius,
            prefixes,
        ):
            shifted_source_hits.append(n)
    require(
        tuple(shifted_source_hits) == (0, 169, 338, 507),
        "changed-source positive control changed",
    )

    print("LRC14 SEMANTIC INNER TRIANGLE AMPLITUDE AUDIT")
    print(f"p={p} R={R} T={T} semantic_nodes={len(semantic_nodes)}")
    print(f"triangle={triangle} gains={gains} gain_sum={sum(gains)}")
    print(f"following_grid_shift={shift}")
    print(f"following_pieces={following_pieces}")
    print(f"following_piece_length={following_pieces[0][1]-following_pieces[0][0]}")
    print(f"current_piece={current_piece}")
    print(f"raw_density={raw_density} raw_integral={raw_integral}")
    print(f"abstract_C3_nontrivial_DFT={dft1},{dft2}")
    print("following_support_intersect_current_rail2=empty")
    print("following_support_intersect_current_atom=empty")
    print(
        "private_half_no_go="
        f"following={following_half} h2_halves={h2_halves}"
    )
    print(
        "endpoint_current_reanchor="
        f"midpoint_hits={len(endpoint_current_midpoint_hits)} "
        f"whole_I_hits={len(endpoint_current_whole_hits)}"
    )
    print(
        "second_generation="
        f"labels={tuple(second_generation_labels)} "
        f"successor_half={successor_half} "
        f"rail_future_placements={second_generation_placements} hits=0"
    )
    print(f"changed_source_whole_I_hits={tuple(shifted_source_hits)}")
    print("scope=equal_raw_parallel_arms_not_chronological_endpoint_current")
    print("ALL CHECKS PASSED")


if __name__ == "__main__":
    main()
