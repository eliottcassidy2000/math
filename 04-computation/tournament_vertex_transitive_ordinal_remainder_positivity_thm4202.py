#!/usr/bin/env python3
"""Primary exact audit for THM-4202 VT ordinal-remainder positivity.

This imports the current THM-4184 exact subset-DP/ordinal-transfer engine.  The
all-order theorem is proved symbolically in the companion markdown;
the finite rows here are hostile controls, not the source of the proof.
"""

from __future__ import annotations

import hashlib
import itertools
import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "04-computation"))

import lrc_a000568_iso_analogy_s509 as lrc  # noqa: E402
import tournament_ordinal_cocycle_parity_thm4184 as base  # noqa: E402


def need(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def circulant(order: int, connection: set[int]) -> tuple[int, ...]:
    return tuple(
        sum(
            1 << right
            for right in range(order)
            if right != left and (right - left) % order in connection
        )
        for left in range(order)
    )


def circulant_bank(order: int) -> tuple[tuple[int, ...], ...]:
    if order == 1:
        return (base.transitive(1),)
    pairs = tuple((step, order - step) for step in range(1, (order + 1) // 2))
    answer = tuple(
        circulant(order, {pairs[i][bit] for i, bit in enumerate(choice)})
        for choice in itertools.product((0, 1), repeat=len(pairs))
    )
    need(len(set(answer)) == 2 ** ((order - 1) // 2), "circulant bank size")
    return answer


def vt_formula(left: base.TournamentData, right: base.TournamentData) -> Fraction:
    m, n = len(left.out), len(right.out)
    h, s = left.hamilton, right.hamilton
    w, v = left.mass, right.mass
    cross_mass = w * v + w * s + h * v + 2 * h * s
    return (
        h * s * w * v
        + cross_mass * s * w * Fraction(m + 2, m)
        + cross_mass * h * v * Fraction(n - 6, n)
        + Fraction(cross_mass * cross_mass, 2)
        * (1 + Fraction(1, m * n) + Fraction(3, m) - Fraction(5, n))
    )


def incoming_masses(data: base.TournamentData) -> tuple[int, ...]:
    order = len(data.out)
    return tuple(
        sum(
            data.capacities[i][j]
            for j in range(order)
            if data.out[j] & (1 << i)
        )
        for i in range(order)
    )


def one_sided_vt_formula(
    left: base.TournamentData, right: base.TournamentData
) -> tuple[Fraction, tuple[int, ...]]:
    """Exact formula when only the left factor is vertex-transitive."""
    m = len(left.out)
    h, s = left.hamilton, right.hamilton
    w, v = left.mass, right.mass
    cross = tuple(
        2
        * (
            left.starts[0][0] * right.ends[j][0]
            + left.starts[0][1] * right.ends[j][1]
        )
        for j in range(len(right.out))
    )
    total = sum(cross)
    squares = sum(value * value for value in cross)
    degrees = tuple(sum(row) for row in right.capacities)
    incoming = incoming_masses(right)
    linear = sum(
        cross[j] * (v - degrees[j] - 4 * incoming[j])
        for j in range(len(right.out))
    )
    value = (
        h * s * w * v
        + (m + 2) * s * w * total
        + m * h * linear
        + Fraction(m, 2)
        * ((m + 3) * total * total - (5 * m - 1) * squares)
    )
    return value, cross


def odd_symmetry_defect(
    left: base.TournamentData,
    right: base.TournamentData,
    cross: tuple[int, ...],
) -> tuple[Fraction, Fraction, Fraction]:
    """Return variance, covariance, and the exact uniformization penalty."""
    m, n = len(left.out), len(right.out)
    h = left.hamilton
    degrees = tuple(sum(row) for row in right.capacities)
    incoming = incoming_masses(right)
    mean_cross = Fraction(sum(cross), n)
    mean_degree = Fraction(sum(degrees), n)
    mean_incoming = Fraction(sum(incoming), n)
    variance = sum((Fraction(value) - mean_cross) ** 2 for value in cross)
    covariance = sum(
        (Fraction(cross[j]) - mean_cross)
        * (
            Fraction(degrees[j])
            - mean_degree
            + 4 * (Fraction(incoming[j]) - mean_incoming)
        )
        for j in range(n)
    )
    penalty = m * h * covariance + Fraction(m * (5 * m - 1), 2) * variance
    return variance, covariance, penalty


def theta(
    left: base.TournamentData,
    middle: base.TournamentData,
    right: base.TournamentData,
) -> int:
    return (
        base.remainder(base.ordinal_data(left, middle), right)
        - left.hamilton**2 * base.remainder(middle, right)
    )


def reverse_triangle(
    out: tuple[int, ...], triangle: tuple[int, int, int]
) -> tuple[int, ...]:
    answer = list(out)
    for left, right in itertools.combinations(triangle, 2):
        if answer[left] & (1 << right):
            answer[left] ^= 1 << right
            answer[right] |= 1 << left
        else:
            answer[right] ^= 1 << left
            answer[left] |= 1 << right
    return tuple(answer)


def is_directed_triangle(
    out: tuple[int, ...], triangle: tuple[int, int, int]
) -> bool:
    return all(
        sum(
            bool(out[vertex] & (1 << other))
            for other in triangle
            if other != vertex
        )
        == 1
        for vertex in triangle
    )


def main() -> None:
    banks = {
        order: tuple(base.tournament_data(out) for out in circulant_bank(order))
        for order in (1, 3, 5, 7, 9)
    }
    lefts = tuple(data for bank in banks.values() for data in bank)
    rights = tuple(data for order, bank in banks.items() if order >= 3 for data in bank)
    digest = hashlib.sha256()

    checks = 0
    minimum = None
    for left in lefts:
        for right in rights:
            direct = base.remainder(left, right)
            formula = vt_formula(left, right)
            need(formula.denominator == 1, "nonintegral VT formula")
            need(formula.numerator == direct, "VT formula mismatch")
            need(direct > 0, "nonpositive VT remainder")
            row = (
                direct,
                len(left.out),
                len(right.out),
                base.label(left.out),
                base.label(right.out),
            )
            if minimum is None or row < minimum:
                minimum = row
            digest.update(("vt|" + "|".join(map(str, row)) + "\n").encode("ascii"))
            checks += 1

    one = banks[1][0]
    cycle = banks[3][0]
    prefixes = [one]
    for _ in range(2):
        prefixes.append(base.ordinal_data(one, prefixes[-1]))
    prefix_checks = 0
    prefix_minimum = None
    for prefix in prefixes:
        prefixed_cycle = base.ordinal_data(prefix, cycle)
        for right in rights:
            direct = base.remainder(prefixed_cycle, right)
            base_value = base.remainder(cycle, right)
            increment = theta(prefix, cycle, right)
            need(direct == base_value + increment, "prefix cocycle mismatch")
            need(direct > base_value > 0 and increment > 0, "prefix strictness")
            row = (direct, len(prefix.out), len(right.out), base.label(right.out))
            if prefix_minimum is None or row < prefix_minimum:
                prefix_minimum = row
            digest.update(("prefix|" + "|".join(map(str, row)) + "\n").encode("ascii"))
            prefix_checks += 1

    # One-sided symmetry/defect identity on every labelled odd right factor
    # at orders three and five, with C3 as the VT left factor.
    one_sided_checks = 0
    for order in (3, 5):
        pair_count = order * (order - 1) // 2
        for bits in range(1 << pair_count):
            out = base.parse(format(bits, f"0{pair_count}b"), order)
            right = base.tournament_data(out)
            direct = base.remainder(cycle, right)
            formula, cross = one_sided_vt_formula(cycle, right)
            need(formula.denominator == 1 and formula.numerator == direct, "one-sided formula")
            variance, covariance, penalty = odd_symmetry_defect(cycle, right, cross)
            uniform = vt_formula(cycle, right)
            need(uniform - penalty == direct, "odd symmetry-defect identity")
            digest.update(
                (
                    "defect|"
                    + "|".join(
                        map(
                            str,
                            (
                                order,
                                base.label(out),
                                direct,
                                uniform,
                                variance,
                                covariance,
                                penalty,
                            ),
                        )
                    )
                    + "\n"
                ).encode("ascii")
            )
            one_sided_checks += 1

    # Score-regular hostile: a directed-triangle reversal preserves every
    # outdegree but destroys the uniform rooted/capacity sidecars.
    order = 9
    seed = circulant(order, set(range(1, 5)))
    hostile = None
    for triangle in itertools.combinations(range(order), 3):
        if not is_directed_triangle(seed, triangle):
            continue
        out = reverse_triangle(seed, triangle)
        need(len({row.bit_count() for row in out}) == 1, "score switch failed")
        data = base.tournament_data(out)
        states = tuple(data.ends)
        degrees = tuple(sum(row) for row in data.capacities)
        incoming = tuple(
            sum(
                data.capacities[i][j]
                for j in range(order)
                if out[j] & (1 << i)
            )
            for i in range(order)
        )
        if len(set(states)) > 1 or len(set(degrees)) > 1 or len(set(incoming)) > 1:
            direct = base.remainder(cycle, data)
            false_uniform = vt_formula(cycle, data)
            need(false_uniform.denominator == 1, "hostile false formula integral")
            need(direct != false_uniform, "hostile did not separate")
            hostile = (
                triangle,
                base.label(out),
                data.hamilton,
                states,
                degrees,
                incoming,
                direct,
                false_uniform.numerator,
                false_uniform.numerator - direct,
            )
            break
    need(hostile is not None, "regular hostile not found")
    hostile_cross = tuple(2 * (row[0] + 2 * row[1]) for row in hostile[3])
    hostile_variance, hostile_covariance, hostile_penalty = odd_symmetry_defect(
        cycle, data, hostile_cross
    )
    need(hostile_penalty == hostile[8], "hostile defect penalty")
    digest.update(("hostile|" + repr(hostile) + "\n").encode("ascii"))

    # Existing LRC projection hostile, replayed directly at its two rational
    # times.  The pointed tournament class agrees but the exact safe bit flips.
    speeds = (0, 1, 3, 4)
    lrc_rows = []
    for time in (Fraction(41, 96), Fraction(7, 48)):
        positions = lrc.positions(speeds, time)
        adjacency = lrc.phase_tournament(positions)
        distances = tuple(lrc.dist0(point) for point in positions[1:])
        safe = all(distance >= Fraction(1, 4) for distance in distances)
        lrc_rows.append(
            (
                time,
                positions,
                distances,
                safe,
                lrc.canonical(adjacency, fixed_zero=False),
                lrc.canonical(adjacency, fixed_zero=True),
                lrc.count_hp(adjacency),
                tuple(sorted(map(sum, adjacency))),
                sum(adjacency[0]),
            )
        )
    need(lrc_rows[0][3] and not lrc_rows[1][3], "LRC safe flip")
    need(lrc_rows[0][4] == lrc_rows[1][4], "LRC unmarked class mismatch")
    need(lrc_rows[0][5] == lrc_rows[1][5], "LRC pointed class mismatch")
    digest.update(("lrc|" + repr(lrc_rows) + "\n").encode("ascii"))

    print("VT_ORDINAL_REMAINDER_PRIMARY_ACCEPT")
    print("circulant_banks", " ".join(f"n{n}={len(banks[n])}" for n in banks))
    print("vt_pair_checks", checks)
    print("vt_minimum", minimum)
    print("prefix_cycle_checks", prefix_checks)
    print("prefix_cycle_minimum", prefix_minimum)
    print("one_sided_odd_defect_checks", one_sided_checks)
    print("regular_hostile_triangle", hostile[0])
    print("regular_hostile_label", hostile[1])
    print("regular_hostile_hamilton", hostile[2])
    print("regular_hostile_end_states", hostile[3])
    print("regular_hostile_capacity_degrees", hostile[4])
    print("regular_hostile_incoming_masses", hostile[5])
    print(
        "regular_hostile_remainder",
        hostile[6],
        "false_uniform",
        hostile[7],
        "difference",
        hostile[8],
    )
    print(
        "regular_hostile_symmetry_defect",
        "cross",
        hostile_cross,
        "variance",
        hostile_variance,
        "covariance",
        hostile_covariance,
        "penalty",
        hostile_penalty,
    )
    for row in lrc_rows:
        print(
            "lrc_projection_row",
            row[0],
            "positions",
            row[1],
            "distances",
            row[2],
            "safe",
            row[3],
            "H",
            row[6],
            "score",
            row[7],
            "stationary",
            row[8],
        )
    print("lrc_same_unmarked", lrc_rows[0][4] == lrc_rows[1][4])
    print("lrc_same_pointed", lrc_rows[0][5] == lrc_rows[1][5])
    print("semantic_sha256", digest.hexdigest())


if __name__ == "__main__":
    main()
