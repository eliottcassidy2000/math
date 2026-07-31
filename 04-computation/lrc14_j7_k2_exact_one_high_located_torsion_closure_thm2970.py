#!/usr/bin/env python3
"""Exact forced-high closure for the projected k=2 band 1680..1742.

It reads the immutable scalar-atlas output, reconstructs the danger-comb
carrier and every unit ray independently, and then applies the located
torsion-cell lemma to every remaining literal upper class.  The 39 rows whose
first label already lies above the projected wall are reported but are outside
the one-high/three-low reduction.
"""

from __future__ import annotations

import argparse
import hashlib
from collections import Counter, defaultdict
from fractions import Fraction as F
from math import gcd, lcm
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
ATLAS_SOURCE = ROOT / "04-computation/lrc14_j7_k2_scalar_band_1680_1742_thm2941.py"
ATLAS_OUTPUT = ROOT / "05-knowledge/results/lrc14_j7_k2_scalar_band_1680_1742_thm2941.out"
DEFAULT_OUTPUT = ROOT / "05-knowledge/results/lrc14_j7_k2_exact_one_high_located_torsion_closure_thm2970.out"

EXPECTED_HASHES = {
    ATLAS_SOURCE: "1224d5594571f21c91c55fe3ab165c4fc34ba7968719862d12660d24efac919d",
    ATLAS_OUTPUT: "c20607cb478ed491d7000f2b8a49213f57d1606a5152700ac3b50c836e2dc66c",
}

RULER = 14 * lcm(*range(1, 15))
ONE_FOURTEENTH = RULER // 14
THIRTEEN_FOURTEENTHS = 13 * ONE_FOURTEENTH
EXPECTED_BEFORE = ((1683, 1), (1694, 10), (1702, 3), (1708, 14),
                   (1722, 11), (1724, 2), (1732, 2), (1736, 15))
EXPECTED_AFTER = ((1694, 7), (1702, 3), (1708, 9), (1722, 7),
                  (1724, 2), (1732, 2), (1736, 9))
EXPECTED_TWO_HIGH_MIN = F(121_670_443, 941_781_292_179)
EXPECTED_LOW2_HIGH_MARGIN = F(3_131, 11_953_620)
EXPECTED_CLOSEST_EMPTY_GAP = F(1_323_977_953, 569_777_681_768_295)
EXPECTED_CLASS_COUNT = 13
EXPECTED_CLASS_SEMANTIC_SHA256 = "1642cd1b5f7bcca7c79abf37850e7eee3c7b4aede31a17dada82dc1cac6dc427"


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def sha256(path: Path) -> str:
    payload = path.read_bytes().replace(b"\r\n", b"\n").replace(b"\r", b"\n")
    return hashlib.sha256(payload).hexdigest()


def ftext(value: F) -> str:
    return f"{value.numerator}/{value.denominator}"


def parse_fraction(text: str) -> F:
    return F(text)


def parse_atlas():
    rows = []
    for line in ATLAS_OUTPUT.read_text(encoding="utf-8").splitlines():
        if not line.startswith("SURVIVOR;"):
            continue
        fields = dict(
            field.split("=", 1)
            for field in line.split(";")[1:]
            if "=" in field
        )
        suffix = []
        for item in fields["suffix"].split(","):
            kind, label, value = item.split(":")
            suffix.append(
                (
                    kind,
                    None if label == "None" else int(label),
                    parse_fraction(value),
                )
            )
        rows.append(
            {
                "body": tuple(map(int, fields["E"].split(","))),
                "h": parse_fraction(fields["h"]),
                "components": int(fields["r"]),
                "L": int(fields["L"]),
                "floor": int(fields["largest_floor"]),
                "first": int(fields["z1"]),
                "first_delta": parse_fraction(fields["delta1"]),
                "suffix": tuple(suffix),
                "lower": parse_fraction(fields["lower"]),
                "upper": parse_fraction(fields["upper"]),
                "gap": parse_fraction(fields["gap"]),
                "high_tail": any(kind == "HIGH-TAIL" for kind, _label, _value in suffix),
            }
        )
    return tuple(rows)


def merge_intervals(intervals):
    rows = sorted((left, right) for left, right in intervals if left < right)
    if not rows:
        return ()
    merged = [[rows[0][0], rows[0][1]]]
    for left, right in rows[1:]:
        if left <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], right)
        else:
            merged.append([left, right])
    return tuple((left, right) for left, right in merged)


def danger_intervals(speed: int):
    require(RULER % (14 * speed) == 0, ("inexact body ruler", speed))
    step = RULER // speed
    radius = RULER // (14 * speed)
    return (
        ((0, radius),)
        + tuple(
            (index * step - radius, index * step + radius)
            for index in range(1, speed)
        )
        + ((RULER - radius, RULER),)
    )


def carrier_for(body):
    danger = merge_intervals(
        interval
        for speed in body
        for interval in danger_intervals(speed)
    )
    carrier = []
    cursor = 0
    for left, right in danger:
        if cursor < left:
            carrier.append((cursor, left))
        cursor = max(cursor, right)
    if cursor < RULER:
        carrier.append((cursor, RULER))
    return tuple(carrier)


def primitive_vector(numerators):
    cycles = numerators // RULER
    remainders = numerators % RULER
    return (
        cycles * (RULER // 7)
        + np.minimum(remainders, ONE_FOURTEENTH)
        + np.maximum(0, remainders - THIRTEEN_FOURTEENTHS)
    )


def amplitude_data(body, modulus):
    carrier = carrier_for(body)
    mass_numerator = sum(right - left for left, right in carrier)
    labels = np.arange(1, modulus, dtype=np.int64)
    coverage = np.zeros(modulus - 1, dtype=np.int64)
    for left, right in carrier:
        coverage += primitive_vector(labels * right) - primitive_vector(labels * left)
    amplitudes = 7 * coverage - mass_numerator * labels
    require(
        int(np.max(np.abs(amplitudes))) * modulus < 2**63,
        ("int64 comparison risk", body, modulus),
    )
    return carrier, mass_numerator, amplitudes


def delta(amplitudes, label):
    return F(int(amplitudes[label - 1]), 7 * RULER * label)


def high_ray_maxima(amplitudes, modulus, high_floor):
    """One attained positive maximum, or the safe zero supremum, for each d."""
    best = {}
    for residue, raw in enumerate(amplitudes, 1):
        amplitude = int(raw)
        denominator = modulus // gcd(modulus, residue)
        if denominator == 1:
            continue
        incumbent = best.get(denominator, (0, None, None))
        if amplitude > 0:
            label = residue if residue >= high_floor else residue + modulus
            old_amplitude, old_label, _old_residue = incumbent
            if old_label is None or amplitude * old_label > old_amplitude * label:
                incumbent = (amplitude, label, residue)
        best[denominator] = incumbent
    return {
        denominator: (
            F(amplitude, 7 * RULER * label) if label is not None else F(0),
            label,
            residue,
            amplitude,
        )
        for denominator, (amplitude, label, residue) in best.items()
    }


def finite_triples(lows, threshold):
    """All distinct low triples reaching threshold; lows are value-descending."""
    count = len(lows)
    for i in range(count - 2):
        first_value = lows[i][0]
        if first_value + lows[i + 1][0] + lows[i + 2][0] < threshold:
            break
        for j in range(i + 1, count - 1):
            partial = first_value + lows[j][0]
            if partial + lows[j + 1][0] < threshold:
                break
            for k in range(j + 1, count):
                total = partial + lows[k][0]
                if total < threshold:
                    break
                yield tuple(sorted((lows[i][1], lows[j][1], lows[k][1]))), total


def complete_body_cells(carrier, modulus):
    require(RULER % modulus == 0, ("body modulus does not divide ruler", modulus))
    step = RULER // modulus
    result = []
    for left, right in carrier:
        first = (left + step - 1) // step
        last = right // step - 1
        if first <= last:
            result.extend(range(first, last + 1))
    return tuple(result)


def fixed_safe_cells(cells, modulus, labels):
    result = []
    for cell in cells:
        for label in labels:
            residue = (label * cell) % modulus
            # Weak endpoint inequalities are intentional: danger is strict-open.
            if 14 * residue < modulus or 14 * (residue + label) > 13 * modulus:
                break
        else:
            result.append(cell)
    return tuple(result)


def torsion_pairs(cells, denominator):
    by_residue = defaultdict(list)
    for cell in cells:
        by_residue[cell % denominator].append(cell)
    result = []
    for order in range(2, 8):
        if denominator % order:
            continue
        shift = denominator // order
        witness = None
        for residue, first_cells in sorted(by_residue.items()):
            second_cells = by_residue.get((residue + shift) % denominator)
            if second_cells:
                witness = (first_cells[0], second_cells[0])
                break
        if witness is not None:
            result.append((order, shift, witness[0], witness[1]))
    return tuple(result)


def circle_norm(value: F) -> F:
    residue = value - (value.numerator // value.denominator)
    return min(residue, 1 - residue)


def render():
    for path, expected in EXPECTED_HASHES.items():
        require(sha256(path) == expected, ("dependency hash changed", path, sha256(path)))

    rows = parse_atlas()
    require(len(rows) == 58, ("survivor count changed", len(rows)))
    before = tuple(sorted(Counter(row["first"] for row in rows).items()))
    require(before == EXPECTED_BEFORE, ("input height histogram changed", before))
    forced = tuple(row for row in rows if row["high_tail"])
    ordinary = tuple(row for row in rows if not row["high_tail"])
    require((len(forced), len(ordinary)) == (19, 39), "19/39 split changed")

    ordinary_heights = tuple(sorted(Counter(row["first"] for row in ordinary).items()))
    require(ordinary_heights == EXPECTED_AFTER, ("ordinary histogram changed", ordinary_heights))
    require(
        all(
            row["first"] >= row["floor"]
            and all(kind == "EXACT" and label is not None and label > row["first"]
                    for kind, label, _value in row["suffix"])
            for row in ordinary
        ),
        "ordinary all-four-high typing changed",
    )

    analytic_two_high = []
    for row in forced:
        high = [value for kind, _label, value in row["suffix"] if kind == "HIGH-TAIL"]
        exact = sorted(
            (value for kind, _label, value in row["suffix"] if kind == "EXACT"),
            reverse=True,
        )
        require(len(high) == 1 and len(exact) == 3, ("forced certificate shape", row))
        # Duplicate-permitting: two high slots each cost H, and the other two
        # independently cost the top two arbitrary values.  This dominates
        # every literal packet with at least two high labels.
        analytic_gap = row["lower"] - (
            row["first_delta"] + 2 * high[0] + exact[0] + exact[1]
        )
        require(row["first"] < row["floor"], ("forced-high wall inactive", row))
        analytic_two_high.append((analytic_gap, row["first"], row["body"]))
    require(
        sum(gap > 0 for gap, _first, _body in analytic_two_high) == 17
        and tuple((first, body) for gap, first, body in analytic_two_high if gap <= 0)
        == (
            (1722, (2, 8, 9, 10, 12, 14)),
            (1736, (1, 4, 10, 12, 13, 14)),
        ),
        ("analytic two-high hostile ledger changed", analytic_two_high),
    )

    exact_empty = []
    live_rows = []
    classes = []
    exact_two_high = []
    low2_high_margins = []
    for row in forced:
        body = row["body"]
        modulus = 14 * lcm(*body)
        require(modulus == row["L"], ("body modulus changed", body, modulus, row["L"]))
        carrier, mass_numerator, amplitudes = amplitude_data(body, modulus)
        require(F(mass_numerator, RULER) == row["h"], ("carrier mass changed", body))
        require(row["lower"] == row["h"] / 91, ("k=2 lower changed", body))
        require(delta(amplitudes, row["first"]) == row["first_delta"], ("first delta changed", body))

        low_rows = sorted(
            (
                (delta(amplitudes, label), label)
                for label in range(row["first"] + 1, row["floor"])
            ),
            key=lambda item: (-item[0], item[1]),
        )
        required = row["lower"] - row["first_delta"]
        maxima = high_ray_maxima(amplitudes, modulus, row["floor"])
        best_denominator, best_high = max(maxima.items(), key=lambda item: item[1][0])
        # This is the load-bearing at-most-one-high estimate.  Every positive
        # unit ray decreases from its first label at/above the wall, every
        # nonpositive ray is bounded by zero, and all residues mod L were
        # scanned.  For two, three, or four high labels the universal
        # duplicate-permitting upper invoice is
        #   2H + max(H,l1) + max(H,l2).
        # On this exact 19-row universe l2>H, so it specializes to the older
        # 2H+l1+l2 invoice; we assert that specialization rather than assume it.
        low2_high_margins.append(low_rows[1][0] - best_high[0])
        exact_two_upper = (
            2 * best_high[0]
            + max(best_high[0], low_rows[0][0])
            + max(best_high[0], low_rows[1][0])
        )
        exact_two_gap = required - exact_two_upper
        require(exact_two_gap > 0, ("exact-ray two-high survived", row, exact_two_gap))
        exact_two_high.append((exact_two_gap, row["first"], body, best_denominator, best_high[1]))
        exact_gap = required - sum((value for value, _label in low_rows[:3]), F()) - best_high[0]
        if exact_gap > 0:
            exact_empty.append(
                (
                    exact_gap,
                    row["first"],
                    body,
                    tuple(label for _value, label in low_rows[:3]),
                    best_denominator,
                    best_high[1],
                )
            )
            continue

        row_class_count = 0
        cells = complete_body_cells(carrier, modulus)
        for denominator, high in sorted(maxima.items()):
            high_value, high_label, high_residue, high_amplitude = high
            for low_labels, low_total in finite_triples(low_rows, required - high_value):
                require(high_label is not None and high_amplitude > 0, "unattained zero high entered")
                fixed_labels = (row["first"], *low_labels)
                safe = fixed_safe_cells(cells, modulus, fixed_labels)
                torsions = torsion_pairs(safe, denominator)
                require(torsions, ("located-cell failure", row["first"], body, denominator, low_labels))
                unit = high_residue // (modulus // denominator)
                require(gcd(unit, denominator) == 1, ("ray unit changed", denominator, unit))
                for order, shift, first_cell, second_cell in torsions:
                    require(
                        (second_cell - first_cell) % denominator == shift
                        and order * shift == denominator
                        and gcd(unit, order) == 1,
                        ("torsion typing failure", denominator, unit, torsions),
                    )
                excess = row["first_delta"] + low_total + high_value - row["lower"]
                classes.append(
                    (
                        row["first"],
                        body,
                        denominator,
                        low_labels,
                        high_label,
                        unit,
                        excess,
                        len(safe),
                        torsions,
                    )
                )
                row_class_count += 1
        require(row_class_count > 0, ("nonempty exact row had no class", row["first"], body))
        live_rows.append((row["first"], body, row_class_count))

    require(len(exact_empty) == 16 and len(live_rows) == 3, "16/3 exact-ray split changed")
    require(min(exact_two_high)[0] == EXPECTED_TWO_HIGH_MIN, "exact two-high minimum changed")
    require(
        min(low2_high_margins) == EXPECTED_LOW2_HIGH_MARGIN > 0,
        "two-high invoice specialization changed",
    )
    require(min(exact_empty)[0] == EXPECTED_CLOSEST_EMPTY_GAP, "closest exact-empty gap changed")
    require(
        tuple(live_rows)
        == (
            (1708, (2, 8, 9, 10, 12, 14), 2),
            (1722, (2, 8, 9, 10, 12, 14), 10),
            (1736, (1, 8, 10, 12, 13, 14), 1),
        ),
        ("live row ledger changed", live_rows),
    )
    require(len(classes) == EXPECTED_CLASS_COUNT, ("class count changed", len(classes)))
    order_histogram = Counter(order for row in classes for order, _shift, _j, _k in row[-1])
    require(
        order_histogram == Counter({2: 13, 3: 12, 4: 13, 5: 1, 6: 12, 7: 13}),
        ("torsion histogram changed", order_histogram),
    )

    # Sharp strict-open endpoint and first hostile order.
    q7_first = circle_norm(F(16, 1) * F(13, 224))
    q7_second = circle_norm(F(16, 1) * F(29, 224))
    require(q7_first == q7_second == F(1, 14), "q=7 endpoint control changed")
    q8_first = circle_norm(F(63, 1) * F(5, 336))
    q8_second = circle_norm(F(63, 1) * F(11, 336))
    require(q8_first == q8_second == F(1, 16) < F(1, 14), "q=8 hostile changed")
    overlap_q8 = F(1, 7) - F(1, 8)
    require(overlap_q8 == F(1, 56), "q=8 overlap changed")

    semantic_hash = hashlib.sha256(repr(tuple(classes)).encode()).hexdigest()
    if EXPECTED_CLASS_SEMANTIC_SHA256 is not None:
        require(semantic_hash == EXPECTED_CLASS_SEMANTIC_SHA256, "class digest changed")

    lines = [
        "LRC14 projected k=2 forced-high exact-ray/located-torsion closure",
        *(
            f"dependency={path.relative_to(ROOT).as_posix()};sha256={expected}"
            for path, expected in EXPECTED_HASHES.items()
        ),
        f"input_rows=58;height_histogram={before}",
        f"forced_high_rows=19;height_histogram={tuple(sorted(Counter(row['first'] for row in forced).items()))}",
        f"ordinary_all_four_high_rows=39;height_histogram={ordinary_heights};next_frontier_height=1736;lowest_remaining_height=1694",
        "ordinary_status=38 rows L=11760/floor=1020;1 row L=17640/floor=1529;all suffix certificates EXACT;all four later labels exceed z1>=floor",
        "analytic_two_high_control=17 closed/2 open;open rows=z1722:E(2,8,9,10,12,14),z1736:E(1,4,10,12,13,14);the atlas HIGH-TAIL bound is not used for exact-one-high",
        f"exact_two_high_invoice=2H+max(H,l1)+max(H,l2);min_l2_minus_H={ftext(min(low2_high_margins))}",
        f"exact_ray_two_high_gap_min={ftext(min(exact_two_high)[0])};witness=z{min(exact_two_high)[1]}:E{min(exact_two_high)[2]};exactly_one_high=19/19",
        f"exact_ray_empty_rows=16;closest_gap={ftext(min(exact_empty)[0])};witness=z{min(exact_empty)[1]}:E{min(exact_empty)[2]}",
        f"exact_ray_live_rows={tuple(live_rows)};literal_upper_classes={len(classes)};distinct_row_low_packets={len(set((row[0], row[1], row[3]) for row in classes))}",
    ]
    for row in sorted(ordinary, key=lambda item: (item["first"], item["body"])):
        lines.append(
            f"ordinary=z:{row['first']};E:{','.join(map(str,row['body']))};L:{row['L']};floor:{row['floor']};status:ALL-FOUR-HIGH-EXACT"
        )
    for exact_gap, first, body, low_labels, denominator, high_label in sorted(exact_empty):
        lines.append(
            f"exact_empty=z:{first};E:{','.join(map(str,body))};gap:{ftext(exact_gap)};top_lows:{low_labels};best_high_d:{denominator};best_high_label:{high_label}"
        )
    for index, row in enumerate(classes, 1):
        first, body, denominator, lows, high_label, unit, excess, safe_count, torsions = row
        pair_by_order = {order: (j, k) for order, _shift, j, k in torsions}
        lines.append(
            f"class={index};z:{first};E:{','.join(map(str,body))};d:{denominator};lows:{lows};"
            f"high_rep:{high_label};unit:{unit};scalar_excess:{ftext(excess)};safe_cells:{safe_count};"
            f"orders:{tuple(order for order,_shift,_j,_k in torsions)};q2:{pair_by_order[2]};q7:{pair_by_order[7]}"
        )
    lines.extend(
        (
            f"torsion_order_histogram={dict(sorted(order_histogram.items()))};pair_failures=0;q2_closes=13/13;q7_closes=13/13",
            "q7_strict_open=L=14,d=7,u=1,m=1,z=16,j=0,k=1,theta=13/16;arguments have norm exactly 1/14, so closures meet and both open dangers exclude the seam",
            "q8_hostile=L=56,d=8,u=1,m=1,z=63,j=0,k=1,theta=5/6;arguments have norm 1/16<1/14 and both dangers fire;arc-overlap=1/56",
            "first_failed_global_implication=for 39 rows z1>=high_floor, so ordering makes every later label high and the finite-three-low packet does not exist",
            "second_failure_boundary=the analytic HIGH-TAIL placeholder is only an upper bound, not an attained unit ray; exact ray maxima empty 16 of 19 forced-high rows",
            "reversed_gap_near_miss=atlas gap is upper-minus-lower; treating it as lower-minus-upper falsely closes the two analytic hostiles; exact unit-ray maxima must precede the at-most-one-high conclusion",
            "load_bearing_geometric_hypothesis=both located cells are complete body cells and the entire fixed packet is safe on each whole cell, not merely at a root/sample",
            f"class_semantic_sha256={semantic_hash}",
            "all_exact_controls=PASS",
        )
    )
    return "\n".join(lines) + "\n"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    payload = render()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(payload, encoding="utf-8", newline="\n")
    print(payload, end="")


if __name__ == "__main__":
    main()
