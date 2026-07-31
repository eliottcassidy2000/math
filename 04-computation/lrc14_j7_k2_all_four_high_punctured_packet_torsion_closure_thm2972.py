#!/usr/bin/env python3
"""Canonical exact verifier for every k=2 all-four-high literal ray packet.

It reconstructs the 39 all-high rows
from the scalar atlas, proves a finite candidate cutoff, enumerates every
four-label suffix reaching the scalar lower wall, and punctures each resulting
five-label packet at each possible distinguished label.  The retained four
labels are frozen on complete carrier cells and the located torsion-pair test
is applied at every order q=2,...,7 dividing the distinguished denominator.
"""

from __future__ import annotations

import hashlib
import argparse
from collections import Counter
from fractions import Fraction as F
from math import gcd, lcm
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
ATLAS_SOURCE = ROOT / "04-computation/lrc14_j7_k2_scalar_band_1680_1742_thm2941.py"
ATLAS_OUTPUT = ROOT / "05-knowledge/results/lrc14_j7_k2_scalar_band_1680_1742_thm2941.out"
EXPECTED_DEPENDENCIES = {
    ATLAS_SOURCE: "89016f939c961fa979ec5b30812981456df5bfb2af3066f1f1b38e5a83f1a412",
    ATLAS_OUTPUT: "4a36611b26585964e185bbaa3d583be3f1c67a7b608cca785920266bc217a779",
}

RULER = 14 * lcm(*range(1, 15))
ONE_FOURTEENTH = RULER // 14
THIRTEEN_FOURTEENTHS = 13 * ONE_FOURTEENTH
CLASS_MASK_CACHE = {}
EXPECTED_ROW_DIGEST = "c5bbb617c3f1cd15fce322751a0d4857adf731545a0408cd07fb7f5b97a7cf55"
EXPECTED_PACKET_DIGEST = "774e915486b0284718808af4dfb8ef59ddfdd3a7080aea38b2e1ef92e2a90236"
EXPECTED_WITNESS_DIGEST = "4d8e5ba05dd1f8effb1e61b544f054df13be89a1aed730d402b6b1f2843fb7b2"
EXPECTED_PARTIAL_DIGEST = "ef3d7d97893b939c95dc61a9d1fc94fa257f122b1748c0f92fb3054e0a083fe6"
EXPECTED_EMPTY_DIGEST = "01ba4719c80b6fe911b091a7c05124b64eeece964e09c058ef8f9805daca546b"


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def lf_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


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
            suffix.append((kind, None if label == "None" else int(label), F(value)))
        rows.append(
            {
                "body": tuple(map(int, fields["E"].split(","))),
                "L": int(fields["L"]),
                "floor": int(fields["largest_floor"]),
                "first": int(fields["z1"]),
                "first_delta": F(fields["delta1"]),
                "suffix": tuple(suffix),
                "lower": F(fields["lower"]),
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
        + tuple((index * step - radius, index * step + radius) for index in range(1, speed))
        + ((RULER - radius, RULER),)
    )


def carrier_for(body):
    danger = merge_intervals(
        interval for speed in body for interval in danger_intervals(speed)
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
    require(int(np.max(np.abs(amplitudes))) * modulus < 2**63, (body, modulus))
    return carrier, mass_numerator, amplitudes


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


def circle_norm(value: F) -> F:
    residue = value - (value.numerator // value.denominator)
    return min(residue, 1 - residue)


def first_after(residue: int, first: int, period: int) -> int:
    if residue > first:
        return residue
    return residue + ((first - residue) // period + 1) * period


def ray_value(raw: int, label: int) -> F:
    return F(raw, 7 * RULER * label)


def finite_candidates(row):
    """Return the exhaustive positive-ray candidate set for a four-suffix."""
    period = row["L"]
    _carrier, _mass, amplitudes = amplitude_data(row["body"], period)

    # The global top three occur among the first four representatives of every
    # residue because a positive ray is strictly decreasing along z -> z+L.
    heads = []
    for residue, raw_value in enumerate(amplitudes, 1):
        raw = int(raw_value)
        if raw <= 0:
            continue
        label = first_after(residue, row["first"], period)
        for offset in range(4):
            candidate = label + offset * period
            heads.append((ray_value(raw, candidate), candidate, residue, raw))
    heads.sort(key=lambda item: (-item[0], item[1]))
    require(len(heads) >= 3, (row["first"], row["body"], "fewer than three heads"))
    top3 = sum((item[0] for item in heads[:3]), F())
    required = row["lower"] - row["first_delta"]
    threshold = required - top3
    require(threshold > 0, (row["first"], row["body"], threshold))

    candidates = []
    first_omitted = []
    for residue, raw_value in enumerate(amplitudes, 1):
        raw = int(raw_value)
        if raw <= 0:
            continue
        label = first_after(residue, row["first"], period)
        while ray_value(raw, label) >= threshold:
            candidates.append((ray_value(raw, label), label, residue, raw))
            label += period
        first_omitted.append((ray_value(raw, label), label))
    candidates.sort(key=lambda item: (-item[0], item[1]))
    require(all(value < threshold for value, _label in first_omitted), row)
    require(len({label for _value, label, _residue, _raw in candidates}) == len(candidates), row)
    return required, threshold, tuple(candidates), tuple(heads[:3])


def admissible_suffixes(candidates, required: F):
    """Yield every distinct four-candidate suffix with exact sum >= required."""
    values = tuple(item[0] for item in candidates)
    labels = tuple(item[1] for item in candidates)
    n = len(values)
    for i in range(n - 3):
        if sum(values[i : i + 4], F()) < required:
            break
        for j in range(i + 1, n - 2):
            if values[i] + sum(values[j : j + 3], F()) < required:
                break
            for k in range(j + 1, n - 1):
                partial = values[i] + values[j] + values[k]
                if partial + values[k + 1] < required:
                    break
                for ell in range(k + 1, n):
                    total = partial + values[ell]
                    if total < required:
                        break
                    suffix = tuple(sorted((labels[i], labels[j], labels[k], labels[ell])))
                    require(len(set(suffix)) == 4, suffix)
                    yield suffix, total


def denominator_and_unit(label: int, modulus: int) -> tuple[int, int]:
    residue = label % modulus
    divisor = gcd(modulus, residue)
    denominator = modulus // divisor
    unit = residue // divisor
    require(denominator > 1 and gcd(unit, denominator) == 1, (label, modulus))
    return denominator, unit


def safe_mask(cells, modulus: int, label: int) -> int:
    return sum(1 << cell for cell in fixed_safe_cells(cells, modulus, (label,)))


def residue_presence(mask: int, modulus: int, denominator: int) -> int:
    """Fold a modulus-cell bitmap onto residues modulo denominator."""
    require(modulus % denominator == 0, (modulus, denominator))
    quotient = modulus // denominator
    if denominator <= quotient:
        key = (modulus, denominator)
        classes = CLASS_MASK_CACHE.get(key)
        if classes is None:
            rows = [0] * denominator
            for cell in range(modulus):
                rows[cell % denominator] |= 1 << cell
            classes = tuple(rows)
            CLASS_MASK_CACHE[key] = classes
        result = 0
        for residue, class_mask in enumerate(classes):
            if mask & class_mask:
                result |= 1 << residue
        return result
    chunk_mask = (1 << denominator) - 1
    result = 0
    for offset in range(0, modulus, denominator):
        result |= (mask >> offset) & chunk_mask
    return result


def rotate_right(mask: int, shift: int, width: int) -> int:
    low = mask & ((1 << shift) - 1)
    return (mask >> shift) | (low << (width - shift))


def fast_torsion_orders(mask: int, modulus: int, denominator: int):
    """Return (q,shift,residue) for every located q<=7 torsion edge."""
    residues = residue_presence(mask, modulus, denominator)
    result = []
    for order in range(2, 8):
        if denominator % order:
            continue
        shift = denominator // order
        paired = residues & rotate_right(residues, shift, denominator)
        if paired:
            lowest = paired & -paired
            residue = lowest.bit_length() - 1
            result.append((order, shift, residue))
    return tuple(result)


def first_cell(mask: int, residue: int, denominator: int, modulus: int) -> int:
    for cell in range(residue, modulus, denominator):
        if (mask >> cell) & 1:
            return cell
    raise RuntimeError(("missing residue", residue, denominator))


def expand_torsions(mask: int, modulus: int, denominator: int, torsions):
    result = []
    for order, shift, residue in torsions:
        second = (residue + shift) % denominator
        result.append(
            (
                order,
                shift,
                first_cell(mask, residue, denominator, modulus),
                first_cell(mask, second, denominator, modulus),
            )
        )
    return tuple(result)


def semantic_digest(lines) -> str:
    payload = "\n".join(lines).encode("ascii") + b"\n"
    return hashlib.sha256(payload).hexdigest()


def mask_cells(mask: int):
    result = []
    while mask:
        low = mask & -mask
        result.append(low.bit_length() - 1)
        mask -= low
    return tuple(result)


def slow_torsion_orders(mask: int, denominator: int):
    residues = {cell % denominator for cell in mask_cells(mask)}
    return tuple(
        order
        for order in range(2, 8)
        if denominator % order == 0
        and any((residue + denominator // order) % denominator in residues for residue in residues)
    )


def render() -> str:
    for path, expected in EXPECTED_DEPENDENCIES.items():
        require(lf_sha256(path) == expected, ("dependency changed", path, lf_sha256(path)))
    rows = tuple(row for row in parse_atlas() if not row["high_tail"])
    require(len(rows) == 39, ("ordinary-row universe changed", len(rows)))

    row_records = []
    packet_records = []
    witness_records = []
    failure_records = []
    partial_records = []
    packet_count = 0
    killed_count = 0
    all_five_count = 0
    distinguished_count_hist = Counter()
    chosen_index_hist = Counter()
    chosen_order_hist = Counter()
    all_witness_order_hist = Counter()
    row_packet_hist = Counter()
    row_failure_hist = Counter()
    height_packet_hist = Counter()
    minimum_safe = None
    maximum_safe = 0
    chosen_safe_minimum = None
    chosen_safe_maximum = 0
    total_candidates = 0

    for row in rows:
        z1 = row["first"]
        body = row["body"]
        modulus = row["L"]
        required, threshold, candidates, top3 = finite_candidates(row)
        total_candidates += len(candidates)
        cells = complete_body_cells(carrier_for(body), modulus)
        labels_needed = {z1}
        labels_needed.update(item[1] for item in candidates)
        safe_by_label = {label: safe_mask(cells, modulus, label) for label in labels_needed}

        local_packets = 0
        local_failures = 0
        for suffix, suffix_total in admissible_suffixes(candidates, required):
            packet = (z1, *suffix)
            require(len(set(packet)) == 5, packet)
            require(suffix_total >= required, packet)
            packet_count += 1
            local_packets += 1
            packet_records.append(
                f"z={z1};E={','.join(map(str, body))};packet={','.join(map(str, packet))}"
            )

            witnesses = []
            for index, label in enumerate(packet):
                fixed = packet[:index] + packet[index + 1 :]
                mask = safe_by_label[fixed[0]]
                for item in fixed[1:]:
                    mask &= safe_by_label[item]
                safe_count = mask.bit_count()
                minimum_safe = safe_count if minimum_safe is None else min(minimum_safe, safe_count)
                maximum_safe = max(maximum_safe, safe_count)
                denominator, unit = denominator_and_unit(label, modulus)
                compressed_torsions = fast_torsion_orders(mask, modulus, denominator)
                torsions = expand_torsions(mask, modulus, denominator, compressed_torsions)
                for order, shift, first_cell_value, second_cell_value in torsions:
                    require(
                        ((mask >> first_cell_value) & 1)
                        and ((mask >> second_cell_value) & 1)
                        and (second_cell_value - first_cell_value) % denominator == shift
                        and order * shift == denominator
                        and gcd(unit, order) == 1,
                        ("torsion typing failure", packet, index, torsions),
                    )
                if local_packets == 1:
                    require(
                        tuple(item[0] for item in torsions) == slow_torsion_orders(mask, denominator),
                        ("fast/slow torsion mismatch", packet, index),
                    )
                if torsions:
                    witnesses.append((index, label, denominator, unit, safe_count, torsions))
                    all_witness_order_hist.update(order for order, _shift, _j, _k in torsions)

            distinguished_count_hist[len(witnesses)] += 1
            if len(witnesses) < 5:
                live_indices = {item[0] for item in witnesses}
                partial_records.append(
                    f"z={z1};E={','.join(map(str, body))};packet={','.join(map(str, packet))};"
                    f"working={','.join(map(str, sorted(live_indices)))};"
                    f"missing={','.join(str(index) for index in range(5) if index not in live_indices)}"
                )
            if not witnesses:
                local_failures += 1
                row_failure_hist[(z1, body)] += 1
                failure_records.append(
                    f"z={z1};E={','.join(map(str, body))};packet={','.join(map(str, packet))}"
                )
                continue

            killed_count += 1
            if len(witnesses) == 5:
                all_five_count += 1
            chosen = min(
                witnesses,
                key=lambda item: (item[5][0][0], -item[4], item[0], item[1]),
            )
            index, label, denominator, unit, safe_count, torsions = chosen
            chosen_safe_minimum = safe_count if chosen_safe_minimum is None else min(chosen_safe_minimum, safe_count)
            chosen_safe_maximum = max(chosen_safe_maximum, safe_count)
            chosen_index_hist[index] += 1
            chosen_order_hist[torsions[0][0]] += 1
            witness_records.append(
                "K;"
                f"z={z1};E={','.join(map(str, body))};packet={','.join(map(str, packet))};"
                f"i={index};label={label};d={denominator};u={unit};safe={safe_count};"
                f"torsion="
                + ",".join(f"{q}:{shift}:{j}:{k}" for q, shift, j, k in torsions)
            )

        row_packet_hist[(z1, body)] = local_packets
        height_packet_hist[z1] += local_packets
        row_records.append(
            f"z={z1};E={','.join(map(str, body))};L={modulus};candidates={len(candidates)};"
            f"packets={local_packets};failures={local_failures};threshold={threshold};"
            f"top3={','.join(str(item[1]) for item in top3)}"
        )

    require(packet_count == 19_640, ("packet universe changed", packet_count))
    require(killed_count == packet_count and not failure_records, (killed_count, failure_records))
    require(chosen_order_hist == Counter({2: packet_count}), chosen_order_hist)
    require(total_candidates == 1_401, ("candidate universe changed", total_candidates))
    require(len(partial_records) == 10, ("partial-puncture boundary changed", partial_records))

    row_digest = semantic_digest(row_records)
    packet_digest = semantic_digest(packet_records)
    witness_digest = semantic_digest(witness_records)
    partial_digest = semantic_digest(partial_records)
    failure_digest = semantic_digest(failure_records)
    require(row_digest == EXPECTED_ROW_DIGEST, ("row digest changed", row_digest))
    require(packet_digest == EXPECTED_PACKET_DIGEST, ("packet digest changed", packet_digest))
    require(witness_digest == EXPECTED_WITNESS_DIGEST, ("witness digest changed", witness_digest))
    require(partial_digest == EXPECTED_PARTIAL_DIGEST, ("partial digest changed", partial_digest))
    require(failure_digest == EXPECTED_EMPTY_DIGEST, ("failure digest changed", failure_digest))

    # The strict-open q=7 seam and the first failed universal order q=8.
    q7_first = circle_norm(F(16) * F(13, 224))
    q7_second = circle_norm(F(16) * F(29, 224))
    require(q7_first == q7_second == F(1, 14), (q7_first, q7_second))
    q8_first = circle_norm(F(63) * F(5, 336))
    q8_second = circle_norm(F(63) * F(11, 336))
    q8_overlap = F(1, 7) - F(1, 8)
    require(
        q8_first == q8_second == F(1, 16) < F(1, 14) and q8_overlap == F(1, 56),
        (q8_first, q8_second, q8_overlap),
    )

    lines = [
        "LRC14 k=2 all-four-high literal-ray packet torsion closure",
        *(f"dependency={path.relative_to(ROOT).as_posix()};lf_sha256={expected}" for path, expected in EXPECTED_DEPENDENCIES.items()),
        f"rows={len(rows)};packets={packet_count};killed={killed_count};open={len(failure_records)}",
        f"finite_candidates_total={total_candidates};maximum_row_candidates={max(int(line.split('candidates=')[1].split(';')[0]) for line in row_records)}",
        f"row_height_histogram={dict(sorted(Counter(row['first'] for row in rows).items()))}",
        f"packet_height_histogram={dict(sorted(height_packet_hist.items()))}",
        f"below_1736_rows={sum(1 for row in rows if row['first'] < 1736)};below_1736_packets={sum(count for height,count in height_packet_hist.items() if height < 1736)};all_killed=yes",
        f"all_five_distinguished_work={all_five_count}/{packet_count}",
        f"distinguished_count_histogram={dict(sorted(distinguished_count_hist.items()))}",
        f"chosen_index_histogram={dict(sorted(chosen_index_hist.items()))}",
        f"chosen_first_order_histogram={dict(sorted(chosen_order_hist.items()))}",
        f"all_witness_order_histogram={dict(sorted(all_witness_order_hist.items()))}",
        f"safe_cell_range={minimum_safe}..{maximum_safe}",
        f"chosen_safe_cell_range={chosen_safe_minimum}..{chosen_safe_maximum}",
        f"row_semantic_sha256={row_digest}",
        f"packet_semantic_sha256={packet_digest}",
        f"witness_semantic_sha256={witness_digest}",
        f"partial_semantic_sha256={partial_digest};partial_packets={len(partial_records)}",
        f"failure_semantic_sha256={failure_digest}",
        "q7_strict_open=L:14;d:7;u:1;m:1;z:16;j:0;k:1;theta:13/16;norms:1/14,1/14;open dangers disjoint",
        "q8_hostile=L:56;d:8;u:1;m:1;z:63;j:0;k:1;theta:5/6;norms:1/16,1/16;danger overlap:1/56",
        "row_census_begin",
        *row_records,
        "row_census_end",
    ]
    if partial_records:
        lines.extend(("partial_packets_begin", *partial_records, "partial_packets_end"))
    if failure_records:
        lines.extend(("failures_begin", *failure_records, "failures_end"))
    lines.extend(
        (
            "scope=39 exact-suffix rows of the projected k=2 1680..1742 atlas;all distinct later labels meeting the exact scalar wall;no HIGH-TAIL rows",
            "consequence=every exact-suffix row is empty;combined with a separate closure of the 19 HIGH-TAIL rows, the full 1680..1742 band is empty",
            "all_exact_controls=PASS",
        )
    )
    return "\n".join(lines) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output",
        type=Path,
        default=(
            ROOT
            / "05-knowledge/results"
            / "lrc14_j7_k2_all_four_high_punctured_packet_torsion_closure_thm2972.out"
        ),
    )
    args = parser.parse_args()
    payload = render()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(payload, encoding="utf-8", newline="\n")
    print(payload, end="")


if __name__ == "__main__":
    main()

