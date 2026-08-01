#!/usr/bin/env python3
"""Exact projected-k2 closure on 1600 <= z1 <= 1679.

The verifier separates finite packet cones, zero-high cones, and projective
exactly-one-high rays.  Small torsion closes the ordinary cases.  Every
remaining projective class has full cyclic residue projection, yielding a
lawful unit-dependent high-order two-cell chord.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import sys
from collections import Counter
from fractions import Fraction as F
from math import gcd
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
UTILITIES = ROOT / "04-computation/lrc14_j7_k2_all_four_high_punctured_packet_torsion_closure_thm2972.py"
ATLAS_SOURCE = ROOT / "04-computation/lrc14_j7_k2_scalar_band_1600_1679_thm2980.py"
ATLAS_OUTPUT = ROOT / "05-knowledge/results/lrc14_j7_k2_scalar_band_1600_1679_thm2980.out"
OUTPUT = ROOT / "05-knowledge/results/lrc14_j7_k2_1600_1679_zero_high_and_high_order_phase_closure_thm2980.out"

EXPECTED_UTILITIES_SHA256 = "b92f7c6cf295f29a3ec2c64730f6bd4486073b9197f322202df98df497369d97"
EXPECTED_ATLAS_SOURCE_SHA256 = "ae6610e52bb18057e7c15f8599141c0a738419dc7d364d56b6688d652e466922"
EXPECTED_ATLAS_OUTPUT_SHA256 = "86967e15677757e71d4962f5b369afcd40fc95ea09ff1c5a8585e3c6dbe55964"
EXPECTED_EMPTY_SHA256 = "01ba4719c80b6fe911b091a7c05124b64eeece964e09c058ef8f9805daca546b"
EXPECTED_FINITE_PACKET_SHA256 = "fead04db235721418dc93d36e7e1894a88a40fec36e584037ea089040e64e028"
EXPECTED_FINITE_WITNESS_SHA256 = "e1ecfb3b254e93dc044d28e6baf399f74a3f7a7c2ea4a679c90d225a8a93abc0"

EXCEPTIONAL_KEYS = {
    (1612, (1, 8, 10, 12, 13, 14)),
    (1612, (1, 10, 11, 12, 13, 14)),
    (1650, (1, 10, 11, 12, 13, 14)),
}
EXPECTED_FINITE_PACKETS = 69_602
EXPECTED_ZERO_PACKETS = 1_098_353
EXPECTED_PROJECTIVE_CLASSES = 6_934
EXPECTED_SMALL_TORSION_CLASSES = 6_884
EXPECTED_FULL_PROJECTION_CLASSES = 50
EXPECTED_UNIT_TESTS = 754
EXPECTED_NONPOSITIVE_TRIPLES = 27
EXPECTED_NONPOSITIVE_RESIDUE_RAYS = 496_866
EXPECTED_NONPOSITIVE_UNIT_TESTS = 5_159_799
EXPECTED_NONPOSITIVE_METHODS = Counter({"small_torsion": 5_159_052, "full_projection": 747})
EXPECTED_NONPOSITIVE_MINIMUM_FIXED_CELLS = 24_348


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def lf_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def load_utilities():
    require(lf_sha256(UTILITIES) == EXPECTED_UTILITIES_SHA256, ("utilities changed", lf_sha256(UTILITIES)))
    spec = importlib.util.spec_from_file_location("thm2980_torsion_utilities", UTILITIES)
    require(spec is not None and spec.loader is not None, UTILITIES)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


base = load_utilities()
RULER = base.RULER
AMPLITUDE_CACHE = {}
TYPE_AUDIT = {
    "max_primitive_input": 0,
    "max_amplitude_times_modulus": 0,
    "max_label_times_cell": 0,
    "max_safe_linear_form": 0,
    "safe_full_cell_checks": 0,
    "fold_full_cell_checks": 0,
}


def semantic_digest(lines) -> str:
    return hashlib.sha256(("\n".join(lines) + "\n").encode("ascii")).hexdigest()


def parse_rows():
    require(lf_sha256(ATLAS_SOURCE) == EXPECTED_ATLAS_SOURCE_SHA256, ("atlas source changed", lf_sha256(ATLAS_SOURCE)))
    require(lf_sha256(ATLAS_OUTPUT) == EXPECTED_ATLAS_OUTPUT_SHA256, ("atlas output changed", lf_sha256(ATLAS_OUTPUT)))
    rows = []
    for line in ATLAS_OUTPUT.read_text(encoding="utf-8-sig").splitlines():
        if not line.startswith("SURVIVOR;"):
            continue
        fields = dict(field.split("=", 1) for field in line.split(";")[1:] if "=" in field)
        rows.append(
            {
                "body": tuple(map(int, fields["E"].split(","))),
                "L": int(fields["L"]),
                "floor": int(fields["largest_floor"]),
                "first": int(fields["z1"]),
                "first_delta": F(fields["delta1"]),
                "lower": F(fields["lower"]),
                "high_tail": "HIGH-TAIL" in fields["suffix"],
            }
        )
    rows.sort(key=lambda row: (row["first"], row["body"]))
    require(len(rows) == 68, ("row universe", len(rows)))
    require(sum(not row["high_tail"] for row in rows) == 40, "exact row count")
    require(sum(row["high_tail"] for row in rows) == 28, "high-tail row count")
    return tuple(rows)


def amplitude(body, modulus):
    key = (body, modulus)
    values = AMPLITUDE_CACHE.get(key)
    if values is None:
        carrier, _mass, values = base.amplitude_data(body, modulus)
        primitive_bound = max((modulus - 1) * right for _left, right in carrier)
        amplitude_bound = int(np.max(np.abs(values))) * modulus
        require(primitive_bound < 2**63 and amplitude_bound < 2**63, (key, primitive_bound, amplitude_bound))
        TYPE_AUDIT["max_primitive_input"] = max(TYPE_AUDIT["max_primitive_input"], primitive_bound)
        TYPE_AUDIT["max_amplitude_times_modulus"] = max(
            TYPE_AUDIT["max_amplitude_times_modulus"], amplitude_bound
        )
        AMPLITUDE_CACHE[key] = values
    return values


def ray_value(raw: int, label: int) -> F:
    return F(raw, 7 * RULER * label)


def first_after(residue: int, first: int, period: int) -> int:
    if residue > first:
        return residue
    return residue + ((first - residue) // period + 1) * period


def first_at_or_after(residue: int, start: int, period: int) -> int:
    if residue >= start:
        return residue
    return residue + ((start - residue + period - 1) // period) * period


def finite_candidates(row):
    """Exact positive-ray cutoff, with a permissive global top-three bound."""
    period = row["L"]
    values = amplitude(row["body"], period)
    heads = []

    def insert_head(raw, label, residue):
        item = (raw, label, residue)
        position = 0
        while position < len(heads) and heads[position][0] * label >= raw * heads[position][1]:
            position += 1
        heads.insert(position, item)
        if len(heads) > 3:
            heads.pop()

    for residue, raw_value in enumerate(values, 1):
        raw = int(raw_value)
        if raw <= 0:
            continue
        label = first_after(residue, row["first"], period)
        for offset in range(3):
            insert_head(raw, label + offset * period, residue)
    exact_heads = tuple((ray_value(raw, label), label, residue, raw) for raw, label, residue in heads)
    required = row["lower"] - row["first_delta"]
    threshold = required - sum((item[0] for item in exact_heads), F())
    if threshold <= 0:
        return required, threshold, (), exact_heads
    candidates = []
    for residue, raw_value in enumerate(values, 1):
        raw = int(raw_value)
        if raw <= 0:
            continue
        label = first_after(residue, row["first"], period)
        while raw * threshold.denominator >= threshold.numerator * 7 * RULER * label:
            candidates.append((ray_value(raw, label), label, residue, raw))
            label += period
    candidates.sort(key=lambda item: (-item[0], item[1]))
    return required, threshold, tuple(candidates), exact_heads


def admissible_suffixes(candidates, required):
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


def finite_triples(rows, threshold):
    n = len(rows)
    for i in range(n - 2):
        if rows[i][0] + rows[i + 1][0] + rows[i + 2][0] < threshold:
            break
        for j in range(i + 1, n - 1):
            partial = rows[i][0] + rows[j][0]
            if partial + rows[j + 1][0] < threshold:
                break
            for k in range(j + 1, n):
                total = partial + rows[k][0]
                if total < threshold:
                    break
                yield (rows[i][1], rows[j][1], rows[k][1]), total


def denominator_and_unit(label: int, modulus: int):
    residue = label % modulus
    divisor = gcd(modulus, residue)
    denominator = modulus // divisor
    unit = residue // divisor
    require(denominator > 1 and gcd(unit, denominator) == 1, (label, modulus, denominator, unit))
    return denominator, unit


def vector_safe_mask(cells_array, modulus, label):
    require(cells_array.dtype == np.int64 and label > 0, (cells_array.dtype, label))
    max_cell = int(cells_array[-1]) if len(cells_array) else 0
    product_bound = label * max_cell
    safe_linear_bound = 14 * (modulus - 1 + label)
    require(
        label < 2**63 - modulus and product_bound < 2**63 and safe_linear_bound < 2**63,
        (modulus, label, max_cell, product_bound, safe_linear_bound),
    )
    TYPE_AUDIT["max_label_times_cell"] = max(TYPE_AUDIT["max_label_times_cell"], product_bound)
    TYPE_AUDIT["max_safe_linear_form"] = max(TYPE_AUDIT["max_safe_linear_form"], safe_linear_bound)
    residues = (label * cells_array) % modulus
    safe = (14 * residues >= modulus) & (14 * (residues + label) <= 13 * modulus)
    bitmap = np.zeros(modulus, dtype=np.uint8)
    bitmap[cells_array[safe]] = 1
    return int.from_bytes(np.packbits(bitmap, bitorder="little").tobytes(), "little")


def vector_residue_presence(mask, modulus, denominator):
    require(modulus % denominator == 0, (modulus, denominator))
    byte_count = (modulus + 7) // 8
    raw = np.frombuffer(mask.to_bytes(byte_count, "little"), dtype=np.uint8)
    bits = np.unpackbits(raw, bitorder="little")[:modulus]
    present = np.any(bits.reshape(modulus // denominator, denominator), axis=0)
    return int.from_bytes(np.packbits(present, bitorder="little").tobytes(), "little")


def scalar_residue_presence(mask, cells, denominator):
    result = 0
    for cell in cells:
        if (mask >> cell) & 1:
            result |= 1 << (cell % denominator)
    return result


def torsion_orders(mask, modulus, denominator):
    residues = vector_residue_presence(mask, modulus, denominator)
    result = []
    for order in range(2, 8):
        if denominator % order:
            continue
        shift = denominator // order
        paired = residues & base.rotate_right(residues, shift, denominator)
        if paired:
            lowest = paired & -paired
            result.append((order, shift, lowest.bit_length() - 1))
    return tuple(result)


def audit_safe_mask(cells, cells_array, modulus, label, vector_mask):
    require(vector_mask == base.safe_mask(cells, modulus, label), ("vector/scalar safe", modulus, label))
    TYPE_AUDIT["safe_full_cell_checks"] += 1


def audit_fold(mask, cells, modulus, denominator):
    vector = vector_residue_presence(mask, modulus, denominator)
    scalar = scalar_residue_presence(mask, cells, denominator)
    require(vector == scalar, ("vector/scalar fold", modulus, denominator))
    TYPE_AUDIT["fold_full_cell_checks"] += 1


def cyclic_span(residues, denominator, unit):
    """Return minimal cyclic covering-arc length and the deleted largest gap."""
    images = sorted({(unit * residue) % denominator for residue in residues})
    require(images, (residues, denominator, unit))
    gaps = [images[index + 1] - images[index] for index in range(len(images) - 1)]
    gaps.append(images[0] + denominator - images[-1])
    largest_gap = max(gaps)
    return denominator - largest_gap, largest_gap


def render() -> str:
    for key in TYPE_AUDIT:
        TYPE_AUDIT[key] = 0
    rows = parse_rows()

    finite_rows = []
    exceptional_rows = []
    finite_packet_count = finite_failure_count = 0
    exact_packet_count = finite_high_packet_count = 0
    finite_packet_hash = hashlib.sha256()
    finite_witness_hash = hashlib.sha256()
    finite_failure_lines = []
    height_packets = Counter()

    for row in rows:
        key = (row["first"], row["body"])
        if key in EXCEPTIONAL_KEYS:
            exceptional_rows.append(row)
            continue
        required, cutoff, candidates, top3 = finite_candidates(row)
        require(cutoff > 0 and len(candidates) <= 5000, ("finite classification", key, cutoff, len(candidates)))
        modulus = row["L"]
        cells = base.complete_body_cells(base.carrier_for(row["body"]), modulus)
        cells_array = np.asarray(cells, dtype=np.int64)
        labels = {row["first"]}
        labels.update(item[1] for item in candidates)
        safe_masks = {label: base.safe_mask(cells, modulus, label) for label in sorted(labels)}
        audit_safe_mask(cells, cells_array, modulus, row["first"], vector_safe_mask(cells_array, modulus, row["first"]))
        local_packets = local_failures = local_q2 = 0
        folded = False
        for suffix, total in admissible_suffixes(candidates, required):
            packet = (row["first"], *suffix)
            require(len(set(packet)) == 5 and total >= required, packet)
            local_packets += 1
            finite_packet_count += 1
            height_packets[row["first"]] += 1
            if row["high_tail"]:
                finite_high_packet_count += 1
            else:
                exact_packet_count += 1
            packet_line = f"{row['first']}|{','.join(map(str,row['body']))}|{','.join(map(str,packet))}\n"
            finite_packet_hash.update(packet_line.encode("ascii"))
            witness = None
            for index, label in enumerate(packet):
                fixed = packet[:index] + packet[index + 1 :]
                mask = safe_masks[fixed[0]]
                for other in fixed[1:]:
                    mask &= safe_masks[other]
                denominator, unit = denominator_and_unit(label, modulus)
                if not folded:
                    audit_fold(mask, cells, modulus, denominator)
                    folded = True
                # The audited THM-2972 bitset fold is faster on finite-row
                # moduli; one full real-bank fold per row is cross-checked
                # independently below before this path is trusted.
                orders = base.fast_torsion_orders(mask, modulus, denominator)
                if orders:
                    witness = (index, label, denominator, unit, mask.bit_count(), orders[0])
                    break
            if witness is None:
                local_failures += 1
                finite_failure_count += 1
                finite_failure_lines.append(packet_line.strip())
            else:
                if witness[-1][0] == 2:
                    local_q2 += 1
                finite_witness_hash.update((packet_line.strip() + "|" + repr(witness) + "\n").encode("ascii"))
        require(local_packets > 0 and local_failures == 0 and folded, (key, local_packets, local_failures, folded))
        finite_rows.append(
            f"z={row['first']};E={','.join(map(str,row['body']))};high={row['high_tail']};"
            f"candidates={len(candidates)};packets={local_packets};failures={local_failures};q2={local_q2};"
            f"cutoff={cutoff};top3={','.join(str(item[1]) for item in top3)}"
        )

    require(len(finite_rows) == 65 and len(exceptional_rows) == 3, (len(finite_rows), len(exceptional_rows)))
    require(
        finite_packet_count == EXPECTED_FINITE_PACKETS
        and finite_failure_count == 0
        and exact_packet_count == 16_809
        and finite_high_packet_count == 52_793,
        (finite_packet_count, finite_failure_count, exact_packet_count, finite_high_packet_count),
    )
    require(finite_packet_hash.hexdigest() == EXPECTED_FINITE_PACKET_SHA256, finite_packet_hash.hexdigest())
    require(finite_witness_hash.hexdigest() == EXPECTED_FINITE_WITNESS_SHA256, finite_witness_hash.hexdigest())
    require(semantic_digest(finite_failure_lines) == EXPECTED_EMPTY_SHA256, finite_failure_lines)

    zero_packet_count = zero_failure_count = 0
    projective_class_count = small_torsion_class_count = 0
    full_projection_class_count = unit_test_count = 0
    zero_packet_hash = hashlib.sha256()
    zero_witness_hash = hashlib.sha256()
    small_class_hash = hashlib.sha256()
    phase_records = []
    exceptional_records = []
    denominator_failures = Counter()
    minimum_full_projection_margin = None
    nonpositive_triple_count = 0
    nonpositive_residue_ray_count = 0
    nonpositive_unit_tests = 0
    nonpositive_failures = []
    nonpositive_method_hist = Counter()
    nonpositive_records = []
    minimum_nonpositive_fixed_cells = None

    for row in exceptional_rows:
        body, modulus, first, floor = row["body"], row["L"], row["first"], row["floor"]
        required = row["lower"] - row["first_delta"]
        values = amplitude(body, modulus)
        lows = []
        for label in range(first + 1, floor):
            raw = int(values[label - 1])
            if raw > 0:
                lows.append((ray_value(raw, label), label, raw))
        lows.sort(key=lambda item: (-item[0], item[1]))
        require(len(lows) >= 4, (first, body))

        maxima = {}
        rays = {}
        high_start = max(first + 1, floor)
        for residue, raw_value in enumerate(values, 1):
            raw = int(raw_value)
            if raw <= 0:
                continue
            denominator = modulus // gcd(modulus, residue)
            divisor = modulus // denominator
            unit = residue // divisor
            require(gcd(unit, denominator) == 1, (denominator, unit))
            label = first_at_or_after(residue, high_start, modulus)
            value = ray_value(raw, label)
            ray = (value, label, unit, residue, raw)
            rays.setdefault(denominator, []).append(ray)
            old = maxima.get(denominator)
            if old is None or value > old[0] or (value == old[0] and label < old[1]):
                maxima[denominator] = (value, label, residue, raw)

        best_d, best_high = max(maxima.items(), key=lambda item: (item[1][0], -item[0]))
        require(best_high[0] <= lows[1][0], (first, body, best_high[0], lows[1][0]))
        two_positive_gap = required - lows[0][0] - lows[1][0]
        three_gap = required - sum((item[0] for item in lows[:3]), F())
        mixed_gap = required - best_high[0] - lows[0][0] - lows[1][0]
        zero_gap = required - sum((item[0] for item in lows[:4]), F())
        two_gap = required - 2 * best_high[0] - lows[0][0] - lows[1][0]
        require(
            zero_gap < 0 and two_positive_gap > 0 and mixed_gap > 0 and two_gap > 0,
            (first, body, zero_gap, two_positive_gap, mixed_gap, two_gap),
        )

        cells = base.complete_body_cells(base.carrier_for(body), modulus)
        cells_array = np.asarray(cells, dtype=np.int64)
        safe_cache = {first: vector_safe_mask(cells_array, modulus, first)}
        audit_safe_mask(cells, cells_array, modulus, first, safe_cache[first])

        # A nonpositive suffix needs exactly three positive low companions.
        # Positive-high + two-low is excluded by mixed_gap>0, and two or more
        # nonpositive entries are excluded because even the best two positive
        # lows are insufficient.  The first two rows have 18 and 9 such low
        # triples; the third has none.
        nonpositive_triples = tuple(finite_triples(lows, required))
        expected_triples = 18 if body == (1, 8, 10, 12, 13, 14) else 9 if first == 1612 else 0
        require(len(nonpositive_triples) == expected_triples, (first, body, len(nonpositive_triples)))
        if nonpositive_triples:
            units_by_denominator = {}
            for residue, raw_value in enumerate(values, 1):
                if int(raw_value) > 0:
                    continue
                denominator = modulus // gcd(modulus, residue)
                divisor = modulus // denominator
                unit = residue // divisor
                require(denominator > 1 and gcd(unit, denominator) == 1, (residue, denominator, unit))
                units_by_denominator.setdefault(denominator, set()).add(unit)
            nonpositive_residue_ray_count += sum(len(units) for units in units_by_denominator.values())
            local_nonpositive_tests = local_nonpositive_failures = 0
            local_nonpositive_methods = Counter()
            local_minimum_cells = None
            for labels, low_total in nonpositive_triples:
                nonpositive_triple_count += 1
                mask = safe_cache[first]
                for label in labels:
                    if label not in safe_cache:
                        safe_cache[label] = vector_safe_mask(cells_array, modulus, label)
                    mask &= safe_cache[label]
                fixed_cells = mask.bit_count()
                require(fixed_cells > 0 and low_total >= required, (first, body, labels, fixed_cells, low_total - required))
                local_minimum_cells = fixed_cells if local_minimum_cells is None else min(local_minimum_cells, fixed_cells)
                minimum_nonpositive_fixed_cells = (
                    fixed_cells
                    if minimum_nonpositive_fixed_cells is None
                    else min(minimum_nonpositive_fixed_cells, fixed_cells)
                )
                for denominator, units in sorted(units_by_denominator.items()):
                    unit_count = len(units)
                    local_nonpositive_tests += unit_count
                    nonpositive_unit_tests += unit_count
                    residue_mask = vector_residue_presence(mask, modulus, denominator)
                    has_small_torsion = any(
                        denominator % order == 0
                        and residue_mask
                        & base.rotate_right(residue_mask, denominator // order, denominator)
                        for order in range(2, 8)
                    )
                    if has_small_torsion:
                        local_nonpositive_methods["small_torsion"] += unit_count
                        nonpositive_method_hist["small_torsion"] += unit_count
                        continue
                    if residue_mask == (1 << denominator) - 1:
                        local_nonpositive_methods["full_projection"] += unit_count
                        nonpositive_method_hist["full_projection"] += unit_count
                        continue
                    residues = tuple(index for index in range(denominator) if (residue_mask >> index) & 1)
                    for unit in sorted(units):
                        span, largest_gap = cyclic_span(residues, denominator, unit)
                        if 7 * span >= denominator:
                            local_nonpositive_methods["diameter"] += 1
                            nonpositive_method_hist["diameter"] += 1
                        else:
                            local_nonpositive_failures += 1
                            nonpositive_failures.append(
                                (first, body, labels, denominator, unit, fixed_cells, span, largest_gap)
                            )
            require(local_nonpositive_failures == 0, (first, body, local_nonpositive_failures))
            nonpositive_records.append(
                f"z={first};E={','.join(map(str,body))};triples={len(nonpositive_triples)};"
                f"nonpositive_residue_rays={sum(len(units) for units in units_by_denominator.values())};"
                f"denominators={len(units_by_denominator)};unit_tests={local_nonpositive_tests};"
                f"min_fixed_cells={local_minimum_cells};methods={dict(sorted(local_nonpositive_methods.items()))};"
                "first_representative_failures=0;later_height_safe_mass_lower=3/7"
            )

        low_candidates = tuple((value, label, label % modulus, raw) for value, label, raw in lows)
        local_zero = local_zero_failures = 0
        folded = False
        for suffix, total in admissible_suffixes(low_candidates, required):
            packet = (first, *suffix)
            local_zero += 1
            zero_packet_count += 1
            packet_line = f"{first}|{','.join(map(str,body))}|{','.join(map(str,packet))}\n"
            zero_packet_hash.update(packet_line.encode("ascii"))
            witness = None
            for index, label in enumerate(packet):
                fixed = packet[:index] + packet[index + 1 :]
                mask = None
                for other in fixed:
                    if other not in safe_cache:
                        safe_cache[other] = vector_safe_mask(cells_array, modulus, other)
                    mask = safe_cache[other] if mask is None else mask & safe_cache[other]
                denominator, unit = denominator_and_unit(label, modulus)
                if not folded:
                    audit_fold(mask, cells, modulus, denominator)
                    folded = True
                orders = base.fast_torsion_orders(mask, modulus, denominator)
                if orders:
                    witness = (index, label, denominator, unit, mask.bit_count(), orders[0])
                    break
            if witness is None:
                local_zero_failures += 1
                zero_failure_count += 1
            else:
                zero_witness_hash.update((packet_line.strip() + "|" + repr(witness) + "\n").encode("ascii"))
        require(local_zero > 0 and local_zero_failures == 0 and folded, (first, body, local_zero, local_zero_failures))

        local_classes = local_small = local_full = local_units = 0
        for denominator, high in sorted(maxima.items()):
            high_value = high[0]
            threshold = required - high_value
            for labels, low_total in finite_triples(lows, threshold):
                local_classes += 1
                projective_class_count += 1
                mask = safe_cache[first]
                for label in labels:
                    if label not in safe_cache:
                        safe_cache[label] = vector_safe_mask(cells_array, modulus, label)
                    mask &= safe_cache[label]
                orders = base.fast_torsion_orders(mask, modulus, denominator)
                if orders:
                    local_small += 1
                    small_torsion_class_count += 1
                    record = (
                        f"z={first};E={','.join(map(str,body))};d={denominator};"
                        f"lows={','.join(map(str,labels))};orders={','.join(str(item[0]) for item in orders)}\n"
                    )
                    small_class_hash.update(record.encode("ascii"))
                    continue

                local_full += 1
                full_projection_class_count += 1
                denominator_failures[denominator] += 1
                residues_mask = vector_residue_presence(mask, modulus, denominator)
                require(residues_mask == (1 << denominator) - 1, (first, body, denominator, labels))
                remaining = required - low_total
                relevant = tuple(ray for ray in rays[denominator] if ray[0] >= remaining)
                require(relevant, (first, body, denominator, labels, remaining))
                class_full_margin = None
                class_pair_margin = None
                for _value, _high_label, unit, _residue, _raw in relevant:
                    local_units += 1
                    unit_test_count += 1
                    span, largest_gap = cyclic_span(range(denominator), denominator, unit)
                    full_margin = 7 * span - denominator
                    require(largest_gap == 1 and full_margin == 6 * denominator - 7, (denominator, unit, span, largest_gap))
                    t = (denominator + 6) // 7
                    second = (pow(unit, -1, denominator) * t) % denominator
                    pair_span, pair_gap = cyclic_span((0, second), denominator, unit)
                    pair_margin = 7 * t - denominator
                    require(
                        (unit * second) % denominator == t
                        and pair_span == min(t, denominator - t) == t
                        and pair_gap == denominator - t
                        and pair_margin >= 0,
                        (denominator, unit, second, t, pair_span, pair_gap),
                    )
                    class_full_margin = full_margin if class_full_margin is None else min(class_full_margin, full_margin)
                    class_pair_margin = pair_margin if class_pair_margin is None else min(class_pair_margin, pair_margin)
                    minimum_full_projection_margin = (
                        full_margin
                        if minimum_full_projection_margin is None
                        else min(minimum_full_projection_margin, full_margin)
                    )
                phase_records.append(
                    f"z={first};E={','.join(map(str,body))};d={denominator};"
                    f"lows={','.join(map(str,labels))};safe={mask.bit_count()};surjective=1;"
                    f"remaining={remaining};relevant_units={len(relevant)};"
                    f"full_arc_margin={class_full_margin};selected_pair_margin={class_pair_margin}"
                )

        exceptional_records.append(
            f"z={first};E={','.join(map(str,body))};L={modulus};floor={floor};lows={len(lows)};"
            f"denominators={len(maxima)};two_positive_gap={two_positive_gap};three_gap={three_gap};mixed_gap={mixed_gap};"
            f"zero_gap={zero_gap};two_gap={two_gap};best_high_le_second_low=1;"
            f"zero_packets={local_zero};zero_failures={local_zero_failures};"
            f"one_high_classes={local_classes};small_torsion={local_small};full_projection={local_full};"
            f"unit_tests={local_units};best_high_d={best_d};best_high_label={best_high[1]}"
        )

    require(zero_packet_count == EXPECTED_ZERO_PACKETS and zero_failure_count == 0, (zero_packet_count, zero_failure_count))
    require(projective_class_count == EXPECTED_PROJECTIVE_CLASSES, projective_class_count)
    require(small_torsion_class_count == EXPECTED_SMALL_TORSION_CLASSES, small_torsion_class_count)
    require(full_projection_class_count == EXPECTED_FULL_PROJECTION_CLASSES, full_projection_class_count)
    require(unit_test_count == EXPECTED_UNIT_TESTS and minimum_full_projection_margin == 59, (unit_test_count, minimum_full_projection_margin))
    require(denominator_failures == Counter({13: 30, 143: 11, 11: 9}), denominator_failures)
    require(
        nonpositive_triple_count == EXPECTED_NONPOSITIVE_TRIPLES
        and nonpositive_residue_ray_count == EXPECTED_NONPOSITIVE_RESIDUE_RAYS
        and nonpositive_unit_tests == EXPECTED_NONPOSITIVE_UNIT_TESTS
        and nonpositive_method_hist == EXPECTED_NONPOSITIVE_METHODS
        and minimum_nonpositive_fixed_cells == EXPECTED_NONPOSITIVE_MINIMUM_FIXED_CELLS
        and not nonpositive_failures,
        (
            nonpositive_triple_count,
            nonpositive_residue_ray_count,
            nonpositive_unit_tests,
            nonpositive_method_hist,
            minimum_nonpositive_fixed_cells,
            nonpositive_failures[:3],
        ),
    )
    require(TYPE_AUDIT["safe_full_cell_checks"] == 68 and TYPE_AUDIT["fold_full_cell_checks"] == 68, TYPE_AUDIT)

    # If z/L=a>=1, one fixed complete L-cell contains floor(a) full
    # punctured phase turns.  Their safe part has normalized mass at least
    # (6/7)floor(a)/a >=3/7, already above the k=2 cap 25/91.
    require(F(3, 7) > F(25, 91) and F(6, 7) > F(25, 91), "later-height safe-mass boundary")

    q7_first = base.circle_norm(F(16) * F(13, 224))
    q7_second = base.circle_norm(F(16) * F(29, 224))
    require(q7_first == q7_second == F(1, 14), (q7_first, q7_second))
    q8_first = base.circle_norm(F(63) * F(5, 336))
    q8_second = base.circle_norm(F(63) * F(11, 336))
    q8_overlap = F(1, 7) - F(1, 8)
    require(q8_first == q8_second == F(1, 16) and q8_overlap == F(1, 56), (q8_first, q8_second, q8_overlap))
    high_order_controls = []
    for denominator, expected_order in ((11, 11), (13, 13), (143, 143)):
        t = (denominator + 6) // 7
        order = denominator // gcd(denominator, t)
        require(order == expected_order and 7 * t >= denominator and 7 * (t - 1) < denominator, (denominator, t, order))
        high_order_controls.append((denominator, t, order, 7 * t - denominator, 6 * denominator - 7))
    for delta in range(1, 13):
        require(any((unit * delta) % 13 == 1 for unit in range(1, 13)), delta)

    phase_digest = semantic_digest(phase_records)
    lines = [
        "LRC14 projected k=2 1600..1679 finite, nonpositive-ray, and high-order phase closure",
        f"dependency={UTILITIES.relative_to(ROOT).as_posix()};lf_sha256={EXPECTED_UTILITIES_SHA256}",
        f"dependency={ATLAS_SOURCE.relative_to(ROOT).as_posix()};lf_sha256={EXPECTED_ATLAS_SOURCE_SHA256}",
        f"dependency={ATLAS_OUTPUT.relative_to(ROOT).as_posix()};lf_sha256={EXPECTED_ATLAS_OUTPUT_SHA256}",
        f"rows={len(rows)};exact_rows=40;high_tail_rows=28;finite_rows={len(finite_rows)};exceptional_rows={len(exceptional_rows)}",
        f"finite_packets={finite_packet_count};exact_packets={exact_packet_count};finite_high_packets={finite_high_packet_count};failures={finite_failure_count}",
        f"finite_packet_sha256={finite_packet_hash.hexdigest()}",
        f"finite_witness_sha256={finite_witness_hash.hexdigest()}",
        f"finite_failure_sha256={semantic_digest(finite_failure_lines)}",
        f"zero_high_packets={zero_packet_count};zero_high_failures={zero_failure_count}",
        f"zero_packet_sha256={zero_packet_hash.hexdigest()}",
        f"zero_witness_sha256={zero_witness_hash.hexdigest()}",
        f"projective_classes={projective_class_count};small_torsion_classes={small_torsion_class_count};full_projection_classes={full_projection_class_count}",
        f"full_projection_denominators={dict(sorted(denominator_failures.items()))};unit_tests={unit_test_count};"
        f"minimum_full_projection_arc_margin={minimum_full_projection_margin};minimum_selected_pair_margin=1",
        f"small_torsion_class_sha256={small_class_hash.hexdigest()}",
        f"phase_record_sha256={phase_digest};phase_unresolved=0",
        f"nonpositive_triples={nonpositive_triple_count};nonpositive_residue_rays={nonpositive_residue_ray_count};"
        f"nonpositive_unit_tests={nonpositive_unit_tests};minimum_fixed_cells={minimum_nonpositive_fixed_cells};"
        f"methods={dict(sorted(nonpositive_method_hist.items()))};failures={len(nonpositive_failures)}",
        f"nonpositive_failure_sha256={semantic_digest(tuple(map(repr, nonpositive_failures)))};"
        "later_height_safe_mass_lower=3/7;completion_cap=25/91",
        f"finite_height_histogram={dict(sorted(height_packets.items()))}",
        f"type_audit={dict(TYPE_AUDIT)};int64_limit={2**63-1};all_integer_bounds_pass=1;floats_used=0;"
        "block_or_equals_scalar_boolean=1;deterministic_order=rows,candidates,packets,punctures,orders,units",
        f"high_order_controls={tuple(high_order_controls)};fixed_pair_all_units_d13=0;unit_dependent_rule=1",
        "cyclic_carrier_lemma=minimal cyclic covering span >=1/7 closes;full projection and q<=7 torsion are sufficient special cases",
        "strict_open_controls=q7 equality disjoint;q8 overlap 1/56;phase ceiling equality disjoint",
        "exceptional_census_begin",
        *exceptional_records,
        "exceptional_census_end",
        "finite_row_census_begin",
        *finite_rows,
        "finite_row_census_end",
        "phase_records_begin",
        *phase_records,
        "phase_records_end",
        "nonpositive_records_begin",
        *nonpositive_records,
        "nonpositive_records_end",
        "scope=all 68 rows of the exact projected-k2 scalar atlas on 1600<=z1<=1679;positive finite/projective packets and every nonpositive ray height after the 27 fixed positive-low triples",
        "consequence=the projected scalar band 1600..1679 is empty;with the separately proved higher-band closure,z1<=1599",
        "all_exact_controls=PASS",
    ]
    return "\n".join(lines) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()
    payload = render()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(payload, encoding="utf-8", newline="\n")
    print(payload, end="")


if __name__ == "__main__":
    main()
