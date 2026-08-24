#!/usr/bin/env python3
"""No-import independent referee for THM-4013's extended Rule 30 census.

This file deliberately does not import, execute, or read either the primary
``probe_extended.py`` or its THM-3511 helper.  It reconstructs the finite
objects from the integer Rule 30 map and finite binary-tree sections.
"""

from __future__ import annotations

from collections import defaultdict, namedtuple
import hashlib
import json
import sys


sys.stdout.reconfigure(newline="\n")

LAST_SCALE = 15
TREE_HEADROOM = 4
UNIT_MODULUS_BITS = 4
HISTORY_LIMIT = 6
WORD_BITS = 64


def ensure(ok: bool, label: object) -> None:
    if not ok:
        raise RuntimeError(f"independent referee failure: {label!r}")


def step_rule30_integer(value: int, mask: int | None = None) -> int:
    result = value ^ ((value << 1) | (value << 2))
    return result if mask is None else result & mask


def two_order(value: int) -> int:
    ensure(value != 0, "2-adic order of zero")
    return (value & -value).bit_length() - 1


def finite_orbit(stop: int, width: int | None) -> list[int]:
    mask = None if width is None else (1 << width) - 1
    orbit = [1]
    value = 1
    for _ in range(stop):
        value = step_rule30_integer(value, mask)
        orbit.append(value)
    return orbit


def direct_marked_permutation(scale: int, valuation: int, depth: int) -> tuple[int, ...]:
    """Compute the marked normalized power map directly, without a word."""
    q = 1 << scale
    width = valuation + depth
    mask = (1 << width) - 1
    low_mask = (1 << valuation) - 1
    values = [1 | (tail << valuation) for tail in range(1 << depth)]
    for _ in range(q):
        values = [step_rule30_integer(value, mask) for value in values]
    ensure(all((value & low_mask) == 1 for value in values),
           ("marked prefix", scale, valuation, depth))
    answer = tuple(value >> valuation for value in values)
    ensure(sorted(answer) == list(range(1 << depth)), ("marked bijection", scale))
    return answer


def one_step_section(prefix: int, valuation: int, depth: int) -> tuple[int, ...]:
    """Section of one Rule 30 integer step over a fixed low-bit prefix."""
    width = valuation + depth
    mask = (1 << width) - 1
    low_mask = (1 << valuation) - 1
    outputs = [step_rule30_integer(prefix | (tail << valuation), mask)
               for tail in range(1 << depth)]
    expected_prefix = outputs[0] & low_mask
    ensure(all((value & low_mask) == expected_prefix for value in outputs),
           ("one-step common output prefix", prefix, valuation, depth))
    answer = tuple(value >> valuation for value in outputs)
    ensure(sorted(answer) == list(range(1 << depth)),
           ("one-step section bijection", prefix, valuation, depth))
    return answer


def inverse_permutation(p: tuple[int, ...]) -> tuple[int, ...]:
    result = [0] * len(p)
    for source, destination in enumerate(p):
        result[destination] = source
    return tuple(result)


def normalized_unit_mod(orbit: list[int], valuations: tuple[int, ...],
                        scale: int, phase: int, precision: int) -> int:
    q = 1 << scale
    modulus = 1 << WORD_BITS
    numerator = (orbit[phase + q] - orbit[phase]) % modulus
    valuation = valuations[scale]
    ensure(numerator % (1 << valuation) == 0,
           ("unit divisibility", scale, phase, precision))
    unit = (numerator >> valuation) & ((1 << precision) - 1)
    ensure(unit & 1 == 1, ("unit parity", scale, phase, precision))
    return unit


def arithmetic_history(orbit: list[int], valuations: tuple[int, ...],
                       gaps: tuple[int, ...], scale: int,
                       phase: int) -> tuple[object, ...]:
    result: list[object] = [None] * max(0, HISTORY_LIMIT - scale - 1)
    start = max(0, scale - HISTORY_LIMIT + 1)
    for level in range(start, scale + 1):
        d = gaps[level]
        precision = UNIT_MODULUS_BITS + d
        modulus = 1 << precision
        q = 1 << level
        left = normalized_unit_mod(orbit, valuations, level, phase, precision)
        right = normalized_unit_mod(orbit, valuations, level, phase + q, precision)
        ratio = (-right * pow(left, -1, modulus)) % modulus
        ensure((ratio - 1) % (1 << d) == 0,
               ("projective gap", scale, phase, level))
        coordinate = ((1 - ratio) >> d) & ((1 << UNIT_MODULUS_BITS) - 1)
        result.append((d, ratio & 15, coordinate))
    ensure(len(result) == HISTORY_LIMIT, ("history length", scale, phase))
    return tuple(result)


def carry_observation(exact: list[int], valuations: tuple[int, ...], n: int):
    alpha = two_order(n)
    base_valuation = valuations[alpha]
    center = (exact[n] >> n) & 1
    if base_valuation > n:
        ensure(center == 0, ("invisible center", n))
        return ("pre",), center

    level = n - base_valuation
    modulus = 1 << (level + 1)
    high_xor = 0
    lower_total = 0
    for scale in range(n.bit_length()):
        if n & (1 << scale) == 0:
            continue
        phase = n & ((1 << scale) - 1)
        valuation = valuations[scale]
        numerator = exact[phase + (1 << scale)] - exact[phase]
        ensure(numerator % (1 << valuation) == 0,
               ("exact shell divisibility", n, scale))
        unit = numerator >> valuation
        ensure(unit & 1 == 1, ("exact shell oddness", n, scale))
        residue = ((1 << (valuation - base_valuation)) * unit) % modulus
        high_xor ^= (residue >> level) & 1
        lower_total += residue & ((1 << level) - 1)
    carry = (lower_total >> level) & 1 if level else 0
    ensure(center == (high_xor ^ carry), ("center carry identity", n))
    return ("visible", high_xor, carry), center


Snapshot = namedtuple(
    "Snapshot",
    "n scale phase gap extension bank portrait units carry history chains center next_bank",
)


def target(record: Snapshot):
    return record.center, record.next_bank


def compact(record: Snapshot):
    return (
        record.n, record.scale, record.phase, record.gap, record.extension,
        record.bank, record.units, record.carry, record.history[-3:],
        record.chains, target(record),
    )


def construct_snapshots(exact: list[int], low: list[int],
                        valuations: tuple[int, ...], gaps: tuple[int, ...]):
    snapshots: list[Snapshot] = []
    stream = hashlib.sha256()
    tail_checks = conjugacy_entries = routed_checks = 0

    for scale in range(LAST_SCALE + 1):
        q = 1 << scale
        valuation = valuations[scale]
        gap = gaps[scale]
        depth = gap + TREE_HEADROOM
        size = 1 << depth
        marked = direct_marked_permutation(scale, valuation, depth)

        # T carries the marked chart to the current physical phase chart.
        transport = tuple(range(size))
        transport_inverse = tuple(range(size))
        low_prefix_mask = (1 << valuation) - 1

        for phase in range(q):
            owner = tuple(
                transport[marked[transport_inverse[value]]]
                for value in range(size)
            )
            ensure(sorted(owner) == list(range(size)), ("phase owner", scale, phase))
            conjugacy_entries += size

            tail = (low[phase] >> valuation) & (size - 1)
            first_tail = (low[phase + q] >> valuation) & (size - 1)
            second_tail = (low[phase + 2 * q] >> valuation) & (size - 1)
            ensure(owner[tail] == first_tail, ("first physical route", scale, phase))
            ensure(owner[owner[tail]] == second_tail,
                   ("second physical route", scale, phase))
            tail_checks += 2

            extension = tail & ((1 << gap) - 1)
            chains = []
            recovered_next = []
            for ray in range(3):
                initial = extension | (ray << gap)
                middle = owner[initial]
                final = owner[middle]
                ensure((final & ((1 << gap) - 1)) == extension,
                       ("routed base fibre", scale, phase, ray))
                chains.append((initial, middle, final))
                recovered_next.append((final >> gap) & 15)
                routed_checks += 1

            carry, center = carry_observation(exact, valuations, q + phase)
            units = tuple(
                normalized_unit_mod(low, valuations, scale, phase + j * q,
                                    UNIT_MODULUS_BITS)
                for j in range(4)
            )
            history = arithmetic_history(low, valuations, gaps, scale, phase)
            record = Snapshot(
                q + phase, scale, phase, gap, extension,
                tuple(owner[index] & 15 for index in range(3)),
                tuple(owner[index] & 15 for index in range(16)),
                units, carry, history, tuple(chains), center,
                tuple(recovered_next),
            )
            snapshots.append(record)
            encoded = json.dumps(record, separators=(",", ":")).encode("ascii")
            stream.update(len(encoded).to_bytes(4, "little"))
            stream.update(encoded)

            prefix = low[phase] & low_prefix_mask
            increment = one_step_section(prefix, valuation, depth)
            increment_inverse = inverse_permutation(increment)
            transport = tuple(increment[transport[value]] for value in range(size))
            transport_inverse = tuple(
                transport_inverse[increment_inverse[value]] for value in range(size)
            )
            ensure(all(transport_inverse[transport[value]] == value
                       for value in range(size)),
                   ("transport inverse", scale, phase))

    ensure(tuple(record.n for record in snapshots) == tuple(range(1, 1 << 16)),
           "complete chronological universe")
    return snapshots, stream.hexdigest(), (tail_checks, conjugacy_entries, routed_checks)


def base_key(record: Snapshot, history_depth: int = 3):
    ensure(0 <= history_depth <= HISTORY_LIMIT, ("history depth", history_depth))
    history = record.history[-history_depth:] if history_depth else ()
    return (record.gap, record.bank, record.portrait, record.units,
            record.carry, history)


def observation_key(record: Snapshot, label: str):
    if label.startswith("history"):
        return base_key(record, int(label[7:]))
    parts = list(base_key(record, 3))
    parts.append(record.extension)
    if label in {"c0", "c0_c1", "c0_c2", "c0_both",
                 "c0_c1_first", "c0_c2_first", "c0_both_first"}:
        parts.extend(record.chains[0][1:])
    if label in {"c0_c1_first", "c0_both_first"}:
        parts.append(record.chains[1][1])
    if label in {"c0_c2_first", "c0_both_first"}:
        parts.append(record.chains[2][1])
    if label in {"c0_c1", "c0_both"}:
        parts.append(record.chains[1])
    if label in {"c0_c2", "c0_both"}:
        parts.append(record.chains[2])
    return tuple(parts)


def mismatch_census(records: list[Snapshot], label: str):
    fibres: dict[tuple[object, ...], list[Snapshot]] = defaultdict(list)
    for record in records:
        fibres[observation_key(record, label)].append(record)

    mismatches = same_scale = cross_scale = 0
    first_pair = None
    first_cross_pair = None
    for key, fibre in fibres.items():
        if len({target(record) for record in fibre}) <= 1:
            continue
        mismatches += 1
        targets_by_scale: dict[int, set[object]] = defaultdict(set)
        for record in fibre:
            targets_by_scale[record.scale].add(target(record))
        if any(len(values) > 1 for values in targets_by_scale.values()):
            same_scale += 1
        scale_list = sorted(targets_by_scale)
        has_cross = any(
            left_target != right_target
            for i, left_scale in enumerate(scale_list)
            for right_scale in scale_list[i + 1:]
            for left_target in targets_by_scale[left_scale]
            for right_target in targets_by_scale[right_scale]
        )
        cross_scale += int(has_cross)

        for right_index, right in enumerate(fibre):
            for left in fibre[:right_index]:
                if target(left) == target(right):
                    continue
                candidate = (right.n, left.n, key, left, right)
                if first_pair is None or candidate[:2] < first_pair[:2]:
                    first_pair = candidate
                if left.scale != right.scale:
                    if first_cross_pair is None or candidate[:2] < first_cross_pair[:2]:
                        first_cross_pair = candidate

    nontrivial = sum(len(fibre) > 1 for fibre in fibres.values())
    counts = (len(fibres), nontrivial, mismatches, same_scale, cross_scale)
    return counts, first_pair, first_cross_pair


def adaptive_census(records: list[Snapshot], first_address_only: bool):
    base_fibres: dict[tuple[object, ...], list[Snapshot]] = defaultdict(list)
    for record in records:
        base_fibres[observation_key(record, "c0")].append(record)

    categories = ("base_enough", "either_one", "c1_only", "c2_only",
                  "both_needed", "pair_unresolved")
    fibre_counts = dict.fromkeys(categories, 0)
    record_counts = dict.fromkeys(categories, 0)
    examples: dict[str, tuple[int, ...]] = {}

    def view(record: Snapshot, ray: int):
        return record.chains[ray][1] if first_address_only else record.chains[ray]

    def resolves(fibre: list[Snapshot], rays: tuple[int, ...]) -> bool:
        quotient: dict[tuple[object, ...], set[object]] = defaultdict(set)
        for record in fibre:
            quotient[tuple(view(record, ray) for ray in rays)].add(target(record))
        return all(len(values) == 1 for values in quotient.values())

    for fibre in base_fibres.values():
        if len({target(record) for record in fibre}) == 1:
            category = "base_enough"
        else:
            c1 = resolves(fibre, (1,))
            c2 = resolves(fibre, (2,))
            pair = resolves(fibre, (1, 2))
            if not pair:
                category = "pair_unresolved"
            elif c1 and c2:
                category = "either_one"
            elif c1:
                category = "c1_only"
            elif c2:
                category = "c2_only"
            else:
                category = "both_needed"
            examples.setdefault(category, tuple(record.n for record in fibre[:3]))
        fibre_counts[category] += 1
        record_counts[category] += len(fibre)

    table = tuple((category, fibre_counts[category], record_counts[category])
                  for category in categories)
    return table, tuple((category, examples.get(category)) for category in categories)


def direct_next_control(records: list[Snapshot]) -> int:
    checks = 0
    for record in records:
        if record.scale == LAST_SCALE:
            continue
        next_n = (1 << (record.scale + 1)) + record.phase
        directly_built_bank = records[next_n - 1].bank
        ensure(record.next_bank == directly_built_bank,
               ("direct next-scale bank", record.n, next_n))
        checks += 1
    return checks


def same_scale_control() -> tuple[int, tuple[tuple[int, int, int, int], ...]]:
    value = 1
    checks = 0
    stop = 3 * (1 << LAST_SCALE)
    for time in range(1, stop + 1):
        value = step_rule30_integer(value)
        ensure(3 * (1 << (2 * time)) < 2 * value < 4 * (1 << (2 * time)),
               ("top-bit bound", time))
        checks += 1
    margins = []
    for scale in range(LAST_SCALE + 1):
        q = 1 << scale
        c = 1 << (2 * q)
        margin = 4 * c * c - 13 * c + 4
        ensure(margin > 0, ("positive exact defect margin", scale))
        margins.append((scale, q, margin.bit_length(), margin & 65535))
    return checks, tuple(margins)


def main() -> None:
    max_n = (1 << (LAST_SCALE + 1)) - 1
    low = finite_orbit(5 * (1 << LAST_SCALE), WORD_BITS)
    exact = finite_orbit(max_n, None)
    ensure(all((exact[index] & ((1 << WORD_BITS) - 1)) == low[index]
               for index in range(max_n + 1)), "exact/64-bit orbit agreement")

    valuations = tuple(two_order(low[1 << scale] - 1)
                       for scale in range(LAST_SCALE + 2))
    gaps = tuple(valuations[index + 1] - valuations[index]
                 for index in range(LAST_SCALE + 1))
    ensure(all(gap > 0 for gap in gaps), "positive valuation gaps")

    records, stream_hash, routing_controls = construct_snapshots(
        exact, low, valuations, gaps
    )
    transition_checks = direct_next_control(records)
    top_bit_checks, margins = same_scale_control()

    labels = ("history0", "history1", "history2", "history3", "history4",
              "history5", "history6", "c0", "c0_c1_first",
              "c0_c2_first", "c0_both_first", "c0_c1", "c0_c2", "c0_both")
    tables = []
    witnesses = {}
    cross_witnesses = {}
    for label in labels:
        counts, first, first_cross = mismatch_census(records, label)
        tables.append((label,) + counts)
        if first is not None:
            witnesses[label] = (first[2], compact(first[3]), compact(first[4]))
        if first_cross is not None:
            cross_witnesses[label] = (
                first_cross[2], compact(first_cross[3]), compact(first_cross[4])
            )

    adaptive_first = adaptive_census(records, True)
    adaptive_full = adaptive_census(records, False)

    n6 = records[5]
    ensure(n6.n == 6, "n=6 hostile lookup")
    # Re-evaluate n=6 from its stored physical owner chains is insufficient
    # for the marked-origin hostile, so reconstruct that phase owner once.
    scale, phase, gap = n6.scale, n6.phase, n6.gap
    valuation = valuations[scale]
    depth = gap + TREE_HEADROOM
    marked = direct_marked_permutation(scale, valuation, depth)
    size = 1 << depth
    transport = tuple(range(size))
    transport_inverse = tuple(range(size))
    for old_phase in range(phase):
        prefix = low[old_phase] & ((1 << valuation) - 1)
        inc = one_step_section(prefix, valuation, depth)
        inc_inv = inverse_permutation(inc)
        transport = tuple(inc[transport[x]] for x in range(size))
        transport_inverse = tuple(transport_inverse[inc_inv[x]] for x in range(size))
    owner6 = tuple(transport[marked[transport_inverse[x]]] for x in range(size))
    zero_ray_bank = tuple((owner6[owner6[ray << gap]] >> gap) & 15
                          for ray in range(3))
    ensure(zero_ray_bank == (11, 12, 9), "n=6 marked-origin hostile prediction")
    ensure(n6.next_bank == (15, 8, 5), "n=6 physical next bank")

    semantic = (
        max_n, valuations, gaps, routing_controls, transition_checks,
        top_bit_checks, margins, tuple(tables), adaptive_first, adaptive_full,
        stream_hash, witnesses["history3"], cross_witnesses["history3"],
        zero_ray_bank, n6.next_bank,
    )
    semantic_hash = hashlib.sha256(
        json.dumps(semantic, separators=(",", ":")).encode("ascii")
    ).hexdigest()

    print("RULE30_DEPTH15_HISTORY_ADAPTIVE_ROUTE_THM4013_INDEPENDENT_AUDIT")
    print("status=FINITE-EXACT_CANONICAL;scope=n=1..65535;no_all-scale_claim;no_prize_claim")
    print("implementation_firewall=no_import_no_exec_no_read_of_primary_or_thm3511_module")
    print(f"universe=1..{max_n};records={len(records)};scales=0..{LAST_SCALE}")
    print(f"valuations={valuations}")
    print(f"gaps={gaps}")
    print(f"record_stream_sha256={stream_hash}")
    print(f"physical_tail_conjugacy_route_controls={routing_controls}")
    print(f"direct_next_bank_controls={transition_checks}")
    print(f"same_scale_top_bit_checks={top_bit_checks}")
    print("same_scale_margins=" + repr(margins).replace(" ", ""))
    print("censuses=" + repr(tuple(tables)).replace(" ", ""))
    print("adaptive_first_address=" + repr(adaptive_first).replace(" ", ""))
    print("adaptive_full_chain=" + repr(adaptive_full).replace(" ", ""))
    print("history3_first_mismatch=" + repr(witnesses["history3"]).replace(" ", ""))
    print("history3_first_cross_scale=" + repr(cross_witnesses["history3"]).replace(" ", ""))
    print(f"n6_zero_ray_vs_physical={zero_ray_bank},{n6.next_bank}")
    print(f"semantic_sha256={semantic_hash}")
    print("ALL INDEPENDENT CHECKS PASSED")


if __name__ == "__main__":
    main()
