#!/usr/bin/env python3
"""Fail-closed exact closure of all fifteen projected-k2 rows at z1=1736.

Nine ordinary rows use the all-label exact status/projected engine.  Five of
the six HIGH-TAIL rows close by a strictly negative exact-ray gap.  The sixth
has a positive high-ray gap and therefore closes only after dropping the
forced-high constraint, which enlarges the universe, and applying the same
all-label PROJECT engine.  No implication in the opposite direction is used.
"""

from __future__ import annotations

import argparse
import concurrent.futures as cf
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations, combinations_with_replacement, product
from math import comb, gcd, lcm
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
FRONTIER_ENGINE = ROOT / "04-computation" / "lrc14_j7_k2_frontier_ray_status_closure_thm2941.py"
PROJECTED_ENGINE = ROOT / "04-computation" / "lrc14_j7_five_aligned_two_drift_projected_closure_thm2941.py"
UNIFORM_ENGINE = ROOT / "04-computation" / "lrc14_j7_k3_uniform_ray_status_closure_thm2941.py"
SCALAR_SOURCE = ROOT / "04-computation" / "lrc14_j7_k2_scalar_band_1680_1742_thm2941.py"
SCALAR_OUTPUT = ROOT / "05-knowledge" / "results" / "lrc14_j7_k2_scalar_band_1680_1742_thm2941.out"
OUTPUT_PATH = ROOT / "05-knowledge" / "results" / "lrc14_j7_k2_z1736_exact_ray_status_projected_closure_thm2941.out"

EXPECTED_FRONTIER_ENGINE_SHA256 = "bc4338935603b8971a99905033753458c880845213fddb5c4c19d8d53d6bc95b"
EXPECTED_PROJECTED_ENGINE_SHA256 = "76f891edfcc029a08202481304a809e03e8bd81f247afaeabab685825c4d3662"
EXPECTED_UNIFORM_ENGINE_SHA256 = "34ab29162ed33d90093e6d2bf781def36c420a1cd6596158b5d6579a3a8f3f46"
EXPECTED_SCALAR_SOURCE_SHA256 = "89016f939c961fa979ec5b30812981456df5bfb2af3066f1f1b38e5a83f1a412"
EXPECTED_SCALAR_OUTPUT_SHA256 = "4a36611b26585964e185bbaa3d583be3f1c67a7b608cca785920266bc217a779"
EXPECTED_PROFILE_SHA256 = "1c697481b42a05547b84e7fd8a45114faf6df2f923185de823e08ec5bf07643d"
EXPECTED_SEMANTIC_SHA256 = "a950a1a6d0d1fc887dde710778e87181b8de76db64688fd3d9b713fe03e72cd4"

EXACT_CASES = tuple((1736, body, "PROJECT") for body in (
    (1, 2, 8, 10, 12, 14),
    (1, 4, 8, 10, 12, 14),
    (1, 6, 8, 10, 12, 14),
    (2, 3, 8, 10, 12, 14),
    (2, 4, 6, 8, 10, 14),
    (2, 4, 8, 10, 12, 14),
    (2, 6, 8, 10, 12, 14),
    (3, 4, 8, 10, 12, 14),
    (4, 6, 8, 10, 12, 14),
))
HIGH_CASES = tuple((1736, body) for body in (
    (1, 2, 10, 11, 12, 14),
    (1, 2, 10, 12, 13, 14),
    (1, 4, 10, 12, 13, 14),
    (1, 8, 10, 12, 13, 14),
    (2, 3, 10, 11, 12, 14),
    (2, 8, 9, 10, 12, 14),
))
POSITIVE_HIGH_CASE = (1736, (1, 8, 10, 12, 13, 14))
NEXT_OCCUPIED = (
    (1732, (1, 4, 8, 10, 12, 14)),
    (1732, (2, 4, 8, 10, 12, 14)),
)
PRIOR_LEDGER = 2_239_789
CURRENT_LEDGER = 2_239_774

EXPECTED_EXACT = {
    (1736, (1, 2, 8, 10, 12, 14)): ((46, 0, 16, 30, 43), 4196, 43, 43, F(23897, 1434069), 4, "77491cb5de95588fee59c489e0f8427b9365799aa40d02350f8c68765e35e5e3"),
    (1736, (1, 4, 8, 10, 12, 14)): ((617, 0, 611, 6, 8), 4196, 8, 8, F(161107, 6585579), 3, "9cf36665aa9b9aa3b5d1845bb8abc3d8b23d1c917610ba951cb2c1cd651cc490"),
    (1736, (1, 6, 8, 10, 12, 14)): ((30, 0, 23, 7, 10), 4196, 10, 10, F(590, 2821), 2, "ae15f835044caddbfe3a81d637ce360b8d7d3a331c7e9906d2bbbbf86405a2d1"),
    (1736, (2, 3, 8, 10, 12, 14)): ((3, 0, 1, 2, 2), 4196, 2, 2, F(681, 2821), 1, "bb6d1f3147315ebda37641915eb8a385283bfad504e12c5e6fd4e27162958439"),
    (1736, (2, 4, 6, 8, 10, 14)): ((1, 0, 1, 0, 0), 4496, 0, 0, None, 0, "2e38e77b22c314a449e91fafed92a43826ac6aa403ae6a8acb6cf58239fbaf5d"),
    (1736, (2, 4, 8, 10, 12, 14)): ((198, 0, 109, 89, 194), 4496, 194, 194, F(681, 2821), 2, "d3484bf4398240d2a64b78dd5f4075037bfde530a599d8f9fea77f9123140748"),
    (1736, (2, 6, 8, 10, 12, 14)): ((30, 0, 9, 21, 29), 4356, 29, 29, F(49915, 510601), 2, "34c8a735ae7198acdef49ff96ceff9b151226e97cf0880c98abbba038f4b6e1c"),
    (1736, (3, 4, 8, 10, 12, 14)): ((3, 0, 3, 0, 0), 4476, 0, 0, None, 0, "2e38e77b22c314a449e91fafed92a43826ac6aa403ae6a8acb6cf58239fbaf5d"),
    (1736, (4, 6, 8, 10, 12, 14)): ((1, 0, 1, 0, 0), 4776, 0, 0, None, 0, "2e38e77b22c314a449e91fafed92a43826ac6aa403ae6a8acb6cf58239fbaf5d"),
}
EXPECTED_HIGH_GAPS = {
    (1736, (1, 2, 10, 11, 12, 14)): F(-17704278379, 151099294846500),
    (1736, (1, 2, 10, 12, 13, 14)): F(-2122333151, 17733467524500),
    (1736, (1, 4, 10, 12, 13, 14)): F(-7353136579, 118459563063660),
    POSITIVE_HIGH_CASE: F(91785, 20406470812),
    (1736, (2, 3, 10, 11, 12, 14)): F(-47509753883, 178751453380500),
    (1736, (2, 8, 9, 10, 12, 14)): F(-7891025249, 164964861300785),
}
EXPECTED_HIGH_PROJECT = (
    (827, 2, 607, 218, 749),
    48660,
    749,
    749,
    F(121, 18109),
    31,
    "aa644ef827e84923f548a30c13d415d02f8633b949787b115f9d7b5925774ae5",
)
SUFFIX_SLOTS = 4
HIGH_WALL_RATIO = F(13, 150)
SCALAR_ETA = F(1, 91)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_sha256(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


require(file_sha256(FRONTIER_ENGINE) == EXPECTED_FRONTIER_ENGINE_SHA256, "frontier engine changed")
require(file_sha256(PROJECTED_ENGINE) == EXPECTED_PROJECTED_ENGINE_SHA256, "projected engine changed")
require(file_sha256(UNIFORM_ENGINE) == EXPECTED_UNIFORM_ENGINE_SHA256, "uniform engine changed")
require(file_sha256(SCALAR_SOURCE) == EXPECTED_SCALAR_SOURCE_SHA256, "scalar source changed")
require(file_sha256(SCALAR_OUTPUT) == EXPECTED_SCALAR_OUTPUT_SHA256, "scalar output changed")
require(PRIOR_LEDGER - len(EXACT_CASES) - len(HIGH_CASES) == CURRENT_LEDGER, "ledger arithmetic changed")


def load(name: str, path: Path):
    spec = spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, "cannot load inherited closure engine")
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def projected_safe_lower_bound(projected, cells, ruler, labels, stop_at_cap=True):
    common_danger = ((F(0), F(1)),)
    cells_used = 0
    for cell in cells:
        local_union = projected.merge_fraction(
            [
                interval
                for label in labels
                for interval in projected.phase_danger(cell, label, ruler)
            ]
        )
        common_danger = projected.intersect_fraction(common_danger, local_union)
        cells_used += 1
        if stop_at_cap and projected.interval_mass(common_danger) <= 1 - F(25, 91):
            break
    return 1 - projected.interval_mass(common_danger), cells_used, common_danger


def direct_projected_mass(projected, body, ruler, labels):
    carrier = tuple(
        (F(left, projected.A.RULER), F(right, projected.A.RULER))
        for left, right in projected.A.carrier_for(body)
    )
    removed = projected.merge_fraction(
        [interval for label in labels for interval in projected.danger_fraction(label)]
    )
    residual = projected.subtract_fraction(carrier, removed)
    image = []
    for left, right in residual:
        scaled_left = ruler * left
        scaled_right = ruler * right
        for integer in range(
            projected.floor_fraction(scaled_left), projected.ceil_fraction(scaled_right)
        ):
            piece_left = max(scaled_left, F(integer)) - integer
            piece_right = min(scaled_right, F(integer + 1)) - integer
            if piece_left < piece_right:
                image.append((piece_left, piece_right))
    return projected.interval_mass(projected.merge_fraction(image))


def ray_and_status(base, first, body, carrier, h, lower, ruler):
    uniform = base.U
    first_delta = uniform.delta(carrier, h, first)
    first_d = ruler // gcd(ruler, first)
    amplitudes = [F(0)]
    signs = Counter()
    for residue in range(1, ruler):
        amplitude = residue * uniform.delta(carrier, h, residue)
        require(
            (residue + ruler) * uniform.delta(carrier, h, residue + ruler) == amplitude,
            (first, body, "ray recurrence", residue),
        )
        amplitudes.append(amplitude)
        signs[(amplitude > 0) - (amplitude < 0)] += 1
    require(
        all(
            amplitudes[ruler - residue] == -amplitudes[residue]
            for residue in range(1, ruler)
        ),
        (first, body, "ray antipode"),
    )
    require(
        ruler * uniform.delta(carrier, h, ruler) == 0,
        (first, body, "aligned ray is nonzero"),
    )
    require(
        signs[1] == signs[-1] and sum(signs.values()) == ruler - 1,
        (first, body, "ray sign census", signs),
    )
    ray_digest = sha256(repr(tuple(amplitudes)).encode()).hexdigest()
    divisors = tuple(d for d in uniform.support.divisors(ruler) if d > 1)

    @lru_cache(maxsize=None)
    def top_values(d, multiplicity):
        candidates = []
        for direction in range(1, d):
            if gcd(direction, d) != 1:
                continue
            residue = (ruler // d) * direction
            amplitude = amplitudes[residue]
            if amplitude < 0:
                continue
            label = residue
            if label <= first:
                label += ((first + 1 - label + ruler - 1) // ruler) * ruler
            candidates.extend(
                (amplitude / (label + offset * ruler), label + offset * ruler, residue)
                for offset in range(multiplicity)
            )
        candidates.sort(key=lambda row: (-row[0], row[1], row[2]))
        require(len(candidates) >= multiplicity, (first, body, "missing ray", d))
        return tuple(candidates[:multiplicity])

    trials = comb(len(divisors) + SUFFIX_SLOTS - 1, SUFFIX_SLOTS)
    scalar = []
    for tail in combinations_with_replacement(divisors, SUFFIX_SLOTS):
        upper = first_delta
        labels = []
        for d, multiplicity in Counter(tail).items():
            chosen = top_values(d, multiplicity)
            upper += sum((value for value, _label, _residue in chosen), F())
            labels.extend(label for _value, label, _residue in chosen)
        if upper >= lower:
            scalar.append(
                (tuple(sorted((first_d, *tail))), upper, (first, *sorted(labels)))
            )
    require(
        len({ds for ds, _upper, _labels in scalar}) == len(scalar),
        (first, body, "duplicate denominator state"),
    )

    actual_ruler, ranges = uniform.support.safe_cell_ranges(body)
    require(actual_ruler == ruler, (first, body, "safe-cell ruler"))
    arcs_cache = {}
    crude_kills = []
    crude_survivors = []
    for ds, upper, labels in scalar:
        common = lcm(*ds)
        if common not in arcs_cache:
            arcs_cache[common] = uniform.fibre.projected_support_arcs(common, ranges)
        arcs = arcs_cache[common]
        witness = None
        for q in uniform.support.divisors(common):
            histogram = uniform.fibre.residue_load_histogram(arcs, q)
            target = max(load for load, count in histogram if count)
            capacity = sum(uniform.fibre_cap(common, d, q) for d in ds)
            if target > capacity:
                witness = (q, common // q, target, capacity)
                break
        row = (ds, upper, labels, witness)
        (crude_survivors if witness is None else crude_kills).append(row)

    status_kills = []
    status_survivors = []
    for ds, upper, labels, _crude_witness in crude_survivors:
        common = lcm(*ds)
        arcs = arcs_cache[common]
        witness = None
        for modulus in uniform.support.divisors(common):
            q = common // modulus
            marginals, capacities = base.hunter_status_data5(common, ds, q)
            histogram = uniform.fibre.residue_load_histogram(arcs, q)
            feasible, certificate = base.common_status_feasible5(
                q, marginals, capacities, histogram
            )
            if not feasible:
                require(certificate is not None, (first, body, ds, "missing witness"))
                witness = (
                    q,
                    modulus,
                    marginals,
                    tuple(sorted(set(capacities))),
                    histogram,
                    certificate,
                )
                break
        row = (ds, upper, labels, witness)
        (status_survivors if witness is None else status_kills).append(row)

    stages = (scalar, crude_kills, status_kills, status_survivors)
    stage_digests = tuple(
        sha256(repr(tuple(rows)).encode()).hexdigest() for rows in stages
    )
    return (
        tuple(amplitudes),
        ray_digest,
        len(divisors),
        trials,
        first_delta,
        first_d,
        tuple(scalar),
        tuple(crude_kills),
        tuple(status_kills),
        tuple(status_survivors),
        stage_digests,
    )


def projected_packets(
    projected, first, body, carrier, h, lower, ruler, amplitudes, first_delta, states
):
    del carrier, h
    cells = projected.body_cells(projected.A.carrier_for(body), ruler)

    @lru_cache(maxsize=None)
    def ray_heads(d, multiplicity):
        candidates = []
        for direction in range(1, d):
            if gcd(direction, d) != 1:
                continue
            residue = (ruler // d) * direction
            amplitude = amplitudes[residue]
            if amplitude < 0:
                continue
            label = residue
            if label <= first:
                label += ((first + 1 - label + ruler - 1) // ruler) * ruler
            candidates.extend(
                (amplitude / (label + offset * ruler), label + offset * ruler, residue)
                for offset in range(multiplicity)
            )
        candidates.sort(key=lambda row: (-row[0], row[1], row[2]))
        return tuple(candidates[:multiplicity])

    @lru_cache(maxsize=None)
    def eligible_group(d, multiplicity, slack):
        heads = ray_heads(d, multiplicity)
        top_sum = sum((row[0] for row in heads), F())
        threshold = heads[-1][0] - slack
        require(
            threshold > 0,
            (first, body, "scalar slack did not make ray finite", d, multiplicity, slack),
        )
        candidates = []
        for direction in range(1, d):
            if gcd(direction, d) != 1:
                continue
            residue = (ruler // d) * direction
            amplitude = amplitudes[residue]
            if amplitude <= 0:
                continue
            label = residue
            if label <= first:
                label += ((first + 1 - label + ruler - 1) // ruler) * ruler
            while amplitude / label >= threshold:
                candidates.append((amplitude / label, label, residue))
                label += ruler
        candidates.sort(key=lambda row: (-row[0], row[1], row[2]))
        choices = []
        for chosen in combinations(candidates, multiplicity):
            labels = tuple(sorted(row[1] for row in chosen))
            if len(set(labels)) != multiplicity:
                continue
            value = sum((row[0] for row in chosen), F())
            deficit = top_sum - value
            if deficit <= slack:
                choices.append((deficit, labels, value))
        choices.sort()
        require(
            choices and choices[0][0] == 0,
            (first, body, "maximizing group choice missing", d, multiplicity),
        )
        return top_sum, threshold, tuple(candidates), tuple(choices)

    packet_count = 0
    killed_count = 0
    survivor_packets = []
    minimum_margin = None
    maximum_cells_used = 0
    minimum_row = None
    state_audit = []
    first_d = ruler // gcd(ruler, first)
    for state_index, (ds, upper, maximizing_labels, witness) in enumerate(states, 1):
        require(witness is None, (first, body, "status survivor has witness", ds))
        slack = upper - lower
        tail = list(ds)
        tail.remove(first_d)
        groups = []
        for d, multiplicity in sorted(Counter(tail).items()):
            top_sum, threshold, candidates, choices = eligible_group(d, multiplicity, slack)
            groups.append((d, multiplicity, top_sum, threshold, len(candidates), choices))
        expected_upper = first_delta + sum((group[2] for group in groups), F())
        require(
            expected_upper == upper,
            (first, body, "state upper changed", state_index, expected_upper, upper),
        )
        state_packets = 0
        state_kills = 0
        state_minimum = None
        for selection in product(*(group[5] for group in groups)):
            total_deficit = sum((choice[0] for choice in selection), F())
            if total_deficit > slack:
                continue
            labels = (first, *sorted(label for choice in selection for label in choice[1]))
            require(
                len(labels) == 5 and len(set(labels)) == 5,
                (first, body, "literal labels not distinct", labels),
            )
            actual_upper = first_delta + sum((choice[2] for choice in selection), F())
            require(
                actual_upper == upper - total_deficit and actual_upper >= lower,
                (first, body, "literal scalar arithmetic", labels),
            )
            projected_lower, cells_used, _common = projected_safe_lower_bound(
                projected, cells, ruler, labels
            )
            margin = projected_lower - F(25, 91)
            packet_count += 1
            state_packets += 1
            maximum_cells_used = max(maximum_cells_used, cells_used)
            minimum_margin = margin if minimum_margin is None else min(minimum_margin, margin)
            state_minimum = margin if state_minimum is None else min(state_minimum, margin)
            if minimum_row is None or (margin, labels) < (minimum_row[0], minimum_row[1]):
                minimum_row = (margin, labels, projected_lower, cells_used)
            if margin > 0:
                killed_count += 1
                state_kills += 1
            else:
                survivor_packets.append((state_index, ds, labels, actual_upper, projected_lower))
        require(state_packets > 0, (first, body, "status state has no literal packet", ds))
        state_audit.append(
            (
                state_index,
                ds,
                upper,
                maximizing_labels,
                slack,
                tuple(
                    (d, multiplicity, top_sum, threshold, candidate_count, len(choices))
                    for d, multiplicity, top_sum, threshold, candidate_count, choices in groups
                ),
                state_packets,
                state_kills,
                state_minimum,
            )
        )
    require(
        packet_count == killed_count and not survivor_packets,
        (first, body, "projected packet survived", survivor_packets),
    )
    direct_mass = None
    if minimum_row is not None:
        _margin, control_labels, prefix_mass, _prefix_cells = minimum_row
        full_cell_mass, full_cells, _common = projected_safe_lower_bound(
            projected, cells, ruler, control_labels, stop_at_cap=False
        )
        direct_mass = direct_projected_mass(projected, body, ruler, control_labels)
        require(
            full_cells == len(cells)
            and direct_mass == full_cell_mass
            and direct_mass >= prefix_mass,
            (first, body, "independent projected control"),
        )
    state_digest = sha256(repr(tuple(state_audit)).encode()).hexdigest()
    return (
        len(cells),
        packet_count,
        killed_count,
        minimum_margin,
        maximum_cells_used,
        minimum_row,
        direct_mass,
        state_digest,
        tuple(state_audit),
    )


def exact_worker(item):
    try:
        base = load("k2_z1736_frontier", FRONTIER_ENGINE)
        projected = load("k2_z1736_projected", PROJECTED_ENGINE)
        first, body, mode = item
        carrier = base.U.suffix.A.carrier_for(body)
        require(projected.A.carrier_for(body) == carrier, "carrier engines disagree")
        h = F(sum(right - left for left, right in carrier), base.U.suffix.A.RULER)
        lower = h * base.U.suffix.ETAS[2]
        ruler = 14 * lcm(*body)
        require(projected.A.RULER % ruler == 0, "unexpected ruler")
        (
            amplitudes,
            ray_digest,
            divisor_count,
            trials,
            first_delta,
            first_d,
            scalar,
            crude_kills,
            status_kills,
            states,
            stage_digests,
        ) = ray_and_status(base, first, body, carrier, h, lower, ruler)
        projected_result = projected_packets(
            projected,
            first,
            body,
            carrier,
            h,
            lower,
            ruler,
            amplitudes,
            first_delta,
            states,
        )
        counts = (
            len(scalar),
            len(crude_kills),
            len(status_kills),
            len(states),
            projected_result[1],
        )
        row = (
            first,
            body,
            mode,
            h,
            len(carrier),
            ruler,
            lower,
            first_delta,
            first_d,
            ray_digest,
            divisor_count,
            trials,
            counts,
            stage_digests,
            *projected_result[:-1],
        )
        return ("OK", row)
    except RuntimeError as exc:
        return ("SURVIVOR", item, repr(exc))


def high_worker(item):
    module = load("k2_z1736_uniform", UNIFORM_ENGINE)
    first, body = item
    carrier = module.suffix.A.carrier_for(body)
    h = F(sum(right - left for left, right in carrier), module.suffix.A.RULER)
    lower = h * SCALAR_ETA
    ruler = 14 * lcm(*body)
    require(module.suffix.A.RULER % ruler == 0, (item, "body ruler left master ruler"))
    wall = HIGH_WALL_RATIO * ruler
    high_floor = max(15, wall.numerator // wall.denominator + 1)
    require(first < high_floor, (item, "case does not force a later high label"))

    amplitudes = [F(0)]
    signs = Counter()
    for residue in range(1, ruler):
        amplitude = residue * module.delta(carrier, h, residue)
        require(
            (residue + ruler) * module.delta(carrier, h, residue + ruler) == amplitude,
            (item, "ray recurrence", residue),
        )
        amplitudes.append(amplitude)
        signs[(amplitude > 0) - (amplitude < 0)] += 1
    require(
        all(amplitudes[ruler - residue] == -amplitudes[residue] for residue in range(1, ruler)),
        (item, "ray antipode"),
    )
    require(ruler * module.delta(carrier, h, ruler) == 0, (item, "aligned ray is nonzero"))
    require(
        signs[1] == signs[-1] and sum(signs.values()) == ruler - 1,
        (item, "ray sign census", signs),
    )

    arbitrary = []
    high = []
    omitted = []
    for residue in range(1, ruler):
        amplitude = amplitudes[residue]
        if amplitude <= 0:
            continue
        first_label = residue
        if first_label <= first:
            first_label += ((first + 1 - first_label + ruler - 1) // ruler) * ruler
        for offset in range(SUFFIX_SLOTS):
            label = first_label + offset * ruler
            arbitrary.append((amplitude / label, label, residue, offset))
        fifth_label = first_label + SUFFIX_SLOTS * ruler
        omitted.append((amplitude / fifth_label, fifth_label, residue))
        high_label = residue
        if high_label < high_floor:
            high_label += ((high_floor - high_label + ruler - 1) // ruler) * ruler
        high.append((amplitude / high_label, high_label, residue))

    rank4 = tuple(sorted(arbitrary, key=lambda row: (-row[0], row[1:]))[:4])
    require(len(rank4) == 4 and rank4[-1][0] > 0, (item, "positive top four missing"))
    omitted_max = min(omitted, key=lambda row: (-row[0], row[1:]))
    require(omitted_max[0] <= rank4[-1][0], (item, "first-four truncation failed"))
    best_high = min(high, key=lambda row: (-row[0], row[1:]))
    if any(label >= high_floor for _value, label, _residue, _offset in rank4):
        constrained = rank4
        branch = "UNRESTRICTED_TOP4_ALREADY_HIGH"
    else:
        constrained = (*rank4[:3], (best_high[0], best_high[1], best_high[2], 0))
        branch = "GLOBAL_TOP3_PLUS_BEST_HIGH"
    require(
        len({row[1] for row in constrained}) == SUFFIX_SLOTS
        and any(row[1] >= high_floor for row in constrained),
        (item, "constrained suffix malformed", constrained),
    )
    first_delta = module.delta(carrier, h, first)
    upper = first_delta + sum((row[0] for row in constrained), F())
    gap = upper - lower
    for value, label, residue, _offset in constrained:
        require(
            module.delta(carrier, h, label) == value
            and label > first
            and label % ruler == residue % ruler,
            (item, "chosen singleton control", label),
        )
    return (
        first,
        body,
        h,
        len(carrier),
        ruler,
        high_floor,
        first_delta,
        lower,
        signs[1],
        signs[-1],
        signs[0],
        sha256(repr(tuple(amplitudes)).encode()).hexdigest(),
        rank4,
        omitted_max,
        best_high,
        branch,
        tuple(constrained),
        upper,
        gap,
    )


def exact_fields(row):
    return (row[12], row[14], row[15], row[16], row[17], row[18], row[21])


def exact_line(row, prefix: str = "EXACT") -> str:
    first, body = row[:2]
    counts = row[12]
    cells, packets, kills, margin, max_prefix, minimum, direct, state_hash = row[14:22]
    return (
        f"{prefix};z1={first};E={','.join(map(str, body))};counts={counts};cells={cells};"
        f"packets={packets};kills={kills};min_margin={margin};max_prefix={max_prefix};"
        f"minimum={minimum};direct={direct};state_sha256={state_hash};closed=1"
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--output", type=Path, default=OUTPUT_PATH)
    args = parser.parse_args()
    require(args.workers >= 1, "workers must be positive")
    with cf.ProcessPoolExecutor(max_workers=args.workers) as pool:
        exact_raw = tuple(pool.map(exact_worker, EXACT_CASES))
        high = tuple(pool.map(high_worker, HIGH_CASES))
        high_project_items = tuple((row[0], row[1], "PROJECT") for row in high if row[-1] >= 0)
        high_project_raw = tuple(pool.map(exact_worker, high_project_items))

    failures = tuple(row for row in (*exact_raw, *high_project_raw) if row[0] != "OK")
    require(not failures, f"a widened all-label packet survived: {failures}")
    exact = tuple(row[1] for row in exact_raw)
    high_project = tuple(row[1] for row in high_project_raw)
    require(tuple((row[0], row[1]) for row in exact) == tuple(item[:2] for item in EXACT_CASES), "ordinary case order changed")
    require(tuple((row[0], row[1]) for row in high) == HIGH_CASES, "high case order changed")
    require(high_project_items == (POSITIVE_HIGH_CASE + ("PROJECT",),), "HIGH-to-PROJECT routing changed")
    require(len(high_project) == 1, "widened project census changed")

    for row in exact:
        key = (row[0], row[1])
        require(exact_fields(row) == EXPECTED_EXACT[key], f"ordinary closure changed: {key}")
        counts, _cells, packets, kills, margin, _prefix, _state = exact_fields(row)
        require(counts[3] == 0 or (packets == kills and margin is not None and margin > 0), f"ordinary packet survived: {key}")
    for row in high:
        key = (row[0], row[1])
        require(row[-1] == EXPECTED_HIGH_GAPS[key], f"high gap changed: {key}")
        require((row[-1] >= 0) == (key == POSITIVE_HIGH_CASE), f"high routing sign changed: {key}")
    require(exact_fields(high_project[0]) == EXPECTED_HIGH_PROJECT, "widened PROJECT closure changed")

    ordinary_totals = tuple(sum(row[12][i] for row in exact) for i in range(5))
    require(ordinary_totals == (929, 0, 774, 155, 286), "ordinary aggregate changed")
    require(sum(row[16] for row in exact) == 286, "ordinary projected kills changed")
    profile_payload = (EXACT_CASES, HIGH_CASES, exact_raw, high, high_project_raw)
    profile_hash = sha256(repr(profile_payload).encode()).hexdigest()
    require(profile_hash == EXPECTED_PROFILE_SHA256, "closure profile changed")
    semantic_payload = (
        EXPECTED_FRONTIER_ENGINE_SHA256,
        EXPECTED_PROJECTED_ENGINE_SHA256,
        EXPECTED_UNIFORM_ENGINE_SHA256,
        EXPECTED_SCALAR_SOURCE_SHA256,
        EXPECTED_SCALAR_OUTPUT_SHA256,
        EXPECTED_PROFILE_SHA256,
        ordinary_totals,
        EXPECTED_HIGH_GAPS,
        EXPECTED_HIGH_PROJECT,
        NEXT_OCCUPIED,
        PRIOR_LEDGER,
        CURRENT_LEDGER,
    )
    semantic_hash = sha256(repr(semantic_payload).encode()).hexdigest()
    require(semantic_hash == EXPECTED_SEMANTIC_SHA256, "semantic digest changed")

    lines = [
        "LRC14 projected k=2 exact ray/status/projected closure at z1=1736",
        f"frontier_engine_sha256={file_sha256(FRONTIER_ENGINE)}",
        f"projected_engine_sha256={file_sha256(PROJECTED_ENGINE)}",
        f"uniform_engine_sha256={file_sha256(UNIFORM_ENGINE)}",
        f"scalar_source_sha256={file_sha256(SCALAR_SOURCE)}",
        f"scalar_output_sha256={file_sha256(SCALAR_OUTPUT)}",
        "scope=all fifteen z1736 scalar rows;all later distinct nonaligned labels;no finite label horizon",
        "routing=nine ordinary PROJECT;five HIGH strict-negative exact rays;one positive-HIGH row widened to all-label PROJECT",
        "direction=dropping forced-HIGH enlarges the universe;the widened closure is sufficient only in that direction",
        f"ordinary_totals=scalar:{ordinary_totals[0]};crude:{ordinary_totals[1]};status:{ordinary_totals[2]};residual:{ordinary_totals[3]};packets:{ordinary_totals[4]};kills:286",
    ]
    lines.extend(exact_line(row) for row in exact)
    for row in high:
        key = (row[0], row[1])
        lines.append(
            f"HIGH;z1={key[0]};E={','.join(map(str, key[1]))};L={row[4]};high_floor={row[5]};"
            f"branch={row[15]};upper={row[-2]};gap={row[-1]};"
            f"route={'PROJECT' if key == POSITIVE_HIGH_CASE else 'EXACT-RAY'};"
            f"ray_closed={int(row[-1] < 0)}"
        )
    lines.append(exact_line(high_project[0], "HIGH-PROJECT"))
    lines.extend((
        "closure=15/15;survivors=0;global_min_projected_margin=121/18109;global_max_prefix=31",
        f"next_occupied={NEXT_OCCUPIED}",
        "consequence=projected k=2 first drift label z1<=1732",
        f"ledger={PRIOR_LEDGER}-15={CURRENT_LEDGER}",
        f"profile_sha256={profile_hash}",
        f"semantic_sha256={semantic_hash}",
        "all_exact_controls=PASS",
    ))
    text = "\n".join(lines) + "\n"
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(text, encoding="utf-8", newline="\n")
    print(text, end="")


if __name__ == "__main__":
    main()
