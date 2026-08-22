#!/usr/bin/env python3
"""Independent audit of the two costly z1=216 gcd-eight projected rows.

This file never imports the lrc_probe scratch wrapper.  It parses the pinned
THM-2941 atlas directly, invokes the promoted THM-3078 screen engine with a
separately spelled exact Farkas checker, and rebuilds every terminal complete
cell on the full Z/LZ grid rather than through the inherited safe-range word.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import math
import multiprocessing as mp
import re
from collections import Counter, defaultdict
from fractions import Fraction as F
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
ATLAS = ROOT / "05-knowledge/results/lrc14_j7_k3_projected_scalar_body_atlas_thm2941.out"
SOURCE_3078 = ROOT / "04-computation/lrc14_j7_k3_z234_direct_farkas_four_two_high_boundary_thm3078.py"
SOURCE_3139 = ROOT / "04-computation/lrc14_j7_k3_z225_terminal_z224_screen_double_layer_descent_thm3139.py"
LEVEL = 216
TARGETS = (238, 370)
BASE_RULER = 5_045_040
EXPECTED_ATLAS_SHA = "cee82237ce1f51729813b9c916edd3353204c18172abe1d71278dee2c5562eda"
EXPECTED_SOURCE_3078_SHA = "2a051babe109f56056fe61476870f8e2e13cfc99b2f9bb7ac122b8780c8fa168"
EXPECTED_SOURCE_3139_SHA = "92e1c22088998d37db89020ac1ebbb7d1ee17e3539152ab205f8b0cd92532e36"
EXPECTED_Z216_ROW_SHA = "53db9e1d3df2cf2b0398847682d909da81705e43a53ae2553d102fd152337649"


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def lf_sha(path: Path) -> str:
    payload = path.read_bytes().replace(b"\r\n", b"\n")
    require(b"\r" not in payload, (path, "bare CR"))
    return hashlib.sha256(payload).hexdigest()


def load(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def atlas_rows():
    pattern = re.compile(
        r"^row=E=([0-9,]+);h=([^;]+);r=([0-9]+);"
        r"L=([0-9]+);high=([0-9]+);z1=([0-9]+);"
    )
    all_rows = []
    selected = []
    for line_number, line in enumerate(ATLAS.read_text(encoding="utf-8").splitlines(), 1):
        if not line.startswith("row="):
            continue
        match = pattern.match(line)
        require(match is not None, (line_number, line))
        body = tuple(map(int, match.group(1).split(",")))
        components = int(match.group(3))
        ruler = int(match.group(4))
        high_floor = int(match.group(5))
        first = int(match.group(6))
        require(ruler == 14 * math.lcm(*body), (body, "ruler"))
        require(high_floor == (13 * ruler) // 132 + 1, (body, "high floor"))
        record = (body, ruler, high_floor, first < high_floor, components, line_number)
        all_rows.append(record)
        if first == LEVEL:
            selected.append(record)
    require(len(all_rows) == 6060, len(all_rows))
    require(len(selected) == 480, len(selected))
    canonical = tuple(row[:4] for row in selected)
    require(hashlib.sha256(repr(canonical).encode()).hexdigest() == EXPECTED_Z216_ROW_SHA,
            "z216 row ordering")
    return tuple(selected)


class ExactDualAudit:
    def __init__(self):
        self.count = 0
        self.digest = hashlib.sha256()
        self.first = None

    def __call__(self, q, marginals, capacities, histogram, certificate):
        thresholds, alpha, z = certificate
        tail_rows = []
        tail_rhs = []
        rebuilt_thresholds = []
        for threshold, _count in histogram:
            if threshold <= 0:
                continue
            demand = sum(count for load, count in histogram if load >= threshold)
            good = tuple(int(capacity >= threshold) for capacity in capacities)
            if all(good):
                continue
            rebuilt_thresholds.append(threshold)
            tail_rows.append(good)
            tail_rhs.append(demand)
        require(tuple(rebuilt_thresholds) == tuple(thresholds), "threshold mismatch")
        require(len(alpha) == len(tail_rows) and len(z) == 5, "Farkas shape")
        require(all(value >= 0 for value in alpha), "negative Farkas alpha")
        equality_rows = [
            (1,) * 16,
            *[
                tuple((pattern >> index) & 1 for pattern in range(16))
                for index in range(4)
            ],
        ]
        equality_rhs = (q, *marginals)
        slacks = tuple(
            sum(z[row] * equality_rows[row][pattern] for row in range(5))
            - sum(alpha[row] * tail_rows[row][pattern] for row in range(len(alpha)))
            for pattern in range(16)
        )
        contradiction = sum(z[row] * equality_rhs[row] for row in range(5))
        contradiction -= sum(alpha[row] * tail_rhs[row] for row in range(len(alpha)))
        require(all(value >= 0 for value in slacks), "negative Farkas column")
        require(contradiction < 0, "nonnegative Farkas contradiction")
        instance = (q, marginals, tuple(sorted(set(capacities))), histogram)
        self.digest.update(repr(instance).encode() + b"\n")
        self.count += 1
        if self.first is None:
            self.first = (q, marginals, capacities, histogram, certificate)
        return contradiction

    def hostile_control(self):
        require(self.first is not None, "missing Farkas hostile control")
        q, marginals, capacities, histogram, certificate = self.first
        thresholds, alpha, z = certificate
        bad_alpha = (-F(1), *tuple(alpha[1:]))
        try:
            self(q, marginals, capacities, histogram, (thresholds, bad_alpha, z))
        except RuntimeError as error:
            require("negative Farkas alpha" in str(error), error)
        else:
            raise RuntimeError("negative-alpha mutation accepted")


def screen_worker(item):
    index, atlas_row = item
    body, ruler, high_floor, wall, components, _line_number = atlas_row
    module = load(f"independent_g8_screen_{index}", SOURCE_3078)
    audit = ExactDualAudit()
    module.eng.independent_farkas_check = audit
    task = (LEVEL, body, ruler, high_floor, wall)
    row = tuple(module.screen_engine.evaluate(task))
    require(row[:6] == (LEVEL, body, ruler, high_floor, ruler // 8, wall),
            (index, row[:6]))
    require(row[16] == row[11] == audit.count, (index, row[11], row[16], audit.count))
    audit.hostile_control()
    return index, row, audit.count, audit.digest.hexdigest(), components


def merge_intervals(rows):
    merged = []
    for left, right in sorted((left, right) for left, right in rows if left < right):
        if merged and left <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], right))
        else:
            merged.append((left, right))
    return tuple(merged)


def base_carrier(body):
    danger = []
    for speed in body:
        step = BASE_RULER // speed
        radius = BASE_RULER // (14 * speed)
        danger.append((0, radius))
        danger.extend((k * step - radius, k * step + radius) for k in range(1, speed))
        danger.append((BASE_RULER - radius, BASE_RULER))
    merged = merge_intervals(danger)
    safe = []
    cursor = 0
    for left, right in merged:
        if cursor < left:
            safe.append((cursor, left))
        cursor = max(cursor, right)
    if cursor < BASE_RULER:
        safe.append((cursor, BASE_RULER))
    return tuple(safe)


def danger_primitive(numerator):
    cycles, remainder = divmod(numerator, BASE_RULER)
    return (
        cycles * (BASE_RULER // 7)
        + min(remainder, BASE_RULER // 14)
        + max(0, remainder - 13 * BASE_RULER // 14)
    )


def singleton_coverage(carrier, label):
    numerator = sum(
        danger_primitive(label * right) - danger_primitive(label * left)
        for left, right in carrier
    )
    return F(numerator, BASE_RULER * label)


def independent_literal_tables(body, first, ruler, high_floor, needed):
    carrier = base_carrier(body)
    h = F(sum(right - left for left, right in carrier), BASE_RULER)

    def delta(label):
        return singleton_coverage(carrier, label) - h / 7

    low = {d: [] for d in needed}
    for label in range(first + 1, high_floor):
        d = ruler // math.gcd(ruler, label)
        if d in low:
            low[d].append((delta(label), label))
    signs = Counter()
    for d in low:
        signs.update((value > 0) - (value < 0) for value, _label in low[d])
        low[d].sort(key=lambda item: (-item[0], item[1]))
        low[d] = tuple(low[d])

    high = {}
    checks = 0
    for d in needed:
        step = ruler // d
        best = (F(0), None, None, None)
        for unit in range(1, d):
            if math.gcd(unit, d) != 1:
                continue
            residue = step * unit
            amplitude = residue * delta(residue)
            require((residue + ruler) * delta(residue + ruler) == amplitude,
                    (d, unit, "ray recurrence"))
            checks += 1
            if residue >= high_floor:
                label = residue
            else:
                label = residue + ((high_floor - residue + ruler - 1) // ruler) * ruler
            value = amplitude / label if amplitude > 0 else F(0)
            candidate = (value, label, unit, amplitude)
            if candidate[0] > best[0] or (
                candidate[0] == best[0]
                and candidate[1] is not None
                and (best[1] is None or candidate[1] < best[1])
            ):
                best = candidate
        require(best[1] is not None, (d, "empty high class"))
        high[d] = best
    return low, high, tuple(sorted(signs.items())), checks


def suffix_slots(mask, first_d):
    slots = list(mask)
    slots.remove(first_d)
    require(len(slots) == 3, (mask, first_d))
    return tuple(slots)


def independent_gap(required, residual, first_d, low, high):
    best = None
    for mask in residual:
        slots = suffix_slots(mask, first_d)
        for choice in range(8):
            if choice.bit_count() < 2:
                continue
            value = F(0)
            possible = True
            for pos, d in enumerate(slots):
                if (choice >> pos) & 1:
                    value += high[d][0]
                elif low[d]:
                    value += low[d][0][0]
                else:
                    possible = False
                    break
            if possible and (best is None or value > best):
                best = value
    require(best is not None, "no two-high hostile")
    return required - best


def finite_low_pairs(first, second, same, threshold):
    rows = []
    if same:
        for index, left in enumerate(first):
            if index + 1 >= len(first) or left[0] + first[index + 1][0] < threshold:
                break
            for right in first[index + 1:]:
                if left[0] + right[0] < threshold:
                    break
                if left[1] != right[1]:
                    rows.append((left, right))
    else:
        for left in first:
            if not second or left[0] + second[0][0] < threshold:
                break
            for right in second:
                if left[0] + right[0] < threshold:
                    break
                if left[1] != right[1]:
                    rows.append((left, right))
    return tuple(rows)


def independent_zero_and_cases(required, residual, first_d, low, high):
    zero = []
    cases = set()
    for mask in residual:
        slots = suffix_slots(mask, first_d)
        if all(low[d] for d in slots):
            value = sum((low[d][0][0] for d in slots), F(0))
            if value >= required:
                zero.append((mask, tuple(low[d][0][1] for d in slots), value - required))
        for high_index, high_d in enumerate(slots):
            low_ds = tuple(d for pos, d in enumerate(slots) if pos != high_index)
            if not low[low_ds[0]] or not low[low_ds[1]]:
                continue
            threshold = required - high[high_d][0]
            for left, right in finite_low_pairs(
                low[low_ds[0]], low[low_ds[1]], low_ds[0] == low_ds[1], threshold
            ):
                low_rows = tuple(sorted(((low_ds[0], left[1]), (low_ds[1], right[1]))))
                excess = left[0] + right[0] + high[high_d][0] - required
                cases.add((mask, high_d, low_rows, excess))
    return tuple(zero), tuple(sorted(cases))


def direct_complete_cells(body, first, low_labels, ruler):
    cells = np.arange(ruler, dtype=np.int64)
    clean = np.ones(ruler, dtype=np.bool_)
    for label in (*body, first, *low_labels):
        residues = (cells * int(label)) % int(ruler)
        clean &= 14 * residues >= ruler
        clean &= 14 * (residues + int(label)) <= 13 * ruler
    return cells[clean]


def terminal_audit(index, screen):
    module = load(f"independent_g8_terminal_{index}", SOURCE_3078)
    body, ruler, high_floor, first_d = screen[1], screen[2], screen[3], screen[4]
    residual = screen[13]
    module.eng.FIRST = LEVEL
    module.eng.ray.FIRST = LEVEL
    stream = module.eng.ray.Stream(body)
    require((stream.L, stream.high_floor, stream.first_d) ==
            (ruler, high_floor, first_d), (index, "stream metadata"))
    needed = {d for mask in residual for d in suffix_slots(mask, first_d)}
    inherited_low, inherited_high, inherited_signs, inherited_checks = (
        module.eng.build_literal_tables(stream, needed)
    )
    low, high, signs, checks = independent_literal_tables(
        body, LEVEL, ruler, high_floor, needed
    )
    require(low == inherited_low, (index, "low literal table"))
    require(high == inherited_high, (index, "high literal table"))
    require((signs, checks) == (inherited_signs, inherited_checks),
            (index, "ray sidecars"))
    required = stream.lower - stream.first_delta
    gap = independent_gap(required, residual, first_d, low, high)
    inherited_gap, _witness = module.eng.duplicate_two_high_gap(
        stream, residual, inherited_low, inherited_high
    )
    require(gap == inherited_gap and gap > 0, (index, gap, inherited_gap))
    zero, cases = independent_zero_and_cases(required, residual, first_d, low, high)
    inherited_zero = module.eng.zero_high_scalar_passes(stream, residual, inherited_low)
    inherited_cases = module.eng.one_high_cases(stream, residual, inherited_low, inherited_high)
    require(zero == inherited_zero, (index, "zero-high bank"))
    require(cases == inherited_cases, (index, "one-high bank"))

    by_labels = defaultdict(set)
    for _mask, high_d, low_rows, _excess in cases:
        labels = tuple(sorted(label for _d, label in low_rows))
        by_labels[labels].add(high_d)
    certificates = {}
    packet = []
    for labels in sorted(by_labels):
        cells = direct_complete_cells(body, LEVEL, labels, ruler)
        inherited_cells = module.vector_fixed_safe_cells(stream, labels)
        require(np.array_equal(cells, inherited_cells), (index, labels, "cell set"))
        for high_d in sorted(by_labels[labels]):
            per_residue = ruler // high_d
            lower_bound = (len(cells) + per_residue - 1) // per_residue
            kappa = (high_d + 6) // 7
            residues = np.unique(cells % high_d)
            exact_count = len(residues)
            require(exact_count > kappa, (index, labels, high_d, exact_count, kappa))
            if lower_bound > kappa:
                method = "coarse"
                slack = lower_bound - kappa
            else:
                method = "exact"
                slack = exact_count - kappa
            certificates[(labels, high_d)] = (method, slack)
            packet.append((labels, high_d, len(cells), lower_bound, exact_count, kappa,
                           slack, method,
                           hashlib.sha256(repr(tuple(map(int, residues))).encode()).hexdigest()))
    coarse = exact = 0
    minimum_slack = None
    for _mask, high_d, low_rows, _excess in cases:
        labels = tuple(sorted(label for _d, label in low_rows))
        method, slack = certificates[(labels, high_d)]
        if method == "coarse":
            coarse += 1
        else:
            exact += 1
        minimum_slack = slack if minimum_slack is None else min(minimum_slack, slack)
    require(coarse + exact == len(cases), (index, coarse, exact, len(cases)))
    terminal = module.terminal_probe((LEVEL, body, residual))
    require((len(zero), len(cases), len(by_labels), coarse, exact, minimum_slack) ==
            (terminal[7], terminal[8], terminal[9], terminal[10], terminal[11], terminal[16]),
            (index, "terminal aggregate"))
    require(terminal[12] == terminal[14] == 0 and terminal[19] and terminal[20],
            (index, "terminal closure"))
    semantic = (*terminal[:6], *terminal[7:])
    return {
        "gap": gap,
        "zero": len(zero),
        "cases": len(cases),
        "label_sets": len(by_labels),
        "coarse": coarse,
        "exact": exact,
        "minimum_slack": minimum_slack,
        "terminal_sha": hashlib.sha256(repr(semantic).encode()).hexdigest(),
        "direct_packet_sha": hashlib.sha256(repr(tuple(packet)).encode()).hexdigest(),
        "ray_checks": checks,
    }


def ftext(value):
    return f"{value.numerator}/{value.denominator}"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--processes", type=int, default=2)
    args = parser.parse_args()
    require(args.processes >= 1, args.processes)
    require(lf_sha(ATLAS) == EXPECTED_ATLAS_SHA, "atlas hash")
    require(lf_sha(SOURCE_3078) == EXPECTED_SOURCE_3078_SHA, "THM3078 source hash")
    require(lf_sha(SOURCE_3139) == EXPECTED_SOURCE_3139_SHA, "THM3139 source hash")
    require(BASE_RULER == 14 * math.lcm(*range(1, 15)), BASE_RULER)
    rows = atlas_rows()
    g8 = tuple(index for index, row in enumerate(rows) if math.gcd(LEVEL, row[1]) == 8)
    require(g8 == (8, 23, 39, 57, 64, 66, 115, 142, 150, 152, 197, 238,
                   272, 277, 300, 304, 306, 338, 370), g8)
    low_cost = tuple(index for index in g8 if rows[index][1] * rows[index][4] <= 2_000_000)
    require(set(g8) - set(low_cost) == set(TARGETS), (low_cost, TARGETS))
    order = {index for index, row in enumerate(rows) if not row[3]}
    natural = {
        index for index, row in enumerate(rows)
        if row[3] and (math.gcd(LEVEL, row[1]), row[1])
        in {(24, 2352), (36, 8820), (72, 2520)}
    }
    require(len(low_cost) == 17 and len(order) == 33 and len(natural) == 47,
            (len(low_cost), len(order), len(natural)))
    require(not (set(TARGETS) & (set(low_cost) | order | natural)), "composition overlap")

    jobs = tuple((index, rows[index]) for index in TARGETS)
    if args.processes == 1:
        screens = tuple(screen_worker(job) for job in jobs)
    else:
        with mp.get_context("spawn").Pool(min(args.processes, len(jobs))) as pool:
            screens = tuple(pool.map(screen_worker, jobs))
    screens = tuple(sorted(screens))

    print("LRC14 z216 costly gcd8 independent audit")
    print(f"atlas_rows={len(rows)};z216_row_sha256={EXPECTED_Z216_ROW_SHA};targets={TARGETS}")
    print("composition=low_cost_g8:17;order:33;natural_wall:47;target_overlap:0")
    for index, screen, audits, instance_sha, components in screens:
        body, ruler, high_floor, first_d, wall = (
            screen[1], screen[2], screen[3], screen[4], screen[5]
        )
        quotient = Counter(math.lcm(*mask) // first_d for mask in screen[13])
        terminal = terminal_audit(index, screen)
        print(
            f"row={index};body={body};L={ruler};components={components};"
            f"high_floor={high_floor};cost={ruler * components};first_d={first_d};wall={wall}"
        )
        print(
            f"screen=states:{screen[9]};crude:{screen[10]};status:{screen[11]};"
            f"residual:{screen[12]};exact_duals:{audits};quotient:{tuple(sorted(quotient.items()))};"
            f"canonical_sha256:{hashlib.sha256(repr(screen[:19]).encode()).hexdigest()};"
            f"mask_sha256:{hashlib.sha256(repr(screen[13]).encode()).hexdigest()};"
            f"instance_sha256:{instance_sha}"
        )
        print(
            f"terminal=gap:{ftext(terminal['gap'])};zero:{terminal['zero']};"
            f"cases:{terminal['cases']};label_sets:{terminal['label_sets']};"
            f"coarse:{terminal['coarse']};exact:{terminal['exact']};maxgap:0;failures:0;"
            f"minimum_slack:{terminal['minimum_slack']};ray_checks:{terminal['ray_checks']};"
            f"terminal_sha256:{terminal['terminal_sha']};"
            f"direct_packet_sha256:{terminal['direct_packet_sha']}"
        )
    require((28 + 6) // 7 == 4 and 2 * ((28 - 1) // 14) + 1 == 3,
            "translated/centered hostile")
    print("hostiles=negative_Farkas_alpha_rejected;zero_high_not_used;two_high_gap_strict;"
          "full_grid_complete_cells;translated_kappa28_4_not_centered_beta28_3")
    print("consequence=373186-2=373184;z216:382wall_to_380wall;cap_stays_216")
    print("all_exact_controls=PASS")


if __name__ == "__main__":
    main()
