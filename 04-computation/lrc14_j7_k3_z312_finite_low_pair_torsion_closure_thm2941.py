#!/usr/bin/env python3
"""Close the 70 projected k=3,z1=312 states by low-pair p-torsion.

The exact ray/status frontier leaves 70 four-denominator multisets on four
bodies.  The inherited projected high-label obligation requires at least one
label above the high wall.  A deliberately duplicate-permitting scalar
relaxation shows that a packet with two or three such labels cannot attain the
scalar lower wall.  Hence every possible packet in the inherited high-wall
universe consists of the fixed first label 312, two labels below the wall, and
one high label.  An unrestricted zero-high hostile control actually passes
the scalar wall in 69 of 70 states, so the inherited obligation is essential.

The two low slots are finite.  Exact residue-ray upper bounds reduce all 70
states to 80 (state, high denominator, low pair) cases, involving only 20
body-distinct low pairs.  For every case there are two complete body-carrier
cells missed by 312 and both low labels whose indices differ by d/p modulo the
exact high denominator d, with p=2 in 78 cases and p=3 in two cases.  Writing
the high label as (L/d)u+mL, gcd(u,d)=1, its phases on those cells differ by
u/p modulo one.  This is a half turn for p=2 and a third turn for p=3, so the
two strict radius-1/14 danger arcs are disjoint for every unit u and every
height m.  The projected drift-safe residual therefore has full mass one,
strictly above the three-aligned completion cap 36/91.

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


ROOT = Path(__file__).resolve().parents[1]
FRONTIER_PATH = (
    ROOT / "04-computation/lrc14_j7_k3_z312_ray_status_frontier_thm2941.py"
)
FRONTIER_OUTPUT_PATH = (
    ROOT / "05-knowledge/results/lrc14_j7_k3_z312_ray_status_frontier_thm2941.out"
)
PROJECTION_PATH = (
    ROOT / "04-computation/lrc14_j7_five_aligned_two_drift_projected_closure_thm2941.py"
)
BODY_ATLAS_PATH = (
    ROOT / "04-computation/lrc14_j7_k3_projected_scalar_body_atlas_thm2941.py"
)
BODY_ATLAS_OUTPUT_PATH = (
    ROOT / "05-knowledge/results/lrc14_j7_k3_projected_scalar_body_atlas_thm2941.out"
)
DEFAULT_OUTPUT = (
    ROOT / "05-knowledge/results/"
    "lrc14_j7_k3_z312_finite_low_pair_torsion_closure_thm2941.out"
)

EXPECTED_FRONTIER_SOURCE_SHA256 = (
    "d6a80c15c5c4d8ef2ea8be9fc886c40e70e3189123b5d0b3fce48765fa301977"
)
EXPECTED_FRONTIER_OUTPUT_SHA256 = (
    "d03fc39ed1f5f64cd2be4e7e28f5cf23e8d7adc0c6a737abc6944bdb7672515f"
)
EXPECTED_PROJECTION_SHA256 = (
    "76f891edfcc029a08202481304a809e03e8bd81f247afaeabab685825c4d3662"
)
EXPECTED_BODY_ATLAS_SHA256 = (
    "2af6d96882f336a409a8657070ed76a75c09a53b3789101b83103b051e864ded"
)
EXPECTED_BODY_ATLAS_OUTPUT_SHA256 = (
    "cee82237ce1f51729813b9c916edd3353204c18172abe1d71278dee2c5562eda"
)
EXPECTED_SEMANTIC_SHA256 = (
    "95f9da4ccf85145f974df62d15c79acfa7000e0917e877bb5b32926e4ff6c3e8"
)

FIRST = 312
ALIGNED_THREE_UNION_CAP = F(36, 91)
EXPECTED_BODY_CASES = {
    (1, 8, 10, 11, 12, 14): (26, 6, F(271403663, 168333225060)),
    (1, 8, 10, 12, 13, 14): (50, 12, F(22539649297, 15003760917600)),
    (1, 8, 11, 12, 13, 14): (2, 1, F(295936150144, 96567353117235)),
    (2, 8, 10, 11, 12, 14): (2, 1, F(41681149, 16799037255)),
}
EXPECTED_TORSION_HISTOGRAM = {2: 78, 3: 2}
EXPECTED_ZERO_HIGH_PASS_COUNTS = {
    (1, 8, 10, 11, 12, 14): 24,
    (1, 8, 10, 12, 13, 14): 41,
    (1, 8, 11, 12, 13, 14): 2,
    (2, 8, 10, 11, 12, 14): 2,
}


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def file_sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, f"cannot load {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


for path, expected in (
    (FRONTIER_PATH, EXPECTED_FRONTIER_SOURCE_SHA256),
    (FRONTIER_OUTPUT_PATH, EXPECTED_FRONTIER_OUTPUT_SHA256),
    (PROJECTION_PATH, EXPECTED_PROJECTION_SHA256),
    (BODY_ATLAS_PATH, EXPECTED_BODY_ATLAS_SHA256),
    (BODY_ATLAS_OUTPUT_PATH, EXPECTED_BODY_ATLAS_OUTPUT_SHA256),
):
    if expected is not None:
        require(file_sha256(path) == expected, ("dependency changed", path))

frontier = load_module("z312_torsion_frontier", FRONTIER_PATH)
projection = load_module("z312_torsion_projection", PROJECTION_PATH)
ray = frontier.ray


def first_on_ray(residue, modulus, threshold):
    if residue >= threshold:
        return residue
    return residue + ((threshold - residue + modulus - 1) // modulus) * modulus


def ray_rows(stream, denominator):
    """All finite low labels and one rigorous high supremum per unit ray."""
    low = []
    high = []
    step = stream.L // denominator
    for unit in range(1, denominator):
        if math.gcd(unit, denominator) != 1:
            continue
        residue = step * unit
        amplitude = residue * ray.local.delta(stream.carrier, stream.h, residue)
        require(
            (residue + stream.L)
            * ray.local.delta(stream.carrier, stream.h, residue + stream.L)
            == amplitude,
            ("ray recurrence failed", stream.body, denominator, unit),
        )
        label = first_on_ray(residue, stream.L, FIRST + 1)
        while label < stream.high_floor:
            require(
                stream.L // math.gcd(stream.L, label) == denominator,
                ("low denominator changed", label, denominator),
            )
            low.append((amplitude / label, label, unit))
            label += stream.L
        high_label = first_on_ray(residue, stream.L, stream.high_floor)
        require(
            stream.L // math.gcd(stream.L, high_label) == denominator,
            ("high denominator changed", high_label, denominator),
        )
        # Positive rays decrease and attain this maximum.  Zero rays attain
        # zero.  Negative rays increase to the strict supremum zero.  Using
        # zero is therefore a valid (occasionally deliberately unattained)
        # upper bound.
        high_upper = amplitude / high_label if amplitude > 0 else F(0)
        high.append((high_upper, high_label, unit, amplitude))
    require(high, ("empty exact denominator class", denominator))
    return (
        tuple(sorted(low, key=lambda row: (row[0], -row[1], -row[2]), reverse=True)),
        tuple(sorted(high, key=lambda row: (row[0], -row[1], -row[2]), reverse=True)),
    )


def relaxed_two_high_bound(stream, residuals, tables):
    """Duplicate-permitting upper bound for packets with >=2 high labels."""
    best = None
    for ds in residuals:
        suffix = list(ds)
        suffix.remove(stream.first_d)
        for high_indices in ((0, 1), (0, 2), (1, 2), (0, 1, 2)):
            total = F(0)
            valid = True
            for index, denominator in enumerate(suffix):
                rows = tables[denominator][1 if index in high_indices else 0]
                if not rows:
                    valid = False
                    break
                # Slots are optimized independently, even when this repeats a
                # literal label.  Thus this is an upper bound on real packets.
                total += rows[0][0]
            if valid:
                candidate = (total, ds, tuple(suffix), high_indices)
                if best is None or candidate > best:
                    best = candidate
    require(best is not None, ("missing two-high relaxation", stream.body))
    return best


def relaxed_zero_high_bound(stream, residuals, tables):
    """Unrestricted finite-low upper, retained as a hostile scope control."""
    best = None
    rows = []
    required = stream.lower - stream.first_delta
    for ds in residuals:
        suffix = list(ds)
        suffix.remove(stream.first_d)
        if not all(tables[denominator][0] for denominator in suffix):
            continue
        total = sum((tables[denominator][0][0][0] for denominator in suffix), F(0))
        candidate = (
            total,
            ds,
            tuple((denominator, tables[denominator][0][0]) for denominator in suffix),
        )
        rows.append((ds, total, total - required, candidate[2]))
        if best is None or candidate > best:
            best = candidate
    require(best is not None, ("missing zero-high hostile control", stream.body))
    rows = tuple(sorted(rows))
    return best, rows, sum(excess >= 0 for _ds, _total, excess, _labels in rows)


def finite_low_pairs(first, second, *, same, threshold):
    """Every distinct finite low pair whose sum reaches ``threshold``."""
    result = []
    if same:
        for index, left in enumerate(first):
            if index + 1 >= len(first) or left[0] + first[index + 1][0] < threshold:
                break
            for right in first[index + 1 :]:
                if left[0] + right[0] < threshold:
                    break
                if left[1] != right[1]:
                    result.append((left, right))
    else:
        for left in first:
            if not second or left[0] + second[0][0] < threshold:
                break
            for right in second:
                if left[0] + right[0] < threshold:
                    break
                if left[1] != right[1]:
                    result.append((left, right))
    return tuple(result)


def fixed_safe_cells(cells, modulus, labels):
    """Whole body cells missing all strict danger combs in ``labels``."""
    require(all(0 < label < modulus for label in labels), ("nonlocal fixed label", labels))
    safe = []
    for cell in cells:
        for label in labels:
            residue = (label * cell) % modulus
            # The phase segment is [residue,residue+label]/L.  It misses the
            # open radius-1/14 neighborhood of every integer exactly under
            # these two weak endpoint inequalities.
            if 14 * residue < modulus or 14 * (residue + label) > 13 * modulus:
                break
        else:
            safe.append(cell)
    return tuple(safe)


def torsion_pair(safe_cells, denominator):
    """Find a deterministic 2- or 3-torsion pair modulo ``denominator``."""
    by_residue = defaultdict(list)
    for cell in safe_cells:
        by_residue[cell % denominator].append(cell)
    for prime in (2, 3):
        if denominator % prime:
            continue
        shift = denominator // prime
        for residue, first_cells in by_residue.items():
            second_cells = by_residue.get((residue + shift) % denominator)
            if second_cells:
                return first_cells[0], second_cells[0], shift, prime
    return None


def atlas_below_frontier():
    counts = Counter()
    pattern = re.compile(r"^row=E=[0-9,]+;.*;z1=([0-9]+);")
    for line in BODY_ATLAS_OUTPUT_PATH.read_text().splitlines():
        match = pattern.match(line)
        if match:
            counts[int(match.group(1))] += 1
    require(counts[FIRST] == 80, ("z312 atlas count changed", counts[FIRST]))
    require(all(counts[value] == 0 for value in range(307, FIRST)), "atlas gap changed")
    require(counts[306] == 2, ("next occupied height changed", counts[306]))
    return tuple(sorted(counts.items())), 306, counts[306]


def analyze_body(body, residuals):
    stream = ray.Stream(body)
    require(
        stream.first == FIRST < stream.high_floor,
        ("projected high-wall gate changed", stream.first, stream.high_floor),
    )
    required = stream.lower - stream.first_delta
    denominators = sorted({denominator for ds in residuals for denominator in ds})
    tables = {denominator: ray_rows(stream, denominator) for denominator in denominators}

    zero_relaxed = relaxed_zero_high_bound(stream, residuals, tables)
    zero_high_gap = required - zero_relaxed[0][0]
    two_relaxed = relaxed_two_high_bound(stream, residuals, tables)
    two_high_gap = required - two_relaxed[0]
    require(two_high_gap > 0, ("two-high packet survived", body, two_relaxed))

    cells = projection.body_cells(stream.carrier, stream.L)
    safe_cache = {}
    cases = []
    for ds in residuals:
        suffix = list(ds)
        suffix.remove(stream.first_d)
        for high_denominator in sorted(set(suffix)):
            low_denominators = list(suffix)
            low_denominators.remove(high_denominator)
            low_first, low_second = sorted(low_denominators)
            high_row = tables[high_denominator][1][0]
            high_upper, high_representative, high_unit, high_amplitude = high_row
            pairs = finite_low_pairs(
                tables[low_first][0],
                tables[low_second][0],
                same=low_first == low_second,
                threshold=required - high_upper,
            )
            for left, right in pairs:
                low_labels = tuple(sorted((left[1], right[1])))
                cache_key = low_labels
                if cache_key not in safe_cache:
                    safe_cache[cache_key] = fixed_safe_cells(
                        cells, stream.L, (FIRST, *low_labels)
                    )
                certificate = torsion_pair(safe_cache[cache_key], high_denominator)
                require(certificate is not None, ("no torsion pair", body, ds, high_denominator, low_labels))
                cell_j, cell_k, shift, prime = certificate
                require(
                    (cell_k - cell_j) % high_denominator == shift
                    and prime * shift == high_denominator,
                    ("torsion congruence changed", certificate, high_denominator),
                )
                fixed_labels = (FIRST, *low_labels)
                fixed_dangers = tuple(
                    (cell, label, projection.phase_danger(cell, label, stream.L))
                    for cell in (cell_j, cell_k)
                    for label in fixed_labels
                )
                require(
                    all(not intervals for _cell, _label, intervals in fixed_dangers),
                    ("fixed label entered located cell", body, ds, fixed_dangers),
                )
                require(
                    math.gcd(high_unit, high_denominator) == 1
                    and high_denominator % prime == 0
                    and high_unit % prime != 0,
                    ("representative lost unit torsion", high_row, certificate),
                )
                high_j = projection.phase_danger(cell_j, high_representative, stream.L)
                high_k = projection.phase_danger(cell_k, high_representative, stream.L)
                require(
                    not projection.intersect_fraction(high_j, high_k),
                    ("representative high dangers intersect", body, ds, high_row, certificate),
                )
                require(F(1, prime) > F(1, 7), ("torsion gap too small", prime))
                cases.append(
                    (
                        ds,
                        high_denominator,
                        low_labels,
                        left[0] + right[0],
                        high_upper,
                        high_representative,
                        high_unit,
                        high_amplitude,
                        len(safe_cache[cache_key]),
                        certificate,
                    )
                )

    cases = tuple(sorted(cases))
    distinct_low_pairs = tuple(sorted(safe_cache))
    return (
        body,
        stream.L,
        stream.high_floor,
        required,
        zero_relaxed,
        zero_high_gap,
        two_relaxed,
        two_high_gap,
        len(cells),
        distinct_low_pairs,
        cases,
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--processes", type=int, default=4)
    args = parser.parse_args()
    require(args.processes >= 1, "process count must be positive")

    require(
        set(frontier.EXPECTED_RESIDUALS) == set(EXPECTED_BODY_CASES),
        "frontier residual bodies changed",
    )
    tasks = tuple(frontier.EXPECTED_RESIDUALS.items())
    if args.processes == 1:
        records = tuple(analyze_body(body, residuals) for body, residuals in tasks)
    else:
        with mp.get_context("spawn").Pool(min(args.processes, len(tasks))) as pool:
            records = tuple(pool.starmap(analyze_body, tasks))
    for record in records:
        (
            body, _L, _high, _required, _zero_relaxed, _zero_gap,
            _two_relaxed, gap, _cell_count, low_pairs, cases,
        ) = record
        expected_cases, expected_pairs, expected_gap = EXPECTED_BODY_CASES[body]
        require(
            (len(cases), len(low_pairs), gap)
            == (expected_cases, expected_pairs, expected_gap),
            ("body closure ledger changed", body, len(cases), len(low_pairs), gap),
        )
        require(
            _zero_relaxed[2] == EXPECTED_ZERO_HIGH_PASS_COUNTS[body],
            ("zero-high hostile ledger changed", body, _zero_relaxed[2]),
        )
    all_cases = tuple((record[0], case) for record in records for case in record[10])
    require(len(all_cases) == 80, ("global case count changed", len(all_cases)))
    zero_high_passes = sum(record[4][2] for record in records)
    require(zero_high_passes == 69, ("zero-high global hostile ledger changed", zero_high_passes))
    torsion_histogram = Counter(case[-1][-1] for _body, case in all_cases)
    require(
        dict(sorted(torsion_histogram.items())) == EXPECTED_TORSION_HISTOGRAM,
        ("torsion histogram changed", torsion_histogram),
    )
    atlas_counts, next_height, next_count = atlas_below_frontier()
    require(ALIGNED_THREE_UNION_CAP < 1, "aligned cap ceased to be proper")

    semantic_payload = (
        FIRST,
        records,
        atlas_counts,
        next_height,
        next_count,
        ALIGNED_THREE_UNION_CAP,
        tuple(sorted(torsion_histogram.items())),
    )
    semantic_hash = hashlib.sha256(repr(semantic_payload).encode()).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(semantic_hash == EXPECTED_SEMANTIC_SHA256, "semantic digest changed")

    lines = [
        "LRC14 projected k=3 z1=312 finite-low-pair p-torsion closure",
        f"frontier_source_sha256={file_sha256(FRONTIER_PATH)}",
        f"frontier_output_sha256={file_sha256(FRONTIER_OUTPUT_PATH)}",
        f"projection_source_sha256={file_sha256(PROJECTION_PATH)}",
        f"body_atlas_source_sha256={file_sha256(BODY_ATLAS_PATH)}",
        f"body_atlas_output_sha256={file_sha256(BODY_ATLAS_OUTPUT_PATH)}",
        "inherited_frontier=80 bodies;18249 states;7481 crude kills;10698 status kills;70 residuals on 4 bodies",
        "high_wall_gate=the pinned projected scalar atlas/frontier requires at least one label >= high_floor because first=312<high_floor on every residual body",
        f"hostile_zero_high_control=the unrestricted zero-high scalar relaxation passes {zero_high_passes}/70 states;zero-high packets are excluded only by the inherited projected high-label obligation",
        "scalar_split=within that high-wall universe,the duplicate-permitting >=2-high upper is strictly below the wall on every residual body;therefore exactly one high label",
    ]
    for record in records:
        (
            body, L, high_floor, required, zero_relaxed, zero_gap,
            two_relaxed, gap, cell_count, low_pairs, cases,
        ) = record
        lines.append(
            f"body={body};L={L};high_floor={high_floor};required_suffix={required};"
            f"zero_high_unrestricted_upper={zero_relaxed[0][0]};zero_high_gap:{zero_gap};"
            f"zero_high_witness=ds:{zero_relaxed[0][1]};labels:{zero_relaxed[0][2]};"
            f"zero_high_scalar_passes:{zero_relaxed[2]}/{len(zero_relaxed[1])};"
            f"two_high_relaxed_upper={two_relaxed[0]};two_high_gap={gap};"
            f"two_high_witness=ds:{two_relaxed[1]};suffix:{two_relaxed[2]};high_slots:{two_relaxed[3]};"
            f"body_cells={cell_count};finite_cases={len(cases)};distinct_low_pairs={len(low_pairs)}"
        )
    lines.extend(
        (
            "finite_reduction=80 necessary (state,high-denominator,low-pair) cases;20 body-distinct low pairs",
            f"torsion_histogram={dict(sorted(torsion_histogram.items()))}",
            "symbolic_law=z=(L/d)u+mL,gcd(u,d)=1,k-j=d/p mod d => phase_k-phase_j=u/p mod 1",
            "topology=p=2 gives a half turn;p=3 gives a one-third/two-thirds turn;strict radius-1/14 danger arcs are disjoint",
            f"aligned_completion=projected drift-safe mass 1 > three-aligned cap {ALIGNED_THREE_UNION_CAP}",
        )
    )
    for body, case in all_cases:
        (
            ds, high_d, lows, low_sum, high_upper, high_label, high_unit,
            high_amplitude, safe_count, certificate,
        ) = case
        lines.append(
            f"case=E:{body};ds:{ds};high_d:{high_d};lows:{lows};"
            f"low_sum:{low_sum};high_upper:{high_upper};high_rep:{high_label};"
            f"high_unit:{high_unit};high_amplitude:{high_amplitude};"
            f"fixed_safe_cells:{safe_count};certificate:{certificate}"
        )
    lines.extend(
        (
            f"atlas_gap=307..311 empty;next_occupied_height:{next_height};bodies:{next_count}",
            "conclusion=all 70 z1=312 residual denominator states are empty uniformly;all 80 z1=312 scalar-atlas rows close;projected k3 cap<=306;next exact occupied frontier=306 with 2 bodies",
            f"semantic_sha256={semantic_hash}",
            "all_exact_controls=PASS",
        )
    )
    payload = "\n".join(lines) + "\n"
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(payload, encoding="utf-8", newline="\n")
    print(payload, end="")


if __name__ == "__main__":
    main()
