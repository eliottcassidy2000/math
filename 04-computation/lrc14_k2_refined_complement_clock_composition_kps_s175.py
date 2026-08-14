#!/usr/bin/env python3
"""Compose THM-3366's pool-14 terminals with the refined k=2 ledger.

The raw support census and the later septimal/d=6 screens are not disjoint
counts, so their numerical deletions cannot be subtracted independently.
This referee reconstructs the exact post-d=6 occurrence count on every
(body,D) row, independently constructs the pool-14 complement-clock terminal
set, and intersects the two ledgers before reporting a new residual.
"""

from collections import Counter
from fractions import Fraction as Q
from functools import lru_cache
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
from shutil import which
from subprocess import PIPE, run
from tempfile import TemporaryDirectory


ROOT = Path(__file__).resolve().parents[1]
D6_PATH = ROOT / "04-computation/lrc14_k2_d6_located_phase_closure_thm2928.py"
COMPLEMENT_BASE_PATH = (
    ROOT / "04-computation/lrc14_k1_body_complement_clock_scan_kps_s172.py"
)
EXPECTED_D6_SHA256 = "9f300459b273ad1825d3fe3e9274c6afe609f2d581e9df3d2be1780d347e541b"
EXPECTED_COMPLEMENT_BASE_LF_SHA256 = (
    "bdb2001cf22f7e92884e895b0095021e42e8f1febd9adbf779b250a2f6c53507"
)
POOL = tuple(range(1, 15))
EXPECTED_COMPLEMENT_ROWS = 19_198
EXPECTED_COMPLEMENT_SEMANTIC = (
    "7ea43d80c994de6f645edebd9872024b3106e4cf902c1518c03dd8736e05ef19"
)
EXPECTED_PRE = (26_899_164_786, 200_141_092_521, 4_354, 2_966, 147)
EXPECTED_COMPLETION_HISTOGRAM = ((1, 12_662), (2, 2_764), (3, 998), (4, 1_106), (5, 1_668))
EXPECTED_REFINED_DELETION = (298, 71_575_318)
EXPECTED_POST = (26_899_164_065, 200_069_517_203, 4_056, 2_966, 144)
EXPECTED_POST_BY_C = {
    1: 3_524_756,
    2: 46_320_209_782,
    3: 112_778_172_674,
    4: 39_482_116_896,
    5: 1_485_493_095,
}
EXPECTED_PRE_SEMANTIC = "a0835d54b3ac720e57d5788945bed44a2fda134897918415d8df610634f7a973"
EXPECTED_POST_SEMANTIC = "032b96bf1d9877a68c395e9508f51c20e5476cf7500891d9c5b6920bc0697238"


def require(condition, detail):
    if not condition:
        raise RuntimeError(detail)


def file_hash(path):
    return sha256(path.read_bytes()).hexdigest()


def normalized_hash(path):
    return sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def load_module(name, path):
    spec = spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, ("module import", path))
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def target_mask(base, D, gaps, points, atoms):
    target = 0
    for index, x in enumerate(points):
        if (D * x).denominator != 1 and any(left < x < right for left, right in gaps):
            target |= 1 << index
    offset = len(points)
    for index, (left, right) in enumerate(atoms):
        if any(
            max(left, gap_left) < min(right, gap_right)
            for gap_left, gap_right in gaps
        ):
            target |= 1 << (offset + index)
    return target


def build_solver(base):
    points = base.arrangement_points(POOL)
    atoms = tuple(zip(points, points[1:]))
    samples = points + tuple((left + right) / 2 for left, right in atoms)
    masks = tuple(
        sum((1 << i for i, x in enumerate(samples) if base.danger(c, x)), 0)
        for c in POOL
    )
    union = 0
    for mask in masks:
        union |= mask
    by_bit = tuple(
        tuple(i for i, mask in enumerate(masks) if mask & (1 << bit))
        for bit in range(len(samples))
    )

    @lru_cache(maxsize=None)
    def exact(remaining, depth):
        if remaining == 0:
            return ()
        if depth == 0:
            return None
        bits = tuple(i for i in range(len(samples)) if remaining & (1 << i))
        pivot = min(bits, key=lambda i: len(by_bit[i]))
        for candidate in by_bit[pivot]:
            reduced = remaining & ~masks[candidate]
            if reduced == remaining:
                continue
            suffix = exact(reduced, depth - 1)
            if suffix is not None:
                return (POOL[candidate],) + suffix
        return None

    def solve(target):
        if target & ~union:
            return None
        for depth in range(6):
            result = exact(target, depth)
            if result is not None:
                return tuple(sorted(result))
        return None

    return points, atoms, solve, exact


def run_engine_portable(base, queries):
    """Replay the frozen C++ engine with the available host compiler."""
    compiler = which("g++") or which("clang++")
    require(compiler is not None, "no C++ compiler")
    input_lines = [str(len(queries))]
    for D in sorted(queries):
        input_lines.append(f"{D} {len(queries[D])}")
        input_lines.extend(f"{c} {need}" for c, need in queries[D])
    input_text = "\n".join(input_lines) + "\n"
    outputs = []
    with TemporaryDirectory(prefix="lrc14-k2-composition-") as temporary:
        for optimization in ("-O2", "-O3"):
            executable = Path(temporary) / f"engine-{optimization[1:]}.exe"
            compiled = run(
                [
                    compiler,
                    "-std=c++17",
                    optimization,
                    "-DNDEBUG",
                    str(base.ENGINE_PATH),
                    "-o",
                    str(executable),
                ],
                stdout=PIPE,
                stderr=PIPE,
                text=True,
                check=False,
            )
            require(compiled.returncode == 0, ("engine compile", compiled.stderr))
            result = run(
                [str(executable)],
                input=input_text,
                stdout=PIPE,
                stderr=PIPE,
                text=True,
                check=False,
            )
            require(result.returncode == 0, ("engine run", result.stderr))
            outputs.append(result.stdout)
    require(outputs[0] == outputs[1], "O2/O3 engine mismatch")
    output = outputs[0]
    require(
        sha256(output.encode()).hexdigest() == base.EXPECTED_ENGINE_OUTPUT_SHA256,
        "engine transcript changed",
    )
    raw = {}
    answers = {}
    for line in output.splitlines():
        fields = line.split()
        if not fields:
            continue
        if fields[0] == "C":
            _tag, D, c, count, state_count, minimum, maximum = fields
            raw[(int(D), int(c))] = (
                int(count),
                int(state_count),
                int(minimum),
                int(maximum),
            )
        elif fields[0] == "Q":
            _tag, D, c, need, count, first = fields
            answers[(int(D), int(c), int(need))] = (int(count), int(first))
    require(len(raw) == 5 * len(queries), "engine raw parse")
    require(len(answers) == sum(map(len, queries.values())), "engine query parse")
    return raw, answers, output


def aggregate(d6, rows_by_D, answers, allowed_keys=None):
    occurrences = 0
    rows = set()
    bodies = set()
    divisors = set()
    occurrence_by_c = Counter()
    row_counts = {}
    semantic = sha256()

    for D in sorted(rows_by_D):
        q = D // 7
        coefficient = d6.fixed_d6_shape_count(D)
        for body, L, support_count, histogram, N_values in rows_by_D[D]:
            key = (body, D)
            if allowed_keys is not None and key not in allowed_keys:
                continue
            row_total = 0
            N4 = N_values[3]
            decrement = coefficient if coefficient and 1 <= N4 <= q // 6 else 0
            for c, N in enumerate(N_values, 1):
                count = answers[(D, c, N)][0]
                if c == 4:
                    count -= decrement
                require(count >= 0, ("negative occurrence", D, body, c, count))
                row_total += count
                occurrence_by_c[c] += count
            if not row_total:
                continue
            row_counts[key] = row_total
            occurrences += row_total
            rows.add(key)
            bodies.add(body)
            divisors.add(D)
            semantic.update(
                f"{D}|{body}|{L}|{support_count}|{histogram}|{N_values}|{row_total}\n".encode()
            )

    shapes = 0
    for D in sorted({D for _body, D in rows}):
        q = D // 7
        coefficient = d6.fixed_d6_shape_count(D)
        local_rows = [
            row
            for row in rows_by_D[D]
            if (row[0], D) in rows
        ]
        for c in range(1, 6):
            minimum_N = min(row[4][c - 1] for row in local_rows)
            count = answers[(D, c, minimum_N)][0]
            if c == 4 and coefficient and 1 <= minimum_N <= q // 6:
                count -= coefficient
            require(count >= 0, ("negative shape", D, c, count))
            shapes += count

    return {
        "summary": (shapes, occurrences, len(rows), len(bodies), len(divisors)),
        "rows": rows,
        "row_counts": row_counts,
        "occurrence_by_c": occurrence_by_c,
        "semantic": semantic.hexdigest(),
    }


def main():
    require(file_hash(D6_PATH) == EXPECTED_D6_SHA256, "d6 dependency changed")
    require(
        normalized_hash(COMPLEMENT_BASE_PATH) == EXPECTED_COMPLEMENT_BASE_LF_SHA256,
        "complement base changed",
    )
    d6 = load_module("kps_s175_d6", D6_PATH)
    complement = load_module("kps_s175_complement", COMPLEMENT_BASE_PATH)

    rows_by_D, _nonseptimal = d6.base.build_rows()
    queries = d6.base.engine_queries(rows_by_D)
    _raw, answers, engine_output = run_engine_portable(d6.base, queries)

    points, atoms, solve, exact = build_solver(complement)
    closed = set()
    completion_semantic = sha256()
    completion_histogram = Counter()
    for D in sorted(rows_by_D):
        for body, L, support_count, _histogram, _N_values in rows_by_D[D]:
            check_L, ranges = d6.base.support.safe_cell_ranges(body)
            require(check_L == L, ("body ruler", body, L, check_L))
            arcs = complement.residue_arcs(D, ranges)
            require(
                sum(right - left for left, right in arcs) == support_count,
                ("support arcs", body, D),
            )
            gaps = complement.unsupported_gaps(D, arcs)
            target = target_mask(complement, D, gaps, points, atoms)
            completion = solve(target)
            if completion is None:
                continue
            key = (body, D)
            closed.add(key)
            completion_histogram[len(completion)] += 1
            completion_semantic.update(
                f"{body}|{L}|{D}|{support_count}|{completion}\n".encode()
            )

    require(len(closed) == EXPECTED_COMPLEMENT_ROWS, len(closed))
    require(
        tuple(sorted(completion_histogram.items())) == EXPECTED_COMPLETION_HISTOGRAM,
        ("completion histogram", tuple(sorted(completion_histogram.items()))),
    )
    require(
        completion_semantic.hexdigest() == EXPECTED_COMPLEMENT_SEMANTIC,
        (
            "complement semantic changed",
            completion_semantic.hexdigest(),
            EXPECTED_COMPLEMENT_SEMANTIC,
        ),
    )

    before = aggregate(d6, rows_by_D, answers)
    require(before["summary"] == EXPECTED_PRE, ("pre-composition", before["summary"]))
    allowed = before["rows"] - closed
    after = aggregate(d6, rows_by_D, answers, allowed)
    killed_rows = before["rows"] & closed
    killed_occurrences = sum(before["row_counts"][key] for key in killed_rows)
    require(
        before["summary"][1] - after["summary"][1] == killed_occurrences,
        "occurrence subtraction",
    )
    require(before["semantic"] == EXPECTED_PRE_SEMANTIC, "pre semantic changed")
    require(
        (len(killed_rows), killed_occurrences) == EXPECTED_REFINED_DELETION,
        ("refined deletion", len(killed_rows), killed_occurrences),
    )
    require(after["summary"] == EXPECTED_POST, ("post summary", after["summary"]))
    require(
        dict(after["occurrence_by_c"]) == EXPECTED_POST_BY_C,
        ("post by-c", dict(after["occurrence_by_c"])),
    )
    require(after["semantic"] == EXPECTED_POST_SEMANTIC, "post semantic changed")

    print("LRC14 K=2 REFINED COMPLEMENT-CLOCK COMPOSITION")
    print("status=FINITE-EXACT intersection of post-d6 and pool-14 ledgers")
    print(
        f"arrangement_points={len(points)};atoms={len(atoms)};"
        f"support_closed_rows={len(closed)};"
        f"support_completion_histogram={tuple(sorted(completion_histogram.items()))}"
    )
    print(f"pre_summary={before['summary']}")
    print(f"pre_occurrence_by_c={before['occurrence_by_c']}")
    print(
        f"refined_closed_rows={len(killed_rows)};"
        f"refined_closed_occurrences={killed_occurrences}"
    )
    print(f"post_summary={after['summary']}")
    print(f"post_occurrence_by_c={after['occurrence_by_c']}")
    print(f"completion_semantic_sha256={completion_semantic.hexdigest()}")
    print(f"pre_semantic_sha256={before['semantic']}")
    print(f"post_semantic_sha256={after['semantic']}")
    print(f"engine_output_sha256={sha256(engine_output.encode()).hexdigest()}")
    print(f"solver_cache={exact.cache_info()}")
    print("verdict=PASS")


if __name__ == "__main__":
    main()
