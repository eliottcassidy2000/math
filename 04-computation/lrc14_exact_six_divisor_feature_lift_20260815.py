#!/usr/bin/env python3
"""Exact divisor and feature lifts of the refined six-clock mutation graph.

The body-only exact-six relation forgets the divisor row which produced an
edge.  This companion first reconstructs the complete current k=3 residual,
enumerates every pool-14 six-completion on its 1,897 body/divisor keys, and
then retains the divisor while transporting the current denominator-feature
class.  It also revisits the sole refined k=2 body three-cycle on the smallest
same-divisor state space needed to decide whether that cycle lifts.

The resulting graphs are finite exact bookkeeping objects.  A surviving
feature class is still an aggregate class, not a literal denominator tuple,
common phase, next-sector row, ancestry packet, or physical iteration.
"""

from __future__ import annotations

import ast
from collections import Counter, defaultdict
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations
from pathlib import Path
from shutil import which
from subprocess import PIPE, run
from tempfile import TemporaryDirectory


ROOT = Path(__file__).resolve().parents[1]
K3_COMPOSITION_PATH = ROOT / "04-computation/lrc14_k3_refined_complement_clock_composition_kps_s176.py"
K3_LEDGER_PATH = ROOT / "04-computation/lrc14_k3_four_drift_expected_spike_gf_thm2928.py"
COMPLEMENT_PATH = ROOT / "04-computation/lrc14_k1_body_complement_clock_scan_kps_s172.py"
K2_D6_PATH = ROOT / "04-computation/lrc14_k2_d6_located_phase_closure_thm2928.py"

EXPECTED_DEPENDENCIES = {
    K3_COMPOSITION_PATH: "27e4bff52705189bf8ff73db42d76d4e2fc94c44330d295f166f8d4217cb1804",
    K3_LEDGER_PATH: "05e365a654b32e66b814dcbce9385a2d13c22a2c84a5474e0855dcab6262b055",
    COMPLEMENT_PATH: "bdb2001cf22f7e92884e895b0095021e42e8f1febd9adbf779b250a2f6c53507",
    K2_D6_PATH: "9f300459b273ad1825d3fe3e9274c6afe609f2d581e9df3d2be1780d347e541b",
}

POOL = tuple(range(1, 15))
A = (1, 3, 7, 8, 9, 10)
B = (2, 4, 5, 7, 9, 12)
C = (4, 5, 7, 8, 9, 10)
K2_KEYS = ((A, 17_640), (C, 17_640), (B, 17_640), (B, 4_410), (A, 4_410))
EXPECTED_K3_POST = (398_241_574, 2_548_893_834, 1_897, 1_823, 107)
EXPECTED_K3_CLASSES = (1_823, 20, 54)
EXPECTED_K3_STATUS_HISTOGRAM = (("absent_current_key", 8), ("live_self", 12))
EXPECTED_K3_FEATURE_SUMMARY = (12, 38, 298_227)
EXPECTED_LONGEST_BODY_PATH = (
    (1, 2, 3, 7, 10, 12),
    (1, 5, 6, 7, 8, 10),
    (3, 4, 5, 7, 8, 12),
)
EXPECTED_K2_ENGINE_OUTPUT_SHA256 = "00e833b13be05d6708e86e9426e158bf05b63641354f875ca13fbc87a3c3dee1"
EXPECTED_SEMANTIC_SHA256 = "c423bca65a43fdbfddc25b972a373002890f2c0a6f42c6aea96a8c2e2dacb8e5"


def require(condition, detail):
    if not condition:
        raise RuntimeError(detail)


def lf_sha256(path):
    return sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def load_module(name, path):
    spec = spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, ("module", path))
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def arrangement(complement):
    points = complement.arrangement_points(POOL)
    atoms = tuple(zip(points, points[1:]))
    samples = points + tuple((left + right) / 2 for left, right in atoms)
    masks = {
        clock: sum(
            (1 << index for index, value in enumerate(samples)
             if complement.danger(clock, value)),
            0,
        )
        for clock in POOL
    }
    six_unions = tuple(
        (
            body,
            combine_masks(masks[clock] for clock in body),
        )
        for body in combinations(POOL, 6)
    )
    return points, atoms, six_unions


def combine_masks(values):
    result = 0
    for value in values:
        result |= value
    return result


def target_for_row(composition, complement, divisor, arcs, points, atoms):
    gaps = complement.unsupported_gaps(divisor, arcs)
    return composition.target_mask(complement, divisor, gaps, points, atoms)


def six_completions(target, six_unions):
    return tuple(body for body, mask in six_unions if target & ~mask == 0)


def completion_summary(values):
    return (
        len(values),
        values[:3],
        sha256(repr(values).encode("ascii")).hexdigest(),
    )


def nontrivial_sccs(vertices, adjacency):
    reverse = {vertex: set() for vertex in vertices}
    for source, targets in adjacency.items():
        for target in targets:
            reverse[target].add(source)

    def reach(start, graph):
        seen = {start}
        stack = [start]
        while stack:
            source = stack.pop()
            for target in graph[source]:
                if target not in seen:
                    seen.add(target)
                    stack.append(target)
        return seen

    remaining = set(vertices)
    result = []
    while remaining:
        vertex = min(remaining)
        component = tuple(sorted(reach(vertex, adjacency) & reach(vertex, reverse)))
        require(component, ("empty SCC", vertex))
        remaining.difference_update(component)
        if len(component) > 1:
            result.append(component)
    return tuple(sorted(result))


def longest_path(vertices, edges, potential):
    incoming = defaultdict(list)
    for source, target in edges:
        incoming[target].append(source)
    best = {vertex: (vertex,) for vertex in vertices}
    for vertex in sorted(vertices, key=lambda item: (potential(item), item)):
        candidates = tuple(best[source] + (vertex,) for source in incoming[vertex])
        if candidates:
            best[vertex] = max(candidates, key=lambda path: (len(path), path))
    return max(best.values(), key=lambda path: (len(path), path))


def k3_feature_counts(k3, row, divisor):
    support_count, _body, _ruler, arcs = row
    q = divisor // 7
    distribution = k3.expected_distribution(divisor)
    c3_distribution = k3.c3_denominator_distribution(divisor)
    histogram_q = k3.combined.residue_load_histogram(arcs, q)
    n_by_c = {
        c: sum(count for load, count in histogram_q if load > c)
        for c in range(5)
    }
    small_loads = {}
    for small_divisor in (2, 3, 4):
        if divisor % small_divisor == 0:
            histogram = k3.combined.residue_load_histogram(arcs, small_divisor)
            small_loads[small_divisor] = k3.combined.top_class_load(histogram, 1)

    patterns = {pattern for pattern, _c, _capacity in distribution}
    support_limits = {}
    for pattern in patterns:
        weights, marginals = k3.base.small_vectors(pattern, divisor, small_loads)
        support_limits[pattern] = k3.base.status_need_limit(weights, marginals)

    result = Counter()
    for (pattern, c, capacity), multiplicity in distribution.items():
        if c == 3:
            continue
        if not k3.mean_passes(n_by_c[c], c, q):
            continue
        if capacity + support_limits[pattern] < support_count:
            continue
        result[("generic", pattern, c, capacity)] += multiplicity

    if k3.mean_passes(n_by_c[3], 3, q):
        for feature, multiplicity in c3_distribution.items():
            pattern, spike_divisor, capacity, allowance = feature
            if capacity + support_limits[pattern] < support_count:
                continue
            if allowance < n_by_c[3]:
                continue
            result[("one_spike", pattern, 3, capacity, spike_divisor, allowance)] += multiplicity
    return result


def compile_k2_queries(d6, queries):
    compiler = which("g++") or which("clang++")
    require(compiler is not None, "no C++ compiler")
    lines = [str(len(queries))]
    for divisor in sorted(queries):
        lines.append(f"{divisor} {len(queries[divisor])}")
        lines.extend(f"{c} {need}" for c, need in queries[divisor])
    input_text = "\n".join(lines) + "\n"
    outputs = []
    with TemporaryDirectory(prefix="lrc14-six-divisor-lift-") as temporary:
        for optimization in ("-O2", "-O3"):
            executable = Path(temporary) / f"engine-{optimization[1:]}.exe"
            compiled = run(
                [compiler, "-std=c++17", optimization, "-DNDEBUG",
                 str(d6.base.ENGINE_PATH), "-o", str(executable)],
                stdout=PIPE,
                stderr=PIPE,
                text=True,
                check=False,
            )
            require(compiled.returncode == 0, ("compile", optimization, compiled.stderr))
            result = run(
                [str(executable)],
                input=input_text,
                stdout=PIPE,
                stderr=PIPE,
                text=True,
                check=False,
            )
            require(result.returncode == 0, ("engine", optimization, result.stderr))
            outputs.append(result.stdout)
    require(outputs[0] == outputs[1], "k2 O2/O3 mismatch")
    answers = {}
    for line in outputs[0].splitlines():
        fields = line.split()
        if fields and fields[0] == "Q":
            _tag, divisor, c, need, count, first = fields
            answers[(int(divisor), int(c), int(need))] = (int(count), int(first))
    require(len(answers) == sum(map(len, queries.values())), ("k2 answers", len(answers)))
    return answers, sha256(outputs[0].encode()).hexdigest()


def k2_target_rows(d6):
    rows = defaultdict(list)
    unavailable = []
    for body, divisor in K2_KEYS:
        ruler, ranges = d6.base.support.safe_cell_ranges(body)
        if ruler % divisor:
            unavailable.append((body, divisor, "not_divisor"))
            continue
        support_count = d6.base.support.support_size_bitset(divisor, ranges)
        if support_count * d6.base.SUPPORT_CUTOFF.denominator > divisor * d6.base.SUPPORT_CUTOFF.numerator:
            unavailable.append((body, divisor, "support_cutoff"))
            continue
        q = divisor // 7
        arcs = d6.base.combined.projected_support_arcs(divisor, ranges)
        require(sum(right - left for left, right in arcs) == support_count,
                ("k2 arcs", body, divisor))
        histogram = d6.base.combined.residue_load_histogram(arcs, q)
        needs = tuple(
            sum(count for load, count in histogram if load > c)
            for c in range(1, 6)
        )
        rows[divisor].append((body, ruler, support_count, histogram, needs))
    return dict(rows), tuple(unavailable)


def k2_row_counts(d6, rows, answers):
    counts = {}
    for divisor, records in rows.items():
        q = divisor // 7
        coefficient = d6.fixed_d6_shape_count(divisor)
        for body, _ruler, _support, _histogram, needs in records:
            total = 0
            for c, need in enumerate(needs, 1):
                count = answers[(divisor, c, need)][0]
                if c == 4 and coefficient and 1 <= needs[3] <= q // 6:
                    count -= coefficient
                require(count >= 0, ("k2 negative", body, divisor, c, count))
                total += count
            counts[(body, divisor)] = total
    return counts


def main():
    source = Path(__file__)
    tree = ast.parse(source.read_text(encoding="utf-8"))
    require(not any(isinstance(node, ast.Assert) for node in ast.walk(tree)), "assert found")
    for path, expected in EXPECTED_DEPENDENCIES.items():
        require(lf_sha256(path) == expected, ("dependency", path, lf_sha256(path), expected))

    composition = load_module("six_lift_k3_composition", K3_COMPOSITION_PATH)
    k3 = load_module("six_lift_k3_ledger", K3_LEDGER_PATH)
    complement = load_module("six_lift_complement", COMPLEMENT_PATH)
    points, atoms, six_unions = arrangement(complement)
    _solver_points, _solver_atoms, solve_small, _exact = composition.build_solver(complement)
    require((points, atoms) == (_solver_points, _solver_atoms), "arrangement mismatch")

    by_divisor, body_count, body_divisor_rows = k3.q7.build_rows()
    require((body_count, body_divisor_rows) == (3_003, 251_536), "k3 universe")
    before = composition.aggregate(k3, by_divisor)
    require(before["summary"] == composition.EXPECTED_PRE, ("k3 pre", before["summary"]))
    records = {
        (body, divisor): row
        for divisor, divisor_rows in by_divisor.items()
        for row in divisor_rows
        for _support, body, _ruler, _arcs in (row,)
        if (body, divisor) in before["rows"]
    }
    require(len(records) == before["summary"][2], ("k3 records", len(records)))
    targets = {}
    closed = set()
    for key, row in records.items():
        body, divisor = key
        _support, check_body, _ruler, arcs = row
        require(body == check_body, ("k3 body", key, row))
        target = target_for_row(composition, complement, divisor, arcs, points, atoms)
        targets[key] = target
        if solve_small(target) is not None:
            closed.add(key)
    require(len(closed) == 7, ("k3 current <=5", len(closed)))
    allowed = before["rows"] - closed
    after = composition.aggregate(k3, by_divisor, allowed)
    require(after["summary"] == EXPECTED_K3_POST, ("k3 post", after["summary"]))

    completions = {
        key: six_completions(targets[key], six_unions)
        for key in sorted(after["rows"])
    }
    self_keys = tuple(sorted(
        key for key, values in completions.items() if key[0] in values
    ))
    no_six = tuple(sorted(key for key, values in completions.items() if not values))
    nonself = tuple(sorted(
        (body, target, divisor)
        for (body, divisor), values in completions.items()
        for target in values
        if target != body
    ))
    require((len(self_keys), len(nonself), len(no_six)) == EXPECTED_K3_CLASSES,
            (len(self_keys), len(nonself), len(no_six)))
    require(all(len(completions[key]) == 1 for key in after["rows"] if completions[key]),
            "k3 exact-six multiplicity")
    require(all(divisor == records[(body, divisor)][2] for body, divisor in self_keys),
            "k3 self not full-period")

    body_edges = tuple(sorted({(source_body, target_body) for source_body, target_body, _divisor in nonself}))
    require(len(body_edges) == 20, ("k3 body edges", len(body_edges)))
    phi = lambda body: int(1 in body) + 2 * int(5 in body) - 2 * int(10 in body)
    deltas = Counter(phi(target) - phi(source_body) for source_body, target in body_edges)
    require(all(delta > 0 for delta in deltas), ("k3 potential", deltas))
    body_vertices = tuple(sorted({body for edge in body_edges for body in edge}))
    longest_body = longest_path(body_vertices, body_edges, phi)
    require(longest_body == EXPECTED_LONGEST_BODY_PATH, ("k3 longest", longest_body))

    exact_six_keys = set(self_keys) | {(source_body, divisor) for source_body, _target, divisor in nonself}
    row_adjacency = {key: set() for key in exact_six_keys}
    row_status = []
    feature_join = []
    feature_cache = {}
    for source_body, target_body, divisor in nonself:
        source_key = (source_body, divisor)
        target_key = (target_body, divisor)
        if target_key not in after["rows"]:
            status = "absent_current_key"
        elif target_key in no_six:
            status = "live_no_six"
        elif target_key in self_keys:
            status = "live_self"
        else:
            status = "live_nonself"
        if target_key in exact_six_keys:
            row_adjacency[source_key].add(target_key)
        row_status.append((source_body, target_body, divisor, status))

        source_features = feature_cache.setdefault(
            source_key, k3_feature_counts(k3, records[source_key], divisor)
        )
        target_features = Counter()
        if target_key in after["rows"]:
            target_features = feature_cache.setdefault(
                target_key, k3_feature_counts(k3, records[target_key], divisor)
            )
            require(sum(target_features.values()) == after["row_counts"][target_key],
                    ("k3 target feature count", target_key))
        require(sum(source_features.values()) == after["row_counts"][source_key],
                ("k3 source feature count", source_key))
        common = set(source_features) & set(target_features)
        require(all(source_features[item] == target_features[item] for item in common),
                ("k3 feature multiplicity", source_key, target_key))
        feature_join.append((
            source_body,
            target_body,
            divisor,
            len(source_features),
            sum(source_features.values()),
            len(common),
            sum(source_features[item] for item in common),
            status,
        ))
    require(not nontrivial_sccs(tuple(row_adjacency), row_adjacency), "k3 row SCC")

    d6 = load_module("six_lift_k2_d6", K2_D6_PATH)
    k2_rows, unavailable = k2_target_rows(d6)
    queries = d6.base.engine_queries(k2_rows)
    answers, k2_engine_hash = compile_k2_queries(d6, queries)
    require(k2_engine_hash == EXPECTED_K2_ENGINE_OUTPUT_SHA256,
            ("k2 engine hash", k2_engine_hash))
    k2_counts = k2_row_counts(d6, k2_rows, answers)
    for key in K2_KEYS:
        k2_counts.setdefault(key, 0)

    k2_completions = {}
    for body, divisor in K2_KEYS:
        ruler, ranges = d6.base.support.safe_cell_ranges(body)
        require(ruler % divisor == 0, ("k2 divisor", body, divisor))
        arcs = complement.residue_arcs(divisor, ranges)
        target = target_for_row(composition, complement, divisor, arcs, points, atoms)
        k2_completions[(body, divisor)] = six_completions(target, six_unions)

    require(C in k2_completions[(A, 17_640)], "k2 A->C")
    require(B in k2_completions[(C, 17_640)], "k2 C->B")
    require(A in k2_completions[(B, 4_410)], "k2 B->A")
    require(k2_completions[(B, 17_640)] == (B,), ("k2 B self", k2_completions[(B, 17_640)]))
    require(all(k2_counts[key] > 0 for key in ((A, 17_640), (C, 17_640), (B, 17_640), (B, 4_410))),
            ("k2 positive chain", k2_counts))
    require(k2_counts[(A, 4_410)] == 0, ("k2 A@4410 terminal", k2_counts[(A, 4_410)]))

    k2_lift_edges = []
    for source_body, target_body, divisor in (
        (A, C, 17_640),
        (C, B, 17_640),
        (B, A, 4_410),
    ):
        if k2_counts[(target_body, divisor)] > 0:
            k2_lift_edges.append(((source_body, divisor), (target_body, divisor)))
    k2_vertices = tuple(sorted(set(K2_KEYS)))
    k2_adjacency = {key: set() for key in k2_vertices}
    for source_key, target_key in k2_lift_edges:
        k2_adjacency[source_key].add(target_key)
    require(not nontrivial_sccs(k2_vertices, k2_adjacency), "k2 divisor-lift SCC")
    require({17_640} & {17_640} & {4_410} == set(), "k2 label intersection")

    feature_summary = (
        sum(row[5] > 0 for row in feature_join),
        sum(row[5] for row in feature_join),
        sum(row[6] for row in feature_join),
    )
    status_histogram = tuple(sorted(Counter(row[3] for row in row_status).items()))
    require(status_histogram == EXPECTED_K3_STATUS_HISTOGRAM,
            ("k3 target statuses", status_histogram))
    require(feature_summary == EXPECTED_K3_FEATURE_SUMMARY,
            ("k3 feature summary", feature_summary))
    require(all((row[5] > 0) == (row[7] == "live_self") for row in feature_join),
            "k3 feature transport does not match self targets")
    k2_summary = tuple(sorted(
        (key, k2_counts[key], completion_summary(k2_completions[key]))
        for key in K2_KEYS
    ))
    semantic_payload = (
        tuple(sorted((key, values) for key, values in completions.items())),
        nonself,
        tuple(sorted(deltas.items())),
        longest_body,
        tuple(row_status),
        tuple(feature_join),
        k2_summary,
        tuple(k2_lift_edges),
    )
    semantic = sha256(repr(semantic_payload).encode("ascii")).hexdigest()
    require(semantic == EXPECTED_SEMANTIC_SHA256,
            ("semantic", semantic, EXPECTED_SEMANTIC_SHA256))

    print("LRC14 exact-six divisor and feature lift")
    print("status=FINITE-EXACT unnumbered structural sidecar;not an LRC14 terminal or physical iteration theorem")
    print(f"k3_current_keys={after['summary'][2]};exact6_self={len(self_keys)};exact6_nonself_edges={len(nonself)};no_six={len(no_six)}")
    print(f"k3_nonself_edges={nonself}")
    print(f"k3_boolean_potential=1_[1inF]+2_[5inF]-2_[10inF];positive_delta_histogram={tuple(sorted(deltas.items()))}")
    print(f"k3_longest_nonself_path={longest_body};length={len(longest_body)-1}")
    print(f"k3_same_divisor_target_status_histogram={status_histogram};rows={tuple(row_status)}")
    print(f"k3_feature_join_rows={tuple(feature_join)}")
    print(f"k3_feature_join_summary=positive_edges,total_classes,total_occurrences={feature_summary}")
    print(f"k2_body_projection_cycle={(A, C, B, A)};edge_divisors={(17_640, 17_640, 4_410)};common_divisor=none")
    print(f"k2_key_counts_and_completion_summaries={k2_summary}")
    print(f"k2_same_divisor_lift_edges={tuple(k2_lift_edges)};nontrivial_sccs=()")
    print(f"k2_targeted_engine_output_sha256={k2_engine_hash};unavailable={unavailable}")
    print("projection_lemma=a projected directed cycle lifts under a conserved sidecar only if one sidecar label is compatible around the entire cycle")
    print("typed_boundary=feature classes aggregate literal denominator multisets and retain no common phase,next sector,tail labels,ancestry,owners,or physical time")
    print(f"semantic_sha256={semantic}")
    print(f"script_sha256_lf={lf_sha256(source)}")
    print("reproducibility=exact Fraction/bitsets/GF counts;targeted O2/O3 engine;no float/no assert/no network")
    print("verdict=PASS;k2_C3_is_divisor_splicing;k3_body_DAG_has_strict_Boolean_potential;feature_lift_acyclic")


if __name__ == "__main__":
    main()
