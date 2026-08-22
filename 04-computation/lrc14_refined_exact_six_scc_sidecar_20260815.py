#!/usr/bin/env python3
"""Targeted exact SCC audit for the refined six-clock mutation relation.

The frozen raw exact-six relation has only two nontrivial SCCs, both directed
three-cycles.  A refined current-row relation is obtained only by deleting raw
body/divisor rows, so it cannot create or merge SCCs.  This sidecar therefore
reconstructs every raw nonself edge internal to those two SCCs, then evaluates
the current k=2 and k=3 occurrence ledgers on exactly the six source rows.

No global refined ledger is replayed.  The k=2 finite-state engine receives
only the targeted queries and is compiled under both -O2 and -O3.  The k=3
composition evaluates only the three targeted divisors.  Strict endpoints and
open atoms are retained when reconstructing the raw edge rows.

This is a FINITE-EXACT structural sidecar, not an LRC(14) terminal and not a
physical iteration theorem.
"""

from __future__ import annotations

import ast
from collections import defaultdict
from fractions import Fraction as Q
from functools import reduce
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations
from pathlib import Path
from shutil import which
from subprocess import PIPE, run
from tempfile import TemporaryDirectory


ROOT = Path(__file__).resolve().parents[1]
RAW_SCRIPT = ROOT / "04-computation/lrc14_exact_six_complement_mutation_graph_probe_20260815.py"
RAW_OUTPUT = ROOT / "05-knowledge/results/lrc14_exact_six_complement_mutation_graph_probe_20260815.out"
BASE_PATH = ROOT / "04-computation/lrc14_k1_body_complement_clock_scan_kps_s172.py"
SUPPORT_PATH = ROOT / "04-computation/lrc14_two_drift_body_projection_support_thm2928.py"
K2_PATH = ROOT / "04-computation/lrc14_k2_refined_complement_clock_composition_kps_s175.py"
D6_PATH = ROOT / "04-computation/lrc14_k2_d6_located_phase_closure_thm2928.py"
K2_BASE_PATH = ROOT / "04-computation/lrc14_k2_septimal_floor_exception_gf_thm2928.py"
K2_ENGINE_PATH = ROOT / "04-computation/lrc14_k2_septimal_floor_exception_engine_thm2928.cpp"
K3_PATH = ROOT / "04-computation/lrc14_k3_refined_complement_clock_composition_kps_s176.py"
K3_LEDGER_PATH = ROOT / "04-computation/lrc14_k3_four_drift_expected_spike_gf_thm2928.py"
K3_Q7_PATH = ROOT / "04-computation/lrc14_k3_four_drift_q7_all_D_gf_thm2928.py"

EXPECTED_DEPENDENCIES = {
    RAW_SCRIPT: "a799d77af7d930e6a46ab7f22544d78feda327a2a9221025a2c28cf9ed5c85ea",
    RAW_OUTPUT: "3627ec4fc55ecbca4a571d76dac1e87c9c9f20c4ff6ba9d1346642ccae5117c4",
    BASE_PATH: "bdb2001cf22f7e92884e895b0095021e42e8f1febd9adbf779b250a2f6c53507",
    SUPPORT_PATH: "778842c0e8e7172835ca6ae673fb6156f212d4296e672bce4e7cc2815195bf1a",
    K2_PATH: "414b3777cc44f2e059d5bc4258ab555fc055de8c0f8b3ebfd99ce0c90c7da14c",
    D6_PATH: "9f300459b273ad1825d3fe3e9274c6afe609f2d581e9df3d2be1780d347e541b",
    K2_BASE_PATH: "085b4e2747a48bdbc1125e894af7d4f647dfdd7be86a00cf02dea2a8667e26dc",
    K2_ENGINE_PATH: "664d0df36d104d959279605c8ea8539d61ab595b155e5157fa7d0433f1b7944c",
    K3_PATH: "27e4bff52705189bf8ff73db42d76d4e2fc94c44330d295f166f8d4217cb1804",
    K3_LEDGER_PATH: "05e365a654b32e66b814dcbce9385a2d13c22a2c84a5474e0855dcab6262b055",
    K3_Q7_PATH: "3cc07195d580c5c5c01457ea95b58837a25c2176d326d12feaccc8e0bfa28dcc",
}

POOL = tuple(range(1, 15))
CYCLE_A = (
    (1, 3, 7, 8, 9, 10),
    (2, 4, 5, 7, 9, 12),
    (4, 5, 7, 8, 9, 10),
)
CYCLE_B = (
    (2, 3, 5, 7, 9, 12),
    (3, 7, 8, 9, 10, 12),
    (4, 5, 6, 7, 9, 10),
)
EXPECTED_RAW_SCCS = (CYCLE_A, CYCLE_B)
EXPECTED_RAW_EDGE_ROWS = (
    ((1, 3, 7, 8, 9, 10), (4, 5, 7, 8, 9, 10), 35_280, 17_640, 9_984),
    ((2, 3, 5, 7, 9, 12), (3, 7, 8, 9, 10, 12), 17_640, 4_410, 3_280),
    ((2, 4, 5, 7, 9, 12), (1, 3, 7, 8, 9, 10), 17_640, 4_410, 3_000),
    ((3, 7, 8, 9, 10, 12), (4, 5, 6, 7, 9, 10), 35_280, 17_640, 9_984),
    ((4, 5, 6, 7, 9, 10), (2, 3, 5, 7, 9, 12), 17_640, 8_820, 4_848),
    ((4, 5, 7, 8, 9, 10), (2, 4, 5, 7, 9, 12), 35_280, 17_640, 10_116),
)
EXPECTED_K2_WEIGHTS = (
    ((1, 3, 7, 8, 9, 10), (4, 5, 7, 8, 9, 10), 17_640, 58_269),
    ((2, 3, 5, 7, 9, 12), (3, 7, 8, 9, 10, 12), 4_410, 0),
    ((2, 4, 5, 7, 9, 12), (1, 3, 7, 8, 9, 10), 4_410, 3_199),
    ((3, 7, 8, 9, 10, 12), (4, 5, 6, 7, 9, 10), 17_640, 58_269),
    ((4, 5, 6, 7, 9, 10), (2, 3, 5, 7, 9, 12), 8_820, 145_541),
    ((4, 5, 7, 8, 9, 10), (2, 4, 5, 7, 9, 12), 17_640, 522_861),
)
EXPECTED_K3_WEIGHTS = tuple((source, target, divisor, 0) for source, target, _ruler, divisor, _support in EXPECTED_RAW_EDGE_ROWS)
EXPECTED_K2_ENGINE_OUTPUT_SHA256 = "859d75a7b7870c777438d6abee483c3458f761aca1e1e9a54179574385a399ef"
EXPECTED_SEMANTIC_SHA256 = "d3be3507d47946f2fa688f70faa572aad1312b8e2ecee30330af9e684bbbed31"


def require(condition: bool, detail: object) -> None:
    if not condition:
        raise RuntimeError(detail)


def lf_sha256(path: Path) -> str:
    return sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def load_module(name: str, path: Path):
    spec = spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, ("module", path))
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def frozen_raw_sccs() -> tuple:
    text = RAW_OUTPUT.read_text(encoding="utf-8")
    observed = []
    for sector in (2, 3):
        prefix = f"k{sector}_nontrivial_sccs="
        lines = [line for line in text.splitlines() if line.startswith(prefix)]
        require(len(lines) == 1, ("raw SCC line", sector, lines))
        observed.append(ast.literal_eval(lines[0][len(prefix):]))
    require(tuple(observed) == (EXPECTED_RAW_SCCS,) * 2, ("raw SCCs", observed))
    return tuple(observed)


def arrangement(base):
    points = base.arrangement_points(POOL)
    atoms = tuple(zip(points, points[1:]))
    samples = points + tuple((left + right) / 2 for left, right in atoms)
    masks = {
        clock: sum(
            (1 << index for index, value in enumerate(samples) if base.danger(clock, value)),
            0,
        )
        for clock in POOL
    }
    return points, atoms, masks


def target_mask(base, divisor, gaps, points, atoms):
    target = 0
    for index, value in enumerate(points):
        if (divisor * value).denominator != 1 and any(
            left < value < right for left, right in gaps
        ):
            target |= 1 << index
    offset = len(points)
    for index, (left, right) in enumerate(atoms):
        if any(
            max(left, gap_left) < min(right, gap_right)
            for gap_left, gap_right in gaps
        ):
            target |= 1 << (offset + index)
    return target


def raw_internal_edge_rows(base):
    support = base.load_support()
    points, atoms, masks = arrangement(base)
    cycle_by_body = {
        body: cycle for cycle in EXPECTED_RAW_SCCS for body in cycle
    }
    cycle_masks = {
        body: reduce(int.__or__, (masks[clock] for clock in body), 0)
        for body in cycle_by_body
    }
    small_unions = tuple(
        reduce(int.__or__, (masks[clock] for clock in clocks), 0)
        for size in range(6)
        for clocks in combinations(POOL, size)
    )
    cutoffs = {2: Q(887, 990), 3: Q(125, 143)}
    observed = {2: [], 3: []}

    for body in sorted(cycle_by_body):
        ruler, ranges = support.safe_cell_ranges(body)
        for divisor in support.divisors(ruler):
            support_count = support.support_size_bitset(divisor, ranges)
            possible_sectors = tuple(
                sector
                for sector, cutoff in cutoffs.items()
                if Q(support_count, divisor) <= cutoff
            )
            if not possible_sectors:
                continue
            arcs = base.residue_arcs(divisor, ranges)
            gaps = base.unsupported_gaps(divisor, arcs)
            target = target_mask(base, divisor, gaps, points, atoms)
            targets = tuple(
                completion
                for completion in cycle_by_body[body]
                if completion != body and target & ~cycle_masks[completion] == 0
            )
            if not targets:
                continue
            if any(target & ~union == 0 for union in small_unions):
                continue
            for sector in possible_sectors:
                observed[sector].extend(
                    (body, completion, ruler, divisor, support_count)
                    for completion in targets
                )

    result = {
        sector: tuple(sorted(rows)) for sector, rows in observed.items()
    }
    require(result == {2: EXPECTED_RAW_EDGE_ROWS, 3: EXPECTED_RAW_EDGE_ROWS}, result)
    return result


def targeted_k2_rows(d6):
    keys = {(source, divisor) for source, _target, _ruler, divisor, _support in EXPECTED_RAW_EDGE_ROWS}
    rows = defaultdict(list)
    for body, divisor in sorted(keys, key=lambda key: (key[1], key[0])):
        ruler, ranges = d6.base.support.safe_cell_ranges(body)
        require(ruler % divisor == 0, ("k2 divisor", body, ruler, divisor))
        support_count = d6.base.support.support_size_bitset(divisor, ranges)
        require(Q(support_count, divisor) <= d6.base.SUPPORT_CUTOFF, ("k2 cutoff", body, divisor))
        q = divisor // 7
        arcs = d6.base.combined.projected_support_arcs(divisor, ranges)
        require(sum(right - left for left, right in arcs) == support_count, ("k2 support", body, divisor))
        histogram = d6.base.combined.residue_load_histogram(arcs, q)
        needs = tuple(
            sum(count for load, count in histogram if load > c)
            for c in range(1, 6)
        )
        rows[divisor].append((body, ruler, support_count, histogram, needs))
    return dict(rows)


def targeted_engine(base, queries):
    compiler = which("g++") or which("clang++")
    require(compiler is not None, "no C++ compiler")
    lines = [str(len(queries))]
    for divisor in sorted(queries):
        lines.append(f"{divisor} {len(queries[divisor])}")
        lines.extend(f"{c} {need}" for c, need in queries[divisor])
    input_text = "\n".join(lines) + "\n"
    outputs = []
    with TemporaryDirectory(prefix="lrc14-refined-six-scc-") as temporary:
        for optimization in ("-O2", "-O3"):
            executable = Path(temporary) / f"engine-{optimization[1:]}"
            compiled = run(
                [compiler, "-std=c++17", optimization, "-DNDEBUG", str(base.ENGINE_PATH), "-o", str(executable)],
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
    require(outputs[0] == outputs[1], "targeted O2/O3 mismatch")
    output_hash = sha256(outputs[0].encode()).hexdigest()
    require(output_hash == EXPECTED_K2_ENGINE_OUTPUT_SHA256, ("engine hash", output_hash))
    answers = {}
    for line in outputs[0].splitlines():
        fields = line.split()
        if fields and fields[0] == "Q":
            _tag, divisor, c, need, count, first = fields
            answers[(int(divisor), int(c), int(need))] = (int(count), int(first))
    require(len(answers) == sum(map(len, queries.values())), ("answers", len(answers)))
    return answers, output_hash


def targeted_k3_rows(k3):
    keys = {(source, divisor) for source, _target, _ruler, divisor, _support in EXPECTED_RAW_EDGE_ROWS}
    rows = defaultdict(list)
    for body, divisor in sorted(keys, key=lambda key: (key[1], key[0])):
        ruler, ranges = k3.q7.support.safe_cell_ranges(body)
        require(ruler % divisor == 0, ("k3 divisor", body, ruler, divisor))
        support_count = k3.q7.support.support_size_bitset(divisor, ranges)
        require(Q(support_count, divisor) <= k3.q7.base.SUPPORT_CUTOFF, ("k3 cutoff", body, divisor))
        arcs = k3.q7.combined.projected_support_arcs(divisor, ranges)
        require(sum(right - left for left, right in arcs) == support_count, ("k3 support", body, divisor))
        rows[divisor].append((support_count, body, ruler, arcs))
    for records in rows.values():
        records.sort(key=lambda row: (row[0], row[1], row[2]))
    return dict(rows)


def weighted_edges(row_counts):
    return tuple(
        (source, target, divisor, row_counts.get((source, divisor), 0))
        for source, target, _ruler, divisor, _support in EXPECTED_RAW_EDGE_ROWS
    )


def nontrivial_sccs(weighted):
    vertices = tuple(sorted({body for cycle in EXPECTED_RAW_SCCS for body in cycle}))
    adjacency = {vertex: set() for vertex in vertices}
    reverse = {vertex: set() for vertex in vertices}
    for source, target, _divisor, weight in weighted:
        if weight:
            adjacency[source].add(target)
            reverse[target].add(source)

    def reach(start, graph):
        seen = {start}
        stack = [start]
        while stack:
            vertex = stack.pop()
            for target in graph[vertex]:
                if target not in seen:
                    seen.add(target)
                    stack.append(target)
        return seen

    remaining = set(vertices)
    components = []
    while remaining:
        vertex = min(remaining)
        component = tuple(sorted(reach(vertex, adjacency) & reach(vertex, reverse)))
        require(component, ("empty component", vertex))
        remaining.difference_update(component)
        if len(component) > 1:
            components.append(component)
    return tuple(sorted(components))


def main() -> None:
    source = Path(__file__)
    tree = ast.parse(source.read_text(encoding="utf-8"))
    require(not any(isinstance(node, ast.Assert) for node in ast.walk(tree)), "assert found")
    for path, expected in EXPECTED_DEPENDENCIES.items():
        require(lf_sha256(path) == expected, ("dependency", path, lf_sha256(path), expected))

    raw_sccs = frozen_raw_sccs()
    base = load_module("refined_scc_base", BASE_PATH)
    raw_rows = raw_internal_edge_rows(base)

    c2 = load_module("refined_scc_c2", K2_PATH)
    d6 = load_module("refined_scc_d6", D6_PATH)
    require(c2.D6_PATH.resolve() == D6_PATH.resolve(), ("k2 current input", c2.D6_PATH))
    require(d6.base.ENGINE_PATH.resolve() == K2_ENGINE_PATH.resolve(), ("k2 engine", d6.base.ENGINE_PATH))
    require(d6.base.SUPPORT_CUTOFF == Q(887, 990), ("k2 cutoff", d6.base.SUPPORT_CUTOFF))
    rows2 = targeted_k2_rows(d6)
    queries2 = d6.base.engine_queries(rows2)
    answers2, engine_hash = targeted_engine(d6.base, queries2)
    current2 = c2.aggregate(d6, rows2, answers2)
    weights2 = weighted_edges(current2["row_counts"])
    require(weights2 == EXPECTED_K2_WEIGHTS, ("k2 weights", weights2))

    c3 = load_module("refined_scc_c3", K3_PATH)
    k3 = load_module("refined_scc_k3", K3_LEDGER_PATH)
    require(c3.K3_PATH.resolve() == K3_LEDGER_PATH.resolve(), ("k3 current input", c3.K3_PATH))
    require(k3.Q7_PATH.resolve() == K3_Q7_PATH.resolve(), ("k3 q7 input", k3.Q7_PATH))
    require(k3.q7.base.SUPPORT_CUTOFF == Q(125, 143), ("k3 cutoff", k3.q7.base.SUPPORT_CUTOFF))
    rows3 = targeted_k3_rows(k3)
    current3 = c3.aggregate(k3, rows3)
    weights3 = weighted_edges(current3["row_counts"])
    require(weights3 == EXPECTED_K3_WEIGHTS, ("k3 weights", weights3))

    sccs2 = nontrivial_sccs(weights2)
    sccs3 = nontrivial_sccs(weights3)
    require(sccs2 == (CYCLE_A,), ("k2 SCCs", sccs2))
    require(not sccs3, ("k3 SCCs", sccs3))

    semantic = sha256(
        repr((raw_sccs, raw_rows, weights2, weights3, sccs2, sccs3)).encode("ascii")
    ).hexdigest()
    require(semantic == EXPECTED_SEMANTIC_SHA256, ("semantic", semantic))

    print("LRC14 refined exact-six SCC targeted sidecar")
    print("status=FINITE-EXACT unnumbered structural sidecar;not an LRC14 terminal or physical iteration theorem")
    print("subgraph_gate=refined relation deletes raw body-divisor rows;it cannot create or merge raw SCCs")
    print(f"raw_nontrivial_sccs_k2_k3={raw_sccs}")
    print(f"raw_internal_nonself_edge_rows_k2_k3={raw_rows[2]}")
    print(f"k2_refined_internal_edge_weights={weights2}")
    print(f"k2_refined_nontrivial_sccs={sccs2}")
    print(f"k3_refined_internal_edge_weights={weights3}")
    print("k3_refined_cycle_source_rows_positive=0;k3_refined_nontrivial_sccs=()")
    print("scope=all possible refined nontrivial SCCs exhausted through the two frozen raw SCCs")
    print(f"targeted_k2_engine_output_sha256={engine_hash}")
    print(f"semantic_sha256={semantic}")
    print("dependency_lf_sha256=" + repr(tuple((path.relative_to(ROOT).as_posix(), digest) for path, digest in EXPECTED_DEPENDENCIES.items())))
    print(f"script_sha256_lf={lf_sha256(source)}")
    print("reproducibility=exact Fraction/bitsets/current-row ledgers;targeted O2/O3 engine;no assert/no float/no network")
    print("verdict=PASS;k2_body_quotient_exactly_one_three_cycle;k3_acyclic")


if __name__ == "__main__":
    main()
