#!/usr/bin/env python3
"""Portable exact pruner/evaluator for THM-741's 21 j=4 flood bodies.

The historical THM-741 driver hard-codes a Windows scratch directory and its
resume ledger is not present in this checkout.  This companion imports only
its exact interval kernel (hash guarded), accepts portable CLI paths, and adds
one proved pruning invariant before the expensive E3 interval construction.

Let G=G(E2) have r components and measure m, let a=v3, and put
m_a=|G\D_a|.  For every b>a, THM-732's one-speed discrepancy bound gives

    |G intersect D_b| <= m/7 + sqrt(2) r/(7b).

Consequently

    |G \ (D_a union D_b)|
      >= m_a - m/7 - sqrt(2) r/(7b).                 (*).

Using S2=99/70>sqrt(2), every integer

    b > U(G,a) = S2*r / (7*(m_a-m/7))

is closed whenever the denominator is positive.  This keeps the fixed E2
carrier instead of paying for the fragmented G(E2 union {a}).  Before even
computing m_a, THM-733/P2 gives

    m_a-m/7 >= 5m/7 - 8r/(49a),

which supplies a cheap all-b pre-screen.  Remaining integer b are swept by
the same exact Fraction interval oracle as THM-741.  No Fano symmetry is used.

Examples:

  python3 04-computation/lrc14_j4_flood_portable_pruner_codex_S14.py --edge 5 7
  python3 ... --edge 5 7 --frontier-only
  python3 ... --all --workers 8 --state /tmp/thm741_flood_S14.jsonl

The body-level JSONL state is optional, append-only, and platform-neutral.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import multiprocessing as mp
import os
import time
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
THM741_PATH = ROOT / "04-computation/lrc14_thm741_2002_body_j4_tree_kps_S128c5.py"
THM741_SHA256 = "5aa81d9d78273c8f9e3e7a6574091a3bc3f64ab6086c7024c15f9420c99dac96"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_thm741():
    digest = sha256(THM741_PATH)
    if digest != THM741_SHA256:
        raise RuntimeError(f"THM-741 dependency changed: expected {THM741_SHA256}, got {digest}")
    spec = importlib.util.spec_from_file_location("thm741_portable_dependency", THM741_PATH)
    if spec is None or spec.loader is None:
        raise RuntimeError("cannot load THM-741 dependency")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


THM = load_thm741()


def flood_body(edge: tuple[int, int]) -> tuple[int, ...]:
    a, b = edge
    if not (1 <= a < b <= 7):
        raise ValueError("a flood edge must satisfy 1<=a<b<=7")
    return tuple(sorted((*range(8, 15), a, b)))


def floor_fraction(value: F) -> int:
    return value.numerator // value.denominator


def missing_obligations(speeds: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(q for q in range(2, 15) if not any(speed % q == 0 for speed in speeds))


def p2_precloses(m: F, r: int, a: int) -> bool:
    """Whether P2 plus (*) closes every b>a against the fixed E2 body."""
    denominator_lower = F(5, 7) * m - F(8 * r, 49 * a)
    return (
        denominator_lower > 0
        and denominator_lower - THM.S2 * r / (7 * (a + 1)) > 0
    )


def first_p2_preclosed(m: F, r: int, start: int, stop: int) -> int:
    """First a in [start,stop) closed by the monotone P2 screen, else stop."""
    if start >= stop or not p2_precloses(m, r, stop - 1):
        return stop
    low, high = start, stop - 1
    while low < high:
        middle = (low + high) // 2
        if p2_precloses(m, r, middle):
            high = middle
        else:
            low = middle + 1
    assert p2_precloses(m, r, low)
    assert low == start or not p2_precloses(m, r, low - 1)
    return low


COUNTER_KEYS = (
    "E1",
    "E2",
    "v3_nodes",
    "p2_preclosed_nodes",
    "exact_m3_nodes",
    "exact_m3_closed_nodes",
    "fallback_nodes",
    "candidate_nodes",
    "candidate_sweeps",
    "positive_sweeps",
    "zero_sweeps",
)


def body_work(payload: tuple[tuple[int, int], bool, int | None]) -> dict[str, object]:
    edge, frontier_only, v1_only = payload
    started = time.time()
    E = flood_body(edge)
    Eset = set(E)
    assert missing_obligations(E) == ()
    GE, rE, mE = THM.good_norm(E)
    threshold = 3 * mE / (THM.S2 * rE)
    V1 = THM.minV(4, threshold.numerator, threshold.denominator)

    counters = {key: 0 for key in COUNTER_KEYS}
    minimum_positive: F | None = None
    minimum_family: tuple[int, ...] | None = None
    failures: list[tuple[int, ...]] = []

    v1_values = range(1, V1) if v1_only is None else (v1_only,)
    for v1 in v1_values:
        if not 1 <= v1 < V1:
            raise ValueError(f"v1 branch {v1} outside [1,{V1})")
        if v1 in Eset:
            continue
        r1, m1, G1 = THM.subtract(GE, v1)
        assert m1 > 0
        counters["E1"] += 1
        threshold1 = 4 * m1 / (THM.S2 * r1)
        V2 = THM.minV(3, threshold1.numerator, threshold1.denominator)
        for v2 in range(v1 + 1, V2):
            if v2 in Eset:
                continue
            r2, m2, G2 = THM.subtract(G1, v2)
            assert m2 > 0
            counters["E2"] += 1
            threshold2 = 5 * m2 / (THM.S2 * r2)
            V3 = THM.minV(2, threshold2.numerator, threshold2.denominator)
            v3_start = v2 + 1
            preclosed_start = first_p2_preclosed(m2, r2, v3_start, V3)
            total_v3 = max(0, V3 - v3_start) - sum(
                v3_start <= speed < V3 for speed in Eset
            )
            bulk_preclosed = max(0, V3 - preclosed_start) - sum(
                preclosed_start <= speed < V3 for speed in Eset
            )
            counters["v3_nodes"] += total_v3
            counters["p2_preclosed_nodes"] += bulk_preclosed
            for v3 in range(v3_start, preclosed_start):
                if v3 in Eset:
                    continue
                assert not p2_precloses(m2, r2, v3)

                # Exact measure only; no E3 interval list yet.
                m3 = THM.subtract_sparse(G2, v3)
                assert m3 > 0
                counters["exact_m3_nodes"] += 1
                denominator = m3 - m2 / 7
                if denominator > 0:
                    upper = THM.S2 * r2 / (7 * denominator)
                    bmax = floor_fraction(upper)
                    if bmax <= v3:
                        counters["exact_m3_closed_nodes"] += 1
                        continue
                else:
                    # The fixed-E2 inequality is silent.  Fall back to the
                    # original exact E3 one-peel tail; this branch is retained
                    # for generality even though the audited flood run need not
                    # use it.
                    r3, exact_m3, G3 = THM.subtract(G2, v3)
                    assert exact_m3 == m3
                    upper = THM.S2 * r3 / (6 * m3)
                    bmax = floor_fraction(upper)
                    counters["fallback_nodes"] += 1

                candidates = [
                    v4 for v4 in range(v3 + 1, bmax + 1) if v4 not in Eset
                ]
                if not candidates:
                    counters["exact_m3_closed_nodes"] += 1
                    continue
                counters["candidate_nodes"] += 1
                counters["candidate_sweeps"] += len(candidates)
                if frontier_only:
                    continue

                if denominator > 0:
                    r3, exact_m3, G3 = THM.subtract(G2, v3)
                    assert exact_m3 == m3 and r3 == len(G3)
                for v4 in candidates:
                    lonely_measure = THM.subtract_sparse(G3, v4)
                    family = tuple(sorted(E + (v1, v2, v3, v4)))
                    if lonely_measure > 0:
                        counters["positive_sweeps"] += 1
                        if minimum_positive is None or lonely_measure < minimum_positive:
                            minimum_positive = lonely_measure
                            minimum_family = family
                    else:
                        counters["zero_sweeps"] += 1
                        if not missing_obligations(family):
                            failures.append(family)

    result: dict[str, object] = {
        "edge": list(edge),
        "E": list(E),
        "r": rE,
        "m": str(mE),
        "V1": V1,
        **counters,
        "minimum_positive": None if minimum_positive is None else str(minimum_positive),
        "minimum_family": None if minimum_family is None else list(minimum_family),
        "failures": [list(family) for family in failures],
        "frontier_only": frontier_only,
        "v1_only": v1_only,
        "seconds": round(time.time() - started, 3),
    }
    return result


def aggregate_branches(
    edge: tuple[int, int], frontier_only: bool, branches: list[dict[str, object]], wall: float
) -> dict[str, object]:
    if not branches:
        raise RuntimeError("no v1 branches to aggregate")
    first = branches[0]
    result: dict[str, object] = {
        "edge": list(edge),
        "E": first["E"],
        "r": first["r"],
        "m": first["m"],
        "V1": first["V1"],
        **{key: sum(int(branch[key]) for branch in branches) for key in COUNTER_KEYS},
        "failures": [failure for branch in branches for failure in branch["failures"]],
        "frontier_only": frontier_only,
        "v1_only": None,
        "seconds": round(sum(float(branch["seconds"]) for branch in branches), 3),
        "wall_seconds": round(wall, 3),
    }
    positive = [
        (F(branch["minimum_positive"]), branch["minimum_family"])
        for branch in branches
        if branch["minimum_positive"] is not None
    ]
    if positive:
        minimum, family = min(positive)
        result["minimum_positive"] = str(minimum)
        result["minimum_family"] = family
    else:
        result["minimum_positive"] = None
        result["minimum_family"] = None
    return result


def render(result: dict[str, object]) -> str:
    lines = [
        f"edge={tuple(result['edge'])} E={tuple(result['E'])}",
        f"root r={result['r']} m={result['m']} V1={result['V1']}",
        (
            "nodes E1/E2/v3="
            f"{result['E1']}/{result['E2']}/{result['v3_nodes']} ; "
            "P2-preclosed/exact-m3/exact-m3-closed/fallback="
            f"{result['p2_preclosed_nodes']}/{result['exact_m3_nodes']}/"
            f"{result['exact_m3_closed_nodes']}/{result['fallback_nodes']}"
        ),
        (
            "bottom candidate nodes/sweeps="
            f"{result['candidate_nodes']}/{result['candidate_sweeps']} ; "
            f"positive/zero={result['positive_sweeps']}/{result['zero_sweeps']}"
        ),
        (
            f"minimum positive={result['minimum_positive']} at {result['minimum_family']} ; "
            f"covering failures={len(result['failures'])}"
        ),
        f"frontier_only={result['frontier_only']}",
    ]
    return "\n".join(lines)


def read_done(state: Path | None) -> dict[tuple[int, int], dict[str, object]]:
    done: dict[tuple[int, int], dict[str, object]] = {}
    if state is None or not state.exists():
        return done
    for line in state.read_text().splitlines():
        try:
            row = json.loads(line)
            done[tuple(row["edge"])] = row
        except (json.JSONDecodeError, KeyError, TypeError):
            continue
    return done


def main() -> None:
    if not __debug__:
        raise RuntimeError("exact flood audit requires assertions; do not run with python -O")
    parser = argparse.ArgumentParser()
    choice = parser.add_mutually_exclusive_group(required=True)
    choice.add_argument("--edge", nargs=2, type=int, metavar=("A", "B"))
    choice.add_argument("--all", action="store_true")
    parser.add_argument("--frontier-only", action="store_true")
    parser.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 2) - 2))
    parser.add_argument("--state", type=Path)
    args = parser.parse_args()

    edges = (
        [tuple(args.edge)]
        if args.edge is not None
        else list(combinations(range(1, 8), 2))
    )
    edges = [tuple(sorted(edge)) for edge in edges]
    done = read_done(args.state)
    pending = [edge for edge in edges if edge not in done]

    print("THM-741 J=4 FLOOD PORTABLE FIXED-E2 PRUNER")
    print("=" * 86)
    print(f"dependency_sha256={THM741_SHA256}")
    print("floods are K7 edges only as an index; no Fano symmetry quotient is used")
    print("new invariant: L(E2+a+b)>=m(E2+a)-m(E2)/7-S2*r(E2)/(7b)")
    print(
        f"requested={len(edges)} resumed={len(done)} pending={len(pending)} "
        f"frontier_only={args.frontier_only} workers={args.workers}"
    )
    results = list(done.values())
    pool = None
    if len(pending) == 1 and args.workers > 1:
        edge = pending[0]
        E = flood_body(edge)
        _, rE, mE = THM.good_norm(E)
        threshold = 3 * mE / (THM.S2 * rE)
        V1 = THM.minV(4, threshold.numerator, threshold.denominator)
        branch_payloads = [
            (edge, args.frontier_only, v1)
            # Large v1 branches are empirically the expensive ones.  Submit
            # them first so the portable intra-body pool does not finish with
            # one long straggler; this changes scheduling only, never logic.
            for v1 in range(V1 - 1, 0, -1)
            if v1 not in set(E)
        ]
        branch_started = time.time()
        with mp.Pool(min(args.workers, len(branch_payloads))) as branch_pool:
            branches = []
            for index, branch in enumerate(
                branch_pool.imap_unordered(body_work, branch_payloads, chunksize=1), start=1
            ):
                branches.append(branch)
                if index % 20 == 0 or index == len(branch_payloads):
                    print(f"v1 branches complete={index}/{len(branch_payloads)}", flush=True)
        iterator = iter(
            [
                aggregate_branches(
                    edge, args.frontier_only, branches, time.time() - branch_started
                )
            ]
        )
    else:
        payloads = [(edge, args.frontier_only, None) for edge in pending]
        if len(payloads) <= 1 or args.workers == 1:
            iterator = map(body_work, payloads)
        else:
            pool = mp.Pool(min(args.workers, len(payloads)))
            iterator = pool.imap_unordered(body_work, payloads, chunksize=1)
    try:
        for result in iterator:
            results.append(result)
            print(render(result), flush=True)
            if args.state is not None:
                args.state.parent.mkdir(parents=True, exist_ok=True)
                with args.state.open("a", encoding="utf-8") as state_file:
                    state_file.write(json.dumps(result, sort_keys=True) + "\n")
    finally:
        if pool is not None:
            pool.close()
            pool.join()

    selected = [row for row in results if tuple(row["edge"]) in set(edges)]
    print("-" * 86)
    print(
        "aggregate completed=%d/21 v3=%d preclosed=%d exact_m3=%d candidate_sweeps=%d failures=%d"
        % (
            len(selected),
            sum(int(row["v3_nodes"]) for row in selected),
            sum(int(row["p2_preclosed_nodes"]) for row in selected),
            sum(int(row["exact_m3_nodes"]) for row in selected),
            sum(int(row["candidate_sweeps"]) for row in selected),
            sum(len(row["failures"]) for row in selected),
        )
    )
    if not args.frontier_only and len(selected) == len(edges) and not any(row["failures"] for row in selected):
        print(f"EXACTLY CLOSED REQUESTED FLOOD BODIES={len(edges)}")
    elif args.frontier_only:
        print("PORTABLE FRONTIER ONLY: exact bottom sweeps deliberately not evaluated")
    else:
        print("REQUEST NOT CLOSED: incomplete rows or covering zero-measure failures remain")
    print(f"source_sha256={sha256(Path(__file__).resolve())}")
    print("ALL ASSERTIONS PASSED")


if __name__ == "__main__":
    main()
