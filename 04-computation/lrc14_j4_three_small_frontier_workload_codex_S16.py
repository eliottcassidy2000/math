#!/usr/bin/env python3
r"""Exact positivity ledger for the residual three-small THM-741 tail.

The proved four-small shadow leaves only completed families with at most three
labels in {1,...,7}.  Thirteen of the 35 three-small bases already contain one
of the exact anchor edges 56,57,67.  This script works on the other 22 bases
K.  Starting from P={8,...,14} union K, it applies the proved THM-735/732
threshold ladder for three ordered external speeds 15<=a<b<c:

  all a,b,c >= Va:       3/Va < 4m(P)/(S2 r(P));
  after exact a, b,c>=Vb: 2/Vb < 5m(P+a)/(S2 r(P+a));
  after exact a,b, c>cap: cap=S2 r(P+a+b)/(6m(P+a+b)).

Only triples below these strict caps remain as exact interval obligations.
They are sharded by completed-family base K and evaluated with the exact
sparse interval kernel.  Five flat-index quantiles per K are independently
replayed with full interval subtraction and serialized into a manifest.  A
zero at any level is an actual covering alarm and is recorded, never pruned.

Tournament Analysis uses the 22 completed-family small bases as vertices.
The pair observable compares exact finite triple counts, then b-node counts,
then a-node counts; the gauge points toward smaller work, with lexicographic K
as the tie path.  This preserves proof-job size but destroys literal interval
geometry and margins.  Runners, Fano flags, and root edges are worse vertices
here because they do not preserve the completed-family base being decided.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import multiprocessing as mp
import os
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CORE_PATH = ROOT / "04-computation/lrc14_thm741_2002_body_j4_tree_kps_S128c5.py"
CORE_SHA256 = "5aa81d9d78273c8f9e3e7a6574091a3bc3f64ab6086c7024c15f9420c99dac96"
FOUR_SMALL_OUTPUT = ROOT / "05-knowledge/results/lrc14_j4_next_anchor_shadow_frontier_codex_S16.out"
FOUR_SMALL_SHA256 = "2b18b991d47c221b39d50966e5a0207a81b31ef7884755f702672d980137a555"
H = frozenset(range(8, 15))
ANCHORS = (frozenset((5, 6)), frozenset((5, 7)), frozenset((6, 7)))
CHECKPOINT_SCHEMA = "lrc14-three-small-exact-v1"
DEFAULT_CHECKPOINT = Path(os.environ.get("TMPDIR", "/tmp")) / (
    "lrc14_j4_three_small_frontier_workload_codex_S16.jsonl"
)
_WORKER_CORE = None


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_core():
    require(sha256(CORE_PATH) == CORE_SHA256, "THM-741 core hash changed")
    require(sha256(FOUR_SMALL_OUTPUT) == FOUR_SMALL_SHA256, "four-small theorem output changed")
    spec = importlib.util.spec_from_file_location("thm741_three_small_dependency", CORE_PATH)
    require(spec is not None and spec.loader is not None, "cannot load THM-741 core")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def floor_fraction(value: F) -> int:
    return value.numerator // value.denominator


def init_worker() -> None:
    global _WORKER_CORE
    _WORKER_CORE = load_core()


def bank_work_worker(K: tuple[int, int, int]) -> dict[str, object]:
    require(_WORKER_CORE is not None, "worker core not initialized")
    return bank_work(_WORKER_CORE, K)


def bank_work(core, K: tuple[int, int, int]) -> dict[str, object]:
    base = tuple(sorted(H | frozenset(K)))
    require(len(base) == 10, f"bad three-small base K={K}")
    good, r, m = core.good_norm(base)
    require(m > 0, f"covering base K={K}")
    Va = core.minV(3, *(4 * m / (core.S2 * r)).as_integer_ratio())

    a_nodes = b_nodes = triple_candidates = 0
    zero_after_a = zero_after_b = 0
    max_Vb = max_cap = 0
    min_after_a: tuple[F, int] | None = None
    min_after_b: tuple[F, int, int] | None = None
    for a in range(15, Va):
        r1, m1, good1 = core.subtract(good, a)
        a_nodes += 1
        if m1 == 0:
            zero_after_a += 1
            continue
        candidate1 = (m1, a)
        if min_after_a is None or candidate1 < min_after_a:
            min_after_a = candidate1
        Vb = core.minV(2, *(5 * m1 / (core.S2 * r1)).as_integer_ratio())
        max_Vb = max(max_Vb, Vb)
        for b in range(a + 1, Vb):
            r2, m2, _ = core.subtract(good1, b)
            b_nodes += 1
            if m2 == 0:
                zero_after_b += 1
                continue
            candidate2 = (m2, a, b)
            if min_after_b is None or candidate2 < min_after_b:
                min_after_b = candidate2
            cap = floor_fraction(core.S2 * r2 / (6 * m2))
            max_cap = max(max_cap, cap)
            triple_candidates += max(0, cap - b)

    require(zero_after_a == zero_after_b == 0, f"intermediate covering alarm K={K}")
    require(triple_candidates > 0, f"empty exact triple bank K={K}")
    sample_indices = {
        0,
        (triple_candidates - 1) // 4,
        (triple_candidates - 1) // 2,
        3 * (triple_candidates - 1) // 4,
        triple_candidates - 1,
    }
    require(len(sample_indices) == 5, f"collapsed cross-check sample K={K}")

    positive = zero = flat_index = 0
    minimum: tuple[F, int, int, int] | None = None
    failures = []
    sample_records = []
    for a in range(15, Va):
        _, m1, good1 = core.subtract(good, a)
        require(m1 > 0, f"one-external replay alarm K={K}, a={a}")
        r1 = len(good1)
        Vb = core.minV(2, *(5 * m1 / (core.S2 * r1)).as_integer_ratio())
        for b in range(a + 1, Vb):
            _, m2, good2 = core.subtract(good1, b)
            require(m2 > 0, f"two-external replay alarm K={K}, a={a}, b={b}")
            r2 = len(good2)
            cap = floor_fraction(core.S2 * r2 / (6 * m2))
            for c in range(b + 1, cap + 1):
                value = core.subtract_sparse(good2, c)
                if value > 0:
                    positive += 1
                    candidate = (value, a, b, c)
                    if minimum is None or candidate < minimum:
                        minimum = candidate
                else:
                    zero += 1
                    failures.append(tuple(sorted(H | frozenset(K) | {a, b, c})))
                if flat_index in sample_indices:
                    full_r, full_value, _ = core.subtract(good2, c)
                    require(
                        full_value == value,
                        f"sparse/full mismatch K={K}, a={a}, b={b}, c={c}",
                    )
                    sample_records.append(
                        f"K={','.join(map(str, K))};index={flat_index}/{triple_candidates};"
                        f"a={a};b={b};c={c};value={value};full_r={full_r}"
                    )
                flat_index += 1
    require(flat_index == triple_candidates, f"candidate replay count changed K={K}")
    require(positive + zero == triple_candidates, f"bad exact sign ledger K={K}")
    require(len(sample_records) == 5, f"bad cross-check record count K={K}")

    return {
        "K": K,
        "r": r,
        "m": m,
        "Va": Va,
        "a_nodes": a_nodes,
        "b_nodes": b_nodes,
        "triple_candidates": triple_candidates,
        "zero_after_a": zero_after_a,
        "zero_after_b": zero_after_b,
        "max_Vb": max_Vb,
        "max_cap": max_cap,
        "min_after_a": min_after_a,
        "min_after_b": min_after_b,
        "positive": positive,
        "zero": zero,
        "minimum": minimum,
        "failures": tuple(failures),
        "sample_records": tuple(sample_records),
    }


def canonical_json(payload: object) -> str:
    return json.dumps(payload, sort_keys=True, separators=(",", ":"))


def fraction_payload(value: F) -> list[int]:
    return [value.numerator, value.denominator]


def fraction_from_payload(payload: object, label: str) -> F:
    require(isinstance(payload, list) and len(payload) == 2, f"bad fraction payload for {label}")
    require(all(isinstance(item, int) for item in payload), f"noninteger fraction payload for {label}")
    numerator, denominator = payload
    require(denominator > 0, f"nonpositive fraction denominator for {label}")
    return F(numerator, denominator)


def optional_minimum_payload(value: tuple | None) -> list[object] | None:
    if value is None:
        return None
    return [fraction_payload(value[0]), *value[1:]]


def optional_minimum_from_payload(
    payload: object,
    length: int,
    label: str,
) -> tuple | None:
    if payload is None:
        return None
    require(isinstance(payload, list) and len(payload) == length, f"bad {label} payload")
    require(all(isinstance(item, int) for item in payload[1:]), f"noninteger {label} coordinates")
    return (fraction_from_payload(payload[0], label), *payload[1:])


def row_payload(row: dict[str, object]) -> dict[str, object]:
    return {
        "K": list(row["K"]),
        "r": int(row["r"]),
        "m": fraction_payload(row["m"]),
        "Va": int(row["Va"]),
        "a_nodes": int(row["a_nodes"]),
        "b_nodes": int(row["b_nodes"]),
        "triple_candidates": int(row["triple_candidates"]),
        "zero_after_a": int(row["zero_after_a"]),
        "zero_after_b": int(row["zero_after_b"]),
        "max_Vb": int(row["max_Vb"]),
        "max_cap": int(row["max_cap"]),
        "min_after_a": optional_minimum_payload(row["min_after_a"]),
        "min_after_b": optional_minimum_payload(row["min_after_b"]),
        "positive": int(row["positive"]),
        "zero": int(row["zero"]),
        "minimum": optional_minimum_payload(row["minimum"]),
        "failures": [list(family) for family in row["failures"]],
        "sample_records": list(row["sample_records"]),
    }


def row_from_payload(payload: object) -> dict[str, object]:
    require(isinstance(payload, dict), "checkpoint row is not an object")
    expected_keys = {
        "K",
        "r",
        "m",
        "Va",
        "a_nodes",
        "b_nodes",
        "triple_candidates",
        "zero_after_a",
        "zero_after_b",
        "max_Vb",
        "max_cap",
        "min_after_a",
        "min_after_b",
        "positive",
        "zero",
        "minimum",
        "failures",
        "sample_records",
    }
    require(set(payload) == expected_keys, "checkpoint row keys changed")
    K_payload = payload["K"]
    require(isinstance(K_payload, list) and len(K_payload) == 3, "bad checkpoint K")
    require(all(isinstance(item, int) for item in K_payload), "noninteger checkpoint K")
    K = tuple(K_payload)
    integer_keys = (
        "r",
        "Va",
        "a_nodes",
        "b_nodes",
        "triple_candidates",
        "zero_after_a",
        "zero_after_b",
        "max_Vb",
        "max_cap",
        "positive",
        "zero",
    )
    require(all(isinstance(payload[key], int) for key in integer_keys), f"noninteger row field K={K}")
    failures_payload = payload["failures"]
    require(isinstance(failures_payload, list), f"bad failure ledger K={K}")
    failures = []
    for family in failures_payload:
        require(isinstance(family, list), f"bad failure family K={K}")
        require(all(isinstance(item, int) for item in family), f"noninteger failure family K={K}")
        failures.append(tuple(family))
    samples_payload = payload["sample_records"]
    require(isinstance(samples_payload, list), f"bad samples K={K}")
    require(all(isinstance(item, str) for item in samples_payload), f"nonstr sample K={K}")
    row = {
        "K": K,
        "r": payload["r"],
        "m": fraction_from_payload(payload["m"], f"m K={K}"),
        "Va": payload["Va"],
        "a_nodes": payload["a_nodes"],
        "b_nodes": payload["b_nodes"],
        "triple_candidates": payload["triple_candidates"],
        "zero_after_a": payload["zero_after_a"],
        "zero_after_b": payload["zero_after_b"],
        "max_Vb": payload["max_Vb"],
        "max_cap": payload["max_cap"],
        "min_after_a": optional_minimum_from_payload(payload["min_after_a"], 2, f"min_after_a K={K}"),
        "min_after_b": optional_minimum_from_payload(payload["min_after_b"], 3, f"min_after_b K={K}"),
        "positive": payload["positive"],
        "zero": payload["zero"],
        "minimum": optional_minimum_from_payload(payload["minimum"], 4, f"minimum K={K}"),
        "failures": tuple(failures),
        "sample_records": tuple(samples_payload),
    }
    require(row["positive"] + row["zero"] == row["triple_candidates"], f"bad sign ledger K={K}")
    require(len(row["sample_records"]) == 5, f"bad cached sample count K={K}")
    require(row["zero"] == len(row["failures"]), f"bad cached failure count K={K}")
    return row


def row_sha256(row: dict[str, object]) -> str:
    return hashlib.sha256(canonical_json(row_payload(row)).encode()).hexdigest()


def checkpoint_record(row: dict[str, object], script_sha256: str) -> dict[str, object]:
    payload = row_payload(row)
    return {
        "schema": CHECKPOINT_SCHEMA,
        "core_sha256": CORE_SHA256,
        "four_small_sha256": FOUR_SMALL_SHA256,
        "script_sha256": script_sha256,
        "K": payload["K"],
        "row_sha256": hashlib.sha256(canonical_json(payload).encode()).hexdigest(),
        "row": payload,
    }


def append_checkpoint(path: Path, row: dict[str, object], script_sha256: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    line = canonical_json(checkpoint_record(row, script_sha256)) + "\n"
    with path.open("a", encoding="utf-8") as handle:
        handle.write(line)
        handle.flush()
        os.fsync(handle.fileno())


def read_checkpoint(
    path: Path,
    script_sha256: str,
    residual: tuple[tuple[int, int, int], ...],
) -> dict[tuple[int, int, int], dict[str, object]]:
    if not path.exists():
        return {}
    raw_lines = path.read_bytes().splitlines()
    rows: dict[tuple[int, int, int], dict[str, object]] = {}
    for index, raw_line in enumerate(raw_lines):
        if not raw_line.strip():
            continue
        try:
            record = json.loads(raw_line)
        except (UnicodeDecodeError, json.JSONDecodeError) as error:
            require(index == len(raw_lines) - 1, f"malformed checkpoint line {index + 1}: {error}")
            continue
        require(isinstance(record, dict), f"checkpoint line {index + 1} is not an object")
        require(record.get("schema") == CHECKPOINT_SCHEMA, "checkpoint schema changed")
        require(record.get("core_sha256") == CORE_SHA256, "checkpoint core hash changed")
        require(record.get("four_small_sha256") == FOUR_SMALL_SHA256, "checkpoint premise hash changed")
        require(record.get("script_sha256") == script_sha256, "checkpoint script hash changed; use --fresh")
        row = row_from_payload(record.get("row"))
        K = row["K"]
        require(K in residual, f"checkpoint contains nonresidual K={K}")
        require(record.get("K") == list(K), f"checkpoint K wrapper mismatch K={K}")
        digest = row_sha256(row)
        require(record.get("row_sha256") == digest, f"checkpoint row hash mismatch K={K}")
        if K in rows:
            require(row_sha256(rows[K]) == digest, f"conflicting duplicate checkpoint K={K}")
        rows[K] = row
    return rows


def render_report(
    rows: list[dict[str, object]],
    residual: tuple[tuple[int, int, int], ...],
    script_sha256: str,
) -> str:
    rows = sorted(rows, key=lambda row: row["K"])
    require(tuple(row["K"] for row in rows) == residual, "row K coverage/order changed")
    require(sum(int(row["zero_after_a"]) for row in rows) == 0, "covering one-external alarm")
    require(sum(int(row["zero_after_b"]) for row in rows) == 0, "covering two-external alarm")
    total_candidates = sum(int(row["triple_candidates"]) for row in rows)
    total_positive = sum(int(row["positive"]) for row in rows)
    total_zero = sum(int(row["zero"]) for row in rows)
    failures = tuple(family for row in rows for family in row["failures"])
    sample_records = tuple(
        record
        for row in rows
        for record in row["sample_records"]
    )
    require(total_candidates == 1_357_920, f"candidate total changed: {total_candidates}")
    require(len(sample_records) == 110, f"bad manifest sample count {len(sample_records)}")
    manifest_sha256 = hashlib.sha256("\n".join(sample_records).encode()).hexdigest()
    row_certificates = tuple((row["K"], row_sha256(row)) for row in rows)
    row_manifest_sha256 = hashlib.sha256(
        "\n".join(canonical_json(row_payload(row)) for row in rows).encode()
    ).hexdigest()
    positive_minima = [
        (row["minimum"][0], row["K"], *row["minimum"][1:])
        for row in rows
        if row["minimum"] is not None
    ]
    global_minimum = min(positive_minima) if positive_minima else None
    order = sorted(
        rows,
        key=lambda row: (
            int(row["triple_candidates"]),
            int(row["b_nodes"]),
            int(row["a_nodes"]),
            row["K"],
        ),
    )
    lex = tuple(sorted(residual))
    rank = {row["K"]: index for index, row in enumerate(order)}
    flips = sum(
        rank[left] > rank[right]
        for index, left in enumerate(lex)
        for right in lex[index + 1 :]
    )
    global_payload = {
        "row_manifest_sha256": row_manifest_sha256,
        "crosscheck_manifest_sha256": manifest_sha256,
        "total_candidates": total_candidates,
        "total_positive": total_positive,
        "total_zero": total_zero,
        "failures": [list(family) for family in failures],
        "global_minimum": (
            None
            if global_minimum is None
            else [
                fraction_payload(global_minimum[0]),
                list(global_minimum[1]),
                global_minimum[2],
                global_minimum[3],
                global_minimum[4],
            ]
        ),
    }
    global_certificate_sha256 = hashlib.sha256(canonical_json(global_payload).encode()).hexdigest()

    lines = [
        "THM-741 THREE-SMALL EXACT POSITIVITY SWEEP",
        "=" * 92,
        f"dependency_sha256={CORE_SHA256}",
        f"four_small_output_sha256={FOUR_SMALL_SHA256}",
        "premise: anchor shadow plus exact four-small theorem; whole flood bodies remain 3/21",
        "three-small completed bases: total=35 anchor-shadowed=13 residual=22",
        "K ; r ; m ; Va ; a-nodes ; b-nodes ; triples ; positive/zero ; "
        "minimum ; maxVb ; maxcap ; min-m1 ; min-m2",
    ]
    for row in order:
        min1 = row["min_after_a"]
        min2 = row["min_after_b"]
        require(min1 is not None and min2 is not None, f"missing intermediate minimum K={row['K']}")
        minimum = row["minimum"]
        minimum_text = (
            "none"
            if minimum is None
            else f"{minimum[0]}@{minimum[1]},{minimum[2]},{minimum[3]}"
        )
        lines.append(
            f"  {''.join(map(str, row['K']))} ; {row['r']} ; {row['m']} ; {row['Va']} ; "
            f"{row['a_nodes']} ; {row['b_nodes']} ; {row['triple_candidates']} ; "
            f"{row['positive']}/{row['zero']} ; {minimum_text} ; "
            f"{row['max_Vb']} ; {row['max_cap']} ; "
            f"{min1[0]}@{min1[1]} ; {min2[0]}@{min2[1]},{min2[2]}"
        )
    lines.append("row certificates (canonical complete per-K payloads):")
    for K, digest in row_certificates:
        lines.append(f"  K={''.join(map(str, K))};sha256={digest}")
    lines.append(f"row_manifest_sha256={row_manifest_sha256}")
    lines.append(
        "TOTALS: a_nodes={} b_nodes={} finite_exact_triples={} zero_after_a={} zero_after_b={}".format(
            sum(int(row["a_nodes"]) for row in rows),
            sum(int(row["b_nodes"]) for row in rows),
            sum(int(row["triple_candidates"]) for row in rows),
            sum(int(row["zero_after_a"]) for row in rows),
            sum(int(row["zero_after_b"]) for row in rows),
        )
    )
    if global_minimum is None:
        lines.append(
            f"EXACT TOTALS: positive={total_positive} zero={total_zero} "
            f"covering_failures={len(failures)} global_minimum=none"
        )
    else:
        lines.append(
            f"EXACT TOTALS: positive={total_positive} zero={total_zero} "
            f"covering_failures={len(failures)} global_minimum={global_minimum[0]} "
            f"at K={global_minimum[1]},a={global_minimum[2]},b={global_minimum[3]},c={global_minimum[4]}"
        )
    if failures:
        lines.append("ZERO/FAILURE FAMILIES:")
        for family in sorted(failures):
            lines.append(f"  {family}")
    else:
        lines.append("ZERO/FAILURE FAMILIES: none")
    lines.append("sparse/full cross-check rule: flat indices 0,1/4,1/2,3/4,last in every K-bank")
    lines.append("sparse/full cross-check manifest records:")
    for record in sample_records:
        lines.append(f"  {record}")
    lines.append(f"crosscheck_samples={len(sample_records)} mismatches=0 manifest_sha256={manifest_sha256}")
    lines.append(f"global_certificate_sha256={global_certificate_sha256}")
    if total_zero == 0 and not failures and total_positive == total_candidates:
        lines.append("VERDICT: PROVED THREE-SMALL SHADOW; every exact final triple is positive")
    else:
        lines.append("VERDICT: THREE-SMALL SHADOW NOT CLOSED; zero/failure ledger above")
    lines.append("Tournament Analysis vertices: 22 residual completed-family K-bases")
    lines.append("pair observable: lex(-triple_count,-b_nodes,-a_nodes); gauge points toward less work; tie=lex K")
    lines.append(
        f"fingerprint: score_hist={{0..21:1}}, directed_3cycles=0, SCC_sizes=22x1, "
        f"edge_flips_vs_lex={flips}, Hamiltonian_paths=1"
    )
    lines.append("tie Hamiltonian path=" + " -> ".join("".join(map(str, row["K"])) for row in order))
    lines.append("kept: completed-family base and exact finite workload; destroyed: interval geometry and margins")
    lines.append("challenged vertices: runners/Fano flags/root edges do not preserve the completed-family bank")
    lines.append(f"source_sha256={script_sha256}")
    lines.append("ALL EXACT WORKLOAD CHECKS PASSED")
    return "\n".join(lines) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--workers",
        type=int,
        default=min(4, max(1, (os.cpu_count() or 2) // 2)),
        help="completed-family K shards evaluated concurrently",
    )
    parser.add_argument(
        "--checkpoint",
        type=Path,
        default=DEFAULT_CHECKPOINT,
        help="resume-safe per-K JSONL checkpoint (default: TMPDIR)",
    )
    parser.add_argument(
        "--fresh",
        action="store_true",
        help="remove any existing checkpoint before computing",
    )
    parser.add_argument(
        "--render-only",
        action="store_true",
        help="require a complete checkpoint and render without recomputing",
    )
    args = parser.parse_args()
    require(args.workers >= 1, "workers must be positive")
    require(not (args.fresh and args.render_only), "--fresh and --render-only are incompatible")

    core = load_core()
    require(core.S2 == F(99, 70) and core.S2 * core.S2 > 2, "bad sqrt(2) majorant")
    premise_payload = FOUR_SMALL_OUTPUT.read_text()
    require("PROVED FOUR-SMALL SHADOW" in premise_payload, "missing four-small premise")

    all_K = tuple(combinations(range(1, 8), 3))
    shadowed = tuple(K for K in all_K if any(anchor <= frozenset(K) for anchor in ANCHORS))
    residual = tuple(K for K in all_K if K not in shadowed)
    require(len(shadowed) == 13 and len(residual) == 22, "bad three-small shadow split")
    require(all(len(frozenset(K) & frozenset((5, 6, 7))) <= 1 for K in residual), "bad residual K")

    script_sha256 = sha256(Path(__file__).resolve())
    checkpoint_path = args.checkpoint.expanduser().resolve()
    if args.fresh and checkpoint_path.exists():
        checkpoint_path.unlink()
    rows_by_K = read_checkpoint(checkpoint_path, script_sha256, residual)
    missing = tuple(K for K in residual if K not in rows_by_K)
    if args.render_only:
        require(not missing, f"incomplete checkpoint: missing {len(missing)} K-banks")
    elif missing and args.workers == 1:
        for K in missing:
            row = bank_work(core, K)
            append_checkpoint(checkpoint_path, row, script_sha256)
            rows_by_K[K] = row
    elif missing:
        context = mp.get_context("spawn")
        with context.Pool(min(args.workers, len(missing)), initializer=init_worker) as pool:
            for row in pool.imap_unordered(bank_work_worker, missing, chunksize=1):
                K = row["K"]
                require(K not in rows_by_K, f"duplicate computed row K={K}")
                append_checkpoint(checkpoint_path, row, script_sha256)
                rows_by_K[K] = row
    require(tuple(sorted(rows_by_K)) == residual, "checkpoint/computation did not cover residual K-bases")
    rows = [rows_by_K[K] for K in residual]
    print(render_report(rows, residual, script_sha256), end="")


if __name__ == "__main__":
    main()
