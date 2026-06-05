#!/usr/bin/env python3
"""Exact rational LRC tools for small candidate and shard checks.

Conventions:
  * ``moving`` is the number of non-stationary speeds in the set.
  * ``n_total = moving + 1`` includes the stationary observer.
  * The LRC threshold is ``delta = 1/n_total`` unless overridden.

The exact maximin scan enumerates all breakpoints where the lower envelope of
``||v t||`` can change slope: pair-sum, pair-difference, antipode, and endpoints.
This is intentionally small and auditable, not optimized for very large boxes.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import time
from collections import Counter, defaultdict, deque
from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations
from pathlib import Path
from typing import Iterable


ONE = Fraction(1, 1)
VERIFIER_VERSION = "lrc_exact_tools_v1"
CANDIDATE_FORMAT = "0,1,(2m+1)/(2v_i),m/(v_i+v_j),m/|v_i-v_j|"


def parse_fraction(text: str) -> Fraction:
    if "/" in text:
        a, b = text.split("/", 1)
        return Fraction(int(a), int(b))
    return Fraction(text)


def fmt_fraction(x: Fraction) -> str:
    return str(x.numerator) if x.denominator == 1 else f"{x.numerator}/{x.denominator}"


def verifier_fingerprint() -> dict[str, str]:
    path = Path(__file__)
    digest = hashlib.blake2b(path.read_bytes(), digest_size=12).hexdigest()
    return {
        "tool": path.name,
        "version": VERIFIER_VERSION,
        "hash_blake2b_12": digest,
    }


def parse_speeds(text: str) -> tuple[int, ...]:
    speeds = tuple(sorted(int(x) for x in text.replace(",", " ").split()))
    if not speeds:
        raise ValueError("empty speed set")
    if speeds[0] <= 0:
        raise ValueError("speeds must be positive nonzero integers")
    if len(set(speeds)) != len(speeds):
        raise ValueError("speeds must be distinct")
    return speeds


def primitive_gcd(speeds: Iterable[int]) -> int:
    g = 0
    for v in speeds:
        g = math.gcd(g, abs(v))
    return g


def normalize_speeds(speeds: Iterable[int]) -> tuple[int, ...]:
    raw = tuple(sorted(abs(int(v)) for v in speeds if int(v) != 0))
    if len(set(raw)) != len(raw):
        raise ValueError(f"non-distinct normalized speeds: {raw}")
    g = primitive_gcd(raw)
    if g <= 0:
        raise ValueError("empty speed set")
    return tuple(v // g for v in raw)


def frac_part(x: Fraction) -> Fraction:
    return x - math.floor(x)


def norm_to_integer(x: Fraction) -> Fraction:
    r = frac_part(x)
    return min(r, ONE - r)


def score_time(speeds: tuple[int, ...], t: Fraction) -> Fraction:
    return min(norm_to_integer(v * t) for v in speeds)


def active_records(
    speeds: tuple[int, ...], t: Fraction, score: Fraction | None = None
) -> tuple[dict[str, object], ...]:
    if score is None:
        score = score_time(speeds, t)
    out: list[dict[str, object]] = []
    for idx, v in enumerate(speeds):
        phase = frac_part(v * t)
        dist = min(phase, ONE - phase)
        if dist == score:
            sign = "-" if phase <= Fraction(1, 2) else "+"
            out.append(
                {
                    "index": idx,
                    "speed": v,
                    "phase": fmt_fraction(phase),
                    "side": sign,
                }
            )
    return tuple(out)


def candidate_times(speeds: tuple[int, ...]) -> set[Fraction]:
    out: set[Fraction] = {Fraction(0), Fraction(1)}
    for v in speeds:
        for k in range(v):
            out.add(Fraction(2 * k + 1, 2 * v))
        for k in range(v + 1):
            out.add(Fraction(k, v))
    for a, b in combinations(speeds, 2):
        for den in (a + b, abs(a - b)):
            if den <= 0:
                continue
            for m in range(den + 1):
                out.add(Fraction(m, den))
    return {t for t in out if 0 <= t <= 1}


@dataclass(frozen=True)
class ExactResult:
    speeds: tuple[int, ...]
    threshold: Fraction
    maximin: Fraction
    witnesses: tuple[Fraction, ...]
    candidate_count: int

    @property
    def best_time(self) -> Fraction:
        return self.witnesses[0]

    @property
    def margin(self) -> Fraction:
        return self.maximin - self.threshold


def exact_maximin(
    speeds: tuple[int, ...], threshold: Fraction | None = None
) -> ExactResult:
    if threshold is None:
        threshold = Fraction(1, len(speeds) + 1)
    candidates = candidate_times(speeds)
    best = Fraction(-1)
    witnesses: list[Fraction] = []
    for t in sorted(candidates):
        score = score_time(speeds, t)
        if score > best:
            best = score
            witnesses = [t]
        elif score == best:
            witnesses.append(t)
    return ExactResult(
        speeds=speeds,
        threshold=threshold,
        maximin=best,
        witnesses=tuple(witnesses),
        candidate_count=len(candidates),
    )


def contact_profile(result: ExactResult) -> dict[str, object]:
    contacts = []
    denominators: set[int] = set()
    for t in result.witnesses:
        recs = active_records(result.speeds, t, result.maximin)
        denominators.add(t.denominator)
        contacts.append(
            {
                "time": fmt_fraction(t),
                "denominator": t.denominator,
                "active": list(recs),
            }
        )
    classes: dict[int, list[tuple[tuple[int, str], ...]]] = defaultdict(list)
    for item in contacts:
        cls = tuple((int(r["index"]), str(r["side"])) for r in item["active"])
        classes[int(item["denominator"])].append(cls)
    return {
        "D": sorted(denominators),
        "contacts": contacts,
        "denominator_classes": {
            str(q): [list(map(list, cls)) for cls in sorted(set(v))]
            for q, v in sorted(classes.items())
        },
    }


def gf2_rank(rows: list[int]) -> int:
    basis: dict[int, int] = {}
    rank = 0
    for x in rows:
        y = x
        while y:
            bit = y.bit_length() - 1
            if bit not in basis:
                basis[bit] = y
                rank += 1
                break
            y ^= basis[bit]
    return rank


def contact_graph_stats(result: ExactResult) -> dict[str, object]:
    k = len(result.speeds)
    contact_nodes = []
    edges: list[tuple[int, int]] = []
    incidence_rows: list[int] = []
    for c_idx, t in enumerate(result.witnesses):
        active = active_records(result.speeds, t, result.maximin)
        mask = 0
        for rec in active:
            idx = int(rec["index"])
            mask |= 1 << idx
            edges.append((idx, c_idx))
        incidence_rows.append(mask)
        contact_nodes.append(fmt_fraction(t))

    # Undirected bipartite connectivity/cyclomatic rank.
    total_nodes = k + len(contact_nodes)
    adj = [[] for _ in range(total_nodes)]
    for runner, c_idx in edges:
        u = runner
        v = k + c_idx
        adj[u].append(v)
        adj[v].append(u)
    seen = [False] * total_nodes
    components = 0
    for i in range(total_nodes):
        if seen[i]:
            continue
        components += 1
        q = deque([i])
        seen[i] = True
        while q:
            u = q.popleft()
            for v in adj[u]:
                if not seen[v]:
                    seen[v] = True
                    q.append(v)
    cyclomatic = len(edges) - total_nodes + components
    return {
        "runner_count": k,
        "contact_count": len(contact_nodes),
        "edge_count": len(edges),
        "components": components,
        "connected": components == 1,
        "f2_rank": gf2_rank(incidence_rows),
        "cyclomatic_rank": cyclomatic,
    }


def half_turn_tournament(speeds: tuple[int, ...], t: Fraction) -> list[int]:
    n = len(speeds)
    out = [0] * n
    phases = [frac_part(v * t) for v in speeds]
    for i in range(n):
        for j in range(i + 1, n):
            diff = frac_part(phases[j] - phases[i])
            if diff == Fraction(1, 2):
                i_to_j = i < j
            else:
                i_to_j = diff < Fraction(1, 2)
            if i_to_j:
                out[i] |= 1 << j
            else:
                out[j] |= 1 << i
    return out


def scc_sizes(adj_bits: list[int]) -> list[int]:
    n = len(adj_bits)
    rev = [0] * n
    for i, bits in enumerate(adj_bits):
        for j in range(n):
            if bits & (1 << j):
                rev[j] |= 1 << i

    seen = [False] * n
    order: list[int] = []

    def dfs1(u: int) -> None:
        seen[u] = True
        bits = adj_bits[u]
        for v in range(n):
            if bits & (1 << v) and not seen[v]:
                dfs1(v)
        order.append(u)

    def dfs2(u: int, comp: list[int]) -> None:
        seen[u] = True
        comp.append(u)
        bits = rev[u]
        for v in range(n):
            if bits & (1 << v) and not seen[v]:
                dfs2(v, comp)

    for u in range(n):
        if not seen[u]:
            dfs1(u)
    seen = [False] * n
    sizes = []
    for u in reversed(order):
        if not seen[u]:
            comp: list[int] = []
            dfs2(u, comp)
            sizes.append(len(comp))
    return sorted(sizes, reverse=True)


def hamiltonian_paths(adj_bits: list[int]) -> int:
    n = len(adj_bits)
    size = 1 << n
    dp = [[0] * n for _ in range(size)]
    for i in range(n):
        dp[1 << i][i] = 1
    for mask in range(size):
        for last in range(n):
            val = dp[mask][last]
            if not val:
                continue
            avail = adj_bits[last] & ~mask
            while avail:
                bit = avail & -avail
                nxt = bit.bit_length() - 1
                dp[mask | bit][nxt] += val
                avail ^= bit
    return sum(dp[-1])


def tournament_stats(speeds: tuple[int, ...], t: Fraction) -> dict[str, object]:
    adj = half_turn_tournament(speeds, t)
    H = hamiltonian_paths(adj)
    return {
        "scc_sizes": scc_sizes(adj),
        "hamiltonian_paths": H,
        "entropy_log2_H": math.log2(H) if H > 0 else float("-inf"),
    }


def certificate_for_speeds(
    speeds: tuple[int, ...], threshold: Fraction | None = None
) -> dict[str, object]:
    speeds = normalize_speeds(speeds)
    result = exact_maximin(speeds, threshold)
    contacts = contact_profile(result)
    graph = contact_graph_stats(result)
    tstats = tournament_stats(speeds, result.best_time)
    return {
        "schema": "lrc_exact_certificate_v1",
        "speeds": list(speeds),
        "moving_runners": len(speeds),
        "n_total": result.threshold.denominator
        if result.threshold.numerator == 1
        else len(speeds) + 1,
        "threshold": fmt_fraction(result.threshold),
        "maximin": fmt_fraction(result.maximin),
        "margin": fmt_fraction(result.margin),
        "witnesses": [fmt_fraction(t) for t in result.witnesses],
        "candidate_count": result.candidate_count,
        "contact_profile": contacts,
        "contact_graph": graph,
        "tournament_at_first_witness": tstats,
    }


def verify_certificate(path: Path) -> dict[str, object]:
    cert = json.loads(path.read_text())
    if cert.get("schema") != "lrc_exact_certificate_v1":
        return {"ok": False, "reason": "unsupported_schema", "schema": cert.get("schema")}
    speeds = tuple(int(v) for v in cert["speeds"])
    threshold = parse_fraction(str(cert["threshold"]))
    result = exact_maximin(speeds, threshold)
    expected_max = parse_fraction(str(cert["maximin"]))
    expected_witnesses = tuple(parse_fraction(str(t)) for t in cert["witnesses"])
    ok = result.maximin == expected_max and result.witnesses == expected_witnesses
    first_bad = None
    if result.maximin != expected_max:
        first_bad = {
            "field": "maximin",
            "expected": fmt_fraction(expected_max),
            "actual": fmt_fraction(result.maximin),
        }
    elif result.witnesses != expected_witnesses:
        first_bad = {
            "field": "witnesses",
            "expected": [fmt_fraction(t) for t in expected_witnesses],
            "actual": [fmt_fraction(t) for t in result.witnesses],
        }
    return {
        "ok": ok,
        "path": str(path),
        "speeds": list(speeds),
        "threshold": fmt_fraction(threshold),
        "maximin": fmt_fraction(result.maximin),
        "witness_count": len(result.witnesses),
        "candidate_count": result.candidate_count,
        "first_bad": first_bad,
    }


def cover_scan(
    speeds: tuple[int, ...],
    delta: Fraction,
    qmax: int,
    near_limit: int = 12,
) -> dict[str, object]:
    speeds = normalize_speeds(speeds)
    last_uncovered: list[dict[str, object]] = []
    minimal_q = None
    for q in range(1, qmax + 1):
        uncovered = []
        for r in range(q):
            left = Fraction(r, q)
            right = Fraction(r + 1, q)
            best_v = None
            best_value = None
            for v in speeds:
                endpoint_max = max(norm_to_integer(v * left), norm_to_integer(v * right))
                cert_value = endpoint_max + Fraction(abs(v), 2 * q)
                if best_value is None or cert_value < best_value:
                    best_value = cert_value
                    best_v = v
            assert best_value is not None and best_v is not None
            slack = delta - best_value
            if not (best_value < delta):
                mid = Fraction(2 * r + 1, 2 * q)
                uncovered.append(
                    {
                        "interval": [fmt_fraction(left), fmt_fraction(right)],
                        "midpoint": fmt_fraction(mid),
                        "best_speed": best_v,
                        "cert_value": fmt_fraction(best_value),
                        "slack": fmt_fraction(slack),
                    }
                )
        last_uncovered = uncovered
        if not uncovered:
            minimal_q = q
            break
    near = sorted(
        last_uncovered,
        key=lambda x: (parse_fraction(str(x["slack"])), x["interval"][0]),
        reverse=True,
    )[:near_limit]
    return {
        "schema": "lrc_modular_cover_scan_v1",
        "speeds": list(speeds),
        "delta": fmt_fraction(delta),
        "qmax": qmax,
        "minimal_certified_q": minimal_q,
        "last_q": minimal_q if minimal_q is not None else qmax,
        "uncovered_count_at_last_q": len(last_uncovered),
        "near_misses": near,
    }


def stable_shard(speeds: tuple[int, ...], shards: int) -> int:
    msg = ",".join(map(str, speeds)).encode()
    digest = hashlib.blake2b(msg, digest_size=8).digest()
    return int.from_bytes(digest, "big") % shards


def bit_sizes_for_result(result: ExactResult) -> dict[str, int]:
    nums = [result.maximin.numerator, result.best_time.numerator]
    dens = [result.maximin.denominator, result.best_time.denominator]
    return {
        "max_numerator_bits": max(abs(x).bit_length() for x in nums),
        "max_denominator_bits": max(abs(x).bit_length() for x in dens),
    }


def run_pilot(args: argparse.Namespace) -> dict[str, object]:
    run_id = args.run_id or time.strftime("%Y%m%dT%H%M%SZ", time.gmtime())
    outdir = Path("runs") / run_id / "primitive_n14"
    outdir.mkdir(parents=True, exist_ok=True)
    jsonl_path = outdir / f"{args.shard_index:04d}.jsonl"
    summary_path = outdir / f"{args.shard_index:04d}.summary.json"
    threshold = parse_fraction(args.delta)
    start = time.perf_counter()
    deadline = start + args.time_budget_s if args.time_budget_s > 0 else None

    total_seen = 0
    shard_seen = 0
    processed = 0
    primitive_skipped = 0
    failures = 0
    tight = 0
    overflow_count = 0
    timeout_count = 0
    max_num_bits = 0
    max_den_bits = 0
    best_margin: Fraction | None = None
    best_rows: list[dict[str, object]] = []
    failure_mode = "completed"

    with jsonl_path.open("w") as fh:
        for combo in combinations(range(1, args.max_speed + 1), args.moving):
            total_seen += 1
            if stable_shard(combo, args.shards) != args.shard_index:
                continue
            shard_seen += 1
            if primitive_gcd(combo) != 1:
                primitive_skipped += 1
                continue
            if deadline is not None and time.perf_counter() >= deadline:
                timeout_count += 1
                failure_mode = "time_budget_exhausted"
                break
            result = exact_maximin(combo, threshold)
            bits = bit_sizes_for_result(result)
            max_num_bits = max(max_num_bits, bits["max_numerator_bits"])
            max_den_bits = max(max_den_bits, bits["max_denominator_bits"])
            contacts = contact_profile(result)
            graph = contact_graph_stats(result)
            status = "strict"
            if result.maximin < threshold:
                failures += 1
                status = "counterexample"
            elif result.maximin == threshold:
                tight += 1
                status = "tight"
            margin = result.margin
            row = {
                "schema": "lrc_primitive_pilot_row_v1",
                "run_id": run_id,
                "shard_index": args.shard_index,
                "shards": args.shards,
                "moving_runners": args.moving,
                "n_total": args.n_total,
                "threshold": fmt_fraction(threshold),
                "max_speed": args.max_speed,
                "speeds": list(combo),
                "status": status,
                "maximin": fmt_fraction(result.maximin),
                "margin": fmt_fraction(margin),
                "first_witness": fmt_fraction(result.best_time),
                "witness_count": len(result.witnesses),
                "candidate_count": result.candidate_count,
                "active_count": len(active_records(combo, result.best_time, result.maximin)),
                "D": contacts["D"],
                "contact_graph": graph,
                "contact_matrix": contacts["contacts"],
                "rank_signature": {
                    "F2": graph["f2_rank"],
                    "beta1": graph["cyclomatic_rank"],
                    "contacts": graph["contact_count"],
                    "edges": graph["edge_count"],
                },
                "exact_candidate_format": CANDIDATE_FORMAT,
                "verifier": verifier_fingerprint(),
                **bits,
            }
            fh.write(json.dumps(row, sort_keys=True) + "\n")
            processed += 1
            if best_margin is None or margin < best_margin:
                best_margin = margin
                best_rows = [row]
            elif margin == best_margin and len(best_rows) < 8:
                best_rows.append(row)

    wall = time.perf_counter() - start
    summary = {
        "schema": "lrc_primitive_pilot_summary_v1",
        "run_id": run_id,
        "parameter_slice": {
            "moving_runners": args.moving,
            "n_total": args.n_total,
            "delta": fmt_fraction(threshold),
            "max_speed": args.max_speed,
            "shard_index": args.shard_index,
            "shards": args.shards,
            "time_budget_s": args.time_budget_s,
            "canonicalization": "sorted primitive gcd-normal form",
            "pruning": "stable blake2b shard filter; gcd primitive filter",
        },
        "certificate_format": "lrc_exact_certificate_v1",
        "exact_candidate_format": CANDIDATE_FORMAT,
        "verifier": verifier_fingerprint(),
        "output_jsonl": str(jsonl_path),
        "walltime_s": wall,
        "failure_mode": failure_mode,
        "total_combinations_seen": total_seen,
        "shard_combinations_seen": shard_seen,
        "primitive_skipped": primitive_skipped,
        "processed_rows": processed,
        "throughput_rows_per_s": processed / wall if wall else None,
        "counterexample_count": failures,
        "tight_count": tight,
        "timeout_count": timeout_count,
        "overflow_count": overflow_count,
        "max_numerator_bits": max_num_bits,
        "max_denominator_bits": max_den_bits,
        "best_margin": fmt_fraction(best_margin) if best_margin is not None else None,
        "best_rows": best_rows,
    }
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    return summary


def validation_table(candidates: list[tuple[str, tuple[int, ...]]], threshold: Fraction) -> list[dict[str, object]]:
    rows = []
    for label, speeds in candidates:
        result = exact_maximin(normalize_speeds(speeds), threshold)
        contacts = contact_profile(result)
        graph = contact_graph_stats(result)
        tstats = tournament_stats(result.speeds, result.best_time)
        rows.append(
            {
                "label": label,
                "speeds": list(result.speeds),
                "M": fmt_fraction(result.maximin),
                "margin": fmt_fraction(result.margin),
                "witnesses": [fmt_fraction(t) for t in result.witnesses],
                "D(A)": contacts["D"],
                "contact_denominators": len(contacts["D"]),
                "SCC": tstats["scc_sizes"],
                "F2_rank": graph["f2_rank"],
                "cyclomatic_rank": graph["cyclomatic_rank"],
                "tournament_entropy_log2H": round(float(tstats["entropy_log2_H"]), 6),
            }
        )
    return rows


def cmd_cert(args: argparse.Namespace) -> None:
    speeds = parse_speeds(args.speeds)
    threshold = parse_fraction(args.delta)
    cert = certificate_for_speeds(speeds, threshold)
    text = json.dumps(cert, indent=2, sort_keys=True) + "\n"
    if args.out:
        Path(args.out).parent.mkdir(parents=True, exist_ok=True)
        Path(args.out).write_text(text)
    print(text, end="")


def cmd_verify(args: argparse.Namespace) -> None:
    print(json.dumps(verify_certificate(Path(args.certificate)), sort_keys=True))


def cmd_cover(args: argparse.Namespace) -> None:
    speeds = parse_speeds(args.speeds)
    delta = parse_fraction(args.delta)
    print(json.dumps(cover_scan(speeds, delta, args.qmax, args.near), indent=2, sort_keys=True))


def cmd_pilot(args: argparse.Namespace) -> None:
    print(json.dumps(run_pilot(args), indent=2, sort_keys=True))


def cmd_table(args: argparse.Namespace) -> None:
    threshold = parse_fraction(args.delta)
    candidates = [
        ("AP13", tuple(range(1, 14))),
        ("Vstar_AP12_to_24", (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 24)),
        ("AP13_11_to_22", (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 13, 22)),
        ("AP13_13_to_26", (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 26)),
        ("seven_ladder", (7, 14, 21, 28, 35, 42, 49, 1, 2, 3, 4, 5, 6)),
    ]
    rows = validation_table(candidates, threshold)
    if args.out:
        Path(args.out).parent.mkdir(parents=True, exist_ok=True)
        Path(args.out).write_text(json.dumps(rows, indent=2, sort_keys=True) + "\n")
    print(json.dumps(rows, indent=2, sort_keys=True))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    sub = ap.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("cert")
    p.add_argument("--speeds", required=True)
    p.add_argument("--delta", default=None)
    p.add_argument("--out")
    p.set_defaults(func=cmd_cert)

    p = sub.add_parser("verify")
    p.add_argument("certificate")
    p.set_defaults(func=cmd_verify)

    p = sub.add_parser("cover")
    p.add_argument("--speeds", required=True)
    p.add_argument("--delta", required=True)
    p.add_argument("--qmax", type=int, default=256)
    p.add_argument("--near", type=int, default=12)
    p.set_defaults(func=cmd_cover)

    p = sub.add_parser("pilot")
    p.add_argument("--run-id")
    p.add_argument("--moving", type=int, default=14)
    p.add_argument("--n-total", type=int, default=15)
    p.add_argument("--delta", default="1/15")
    p.add_argument("--max-speed", type=int, default=24)
    p.add_argument("--shard-index", type=int, default=0)
    p.add_argument("--shards", type=int, default=256)
    p.add_argument("--time-budget-s", type=float, default=30.0)
    p.set_defaults(func=cmd_pilot)

    p = sub.add_parser("table")
    p.add_argument("--delta", default="1/14")
    p.add_argument("--out")
    p.set_defaults(func=cmd_table)

    args = ap.parse_args()
    if getattr(args, "delta", None) is None:
        args.delta = f"1/{len(parse_speeds(args.speeds)) + 1}"
    args.func(args)


if __name__ == "__main__":
    main()
