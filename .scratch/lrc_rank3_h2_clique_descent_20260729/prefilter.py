#!/usr/bin/env python3
"""Cheap exact-cap prefilters for the K>=20 strongest-r3 H2 graph pilot."""

from __future__ import annotations

import hashlib
import importlib.util
import multiprocessing as mp
import os
from collections import Counter
from fractions import Fraction as F
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
PILOT_PATH = ROOT / ".scratch/lrc_rank3_h2_clique_descent_20260729/pilot.py"
PAIR_PATH = (
    ROOT
    / ".scratch/lrc_sparse_j6_gate_audit_20260729/"
    / "all_root_q5_pair_suffix_census.py"
)
PAIR_SHA256 = (
    "99a93add65caf8acb1658a2dd363fd04f98707c2d3c853ca0f56f4f577ac0a4c"
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def load(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, f"cannot load {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


require(hashlib.sha256(PAIR_PATH.read_bytes()).hexdigest() == PAIR_SHA256, "pair source changed")
P = load(PILOT_PATH, "rank3_h2_prefilter_pilot")
C = load(PAIR_PATH, "rank3_h2_prefilter_pair")


def ftext(value: F) -> str:
    return f"{value.numerator}/{value.denominator}"


def hunter_star_cap(q: tuple[F, ...], pair_cap: F) -> tuple[F, tuple[F, ...]]:
    """Evaluate THM-2897's exact piecewise-linear G_5 cap."""

    require(len(q) == 5 and all(q[i] >= q[i + 1] for i in range(4)), "bad q ranks")
    upper = min(q[0], pair_cap)

    def clip(value: F) -> F:
        return max(F(0), min(upper, value))

    candidates = {
        F(0),
        upper,
        clip(pair_cap / 2),
        *(
            value
            for rank in q[1:]
            for value in (clip(rank), clip(pair_cap - rank))
        ),
    }

    def invoice(a: F) -> F:
        return a + sum(
            (min(a, rank, pair_cap - a) for rank in q[1:]),
            F(0),
        )

    decorated = tuple(sorted((invoice(a), a) for a in candidates))
    maximum = max(decorated)
    require(
        all(value <= maximum[0] for value, _ in decorated),
        "Hunter breakpoint maximum failed",
    )
    return maximum[0], tuple(a for value, a in decorated if value == maximum[0])


def profile(row: dict[str, object]) -> dict[str, object]:
    body = row["body"]
    root_good, _, _ = P.S.G.T.CORE.good_norm(body)
    carrier = P.S.R.subtract_local(root_good, row["apex"])
    require(P.mass(carrier) == row["h"], "carrier mass changed")
    excluded = set(row["prefix"])
    ranked, tail_single, tail_first = C.sealed_ranking(carrier, excluded)
    require(tuple(ranked[:5]) == row["top5"], "top-five ranking changed")
    pair = C.pair_cap_from_ranking(carrier, ranked, tail_single)
    b2 = pair["cap"]
    q = tuple(value for value, _ in row["top5"])
    r3 = sum(q[:3], F(0))
    theta = row["h"] - r3
    hunter, maximizers = hunter_star_cap(q, b2)
    margins = {
        "no_edge": theta - b2,
        "partition": row["h"] - r3 - b2,
        "q5_pair": row["h"] - q[4] - 2 * b2,
        "hunter": row["h"] - hunter,
    }
    require(margins["no_edge"] == margins["partition"], "margin alias changed")
    return {
        **row,
        "b2": b2,
        "pair_head": pair["head"],
        "pair_tail": pair["tail"],
        "pair_maximizer": pair["maximizer"],
        "pair_paid": pair["paid"],
        "pair_paid_digest": pair["paid_digest"],
        "pair_tail_first": tail_first,
        "hunter": hunter,
        "hunter_maximizers": maximizers,
        "margins": margins,
    }


def ledger_line(row: dict[str, object]) -> str:
    return (
        f"E={row['body']};K={row['K']};rank={row['rank']};a={row['apex']};"
        f"P={row['prefix']};h={ftext(row['h'])};"
        f"B2={ftext(row['b2'])};Bhead={ftext(row['pair_head'])};"
        f"Btail={ftext(row['pair_tail'])};Bmax={row['pair_maximizer']};"
        f"paid={row['pair_paid']};paid_digest={row['pair_paid_digest']};"
        f"G5={ftext(row['hunter'])};Garg={tuple(map(ftext,row['hunter_maximizers']))};"
        + ";".join(f"{name}={ftext(value)}" for name, value in row["margins"].items())
        + "\n"
    )


def main() -> None:
    workers = min(os.cpu_count() or 1, 8)
    with mp.get_context("spawn").Pool(workers) as pool:
        rows = list(pool.imap(profile, P.TARGET_ROWS, chunksize=1))
    require(
        tuple((row["body"], row["rank"]) for row in rows)
        == tuple((row["body"], row["rank"]) for row in P.TARGET_ROWS),
        "worker order changed",
    )
    digest = hashlib.sha256(b"LRC14/j6/highK-r3-cap-prefilter/v1\n")
    for row in rows:
        line = ledger_line(row)
        print("ROW;" + line.rstrip())
        digest.update(line.encode())
    bits = Counter(
        tuple(row["margins"][name] > 0 for name in ("partition", "q5_pair", "hunter"))
        for row in rows
    )
    union = [
        row
        for row in rows
        if any(row["margins"][name] > 0 for name in ("partition", "q5_pair", "hunter"))
    ]
    survivors = [
        row
        for row in rows
        if all(row["margins"][name] <= 0 for name in ("partition", "q5_pair", "hunter"))
    ]
    print("LRC14 j6 K>=20 strongest-r3 cap prefilter")
    print(
        f"branches={len(rows)};partition_no_edge="
        f"{sum(row['margins']['partition']>0 for row in rows)};"
        f"q5_pair={sum(row['margins']['q5_pair']>0 for row in rows)};"
        f"hunter={sum(row['margins']['hunter']>0 for row in rows)};"
        f"union={len(union)};survivors={len(survivors)}"
    )
    print(f"bit_histogram={tuple(sorted(bits.items()))}")
    print(
        "survivor_keys="
        + repr(
            tuple(
                (
                    row["body"],
                    row["rank"],
                    row["apex"],
                    tuple(ftext(row["margins"][name]) for name in ("partition", "q5_pair", "hunter")),
                )
                for row in survivors
            )
        )
    )
    print(f"pair_paid={sum(row['pair_paid'] for row in rows)}")
    print(f"canonical_ledger_sha256={digest.hexdigest()}")
    print("scope=exact upper-cap prefilter on K>=20 strongest-r3 branches")
    print("all_exact_controls=PASS")


if __name__ == "__main__":
    main()
