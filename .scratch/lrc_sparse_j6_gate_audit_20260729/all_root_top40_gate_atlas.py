#!/usr/bin/env python3
"""Discovery atlas for top-forty adaptive gates on all 3432 j=6 roots."""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import multiprocessing as mp
from collections import Counter
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
KERNEL = (
    ROOT
    / "04-computation/lrc14_thm2885_eight_body_top15_hitting_gate_codex_20260729.py"
)
KERNEL_SHA = "dff97f67b1104c25589802a6a2f216b6e7bfedd58eebfa1bcce615d59c1e872f"
FIRST = 15
BASE = 1600
COUNT = 40
S2 = F(99, 70)
BODIES = tuple(combinations(range(1, 15), 7))


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def load():
    require(
        hashlib.sha256(KERNEL.read_bytes()).hexdigest() == KERNEL_SHA,
        "single-comb kernel changed",
    )
    spec = importlib.util.spec_from_file_location("j6_top40_kernel", KERNEL)
    require(spec is not None and spec.loader is not None, "cannot load kernel")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


T = load()


def ftext(value: F) -> str:
    return f"{value.numerator}/{value.denominator}"


def ceiling(value: F) -> int:
    return (value.numerator + value.denominator - 1) // value.denominator


def stratum(body: tuple[int, ...]) -> str:
    return ("low", "one", "both")[int(13 in body) + int(14 in body)]


def profile(body: tuple[int, ...]) -> dict[str, object]:
    good, components, mass = T.CORE.good_norm(body)
    rows = T.coverages_many(good, range(FIRST, BASE + 1))
    ranked_base = sorted(rows, key=lambda item: (-item[0], item[1]))
    q40 = ranked_base[COUNT - 1][0]
    if q40 <= mass / 7:
        return {
            "body": body,
            "stratum": stratum(body),
            "status": "unsealable",
            "m": mass,
            "r": components,
            "q40": q40,
        }
    threshold = S2 * components / (7 * (q40 - mass / 7))
    tail_first = max(BASE + 1, ceiling(threshold))
    rows.extend(T.coverages_many(good, range(BASE + 1, tail_first)))
    require(
        mass / 7 + S2 * components / (7 * tail_first) <= q40,
        f"rank-forty tail did not seal: {body}",
    )
    top = tuple(sorted(rows, key=lambda item: (-item[0], item[1]))[:COUNT])
    margins = tuple(
        mass
        - sum((value for value, _ in top[k : k + 6]), F(0))
        for k in range(COUNT - 5)
    )
    least = next((k for k, margin in enumerate(margins) if margin > 0), None)
    if least is None:
        return {
            "body": body,
            "stratum": stratum(body),
            "status": "no_gate",
            "m": mass,
            "r": components,
            "top": top,
            "tail_first": tail_first,
            "threshold": threshold,
            "margins": margins,
        }
    by_speed = {speed: value for value, speed in rows}
    controls = tuple(
        dict.fromkeys(
            (
                FIRST,
                top[0][1],
                top[-1][1],
                tail_first - 1,
                BASE,
                top[COUNT // 2][1],
                FIRST + 1,
                BASE - 1,
            )
        )
    )[:4]
    require(
        len(controls) == 4
        and all(
            by_speed[speed] == T.coverage(good, speed)
            for speed in controls
        ),
        f"scalar controls changed: {body}",
    )
    return {
        "body": body,
        "stratum": stratum(body),
        "status": "gate",
        "m": mass,
        "r": components,
        "top": top,
        "K": least,
        "margin": margins[least],
        "K24": margins[24],
        "threshold": threshold,
        "tail_first": tail_first,
        "controls": len(controls),
    }


def ledger_line(row: dict[str, object]) -> str:
    if row["status"] != "gate":
        return f"E={row['body']};status={row['status']}\n"
    return (
        f"E={','.join(map(str, row['body']))};S={row['stratum']};"
        f"h={ftext(row['m'])};r={row['r']};K={row['K']};"
        f"margin={ftext(row['margin'])};K24={ftext(row['K24'])};"
        f"T={ftext(row['threshold'])};tail={row['tail_first']};top="
        + ",".join(
            f"{speed}:{ftext(value)}" for value, speed in row["top"]
        )
        + "\n"
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=4)
    args = parser.parse_args()
    require(len(BODIES) == 3432 and args.workers >= 1, "bad atlas universe")
    if args.workers == 1:
        rows = [profile(body) for body in BODIES]
    else:
        with mp.get_context("spawn").Pool(args.workers) as pool:
            rows = list(pool.map(profile, BODIES, chunksize=4))
    require(
        tuple(row["body"] for row in rows) == BODIES,
        "worker order changed",
    )
    failures = [row for row in rows if row["status"] != "gate"]
    gates = [row for row in rows if row["status"] == "gate"]
    distribution = tuple(sorted(Counter(row["K"] for row in gates).items()))
    stratified = tuple(
        (
            name,
            tuple(
                sorted(
                    Counter(
                        row["K"] for row in gates if row["stratum"] == name
                    ).items()
                )
            ),
        )
        for name in ("low", "one", "both")
    )
    sharp = min(
        gates,
        key=lambda row: (row["margin"], row["body"]),
    )
    maximum_k = max(
        gates,
        key=lambda row: (row["K"], tuple(-x for x in row["body"])),
    )
    maximum_tail = max(
        gates,
        key=lambda row: (row["tail_first"], tuple(-x for x in row["body"])),
    )
    digest = hashlib.sha256(b"LRC14/j6/all-root-top40-gates/v1\n")
    for row in rows:
        digest.update(ledger_line(row).encode())
    print("LRC14 J6 ALL-ROOT TOP-FORTY GATE DISCOVERY")
    print(
        f"universe={len(rows)};gates={len(gates)};failures={len(failures)};"
        f"controls={sum(row.get('controls', 0) for row in gates)}"
    )
    print(f"K_distribution={distribution}")
    print(f"stratified_K_distribution={stratified}")
    print(
        f"maximum_K={maximum_k['K']};E={maximum_k['body']};"
        f"margin={ftext(maximum_k['margin'])}"
    )
    print(
        f"minimum_margin={ftext(sharp['margin'])};E={sharp['body']};"
        f"K={sharp['K']}"
    )
    print(
        f"maximum_tail_first={maximum_tail['tail_first']};"
        f"E={maximum_tail['body']};T={ftext(maximum_tail['threshold'])}"
    )
    if failures:
        print(
            "failure_rows="
            + "|".join(
                f"{row['body']}:{row['status']}" for row in failures[:20]
            )
        )
    print(f"ledger_sha256={digest.hexdigest()}")
    print("scope=root hitting gates only;no suffix closure;not LRC14")


if __name__ == "__main__":
    main()
