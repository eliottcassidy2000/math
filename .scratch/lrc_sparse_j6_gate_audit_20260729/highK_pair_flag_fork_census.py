#!/usr/bin/env python3
"""Census THM-2897's three pair-cap flag forks on the 394 high-K core."""

from __future__ import annotations

import argparse
import ast
import hashlib
import importlib.util
import multiprocessing as mp
import os
from collections import Counter, defaultdict
from fractions import Fraction as F
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
CENSUS_PATH = (
    ROOT / ".scratch/lrc_sparse_j6_gate_audit_20260729/"
    "all_root_q5_pair_suffix_census.py"
)
CENSUS_SHA = "02e8802eb60f6d598903b84e47f454592a78001f5143c02450533c0899e06c07"
HIGH_OUTPUT = (
    ROOT / ".scratch/lrc_sparse_j6_gate_audit_20260729/"
    "all_root_q5_pair_suffix_highK.out"
)
HIGH_OUTPUT_SHA = "a399fe314d6ad37ee0855537ab16b4590f4101210fcf2f7295fcbba36fb40b39"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load(path: Path, expected: str, name: str):
    require(file_sha256(path) == expected, f"{path.name} changed")
    spec = importlib.util.spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, f"cannot load {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


C = load(CENSUS_PATH, CENSUS_SHA, "j6_highK_pair_flag_census")
require(file_sha256(HIGH_OUTPUT) == HIGH_OUTPUT_SHA, "high-K output changed")
OPEN_ROOT_ROWS = ast.literal_eval(
    next(
        line.split("=", 1)[1]
        for line in HIGH_OUTPUT.read_text().splitlines()
        if line.startswith("open_roots=")
    )
)
HIGH_BODIES = tuple(row[0] for row in OPEN_ROOT_ROWS)
require(len(HIGH_BODIES) == 62 and len(set(HIGH_BODIES)) == 62, "root set changed")
EXPECTED_MARGINAL_COUNTS = (392, 308, 158)
EXPECTED_BIT_HISTOGRAM = (
    ((False, False, False), 2),
    ((True, False, False), 84),
    ((True, True, False), 150),
    ((True, True, True), 158),
)
EXPECTED_LEDGER_DIGEST = (
    "3575f3b5c6850b2a72a791968526f6748d9bcc7aed3c083938ce6f79cde25d10"
)
EXPECTED_SURVIVORS: tuple[tuple[object, ...], ...] | None = (
    (
        (1, 3, 7, 8, 10, 11, 13),
        21,
        1,
        18,
        (
            "-1848989/5703417720",
            "-1265297389/65589303780",
            "-125320907/2851708860",
        ),
        "21237/140140",
        "5179/127512",
        "70819523/814773960",
    ),
    (
        (1, 4, 9, 10, 12, 13, 14),
        20,
        1,
        22,
        (
            "-157867/116035920",
            "-58997/2109744",
            "-3147959/58017960",
        ),
        "4643/25740",
        "1201/22932",
        "2423647/23207184",
    ),
)


def ftext(value: F) -> str:
    return f"{value.numerator}/{value.denominator}"


def profile_body(body: tuple[int, ...]) -> tuple[dict[str, object], ...]:
    root = C.A.profile_body(body)
    require(root["adaptive_k"] >= 20, f"root left high-K universe: {body}")
    root_good, root_r, root_h = C.A.T.CORE.good_norm(body)
    require(
        root_r == root["r"] and root_h == root["m"],
        f"root carrier changed: {body}",
    )
    open_rows: list[dict[str, object]] = []
    for rank in range(1, root["adaptive_k"] + 1):
        apex = root["top"][rank - 1][1]
        prefix = tuple(speed for _, speed in root["top"][:rank])
        excluded = set(prefix)
        carrier = C.H.R.subtract_local(root_good, apex)
        direct, direct_r, direct_h = C.H.T.CORE.good_norm(
            tuple(sorted((*body, apex)))
        )
        h = C.mass(carrier)
        require(
            carrier == direct and len(carrier) == direct_r and h == direct_h > 0,
            f"literal/direct carrier mismatch: {body}, rank {rank}",
        )
        ranked, tail_single, tail_first = C.sealed_ranking(carrier, excluded)
        top5 = tuple(ranked[:5])
        q3 = top5[2][0]
        q5 = top5[4][0]
        scalar_margin = h - sum((value for value, _ in top5), F(0))
        pair = C.pair_cap_from_ranking(carrier, ranked, tail_single)
        b2 = pair["cap"]
        q5_pair_margin = h - q5 - 2 * b2
        if scalar_margin > 0 or q5_pair_margin > 0:
            continue
        margins = (
            F(4, 7) * h - b2,
            F(5, 7) * h - q3 - b2,
            F(6, 7) * h - 2 * b2,
        )
        bits = tuple(margin > 0 for margin in margins)
        open_rows.append(
            {
                "body": body,
                "K": root["adaptive_k"],
                "rank": rank,
                "apex": apex,
                "prefix": prefix,
                "h": h,
                "r": direct_r,
                "q3": q3,
                "q5": q5,
                "b2": b2,
                "tail_first": tail_first,
                "scalar_margin": scalar_margin,
                "q5_pair_margin": q5_pair_margin,
                "margins": margins,
                "bits": bits,
                "paid": pair["paid"],
                "paid_digest": pair["paid_digest"],
            }
        )
    return tuple(open_rows)


def ledger_line(row: dict[str, object]) -> str:
    return (
        f"E={row['body']};K={row['K']};rank={row['rank']};a={row['apex']};"
        f"P={row['prefix']};h={ftext(row['h'])};r={row['r']};"
        f"q3={ftext(row['q3'])};q5={ftext(row['q5'])};B2={ftext(row['b2'])};"
        f"S={ftext(row['scalar_margin'])};Q5={ftext(row['q5_pair_margin'])};"
        f"H3={ftext(row['margins'][0])};H2={ftext(row['margins'][1])};"
        f"H1={ftext(row['margins'][2])};bits={row['bits']};"
        f"tail={row['tail_first']};paid={row['paid']};"
        f"paid_digest={row['paid_digest']}\n"
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--workers", type=int, default=min(os.cpu_count() or 1, 4)
    )
    args = parser.parse_args()
    require(args.workers >= 1, "workers must be positive")
    if args.workers == 1:
        local_rows = [profile_body(body) for body in HIGH_BODIES]
    else:
        context = mp.get_context("spawn")
        with context.Pool(args.workers) as pool:
            local_rows = list(
                pool.imap(profile_body, HIGH_BODIES, chunksize=2)
            )
    rows = [row for local in local_rows for row in local]
    require(len(rows) == 394, "high-K q5/top5 core changed")
    bit_histogram = tuple(sorted(Counter(row["bits"] for row in rows).items()))
    marginal_counts = tuple(
        sum(row["bits"][index] for row in rows) for index in range(3)
    )
    union_count = sum(any(row["bits"]) for row in rows)
    per_k: dict[int, Counter[tuple[bool, bool, bool]]] = defaultdict(Counter)
    for row in rows:
        per_k[row["K"]][row["bits"]] += 1
    per_k_histogram = tuple(
        (k, tuple(sorted(counter.items())))
        for k, counter in sorted(per_k.items())
    )
    minima = tuple(
        min(
            rows,
            key=lambda row: (
                row["margins"][index],
                row["body"],
                row["rank"],
            ),
        )
        for index in range(3)
    )
    digest = hashlib.sha256(b"LRC14/j6/highK-pair-flag-fork/v1\n")
    for row in rows:
        digest.update(ledger_line(row).encode())
    digest_text = digest.hexdigest()
    survivors = tuple(
        (
            row["body"],
            row["K"],
            row["rank"],
            row["apex"],
            tuple(ftext(margin) for margin in row["margins"]),
            ftext(row["h"]),
            ftext(row["q3"]),
            ftext(row["b2"]),
        )
        for row in rows
        if not any(row["bits"])
    )
    require(
        marginal_counts == EXPECTED_MARGINAL_COUNTS,
        "marginal eligibility counts changed",
    )
    require(
        bit_histogram == EXPECTED_BIT_HISTOGRAM,
        "eligibility histogram changed",
    )
    require(digest_text == EXPECTED_LEDGER_DIGEST, "flag-fork ledger changed")
    if EXPECTED_SURVIVORS is not None:
        require(survivors == EXPECTED_SURVIVORS, "survivor rows changed")

    print("LRC14 j6 K>=20 pair-cap flag-fork census")
    print(f"universe=roots:{len(HIGH_BODIES)},branches:{len(rows)}")
    print(f"marginal_counts_H3_H2_H1={marginal_counts}")
    print(f"eligibility_union={union_count};survivors={len(rows)-union_count}")
    print(f"bit_histogram={bit_histogram}")
    print(f"per_K_histogram={per_k_histogram}")
    print(f"survivor_rows={survivors}")
    print(
        "minimum_margins="
        + repr(
            tuple(
                (
                    ftext(row["margins"][index]),
                    row["body"],
                    row["rank"],
                    row["apex"],
                )
                for index, row in enumerate(minima)
            )
        )
    )
    print(f"pair_paid={sum(row['paid'] for row in rows)}")
    print(f"canonical_ledger_sha256={digest_text}")
    print("scope=eligibility only;flag children remain proof obligations")
    print("all_exact_controls=PASS")


if __name__ == "__main__":
    main()
