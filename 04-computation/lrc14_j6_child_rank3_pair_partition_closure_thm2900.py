#!/usr/bin/env python3
"""Locked THM-2900 census of the THM-2897 q3+B2 child certificate."""

from __future__ import annotations

import hashlib
import importlib.util
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
ENGINE = (
    ROOT
    / "04-computation/lrc14_j6_suffix_parity_flag_closure_thm2895_independent_audit.py"
)
ENGINE_SHA = "9023a4042dc8def3f8781668e721549972fb1458d07f2fab89bf7ac3a6f745cc"
EXPECTED_SHARP = (
    "1720879997/753665913000",
    (2, 8, 9, 10, 11, 13, 14),
    1,
    19,
    (37, 125),
    "1113566829/28324045750",
    "1335801387779/17334315999000",
)
EXPECTED_LEDGER_DIGEST: str | None = (
    "60ad47208695c5aa22ec79ca428eccfdd5e214cd7cf4a20a85629aa5023a6fa2"
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def load():
    require(
        hashlib.sha256(ENGINE.read_bytes()).hexdigest() == ENGINE_SHA,
        "independent engine changed",
    )
    spec = importlib.util.spec_from_file_location("thm2897_child_engine", ENGINE)
    require(spec is not None and spec.loader is not None, "cannot load engine")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


E = load()


def ftext(value: F) -> str:
    return f"{value.numerator}/{value.denominator}"


def main() -> None:
    roots = [E.root_profile(body) for body in E.ROOTS]
    branches = [
        E.suffix_profile(root, rank)
        for root in roots
        for rank in range(1, root["K"] + 1)
    ]
    open_branches = [
        E.h4_profile(branch)
        for branch in branches
        if not branch["closed"]
    ]
    rows = [
        (branch, E.pair_profile(branch, pair))
        for branch in open_branches
        for pair in combinations(branch["H"], 2)
    ]
    margins = [
        (pair["m"] - pair["top3"][2][0] - pair["B2"], branch, pair)
        for branch, pair in rows
    ]
    sharp = min(
        margins,
        key=lambda item: (
            item[0],
            item[1]["body"],
            item[1]["rank"],
            item[2]["hpair"],
        ),
    )
    require(len(rows) == 784, "pair universe changed")
    require(all(margin > 0 for margin, _, _ in margins), "q3+B2 failure")
    sharp_result = (
        ftext(sharp[0]),
        sharp[1]["body"],
        sharp[1]["rank"],
        sharp[1]["apex"],
        sharp[2]["hpair"],
        ftext(sharp[2]["top3"][2][0]),
        ftext(sharp[2]["B2"]),
    )
    require(sharp_result == EXPECTED_SHARP, "sharp q3+B2 row changed")
    ledger = hashlib.sha256(b"LRC14/j6/child-q3-plus-B2/v1\n")
    for branch, pair in rows:
        margin = pair["m"] - pair["top3"][2][0] - pair["B2"]
        ledger.update(
            (
                f"E={branch['body']};rank={branch['rank']};"
                f"a={branch['apex']};L={pair['hpair']};"
                f"h={ftext(pair['m'])};"
                f"q3={ftext(pair['top3'][2][0])};"
                f"B2={ftext(pair['B2'])};"
                f"margin={ftext(margin)}\n"
            ).encode()
        )
    ledger_digest = ledger.hexdigest()
    if EXPECTED_LEDGER_DIGEST is not None:
        require(
            ledger_digest == EXPECTED_LEDGER_DIGEST,
            "child q3+B2 ledger changed",
        )
    print("THM-2900 CHILD q3+B2 INDEPENDENT AUDIT")
    print(f"pairs={len(rows)};closed={sum(margin > 0 for margin,_,_ in margins)}")
    print(
        f"minimum_margin={ftext(sharp[0])};"
        f"E={sharp[1]['body']};rank={sharp[1]['rank']};"
        f"apex={sharp[1]['apex']};H4pair={sharp[2]['hpair']};"
        f"q3={ftext(sharp[2]['top3'][2][0])};"
        f"B2={ftext(sharp[2]['B2'])}"
    )
    print("recursive_H2_rows_needed=0")
    print(f"canonical_ledger_sha256={ledger_digest}")
    print(
        "mode="
        + ("DISCOVERY" if EXPECTED_LEDGER_DIGEST is None else "LOCKED")
    )
    print("all_independent_controls=PASS")


if __name__ == "__main__":
    main()
