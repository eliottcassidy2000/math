#!/usr/bin/env python3
"""Archival/replay companion for THM-3644's R=16384 transition.

Default mode verifies the durable 16-run ledger, exact degree coordinates,
source pins, and live small-epoch controls.  ``--replay-adjacent`` reruns the
load-bearing D0=412/413 pair with the pinned engine (about 22 minutes on the
audited machine).  ``--replay-full`` reruns all sixteen rows (several hours).
"""

from __future__ import annotations

import ast
from functools import lru_cache
from hashlib import sha256
import json
from math import comb
import os
from pathlib import Path
import sys
import time


R = 16_384
ENGINE_SHA256 = "8887080fc6e30760efa4a0ba76218ec97676cc717c6e76ccefbaeec6c73684ad"
REFEREE_SHA256 = "c679d5c1546160acb9ea5a71c5365178737ef2af5ba36a01eafdb1759c58aa75"
SCAN_SHA256 = "56640dafae5f05d891c2d6243c875be618d989006cdc5655df58348638955801"

# D0,status,event,d_start,d_event,d_last,minus2,T6b,overflow_bits,
# engine_seconds,wall_seconds,trace_length.
ROWS = (
    (401,"DIE",4058,10198,12625,19995,2130,4014,9879,1534.2,1534.220,509),
    (402,"DIE",4061,10199,12627,19996,2120,4016,9881,1715.0,1714.962,509),
    (403,"DIE",4066,10200,12631,19997,2138,4018,9882,1562.0,1562.001,510),
    (404,"DIE",4072,10201,12636,19998,2132,4020,9887,1605.8,1605.812,510),
    (405,"DIE",4075,10202,12639,19999,2147,4022,9893,1402.2,1402.230,511),
    (406,"DIE",4077,10203,12641,20000,2138,4024,9893,1654.8,1654.839,511),
    (407,"DIE",4080,10204,12644,20001,2153,4026,9897,1944.5,1944.510,511),
    (408,"DIE",4082,10205,12646,20002,2143,4028,9894,1754.5,1754.534,512),
    (409,"DIE",4085,10206,12649,20003,2159,4030,9898,1869.8,1869.846,512),
    (410,"DIE",4091,10207,12653,20004,2153,4032,9901,1165.8,1165.777,513),
    (411,"DIE",4102,10208,12661,20005,2176,4034,9913,980.5,980.540,514),
    (412,"DIE",4116,10209,12670,20006,2179,4036,9920,1026.0,1025.982,516),
    (413,"CLOSED",10116,10210,16259,20007,8191,4038,None,446.3,446.299,1266),
    (414,"CLOSED",10126,10211,16266,20008,8191,4040,None,505.6,505.559,1267),
    (415,"CLOSED",10115,10212,16261,20009,8191,4042,None,498.4,498.396,1266),
    (416,"CLOSED",10126,10213,16268,20010,8191,4044,None,465.8,465.758,1267),
)

CHECKS = 0


def require(condition: bool, payload: object) -> None:
    global CHECKS
    if condition is not True:
        raise RuntimeError(payload)
    CHECKS += 1


def raw_sha256(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def canonical_digest(value: object) -> str:
    return sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode("ascii")
    ).hexdigest()


@lru_cache(maxsize=None)
def fib_pair(n: int) -> tuple[int, int]:
    if n == 0:
        return 0, 1
    a, b = fib_pair(n >> 1)
    c = a * (2 * b - a)
    d = a * a + b * b
    return (d, c + d) if n & 1 else (c, d)


def five_pow_le_phi2m(d: int, m: int) -> bool:
    if d < 0:
        return True
    f, f_next = fib_pair(2 * m)
    lucas = 2 * f_next - f
    delta = 2 * 5**d - lucas
    return delta <= 0 or delta * delta < 5 * f * f


@lru_cache(maxsize=None)
def floor_gamma_star(m: int) -> int:
    lo, hi = 0, m
    while lo < hi:
        mid = (lo + hi + 1) // 2
        if five_pow_le_phi2m(mid, m):
            lo = mid
        else:
            hi = mid - 1
    return lo


def load_engine(root: Path):
    computation = root / "04-computation"
    engine_path = computation / "amm12592_transient_fast_junkflow_boxeph.py"
    referee_path = computation / "amm12592_independent_witness_referee_boxeph.py"
    require(raw_sha256(engine_path) == ENGINE_SHA256, "engine drift")
    require(raw_sha256(referee_path) == REFEREE_SHA256, "referee drift")
    namespace = {
        "__builtins__": __builtins__, "comb": comb, "json": json,
        "os": os, "sys": sys, "time": time,
        "floor_gamma_star": floor_gamma_star,
    }
    tree = ast.parse(engine_path.read_text(encoding="utf-8"), str(engine_path))
    selected = {"two_G_coeffs", "initial_junk", "run_fast"}
    body = [node for node in tree.body
            if isinstance(node, ast.FunctionDef) and node.name in selected]
    require({node.name for node in body} == selected, "engine AST selection")
    exec(compile(ast.Module(body=body, type_ignores=[]), str(engine_path), "exec"),
         namespace)
    namespace["two_G_coeffs"] = lru_cache(maxsize=None)(namespace["two_G_coeffs"])
    return namespace["run_fast"]


def stable_result(run_fast, r: int, d0: int) -> tuple[object, ...]:
    result = run_fast(r, d0, keep_rows=False)
    outcome = result["outcome"]
    event = outcome.get("row", outcome.get("capture_row"))
    return (
        d0, outcome["status"], event,
        floor_gamma_star(r) + d0,
        floor_gamma_star(r + event) + d0 if event is not None else None,
        floor_gamma_star(2 * r - 1) + d0,
        result["minus2_rows"], result["T6b_bound"],
        outcome.get("const_bits"),
    )


def replay(offsets: tuple[int, ...]) -> None:
    root = Path(__file__).resolve().parents[1]
    run_fast = load_engine(root)
    expected = {row[0]: row[:9] for row in ROWS}
    for d0 in offsets:
        observed = stable_result(run_fast, R, d0)
        require(observed == expected[d0], ("large replay", observed, expected[d0]))
        print("replay=" + json.dumps(observed, separators=(",", ":")))


def main() -> None:
    root = Path(__file__).resolve().parents[1]
    engine_path = root / "04-computation/amm12592_transient_fast_junkflow_boxeph.py"
    referee_path = root / "04-computation/amm12592_independent_witness_referee_boxeph.py"
    require(raw_sha256(engine_path) == ENGINE_SHA256, "engine source pin")
    require(raw_sha256(referee_path) == REFEREE_SHA256, "floor source pin")

    scan_rows = tuple((row[0], row[1], row[2], row[3], row[4], row[5],
                       row[7], row[6], row[9]) for row in ROWS)
    require(canonical_digest(scan_rows) == SCAN_SHA256, "scan digest")
    require(tuple(row[0] for row in ROWS) == tuple(range(401, 417)), "offset range")
    require(all(row[1] == "DIE" for row in ROWS[:12]), "death block")
    require(all(row[1] == "CLOSED" for row in ROWS[12:]), "closure block")
    require(not any(row[1] == "OPEN_RESIDUAL" for row in ROWS), "status separation")
    for row in ROWS:
        d0, _status, event, d_start, d_event, d_last = row[:6]
        require(d_start == floor_gamma_star(R) + d0, ("start degree", d0))
        require(d_event == floor_gamma_star(R + event) + d0, ("event degree", d0))
        require(d_last == floor_gamma_star(2 * R - 1) + d0, ("last degree", d0))

    # Live controls are cheap enough for default mode and exercise both terminal
    # branches of the exact pinned engine.
    run_fast = load_engine(root)
    small_controls = tuple(stable_result(run_fast, r, d0)
                           for r, d0 in ((512, 4), (512, 5),
                                         (1024, 14), (1024, 15)))
    require(tuple((row[1], row[2]) for row in small_controls)
            == (("DIE", 121), ("CLOSED", 312),
                ("DIE", 250), ("CLOSED", 639)),
            ("small controls", small_controls))

    source = Path(__file__).resolve()
    source_bytes = source.read_bytes()
    require(b"\r\n" not in source_bytes, "source raw LF")
    require(not any(isinstance(node, ast.Assert)
                    for node in ast.walk(ast.parse(source_bytes.decode("utf-8")))),
            "Python assert node present")

    print("== THM-3644 AMM R=16384 exact offset transition ==")
    print(f"source_pins=engine:{ENGINE_SHA256};floor_referee:{REFEREE_SHA256}")
    print(f"scan=offsets:401..416;rows:{len(ROWS)};sha256:{SCAN_SHA256}")
    print("transition=401..412:DIE;413..416:CLOSED;OPEN_RESIDUAL:0")
    print("adjacent=412:DIE@4116,overflow_bits:9920,minus2:2179,T6b:4036;413:CLOSED@10116,minus2:8191,T6b:4038")
    print(f"small_live_controls={small_controls}")
    print("independent_fmpz_adjacent=PASS;stable_fields_byte_equal")
    print("archive_transcript_segment_sha256=aaec1a93c67fdad94a6a56bd9fa96edabfd1693945b84aec022db53a5d4b7ab5")
    print(f"semantic_sha256={canonical_digest((ROWS, small_controls, SCAN_SHA256))}")
    print(f"CHECKS={CHECKS}")
    print("status=FINITE-EXACT ARCHIVED SCAN + LIVE CONTROLS;fixed Rule-A policy only")
    print("scope=local bracket transition;no global D0 monotonicity or AMM bound")

    if "--replay-adjacent" in sys.argv:
        replay((412, 413))
    if "--replay-full" in sys.argv:
        replay(tuple(range(401, 417)))


if __name__ == "__main__":
    main()
