#!/usr/bin/env python3
"""Exact THM-2941 closure for the unique 1770..1779 scalar survivor."""

from __future__ import annotations

import argparse
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "04-computation" / "lrc14_j7_k2_1780_1789_exact_descent_closure_thm2941.py"
CASE = (1776, (1, 4, 8, 10, 12, 14), "PROJECT")
EXPECTED_SOURCE_SHA256 = "b57c68539a5da6fbe456e29f07bc7ff55778702e70990bcfa54e7409ad56485b"
EXPECTED_PROFILE_SHA256 = "1644442c7a407e64c45631718442e487c8272b5aab75f9ef8e18a42b06cbec9d"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def load():
    spec = spec_from_file_location("z1776_inherited_exact_worker", SOURCE)
    if spec is None or spec.loader is None:
        raise RuntimeError("cannot load inherited exact worker")
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    source_hash = sha256(SOURCE.read_bytes()).hexdigest()
    require(source_hash == EXPECTED_SOURCE_SHA256, "inherited closure source changed")
    row = load().exact_worker(CASE)
    first, body = row[:2]
    counts = row[12]
    cells, packets, kills, margin, prefix, minimum, direct, state_hash = row[14:22]
    closed = counts[3] == 0 or (packets == kills and margin is not None and margin > 0)
    require(closed, "z1776 row survived")
    require(counts == (77, 29, 48, 0, 0), "z1776 status census changed")
    payload = (CASE, row)
    profile_hash = sha256(repr(payload).encode()).hexdigest()
    require(profile_hash == EXPECTED_PROFILE_SHA256, "z1776 profile digest changed")
    lines = (
        "LRC14 projected k=2 z1776 exact descent closure",
        f"source_sha256={source_hash}",
        f"z1={first};E={','.join(map(str, body))};counts={counts};cells={cells};"
        f"packets={packets};kills={kills};min_margin={margin};max_prefix={prefix};"
        f"minimum={minimum};direct={direct};state_sha256={state_hash}",
        "closure=1/1",
        f"profile_sha256={profile_hash}",
        "all_exact_controls=PASS",
    )
    text = "\n".join(lines) + "\n"
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(text, encoding="utf-8", newline="\n")
    print(text, end="")


if __name__ == "__main__":
    main()
