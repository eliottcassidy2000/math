#!/usr/bin/env python3
"""Exact THM-2941 closure for the twelve scalar survivors on 1780..1789."""

from __future__ import annotations

import argparse
import concurrent.futures as cf
import os
from fractions import Fraction as F
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from math import lcm
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
EXACT_PATH = ROOT / "04-computation" / "lrc14_j7_k2_exact_descent_1800_1824_closure_thm2941.py"
HIGH_PATH = ROOT / "04-computation" / "lrc14_j7_k2_high_wall_descent_1800_1810_closure_thm2941.py"
EXPECTED_EXACT_SHA256 = "1bc6674fbb9b6f4c8979c229c164d267e60911ed582fb3184813d45c21da2adf"
EXPECTED_HIGH_SHA256 = "c12c17297aa8a96cbdd7d9d529838c776b160ace86b92d0030f5df447fe6877b"

EXACT_CASES = (
    (1780, (1, 4, 8, 10, 12, 14), "PROJECT"),
    (1784, (1, 4, 8, 10, 12, 14), "PROJECT"),
    (1788, (1, 4, 8, 10, 12, 14), "PROJECT"),
    (1784, (1, 6, 8, 10, 12, 14), "PROJECT"),
    (1780, (2, 4, 8, 10, 12, 14), "PROJECT"),
    (1784, (2, 4, 8, 10, 12, 14), "PROJECT"),
    (1780, (2, 6, 8, 10, 12, 14), "PROJECT"),
    (1784, (2, 6, 8, 10, 12, 14), "PROJECT"),
    (1780, (4, 6, 8, 10, 12, 14), "PROJECT"),
    (1784, (4, 6, 8, 10, 12, 14), "PROJECT"),
)
HIGH_CASES = (
    (1780, (1, 8, 10, 12, 13, 14)),
    (1784, (2, 8, 9, 10, 12, 14)),
)
EXPECTED_PROFILE_SHA256 = "1196675a616cb54d49343773152376302570a6e1601cbadf1f56c0bcafff37de"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def load(name: str, path: Path):
    spec = spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(("cannot load", path))
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def exact_worker(item):
    module = load("k2_1780_exact_worker", EXACT_PATH)
    first, body, mode = item
    carrier = module.U.suffix.A.carrier_for(body)
    require(module.P.A.carrier_for(body) == carrier, "carrier engines disagree")
    h = F(sum(right - left for left, right in carrier), module.U.suffix.A.RULER)
    lower = h * module.U.suffix.ETAS[2]
    L = 14 * lcm(*body)
    require(module.P.A.RULER % L == 0, "unexpected ruler")
    (
        amplitudes,
        ray_digest,
        divisor_count,
        trials,
        first_delta,
        first_d,
        scalar,
        crude_kills,
        status_kills,
        states,
        stage_digests,
    ) = module.ray_and_status(first, body, carrier, h, lower, L)
    projected = module.projected_packets(
        first, body, carrier, h, lower, L, amplitudes, first_delta, states
    )
    counts = (
        len(scalar),
        len(crude_kills),
        len(status_kills),
        len(states),
        projected[1],
    )
    return (
        first,
        body,
        mode,
        h,
        len(carrier),
        L,
        lower,
        first_delta,
        first_d,
        ray_digest,
        divisor_count,
        trials,
        counts,
        stage_digests,
        *projected[:-1],
    )


def high_worker(item):
    module = load("k2_1780_high_worker", HIGH_PATH)
    module.EXPECTED_GAPS = {item: None}
    strict_require = module.require

    def allow_surviving_gap(condition, message):
        if (
            not condition
            and isinstance(message, tuple)
            and len(message) >= 2
            and message[1] == "forced-high scalar row survived"
        ):
            return
        strict_require(condition, message)

    module.require = allow_surviving_gap
    return module.profile(item)


def exact_summary(row):
    first, body = row[:2]
    counts = row[12]
    cells, packets, kills, margin, prefix, minimum, direct, state_hash = row[14:22]
    closed = counts[3] == 0 or (packets == kills and margin is not None and margin > 0)
    return closed, (
        f"EXACT;z1={first};E={','.join(map(str, body))};counts={counts};cells={cells};"
        f"packets={packets};kills={kills};min_margin={margin};max_prefix={prefix};"
        f"minimum={minimum};direct={direct};state_sha256={state_hash};closed={int(closed)}"
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--hash-seed", type=int, default=0)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    require(args.workers >= 1, "worker count must be positive")
    require(args.hash_seed >= 0, "hash seed must be nonnegative")
    os.environ["PYTHONHASHSEED"] = str(args.hash_seed)
    require(sha256(EXACT_PATH.read_bytes()).hexdigest() == EXPECTED_EXACT_SHA256,
            "exact descent source changed")
    require(sha256(HIGH_PATH.read_bytes()).hexdigest() == EXPECTED_HIGH_SHA256,
            "high-wall source changed")
    with cf.ProcessPoolExecutor(max_workers=args.workers) as pool:
        exact = tuple(pool.map(exact_worker, EXACT_CASES))
        high = tuple(pool.map(high_worker, HIGH_CASES))
        high_project_items = tuple(
            (row[0], row[1], "PROJECT") for row in high if row[-1] >= 0
        )
        high_project = tuple(pool.map(exact_worker, high_project_items))
    exact = tuple(sorted(exact, key=lambda row: (row[0], row[1])))
    high = tuple(sorted(high, key=lambda row: (row[0], row[1])))
    high_project = tuple(sorted(high_project, key=lambda row: (row[0], row[1])))
    lines = [
        "LRC14 projected k=2 exact descent closure 1780..1789",
        f"exact_cases={len(exact)};high_cases={len(high)};high_project={len(high_project)}",
    ]
    closed = 0
    for row in exact:
        verdict, line = exact_summary(row)
        closed += int(verdict)
        lines.append(line)
    for row in high:
        gap = row[-1]
        verdict = gap < 0
        closed += int(verdict)
        lines.append(
            f"HIGH;z1={row[0]};E={','.join(map(str, row[1]))};L={row[4]};"
            f"high_floor={row[5]};branch={row[15]};upper={row[-2]};gap={gap};closed={int(verdict)}"
        )
    for row in high_project:
        verdict, line = exact_summary(row)
        closed += int(verdict)
        lines.append("HIGH-PROJECT;" + line.split(";", 1)[1])
    payload = (exact, high, high_project)
    profile_hash = sha256(repr(payload).encode()).hexdigest()
    require(len(exact) == 10, "ordinary exact-case census changed")
    require(len(high) == 2, "forced-high case census changed")
    require(len(high_project) == 0, "unexpected surviving high-project case")
    require(closed == 12, "1780..1789 closure changed")
    require(
        profile_hash == EXPECTED_PROFILE_SHA256,
        f"closure profile digest changed: expected {EXPECTED_PROFILE_SHA256}; got {profile_hash}",
    )
    lines.extend((
        f"closure={closed}/12",
        f"profile_sha256={profile_hash}",
        "all_exact_controls=PASS",
    ))
    text = "\n".join(lines) + "\n"
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(text, encoding="utf-8", newline="\n")
    print(text, end="")


if __name__ == "__main__":
    main()
