#!/usr/bin/env python3
"""Exact projected k=2 closure of the two ordinary rows below 1780.

The global 1750..1799 scalar band and the independent z1=1790 and z1=1788
z1=1784, and z1=1780 verifiers leave two ordinary exact-suffix rows in this
range.  The remaining HIGH-TAIL row is closed by the companion exact-ray
scalar verifier.  This file closes exactly those two ordinary rows.

The inherited data-driven pipeline is

  residue rays -> denominator multisets -> all-divisor capacity
  -> exact common-K5 status -> finite scalar-slack literal packets
  -> lossless projected residual.

Every residual packet has projected mass at least 25/91 and so cannot fit in
the two remaining aligned danger combs.  The finite label enumeration is
derived from positive exact slack on each hyperbolic residue ray; no label
horizon is imposed.
"""

from __future__ import annotations

import argparse
import multiprocessing as mp
from fractions import Fraction as F
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from math import lcm
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
BASE = (
    ROOT
    / "04-computation"
    / "lrc14_j7_k2_exact_descent_1800_1824_closure_thm2941.py"
)
BAND = ROOT / "05-knowledge/results/lrc14_j7_k2_scalar_band_1750_1799_thm2941.out"
Z1790 = ROOT / "05-knowledge/results/lrc14_j7_k2_z1790_exact_descent_closure_thm2941.out"
Z1788 = ROOT / "05-knowledge/results/lrc14_j7_k2_z1788_ray_status_closure_thm2941.out"
Z1784 = ROOT / "05-knowledge/results/lrc14_j7_k2_z1784_all_label_closure_thm2941.out"
Z1780 = ROOT / "05-knowledge/results/lrc14_j7_k2_z1780_hybrid_closure_thm2941.out"
HIGH = ROOT / "05-knowledge/results/lrc14_j7_k2_high_wall_descent_1768_1790_closure_thm2941.out"
OUTPUT_PATH = (
    ROOT
    / "05-knowledge/results/lrc14_j7_k2_exact_descent_1758_1790_closure_thm2941.out"
)
EXPECTED_BASE_SHA256 = (
    "e9242e697efb0ba452ca0c51ee65da4d02b427557ef627414286661cfea1cf79"
)
EXPECTED_BAND_SHA256 = (
    "2ce806d361d7eb97d9ae2d23e438898c8e1f895a89501c9a1847e51f61ca8009"
)
EXPECTED_Z1790_SHA256 = (
    "b03b46a6c438773ef1c433c435828a5426e18c998999cbee543541904b85f20b"
)
EXPECTED_Z1788_SHA256 = (
    "9656de4784eb936919113e8a4151f60514c21089a3f110994b94654ea0c070ba"
)
EXPECTED_Z1784_SHA256 = (
    "318f0908b2c84d68a5ef75884316f0ebc600d6aeed28d640f7489d6f179fb6bf"
)
EXPECTED_Z1780_SHA256 = (
    "afa96047c84a544b77142e34dfd25ad3923bace4d54f3f5ce47b0910f99ee27b"
)
EXPECTED_HIGH_SHA256 = (
    "ce64016d49a8cc93a8d2f532091d7253bf62fd9f3fa7d81d6141f05740d8886d"
)

HOSTILE = (1, 4, 8, 10, 12, 14)
SLICES = (
    (1776, HOSTILE, "STATUS"),
    (1758, HOSTILE, "STATUS"),
)
HIGH_SCALAR_KEYS = (
    (1768, (1, 8, 10, 12, 13, 14)),
)
ITEMS = SLICES

EXPECTED_COUNTS = {
    (1776, HOSTILE): (77, 29, 48, 0, 0),
    (1758, HOSTILE): (2, 0, 2, 0, 0),
}
EXPECTED_PROFILE_SHA256 = (
    "fc2bfb402a9692f26973091eb88e890937aa2e4bfa45a7a4bcae7f76c9a09af0"
)
EXPECTED_SEMANTIC_SHA256 = (
    "815b87faba8b0216ef01eefcc7eee54a66eaf2066448b2a54802fa4a96a301b3"
)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def file_sha256(path):
    return sha256(path.read_bytes()).hexdigest()


def ftext(value):
    if value is None:
        return "NONE"
    return f"{value.numerator}/{value.denominator}"


require(file_sha256(BASE) == EXPECTED_BASE_SHA256, "exact descent engine changed")
require(file_sha256(BAND) == EXPECTED_BAND_SHA256, "1750..1799 atlas changed")
require(file_sha256(Z1790) == EXPECTED_Z1790_SHA256, "z1790 closure changed")
require(file_sha256(Z1788) == EXPECTED_Z1788_SHA256, "z1788 closure changed")
require(file_sha256(Z1784) == EXPECTED_Z1784_SHA256, "z1784 closure changed")
require(file_sha256(Z1780) == EXPECTED_Z1780_SHA256, "z1780 closure changed")
if EXPECTED_HIGH_SHA256 is not None:
    require(file_sha256(HIGH) == EXPECTED_HIGH_SHA256, "remaining high closure changed")
SPEC = spec_from_file_location("k2_exact_1758_base", BASE)
require(SPEC is not None and SPEC.loader is not None, "cannot load exact descent engine")
D = module_from_spec(SPEC)
SPEC.loader.exec_module(D)
D.EXPECTED_COUNTS = {item[:2]: None for item in SLICES}


def atlas_partition():
    keys = []
    high = []
    for line in BAND.read_text().splitlines():
        if not line.startswith("SURVIVOR;"):
            continue
        fields = dict(field.split("=", 1) for field in line.split(";")[1:] if "=" in field)
        key = (int(fields["z1"]), tuple(map(int, fields["E"].split(","))))
        keys.append(key)
        if "HIGH-TAIL" in line:
            high.append(key)
    return tuple(sorted(keys)), tuple(sorted(high))


ATLAS_KEYS, ATLAS_HIGH_KEYS = atlas_partition()
require(
    tuple(key for key in ATLAS_KEYS if 1751 <= key[0] <= 1779)
    == tuple(sorted((*tuple(item[:2] for item in SLICES), *HIGH_SCALAR_KEYS))),
    "exact/high cases do not partition the 1751..1779 atlas",
)
require(
    tuple(key for key in ATLAS_HIGH_KEYS if 1751 <= key[0] <= 1779)
    == tuple(sorted(HIGH_SCALAR_KEYS)),
    "remaining HIGH-TAIL classification changed",
)


def all_label_profile(item):
    first, body, mode = item
    carrier = D.U.suffix.A.carrier_for(body)
    require(D.P.A.carrier_for(body) == carrier, (first, body, "carrier mismatch"))
    h = F(sum(right - left for left, right in carrier), D.U.suffix.A.RULER)
    lower = h * D.U.suffix.ETAS[2]
    L = 14 * lcm(*body)
    require(D.P.A.RULER % L == 0, (first, body, "ruler changed", L))
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
    ) = D.ray_and_status(first, body, carrier, h, lower, L)
    projected = D.projected_packets(
        first, body, carrier, h, lower, L, amplitudes, first_delta, states
    )
    counts = (
        len(scalar),
        len(crude_kills),
        len(status_kills),
        len(states),
        projected[1],
    )
    expected = EXPECTED_COUNTS[(first, body)]
    if expected is not None:
        require(counts == expected, (first, body, "counts changed", counts))
    require(projected[1] == projected[2], (first, body, "projected packet survived"))
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


def profile(item):
    return all_label_profile(item)


def render(profiles):
    require(tuple((row[0], row[1], row[2]) for row in profiles) == ITEMS, "item universe changed")
    require(all(row[12][3] == 0 or row[15] == row[16] for row in profiles), "closure mismatch")
    totals = tuple(sum(row[12][index] for row in profiles) for index in range(5))
    positive_margins = tuple(row[17] for row in profiles if row[17] is not None)
    require(totals[4] == sum(row[16] for row in profiles), "a projected packet survived")
    require(all(margin > 0 for margin in positive_margins), "nonpositive projected margin")
    actual_counts = tuple((row[0], row[1], row[12]) for row in profiles)
    for first, body, counts in actual_counts:
        expected = EXPECTED_COUNTS[(first, body)]
        if expected is not None:
            require(counts == expected, (first, body, "stage counts changed", counts))
    profile_hash = sha256(repr(tuple(profiles)).encode()).hexdigest()
    if EXPECTED_PROFILE_SHA256 is not None:
        require(profile_hash == EXPECTED_PROFILE_SHA256, "profile digest changed")
    semantic_payload = (
        ITEMS,
        HIGH_SCALAR_KEYS,
        D.QUANTIFIER,
        D.ALIGNED_TWO_UNION_CAP,
        actual_counts,
        totals,
        positive_margins,
        profile_hash,
    )
    semantic_hash = sha256(repr(semantic_payload).encode()).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(semantic_hash == EXPECTED_SEMANTIC_SHA256, "semantic digest changed")
    lines = [
        "LRC14 projected k=2 exact descent closure 1758..1776",
        f"exact_descent_engine_sha256={file_sha256(BASE)}",
        f"scalar_band_sha256={file_sha256(BAND)}",
        f"z1790_closure_sha256={file_sha256(Z1790)}",
        f"z1788_closure_sha256={file_sha256(Z1788)}",
        f"z1784_closure_sha256={file_sha256(Z1784)}",
        f"z1780_closure_sha256={file_sha256(Z1780)}",
        f"high_closure_sha256={file_sha256(HIGH)}",
        "scope=two remaining exact-suffix atlas rows;all later distinct nonaligned labels;no label horizon",
        f"totals=scalar:{totals[0]};crude:{totals[1]};status:{totals[2]};residual:{totals[3]};packets:{totals[4]};kills:{sum(row[16] for row in profiles)}",
    ]
    for row in profiles:
        lines.append(
            f"CASE;z1={row[0]};E={','.join(map(str, row[1]))};mode={row[2]};"
            f"L={row[5]};counts={row[12]};stage_digests={row[13]};"
            f"projected_margin={ftext(row[17])};max_cells={row[18]};"
            f"minimum_row={row[19]};direct_control_mass={ftext(row[20])};"
            f"state_digest={row[21]};conclusion=EMPTY"
        )
    lines.extend(
        (
            "ordinary_rows=2;survivors=0",
            f"profile_sha256={profile_hash}",
            f"semantic_sha256={semantic_hash}",
            "all_exact_controls=PASS",
        )
    )
    return "\n".join(lines) + "\n"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=min(len(ITEMS), mp.cpu_count() or 1))
    parser.add_argument("--output", type=Path, default=OUTPUT_PATH)
    args = parser.parse_args()
    require(args.workers >= 1, "workers must be positive")
    if args.workers == 1:
        profiles = [profile(item) for item in ITEMS]
    else:
        with mp.get_context("spawn").Pool(args.workers) as pool:
            profiles = list(pool.imap(profile, ITEMS, chunksize=1))
    order = {item: index for index, item in enumerate(ITEMS)}
    profiles.sort(key=lambda row: order[(row[0], row[1], row[2])])
    output = render(tuple(profiles))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(output, encoding="utf-8", newline="\n")
    print(output, end="")


if __name__ == "__main__":
    main()
