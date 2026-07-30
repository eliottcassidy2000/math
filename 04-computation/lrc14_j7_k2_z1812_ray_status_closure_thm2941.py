#!/usr/bin/env python3
"""Exact all-nonaligned-label closure of the projected k=2 row at z1=1812.

The guarded global atlas on ``1810 <= z1 <= 1835`` has one row at z1=1812,
on ``E=(1,4,8,10,12,14)`` with ``L=11760``.  The already-proved z1=1824
closure and the same atlas give the inherited cap z1<=1812.  This verifier
therefore attacks the new frontier.

It reconstructs, over every four-tuple of later distinct nonaligned labels,
the exact residue-ray scalar maximum.  The first denominator is 980.  Of all
``C(59+3,4)=557845`` denominator multisets, eleven meet the scalar lower bound.
Exact all-divisor fibre capacity kills four, and exact common K5 status Farkas
certificates kill seven.  Nothing survives, so the atlas's next occupied
height gives the sharper cap z1<=1810.  No finite label horizon or projected
residual is used.
"""

from __future__ import annotations

import argparse
import os
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from math import gcd, lcm
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PARENT = (
    ROOT
    / "04-computation"
    / "lrc14_j7_k2_z1824_ray_status_closure_thm2941.py"
)
PARENT_OUTPUT = (
    ROOT
    / "05-knowledge"
    / "results"
    / "lrc14_j7_k2_z1824_ray_status_closure_thm2941.out"
)
BAND_SOURCE = (
    ROOT
    / "04-computation"
    / "lrc14_j7_k2_scalar_band_1810_1835_thm2941.py"
)
BAND_OUTPUT = (
    ROOT
    / "05-knowledge"
    / "results"
    / "lrc14_j7_k2_scalar_band_1810_1835_thm2941.out"
)
OUTPUT_PATH = (
    ROOT
    / "05-knowledge"
    / "results"
    / "lrc14_j7_k2_z1812_ray_status_closure_thm2941.out"
)

EXPECTED_PARENT_SHA256 = (
    "6d825a99d57c8febedfe212fd4076495204ebefe71e3cb47e0d86cb38de238ab"
)
EXPECTED_PARENT_OUTPUT_SHA256 = (
    "06869b6c718e9f63799904bb5e1fb19c88e8e356f9272e10d3eb14175ac3cebd"
)
EXPECTED_BAND_SOURCE_SHA256 = (
    "058c43d67d0bb110993ec687877edba4f5a7ad472a81455b0a5276b20db7680d"
)
EXPECTED_BAND_OUTPUT_SHA256 = (
    "76f08dc5b70c98dd7c8fa958f5598f5c50e7cb1df26d0be10531ba8185796952"
)
EXPECTED_PROFILE_SHA256 = (
    "e0b8768bc9d0749634712f1b92def48652271f633efbcad573923cef354b4b3a"
)
EXPECTED_SEMANTIC_SHA256 = (
    "251514d6d7f8bea9d993b0140bbdea454fa7b649fa18ee085a7c85cf37635b48"
)

BODY = (1, 4, 8, 10, 12, 14)
FIRST = 1812
NEXT_CAP = 1810
QUANTIFIER = "four later distinct nonaligned labels"
EXPECTED_HEIGHTS = (1810,) * 8 + (1812, 1824)
EXPECTED_PARENT_CAP = (
    "consequence=the proved projected k=2 first-drift cap is z1<=1812"
)
EXPECTED_BAND_ROW = (
    "SURVIVOR;E=1,4,8,10,12,14;h=1049/2940;r=34;L=11760;"
    "largest_floor=1020;z1=1812;delta1=377/1035860;"
    "suffix=EXACT:1836:3583/3148740,EXACT:2060:1867/2119740,"
    "EXACT:2172:263/310415,EXACT:2142:821/1049580;"
    "lower=1049/267540;upper=8890816603/2215998983205;"
    "gap=2101696427/23046389425332"
)
EXPECTED_STAGE = (
    11,
    "3e0859bbe77c95232a3d1096176fe7e1f0ebfd6ddd95f2ded98f1258278ba5ce",
    4,
    "8dcd5e390f85ab0d08979b4b35a80cc354c3c7ebe318d5275ac048a4fa03df93",
    7,
    "cb7ff868ba357f68b21009f640ac5dd4eaaedc92c9210e9b3964c91d8ce3bce0",
    0,
    "2e38e77b22c314a449e91fafed92a43826ac6aa403ae6a8acb6cf58239fbaf5d",
)
EXPECTED_SIGNS = Counter({-1: 5879, 0: 1, 1: 5879})
EXPECTED_CRUDE_WITNESSES = Counter({(420, 7, 6, 5): 4})
EXPECTED_STATUS_MODULI = Counter({(840, 7): 7})


def require(condition: bool, message: object) -> None:
    if not condition:
        raise RuntimeError(message)


def file_sha256(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def ftext(value: F) -> str:
    return f"{value.numerator}/{value.denominator}"


require(file_sha256(PARENT) == EXPECTED_PARENT_SHA256, "z1824 engine changed")
require(
    file_sha256(PARENT_OUTPUT) == EXPECTED_PARENT_OUTPUT_SHA256,
    "z1824 output changed",
)
require(
    EXPECTED_PARENT_CAP in PARENT_OUTPUT.read_text(encoding="utf-8").splitlines(),
    "inherited z1<=1812 cap changed",
)
require(file_sha256(BAND_SOURCE) == EXPECTED_BAND_SOURCE_SHA256, "band source changed")
require(file_sha256(BAND_OUTPUT) == EXPECTED_BAND_OUTPUT_SHA256, "band output changed")

BAND_LINES = BAND_OUTPUT.read_text(encoding="utf-8").splitlines()
require(EXPECTED_BAND_ROW in BAND_LINES, "the unique z1=1812 row changed")
BAND_HEIGHTS = tuple(
    sorted(
        int(field.removeprefix("z1="))
        for line in BAND_LINES
        if line.startswith("SURVIVOR;")
        for field in line.split(";")
        if field.startswith("z1=")
    )
)
require(BAND_HEIGHTS == EXPECTED_HEIGHTS, "global band height profile changed")
require(
    max(height for height in BAND_HEIGHTS if height < FIRST) == NEXT_CAP,
    "next occupied scalar height changed",
)

SPEC = spec_from_file_location("z1812_guarded_parent", PARENT)
require(SPEC is not None and SPEC.loader is not None, "cannot load z1824 parent")
A = module_from_spec(SPEC)
SPEC.loader.exec_module(A)
C = A.C
U = A.U
K = A.K

# Reparameterize the generic, guarded z1836 ray/status engine.  This rebuilds
# the full denominator ledger; no z1824 or z1836 state is carried forward.
C.FIRST = FIRST
C.BODIES = (BODY,)
C.EXPECTED_SIGNS = {BODY: EXPECTED_SIGNS}
C.EXPECTED = {BODY: EXPECTED_STAGE}


def compute_profile():
    carrier = U.suffix.A.carrier_for(BODY)
    h = F(sum(right - left for left, right in carrier), U.suffix.A.RULER)
    lower = h * U.suffix.ETAS[2]
    L = 14 * lcm(*BODY)
    first_delta = U.delta(carrier, h, FIRST)
    first_d = L // gcd(L, FIRST)
    wall = U.suffix.PROJECTED_RATIOS[2] * L
    high_floor = max(15, wall.numerator // wall.denominator + 1)
    require(
        (h, len(carrier), L, lower, first_delta, first_d, high_floor)
        == (F(1049, 2940), 34, 11760, F(1049, 267540), F(377, 1035860), 980, 1020),
        "z1=1812 frontier geometry changed",
    )
    require(FIRST >= high_floor, "z1=1812 fell below the projected wall")

    (
        _amplitudes,
        ray_digest,
        divisor_count,
        trials,
        scalar,
        crude_kills,
        status_kills,
        status_survivors,
        stage_digests,
    ) = C.scalar_status_frontier(BODY, carrier, h, lower, L, first_delta, first_d)
    require(divisor_count == 59 and trials == 557845, "denominator universe changed")
    require(
        (len(scalar), len(crude_kills), len(status_kills), len(status_survivors))
        == (11, 4, 7, 0),
        "z1=1812 closure counts changed",
    )

    actual_L, ranges = U.support.safe_cell_ranges(BODY)
    require(actual_L == L, "safe-cell ruler changed")
    crude_witnesses = Counter()
    for ds, _upper, _labels, witness in crude_kills:
        require(witness is not None, "crude kill lost its witness")
        q, M, target, capacity = witness
        D = lcm(*ds)
        arcs = U.fibre.projected_support_arcs(D, ranges)
        histogram = U.fibre.residue_load_histogram(arcs, q)
        recomputed_target = max(load for load, count in histogram if count)
        recomputed_capacity = sum(U.fibre_cap(D, d, q) for d in ds)
        require(
            D == q * M
            and (target, capacity) == (recomputed_target, recomputed_capacity)
            and target > capacity,
            "crude fibre witness failed exact reconstruction",
        )
        crude_witnesses[(q, M, target, capacity)] += 1
    require(crude_witnesses == EXPECTED_CRUDE_WITNESSES, "crude census changed")

    status_moduli = Counter()
    exact_dual_rows = []
    for ds, upper, labels, witness in status_kills:
        require(witness is not None, "status kill lost its witness")
        q, M, marginals, capacity_values, histogram, certificate = witness
        D = lcm(*ds)
        arcs = U.fibre.projected_support_arcs(D, ranges)
        recomputed_histogram = U.fibre.residue_load_histogram(arcs, q)
        recomputed_marginals, recomputed_capacities = K.hunter_status_data5(D, ds, q)
        require(
            D == q * M
            and recomputed_histogram == histogram
            and recomputed_marginals == marginals
            and tuple(sorted(set(recomputed_capacities))) == capacity_values,
            "status witness failed exact reconstruction",
        )
        dual_gap, minimum_slack = A.exact_farkas_audit(
            q, marginals, recomputed_capacities, histogram, certificate
        )
        status_moduli[(q, M)] += 1
        exact_dual_rows.append(
            (ds, upper, labels, q, M, dual_gap, minimum_slack, certificate)
        )
    require(status_moduli == EXPECTED_STATUS_MODULI, "status modulus census changed")
    require(not status_survivors, "a z1=1812 status survivor appeared")

    minimum_dual_gap = min(row[5] for row in exact_dual_rows)
    minimum_dual_slack = min(row[6] for row in exact_dual_rows)
    profile_payload = (
        BODY,
        FIRST,
        NEXT_CAP,
        QUANTIFIER,
        BAND_HEIGHTS,
        h,
        len(carrier),
        L,
        lower,
        first_delta,
        first_d,
        high_floor,
        ray_digest,
        divisor_count,
        trials,
        tuple(scalar),
        tuple(crude_kills),
        tuple(status_kills),
        stage_digests,
        tuple(sorted(crude_witnesses.items())),
        tuple(sorted(status_moduli.items())),
        tuple(exact_dual_rows),
    )
    profile_hash = sha256(repr(profile_payload).encode()).hexdigest()
    if EXPECTED_PROFILE_SHA256 is not None:
        require(profile_hash == EXPECTED_PROFILE_SHA256, "profile digest changed")
    return (
        h,
        len(carrier),
        L,
        lower,
        first_delta,
        first_d,
        high_floor,
        ray_digest,
        divisor_count,
        trials,
        len(scalar),
        len(crude_kills),
        len(status_kills),
        len(status_survivors),
        stage_digests,
        crude_witnesses,
        status_moduli,
        minimum_dual_gap,
        minimum_dual_slack,
        profile_hash,
    )


def render(profile) -> str:
    (
        h,
        components,
        L,
        lower,
        first_delta,
        first_d,
        high_floor,
        ray_digest,
        divisor_count,
        trials,
        scalar_count,
        crude_count,
        status_count,
        survivor_count,
        stage_digests,
        crude_witnesses,
        status_moduli,
        minimum_dual_gap,
        minimum_dual_slack,
        profile_hash,
    ) = profile
    semantic_payload = (
        EXPECTED_PARENT_SHA256,
        EXPECTED_PARENT_OUTPUT_SHA256,
        EXPECTED_BAND_SOURCE_SHA256,
        EXPECTED_BAND_OUTPUT_SHA256,
        BODY,
        FIRST,
        NEXT_CAP,
        QUANTIFIER,
        BAND_HEIGHTS,
        scalar_count,
        crude_count,
        status_count,
        survivor_count,
        stage_digests,
        tuple(sorted(crude_witnesses.items())),
        tuple(sorted(status_moduli.items())),
        minimum_dual_gap,
        minimum_dual_slack,
        profile_hash,
    )
    semantic_hash = sha256(repr(semantic_payload).encode()).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(semantic_hash == EXPECTED_SEMANTIC_SHA256, "semantic digest changed")
    lines = [
        "LRC14 projected k=2 z1=1812 exact ray/status closure",
        f"z1824_parent_source_sha256={file_sha256(PARENT)}",
        f"z1824_parent_output_sha256={file_sha256(PARENT_OUTPUT)}",
        f"scalar_band_source_sha256={file_sha256(BAND_SOURCE)}",
        f"scalar_band_output_sha256={file_sha256(BAND_OUTPUT)}",
        (
            "scope=unique global scalar row at z1=1812;"
            f"{QUANTIFIER};no finite label horizon"
        ),
        (
            f"E={','.join(map(str, BODY))};h={ftext(h)};r={components};L={L};"
            f"lower={ftext(lower)};first_delta={ftext(first_delta)};"
            f"first_d={first_d};projected_wall={high_floor}"
        ),
        (
            "ray_law=(z+L)delta(z+L)=zdelta(z);"
            f"ray_sha256={ray_digest};denominators={divisor_count};trials={trials}"
        ),
        (
            f"scalar_states={scalar_count};crude_all_divisor_kills={crude_count};"
            f"common_K5_status_kills={status_count};survivors={survivor_count}"
        ),
        f"crude_witness_census={tuple(sorted(crude_witnesses.items()))}",
        f"status_modulus_census={tuple(sorted(status_moduli.items()))}",
        (
            f"exact_Farkas_certificates={status_count};"
            f"minimum_dual_gap={ftext(minimum_dual_gap)};"
            f"minimum_pattern_slack={ftext(minimum_dual_slack)}"
        ),
        f"stage_sha256={stage_digests}",
        (
            "conclusion=the unique projected k=2 z1=1812 row is empty "
            "uniformly over later nonaligned labels"
        ),
        f"consequence=the proved projected k=2 first-drift cap is z1<={NEXT_CAP}",
        f"profile_sha256={profile_hash}",
        f"semantic_sha256={semantic_hash}",
        "all_exact_controls=PASS",
    ]
    return "\n".join(lines) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--hash-seed", type=int, default=0)
    parser.add_argument("--output", type=Path, default=OUTPUT_PATH)
    args = parser.parse_args()
    require(args.hash_seed >= 0, "hash seed must be nonnegative")
    os.environ["PYTHONHASHSEED"] = str(args.hash_seed)
    output = render(compute_profile())
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(output, encoding="utf-8", newline="\n")
    print(output, end="")


if __name__ == "__main__":
    main()
