#!/usr/bin/env python3
"""Exact all-nonaligned-label closure of the projected k=2 row at z1=1824.

The global scalar atlas on ``1810 <= z1 <= 1835`` leaves one row at its
largest occupied height,

    E=(1,4,8,10,12,14),  L=11760,  z1=1824.

This verifier reuses the guarded residue-ray and common K5 Hunter-status
engine that closed the exact-suffix part of the z1=1836 slice.  It does not
reuse that slice's denominator ledger: the first denominator changes from
980 to 245, so all ``C(59+3,4)=557845`` four-suffix denominator multisets are
reconstructed from scratch.

For every later nonaligned label z, endpoint periodicity gives the exact ray law

    delta_E(z+L) = z*delta_E(z)/(z+L).

Thus, for each denominator multiset, the first four positive points on every
primitive residue ray determine its exact nonaligned-label scalar maximum.  Only 38
states meet the scalar lower bound.  Exact all-divisor fibre capacities kill
15; exact rational Farkas certificates for the common 32-cell K5 status table
kill the other 23.  Hence the z1=1824 row is empty uniformly over all later
distinct nonaligned labels.  There is no finite label horizon and no projected-residual
fallback in this closure.
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
ENGINE = (
    ROOT
    / "04-computation"
    / "lrc14_j7_k2_z1836_exact_body_projected_closure_thm2941.py"
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
    / "lrc14_j7_k2_z1824_ray_status_closure_thm2941.out"
)

EXPECTED_ENGINE_SHA256 = (
    "54ab8787af4dd410d3c809e5a15b94fb06e7c6129def4003021bb923f0dfb766"
)
EXPECTED_BAND_SOURCE_SHA256 = (
    "058c43d67d0bb110993ec687877edba4f5a7ad472a81455b0a5276b20db7680d"
)
EXPECTED_BAND_OUTPUT_SHA256 = (
    "76f08dc5b70c98dd7c8fa958f5598f5c50e7cb1df26d0be10531ba8185796952"
)
EXPECTED_PROFILE_SHA256 = (
    "e7e08e818cbb28d8235f6d1b719330ad51024dd5c5c71537b010ebfa030a9ef3"
)
EXPECTED_SEMANTIC_SHA256 = (
    "4a5e22f4a7d5a3d5832226dd8f54d19fadfb3979c5d60e6e0d7eead3af6e070d"
)

BODY = (1, 4, 8, 10, 12, 14)
FIRST = 1824
NEXT_CAP = 1812
QUANTIFIER = "four later distinct nonaligned labels"
EXPECTED_EMPTY_BANDS = "empty_first_bands=1813..1823,1825..1835"
EXPECTED_BAND_ROW = (
    "SURVIVOR;E=1,4,8,10,12,14;h=1049/2940;r=34;L=11760;"
    "largest_floor=1020;z1=1824;delta1=1447/3128160;"
    "suffix=EXACT:1836:3583/3148740,EXACT:2060:1867/2119740,"
    "EXACT:2172:263/310415,EXACT:2142:821/1049580;"
    "lower=1049/267540;upper=36678771989/8922697892640;"
    "gap=4403581813/23199014520864"
)
EXPECTED_STAGE = (
    38,
    "33e9ff8243d644706423af0394d2124d31f9b8661199e057c3daae6945b2e4ed",
    15,
    "90011e470875553d54bf83405af93e4c3a5f959a22ed9256711e3cffcd743069",
    23,
    "9d74b69817c80a2fa41cea5d40702566ff545b58a2db400c0d94938c2e1ffbf1",
    0,
    "2e38e77b22c314a449e91fafed92a43826ac6aa403ae6a8acb6cf58239fbaf5d",
)
EXPECTED_SIGNS = Counter({-1: 5879, 0: 1, 1: 5879})
EXPECTED_CRUDE_WITNESSES = Counter(
    {
        (420, 7, 6, 5): 10,
        (1, 1960, 1504, 1400): 3,
        (1, 980, 840, 700): 2,
    }
)
EXPECTED_STATUS_MODULI = Counter(
    {
        (840, 7): 14,
        (980, 6): 5,
        (980, 3): 1,
        (588, 5): 1,
        (420, 7): 1,
        (490, 6): 1,
    }
)


def require(condition: bool, message: object) -> None:
    if not condition:
        raise RuntimeError(message)


def file_sha256(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def ftext(value: F) -> str:
    return f"{value.numerator}/{value.denominator}"


require(file_sha256(ENGINE) == EXPECTED_ENGINE_SHA256, "closure engine changed")
require(
    file_sha256(BAND_SOURCE) == EXPECTED_BAND_SOURCE_SHA256,
    "1810..1835 scalar source changed",
)
require(
    file_sha256(BAND_OUTPUT) == EXPECTED_BAND_OUTPUT_SHA256,
    "1810..1835 scalar output changed",
)
require(
    EXPECTED_BAND_ROW in BAND_OUTPUT.read_text(encoding="utf-8").splitlines(),
    "the unique z1=1824 scalar row changed",
)
require(
    EXPECTED_EMPTY_BANDS
    in BAND_OUTPUT.read_text(encoding="utf-8").splitlines(),
    "the scalar-empty descent below z1=1824 changed",
)

SPEC = spec_from_file_location("z1824_guarded_ray_status_engine", ENGINE)
require(SPEC is not None and SPEC.loader is not None, "cannot load closure engine")
C = module_from_spec(SPEC)
SPEC.loader.exec_module(C)
U = C.U
K = C.K

# The imported routine is intentionally parameterized by these globals.  Pin
# all of them together before invoking it.  Its transitive dependency guards
# have already run during import.
C.FIRST = FIRST
C.BODIES = (BODY,)
C.EXPECTED_SIGNS = {BODY: EXPECTED_SIGNS}
C.EXPECTED = {BODY: EXPECTED_STAGE}


def exact_farkas_audit(q, marginals, capacities, histogram, certificate):
    """Independently recheck the returned common-status certificate over Q."""

    thresholds, alpha, z = certificate
    require(len(z) == 6, "common-status equality dual has wrong arity")
    tail_rows = []
    tail_rhs = []
    active_thresholds = []
    for threshold, _count in histogram:
        if threshold <= 0:
            continue
        demand = sum(count for load, count in histogram if load >= threshold)
        good = tuple(int(capacity >= threshold) for capacity in capacities)
        if all(good):
            continue
        active_thresholds.append(threshold)
        if not any(good):
            require(
                tuple(thresholds) == (threshold,)
                and tuple(alpha) == (F(1),)
                and tuple(z) == (F(0),) * 6
                and demand > 0,
                "degenerate no-status certificate changed",
            )
            return F(demand), F(0)
        tail_rows.append(good)
        tail_rhs.append(demand)

    require(tuple(active_thresholds) == tuple(thresholds), "dual thresholds changed")
    require(len(alpha) == len(tail_rows), "tail dual has wrong arity")
    equality_rows = [
        (1,) * 32,
        *[
            tuple((pattern >> index) & 1 for pattern in range(32))
            for index in range(5)
        ],
    ]
    equality_rhs = (q, *marginals)
    slacks = tuple(
        sum(z[row] * equality_rows[row][pattern] for row in range(6))
        - sum(alpha[row] * tail_rows[row][pattern] for row in range(len(alpha)))
        for pattern in range(32)
    )
    contradiction = (
        sum(z[row] * equality_rhs[row] for row in range(6))
        - sum(alpha[row] * tail_rhs[row] for row in range(len(alpha)))
    )
    require(
        all(value >= 0 for value in alpha)
        and all(value >= 0 for value in slacks)
        and contradiction < 0,
        "returned Farkas certificate failed exact rational replay",
    )
    return -contradiction, min(slacks)


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
        == (F(1049, 2940), 34, 11760, F(1049, 267540), F(1447, 3128160), 245, 1020),
        "z1=1824 frontier geometry changed",
    )
    require(FIRST >= high_floor, "z1=1824 fell below the projected wall")

    frontier = C.scalar_status_frontier(
        BODY, carrier, h, lower, L, first_delta, first_d
    )
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
    ) = frontier
    require(divisor_count == 59 and trials == 557845, "denominator universe changed")
    require(
        (len(scalar), len(crude_kills), len(status_kills), len(status_survivors))
        == (38, 15, 23, 0),
        "z1=1824 closure counts changed",
    )
    require(not status_survivors, "a z1=1824 common-status survivor appeared")

    crude_witnesses = Counter(
        (row[3][0], row[3][1], row[3][2], row[3][3]) for row in crude_kills
    )
    status_moduli = Counter((row[3][0], row[3][1]) for row in status_kills)
    require(crude_witnesses == EXPECTED_CRUDE_WITNESSES, "crude witness census changed")
    require(status_moduli == EXPECTED_STATUS_MODULI, "status modulus census changed")

    actual_L, ranges = U.support.safe_cell_ranges(BODY)
    require(actual_L == L, "safe-cell ruler changed")
    exact_dual_rows = []
    for ds, upper, labels, witness in status_kills:
        require(witness is not None, "status kill lost its witness")
        q, M, marginals, capacity_values, histogram, certificate = witness
        D = lcm(*ds)
        require(D == q * M, "status factorization changed")
        arcs = U.fibre.projected_support_arcs(D, ranges)
        recomputed_histogram = U.fibre.residue_load_histogram(arcs, q)
        recomputed_marginals, recomputed_capacities = K.hunter_status_data5(D, ds, q)
        require(
            recomputed_histogram == histogram
            and recomputed_marginals == marginals
            and tuple(sorted(set(recomputed_capacities))) == capacity_values,
            "status witness data failed independent reconstruction",
        )
        dual_gap, minimum_slack = exact_farkas_audit(
            q, marginals, recomputed_capacities, histogram, certificate
        )
        exact_dual_rows.append(
            (ds, upper, labels, q, M, dual_gap, minimum_slack, certificate)
        )

    minimum_dual_gap = min(row[5] for row in exact_dual_rows)
    minimum_dual_slack = min(row[6] for row in exact_dual_rows)
    profile_payload = (
        BODY,
        FIRST,
        NEXT_CAP,
        EXPECTED_EMPTY_BANDS,
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
        EXPECTED_ENGINE_SHA256,
        EXPECTED_BAND_SOURCE_SHA256,
        EXPECTED_BAND_OUTPUT_SHA256,
        BODY,
        FIRST,
        QUANTIFIER,
        h,
        components,
        L,
        lower,
        first_delta,
        first_d,
        high_floor,
        divisor_count,
        trials,
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
        "LRC14 projected k=2 z1=1824 exact ray/status closure",
        f"closure_engine_sha256={file_sha256(ENGINE)}",
        f"scalar_band_source_sha256={file_sha256(BAND_SOURCE)}",
        f"scalar_band_output_sha256={file_sha256(BAND_OUTPUT)}",
        (
            "scope=unique global scalar row at z1=1824;"
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
            "conclusion=the unique projected k=2 z1=1824 row is empty "
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
