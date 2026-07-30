#!/usr/bin/env python3
"""All-label ray/status closure of the projected ``k=3, z1=364`` spike.

The canonical one-pass scalar atlas reconstructs exactly 25 six-body rows at
first drift 364 and no scalar rows at any interstitial height from 365 through
377.  This script freezes that complete body tuple and pins both the atlas
source and transcript.  It then applies the exact residue-ray quotient and
the common 16-status Hunter/Farkas referee used at ``z1=378``.  Consequently
the computation has no label horizon and no dependency on ``.scratch``.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
ATLAS_PATH = (
    ROOT
    / "04-computation"
    / "lrc14_j7_k3_projected_scalar_atlas_thm2941.py"
)
ATLAS_OUTPUT_PATH = (
    ROOT
    / "05-knowledge"
    / "results"
    / "lrc14_j7_k3_projected_scalar_atlas_thm2941.out"
)
RAY_ENGINE_PATH = (
    ROOT
    / "04-computation"
    / "lrc14_j7_k3_z378_ray_status_closure_thm2941.py"
)
OUTPUT_PATH = (
    ROOT
    / "05-knowledge"
    / "results"
    / "lrc14_j7_k3_z364_ray_status_closure_thm2941.out"
)
EXPECTED_ATLAS_SHA256 = (
    "ddb2d19c02c4d70cfa74141265ceac585685932e85768baf4cb98aeb3e37935b"
)
EXPECTED_ATLAS_OUTPUT_SHA256 = (
    "ce6807de6d6b7022c97839d0bf9fc8ba3b90e7b97bc5b0d4069e88563e232be6"
)
EXPECTED_RAY_ENGINE_SHA256 = (
    "813d93ffe67cf3bdab06aa27f04680f49b057e387eabc8e713ba9cba35439353"
)
EXPECTED_SEMANTIC_SHA256 = (
    "008294c359301867a83262ebad5696d9ea8887db4161dc4829289c50e9f710da"
)

FIRST = 364
ORIGINAL_BODIES = (
    (1, 2, 6, 8, 12, 14),
    (1, 2, 6, 10, 12, 14),
    (1, 2, 8, 10, 12, 14),
    (1, 4, 6, 8, 12, 14),
    (1, 4, 6, 10, 12, 14),
    (1, 4, 8, 10, 12, 14),
    (1, 6, 8, 9, 12, 14),
    (1, 6, 8, 10, 12, 14),
    (2, 3, 8, 10, 12, 14),
    (2, 3, 8, 11, 12, 14),
    (2, 4, 6, 8, 12, 14),
    (2, 4, 6, 10, 12, 14),
    (2, 4, 8, 10, 12, 14),
    (2, 6, 8, 9, 12, 14),
    (2, 6, 8, 10, 12, 14),
    (2, 6, 8, 10, 13, 14),
    (2, 6, 8, 11, 12, 14),
    (2, 6, 9, 10, 12, 14),
    (2, 6, 10, 11, 12, 14),
    (2, 8, 9, 10, 12, 14),
    (2, 8, 10, 11, 12, 14),
    (2, 8, 10, 12, 13, 14),
    (3, 4, 8, 10, 12, 14),
    (4, 6, 8, 10, 12, 14),
    (4, 8, 10, 11, 12, 14),
)
EXPECTED_COUNTS = dict(
    zip(
        ORIGINAL_BODIES,
        (
            (1, 0, 1),
            (9, 0, 9),
            (18, 1, 17),
            (5, 0, 5),
            (8, 0, 8),
            (331, 111, 220),
            (1, 0, 1),
            (67, 5, 62),
            (9, 1, 8),
            (5, 0, 5),
            (5, 0, 5),
            (5, 0, 5),
            (193, 29, 164),
            (2, 0, 2),
            (151, 14, 137),
            (0, 0, 0),
            (46, 0, 46),
            (1, 0, 1),
            (0, 0, 0),
            (2, 2, 0),
            (0, 0, 0),
            (237, 78, 159),
            (4, 1, 3),
            (9, 0, 9),
            (0, 0, 0),
        ),
        strict=True,
    )
)


def require(condition: bool, message: object) -> None:
    if not condition:
        raise RuntimeError(message)


def file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, ("cannot load", path))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


require(
    file_sha256(ATLAS_PATH) == EXPECTED_ATLAS_SHA256,
    "projected scalar atlas source changed",
)
require(
    file_sha256(ATLAS_OUTPUT_PATH) == EXPECTED_ATLAS_OUTPUT_SHA256,
    "projected scalar atlas transcript changed",
)
require(
    file_sha256(RAY_ENGINE_PATH) == EXPECTED_RAY_ENGINE_SHA256,
    "z1=378 ray/status engine changed",
)
atlas = load_module("k3_projected_scalar_atlas", ATLAS_PATH)
require(
    atlas.EXPECTED_Z364_BODIES == ORIGINAL_BODIES,
    "z1=364 atlas/body handoff changed",
)
ray = load_module("k3_z378_ray_status_engine", RAY_ENGINE_PATH)
ray.FIRST = FIRST


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUTPUT_PATH)
    args = parser.parse_args()

    pair_rows, control_marginals, control_caps, control_certificate = (
        ray.local.controls()
    )
    records = []
    total_states = total_crude = total_status = total_survivors = 0
    certificate_digest = hashlib.sha256()
    for body in ORIGINAL_BODIES:
        stream = ray.Stream(body)
        trials, states, checks, signs = ray.ray_quotient_states(stream)
        crude, status, survivors = ray.common_status_screen(stream, states)
        require(
            (len(states), len(crude), len(status)) == EXPECTED_COUNTS[body],
            (
                "body closure counts changed",
                body,
                len(states),
                len(crude),
                len(status),
            ),
        )
        maximum = max(
            states.items(),
            key=lambda item: (item[1]["excess"], item[0]),
            default=None,
        )
        sign_totals = {
            sign: sum(
                count
                for (_d, candidate_sign), count in signs.items()
                if candidate_sign == sign
            )
            for sign in (-1, 0, 1)
        }
        require(
            sign_totals[-1] == sign_totals[1],
            ("ray sign imbalance", body, sign_totals),
        )
        status_histogram = Counter(witness[1] for witness in status.values())
        for ds, witness in sorted(status.items()):
            certificate_digest.update(f"{body}|{ds}|{witness}\n".encode())
        records.append(
            (
                body,
                stream.h,
                len(stream.carrier),
                stream.L,
                stream.high_floor,
                stream.first_d,
                checks,
                tuple(sorted(sign_totals.items())),
                len(ray.support.divisors(stream.L)) - 1,
                trials,
                tuple(sorted(states.items())),
                tuple(sorted(crude.items())),
                tuple(sorted(status.items())),
                tuple(survivors),
                tuple(sorted(status_histogram.items())),
                maximum,
            )
        )
        total_states += len(states)
        total_crude += len(crude)
        total_status += len(status)
        total_survivors += len(survivors)

    require(
        (total_states, total_crude, total_status, total_survivors)
        == (1_109, 242, 867, 0),
        (
            "z1=364 closure totals changed",
            total_states,
            total_crude,
            total_status,
            total_survivors,
        ),
    )
    semantic_payload = (
        FIRST,
        ORIGINAL_BODIES,
        tuple(records),
        pair_rows,
        control_marginals,
        control_caps,
        control_certificate,
        certificate_digest.hexdigest(),
    )
    semantic_sha256 = hashlib.sha256(repr(semantic_payload).encode()).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(
            semantic_sha256 == EXPECTED_SEMANTIC_SHA256,
            "z1=364 closure semantic digest changed",
        )

    lines = [
        "LRC14 projected k=3 z1=364 all-label ray/status closure",
        f"atlas_source_sha256={file_sha256(ATLAS_PATH)}",
        f"atlas_output_sha256={file_sha256(ATLAS_OUTPUT_PATH)}",
        f"ray_engine_source_sha256={file_sha256(RAY_ENGINE_PATH)}",
        (
            "scope=25 canonical projected scalar body rows at z1=364;"
            "three distinct later nonaligned labels;no finite horizon"
        ),
        (
            "ray_law=(z+L)delta(z+L)=zdelta(z);"
            "A(L-b)=-A(b);all denominator maxima attained"
        ),
        f"pair_overlap_exhaustive_controls={pair_rows}",
        (
            "common_table_controls=loads(9,9):FEASIBLE;"
            "loads(14,8):EXACT_FARKAS_INFEASIBLE"
        ),
    ]
    for record in records:
        (
            body,
            h,
            components,
            L,
            high_floor,
            first_d,
            checks,
            sign_totals,
            denominator_classes,
            trials,
            states,
            crude,
            status,
            survivors,
            status_histogram,
            maximum,
        ) = record
        lines.extend(
            (
                (
                    f"E={body};h={ray.ftext(h)};r={components};L={L};"
                    f"high={high_floor};d1={first_d};ray_checks={checks};"
                    f"ray_signs={dict(sign_totals)};"
                    f"denominator_classes={denominator_classes};"
                    f"denominator_trials={trials};scalar_states={len(states)};"
                    f"crude_kills={len(crude)};status_kills={len(status)};"
                    f"survivors={len(survivors)};"
                    f"status_M_histogram={dict(status_histogram)}"
                ),
                f"  scalar_state_denominators={tuple(ds for ds, _ in states)}",
                f"  crude_kill_denominators={tuple(ds for ds, _ in crude)}",
                f"  status_kill_denominators={tuple(ds for ds, _ in status)}",
                f"  max_state={maximum}",
            )
        )
    lines.extend(
        (
            (
                f"totals=states:{total_states};crude_kills:{total_crude};"
                f"status_kills:{total_status};survivors:{total_survivors}"
            ),
            f"farkas_certificate_sha256={certificate_digest.hexdigest()}",
            "conclusion=all 25 z1=364 candidate body rows are empty",
            f"semantic_sha256={semantic_sha256}",
            "all_exact_controls=PASS",
        )
    )
    payload = "\n".join(lines) + "\n"
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(payload)
    print(payload, end="")


if __name__ == "__main__":
    main()
