#!/usr/bin/env python3
"""Exact THM-3626 companion for the archived R=8192 failed AMM trace."""

from __future__ import annotations

import ast
from hashlib import sha256
import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
PARENT_FILES = (
    (
        "theorem",
        ROOT / "01-canon/theorems/THM-3616-amm-R4096-adjoint-horizon-and-golden-finite-scale-obstruction.md",
        "48116058545de3e8614d7fede8159baaa6d900315d55894d7669c48c4f1f11bf",
    ),
    (
        "script",
        ROOT / "04-computation/amm_binary_sturmian_R4096_adjoint_horizon_thm3616.py",
        "c6753be05870700747aadb3f43e41a24e80f84f37b79e07fe70960f71a9e864a",
    ),
    (
        "output",
        ROOT / "05-knowledge/results/amm_binary_sturmian_R4096_adjoint_horizon_thm3616.out",
        "519e52f14ad8dc8cb42a182d769b00847eed1a19783958184adf36f87c8d679f",
    ),
)

EXPECTED = {
    "R": 8192,
    "D0": 191,
    "death": 2045,
    "fatal_bits": 4937,
    "profile": (5089, 6312, 9987),
    "first_negative": 774,
    "q": 773,
    "negative_count": 1271,
    "negative_endpoints": (774, 2044),
    "boundary_bits": (4939, 4930),
    "boundary_digest": "d42558aa87d6c8e57f739abdb54758ec192d0c6081d3ab850ce04c746cedecc3",
    "active_cells": 808356,
    "active_digest": "dcbbf6f02e3c4b2d01ad0d8e35a6e8b46a73fa8d08c07969306adc544c4a5ba5",
}
EXPECTED_SEMANTIC_SHA256 = (
    "226de0151e7fe6f9ded73b11041bb83d259f8d4d5d44c738de8f9b706d052d59"
)

CHECKS = 0


def require(condition: bool, payload: object) -> None:
    """Optimization-safe exact gate."""
    global CHECKS
    if condition is not True:
        raise RuntimeError(payload)
    CHECKS += 1


def lf_sha256(path: Path) -> str:
    return sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def digest_json(value: object) -> str:
    payload = json.dumps(value, sort_keys=True, separators=(",", ":")).encode("ascii")
    return sha256(payload).hexdigest()


for label, path, expected_hash in PARENT_FILES:
    observed = lf_sha256(path)
    require(observed == expected_hash, ("parent drift", label, observed, expected_hash))

sys.path.insert(0, str(ROOT / "04-computation"))
import amm_binary_sturmian_R4096_adjoint_horizon_thm3616 as A


def main() -> None:
    control_profile = A.degree_profile(4096, 88)
    require(
        (control_profile[0], control_profile[1014], control_profile[-1])
        == A.EXPECTED_4096["profile"],
        "R4096 profile control",
    )
    control_death, _, _ = A.first_death_from_top(4096, control_profile, 4096 // 3)
    require(control_death == A.EXPECTED_4096["death"], "R4096 death control")

    r_value = EXPECTED["R"]
    profile = A.degree_profile(r_value, EXPECTED["D0"])
    require(
        all(profile[index] - profile[index - 1] in (0, 1)
            for index in range(1, len(profile))),
        "binary Sturmian degree word",
    )
    death, scan_fatal, initial_width = A.first_death_from_top(
        r_value, profile, r_value // 3
    )
    require(death == EXPECTED["death"], ("death", death))
    require(initial_width == r_value // 3 + 1, ("initial width", initial_width))
    result = A.adjoint_horizon(r_value, profile, death)
    require(result["fatal"] == scan_fatal < 0, "fatal replay agreement")
    require(abs(scan_fatal).bit_length() == EXPECTED["fatal_bits"], "fatal bits")

    profile_control = (profile[0], profile[death], profile[-1])
    require(profile_control == EXPECTED["profile"], ("profile", profile_control))
    require(result["first_negative"] == EXPECTED["first_negative"], "wall")
    require(result["q"] == EXPECTED["q"], "departure horizon")
    require(death - result["first_negative"] == EXPECTED["negative_count"], "wall size")
    require(
        (result["first_negative"], death - 1) == EXPECTED["negative_endpoints"],
        "wall endpoints",
    )
    require(result["boundary_bits"] == EXPECTED["boundary_bits"], "boundary bits")
    require(result["boundary_digest"] == EXPECTED["boundary_digest"], "boundary digest")
    require(result["active"][0] == EXPECTED["active_cells"], "active cells")
    active_digest = A.digest_json(result["active"])
    require(active_digest == EXPECTED["active_digest"], "active digest")

    # Put all three dyadic golden errors over denominator 8192:
    # |q/R-theta_gold|=(a-1024*sqrt(5))/8192.
    radical_square = 5 * 1024**2
    error_numerators = (2292, 2308, 2299)  # R=2048,4096,8192.
    require(all(value * value > radical_square for value in error_numerators),
            "negative signed golden errors")
    require(error_numerators[1] - error_numerators[2] == 9,
            "R8192 rebound from R4096")
    require(error_numerators[2] - error_numerators[0] == 7,
            "R8192 remains farther than R2048")
    require(result["q"] - 2 * A.EXPECTED_4096["q"] == 9,
            "dyadic q defect")
    require(death - 2 * A.EXPECTED_4096["death"] == 17,
            "dyadic death defect")

    semantic_record = {
        "parent": tuple((label, digest) for label, _, digest in PARENT_FILES),
        "control": (4096, 88, control_death, A.EXPECTED_4096["profile"]),
        "expected": EXPECTED,
        "error_numerators": error_numerators,
        "radical_square": radical_square,
        "dyadic_defects": (9, 17),
    }
    semantic = digest_json(semantic_record)
    if EXPECTED_SEMANTIC_SHA256 != "TO_BE_PINNED":
        require(semantic == EXPECTED_SEMANTIC_SHA256,
                ("semantic", semantic, EXPECTED_SEMANTIC_SHA256))

    source_path = Path(__file__).resolve()
    source_bytes = source_path.read_bytes()
    require(b"\r\n" not in source_bytes, "source is not raw LF")
    require(not any(isinstance(node, ast.Assert)
                    for node in ast.walk(ast.parse(source_bytes.decode("utf-8")))),
            "Python assert node present")

    print("== THM-3626 AMM R8192 exact adjoint horizon ==")
    print("parent_sha256_lf", tuple((label, digest) for label, _, digest in PARENT_FILES))
    print("control_R4096", (88, control_death, A.EXPECTED_4096["profile"]))
    print("epoch", (r_value, EXPECTED["D0"]))
    print("death_fatal_bits_profile", (death, EXPECTED["fatal_bits"], profile_control))
    print("wall", (result["first_negative"], result["q"], EXPECTED["negative_endpoints"]))
    print("boundary", (result["boundary_bits"], result["boundary_digest"]))
    print("active", (result["active"][0], active_digest))
    print("q_over_R", (result["q"], r_value))
    print("golden_error", "(-2299+1024*sqrt(5))/8192<0")
    print("error_rebound", "9/8192 closer than R4096;7/8192 farther than R2048")
    print("dyadic_defects", "q8192-2*q4096=9;death8192-2*death4096=17")
    print("semantic_sha256", semantic)
    print("script_sha256_lf", lf_sha256(source_path))
    print("CHECKS", CHECKS)
    print("SCOPE", "fixed failed Rule-A prefix; no alternative feasibility/asymptotic AMM bound")
    print("PASS")


if __name__ == "__main__":
    main()
