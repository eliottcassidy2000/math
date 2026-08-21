#!/usr/bin/env python3
"""Exact companion for THM-3648's terminal local AMM adjoint wall.

THM-3644 proves that, inside its complete R=16384 bracket, D0=412 dies and
D0=413 closes.  This companion applies THM-3633's checkpointed positive
adjoint to the failed D0=412 trace.  It also compares the resulting phase
with the audited R=8192 failed trace and the nonterminal R=16384,D0=400
trace.  The computation concerns one Rule-A prefix and makes no full-entry
or asymptotic claim.
"""

from __future__ import annotations

import ast
from hashlib import sha256
import json
from pathlib import Path
import sys


if hasattr(sys, "set_int_max_str_digits"):
    sys.set_int_max_str_digits(0)

ROOT = Path(__file__).resolve().parents[1]
COMPUTATION = ROOT / "04-computation"
sys.path.insert(0, str(COMPUTATION))

PARENT_FILES = (
    (
        "THM3633_theorem",
        ROOT / "01-canon/theorems/THM-3633-amm-r16384-fixed-failed-trace-adjoint-phase-shock.md",
        "936afc9acbfff61b9398502d2ea6114ec2b6f852098b06680a33fc9480decd07",
    ),
    (
        "THM3633_script",
        COMPUTATION / "amm_binary_sturmian_R16384_D400_adjoint_phase_shock_thm3633.py",
        "a4888d29f0eb1f59debea6f1545b16928bd9278af77bfa9ae5f5480ee08bdab4",
    ),
    (
        "THM3633_output",
        ROOT / "05-knowledge/results/amm_binary_sturmian_R16384_D400_adjoint_phase_shock_thm3633.out",
        "6754b09009d72ecee1fb5e5e3230adf2a660578cfa34a047dd154ddb28088104",
    ),
    (
        "THM3644_theorem",
        ROOT / "01-canon/theorems/THM-3644-amm12592-exact-offset-threshold-at-R16384.md",
        "1dcde7cf49f8dd757a21f57061fcaf59773edea0f3abc6ad3efda96adc2c6e28",
    ),
    (
        "THM3644_script",
        COMPUTATION / "amm12592_R16384_offset_transition_thm3644.py",
        "ddf248e8939bbdefa7b2544bbd3df1c23e47e53e00aae4013d09149eae2ca59c",
    ),
    (
        "THM3644_output",
        ROOT / "05-knowledge/results/amm12592_R16384_offset_transition_thm3644.out",
        "28cedef2dc0b176b62f0633d93ea23a7c316993ed197743bd54f4466eb860c21",
    ),
)

EXPECTED = {
    "epoch": (16_384, 412),
    "death": 4_116,
    "profile": (10_209, 12_670, 20_006),
    "profile_digest": "b4c95d7379b5f31b0d48e30aeab9d9912b1fead409a51b519aea2f3bfaa5ff43",
    "fatal_bits": 9_919,
    "fatal_digest": "da8aba4bb9ab4ab93723461c1a98380604999fda9e5bc880d1065fd1b1fd40c1",
    "first_negative": 1_565,
    "q": 1_564,
    "negative_endpoints": (1_565, 4_115),
    "boundary_bits": (9_920, 9_919),
    "boundary_digest": "cb782c23c903e5fe244baf206f745259e91122a10bd7d3feb13a6e5be0fb6a58",
    "active_cells": 3_255_076,
    "active_digest": "7012e2dca4351711ac35d9842e46f48be8b2b5ce1aba9f12e68bb0784442ce5b",
    "wall_digest": "b40324f8cb11b874b043bed7a3f3c55757b0a297d4ebac0731a4f0556b0d73d8",
    "phase": {
        "h": 31,
        "edge_bound": 4_036,
        "margin": 80,
        "ell": 2_551,
        "ell_error": (4_599, -2_048),
        "q_error": (-4_580, 2_048),
    },
}
EXPECTED_SEMANTIC_SHA256 = "56e327478603425cc77ae386681b25d7a0ef29778f59183538e35f704609c9e7"
CHECKS = 0


def require(condition: bool, payload: object) -> None:
    global CHECKS
    if condition is not True:
        raise RuntimeError(payload)
    CHECKS += 1


def raw_sha256(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def digest(value: object) -> str:
    return sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode("ascii")
    ).hexdigest()


for label, path, expected_hash in PARENT_FILES:
    require(raw_sha256(path) == expected_hash, (label, "parent drift"))

import amm_binary_sturmian_R16384_D400_adjoint_phase_shock_thm3633 as T
import amm12592_R16384_offset_transition_thm3644 as O


def main() -> None:
    r_value, d0 = EXPECTED["epoch"]
    death = EXPECTED["death"]

    row_412 = next(row for row in O.ROWS if row[0] == 412)
    row_413 = next(row for row in O.ROWS if row[0] == 413)
    require(row_412[:9] == (412, "DIE", 4116, 10209, 12670, 20006,
                            2179, 4036, 9920), "THM3644 failed row")
    require(row_413[:9] == (413, "CLOSED", 10116, 10210, 16259, 20007,
                            8191, 4038, None), "THM3644 adjacent closure")

    profile, last_degree, profile_digest = T.exact_profile_prefix(
        r_value, d0, death
    )
    result = T.checkpointed_adjoint(r_value, profile, death)
    phase = T.phase_invoice(
        r_value, profile[0], death, result["first_negative"], result["q"]
    )

    observed = {
        "epoch": (r_value, d0),
        "death": death,
        "profile": (profile[0], profile[death], last_degree),
        "profile_digest": profile_digest,
        "fatal_bits": abs(result["fatal"]).bit_length(),
        "fatal_digest": result["fatal_digest"],
        "first_negative": result["first_negative"],
        "q": result["q"],
        "negative_endpoints": result["negative_endpoints"],
        "boundary_bits": result["boundary_bits"],
        "boundary_digest": result["boundary_digest"],
        "active_cells": result["active"][0],
        "active_digest": result["active_digest"],
        "wall_digest": result["wall_digest"],
        "phase": phase,
    }
    require(observed == EXPECTED, ("target ledger", observed))

    control_8192 = T.EXPECTED[8192]
    control_400 = T.EXPECTED[16384]
    phase_400 = T.phase_invoice(
        16_384,
        control_400["profile"][0],
        control_400["death"],
        control_400["first_negative"],
        control_400["q"],
    )

    # All three comparisons are exact.  The signed q/R errors are negative:
    # terminal D412: (-1145+512 sqrt(5))/4096;
    # R8192 D191:   (-2299+1024 sqrt(5))/8192.
    require(1145 * 1145 - 5 * 512 * 512 == 305,
            "terminal radical square gap")
    require(2299 * 2299 - 5 * 1024 * 1024 == 42_521,
            "R8192 radical square gap")
    require(EXPECTED["q"] - 2 * control_8192["q"] == 18,
            "terminal dyadic q defect")
    require(2299 - 2290 == 9,
            "common-denominator numerator improvement")
    require(EXPECTED["q"] - control_400["q"] == 48,
            "same-scale q recovery")
    require(EXPECTED["death"] - control_400["death"] == 61,
            "same-scale death shift")
    require(EXPECTED["phase"]["ell"] - phase_400["ell"] == 13,
            "same-scale depth shift")
    require(EXPECTED["phase"]["margin"] - phase_400["margin"] == 37,
            "same-scale margin shift")
    require(EXPECTED["phase"]["h"] - phase_400["h"] == -12,
            "same-scale headroom shift")
    require(37 - 2 * (-12) - 13 == 48, "phase-invoice shift")

    semantic_record = {
        "parents": tuple((label, expected_hash)
                         for label, _path, expected_hash in PARENT_FILES),
        "target": observed,
        "adjacent": (row_412[:9], row_413[:9]),
        "golden": {
            "terminal_error": "(-1145+512*sqrt(5))/4096<0",
            "terminal_square_gap": 305,
            "R8192_error": "(-2299+1024*sqrt(5))/8192<0",
            "R8192_square_gap": 42_521,
            "absolute_error_improvement": "9/8192",
        },
        "same_scale_shift_D400_to_D412": {
            "q": 48,
            "death": 61,
            "h": -12,
            "margin": 37,
            "ell": 13,
            "invoice": "37-2*(-12)-13=48",
        },
        "scope": "terminal failure only within THM3644 bracket; no global threshold or AMM bound",
    }
    semantic = digest(semantic_record)
    if EXPECTED_SEMANTIC_SHA256 != "TO_BE_PINNED":
        require(semantic == EXPECTED_SEMANTIC_SHA256,
                ("semantic digest", semantic))

    source = Path(__file__).resolve()
    source_bytes = source.read_bytes()
    require(b"\r\n" not in source_bytes, "source raw LF")
    require(not any(isinstance(node, ast.Assert)
                    for node in ast.walk(ast.parse(source_bytes.decode("utf-8")))),
            "Python assert node present")

    print("== THM-3648 R16384,D0=412 terminal-local-failure adjoint ==")
    print("parent_sha256_raw=" + repr(semantic_record["parents"]))
    print("adjacent=412:DIE@4116;413:CLOSED@10116;scope:THM3644 bracket")
    print("profile=" + repr(observed["profile"]) + ";profile_digest=" + profile_digest)
    print(f"fatal_bits={observed['fatal_bits']};fatal_digest={observed['fatal_digest']}")
    print("wall=(first_negative:1565,q:1564,endpoints:(1565,4115));wall_digest="
          + observed["wall_digest"])
    print("boundary=bits:(9920,9919);digest=" + observed["boundary_digest"])
    print(f"active_cells={observed['active_cells']};active_digest={observed['active_digest']}")
    print("phase=" + repr(phase))
    print("golden_error=(-1145+512*sqrt(5))/4096<0;square_gap=305")
    print("versus_R8192_D191=absolute_error_improves_by_9/8192;q_defect=+18")
    print("versus_same_scale_D400=(delta_q,delta_death,delta_h,delta_margin,delta_ell)="
          "(48,61,-12,37,13);invoice=37-2*(-12)-13=48")
    print("semantic_sha256=" + semantic)
    print(f"CHECKS={CHECKS};parent_checks={T.CHECKS + O.CHECKS}")
    print("status=FINITE-EXACT TERMINAL-LOCAL-FAILURE PREFIX")
    print("scope=no global offset monotonicity/full-entry infeasibility/asymptotic law/AMM bound")


if __name__ == "__main__":
    main()
