#!/usr/bin/env python3
"""Exact THM-3633 companion for the fixed failed AMM trace (R,D0)=(16384,400).

It hash-pins the current THM-3626 and THM-3616 theorem/script/output triples,
imports their exact
Fibonacci--Lucas and top-distance primitives, reconstructs the complete
R=8192 control, and then certifies the fixed archived R=16384,D0=400 trace.

The adjoint is evaluated with block checkpoints.  This keeps only O(j^2/B)
primal checkpoint cells and O(Bj) temporary charge cells rather than the full
quadratic multiplier triangle.  Every mathematical gate is an unconditional
``require`` or explicit exception, so ``python`` and ``python -O`` have the
same semantics.  The result concerns one known failing prefix; D0=400 is not
known to be the terminal failing offset at R=16384.
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
PARENT_FILES = (
    (
        "THM3626_theorem",
        ROOT / "01-canon/theorems/THM-3626-amm-r8192-adjoint-horizon-and-phase-rebound.md",
        "174ea12b9b3e3679295d31210a77b1cd6ab093f380f8e62c7437d31023606a36",
    ),
    (
        "THM3626_script",
        ROOT / "04-computation/amm_binary_sturmian_R8192_adjoint_horizon_thm3626.py",
        "548e5ddaad0b7220687468552a448f3021ad16b3162a2c11409d139d31ee3be5",
    ),
    (
        "THM3626_output",
        ROOT / "05-knowledge/results/amm_binary_sturmian_R8192_adjoint_horizon_thm3626.out",
        "4319f792b1c69b84ea30496ef0500d3666da73575d4de267a60428d5c757300b",
    ),
    (
        "THM3616_theorem",
        ROOT / "01-canon/theorems/THM-3616-amm-R4096-adjoint-horizon-and-golden-finite-scale-obstruction.md",
        "48116058545de3e8614d7fede8159baaa6d900315d55894d7669c48c4f1f11bf",
    ),
    (
        "THM3616_script",
        ROOT / "04-computation/amm_binary_sturmian_R4096_adjoint_horizon_thm3616.py",
        "c6753be05870700747aadb3f43e41a24e80f84f37b79e07fe70960f71a9e864a",
    ),
    (
        "THM3616_output",
        ROOT / "05-knowledge/results/amm_binary_sturmian_R4096_adjoint_horizon_thm3616.out",
        "519e52f14ad8dc8cb42a182d769b00847eed1a19783958184adf36f87c8d679f",
    ),
)

EXPECTED = {
    8192: {
        "D0": 191,
        "death": 2045,
        "fatal_bits": 4937,
        "profile": (5089, 6312, 9987),
        "first_negative": 774,
        "q": 773,
        "boundary_bits": (4939, 4930),
        "boundary_digest": "d42558aa87d6c8e57f739abdb54758ec192d0c6081d3ab850ce04c746cedecc3",
        "active_cells": 808356,
        "active_digest": "dcbbf6f02e3c4b2d01ad0d8e35a6e8b46a73fa8d08c07969306adc544c4a5ba5",
        "fatal_digest": "1e78264343b58c90f5f0cb24302813b284cba1f9ba4ca60eb7740a74a15a37a9",
        "profile_digest": "a595c188bf04c1397bed449f23bfda56a13a9d6143d4d6f48ff23bb1a10d2df1",
        "wall_digest": "a0ee6805424226ade0ec4ee728689106cb4a44b18df8de19938449fad98171f9",
    },
    16384: {
        "D0": 400,
        "death": 4055,
        "fatal_bits": 9874,
        "profile": (10197, 12622, 19994),
        "first_negative": 1517,
        "q": 1516,
        "boundary_bits": (9870, 9874),
        "boundary_digest": "6ecf448b01a6225777b389d9293084e5c706c7cdf80fdfdf87f859b9e05d787f",
        "active_cells": 3221991,
        "active_digest": "de05e4c1cd3748884a8866af9fb274ae01b3d749c80bc63dcfe63806f73a87e6",
        "fatal_digest": "03c36f2c031cfcb59f07d49f7838f82c80905f4947a53d1ad754d4ef200bd254",
        "profile_digest": "e9bca4cece8530ea385422cc8958b83e5e2e60a25d4c44d93338d523941febed",
        "wall_digest": "3b443de3b7116591e37e6d0a2e8a96e4b716204028b190da4c860f741cc53461",
    },
}
EXPECTED_SEMANTIC_SHA256 = "08450a048b8e1379d3945e9464021b53203bb10f82638daf6bbfca076bd71a00"
BLOCK = 64
CHECKS = 0


def require(condition: bool, payload: object) -> None:
    """Optimization-safe unconditional gate."""
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
    observed_hash = lf_sha256(path)
    require(
        observed_hash == expected_hash,
        ("parent drift", label, observed_hash, expected_hash),
    )

sys.path.insert(0, str(ROOT / "04-computation"))
import amm_binary_sturmian_R8192_adjoint_horizon_thm3626 as P

A = P.A


def forward_step(
    state: tuple[int, ...], degree: int, delta: int, want_charges: bool
) -> tuple[tuple[int, ...], tuple[int, ...] | None]:
    """Apply one exact top-distance Rule-A step, optionally retaining charges."""
    require(delta in (0, 1), ("nonbinary degree step", degree, delta))
    width = len(state) - 1
    next_state = [0] * width
    charges = [0] * width if want_charges else None
    choose_r = 1
    for distance in range(width):
        if delta == 0:
            load = state[distance] + state[distance + 1]
        else:
            load = (
                (state[distance - 1] if distance else 0)
                + 2 * state[distance]
                + state[distance + 1]
            )
        choose_next = choose_r * (degree - 1 - distance) // (distance + 1)
        correction = min(choose_next, max(-choose_r, load))
        next_state[distance] = load - correction
        if charges is not None:
            charge = correction + choose_r
            if charge < 0:
                raise RuntimeError(("negative clamp charge", degree, distance, charge))
            charges[distance] = charge
        choose_r = choose_next
    return tuple(next_state), None if charges is None else tuple(charges)


def exact_profile_prefix(
    r_value: int, offset: int, death: int
) -> tuple[tuple[int, ...], int, str]:
    """Build exactly the mechanical word used through death plus its tail pin."""
    profile = tuple(
        A.golden_floor(r_value + row) + offset for row in range(death + 1)
    )
    require(
        all(
            profile[row] - profile[row - 1] in (0, 1)
            for row in range(1, len(profile))
        ),
        (r_value, "binary Sturmian prefix"),
    )
    last_degree = A.golden_floor(2 * r_value - 1) + offset
    profile_digest = digest_json((profile, last_degree))
    return profile, last_degree, profile_digest


def checkpointed_adjoint(
    r_value: int, profile: tuple[int, ...], death: int
) -> dict[str, object]:
    """Certify the first death and full sign wall with bounded storage."""
    initial_top_load = A.initial_halved_load(r_value, profile[0], profile[0])
    initial_top_residual = initial_top_load - min(1, max(0, initial_top_load))
    require(initial_top_residual == 0, (r_value, "row-zero fatal", initial_top_load))

    state = A.initial_top_state(r_value, profile[0], death)
    checkpoints: dict[int, tuple[int, ...]] = {0: state}
    boundaries = [0]
    for row in range(1, death):
        require(state[0] == 0, (r_value, "premature fatal row", row, state[0]))
        delta = profile[row] - profile[row - 1]
        state, _ = forward_step(state, profile[row], delta, False)
        if row % BLOCK == 0 or row == death - 1:
            checkpoints[row] = state
            boundaries.append(row)

    fatal = state[0]
    require(fatal < 0, (r_value, "fatal sign", fatal))

    level = (1,)
    rhs = fatal
    active_cells = 0
    multiplier_mass = 0
    maximum_multiplier = 0
    last_negative_cut: int | None = None
    last_negative_rhs: int | None = None
    last_negative_active: tuple[int, int, int] | None = None
    wall_hasher = sha256()

    for block_index in range(len(boundaries) - 1, 0, -1):
        start = boundaries[block_index - 1]
        end = boundaries[block_index]
        state = checkpoints[start]
        block_charges: list[tuple[int, ...]] = []
        for row in range(start + 1, end + 1):
            delta = profile[row] - profile[row - 1]
            state, charges = forward_step(state, profile[row], delta, True)
            require(charges is not None, (r_value, row, "missing charges"))
            block_charges.append(charges)

        for row, charges in zip(
            range(end, start, -1), reversed(block_charges), strict=True
        ):
            require(len(level) == len(charges), (r_value, row, "dual width"))
            contribution = sum(
                multiplier * charge
                for multiplier, charge in zip(level, charges, strict=True)
            )
            require(contribution >= 0, (r_value, row, "negative pairing"))
            rhs += contribution

            active_cells += len(level)
            multiplier_mass += sum(level)
            maximum_multiplier = max(maximum_multiplier, max(level))
            wall_hasher.update(f"{row}:{rhs}\n".encode("ascii"))

            if rhs < 0:
                last_negative_cut = row
                last_negative_rhs = rhs
                last_negative_active = (
                    active_cells,
                    multiplier_mass,
                    maximum_multiplier,
                )
            elif rhs > 0:
                require(
                    last_negative_cut == row + 1
                    and last_negative_rhs is not None
                    and last_negative_active is not None,
                    (r_value, row, "nonadjacent sign wall", last_negative_cut),
                )
                first_negative = last_negative_cut
                q_value = row
                boundary = (rhs, last_negative_rhs)
                active = last_negative_active
                ell = death - first_negative
                require(
                    active[0] == ell * (ell + 1) // 2,
                    (r_value, "triangular active count", active[0], ell),
                )
                return {
                    "fatal": fatal,
                    "fatal_digest": digest_json(fatal),
                    "first_negative": first_negative,
                    "q": q_value,
                    "negative_endpoints": (first_negative, death - 1),
                    "negative_count": death - first_negative,
                    "boundary": boundary,
                    "boundary_bits": tuple(abs(value).bit_length() for value in boundary),
                    "boundary_digest": digest_json(boundary),
                    "active": active,
                    "active_digest": digest_json(active),
                    "wall_digest": wall_hasher.hexdigest(),
                }
            else:
                raise RuntimeError((r_value, row, "zero adjoint boundary"))

            delta = profile[row] - profile[row - 1]
            level = A.previous_multiplier_level(level, delta)

    raise RuntimeError((r_value, "no positive/negative adjoint wall"))


def phase_invoice(
    r_value: int,
    d_zero: int,
    death: int,
    first_negative: int,
    q_value: int,
) -> dict[str, object]:
    """Return the exact integer/quadratic phase decomposition."""
    require(r_value % 8 == 0, (r_value, "phase denominator"))
    ell = death - first_negative
    h_value = 5 * r_value // 8 - d_zero
    edge_bound = 2 * d_zero - r_value + 2
    margin = death - edge_bound

    # Store a+b*sqrt(5) as (a,b).  Here
    # beta=(sqrt(5)-1)/8 and theta=(3-sqrt(5))/8.
    beta_r = (-r_value // 8, r_value // 8)
    ell_error = (ell - beta_r[0], -beta_r[1])
    q_error = (q_value - 3 * r_value // 8, r_value // 8)
    invoice_rhs = (
        margin + 1 - 2 * h_value - ell_error[0],
        -ell_error[1],
    )
    require(q_error == invoice_rhs, (r_value, "phase invoice", q_error, invoice_rhs))
    require(q_value == death - ell - 1, (r_value, "q/death/depth identity"))
    return {
        "h": h_value,
        "edge_bound": edge_bound,
        "margin": margin,
        "ell": ell,
        "ell_error": ell_error,
        "q_error": q_error,
    }


def certify_epoch(r_value: int) -> dict[str, object]:
    expected = EXPECTED[r_value]
    profile, last_degree, profile_digest = exact_profile_prefix(
        r_value, expected["D0"], expected["death"]
    )
    result = checkpointed_adjoint(r_value, profile, expected["death"])
    profile_control = (profile[0], profile[expected["death"]], last_degree)
    require(profile_control == expected["profile"], (r_value, "profile", profile_control))
    require(
        abs(result["fatal"]).bit_length() == expected["fatal_bits"],
        (r_value, "fatal bits"),
    )
    require(result["first_negative"] == expected["first_negative"], (r_value, "wall"))
    require(result["q"] == expected["q"], (r_value, "q"))
    require(
        result["negative_endpoints"]
        == (expected["first_negative"], expected["death"] - 1),
        (r_value, "negative endpoints"),
    )
    require(
        result["negative_count"] == expected["death"] - expected["first_negative"],
        (r_value, "negative count"),
    )
    require(result["boundary_bits"] == expected["boundary_bits"], (r_value, "boundary bits"))
    require(
        result["boundary_digest"] == expected["boundary_digest"],
        (r_value, "boundary digest", result["boundary_digest"]),
    )
    require(result["active"][0] == expected["active_cells"], (r_value, "active cells"))

    phase = phase_invoice(
        r_value,
        profile[0],
        expected["death"],
        result["first_negative"],
        result["q"],
    )
    actual_pins = {
        "fatal_digest": result["fatal_digest"],
        "profile_digest": profile_digest,
        "wall_digest": result["wall_digest"],
        "active_digest": result["active_digest"],
    }
    for label, observed in actual_pins.items():
        pinned = expected[label]
        if pinned != "TO_BE_PINNED":
            require(observed == pinned, (r_value, label, observed, pinned))

    return {
        "R": r_value,
        "D0": expected["D0"],
        "death": expected["death"],
        "profile": profile_control,
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


def main() -> None:
    control = certify_epoch(8192)
    target = certify_epoch(16384)

    require(control["q"] == P.EXPECTED["q"], "THM3626 q control")
    require(control["boundary_digest"] == P.EXPECTED["boundary_digest"], "THM3626 boundary control")
    require(control["active_digest"] == P.EXPECTED["active_digest"], "THM3626 active control")

    # Both q/R values lie below theta=(3-sqrt(5))/8.
    require(2299 * 2299 - 5 * 1024 * 1024 == 42521, "R8192 radical square gap")
    require(1157 * 1157 - 5 * 512 * 512 == 27929, "R16384 radical square gap")
    require(2 * control["q"] - target["q"] == 30, "normalized q deterioration")
    require(
        target["death"] - 2 * control["death"] == -35,
        "dyadic death defect",
    )
    require(
        target["phase"]["ell"] - 2 * control["phase"]["ell"] == -4,
        "dyadic depth defect",
    )
    require(
        target["q"] - 2 * control["q"] == -30,
        "dyadic q defect",
    )

    semantic_record = {
        "parents": tuple((label, digest) for label, _, digest in PARENT_FILES),
        "control": control,
        "target": target,
        "radical_square_gaps": (42521, 27929),
        "dyadic_defects": (-30, -35, -4),
        "caveat": "R16384 D0=400 is failing but not known terminal",
    }
    semantic = digest_json(semantic_record)
    if EXPECTED_SEMANTIC_SHA256 != "TO_BE_PINNED":
        require(
            semantic == EXPECTED_SEMANTIC_SHA256,
            ("semantic digest", semantic, EXPECTED_SEMANTIC_SHA256),
        )

    source_path = Path(__file__).resolve()
    source_bytes = source_path.read_bytes()
    require(b"\r\n" not in source_bytes, "source is not raw LF")
    require(
        not any(
            isinstance(node, ast.Assert)
            for node in ast.walk(ast.parse(source_bytes.decode("utf-8")))
        ),
        "Python assert node present",
    )

    print("== AMM fixed R16384,D0=400 checkpointed adjoint sidecar ==")
    print("parent_sha256_lf", tuple((label, digest) for label, _, digest in PARENT_FILES))
    for label, result in (("R8192_control", control), ("R16384_target", target)):
        print(
            label,
            {
                "epoch": (result["R"], result["D0"]),
                "death": result["death"],
                "profile": result["profile"],
                "profile_digest": result["profile_digest"],
                "fatal_bits": result["fatal_bits"],
                "fatal_digest": result["fatal_digest"],
                "wall": (result["first_negative"], result["q"], result["negative_endpoints"]),
                "wall_digest": result["wall_digest"],
                "boundary": (result["boundary_bits"], result["boundary_digest"]),
                "active": (result["active_cells"], result["active_digest"]),
                "phase": result["phase"],
            },
        )
    print("golden_error_R8192", "(-2299+1024*sqrt(5))/8192<0|square_gap:42521")
    print("golden_error_R16384", "(-1157+512*sqrt(5))/4096<0|square_gap:27929")
    print("error_deterioration", "E16384-E8192=15/8192>0")
    print("dyadic_defects", "q:-30|death:-35|depth:-4")
    print("semantic_sha256", semantic)
    print("script_sha256_lf", lf_sha256(source_path))
    print("CHECKS", CHECKS)
    print("SCOPE", "fixed failed R16384,D0=400 prefix; D0=400 not known terminal; no AMM bound")
    print("PASS")


if __name__ == "__main__":
    main()
