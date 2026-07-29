#!/usr/bin/env python3
"""Independent exact audit of the sole depth-three equality branch.

This verifier reads only the emitted third-level ledger, reconstructs
the terminal 13-speed family from its branch key and pair witness, checks the
strict/open endpoint convention directly, and joins the family to THM-2907's
independent endpoint bank.
"""

from __future__ import annotations

import argparse
import ast
import hashlib
import importlib.util
from fractions import Fraction as F
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
BASE_PATH = Path(__file__).with_name(
    "lrc14_j6_seven_body_six_slot_recursive_pair_hunter_thm2923_stage1.py"
)
EXPECTED_BASE_SHA256 = (
    "74e5eb2d0b23dda6366d115376890726ac32dc834a09a71479ea875f05e7615e"
)
DEFAULT_THM2907_ENDPOINT_OUTPUT = (
    ROOT
    / "05-knowledge/results/"
    "lrc14_j6_paircap_exception_endpoint_convention_audit_codex_20260729.out"
)
TARGET = F(1, 14)
TIME = F(1, 28)
EPSILON = F(1, 10**9)
EXPECTED_BODY = (2, 4, 6, 8, 10, 12, 14)
EXPECTED_FAMILY = tuple(range(2, 27, 2))
EXPECTED_PRIMITIVE = tuple(range(1, 14))
EXPECTED_THM2907_ENDPOINT_OUTPUT_SHA256 = (
    "a93a4a724dac6c55806f3358c2f5ab25de8f0261c92906a0161414781a717d20"
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def load_base():
    require(file_sha256(BASE_PATH) == EXPECTED_BASE_SHA256, "stage-one source changed")
    spec = importlib.util.spec_from_file_location("endpoint_equality_base", BASE_PATH)
    require(spec is not None and spec.loader is not None, "cannot load interval base")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


G = load_base()


def fields(line: str) -> dict[str, str]:
    return dict(part.split("=", 1) for part in line.rstrip().split(";")[1:])


def ints(text: str) -> tuple[int, ...]:
    return () if not text else tuple(map(int, text.split(",")))


def distance(speed: int, time: F) -> F:
    residue = (speed * time) % 1
    return min(residue, 1 - residue)


def ftext(value: F) -> str:
    return f"{value.numerator}/{value.denominator}"


def parse_endpoint_bank(path: Path) -> set[tuple[int, ...]]:
    require(
        file_sha256(path) == EXPECTED_THM2907_ENDPOINT_OUTPUT_SHA256,
        "THM-2907 endpoint output changed",
    )
    families: set[tuple[int, ...]] = set()
    for line in path.read_text().splitlines():
        if not line.startswith("family="):
            continue
        family_text = line.split(";", 1)[0].removeprefix("family=")
        family = ast.literal_eval(family_text)
        require(
            isinstance(family, tuple) and all(isinstance(v, int) for v in family),
            "malformed THM-2907 endpoint family",
        )
        families.add(family)
    require(len(families) == 2, "THM-2907 endpoint bank changed")
    return families


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--third-ledger", type=Path, required=True)
    parser.add_argument(
        "--thm2907-endpoint-output",
        type=Path,
        default=DEFAULT_THM2907_ENDPOINT_OUTPUT,
    )
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    third_count = 0
    open_rows: list[dict[str, str]] = []
    for line in args.third_ledger.read_text().splitlines():
        if not line.startswith("THIRD;"):
            continue
        third_count += 1
        data = fields(line)
        if data["closed"] == "0":
            open_rows.append(data)
    require(third_count == 117, "third-pivot universe changed")
    require(len(open_rows) == 1, "sole third-pivot equality changed")
    row = open_rows[0]

    body = ints(row["E"])
    rank = int(row["rank"])
    apex = int(row["a"])
    prefix = ints(row["P"])
    first = int(row["x"])
    first_earlier = ints(row["earlier"])
    second = int(row["y"])
    second_earlier = ints(row["yearlier"])
    third = int(row["z"])
    third_earlier = ints(row["zearlier"])
    witness = ints(row["B2w"])
    margin = F(row["margin"])
    require(
        (
            body,
            rank,
            apex,
            prefix,
            first,
            first_earlier,
            second,
            second_earlier,
            third,
            third_earlier,
            witness,
            margin,
        )
        == (
            EXPECTED_BODY,
            1,
            22,
            (22,),
            18,
            (),
            26,
            (),
            16,
            (),
            (20, 24),
            F(0),
        ),
        "terminal equality branch identity changed",
    )

    family = tuple(sorted((*body, apex, first, second, third, *witness)))
    require(
        len(family) == 13
        and len(set(family)) == 13
        and family == EXPECTED_FAMILY,
        "terminal family is not 2*{1,...,13}",
    )
    primitive = tuple(value // 2 for value in family)
    require(
        all(value % 2 == 0 for value in family)
        and primitive == EXPECTED_PRIMITIVE,
        "terminal scale-down is not {1,...,13}",
    )

    carrier, components, mass = G.direct_carrier(family)
    require(
        carrier == [] and components == 0 and mass == 0,
        "closed-comb engine did not erase the endpoint family",
    )
    distances = {value: distance(value, TIME) for value in family}
    strict = tuple(value for value in family if distances[value] < TARGET)
    boundary = tuple(value for value in family if distances[value] == TARGET)
    require(
        strict == ()
        and boundary == (2, 26)
        and min(distances.values()) == TARGET
        and all(value >= TARGET for value in distances.values()),
        "strict/open endpoint audit failed",
    )
    primitive_residues = tuple((value * TIME * 2) % 1 for value in primitive)
    require(
        primitive_residues == tuple(F(value, 14) for value in range(1, 14)),
        "primitive residue proof changed",
    )
    left = min(distance(value, TIME - EPSILON) for value in family)
    right = min(distance(value, TIME + EPSILON) for value in family)
    require(F(0) < EPSILON < F(1, 364), "perturbation left the local endpoint cell")
    require(
        left == TARGET - 2 * EPSILON
        and right == TARGET - 26 * EPSILON,
        "two-sided local endpoint formulas changed",
    )

    endpoint_bank = parse_endpoint_bank(args.thm2907_endpoint_output)
    require(family in endpoint_bank, "terminal family is absent from THM-2907")

    lines = [
        f"body={body}",
        f"branch=(rank={rank},apex={apex},prefix={prefix},x={first},y={second},z={third})",
        f"pair_witness={witness}",
        f"family={family}",
        f"primitive={primitive}",
        f"time={ftext(TIME)}",
        f"target={ftext(TARGET)}",
        f"strict_danger={strict}",
        f"boundary={boundary}",
        f"closed_carrier=(components={components},mass={ftext(mass)})",
        f"perturbation_clearance=(left={ftext(left)},right={ftext(right)})",
        "thm2907_endpoint_join=PASS",
    ]
    digest = hashlib.sha256(
        (
            "LRC14/THM2915-failure-recursion/endpoint-equality/v1\n"
            + "\n".join(lines)
            + "\n"
        ).encode()
    ).hexdigest()
    rendered = (
        "LRC14 THM2915-failure recursion endpoint equality audit\n"
        + "\n".join(lines)
        + "\n"
        + f"semantic_sha256={digest}\n"
        + "mode=LOCKED\n"
        + "all_exact_controls=PASS\n"
    )
    if args.output is not None:
        args.output.write_text(rendered, encoding="utf-8", newline="\n")
    print(rendered, end="")


if __name__ == "__main__":
    main()
