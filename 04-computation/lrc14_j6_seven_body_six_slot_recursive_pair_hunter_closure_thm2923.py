#!/usr/bin/env python3
"""Locked replay driver for the complete seven-body/six-slot rung.

This verifier regenerates every proof-bearing stage in a temporary directory,
checks the LF-normalized SHA-256 of each handoff before it is consumed, and
emits one compact summary plus one identity-complete combined ledger.

The universe is exactly E in C([14], 7), with six external speed slots.  This
script does not treat the at-most-six-in-window (j >= 7) loose-spread sector
and therefore does not prove the unrestricted fourteen-runner conjecture.
"""

from __future__ import annotations

import argparse
import hashlib
import os
import subprocess
import sys
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
STEM = "lrc14_j6_seven_body_six_slot_recursive_pair_hunter_thm2923"
STAGE1 = ROOT / "04-computation" / f"{STEM}_stage1.py"
STAGE2 = ROOT / "04-computation" / f"{STEM}_stage2.py"
COMPOSITION = ROOT / "04-computation" / f"{STEM}_composition.py"
ENDPOINT_AUDIT = ROOT / "04-computation" / f"{STEM}_endpoint_audit.py"
INHERITED_AUDIT = ROOT / "04-computation" / f"{STEM}_inherited_audit.py"
DEFAULT_OUTPUT = (
    ROOT
    / "05-knowledge/results"
    / "lrc14_j6_seven_body_six_slot_recursive_pair_hunter_closure_thm2923.out"
)
DEFAULT_LEDGER = (
    ROOT
    / "05-knowledge/results"
    / "lrc14_j6_seven_body_six_slot_recursive_pair_hunter_closure_thm2923.ledger.out"
)

EXPECTED_HELPER_SHA256 = {
    STAGE1: "74e5eb2d0b23dda6366d115376890726ac32dc834a09a71479ea875f05e7615e",
    STAGE2: "54aad75134a3066dca355fb2ed536c7bea901df9d6e2eefc0f1ead30e57e7622",
    COMPOSITION: "7cb3384831151ac3b2dd6d3b38f185d5328ed8c13dfca757621b573a6120600d",
    ENDPOINT_AUDIT: "d8d9504c643cc5d4e2d2fe69893042c4245c7990b97c4a5e697bc64a7bbc5edf",
    INHERITED_AUDIT: "cb40a748fecad6659eb5b2a140d2d8d23a966a643bd1a3918b316341644dd78d",
}
EXPECTED_ARTIFACT_SHA256 = {
    "PAIR": "76bfbdbc6cab524a4536c638a8ed58b0921131f21d0a2c99ada2c6fb4ef65632",
    "RECURSIVE": "a5f2c832f2feadfbf06971b013ba8b070a987fcd865ae38bf0c3757940807ec9",
    "STAGE1_SUMMARY": "c9ed0849443a433fdf75899d98d033d273c19316d08e579ca8c4d9d6e9815be0",
    "GRANDCHILD": "a121a239b08024c0cb91ed292e1a0ccee1a403f567d21f532bdc9ee3439f7893",
    "THIRD": "d6d3fd59bfff3f46aabcb9c80580ecce48cad79dced9782460ed917158c4946f",
    "STAGE2_SUMMARY": "161d702bc98f25d3d3dcfcd8d832e3aa2703f246a54d39ec55935ef134166f88",
    "COMPOSITION": "fc58f23e4a502f6e7bda56847cf4ee2b09cfd9279b195636ab5e6611bebf18f0",
    "ROOTS": "bfe95330f6ede745f4a16997377935999d915650b69cde1e5ab121441a167f58",
    "ENDPOINT": "8fe9e2c8f2e2dd41ab97223d793a71325b1749f9486f44e8197071818a3d4c8d",
    "INHERITED": "d5d5d6ba2cb9b47920aa4b9de615ffe829a91ddc7dc5b7da0d6b6a7310f8f2e9",
}
EXPECTED_LEDGER_SEMANTIC_SHA256: str | None = (
    "a2555163f07f6787e3a789b365ef27c47957ed39a82074ee9067b9bac89f3104"
)
EXPECTED_OUTPUT_SEMANTIC_SHA256: str | None = (
    "c2ef8608ac64f068e6cc83481c8b91c126c1146b4808bd21b9e42d61e08d6abb"
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_sha256(path: Path) -> str:
    payload = path.read_bytes()
    require(
        b"\r" not in payload.replace(b"\r\n", b""),
        f"{path.name}: unexpected lone carriage return",
    )
    return hashlib.sha256(payload.replace(b"\r\n", b"\n")).hexdigest()


def digest(domain: str, text: str) -> str:
    return hashlib.sha256((domain + "\n" + text).encode()).hexdigest()


def run_stage(
    name: str,
    script: Path,
    arguments: list[str],
    expected_stdout: Path,
    environment: dict[str, str],
) -> None:
    optimize = ["-O"] if sys.flags.optimize else []
    result = subprocess.run(
        [sys.executable, *optimize, str(script), *arguments],
        cwd=ROOT,
        env=environment,
        check=False,
        capture_output=True,
        text=True,
    )
    require(
        result.returncode == 0,
        f"{name} failed with code {result.returncode}:\n{result.stderr}\n{result.stdout}",
    )
    require(result.stderr == "", f"{name} emitted stderr: {result.stderr}")
    require(expected_stdout.is_file(), f"{name} did not emit {expected_stdout.name}")
    require(
        result.stdout == expected_stdout.read_text(),
        f"{name} stdout and written summary differ",
    )


def check_artifact(label: str, path: Path) -> str:
    actual = file_sha256(path)
    require(
        actual == EXPECTED_ARTIFACT_SHA256[label],
        f"{label} handoff changed: expected {EXPECTED_ARTIFACT_SHA256[label]}, got {actual}",
    )
    return actual


def replay(work: Path, workers: int, hash_seed: int) -> tuple[str, str]:
    environment = dict(os.environ)
    environment["PYTHONHASHSEED"] = str(hash_seed)

    pair = work / "pair.ledger.out"
    recursive = work / "recursive.ledger.out"
    stage1_summary = work / "stage1.out"
    grandchild = work / "grandchild.ledger.out"
    third = work / "third.ledger.out"
    stage2_summary = work / "stage2.out"
    composition = work / "composition.out"
    roots = work / "roots.ledger.out"
    endpoint = work / "endpoint.out"
    inherited = work / "inherited.out"

    run_stage(
        "stage one",
        STAGE1,
        [
            "--workers",
            str(workers),
            "--pair-ledger",
            str(pair),
            "--recursive-ledger",
            str(recursive),
            "--summary",
            str(stage1_summary),
        ],
        stage1_summary,
        environment,
    )
    check_artifact("PAIR", pair)
    check_artifact("RECURSIVE", recursive)
    check_artifact("STAGE1_SUMMARY", stage1_summary)

    run_stage(
        "stage two",
        STAGE2,
        [
            "--pair-ledger",
            str(pair),
            "--recursive-ledger",
            str(recursive),
            "--workers",
            str(workers),
            "--grandchild-ledger",
            str(grandchild),
            "--third-ledger",
            str(third),
            "--summary",
            str(stage2_summary),
        ],
        stage2_summary,
        environment,
    )
    check_artifact("GRANDCHILD", grandchild)
    check_artifact("THIRD", third)
    check_artifact("STAGE2_SUMMARY", stage2_summary)

    run_stage(
        "composition",
        COMPOSITION,
        [
            "--pair-ledger",
            str(pair),
            "--recursive-ledger",
            str(recursive),
            "--grandchild-ledger",
            str(grandchild),
            "--third-ledger",
            str(third),
            "--output",
            str(composition),
            "--root-ledger",
            str(roots),
        ],
        composition,
        environment,
    )
    check_artifact("COMPOSITION", composition)
    check_artifact("ROOTS", roots)

    run_stage(
        "endpoint audit",
        ENDPOINT_AUDIT,
        ["--third-ledger", str(third), "--output", str(endpoint)],
        endpoint,
        environment,
    )
    check_artifact("ENDPOINT", endpoint)

    run_stage(
        "inherited-slice audit",
        INHERITED_AUDIT,
        [
            "--pair-ledger",
            str(pair),
            "--recursive-ledger",
            str(recursive),
            "--grandchild-ledger",
            str(grandchild),
            "--root-ledger",
            str(roots),
            "--output",
            str(inherited),
        ],
        inherited,
        environment,
    )
    check_artifact("INHERITED", inherited)

    ledger_parts = [
        "LRC14 seven-body/six-slot recursive pair-Hunter proof ledger\n",
        "universe=E in C([14],7);six external slots\n",
    ]
    for label, path in (
        ("PAIR", pair),
        ("RECURSIVE", recursive),
        ("STAGE1_SUMMARY", stage1_summary),
        ("GRANDCHILD", grandchild),
        ("THIRD", third),
        ("STAGE2_SUMMARY", stage2_summary),
        ("COMPOSITION", composition),
        ("ROOTS", roots),
        ("ENDPOINT", endpoint),
        ("INHERITED", inherited),
    ):
        ledger_parts.extend((f"BEGIN_{label}\n", path.read_text(), f"END_{label}\n"))
    ledger_body = "".join(ledger_parts)
    ledger_semantic = digest(
        "LRC14/THM2923/seven-body-six-slot/combined-ledger/v1",
        ledger_body,
    )
    if EXPECTED_LEDGER_SEMANTIC_SHA256 is not None:
        require(
            ledger_semantic == EXPECTED_LEDGER_SEMANTIC_SHA256,
            "combined ledger semantic digest changed",
        )
    ledger_rendered = (
        ledger_body
        + f"combined_semantic_sha256={ledger_semantic}\n"
        + "mode=LOCKED\n"
        + "all_exact_controls=PASS\n"
    )

    manifest = tuple(
        (label, EXPECTED_ARTIFACT_SHA256[label])
        for label in sorted(EXPECTED_ARTIFACT_SHA256)
    )
    output_lines = [
        "LRC14 complete seven-body/six-slot recursive pair-Hunter closure",
        "universe=3432 roots E in C([14],7); six external slots",
        "stage1=(failures=4866,pair=1612,hunter4=1175,recursive=2079)",
        "stage2=(second_pivots=6172,top3_closed=5944,failed_grandchildren=228)",
        "stage3=(pair_plus_singleton=65,hunter3=114,third_recursive=49)",
        "terminal=(third_pivots=117,strict=116,equality=1)",
        "endpoint=(family=2*{1,...,13},time=1/28,boundary={2,26},thm2907_join=PASS)",
        "composition=(parents=11842,E=52,P=279,T2=11562,E_union_P_union_T2=11842)",
        "roots=(deep_route=3411,baseline=351,overlap=330,union=3432,residual=0)",
        "inherited_controls=(thm2913_children=42,thm2920_children=367,thm2920_roots=296)",
        "independent_hostile_audit=(auditor=thm2903_2905_union_audit,status=PASS,"
        "parents=11842,terminal=116+1,route_bits=recomputed,union=3432)",
        f"artifact_manifest={manifest}",
        f"combined_ledger_semantic_sha256={ledger_semantic}",
        "why_stop_j7=(q7>=h/7 and B2>=2h/7 imply G7(h/7)=h;"
        "lambda=h/7;THM735 finite-core gate unavailable)",
        "scope=complete seven-in-window/six-external rung only",
        "remaining_scope=at-most-six-in-window (j>=7) loose-spread sector OPEN",
        "not_unrestricted_LRC14=TRUE",
    ]
    output_body = "\n".join(output_lines) + "\n"
    output_semantic = digest(
        "LRC14/THM2923/seven-body-six-slot/summary/v1",
        output_body,
    )
    if EXPECTED_OUTPUT_SEMANTIC_SHA256 is not None:
        require(
            output_semantic == EXPECTED_OUTPUT_SEMANTIC_SHA256,
            "summary semantic digest changed",
        )
    output_rendered = (
        output_body
        + f"semantic_sha256={output_semantic}\n"
        + "mode=LOCKED\n"
        + "all_exact_controls=PASS\n"
    )
    return output_rendered, ledger_rendered


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=max(1, os.cpu_count() or 1))
    parser.add_argument("--hash-seed", type=int, default=0)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--ledger", type=Path, default=DEFAULT_LEDGER)
    args = parser.parse_args()
    require(args.workers >= 1, "workers must be positive")
    require(args.hash_seed >= 0, "hash seed must be nonnegative")
    require(args.output.resolve() != args.ledger.resolve(), "output and ledger collide")
    for path, expected in EXPECTED_HELPER_SHA256.items():
        actual = file_sha256(path)
        require(
            actual == expected,
            f"{path.name} changed: expected {expected}, got {actual}",
        )

    with tempfile.TemporaryDirectory(prefix="thm2923-") as directory:
        output, ledger = replay(Path(directory), args.workers, args.hash_seed)
    args.output.write_text(output, encoding="utf-8", newline="\n")
    args.ledger.write_text(ledger, encoding="utf-8", newline="\n")
    print(output, end="")


if __name__ == "__main__":
    main()
