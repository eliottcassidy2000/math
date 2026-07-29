#!/usr/bin/env python3
"""Exact root-set containment sidecar for THM-2915 and THM-2916.

THM-2916 was promoted while the broader THM-2915 computation was in its
final replay.  This verifier compares the two locked explicit root lists.
It is deliberately a composition sidecar: THM-2915's child computation
does not depend on the later, narrower theorem.
"""

from __future__ import annotations

import ast
import hashlib
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
THM2915_OUT = (
    ROOT
    / "05-knowledge/results/"
    "lrc14_j6_all_open_centre_child_top4_closure_thm2915.out"
)
THM2916_OUT = (
    ROOT
    / "05-knowledge/results/"
    "lrc14_j6_two_h3_row_child_top4_census_codex_20260729.out"
)

EXPECTED_SHA256 = {
    THM2915_OUT: "26d1b492af588dda57f0531b9bb1acd5faca32b12d7cbac195fea31b3d4dd30e",
    THM2916_OUT: "6b31abaadd4f089a9f98a9eea49845c0ed0123a810bb4bf4f3c0309c519e7df9",
}
EXPECTED_THM2916_DIGEST = (
    "c68d09676683f6204df3b04353a3b3107ebbb4285d13a3b6001446372e351e1b"
)
EXPECTED_JOINT_DIGEST = (
    "d6702800de9451cce73801cb0cb9d8566f9e42ad4e63c6b0e77c9668b37cdce9"
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def unique_text(path: Path, prefix: str) -> str:
    matches = [
        line.removeprefix(prefix)
        for line in path.read_text().splitlines()
        if line.startswith(prefix)
    ]
    require(len(matches) == 1, f"{path.name}: expected one {prefix!r} line")
    return matches[0]


def literal(path: Path, prefix: str) -> object:
    return ast.literal_eval(unique_text(path, prefix))


def tuple_digest(roots: set[tuple[int, ...]]) -> str:
    return hashlib.sha256(repr(tuple(sorted(roots))).encode()).hexdigest()


def main() -> None:
    for path, expected in EXPECTED_SHA256.items():
        require(file_sha256(path) == expected, f"{path.name}: hash changed")
        require(unique_text(path, "mode=") == "LOCKED", f"{path.name}: unlocked")
        require(
            unique_text(path, "all_exact_controls=") == "PASS",
            f"{path.name}: exact controls failed",
        )

    broad = set(literal(THM2915_OUT, "current_union_roots="))
    narrow = set(literal(THM2916_OUT, "closed_roots="))
    root_digests = literal(THM2915_OUT, "root_digests=")
    require(
        isinstance(root_digests, dict)
        and root_digests["current_union1315"]
        == "47c0b23646ae5744f5354d6475aa23283c275341e6665ba5532640cddea0c41f",
        "THM2915 current-union digest changed",
    )
    require(
        unique_text(THM2916_OUT, "closed_root_digest=")
        == EXPECTED_THM2916_DIGEST
        and tuple_digest(narrow) == EXPECTED_THM2916_DIGEST,
        "THM2916 closed-root digest changed",
    )

    overlap = broad & narrow
    new_from_narrow = narrow - broad
    joint = broad | narrow
    require(
        len(broad) == 1_315
        and len(narrow) == 394
        and overlap == narrow
        and not new_from_narrow
        and len(joint) == 1_315
        and tuple_digest(joint) == EXPECTED_JOINT_DIGEST,
        "THM2916 is no longer contained in the THM2915 union",
    )

    print("LRC14 THM2915/THM2916 exact root containment")
    print(
        "counts="
        f"(thm2915_union={len(broad)},thm2916_closed={len(narrow)},"
        f"overlap={len(overlap)},new_from_thm2916={len(new_from_narrow)},"
        f"joint_union={len(joint)},residual={3432-len(joint)})"
    )
    print(f"thm2916_overlap_digest={tuple_digest(overlap)}")
    print(f"joint_union_digest={tuple_digest(joint)}")
    print("mode=LOCKED")
    print(
        "scope=exact root-set composition only;"
        "THM2915 remains independent of THM2916;not LRC14"
    )
    print("all_exact_controls=PASS")


if __name__ == "__main__":
    main()
