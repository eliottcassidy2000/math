#!/usr/bin/env python3
"""Exact postcomposition sidecar for THM-2915 and later LRC root banks.

THM-2916 was promoted while the broader THM-2915 computation was in its
final replay; THM-2920 arrived during the clean-checkout replay, and
THM-2919's sharp-H1 promotion arrived during final rebasing.  This verifier
compares all locked explicit root lists.  It is deliberately a composition
sidecar: THM-2915's child computation does not depend on later theorems.
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
THM2920_OUT = (
    ROOT
    / "05-knowledge/results/"
    "lrc14_j6_two_h3_row_pair_hunter_recursive_toothpick_closure_codex_20260729.out"
)
THM2919_OUT = (
    ROOT
    / "05-knowledge/results/"
    "lrc14_j6_sharp_h1_tail_recomposition_codex_20260729.out"
)

EXPECTED_SHA256 = {
    THM2915_OUT: "26d1b492af588dda57f0531b9bb1acd5faca32b12d7cbac195fea31b3d4dd30e",
    THM2916_OUT: "6b31abaadd4f089a9f98a9eea49845c0ed0123a810bb4bf4f3c0309c519e7df9",
    THM2920_OUT: "1a38fd441dfd77a4f5d30d45d3160febc33d2d4eeb6247b223f10a1e31a8aefb",
    THM2919_OUT: "8eb5b8af65466539f070418b1d405b6b8b7c8d7f9595c2a7fa069196f80c0d38",
}
EXPECTED_THM2916_DIGEST = (
    "c68d09676683f6204df3b04353a3b3107ebbb4285d13a3b6001446372e351e1b"
)
EXPECTED_THM2920_DIGEST = (
    "e3045198e08804c78025bd532111377309882911e08bc50604aa7119ac266c71"
)
EXPECTED_COMPLETE_TWO_H3_DIGEST = (
    "772f67e5711cb009012a9c4abeb1b9a288195126f382244231f1f93362b63efc"
)
EXPECTED_THM2915_TUPLE_DIGEST = (
    "d6702800de9451cce73801cb0cb9d8566f9e42ad4e63c6b0e77c9668b37cdce9"
)
EXPECTED_THM2920_NEW_DIGEST = (
    "6c4556f9fd07fafa1fe7f62a4088b7842bd08de32b188eab231313f748f77703"
)
EXPECTED_CURRENT_UNION_DIGEST = (
    "675631dfe7a4bf5924b7f16a100e7eebe60901d25190e493fa4a640325475a50"
)
EXPECTED_REPAIR_OVERLAP = {(1, 2, 3, 6, 9, 11, 12)}
EXPECTED_SHARP_ROOTS = {
    (1, 3, 5, 7, 10, 11, 14),
    (2, 4, 8, 10, 11, 13, 14),
}


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
        if path != THM2919_OUT:
            require(
                unique_text(path, "mode=") == "LOCKED",
                f"{path.name}: unlocked",
            )
        require(
            unique_text(path, "all_exact_controls=") == "PASS",
            f"{path.name}: exact controls failed",
        )

    broad = set(literal(THM2915_OUT, "current_union_roots="))
    narrow = set(literal(THM2916_OUT, "closed_roots="))
    repairs = set(literal(THM2920_OUT, "closed_roots="))
    sharp_roots = set(literal(THM2919_OUT, "additive_route_roots="))
    sharp_current = set(literal(THM2919_OUT, "additive_current_roots="))
    complete_two_h3 = set(
        literal(THM2920_OUT, "complete_two_h3_roots=")
    )
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
    require(
        unique_text(THM2920_OUT, "closed_root_digest=")
        == EXPECTED_THM2920_DIGEST
        and tuple_digest(repairs) == EXPECTED_THM2920_DIGEST
        and unique_text(THM2920_OUT, "complete_two_h3_digest=")
        == EXPECTED_COMPLETE_TWO_H3_DIGEST
        and tuple_digest(complete_two_h3)
        == EXPECTED_COMPLETE_TWO_H3_DIGEST,
        "THM2920 root digests changed",
    )

    narrow_overlap = broad & narrow
    narrow_new = narrow - broad
    repair_overlap = broad & repairs
    repair_new = repairs - broad
    current = broad | repairs
    h_synergies = set(literal(THM2915_OUT, "H_additive_gain_roots="))
    require(
        len(broad) == 1_315
        and len(narrow) == 394
        and narrow_overlap == narrow
        and not narrow_new
        and tuple_digest(broad) == EXPECTED_THM2915_TUPLE_DIGEST
        and len(repairs) == 296
        and not (narrow & repairs)
        and narrow | repairs == complete_two_h3
        and repair_overlap == EXPECTED_REPAIR_OVERLAP
        and repair_overlap <= h_synergies
        and len(repair_new) == 295
        and tuple_digest(repair_new) == EXPECTED_THM2920_NEW_DIGEST
        and sharp_roots == EXPECTED_SHARP_ROOTS
        and len(sharp_current) == 1
        and sharp_current <= sharp_roots
        and sharp_roots <= broad
        and len(current) == 1_610
        and tuple_digest(current) == EXPECTED_CURRENT_UNION_DIGEST,
        "post-THM2915 root composition changed",
    )

    print("LRC14 THM2915 postcomposition with THM2916/2919/2920")
    print(
        "counts="
        f"(thm2915_union={len(broad)},thm2916_closed={len(narrow)},"
        f"thm2916_outside={len(narrow_new)},thm2920_repairs={len(repairs)},"
        f"thm2920_overlap={len(repair_overlap)},"
        f"thm2920_new={len(repair_new)},thm2919_route={len(sharp_roots)},"
        f"thm2919_outside={len(sharp_roots-broad)},current_union={len(current)},"
        f"current_residual={3432-len(current)})"
    )
    print(f"thm2916_overlap_digest={tuple_digest(narrow_overlap)}")
    print(f"thm2920_overlap_roots={tuple(sorted(repair_overlap))}")
    print(f"thm2920_new_digest={tuple_digest(repair_new)}")
    print(f"current_union_digest={tuple_digest(current)}")
    print("mode=LOCKED")
    print(
        "scope=exact root-set composition only;"
        "THM2915 remains independent of THM2916/2919/2920;not LRC14"
    )
    print("all_exact_controls=PASS")


if __name__ == "__main__":
    main()
