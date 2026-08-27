#!/usr/bin/env python3
"""Independent semantic reconstruction of the current post-THM-4254 residual.

The inherited THM-4231/4238/4242 compiler is used only for its frozen
``aggregate_residual`` consequence object.  This script separately removes
the three THM-4252 edges and every inherited edge in the semantic THM-4254
band 755 <= second endpoint <= 768, verifies exact fingerprints, and either
prints a compact audit report or the complete ``q,r`` residual universe.
"""

from __future__ import annotations

import argparse
import contextlib
import hashlib
import io
import os
from pathlib import Path
import runpy


FNV_OFFSET = 0xCBF29CE484222325
FNV_PRIME = 0x100000001B3
MASK64 = (1 << 64) - 1


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def fnv_add(state: int, word: int) -> int:
    for shift in range(0, 64, 8):
        state ^= (word >> shift) & 0xFF
        state = (state * FNV_PRIME) & MASK64
    return state


def edge_fnv(edges: list[tuple[int, int]]) -> int:
    state = FNV_OFFSET
    for q, r in edges:
        state = fnv_add(state, q)
        state = fnv_add(state, r)
    return state


def edge_sha(edges: list[tuple[int, int]]) -> str:
    payload = b"".join(f"{q},{r}\n".encode("ascii") for q, r in edges)
    return hashlib.sha256(payload).hexdigest()


def reconstruct(repo: Path) -> tuple[
    list[tuple[int, int]],
    list[tuple[int, int]],
    list[tuple[int, int]],
]:
    inherited_path = repo / (
        "04-computation/"
        "lrc14_thm4231_4238_4242_cross_residual_postprocess_20260826.py"
    )
    require(inherited_path.is_file(), "missing inherited residual compiler")
    prior_cwd = Path.cwd()
    try:
        os.chdir(repo)
        with contextlib.redirect_stdout(io.StringIO()):
            namespace = runpy.run_path(str(inherited_path))
    finally:
        os.chdir(prior_cwd)
    inherited = list(namespace["aggregate_residual"])
    require(len(inherited) == 181_126, "inherited count changed")
    require(edge_fnv(inherited) == 0xBDF59726990A6C92,
            "inherited FNV changed")
    require(edge_sha(inherited) ==
            "c0e2fe1c69cfe8cfe6e633a1eca0d8d37ca991ecdaa04b98d7c595a99b9be6bf",
            "inherited SHA changed")

    removed_4252 = {(466, 699), (616, 769), (721, 769)}
    require(removed_4252 <= set(inherited), "THM-4252 edge missing")
    post_4252 = [edge for edge in inherited if edge not in removed_4252]
    require(len(post_4252) == 181_123, "post-THM-4252 count changed")
    require(edge_fnv(post_4252) == 0x6EC03ED4C4DC841B,
            "post-THM-4252 FNV changed")
    require(edge_sha(post_4252) ==
            "9a9b6fbe14db00e9d7f8f08ecddaa1e3d263fd063c6b3c003e18c210b3334ef8",
            "post-THM-4252 SHA changed")

    band = [edge for edge in post_4252 if 755 <= edge[1] <= 768]
    require(len(band) == 59, "semantic THM-4254 band count changed")
    require(edge_fnv(band) == 0xB3D54B78BABBCAEC,
            "semantic THM-4254 band FNV changed")
    require(edge_sha(band) ==
            "6b54d8fa3b408325fc309bec3ed769f5e56ce370fa34fa7ad1bb6d7ed4cafc36",
            "semantic THM-4254 band SHA changed")

    band_set = set(band)
    residual = [edge for edge in post_4252 if edge not in band_set]
    require(len(residual) == 181_064, "post-THM-4254 count changed")
    require(edge_fnv(residual) == 0x8F550DCC2E552962,
            "post-THM-4254 FNV changed")
    require(edge_sha(residual) ==
            "0167652b41139bfd00c52236338fdd50e3be604641fe03e71eb66c68ee497d35",
            "post-THM-4254 SHA changed")
    return inherited, band, residual


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("repo", type=Path)
    parser.add_argument("--edges", action="store_true")
    args = parser.parse_args()
    inherited, band, residual = reconstruct(args.repo.resolve())
    if args.edges:
        for q, r in residual:
            print(f"{q},{r}")
        return

    print("POST_THM4254_RESIDUAL_INDEPENDENT_V1")
    print(
        f"INHERITED COUNT {len(inherited)} FNV {edge_fnv(inherited):016x} "
        f"SHA256 {edge_sha(inherited)}"
    )
    print(
        f"THM4254_SEMANTIC_BAND COUNT {len(band)} "
        f"FNV {edge_fnv(band):016x} SHA256 {edge_sha(band)}"
    )
    print(
        f"RESIDUAL COUNT {len(residual)} FNV {edge_fnv(residual):016x} "
        f"SHA256 {edge_sha(residual)}"
    )
    maximum = max(r for _, r in residual)
    top = [edge for edge in residual if edge[1] == maximum]
    print(
        f"TOP ENDPOINT {maximum} ROWS {len(top)} EDGES "
        + " ".join(f"{q},{r}" for q, r in top)
    )
    print("VERDICT PASS")


if __name__ == "__main__":
    main()
