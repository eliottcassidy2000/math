#!/usr/bin/env python3
"""Detached semantic proof-graph audit for reserved THM-4267.

Reads only the canonical frozen post-THM-4256 and post-THM-4266 edge lists,
selects the strict 4:5 ray by the defining integer relation, and recomputes
all removals and fingerprints without importing any proof-graph code.
"""

from __future__ import annotations

import argparse
import hashlib
from pathlib import Path


FNV_OFFSET = 0xCBF29CE484222325
FNV_PRIME = 0x100000001B3
MASK64 = (1 << 64) - 1


def read_edges(path: Path) -> list[tuple[int, int]]:
    edges = []
    for line in path.read_text(encoding="ascii").splitlines():
        q_text, r_text = line.split(",")
        edges.append((int(q_text), int(r_text)))
    assert edges == sorted(set(edges))
    return edges


def fnv_word(state: int, word: int) -> int:
    for shift in range(0, 64, 8):
        state ^= (word >> shift) & 0xFF
        state = state * FNV_PRIME & MASK64
    return state


def fnv(edges: list[tuple[int, int]]) -> str:
    state = FNV_OFFSET
    for q, r in edges:
        state = fnv_word(fnv_word(state, q), r)
    return f"{state:016x}"


def sha(edges: list[tuple[int, int]]) -> str:
    body = b"".join(f"{q},{r}\n".encode("ascii") for q, r in edges)
    return hashlib.sha256(body).hexdigest()


def fp(label: str, edges: list[tuple[int, int]]) -> None:
    print(f"{label} COUNT {len(edges)} FNV {fnv(edges)} SHA256 {sha(edges)}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo", type=Path, default=Path("."))
    args = parser.parse_args()
    root = (
        args.repo.resolve()
        / "05-knowledge/results/lrc14_three_round_learned_carrier_thm4266"
    )
    post4256 = read_edges(root / "post_thm4256_residual.csv")
    post4266 = read_edges(root / "thm4266_final_residual.csv")

    assert len(post4256) == 180_991
    assert fnv(post4256) == "021bf0ed1581657f"
    assert sha(post4256) == "9192c5d73aa5f123ddd10f0115dcaf7231fa518980610042e4cd3f8e73afd44f"
    assert len(post4266) == 177_585
    assert fnv(post4266) == "6ce05d05eb01daed"
    assert sha(post4266) == "009614651bb81e9763b2a9ff4b580497bfb6978a6c69d18cf986346e369374d9"

    ray = [(q, r) for q, r in post4256 if 5 * q == 4 * r and q > 290]
    scales = [q // 4 for q, _ in ray]
    assert all((q, r) == (4 * g, 5 * g) for (q, r), g in zip(ray, scales))
    overlap = [edge for edge in ray if edge not in set(post4266)]
    novel = [edge for edge in ray if edge in set(post4266)]
    final = [edge for edge in post4266 if edge not in set(novel)]

    assert scales == list(range(73, 132)) + [133, 134, 135, 137, 138, 139, 141]
    assert overlap == [(552, 690), (556, 695), (564, 705)]
    assert len(novel) == 63
    assert len(final) == 177_522
    assert max(r for _, r in final) == 688
    assert [edge for edge in final if edge[1] == 688] == [(520, 688)]

    print("THM4267_DETACHED_PROOF_GRAPH_REFEREE_V1")
    fp("POST_THM4256", post4256)
    fp("RAW_STRICT_4_5", ray)
    print("RAW_SCALES " + ",".join(map(str, scales)))
    fp("OVERLAP_THM4266", overlap)
    print("OVERLAP_SCALES " + ",".join(str(q // 4) for q, _ in overlap))
    fp("NOVEL_AFTER_THM4266", novel)
    fp("FINAL_AFTER_THM4267", final)
    print("FINAL_TOP UNIQUE 520,688")
    print("VERDICT PASS")


if __name__ == "__main__":
    main()
