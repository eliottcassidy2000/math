#!/usr/bin/env python3
"""Hash-pinned scalar splice proving the projected k=2 cap z1 <= 1736.

The preceding exact chain proves z1 <= 1742.  A fresh all-3,003-body scalar
atlas exhausts every first label from 1680 through 1742.  Its largest
scalar-eligible height is 1736, so the complete integer interval 1737..1742
is empty by a necessary condition, without interpolation.  This verifier
pins both artifacts and checks that composition.
"""

from __future__ import annotations

import argparse
from collections import Counter
from hashlib import sha256
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CAP_SOURCE = ROOT / "04-computation/lrc14_j7_k2_cap1742_splice_thm2941.py"
CAP_OUTPUT = ROOT / "05-knowledge/results/lrc14_j7_k2_cap1742_splice_thm2941.out"
ATLAS_SOURCE = ROOT / "04-computation/lrc14_j7_k2_scalar_band_1680_1742_thm2941.py"
ATLAS_OUTPUT = ROOT / "05-knowledge/results/lrc14_j7_k2_scalar_band_1680_1742_thm2941.out"
OUTPUT_PATH = ROOT / "05-knowledge/results/lrc14_j7_k2_cap1736_scalar_splice_thm2941.out"

EXPECTED_HASHES = (
    "847f6619cb6aa0039c3335a77b545ec8ad3017820f1a6b3234fe0a0b3ebccb14",
    "61b8cd7ae66a212ea60423a77b7cc6d9b1de67f5af2469e8b13043de17c8dcd6",
    "1224d5594571f21c91c55fe3ab165c4fc34ba7968719862d12660d24efac919d",
    "c20607cb478ed491d7000f2b8a49213f57d1606a5152700ac3b50c836e2dc66c",
)
EXPECTED_PROFILE_SHA256 = (
    "5a17a8f7ed05f7a9cee48e94a980de08a23db99b1ebecc4f79a40d0ff1272a2a"
)
EXPECTED_SEMANTIC_SHA256 = (
    "5dd3a81f7f7ab362dd9c068a7595c067c4891f7d65884744152262b5ee7d47ec"
)

UPPER_CAP = 1742
ATLAS_INTERVAL = tuple(range(1680, 1743))
CLOSED_INTERVAL = tuple(range(1737, 1743))
FINAL_CAP = 1736
EXPECTED_HEIGHT_COUNTS = (
    (1683, 1),
    (1694, 10),
    (1702, 3),
    (1708, 14),
    (1722, 11),
    (1724, 2),
    (1732, 2),
    (1736, 15),
)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def file_sha256(path):
    return sha256(path.read_bytes()).hexdigest()


PATHS = (CAP_SOURCE, CAP_OUTPUT, ATLAS_SOURCE, ATLAS_OUTPUT)
require(
    tuple(file_sha256(path) for path in PATHS) == EXPECTED_HASHES,
    "cap/atlas dependency changed",
)


def atlas_heights(text):
    heights = []
    for line in text.splitlines():
        if not line.startswith("SURVIVOR;"):
            continue
        fields = dict(
            field.split("=", 1) for field in line.split(";")[1:] if "=" in field
        )
        heights.append(int(fields["z1"]))
    return tuple(sorted(Counter(heights).items()))


def render():
    cap = CAP_OUTPUT.read_text(encoding="utf-8")
    atlas = ATLAS_OUTPUT.read_text(encoding="utf-8")
    require(
        "consequence=projected k=2 first drift label z1<=1742" in cap,
        "inherited cap missing",
    )
    require("first_band=1680..1742" in atlas, "atlas universe changed")
    require(
        "first_candidate_rows=189105" in atlas,
        "all-body candidate universe changed",
    )
    height_counts = atlas_heights(atlas)
    require(height_counts == EXPECTED_HEIGHT_COUNTS, "atlas height profile changed")
    require(sum(count for _height, count in height_counts) == 58, "survivor count changed")
    require(max(height for height, _count in height_counts) == FINAL_CAP, "atlas top changed")
    require(
        CLOSED_INTERVAL == tuple(range(FINAL_CAP + 1, UPPER_CAP + 1)),
        "splice is not contiguous",
    )
    require(
        not any(height in CLOSED_INTERVAL for height, _count in height_counts),
        "a scalar row survived in the claimed empty splice",
    )
    require("all_exact_controls=PASS" in atlas, "atlas controls missing")
    profile_payload = (
        EXPECTED_HASHES,
        UPPER_CAP,
        ATLAS_INTERVAL,
        CLOSED_INTERVAL,
        EXPECTED_HEIGHT_COUNTS,
        FINAL_CAP,
    )
    profile_hash = sha256(repr(profile_payload).encode()).hexdigest()
    if EXPECTED_PROFILE_SHA256 is not None:
        require(profile_hash == EXPECTED_PROFILE_SHA256, "profile digest changed")
    semantic_payload = (
        "projected k=2",
        "scalar necessary first-drift atlas",
        profile_payload,
        profile_hash,
    )
    semantic_hash = sha256(repr(semantic_payload).encode()).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(semantic_hash == EXPECTED_SEMANTIC_SHA256, "semantic digest changed")
    lines = [
        "LRC14 projected k=2 cap-1736 scalar splice",
        *(f"dependency_sha256={path.name}:{file_sha256(path)}" for path in PATHS),
        "inherited_chain=all integer first-drift labels z1>=1743 are empty;cap z1<=1742",
        "atlas_universe=all 3003 bodies;every integer first label 1680..1742;scalar necessary condition",
        f"atlas_height_counts={height_counts}",
        f"closed_integer_interval={CLOSED_INTERVAL}",
        "splice_check=no scalar-eligible row at 1737..1742;no interpolation",
        f"consequence=projected k=2 first drift label z1<={FINAL_CAP}",
        f"profile_sha256={profile_hash}",
        f"semantic_sha256={semantic_hash}",
        "all_exact_controls=PASS",
    ]
    return "\n".join(lines) + "\n"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--hash-seed", type=int, default=0)
    parser.add_argument("--output", type=Path, default=OUTPUT_PATH)
    args = parser.parse_args()
    require(args.hash_seed >= 0, "hash seed must be nonnegative")
    output = render()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(output, encoding="utf-8", newline="\n")
    print(output, end="")


if __name__ == "__main__":
    main()
