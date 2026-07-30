#!/usr/bin/env python3
"""Hash-pinned splice proving the projected k=2 cap z1 <= 1742.

The upper descent proves z1 <= 1750.  Two independent all-body packages close
the complete integer slice 1750 and the complete contiguous band 1743..1749.
This verifier checks their exact artifacts and composes the intervals; it does
not substitute interpolation for any missing height.
"""

from __future__ import annotations

import argparse
from hashlib import sha256
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
UPPER_SOURCE = ROOT / "04-computation" / "lrc14_j7_k2_z1758_ray_status_closure_thm2941.py"
UPPER_OUTPUT = ROOT / "05-knowledge" / "results" / "lrc14_j7_k2_z1758_ray_status_closure_thm2941.out"
SLICE_SOURCE = ROOT / "04-computation" / "lrc14_j7_k2_z1750_literal_packet_projected_closure_thm2941.py"
SLICE_OUTPUT = ROOT / "05-knowledge" / "results" / "lrc14_j7_k2_z1750_literal_packet_projected_closure_thm2941.out"
BAND_SOURCE = ROOT / "04-computation" / "lrc14_j7_k2_band_1743_1749_literal_packet_closure_thm2941.py"
BAND_OUTPUT = ROOT / "05-knowledge" / "results" / "lrc14_j7_k2_band_1743_1749_literal_packet_closure_thm2941.out"
OUTPUT_PATH = ROOT / "05-knowledge" / "results" / "lrc14_j7_k2_cap1742_splice_thm2941.out"

EXPECTED_HASHES = (
    "369f0f24fa22e7bd7c6a5611b275b54663eb76dcb04f18d91897e093d02708da",
    "200a74e491008c11c6e012a48a2117e89eac8b803d0553387890d2496beb6dba",
    "70c9faa37a0524673e8178ed82cc6abc040438fff043b944a7ee0227d48c8997",
    "5288ace3fa398125045eeaea8c983c74d0019faa77d426784a612b7d2b4ae700",
    "635c7c0bf79cb524f42ac2b1002ab8b4e10fb67d417443bd21756ed907145865",
    "d046fe640140556b3964d876336ac26c233cc8dbc15a0e10bc5dc721b01f4666",
)
EXPECTED_PROFILE_SHA256 = "73dd88bbc686538240c7ca3aaa53adb075755484b2a120cd92a77f557f1e7d63"
EXPECTED_SEMANTIC_SHA256 = "acc67c59ce2cbf9d84ad4f97942a35dd79db4fc2e93b34614bfe2ba1549fae70"

PREVIOUS_LEDGER = 2_239_804
DECREMENTS = (1, 12, 2)
FINAL_LEDGER = 2_239_789
UPPER_CAP = 1750
CLOSED_INTERVAL = tuple(range(1743, 1751))
FINAL_CAP = 1742


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def file_sha256(path):
    return sha256(path.read_bytes()).hexdigest()


PATHS = (UPPER_SOURCE, UPPER_OUTPUT, SLICE_SOURCE, SLICE_OUTPUT, BAND_SOURCE, BAND_OUTPUT)
require(tuple(file_sha256(path) for path in PATHS) == EXPECTED_HASHES, "splice dependency changed")


def render():
    upper = UPPER_OUTPUT.read_text(encoding="utf-8")
    slice_text = SLICE_OUTPUT.read_text(encoding="utf-8")
    band = BAND_OUTPUT.read_text(encoding="utf-8")
    require("consequence=projected k=2 first drift label z1<=1750" in upper, "upper cap missing")
    require("global_z1758_rows=1;empty=1;survivors=0" in upper, "z1758 ledger missing")
    require("global_scalar_survivors=12" in slice_text, "z1750 survivor count changed")
    require("conclusion=the complete projected k=2 first-drift slice z1=1750 is empty" in slice_text, "z1750 completeness missing")
    require("global_scalar_survivors=2" in band, "lower-band survivor count changed")
    require("conclusion=the complete projected k=2 first-drift integer band 1743..1749 is empty" in band, "lower-band completeness missing")
    require(CLOSED_INTERVAL == (1743, 1744, 1745, 1746, 1747, 1748, 1749, 1750), "integer cover changed")
    require(UPPER_CAP == max(CLOSED_INTERVAL) and FINAL_CAP + 1 == min(CLOSED_INTERVAL), "splice is not contiguous")
    require(PREVIOUS_LEDGER - sum(DECREMENTS) == FINAL_LEDGER, "ledger arithmetic changed")
    profile_payload = (EXPECTED_HASHES, UPPER_CAP, CLOSED_INTERVAL, FINAL_CAP, PREVIOUS_LEDGER, DECREMENTS, FINAL_LEDGER)
    profile_hash = sha256(repr(profile_payload).encode()).hexdigest()
    if EXPECTED_PROFILE_SHA256 is not None:
        require(profile_hash == EXPECTED_PROFILE_SHA256, "profile digest changed")
    semantic_payload = ("projected k=2", "integer first-drift label", profile_payload, profile_hash)
    semantic_hash = sha256(repr(semantic_payload).encode()).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(semantic_hash == EXPECTED_SEMANTIC_SHA256, "semantic digest changed")
    lines = [
        "LRC14 projected k=2 cap-1742 exact splice",
        *(f"dependency_sha256={path.name}:{file_sha256(path)}" for path in PATHS),
        "upper_chain=all integer first-drift labels z1>=1751 are empty;cap z1<=1750",
        "lower_inputs=complete all-body slice z1=1750;complete all-body band z1=1743..1749",
        f"closed_integer_interval={CLOSED_INTERVAL}",
        "splice_check=no missing integer height;no interpolation",
        f"necessary_ledger={PREVIOUS_LEDGER}-1-12-2={FINAL_LEDGER}",
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
