#!/usr/bin/env python3
"""Exact forced-high ray closure for the last remaining k=2 row below 1780.

The global 1750..1799 scalar band leaves ten HIGH-TAIL rows.  The independent
z1=1790 verifier closes its three rows, the z1=1784 verifier closes its one
HIGH-TAIL row, the z1=1780 verifier closes its one, and the independent
z1=1750 verifier closes its four.  This file closes the only remaining row,
at 1768, by replacing the
analytic placeholder with the attained ray
maximum.  It asserts that disjoint partition against the all-body band output.
"""

from __future__ import annotations

import argparse
import multiprocessing as mp
from fractions import Fraction as F
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
BASE = (
    ROOT
    / "04-computation"
    / "lrc14_j7_k2_high_wall_descent_1800_1810_closure_thm2941.py"
)
BAND = ROOT / "05-knowledge/results/lrc14_j7_k2_scalar_band_1750_1799_thm2941.out"
Z1790 = ROOT / "05-knowledge/results/lrc14_j7_k2_z1790_exact_descent_closure_thm2941.out"
Z1784 = ROOT / "05-knowledge/results/lrc14_j7_k2_z1784_all_label_closure_thm2941.out"
Z1780 = ROOT / "05-knowledge/results/lrc14_j7_k2_z1780_hybrid_closure_thm2941.out"
OUTPUT_PATH = (
    ROOT
    / "05-knowledge/results/lrc14_j7_k2_high_wall_descent_1768_1790_closure_thm2941.out"
)
EXPECTED_BASE_SHA256 = (
    "22db2600485332469f6dec9e1356bc121165923a33c1ca6381acddcc89506f9e"
)
EXPECTED_BAND_SHA256 = (
    "2ce806d361d7eb97d9ae2d23e438898c8e1f895a89501c9a1847e51f61ca8009"
)
EXPECTED_Z1790_SHA256 = (
    "b03b46a6c438773ef1c433c435828a5426e18c998999cbee543541904b85f20b"
)
EXPECTED_Z1784_SHA256 = (
    "318f0908b2c84d68a5ef75884316f0ebc600d6aeed28d640f7489d6f179fb6bf"
)
EXPECTED_Z1780_SHA256 = (
    "afa96047c84a544b77142e34dfd25ad3923bace4d54f3f5ce47b0910f99ee27b"
)
CASES = (
    (1768, (1, 8, 10, 12, 13, 14)),
)
EXPECTED_GAPS = {
    CASES[0]: F(-906233, 6582732520),
}
EXPECTED_PROFILE_SHA256 = (
    "24dcf4d1e7d05aa2f69a8e191eb8d57ebed224486dd77c2ecbba034e1e15cd03"
)
EXPECTED_SEMANTIC_SHA256 = (
    "3d62a60ad5841b8190305fb5e60d7632e078466c20517a6930a870b91a9b67d2"
)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def file_sha256(path):
    return sha256(path.read_bytes()).hexdigest()


def ftext(value):
    return f"{value.numerator}/{value.denominator}"


require(file_sha256(BASE) == EXPECTED_BASE_SHA256, "high-wall engine changed")
require(file_sha256(BAND) == EXPECTED_BAND_SHA256, "1750..1799 atlas changed")
require(file_sha256(Z1790) == EXPECTED_Z1790_SHA256, "z1790 closure changed")
require(file_sha256(Z1784) == EXPECTED_Z1784_SHA256, "z1784 closure changed")
require(file_sha256(Z1780) == EXPECTED_Z1780_SHA256, "z1780 closure changed")
SPEC = spec_from_file_location("k2_high_wall_1768_base", BASE)
require(SPEC is not None and SPEC.loader is not None, "cannot load high-wall engine")
H = module_from_spec(SPEC)
SPEC.loader.exec_module(H)
H.EXPECTED_GAPS = EXPECTED_GAPS


def atlas_high_keys():
    keys = []
    for line in BAND.read_text().splitlines():
        if not line.startswith("SURVIVOR;") or "HIGH-TAIL" not in line:
            continue
        fields = dict(field.split("=", 1) for field in line.split(";")[1:] if "=" in field)
        keys.append((int(fields["z1"]), tuple(map(int, fields["E"].split(",")))))
    return tuple(sorted(keys))


require(
    tuple(key for key in atlas_high_keys() if 1751 <= key[0] <= 1779)
    == tuple(sorted(CASES)),
    "remaining 1751..1779 HIGH-TAIL rows changed",
)


def profile(case):
    return H.profile(case)


def render(profiles):
    require(tuple((row[0], row[1]) for row in profiles) == CASES, "case universe changed")
    require(all(row[-1] < 0 for row in profiles), "a high-wall scalar row survived")
    profile_hash = sha256(repr(tuple(profiles)).encode()).hexdigest()
    if EXPECTED_PROFILE_SHA256 is not None:
        require(profile_hash == EXPECTED_PROFILE_SHA256, "profile digest changed")
    semantic_payload = (
        CASES,
        tuple(row[-1] for row in profiles),
        profile_hash,
    )
    semantic_hash = sha256(repr(semantic_payload).encode()).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(semantic_hash == EXPECTED_SEMANTIC_SHA256, "semantic digest changed")
    lines = [
        "LRC14 projected k=2 forced-high exact-ray descent closure at 1768",
        f"high_wall_engine_sha256={file_sha256(BASE)}",
        f"scalar_band_sha256={file_sha256(BAND)}",
        f"z1790_closure_sha256={file_sha256(Z1790)}",
        f"z1784_closure_sha256={file_sha256(Z1784)}",
        f"z1780_closure_sha256={file_sha256(Z1780)}",
        "scope=last remaining scalar-empty HIGH-TAIL row below 1780;all later distinct nonaligned labels;no label horizon",
    ]
    for row in profiles:
        lines.append(
            f"CASE;z1={row[0]};E={','.join(map(str, row[1]))};L={row[4]};"
            f"high_floor={row[5]};ray_signs=+{row[8]}/-{row[9]}/0:{row[10]};"
            f"branch={row[15]};upper={ftext(row[-2])};gap={ftext(row[-1])};"
            "conclusion=SCALAR-EMPTY"
        )
    lines.extend(
        (
            "remaining_high_wall_rows=1;scalar_empty=1;survivors=0",
            f"profile_sha256={profile_hash}",
            f"semantic_sha256={semantic_hash}",
            "all_exact_controls=PASS",
        )
    )
    return "\n".join(lines) + "\n"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=min(len(CASES), mp.cpu_count() or 1))
    parser.add_argument("--output", type=Path, default=OUTPUT_PATH)
    args = parser.parse_args()
    require(args.workers >= 1, "workers must be positive")
    if args.workers == 1:
        profiles = [profile(case) for case in CASES]
    else:
        with mp.get_context("spawn").Pool(args.workers) as pool:
            profiles = list(pool.imap(profile, CASES))
    order = {case: index for index, case in enumerate(CASES)}
    profiles.sort(key=lambda row: order[(row[0], row[1])])
    output = render(tuple(profiles))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(output, encoding="utf-8", newline="\n")
    print(output, end="")


if __name__ == "__main__":
    main()
