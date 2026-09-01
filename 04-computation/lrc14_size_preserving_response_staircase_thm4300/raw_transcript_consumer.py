#!/usr/bin/env python3
"""Cross-check exchanged and augmented THM-4300 raw transcripts."""

from __future__ import annotations

import argparse
import hashlib
from pathlib import Path


def require(ok: bool, message: str) -> None:
    if not ok:
        raise RuntimeError(message)


def lines(path: Path) -> list[str]:
    return path.read_text(encoding="ascii").splitlines()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("exchanged", type=Path)
    parser.add_argument("augmented", type=Path)
    parser.add_argument("exchanged_failures", type=Path)
    parser.add_argument("augmented_failures", type=Path)
    args = parser.parse_args()

    exchanged = lines(args.exchanged)
    augmented = lines(args.augmented)
    require(exchanged[0] ==
            "ENDPOINT597_FULL_EXCHANGED_CARRIER_RAW_VERIFY_V1",
            "exchanged transcript header changed")
    require(exchanged[-1] ==
            "VERDICT PASS EXACT_9019_MASK_EXCHANGED_CARRIER_CLOSES_PREFIX",
            "exchanged transcript verdict changed")
    require("COMPLETED_ROWS 354 TOTAL_BODY_TESTS 5064731100 "
            "TOTAL_EXPOSED 1045532 TOTAL_HIT_INCIDENCES 44773293 "
            "TOTAL_FAILURES 0 LEDGER_FNV f9fe3fc694b6aa98" in exchanged,
            "exchanged raw summary changed")

    exchanged_detail = [line for line in exchanged
                         if line.startswith(("PAIR ", "LAYER "))]
    augmented_detail: list[str] = []
    for line in augmented:
        if line.startswith(("PAIR ", "LAYER ")):
            if line.startswith("PAIR "):
                endpoint = int(line.split(",", 1)[1].split()[0])
                if endpoint < 597:
                    break
            augmented_detail.append(line)
    require(len(exchanged_detail) == 394,
            "exchanged pair/layer line count changed")
    require(exchanged_detail == augmented_detail,
            "exchanged and augmented prefix details differ")
    require(any("PAIR 210,596" in line and "FAILURES 24" in line
                for line in augmented),
            "augmented endpoint-596 boundary changed")

    exchanged_failures = args.exchanged_failures.read_bytes()
    augmented_failures = lines(args.augmented_failures)
    require(exchanged_failures == b"q,r,body_hex\n",
            "exchanged failure ledger is not header-only")
    require(len(augmented_failures) == 25 and
            augmented_failures[0] == "q,r,body_hex" and
            all(line.startswith("210,596,") for line in
                augmented_failures[1:]),
            "augmented endpoint-596 failure ledger changed")

    detail_bytes = ("\n".join(exchanged_detail) + "\n").encode("ascii")
    print("THM4300_RAW_TRANSCRIPT_CONSUMER_V1")
    print(f"EXCHANGED_DETAIL lines={len(exchanged_detail)} "
          f"sha256={hashlib.sha256(detail_bytes).hexdigest()}")
    print("AUGMENTED_PREFIX_MATCH lines=394 diff=0")
    print("BOUNDARY endpoint=596 pair=(210,596) failures=24")
    print("VERDICT PASS RAW_TRANSCRIPTS_AGREE")


if __name__ == "__main__":
    main()
