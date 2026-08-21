#!/usr/bin/env python3
"""Exact no-trace Rule-A probes for the R=32768 offset boundary."""

from pathlib import Path
import sys
import time


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "04-computation"))

import amm12592_R16384_offset_transition_thm3644 as T


EXPECTED_ENGINE_SHA256 = (
    "8887080fc6e30760efa4a0ba76218ec97676cc717c6e76ccefbaeec6c73684ad"
)


def require(condition: bool, payload: object) -> None:
    if condition is not True:
        raise RuntimeError(payload)


def main() -> None:
    offsets = tuple(map(int, sys.argv[1:])) or (854, 855, 856)
    require(T.ENGINE_SHA256 == EXPECTED_ENGINE_SHA256,
            ("engine hash drift", T.ENGINE_SHA256))
    run_fast = T.load_engine(ROOT)
    print(f"R=32768;offsets={offsets};engine={T.ENGINE_SHA256}", flush=True)
    for offset in offsets:
        require(offset >= 0, ("negative offset", offset))
        started = time.perf_counter()
        row = T.stable_result(run_fast, 32768, offset)
        elapsed = time.perf_counter() - started
        print(f"offset={offset};stable_result={row};seconds={elapsed:.3f}",
              flush=True)


if __name__ == "__main__":
    main()
