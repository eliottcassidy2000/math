#!/usr/bin/env python3
"""Freeze the 2,530 small-ruler physical contexts used by the head-263 scan."""

from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
COMPILER = ROOT / "04-computation/lrc14_connected_low_all_heads_universal_forest_thm3352.py"
OUTPUT = ROOT / "04-computation/lrc14_disconnected_head263_contexts_20260812.txt"
EXPECTED_COMPILER = "aebfe98ab72f7eb0dc1718dfb17529e5b3f9288c6ed97d57f048bf3b29281f19"
EXPECTED_CONTEXTS = "efea6bd97522fe1c31a5a88ca9f3223f9e7a8c08e3be85c493e9f62fdfaf06e4"


def require(condition, description):
    if not condition:
        raise RuntimeError(description)


def file_hash(path):
    return sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def load(name, path):
    spec = spec_from_file_location(name, path)
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main():
    require(file_hash(COMPILER) == EXPECTED_COMPILER, ("compiler hash", file_hash(COMPILER)))
    compiler = load("head263_context_compiler", COMPILER)
    rows = tuple(row for row in compiler.feasible_contexts() if row[0] < 4592)
    require(len(rows) == 2530 and rows == tuple(sorted(rows)), ("context universe", len(rows)))
    for L, cell, e, f in rows:
        unit = L // 14
        require(L % 14 == 0, ("ruler", L))
        require(unit <= e * cell % L and e * cell % L + e <= L - unit,
                ("first endpoint", L, cell, e, f))
        require(unit <= f * cell % L and f * cell % L + f <= L - unit,
                ("second endpoint", L, cell, e, f))
    payload = "".join(f"{L} {cell} {e} {f}\n" for L, cell, e, f in rows)
    require(sha256(payload.encode("ascii")).hexdigest() == EXPECTED_CONTEXTS,
            "context semantic hash")
    OUTPUT.write_text(payload, encoding="ascii", newline="\n")
    print("LRC14 DISCONNECTED HEAD263 CONTEXT BANK")
    print("contexts", len(rows), "rulers", len({row[0] for row in rows}),
          "cells", len({row[:2] for row in rows}))
    print("compiler_sha256", file_hash(COMPILER))
    print("context_sha256", file_hash(OUTPUT))


if __name__ == "__main__":
    main()
