"""Bounded deterministic walk in the symplectic orbit of one F16 spread.

This is not an enumeration of all spreads, nor of its complete Sp(8,2) orbit.
Each recorded state is an actual spread obtained by invertible transvections.
"""
from collections import Counter
import hashlib
import importlib.util
import json
from pathlib import Path
import random

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("atlas", HERE / "105_17_probe.py")
atlas = importlib.util.module_from_spec(spec)
spec.loader.exec_module(atlas)


def residual_type(spread):
    pairs = []
    matching_codes = []
    for block in spread:
        duads = []
        for x in block:
            if x.bit_count() == 2:
                duads.append(tuple(i for i in range(10) if x >> i & 1))
            elif x.bit_count() == 8:
                duads.append(tuple(i for i in range(10) if not (x >> i & 1)))
        atlas.need(len(duads) in (1, 5), "length-ten code type")
        if len(duads) == 1:
            pairs.append(duads[0])
        else:
            matching_codes.append(duads)
    atlas.need(len(pairs) == 10 and len(matching_codes) == 7, "17=7+10")
    neighbours = [[] for _ in range(10)]
    for a, b in pairs:
        neighbours[a].append(b)
        neighbours[b].append(a)
    atlas.need(all(len(row) == 2 for row in neighbours), "residual degree two")
    seen = set()
    lengths = []
    for start in range(10):
        if start in seen:
            continue
        previous, current, length = -1, start, 0
        while current not in seen:
            seen.add(current)
            length += 1
            nxt = next(v for v in neighbours[current] if v != previous)
            previous, current = current, nxt
        atlas.need(current == start, "simple residual cycle")
        lengths.append(length)
    return tuple(sorted(lengths))


def main():
    seed_values = (2026102001, 2026102002, 2026102003, 2026102004)
    steps_per_walk = 2500
    vectors = [v for v in range(1, 512) if v.bit_count() % 2 == 0]
    initial, _, _ = atlas.make_spread()
    counts, first_witness = Counter(), {}
    walk_receipts = []
    for seed in seed_values:
        rng = random.Random(seed)
        spread = [list(block) for block in initial]
        word, local = [], Counter()
        for step in range(1, steps_per_walk + 1):
            v = rng.choice(vectors)
            word.append(v)
            spread = [[x ^ (v if (x & v).bit_count() % 2 else 0) for x in block]
                      for block in spread]
            cycle_type = residual_type(spread)
            counts[cycle_type] += 1
            local[cycle_type] += 1
            if cycle_type not in first_witness:
                # A second full checker verifies linearity, isotropy, coverage,
                # duad decomposition, and any even-cycle cubic colourings.
                witness = atlas.audit_spread(spread)
                atlas.need(tuple(witness["cycle_type"]) == cycle_type, "independent type")
                first_witness[cycle_type] = dict(seed=seed, step=step,
                                                transvection_word=word.copy(),
                                                full_witness=witness)
        walk_receipts.append(dict(seed=seed, steps=steps_per_walk,
                                  counts={str(k): v for k, v in sorted(local.items())}))
    result = {
        "universe": "four deterministic pseudorandom walks, 2500 steps each, "
                    "in the Sp(8,2) orbit of the displayed F16 spread; visits "
                    "are counted with multiplicity, not asserted distinct",
        "all_spreads_enumerated": False,
        "full_symplectic_orbit_enumerated": False,
        "total_visits": sum(counts.values()),
        "walks": walk_receipts,
        "cycle_type_counts": {str(k): v for k, v in sorted(counts.items())},
        "first_witnesses": {str(k): v for k, v in sorted(first_witness.items())},
        "source_sha256_LF": hashlib.sha256(
            Path(__file__).read_bytes().replace(b"\r\n", b"\n")).hexdigest(),
        "atlas_source_sha256_LF": hashlib.sha256(
            (HERE / "105_17_probe.py").read_bytes().replace(b"\r\n", b"\n")).hexdigest(),
    }
    Path(__file__).with_suffix(".json").write_text(
        json.dumps(result, indent=2) + "\n", encoding="utf-8", newline="\n")
    print("visits:", result["total_visits"])
    print("cycle types:", result["cycle_type_counts"])
    print("first witnesses:", {str(k): (v["seed"], v["step"])
                               for k, v in sorted(first_witness.items())})
    print("PASS bounded random walk, complete witnesses; no orbit exhaustiveness claim.")


if __name__ == "__main__":
    main()
