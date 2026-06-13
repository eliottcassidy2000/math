#!/usr/bin/env python3
"""
tournament_clock_identify_s24.py

oracle-2026-06-01-S24

Follow-ups to tournament_clock_s24:
 (1) UNIVERSALITY: over many random speed sets, is the set of realizable
     phase-comparator tournament iso-classes (and their H-values) FIXED per n?
 (2) IDENTIFY: what are those tournaments?  (score sequences, regularity,
     rotational/circular structure)
 (3) transitive cell  <=>  all runners inside some open semicircle (a 1/2-gap).
"""

from __future__ import annotations

import random
from fractions import Fraction
from itertools import combinations, permutations
from functools import lru_cache
from collections import Counter

import importlib.util, pathlib
spec = importlib.util.spec_from_file_location(
    "clock", str(pathlib.Path(__file__).resolve().parent / "tournament_clock_s24.py"))
clock = importlib.util.module_from_spec(spec); spec.loader.exec_module(clock)


def is_regular(score):
    return len(set(score)) == 1


def is_rotational_like(adj):
    """Check if adj is iso to a circulant/rotational tournament:
    exists cyclic labelling with i->j depending only on (j-i) mod n."""
    n = len(adj)
    for perm in permutations(range(n)):
        ok = True
        # build relabeled adjacency
        rel = [[adj[perm[i]][perm[j]] for j in range(n)] for i in range(n)]
        pat = None
        for i in range(n):
            row = tuple(rel[i][(i + d) % n] for d in range(1, n))
            if pat is None:
                pat = row
            elif row != pat:
                ok = False
                break
        if ok:
            return True
    return False


def realizable_set(n, trials, seed):
    rng = random.Random(seed)
    classes = {}  # canon -> (H, score, example_adj)
    for _ in range(trials):
        speeds = tuple(sorted({0} | set(rng.sample(range(1, 30), n - 1))))
        if len(speeds) != n:
            continue
        for _, adj in clock.clock_cells(speeds):
            c = clock.canonical_form(adj)
            if c not in classes:
                classes[c] = (clock.hamiltonian_path_count(adj),
                              clock.score_sequence(adj), adj)
    return classes


def semicircle_check(speeds, samples, seed):
    """On random cell midpoints, verify: tournament transitive
    <=> all points lie in some open semicircle."""
    rng = random.Random(seed)
    cells = clock.clock_cells(speeds)
    mism = 0
    for t, adj in cells:
        score = clock.score_sequence(adj)
        transitive = (score == tuple(range(len(speeds))))
        pts = sorted(float(clock.frac(Fraction(s) * t)) for s in speeds)
        # all in a semicircle iff some gap between consecutive (circular) >= 1/2
        gaps = [pts[(i + 1) % len(pts)] - pts[i] for i in range(len(pts) - 1)]
        gaps.append(pts[0] + 1 - pts[-1])
        in_semi = max(gaps) > 0.5
        if transitive != in_semi:
            mism += 1
    return len(cells), mism


def main():
    print("Tournament clock — universality & identification (oracle-S24)\n")

    for n, trials in [(5, 150), (6, 60), (7, 12)]:
        print("=" * 64)
        print(f"n={n}: realizable phase-comparator tournaments over {trials} random speed sets")
        print("=" * 64, flush=True)
        classes = realizable_set(n, trials, seed=n * 100)
        total_iso = {3: 2, 4: 4, 5: 12, 6: 56, 7: 456}.get(n, "?")
        print(f" realizable iso-classes: {len(classes)}  (of {total_iso} total tournaments on {n})")
        Hs = sorted({h for h, _, _ in classes.values()})
        print(f" realizable H-values: {Hs}")
        print(" class details (H, score sequence, regular?, rotational?):")
        for c, (h, score, adj) in sorted(classes.items(), key=lambda kv: kv[1][0]):
            reg = is_regular(score)
            rot = is_rotational_like(adj)
            print(f"   H={h:<4} score={score} regular={reg} rotational={rot}")
        print(flush=True)

    print("=" * 64)
    print("transitive  <=>  all runners in an open semicircle (1/2-gap)")
    print("=" * 64)
    for label, s in {
        "arith 0..4": (0, 1, 2, 3, 4),
        "primes 0,2,3,5,7": (0, 2, 3, 5, 7),
        "spread 0,1,4,9,11": (0, 1, 4, 9, 11),
        "arith 0..5": (0, 1, 2, 3, 4, 5),
        "primes6 0,2,3,5,7,11": (0, 2, 3, 5, 7, 11),
    }.items():
        ncells, mism = semicircle_check(s, 0, seed=1)
        print(f"   [{label}] cells={ncells} mismatches={mism}  "
              f"({'CONFIRMED' if mism == 0 else 'FAILED'})")


if __name__ == "__main__":
    main()
