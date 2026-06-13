#!/usr/bin/env python3
"""
opus-2026-05-29-S10

Probe INV-191: H=63 unlocks at n=8 via complete Omega.

The S8 audit gave one concrete n=8 tournament with H(T)=63 and
Omega(T)=K31.  This script asks whether that mechanism appears stable in
random n=8 samples:

  * sample labeled n=8 tournaments deterministically;
  * collect all H=63 hits;
  * for each hit compute the full directed-odd-cycle conflict graph;
  * canonicalize hit tournaments up to isomorphism (brute force, n=8 only);
  * report whether any H=63 sample has non-complete Omega or alpha_2>0.

All output is intended to be captured with tee into 05-knowledge/results/.
"""

from __future__ import annotations

import itertools
import random
from collections import Counter, defaultdict


N = 8
TRIALS = 100_000
SEED = 2026_0529_10


def rows_from_mask(mask: int, n: int = N) -> tuple[int, ...]:
    rows = [0] * n
    bit = 0
    for i in range(n):
        for j in range(i + 1, n):
            if (mask >> bit) & 1:
                rows[i] |= 1 << j
            else:
                rows[j] |= 1 << i
            bit += 1
    return tuple(rows)


def hamiltonian_paths(rows: tuple[int, ...]) -> int:
    n = len(rows)
    full = (1 << n) - 1
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1 << n):
        missing = full ^ mask
        for v in range(n):
            val = dp[mask][v]
            if not val:
                continue
            nxt = rows[v] & missing
            while nxt:
                lsb = nxt & -nxt
                w = lsb.bit_length() - 1
                dp[mask | lsb][w] += val
                nxt ^= lsb
    return sum(dp[full])


def directed_odd_cycles(rows: tuple[int, ...]) -> list[tuple[int, ...]]:
    n = len(rows)
    cycles: list[tuple[int, ...]] = []
    for length in range(3, n + 1, 2):
        for perm in itertools.permutations(range(n), length):
            if perm[0] != min(perm):
                continue
            ok = True
            for i in range(length):
                if not (rows[perm[i]] >> perm[(i + 1) % length]) & 1:
                    ok = False
                    break
            if ok:
                cycles.append(perm)
    return cycles


def omega_profile(rows: tuple[int, ...]) -> dict[str, object]:
    cycles = directed_odd_cycles(rows)
    supports = [set(c) for c in cycles]
    nonedges = []
    for i, j in itertools.combinations(range(len(cycles)), 2):
        if supports[i].isdisjoint(supports[j]):
            nonedges.append((i, j))
    alpha = [1, len(cycles), len(nonedges)]
    if nonedges:
        # n=8 cannot fit three vertex-disjoint odd cycles, so alpha_3=0.
        ip2 = 1 + 2 * alpha[1] + 4 * alpha[2]
    else:
        ip2 = 1 + 2 * alpha[1]
    return {
        "num_cycles": len(cycles),
        "lengths": Counter(map(len, cycles)),
        "complete": not nonedges,
        "alpha": alpha if nonedges else [1, len(cycles)],
        "ip2": ip2,
        "nonedge_count": len(nonedges),
    }


def score_sequence(rows: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(sorted(row.bit_count() for row in rows))


def canonical(rows: tuple[int, ...]) -> tuple[int, ...]:
    n = len(rows)
    best = None
    for p in itertools.permutations(range(n)):
        inv = [0] * n
        for new, old in enumerate(p):
            inv[old] = new
        bits = []
        for i in range(n):
            row = 0
            old_i = p[i]
            for j in range(n):
                if i != j and ((rows[old_i] >> p[j]) & 1):
                    row |= 1 << j
            bits.append(row)
        tup = tuple(bits)
        if best is None or tup < best:
            best = tup
    assert best is not None
    return best


def automorphism_count(rows: tuple[int, ...]) -> int:
    n = len(rows)
    count = 0
    for p in itertools.permutations(range(n)):
        ok = True
        for i in range(n):
            if not ok:
                break
            for j in range(n):
                if i == j:
                    continue
                if ((rows[i] >> j) & 1) != ((rows[p[i]] >> p[j]) & 1):
                    ok = False
                    break
        if ok:
            count += 1
    return count


def main() -> None:
    rng = random.Random(SEED)
    hits: list[tuple[int, tuple[int, ...], dict[str, object]]] = []
    h_counter: Counter[int] = Counter()
    near_counter: Counter[int] = Counter()

    for trial in range(1, TRIALS + 1):
        mask = rng.getrandbits(N * (N - 1) // 2)
        rows = rows_from_mask(mask)
        h = hamiltonian_paths(rows)
        h_counter[h] += 1
        if 51 <= h <= 75:
            near_counter[h] += 1
        if h == 63:
            hits.append((mask, rows, omega_profile(rows)))

    by_profile: Counter[tuple[tuple[int, int], bool, tuple[int, ...], int]] = Counter()
    by_score: Counter[tuple[int, ...]] = Counter()
    by_canon: dict[tuple[int, ...], list[int]] = defaultdict(list)

    for mask, rows, prof in hits:
        length_key = tuple(sorted((int(k), int(v)) for k, v in prof["lengths"].items()))
        alpha_key = tuple(int(x) for x in prof["alpha"])
        by_profile[(length_key, bool(prof["complete"]), alpha_key, int(prof["ip2"]))] += 1
        by_score[score_sequence(rows)] += 1
        by_canon[canonical(rows)].append(mask)

    print("H=63 complete-Omega probe (opus-2026-05-29-S10)")
    print("=" * 72)
    print(f"n = {N}")
    print(f"trials = {TRIALS}")
    print(f"seed = {SEED}")
    print()
    print(f"H=63 hits = {len(hits)} ({len(hits) / TRIALS:.4%})")
    print(f"near-H distribution [51,75] = {dict(sorted(near_counter.items()))}")
    print()

    print("Omega profiles among H=63 hits:")
    if not by_profile:
        print("  no hits")
    for (length_key, complete, alpha_key, ip2), count in by_profile.most_common():
        print(
            f"  count={count:3d} complete={complete} "
            f"lengths={dict(length_key)} alpha={alpha_key} I2={ip2}"
        )

    print()
    print("Score sequences among H=63 hits:")
    for scores, count in by_score.most_common():
        print(f"  count={count:3d} scores={scores}")

    print()
    print(f"Distinct sampled isomorphism classes among H=63 hits = {len(by_canon)}")
    for idx, (canon, masks) in enumerate(by_canon.items(), 1):
        rows = rows_from_mask(masks[0])
        print(
            f"  class {idx}: labeled_hits={len(masks)}, "
            f"scores={score_sequence(rows)}, aut={automorphism_count(rows)}"
        )

    bad = [
        (mask, prof)
        for mask, _rows, prof in hits
        if not bool(prof["complete"]) or int(prof["ip2"]) != 63
    ]
    print()
    print(f"H=63 hits with non-complete Omega or I2!=63 = {len(bad)}")
    if bad:
        for mask, prof in bad[:5]:
            print(f"  mask={mask} profile={prof}")
    else:
        print("  none in this sample")
        print("  empirical lead: sampled H=63 at n=8 always arrived as Omega=K31.")


if __name__ == "__main__":
    main()
