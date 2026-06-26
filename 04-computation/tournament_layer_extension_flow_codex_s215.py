#!/usr/bin/env python3
"""Tournament layer-extension flow scout.

This script tests the user's diagonal-layer model in two finite ways.

1. Consecutive tiling layers of sizes k and k+1 joined by all k(k+1)
   cross-lines are not generic binary sheets.  If the line value is a pair
   potential x_i XOR y_j, every rectangle has zero parity.  Thus the whole
   K_{k,k+1} sheet is determined by k+(k+1)-1 boundary bits.

2. Growing a tournament from n to n+1 by adding one vertex is a rooted
   extension problem.  For a parent class P, the incident word of the new
   vertex is a binary word on V(P), modulo Aut(P).  Summing these word orbits
   over parent classes gives exactly the number of rooted/perspective
   tournament classes at n+1.

The second point is the clean version of the n=3 -> n=4 -> n=5 coincidence:
rooted perspectives R(3)=A000568(4)=4 and R(4)=A000568(5)=12, but the ladder
breaks at R(5)=48 while A000568(6)=56.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from itertools import permutations
from math import factorial


def pair_index(n: int) -> dict[tuple[int, int], int]:
    return {(i, j): idx for idx, (i, j) in enumerate((i, j) for i in range(n) for j in range(i + 1, n))}


def edge(mask: int, n: int, i: int, j: int, idx: dict[tuple[int, int], int]) -> int:
    """Return 1 iff i beats j."""
    if i == j:
        return 0
    if i < j:
        return (mask >> idx[(i, j)]) & 1
    return 1 - ((mask >> idx[(j, i)]) & 1)


def permute_mask(mask: int, n: int, perm: tuple[int, ...], idx: dict[tuple[int, int], int]) -> int:
    out = 0
    for a in range(n):
        for b in range(a + 1, n):
            # In the permuted labelled tournament, a beats b iff old perm[a]
            # beats old perm[b].
            if edge(mask, n, perm[a], perm[b], idx):
                out |= 1 << idx[(a, b)]
    return out


def canonical(mask: int, n: int, perms: list[tuple[int, ...]], idx: dict[tuple[int, int], int]) -> int:
    return min(permute_mask(mask, n, p, idx) for p in perms)


def automorphisms(mask: int, n: int, perms: list[tuple[int, ...]], idx: dict[tuple[int, int], int]) -> list[tuple[int, ...]]:
    return [p for p in perms if permute_mask(mask, n, p, idx) == mask]


def vertex_orbits(auts: list[tuple[int, ...]], n: int) -> list[tuple[int, ...]]:
    parent = list(range(n))

    def find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a: int, b: int) -> None:
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra

    for p in auts:
        for i, image in enumerate(p):
            union(i, image)
    buckets: dict[int, list[int]] = defaultdict(list)
    for i in range(n):
        buckets[find(i)].append(i)
    return [tuple(v) for v in buckets.values()]


def word_orbits(auts: list[tuple[int, ...]], n: int) -> list[list[int]]:
    """Aut(P)-orbits on incident words, where bit i belongs to old vertex i."""
    seen: set[int] = set()
    orbits: list[list[int]] = []
    for w in range(1 << n):
        if w in seen:
            continue
        orb = set()
        for p in auts:
            moved = 0
            # New word at position i reads old bit p[i].
            for i in range(n):
                if (w >> p[i]) & 1:
                    moved |= 1 << i
            orb.add(moved)
        seen.update(orb)
        orbits.append(sorted(orb))
    return orbits


def burnside_word_orbit_count(auts: list[tuple[int, ...]], n: int) -> int:
    total = 0
    for p in auts:
        seen = [False] * n
        cycles = 0
        for i in range(n):
            if seen[i]:
                continue
            cycles += 1
            j = i
            while not seen[j]:
                seen[j] = True
                j = p[j]
        total += 2**cycles
    return total // len(auts)


def extend_mask(parent: int, n: int, word: int, idx_n: dict[tuple[int, int], int], idx_child: dict[tuple[int, int], int]) -> int:
    """Add vertex n.  Word bit i=1 means new vertex n beats old i."""
    child = 0
    for i in range(n):
        for j in range(i + 1, n):
            if (parent >> idx_n[(i, j)]) & 1:
                child |= 1 << idx_child[(i, j)]
    for i in range(n):
        if (word >> i) & 1:
            # old i loses to new n, so edge i<n is oriented n -> i;
            # bit for (i,n) is 0 under our convention i beats n iff bit=1.
            pass
        else:
            child |= 1 << idx_child[(i, n)]
    return child


def canonical_cached(
    mask: int,
    n: int,
    perms: dict[int, list[tuple[int, ...]]],
    indices: dict[int, dict[tuple[int, int], int]],
    cache: dict[tuple[int, int], int],
) -> int:
    key = (n, mask)
    if key not in cache:
        cache[key] = canonical(mask, n, perms[n], indices[n])
    return cache[key]


def enumerate_classes(max_n: int = 7) -> tuple[
    dict[int, list[int]],
    dict[int, list[tuple[int, ...]]],
    dict[int, dict[tuple[int, int], int]],
    dict[tuple[int, int], int],
]:
    perms = {n: list(permutations(range(n))) for n in range(1, max_n + 1)}
    indices = {n: pair_index(n) for n in range(1, max_n + 1)}
    canon_cache: dict[tuple[int, int], int] = {}
    classes: dict[int, list[int]] = {1: [0]}
    for n in range(1, max_n):
        child_canons = set()
        for parent in classes[n]:
            idx_n = indices[n]
            idx_child = indices[n + 1]
            for w in range(1 << n):
                child = extend_mask(parent, n, w, idx_n, idx_child)
                child_canons.add(canonical_cached(child, n + 1, perms, indices, canon_cache))
        classes[n + 1] = sorted(child_canons)
    return classes, perms, indices, canon_cache


def sheet_rank_table(max_k: int = 8) -> list[dict[str, int]]:
    rows = []
    for k in range(1, max_k + 1):
        generic_bits = k * (k + 1)
        potential_bits = 2 * k  # k+(k+1)-1
        rows.append(
            {
                "k": k,
                "generic_bits": generic_bits,
                "potential_bits": potential_bits,
                "saved_bits": generic_bits - potential_bits,
                "admissible_sheets": 2**potential_bits,
            }
        )
    return rows


def rectangle_law_holds(k: int, x: int, y: int) -> bool:
    vals = {}
    for i in range(k):
        for j in range(k + 1):
            vals[(i, j)] = ((x >> i) & 1) ^ ((y >> j) & 1)
    for i in range(k):
        for ip in range(i + 1, k):
            for j in range(k + 1):
                for jp in range(j + 1, k + 1):
                    if vals[(i, j)] ^ vals[(ip, j)] ^ vals[(i, jp)] ^ vals[(ip, jp)]:
                        return False
    return True


def main() -> None:
    print("Tournament layer-extension flow scout (codex-2026-06-26-S215)")
    print()
    print("LAYER SHEET COMPRESSION")
    print("Line sheet between layer sizes k and k+1 uses e_ij=x_i XOR y_j.")
    print("Every 2x2 rectangle has parity zero; equivalently the sheet has boundary-potential rank.")
    print(" k generic_bits potential_bits saved_bits admissible_sheets")
    for row in sheet_rank_table():
        print(
            f"{row['k']:2d} {row['generic_bits']:12d} {row['potential_bits']:14d} "
            f"{row['saved_bits']:10d} {row['admissible_sheets']:17d}"
        )
    for k in range(1, 7):
        for x in range(1 << k):
            for y in range(1 << (k + 1)):
                assert rectangle_law_holds(k, x, y)
    print("Verified rectangle parity for all boundary words through k=6.")
    print()

    classes, perms, indices, canon_cache = enumerate_classes(7)
    print("UNROOTED AND ROOTED/PERSPECTIVE COUNTS")
    print("A(n)=unrooted classes; R(n)=sum vertex-orbits over classes.")
    print("E(n->n+1)=sum Aut(parent)-word-orbits; theorem check E(n->n+1)=R(n+1).")
    print(" n A(n) R(n) A(n+1)-R(n) E(n->n+1) child_merges")
    rooted: dict[int, int] = {}
    aut_cache: dict[tuple[int, int], list[tuple[int, ...]]] = {}
    for n in range(1, 8):
        r = 0
        orbit_hist = Counter()
        aut_hist = Counter()
        for c in classes[n]:
            auts = automorphisms(c, n, perms[n], indices[n])
            aut_cache[(n, c)] = auts
            vorb = vertex_orbits(auts, n)
            r += len(vorb)
            orbit_hist[len(vorb)] += 1
            aut_hist[len(auts)] += 1
        rooted[n] = r

    for n in range(1, 8):
        orbit_hist = Counter()
        aut_hist = Counter()
        for c in classes[n]:
            auts = aut_cache[(n, c)]
            orbit_hist[len(vertex_orbits(auts, n))] += 1
            aut_hist[len(auts)] += 1
        next_gap = classes.get(n + 1)
        gap_text = str(len(next_gap) - rooted[n]) if next_gap is not None else "-"
        if n < 7:
            e_count = 0
            child_targets = set()
            rooted_states = 0
            for parent in classes[n]:
                auts = aut_cache[(n, parent)]
                orbs = word_orbits(auts, n)
                b_count = burnside_word_orbit_count(auts, n)
                assert b_count == len(orbs)
                e_count += len(orbs)
                rooted_states += len(orbs)
                for orb in orbs:
                    child = extend_mask(parent, n, orb[0], indices[n], indices[n + 1])
                    child_targets.add(canonical_cached(child, n + 1, perms, indices, canon_cache))
            assert e_count == rooted[n + 1]
            child_merges = rooted_states - len(child_targets)
            e_text = str(e_count)
            merge_text = str(child_merges)
        else:
            e_text = "-"
            merge_text = "-"
        print(f"{n:2d} {len(classes[n]):4d} {rooted[n]:4d} {gap_text:>11s} {e_text:11s} {merge_text:12s}")
        print(f"    vertex_orbit_hist={dict(sorted(orbit_hist.items()))} aut_size_hist={dict(sorted(aut_hist.items()))}")
    print()
    print("PERSPECTIVE LADDER READOUT")
    for n in range(1, 7):
        marker = "==" if rooted[n] == len(classes[n + 1]) else "!="
        print(f"  R({n})={rooted[n]} {marker} A({n+1})={len(classes[n+1])}")
    print("The equality R(n)=A(n+1) holds through n=4 and fails at n=5.")
    print()

    print("DUPLICATION LAW")
    print("The line sheet has k(k+1) apparent binary values but only 2k degrees of freedom.")
    print("The one-point extension layer has 2^n apparent incident words per labelled parent,")
    print("but only Burnside word-orbits per parent class.  Both are the same doctrine:")
    print("keep the boundary/address carrier, then quotient by the stabilizer at the last step.")
    print()
    print("TOURNAMENT ANALYSIS")
    print("  vertices=layer carriers: boundary words, rank-one sheets, parent classes, rooted children, unrooted children")
    print("  observable=amount of route/perspective information retained before quotient")
    print("  switch/gauge=more retained information points to less retained information")
    print("  tie Hamiltonian path=boundary_words > rank_one_sheet > parent_word_orbit > rooted_child > unrooted_child")
    print("  challenged assumption=vertices are not original tournament vertices; they are extension carriers/proof states")
    print("  preserves=one-vertex extension and fixed-path layer address")
    print("  destroys=which child root/deletion address produced an unrooted class")


if __name__ == "__main__":
    main()
