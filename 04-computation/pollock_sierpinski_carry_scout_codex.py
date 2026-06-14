#!/usr/bin/env python3
"""Pollock tetrahedral scout through a Sierpinski/carry lens.

Context
-------
Pollock's tetrahedral conjecture says every positive integer is a sum of at
most five tetrahedral numbers Te_k = C(k+2, 3).  The previous Codex scout
rewrote the five-sum problem in terms of D_4, the integers not covered by at
most four tetrahedral numbers:

    n = Te_k + r is missed by both Te_k and Te_{k-1}
    exactly when r and r + tri(k) are both in D_4.

This script asks whether the Sierpinski/Waring dyadic intuition gives a local
obstruction, or whether it only survives as a carry-window/self-correlation
phenomenon.

Tournament Analysis
-------------------
Vertices: dyadic levels e=3..12.
Pairwise observable: tail pair-residue compression
    log2((2^e)^2 / observed_tail_pair_classes_e)
for defect-pairs with triangular gap index k >= 100.
Gauge: orient e -> f when e has larger compression.
Tie Hamiltonian path: increasing dyadic level.
Reported fingerprints: score histogram, directed 3-cycles, SCC sizes,
Hamiltonian-path count, and the champion path.

The point is intentionally asymmetric: single tetrahedral atoms appear
surjective modulo 2^e in this range, so single-residue obstructions are the
wrong vertex set.  Pair residues remain sparse once the triangular-gap
self-correlation condition is imposed.
"""

from __future__ import annotations

from bisect import bisect_left
from collections import Counter, defaultdict
from itertools import combinations
from math import log2


LIMIT = 1_000_000
DYADIC_MAX_E = 12
TAIL_PAIR_CUTOFF = 100
PAIR_CUTOFFS = (0, 50, 100, 300, 500, 800)


def tetra(k: int) -> int:
    return k * (k + 1) * (k + 2) // 6


def tri(k: int) -> int:
    return k * (k + 1) // 2


def v2(n: int) -> int:
    if n == 0:
        return 10**9
    return (n & -n).bit_length() - 1


def atoms_upto(limit: int) -> list[int]:
    atoms: list[int] = []
    k = 1
    while tetra(k) <= limit:
        atoms.append(tetra(k))
        k += 1
    return atoms


def reach_at_most(atoms: list[int], limit: int, terms: int) -> int:
    """Return bitset of integers <= limit reachable by at most `terms` atoms."""
    mask = (1 << (limit + 1)) - 1
    bits = 1
    for _ in range(terms):
        prev = bits
        nxt = prev
        for atom in atoms:
            nxt |= prev << atom
        bits = nxt & mask
    return bits


def bit_positions(bits: int, limit: int) -> list[int]:
    out: list[int] = []
    x = bits & ((1 << (limit + 1)) - 1)
    while x:
        low = x & -x
        pos = low.bit_length() - 1
        out.append(pos)
        x ^= low
    return out


def missing_positive(reach_bits: int, limit: int) -> list[int]:
    all_positive = ((1 << (limit + 1)) - 1) ^ 1
    miss_bits = (~reach_bits) & all_positive
    return bit_positions(miss_bits, limit)


def defect_pairs(defects: list[int], limit: int) -> list[tuple[int, int, int]]:
    defect_set = set(defects)
    pairs: list[tuple[int, int, int]] = []
    k = 1
    while tri(k) <= limit:
        gap = tri(k)
        for left in defects:
            right = left + gap
            if right > limit:
                break
            if right in defect_set:
                pairs.append((k, left, right))
        k += 1
    pairs.sort(key=lambda x: (x[0], x[1], x[2]))
    return pairs


def tetra_residue_stats(e: int, defects: list[int], pairs: list[tuple[int, int, int]]) -> dict[str, object]:
    m = 1 << e
    period_scan = 4 * m + 16
    atom_residues = {tetra(k) % m for k in range(period_scan)}
    defect_residues = {d % m for d in defects}
    pair_classes = {(left % m, right % m) for _, left, right in pairs}
    tail_pair_classes = {
        (left % m, right % m)
        for k, left, right in pairs
        if k >= TAIL_PAIR_CUTOFF
    }
    return {
        "e": e,
        "m": m,
        "period_scan": period_scan,
        "atom_residue_count": len(atom_residues),
        "atom_surjective": len(atom_residues) == m,
        "defect_residue_count": len(defect_residues),
        "pair_class_count": len(pair_classes),
        "tail_pair_class_count": len(tail_pair_classes),
        "possible_pair_classes": m * m,
        "compression": log2((m * m) / max(1, len(tail_pair_classes))),
    }


def adjacency_from_scores(levels: list[int], score: dict[int, float]) -> dict[int, set[int]]:
    adj = {v: set() for v in levels}
    for a, b in combinations(levels, 2):
        if score[a] > score[b]:
            adj[a].add(b)
        elif score[b] > score[a]:
            adj[b].add(a)
        else:
            # Fixed tie path: increasing dyadic level.
            adj[min(a, b)].add(max(a, b))
    return adj


def directed_3_cycles(levels: list[int], adj: dict[int, set[int]]) -> int:
    total = 0
    for a, b, c in combinations(levels, 3):
        if b in adj[a] and c in adj[b] and a in adj[c]:
            total += 1
        elif c in adj[a] and b in adj[c] and a in adj[b]:
            total += 1
    return total


def scc_sizes(levels: list[int], adj: dict[int, set[int]]) -> list[int]:
    radj = {v: set() for v in levels}
    for v in levels:
        for w in adj[v]:
            radj[w].add(v)

    seen: set[int] = set()
    order: list[int] = []

    def dfs(v: int) -> None:
        seen.add(v)
        for w in adj[v]:
            if w not in seen:
                dfs(w)
        order.append(v)

    for v in levels:
        if v not in seen:
            dfs(v)

    seen.clear()
    sizes: list[int] = []

    def rdfs(v: int) -> int:
        seen.add(v)
        size = 1
        for w in radj[v]:
            if w not in seen:
                size += rdfs(w)
        return size

    for v in reversed(order):
        if v not in seen:
            sizes.append(rdfs(v))
    return sorted(sizes, reverse=True)


def hamiltonian_path_count(levels: list[int], adj: dict[int, set[int]]) -> int:
    idx = {v: i for i, v in enumerate(levels)}
    n = len(levels)
    dp: dict[tuple[int, int], int] = {}
    for v in levels:
        dp[(1 << idx[v], v)] = 1
    for mask in range(1 << n):
        for last in levels:
            count = dp.get((mask, last), 0)
            if not count:
                continue
            for nxt in adj[last]:
                bit = 1 << idx[nxt]
                if mask & bit:
                    continue
                dp[(mask | bit, nxt)] = dp.get((mask | bit, nxt), 0) + count
    full = (1 << n) - 1
    return sum(dp.get((full, v), 0) for v in levels)


def greedy_champion_path(levels: list[int], adj: dict[int, set[int]]) -> list[int]:
    remaining = set(levels)
    path: list[int] = []
    while remaining:
        winner = max(remaining, key=lambda v: (len(adj[v] & remaining), v))
        path.append(winner)
        remaining.remove(winner)
    return path


def print_pair_histograms(pairs: list[tuple[int, int, int]]) -> None:
    print("\n=== Triangular-gap pair histograms inside D_4 ===")
    for cutoff in PAIR_CUTOFFS:
        tail = [(k, left, right) for k, left, right in pairs if k >= cutoff]
        v2_hist = Counter(v2(tri(k)) for k, _, _ in tail)
        kmod_hist = Counter(k % 16 for k, _, _ in tail)
        print(
            f"k >= {cutoff:3d}: pairs={len(tail):3d}, "
            f"max_k={max((k for k, _, _ in tail), default=None)}, "
            f"v2(tri(k))={v2_hist.most_common()}, "
            f"k mod 16 top={kmod_hist.most_common(8)}"
        )


def print_dyadic_stats(defects: list[int], pairs: list[tuple[int, int, int]]) -> list[dict[str, object]]:
    print("\n=== Dyadic residue audit ===")
    print("Single tetrahedral atom residues are scanned over k < 4*2^e+16.")
    stats: list[dict[str, object]] = []
    for e in range(1, DYADIC_MAX_E + 1):
        stat = tetra_residue_stats(e, defects, pairs)
        stats.append(stat)
        print(
            f"e={e:2d} m={stat['m']:5d} "
            f"atom_res={stat['atom_residue_count']:5d}/{stat['m']:<5d} "
            f"surj={str(stat['atom_surjective']):5s} "
            f"D4_res={stat['defect_residue_count']:4d} "
            f"pair_classes={stat['pair_class_count']:4d} "
            f"tail(k>={TAIL_PAIR_CUTOFF})={stat['tail_pair_class_count']:4d} "
            f"compression={stat['compression']:.3f}"
        )
    return stats


def print_tournament(stats: list[dict[str, object]]) -> None:
    print("\n=== Tournament Analysis: dyadic compression levels ===")
    levels = [int(s["e"]) for s in stats if int(s["e"]) >= 3]
    score = {int(s["e"]): float(s["compression"]) for s in stats if int(s["e"]) >= 3}
    adj = adjacency_from_scores(levels, score)
    outdegrees = {v: len(adj[v]) for v in levels}
    score_hist = Counter(outdegrees.values())
    cycles = directed_3_cycles(levels, adj)
    scc = scc_sizes(levels, adj)
    hpaths = hamiltonian_path_count(levels, adj)
    path = greedy_champion_path(levels, adj)
    edges = []
    for a, b in combinations(levels, 2):
        if b in adj[a]:
            edges.append(f"{a}>{b}")
        else:
            edges.append(f"{b}>{a}")
    print("vertices:", levels)
    print("observable: log2((2^e)^2 / observed tail pair-residue classes)")
    print("tie_path: increasing e")
    print("scores:", {e: round(score[e], 3) for e in levels})
    print("outdegrees:", outdegrees)
    print("score_hist:", dict(sorted(score_hist.items())))
    print("directed_3_cycles:", cycles)
    print("scc_sizes:", scc)
    print("hamiltonian_paths:", hpaths)
    print("champion_path:", path)
    print("edges:", ", ".join(edges))


def print_lucas_parity() -> None:
    print("\n=== Lucas/Sierpinski parity check ===")
    odd_residues_mod64 = [k for k in range(64) if tetra(k) % 2 == 1]
    verified = all((tetra(k) % 2 == 1) == (k % 4 == 1) for k in range(512))
    print("Te_k = C(k+2,3). Lucas gives Te_k odd iff k == 1 mod 4.")
    print("verified_up_to_k_511:", verified)
    print("odd k residues mod 64:", odd_residues_mod64)


def print_carry_windows(defects: list[int]) -> None:
    print("\n=== Near-next-tetrahedron carry windows for D_4 ===")
    max_d = max(defects)
    tvals = []
    k = 0
    while tetra(k) <= max_d + 20_000:
        tvals.append(tetra(k))
        k += 1
    rows = []
    for d in defects:
        h = bisect_left(tvals, d)
        upper = tvals[h]
        lower = tvals[h - 1] if h else 0
        shell_width = upper - lower
        gap_up = upper - d
        gap_down = d - lower
        rows.append((d, h, gap_up, gap_down, gap_up / shell_width))
    print("last 10 defects with upper-index, gap-up, gap-down, normalized gap-up:")
    for row in rows[-10:]:
        d, h, gap_up, gap_down, norm = row
        print(f"  d={d:6d} h={h:3d} up_gap={gap_up:5d} down_gap={gap_down:5d} norm={norm:.4f}")
    for bound in (1, 5, 10, 50, 100, 500, 1000, 5000, 10000):
        count = sum(1 for _, _, gap_up, _, _ in rows if gap_up <= bound)
        print(f"gap_up <= {bound:5d}: {count:3d}/{len(rows)}")
    tail_close = [(d, h, gap_up) for d, h, gap_up, _, _ in rows if d >= 30_000 and gap_up <= 100]
    print("tail defects within 100 below a tetrahedral number:", tail_close)


def print_defect_residue_samples(defects: list[int], pairs: list[tuple[int, int, int]]) -> None:
    print("\n=== Residue sample tables ===")
    for m in (16, 32, 64):
        defect_hist = Counter(d % m for d in defects)
        print(f"D4 residues mod {m}: classes={len(defect_hist)}, top={defect_hist.most_common(10)}")
    for e in (4, 5, 6):
        m = 1 << e
        tail = [(left % m, right % m) for k, left, right in pairs if k >= 300]
        hist = Counter(tail)
        print(f"tail pair classes k>=300 mod 2^{e}: count={len(hist)}, top={hist.most_common(10)}")


def main() -> None:
    atoms = atoms_upto(LIMIT)
    reach4 = reach_at_most(atoms, LIMIT, 4)
    reach5 = reach_at_most(atoms, LIMIT, 5)
    defects = missing_positive(reach4, LIMIT)
    five_misses = missing_positive(reach5, LIMIT)
    pairs = defect_pairs(defects, LIMIT)

    print("Pollock tetrahedral Sierpinski/carry scout")
    print(f"limit={LIMIT}")
    print(f"tetrahedral_atoms={len(atoms)}, largest_atom={atoms[-1]}")
    print(f"D4_defects={len(defects)}, max_defect={max(defects)}")
    print(f"five_term_misses={len(five_misses)}")
    print(f"triangular_defect_pairs={len(pairs)}")
    print(f"last_pair_by_k={max(pairs, key=lambda x: x[0])}")
    print(f"last_pair_by_right={max(pairs, key=lambda x: x[2])}")

    print_pair_histograms(pairs)
    stats = print_dyadic_stats(defects, pairs)
    print_tournament(stats)
    print_lucas_parity()
    print_carry_windows(defects)
    print_defect_residue_samples(defects, pairs)

    print("\n=== Proof-route interpretation ===")
    print("1. Pure dyadic local obstruction is unlikely: Te_k hits every class mod 2^e in the scan.")
    print("2. Sierpinski survives as Lucas parity of indices and as carry-window structure, not as a missing residue.")
    print("3. The live Pollock target is a no-long triangular self-correlation theorem for D_4.")
    print("4. Pair residues are sparse at high dyadic levels; use them as certificates/pruning, not as single-class obstructions.")


if __name__ == "__main__":
    main()
