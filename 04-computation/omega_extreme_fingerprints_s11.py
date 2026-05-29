#!/usr/bin/env python3
"""
opus-2026-05-29-S11

Fingerprint two opposite extremes of the odd-cycle conflict graph:

1. complete Omega classes, especially the two n=8 H=63 unlocks from THM-344;
2. the n=9 THM-025 real-rootedness counterexample, where Omega is almost
   complete but has one independent triple and non-real roots.

The goal is not another OCF verification.  It is a small structural audit of the
"where have we not looked?" axis: the disjointness graph of directed odd cycles.
"""

from __future__ import annotations

import itertools
import subprocess
from collections import Counter, defaultdict, deque
from math import factorial


def pairs(n: int) -> list[tuple[int, int]]:
    return [(i, j) for i in range(n) for j in range(i + 1, n)]


def gentourng_reps(n: int) -> list[int]:
    m = n * (n - 1) // 2
    result = subprocess.run(["gentourng", str(n)], capture_output=True, text=True, check=False)
    text = result.stdout if result.stdout else result.stderr
    reps = []
    for line in text.splitlines():
        line = line.strip()
        if len(line) == m and all(c in "01" for c in line):
            reps.append(int(line, 2))
    return reps


def rows_from_nauty_bits(n: int, bits: int) -> tuple[int, ...]:
    m = n * (n - 1) // 2
    rows = [0] * n
    for k, (i, j) in enumerate(pairs(n)):
        if bits & (1 << (m - 1 - k)):
            rows[i] |= 1 << j
        else:
            rows[j] |= 1 << i
    return tuple(rows)


def thm025_rows() -> tuple[int, ...]:
    out = {
        0: [1, 3, 6, 7],
        1: [3],
        2: [0, 1, 4, 5, 6, 7],
        3: [2, 5, 7],
        4: [0, 1, 3, 7],
        5: [0, 1, 4, 6, 7, 8],
        6: [1, 3, 4, 7],
        7: [1],
        8: [0, 1, 2, 3, 4, 6, 7],
    }
    rows = []
    for v in range(9):
        mask = 0
        for w in out[v]:
            mask |= 1 << w
        rows.append(mask)
    return tuple(rows)


def canonical_cycle_sequences(n: int) -> list[tuple[tuple[int, ...], int]]:
    seqs = []
    for length in range(3, n + 1, 2):
        for subset in itertools.combinations(range(n), length):
            root = subset[0]
            for tail in itertools.permutations(subset[1:]):
                seq = (root,) + tail
                mask = sum(1 << v for v in seq)
                seqs.append((seq, mask))
    return seqs


SEQ_CACHE: dict[int, list[tuple[tuple[int, ...], int]]] = {}


def odd_cycle_masks(rows: tuple[int, ...]) -> list[tuple[tuple[int, ...], int]]:
    n = len(rows)
    if n not in SEQ_CACHE:
        SEQ_CACHE[n] = canonical_cycle_sequences(n)
    cycles = []
    for seq, mask in SEQ_CACHE[n]:
        ok = True
        for i, v in enumerate(seq):
            if not ((rows[v] >> seq[(i + 1) % len(seq)]) & 1):
                ok = False
                break
        if ok:
            cycles.append((seq, mask))
    return cycles


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


def alpha_vector(cycle_masks: list[int]) -> list[int]:
    alpha = [1, len(cycle_masks)]
    disjoint_pairs = []
    for i, j in itertools.combinations(range(len(cycle_masks)), 2):
        if cycle_masks[i] & cycle_masks[j] == 0:
            disjoint_pairs.append((i, j))
    if not disjoint_pairs:
        return alpha
    alpha.append(len(disjoint_pairs))
    triples = 0
    for i, j, k in itertools.combinations(range(len(cycle_masks)), 3):
        if (
            cycle_masks[i] & cycle_masks[j] == 0
            and cycle_masks[i] & cycle_masks[k] == 0
            and cycle_masks[j] & cycle_masks[k] == 0
        ):
            triples += 1
    if triples:
        alpha.append(triples)
    quads = 0
    if len(cycle_masks) <= 160:
        for combo in itertools.combinations(range(len(cycle_masks)), 4):
            ok = True
            for a, b in itertools.combinations(combo, 2):
                if cycle_masks[a] & cycle_masks[b]:
                    ok = False
                    break
            if ok:
                quads += 1
    if quads:
        alpha.append(quads)
    return alpha


def ip_at_2(alpha: list[int]) -> int:
    return sum((2**k) * a for k, a in enumerate(alpha))


def score_sequence(rows: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(sorted(row.bit_count() for row in rows))


def scc_sizes(rows: tuple[int, ...]) -> tuple[int, ...]:
    n = len(rows)
    rev = [0] * n
    for v, row in enumerate(rows):
        for w in range(n):
            if (row >> w) & 1:
                rev[w] |= 1 << v

    def reach(start: int, graph: tuple[int, ...] | list[int]) -> int:
        seen = 0
        stack = [start]
        while stack:
            v = stack.pop()
            if (seen >> v) & 1:
                continue
            seen |= 1 << v
            nxt = graph[v] & ~seen
            while nxt:
                lsb = nxt & -nxt
                stack.append(lsb.bit_length() - 1)
                nxt ^= lsb
        return seen

    remaining = (1 << n) - 1
    sizes = []
    while remaining:
        v = (remaining & -remaining).bit_length() - 1
        comp = reach(v, rows) & reach(v, rev)
        sizes.append(comp.bit_count())
        remaining &= ~comp
    return tuple(sorted(sizes, reverse=True))


def delete_vertex(rows: tuple[int, ...], remove: int) -> tuple[int, ...]:
    old_to_new = {}
    new = 0
    for v in range(len(rows)):
        if v != remove:
            old_to_new[v] = new
            new += 1
    out = [0] * (len(rows) - 1)
    for v in range(len(rows)):
        if v == remove:
            continue
        nv = old_to_new[v]
        for w in range(len(rows)):
            if w != remove and ((rows[v] >> w) & 1):
                out[nv] |= 1 << old_to_new[w]
    return tuple(out)


def vertex_cycle_counts(n: int, masks: list[int]) -> tuple[int, ...]:
    return tuple(sum(1 for mask in masks if (mask >> v) & 1) for v in range(n))


def transitive_order(rows: tuple[int, ...]) -> tuple[int, ...] | None:
    """Return the unique source-to-sink order if rows form a transitive tournament."""
    n = len(rows)
    order = tuple(sorted(range(n), key=lambda v: rows[v].bit_count(), reverse=True))
    if tuple(rows[v].bit_count() for v in order) != tuple(range(n - 1, -1, -1)):
        return None
    for i, v in enumerate(order):
        for w in order[i + 1 :]:
            if not ((rows[v] >> w) & 1):
                return None
    return order


def weighted_signature_count(bits: str) -> int:
    """Odd-cycle count when a core vertex is inserted into a transitive order.

    bits[i] = 1 means core -> order[i], bits[i] = 0 means order[i] -> core.
    For each 1...0 pair, choose an even number of interior vertices.
    """
    total = 0
    for i, left in enumerate(bits):
        if left != "1":
            continue
        for j in range(i + 1, len(bits)):
            if bits[j] != "0":
                continue
            interior = j - i - 1
            total += 1 if interior == 0 else 1 << (interior - 1)
    return total


def single_core_signature(rows: tuple[int, ...], core: int) -> tuple[str, int] | None:
    sub = delete_vertex(rows, core)
    order = transitive_order(sub)
    if order is None:
        return None
    old_vertices = [v for v in range(len(rows)) if v != core]
    lifted_order = [old_vertices[i] for i in order]
    bits = "".join("1" if ((rows[core] >> v) & 1) else "0" for v in lifted_order)
    return bits, weighted_signature_count(bits)


def signature_spectrum(max_m: int = 8) -> str:
    lines = ["Single-core weighted-signature spectra", "=" * 42]
    for m in range(2, max_m + 1):
        by_value: dict[int, list[str]] = defaultdict(list)
        for mask in range(1 << m):
            bits = format(mask, f"0{m}b")
            by_value[weighted_signature_count(bits)].append(bits)
        support = sorted(by_value)
        lines.append(f"m={m} transitive vertices: support={support}")
        lines.append(f"  missing in range={ [r for r in range(support[0], support[-1] + 1) if r not in by_value] }")
        for target in (3, 10, 31):
            if target in by_value:
                examples = ", ".join(by_value[target][:6])
                lines.append(f"  r={target}: {len(by_value[target])} signatures; examples {examples}")
            else:
                lines.append(f"  r={target}: absent")
    return "\n".join(lines)


def support_tuple(mask: int) -> tuple[int, ...]:
    return tuple(v for v in range(mask.bit_length()) if (mask >> v) & 1)


def fingerprint(label: str, rows: tuple[int, ...]) -> str:
    n = len(rows)
    cycles = odd_cycle_masks(rows)
    masks = [mask for _, mask in cycles]
    alpha = alpha_vector(masks)
    core = (1 << n) - 1
    for mask in masks:
        core &= mask
    pair_intersections = Counter()
    disjoint_pairs = []
    for i, j in itertools.combinations(range(len(masks)), 2):
        inter = (masks[i] & masks[j]).bit_count()
        pair_intersections[inter] += 1
        if inter == 0:
            disjoint_pairs.append((i, j))
    independent_triples = []
    for i, j, k in itertools.combinations(range(len(masks)), 3):
        if masks[i] & masks[j] == masks[i] & masks[k] == masks[j] & masks[k] == 0:
            independent_triples.append((i, j, k))

    lines = []
    lines.append(f"{label}")
    lines.append("-" * len(label))
    lines.append(f"n={n}, H={hamiltonian_paths(rows)}, scores={score_sequence(rows)}, SCC sizes={scc_sizes(rows)}")
    lines.append(f"odd cycles={len(cycles)}, lengths={dict(sorted(Counter(len(seq) for seq, _ in cycles).items()))}")
    lines.append(f"alpha={alpha}, I(Omega,2)={ip_at_2(alpha)}, complete_Omega={len(disjoint_pairs) == 0}")
    lines.append(f"disjoint-pair count={len(disjoint_pairs)}, independent-triples={len(independent_triples)}")
    lines.append(f"cycle-family core={support_tuple(core) if core else ()}")
    if core and core.bit_count() == 1:
        core_vertex = (core & -core).bit_length() - 1
        sig = single_core_signature(rows, core_vertex)
        if sig is not None:
            lines.append(f"single-core signature at vertex {core_vertex}: bits={sig[0]}, weighted_count={sig[1]}")
    lines.append(f"vertex cycle participation={vertex_cycle_counts(n, masks)}")
    lines.append(f"pair intersection sizes={dict(sorted(pair_intersections.items()))}")
    if independent_triples:
        first = independent_triples[0]
        triple_supports = [support_tuple(masks[i]) for i in first]
        lines.append(f"first independent triple supports={triple_supports}")
    if n == 8:
        deletion_rows = []
        for v in range(n):
            sub = delete_vertex(rows, v)
            sub_cycles = odd_cycle_masks(sub)
            sub_masks = [mask for _, mask in sub_cycles]
            sub_alpha = alpha_vector(sub_masks)
            deletion_rows.append(
                (
                    v,
                    rows[v].bit_count(),
                    hamiltonian_paths(sub),
                    len(sub_cycles),
                    dict(sorted(Counter(len(seq) for seq, _ in sub_cycles).items())),
                    sub_alpha,
                    len(sub_alpha) == 2,
                    score_sequence(sub),
                )
            )
        lines.append("deletion profile v:deg,H(T-v),cycles,lengths,alpha,complete,scores")
        for row in deletion_rows:
            lines.append(f"  {row}")
    return "\n".join(lines)


def complete_omega_spectrum(max_n: int = 8) -> str:
    lines = ["Complete-Omega spectrum by isomorphism class", "=" * 52]
    for n in range(3, max_n + 1):
        reps = gentourng_reps(n)
        by_cycles: Counter[int] = Counter()
        examples: dict[int, tuple[int, tuple[int, ...], int]] = {}
        complete_count = 0
        for cid, bits in enumerate(reps):
            rows = rows_from_nauty_bits(n, bits)
            cycles = odd_cycle_masks(rows)
            masks = [mask for _, mask in cycles]
            nonedges = 0
            for a, b in itertools.combinations(masks, 2):
                if a & b == 0:
                    nonedges += 1
                    break
            if nonedges == 0:
                complete_count += 1
                r = len(cycles)
                by_cycles[r] += 1
                examples.setdefault(r, (cid, score_sequence(rows), hamiltonian_paths(rows)))
        lines.append(f"n={n}: complete classes={complete_count}/{len(reps)}")
        if by_cycles:
            support = sorted(by_cycles)
            lines.append(f"  cycle-count support={support}")
            lines.append(f"  missing r in range={ [r for r in range(support[0], support[-1] + 1) if r not in by_cycles] }")
            for r in support:
                cid, scores, h = examples[r]
                lines.append(
                    f"  r={r:2d}: classes={by_cycles[r]:3d}, H=1+2r={1 + 2*r:3d}, "
                    f"example cid={cid}, H={h}, scores={scores}"
                )
    return "\n".join(lines)


def main() -> None:
    print("Omega extremes fingerprint audit (opus-2026-05-29-S11)")
    print("=" * 72)
    print()
    print(complete_omega_spectrum())
    print()
    print(signature_spectrum())
    print()

    n8_reps = gentourng_reps(8)
    for cid in (2519, 3285):
        rows = rows_from_nauty_bits(8, n8_reps[cid])
        aut_labeled = factorial(8)
        print(fingerprint(f"THM-344 n=8 H=63 class cid={cid} (labeled copies={aut_labeled})", rows))
        print()

    print(fingerprint("THM-025 n=9 real-rootedness counterexample", thm025_rows()))


if __name__ == "__main__":
    main()
