#!/usr/bin/env python3
"""S652: THM-410 and square-node blowup enumeration speedups.

This scout tests speedup carriers for tournament isomorphism and H-spectrum
work.  It is intentionally bounded: the goal is to certify reusable formulas
and algorithmic directions, not to run the full A000568(10) census.

Carriers tested here:
  1. THM-410 interval-reversal matchings, where c3 is an additive interval sum.
  2. A broader upset deformation identity around a transitive order.
  3. Lexicographic square substitution Sq(T)=T[T], the user's "each node
     becomes a tournament of the original size" operation.
  4. A module/run-cover dynamic program for H on lexicographic substitutions.
"""

from __future__ import annotations

from collections import Counter
from functools import lru_cache
from itertools import combinations, permutations, product
from math import comb, log2
import random
import time


A000568 = {
    1: 1,
    2: 1,
    3: 2,
    4: 4,
    5: 12,
    6: 56,
    7: 456,
    8: 6880,
    9: 191536,
    10: 9733056,
}


def transitive(n: int) -> list[list[int]]:
    return [[1 if i < j else 0 for j in range(n)] for i in range(n)]


def adj_from_bits(n: int, bits: int) -> list[list[int]]:
    """Upper-pair bit 1 means i->j for i<j; bit 0 means j->i."""
    adj = [[0] * n for _ in range(n)]
    bit = 0
    for i in range(n):
        for j in range(i + 1, n):
            if (bits >> bit) & 1:
                adj[i][j] = 1
            else:
                adj[j][i] = 1
            bit += 1
    return adj


def encode(adj: list[list[int]], perm: tuple[int, ...] | None = None) -> int:
    n = len(adj)
    if perm is None:
        perm = tuple(range(n))
    code = 0
    bit = 0
    for a in range(n):
        u = perm[a]
        for b in range(a + 1, n):
            v = perm[b]
            if adj[u][v]:
                code |= 1 << bit
            bit += 1
    return code


def canon(adj: list[list[int]]) -> int:
    return min(encode(adj, p) for p in permutations(range(len(adj))))


def iso_reps(n: int) -> list[list[list[int]]]:
    reps: dict[int, list[list[int]]] = {}
    for bits in range(1 << comb(n, 2)):
        adj = adj_from_bits(n, bits)
        key = canon(adj)
        reps.setdefault(key, adj_from_bits(n, key))
    return [reps[key] for key in sorted(reps)]


def scores(adj: list[list[int]]) -> list[int]:
    return [sum(row) for row in adj]


def score_hist(adj: list[list[int]]) -> Counter[int]:
    return Counter(scores(adj))


def hist_tuple(hist: Counter[int]) -> tuple[tuple[int, int], ...]:
    return tuple(sorted(hist.items()))


def c3_count(adj: list[list[int]]) -> int:
    n = len(adj)
    total = 0
    for a, b, c in combinations(range(n), 3):
        da = adj[a][b] + adj[a][c]
        db = adj[b][a] + adj[b][c]
        dc = adj[c][a] + adj[c][b]
        if da == db == dc == 1:
            total += 1
    return total


def hamiltonian_paths(adj: list[list[int]]) -> int:
    n = len(adj)
    size = 1 << n
    dp = [[0] * n for _ in range(size)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(size):
        for v in range(n):
            cur = dp[mask][v]
            if not cur:
                continue
            for w in range(n):
                if not ((mask >> w) & 1) and adj[v][w]:
                    dp[mask | (1 << w)][w] += cur
    return sum(dp[-1])


def subset_h_counts(adj: list[list[int]]) -> list[int]:
    """h[mask] = Hamiltonian paths of the induced subtournament on mask."""
    n = len(adj)
    size = 1 << n
    dp = [[0] * n for _ in range(size)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(size):
        for v in range(n):
            cur = dp[mask][v]
            if not cur:
                continue
            for w in range(n):
                if not ((mask >> w) & 1) and adj[v][w]:
                    dp[mask | (1 << w)][w] += cur
    return [sum(dp[mask]) for mask in range(size)]


@lru_cache(maxsize=None)
def ordered_path_cover_counts_by_code(n: int, code: int) -> tuple[int, ...]:
    """cover[r] = ordered r-run directed path covers of all n vertices."""
    adj = adj_from_bits(n, code)
    size = 1 << n
    h_sub = subset_h_counts(adj)
    cover = [[0] * size for _ in range(n + 1)]
    cover[0][0] = 1
    for r in range(1, n + 1):
        for mask in range(size):
            sub = mask
            total = 0
            while sub:
                total += h_sub[sub] * cover[r - 1][mask ^ sub]
                sub = (sub - 1) & mask
            cover[r][mask] = total
    return tuple(cover[r][size - 1] for r in range(n + 1))


def macro_word_counter(outer: list[list[int]]):
    k = len(outer)

    @lru_cache(maxsize=None)
    def rec(last: int, remaining: tuple[int, ...]) -> int:
        if sum(remaining) == 0:
            return 1
        total = 0
        for nxt in range(k):
            if remaining[nxt] == 0:
                continue
            if last != -1 and not outer[last][nxt]:
                continue
            new_remaining = list(remaining)
            new_remaining[nxt] -= 1
            total += rec(nxt, tuple(new_remaining))
        return total

    return rec


def lex_h_module_dp(
    outer: list[list[int]], inners: tuple[list[list[int]], ...]
) -> int:
    """Count H of a lexicographic substitution using module run covers."""
    cover_lists = [
        ordered_path_cover_counts_by_code(len(inner), encode(inner))
        for inner in inners
    ]
    macro = macro_word_counter(outer)
    ranges = [range(1, len(inner) + 1) for inner in inners]
    total = 0
    for run_counts in product(*ranges):
        macro_count = macro(-1, tuple(run_counts))
        if not macro_count:
            continue
        ways = macro_count
        for idx, r in enumerate(run_counts):
            ways *= cover_lists[idx][r]
        total += ways
    return total


def scc_sizes(adj: list[list[int]]) -> list[int]:
    n = len(adj)

    def dfs(start: int, reverse: bool) -> set[int]:
        seen = {start}
        stack = [start]
        while stack:
            v = stack.pop()
            for w in range(n):
                edge = adj[w][v] if reverse else adj[v][w]
                if edge and w not in seen:
                    seen.add(w)
                    stack.append(w)
        return seen

    remaining = set(range(n))
    sizes = []
    while remaining:
        v = next(iter(remaining))
        comp = dfs(v, False) & dfs(v, True)
        sizes.append(len(comp))
        remaining -= comp
    return sorted(sizes, reverse=True)


def all_matchings(n: int) -> list[tuple[tuple[int, int], ...]]:
    @lru_cache(maxsize=None)
    def rec(unused: tuple[int, ...]) -> tuple[tuple[tuple[int, int], ...], ...]:
        if not unused:
            return ((),)
        first = unused[0]
        rest = unused[1:]
        out: list[tuple[tuple[int, int], ...]] = list(rec(rest))
        for idx, other in enumerate(rest):
            tail_unused = rest[:idx] + rest[idx + 1 :]
            for tail in rec(tail_unused):
                out.append(((first, other),) + tail)
        return tuple(out)

    return list(rec(tuple(range(n))))


def reverse_matching(n: int, matching: tuple[tuple[int, int], ...]) -> list[list[int]]:
    adj = transitive(n)
    for a, b in matching:
        adj[a][b] = 0
        adj[b][a] = 1
    return adj


def thm410_formula(matching: tuple[tuple[int, int], ...]) -> int:
    return sum(b - a - 1 for a, b in matching)


def upset_c3_formula(adj: list[list[int]]) -> int:
    """Exact c3 count relative to the natural transitive order.

    For x<y<z, let r_ab be 1 iff the transitive edge a->b is reversed.
    The triple is cyclic exactly when r_xy = r_yz != r_xz.
    """
    total = 0
    for x, y, z in combinations(range(len(adj)), 3):
        rxy = adj[y][x]
        ryz = adj[z][y]
        rxz = adj[z][x]
        if rxy == ryz and rxz != rxy:
            total += 1
    return total


def upset_c3_bitset(adj: list[list[int]]) -> int:
    n = len(adj)
    outmask = [0] * n
    for i in range(n):
        mask = 0
        for j in range(n):
            if adj[i][j]:
                mask |= 1 << j
        outmask[i] = mask

    total = 0
    for y in range(n):
        left0 = left1 = right0 = right1 = 0
        for x in range(y):
            if adj[y][x]:
                left1 |= 1 << x
            else:
                left0 |= 1 << x
        for z in range(y + 1, n):
            if adj[z][y]:
                right1 |= 1 << z
            else:
                right0 |= 1 << z
        lx = left0
        while lx:
            x = (lx & -lx).bit_length() - 1
            total += (right0 & ~outmask[x]).bit_count()
            lx &= lx - 1
        lx = left1
        while lx:
            x = (lx & -lx).bit_length() - 1
            total += (right1 & outmask[x]).bit_count()
            lx &= lx - 1
    return total


def lex_product(
    outer: list[list[int]], inners: tuple[list[list[int]], ...]
) -> list[list[int]]:
    sizes = [len(inner) for inner in inners]
    offsets = [0]
    for size in sizes:
        offsets.append(offsets[-1] + size)
    total = offsets[-1]
    adj = [[0] * total for _ in range(total)]

    for block, inner in enumerate(inners):
        off = offsets[block]
        for a in range(sizes[block]):
            for b in range(sizes[block]):
                if inner[a][b]:
                    adj[off + a][off + b] = 1

    for i in range(len(outer)):
        for j in range(len(outer)):
            if i == j or not outer[i][j]:
                continue
            oi = offsets[i]
            oj = offsets[j]
            for a in range(sizes[i]):
                for b in range(sizes[j]):
                    adj[oi + a][oj + b] = 1
    return adj


def c3_composition_formula(
    outer: list[list[int]], inners: tuple[list[list[int]], ...]
) -> int:
    sizes = [len(inner) for inner in inners]
    total = sum(c3_count(inner) for inner in inners)
    for i, j, k in combinations(range(len(outer)), 3):
        deg_i = outer[i][j] + outer[i][k]
        deg_j = outer[j][i] + outer[j][k]
        deg_k = outer[k][i] + outer[k][j]
        if deg_i == deg_j == deg_k == 1:
            total += sizes[i] * sizes[j] * sizes[k]
    return total


def square_score_hist_formula(adj: list[list[int]]) -> Counter[int]:
    n = len(adj)
    hist = Counter(scores(adj))
    out = Counter()
    for so, mo in hist.items():
        for si, mi in hist.items():
            out[n * so + si] += mo * mi
    return out


def method_tournament() -> tuple[list[str], list[list[int]]]:
    methods = [
        ("raw_labelled_scan", (0, 5, 0, 0, 5)),
        ("iso_canonical_augmentation", (2, 5, 2, 2, 4)),
        ("Burnside_cycle_type_count", (5, 5, 1, 4, 5)),
        ("c3_pruned_extension", (3, 5, 3, 3, 4)),
        ("THM410_interval_matching_cache", (4, 5, 4, 4, 5)),
        ("upset_bitset_deformation_formula", (4, 5, 4, 4, 4)),
        ("modular_substitution_decomposition", (5, 4, 5, 5, 4)),
        ("lexicographic_square_subsequence", (4, 4, 4, 5, 5)),
        ("quotient_word_H_automaton", (4, 3, 5, 5, 3)),
    ]
    n = len(methods)
    adj = [[0] * n for _ in range(n)]
    for i, j in combinations(range(n), 2):
        wi = wj = 0
        for a, b in zip(methods[i][1], methods[j][1]):
            if a > b:
                wi += 1
            elif b > a:
                wj += 1
        if wi == wj:
            wi += methods[i][1][2] >= methods[j][1][2]
            wj += methods[j][1][2] > methods[i][1][2]
        if wi > wj:
            adj[i][j] = 1
        else:
            adj[j][i] = 1
    return [name for name, _ in methods], adj


def sci_int(x: int) -> str:
    return f"{x:.3e}" if x >= 1_000_000 else str(x)


def direct_states(total_vertices: int) -> int:
    return total_vertices * (1 << total_vertices)


def square_module_state_proxy(n: int) -> int:
    # Block path-cover states plus macro states over run-count boxes.
    return n * (3**n) + n * ((n + 1) ** n)


def print_lines(lines: list[str], title: str) -> None:
    lines.append("")
    lines.append(title)
    lines.append("-" * len(title))


def main() -> None:
    lines: list[str] = []
    lines.append("S652 THM-410 / Square-Blowup Enumeration Speedups")
    lines.append("=" * 72)
    lines.append("")
    lines.append("Thesis: enumeration gets faster when the quotient remembers the")
    lines.append("witness carrier.  THM-410 remembers reversed intervals; square-node")
    lines.append("substitution remembers modules.  Forget either carrier and the task")
    lines.append("falls back toward raw edge masks or raw S_N canonicalization.")

    print_lines(lines, "THM-410 Matching Cache")
    for n in range(2, 9):
        dist = Counter()
        mismatches = 0
        for matching in all_matchings(n):
            got = c3_count(reverse_matching(n, matching))
            want = thm410_formula(matching)
            if got != want:
                mismatches += 1
            dist[want] += 1
        preview = ", ".join(f"{k}:{dist[k]}" for k in sorted(dist)[:8])
        lines.append(
            f"n={n}: matchings={sum(dist.values())}, mismatches={mismatches}, "
            f"first c3 bins={preview}"
        )

    print_lines(lines, "Interval-Matching Search Lane")
    lines.append(
        "The cap c3<=B becomes the additive bound sum(b-a-1)<=B over "
        "reversed matching intervals."
    )
    lines.append(" n  matchings  labelled tournaments  A000568(n)  c3-values  <=10  <=20")
    for n in range(5, 11):
        dist = Counter(thm410_formula(m) for m in all_matchings(n))
        lines.append(
            f"{n:2d} {sum(dist.values()):10d} {sci_int(1 << comb(n, 2)):>21} "
            f"{str(A000568.get(n, '?')):>11} {len(dist):10d} "
            f"{sum(v for k, v in dist.items() if k <= 10):5d} "
            f"{sum(v for k, v in dist.items() if k <= 20):5d}"
        )

    print_lines(lines, "Upset Deformation Formula")
    for n in range(3, 7):
        total = 1 << comb(n, 2)
        mismatches = 0
        for bits in range(total):
            adj = adj_from_bits(n, bits)
            c3 = c3_count(adj)
            if c3 != upset_c3_formula(adj) or c3 != upset_c3_bitset(adj):
                mismatches += 1
        lines.append(f"n={n}: labelled={total}, mismatches={mismatches}")
    rng = random.Random(20260605)
    adj96 = transitive(96)
    for i in range(96):
        for j in range(i + 1, 96):
            if rng.random() < 0.22:
                adj96[i][j] = 0
                adj96[j][i] = 1
            else:
                adj96[i][j] = 1
                adj96[j][i] = 0
    t0 = time.perf_counter()
    slow = upset_c3_formula(adj96)
    t1 = time.perf_counter()
    fast = upset_c3_bitset(adj96)
    t2 = time.perf_counter()
    lines.append(
        f"n=96 random upset: triple-scan={t1-t0:.4f}s, "
        f"bitset={t2-t1:.4f}s, equal={slow == fast}"
    )
    lines.append("Formula: for x<y<z, cyclic iff r_xy=r_yz and r_xz is opposite.")

    print_lines(lines, "Lexicographic Square Sq(T)=T[T]")
    for n in range(2, 6):
        total = 1 << comb(n, 2)
        mismatches = 0
        canon_classes = set()
        square_fingerprints = set()
        c3_values = set()
        for bits in range(total):
            adj = adj_from_bits(n, bits)
            canon_classes.add(canon(adj))
            sq = lex_product(adj, tuple(adj for _ in range(n)))
            want_c3 = (n**3 + n) * c3_count(adj)
            want_scores = square_score_hist_formula(adj)
            if c3_count(sq) != want_c3 or score_hist(sq) != want_scores:
                mismatches += 1
            square_fingerprints.add((hist_tuple(want_scores), want_c3))
            c3_values.add(c3_count(adj))
        raw_bits = comb(n * n, 2)
        iso = len(canon_classes)
        lines.append(
            f"n={n}->N={n*n}: labelled={total}, iso={iso}, "
            f"mismatches={mismatches}, c3_multiplier={n**3+n}"
        )
        lines.append(
            f"  distinct base c3={len(c3_values)}, square(score,c3) "
            f"fingerprints={len(square_fingerprints)}, structured images<={iso} "
            f"inside 2^{raw_bits} labelled N-space"
        )
    lines.append("Exact: score(i,a)=n*score_T(i)+score_T(a); c3(Sq(T))=(n^3+n)c3(T).")

    print_lines(lines, "Module H Dynamic Program")
    reps3 = iso_reps(3)
    failures = 0
    checked = 0
    for outer in reps3:
        for inners in product(reps3, repeat=3):
            lex = lex_product(outer, tuple(inners))
            direct_h = hamiltonian_paths(lex)
            module_h = lex_h_module_dp(outer, tuple(inners))
            direct_c3 = c3_count(lex)
            formula_c3 = c3_composition_formula(outer, tuple(inners))
            checked += 1
            if direct_h != module_h or direct_c3 != formula_c3:
                failures += 1
    lines.append(f"n=3 substitution combos checked directly: {checked}, failures={failures}")

    reps4 = iso_reps(4)
    sample4 = []
    for outer in reps4:
        sample4.append((outer, tuple([reps4[0]] * 4)))
        sample4.append((outer, tuple([reps4[-1]] * 4)))
    failures4 = 0
    for outer, inners in sample4:
        lex = lex_product(outer, inners)
        if lex_h_module_dp(outer, inners) != hamiltonian_paths(lex):
            failures4 += 1
        if c3_composition_formula(outer, inners) != c3_count(lex):
            failures4 += 1
    lines.append(f"n=4 total-16-vertex samples checked directly: {len(sample4)}, failures={failures4}")

    c3_outer = [[0, 1, 0], [0, 0, 1], [1, 0, 0]]
    t2 = transitive(2)
    c3_square = lex_product(c3_outer, tuple(c3_outer for _ in range(3)))
    lines.append(
        "C3[C3]: "
        f"H_module={lex_h_module_dp(c3_outer, tuple(c3_outer for _ in range(3)))}, "
        f"H_direct={hamiltonian_paths(c3_square)}, "
        f"naive_multiplicative={hamiltonian_paths(c3_outer) ** 4}, "
        f"c3={c3_count(c3_square)}"
    )
    c3_t2 = lex_product(c3_outer, tuple(t2 for _ in range(3)))
    lines.append(
        "C3[T2,T2,T2]: "
        f"H={hamiltonian_paths(c3_t2)}, "
        f"naive_multiplicative={hamiltonian_paths(c3_outer) * hamiltonian_paths(t2) ** 3}"
    )
    lines.append(
        "Consequence: H multiplies safely over transitive condensation/direct sums, "
        "but strong outer quotients need macro-word interleavings."
    )

    print_lines(lines, "State-Count Proxy")
    lines.append(" n  total vertices  raw Held-Karp states  module-state proxy  raw/module")
    for n in range(3, 8):
        raw = direct_states(n * n)
        mod = square_module_state_proxy(n)
        lines.append(f"{n:2d} {n*n:15d} {sci_int(raw):>22} {sci_int(mod):>19} {raw/mod:11.2f}")

    print_lines(lines, "Tournament Analysis")
    names, mt = method_tournament()
    order = sorted(range(len(names)), key=lambda i: -sum(mt[i]))
    lines.append(f"vertices={', '.join(names)}")
    lines.append(f"score_hist={dict(sorted(score_hist(mt).items()))}")
    lines.append(f"c3={c3_count(mt)}")
    lines.append(f"scc_sizes={scc_sizes(mt)}")
    lines.append(f"H={hamiltonian_paths(mt)}")
    lines.append("ranking:")
    for idx in order:
        lines.append(f"  {sum(mt[idx])}: {names[idx]}")
    lines.append(
        "Challenged assumption: vertices need not be tournaments themselves; "
        "proof routes, interval witnesses, quotient states, modular fibers, "
        "and run-cover packets can be the more faithful vertex set."
    )

    print_lines(lines, "Actionable Synthesis")
    lines.append(
        "1. For low-c3 and near-transitive probes, seed by THM-410 intervals "
        "and then widen to general upsets using r_xy=r_yz!=r_xz bitsets."
    )
    lines.append(
        "2. For A000568-like enumeration, detect modular substitution before "
        "full canonicalization: prime quotient plus decorated fibers should "
        "replace many raw S_N searches."
    )
    lines.append(
        "3. For the user's square operation, Sq(T)=T[T] maps n to n^2 while "
        "keeping score and c3 in base-n arithmetic; H is handled by run-cover "
        "polynomials and macro words."
    )
    lines.append(
        "4. For the n=10 H-spectrum bottleneck, S9 already showed all 157 "
        "high gaps unlock by biased n=10 witnesses.  The next exact win is a "
        "certified structured-witness menu plus module/c3 pruning, not blind "
        "flattening of every side channel."
    )

    print("\n".join(lines))


if __name__ == "__main__":
    main()
