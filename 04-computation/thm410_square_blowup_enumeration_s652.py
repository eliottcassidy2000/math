#!/usr/bin/env python3
"""S652: THM-410 and square-node blowups as enumeration speedup carriers.

Two carriers are tested.

1. THM-410 interval matching kernel:
   start from a transitive tournament and reverse a matching M.  The number of
   directed 3-cycles is the additive interval weight sum(b-a-1).  This gives a
   cheap c3-bounded search layer for near-transitive / H-gap work.

2. Square-node substitution:
   replace every vertex of an n-vertex tournament by an n-vertex tournament.
   This produces n^2 vertices, but the imprimitive skeleton can be counted by
   a Burnside quotient over Aut(base), and H can be evaluated by path-cover
   polynomials of the blocks plus a macro-word DP.

The point is not to enumerate all tournaments on n^2 vertices.  The point is to
extract a reusable carrier for the decomposable slice that ordinary canonical
search sees only after paying for the full label set.
"""

from __future__ import annotations

from collections import Counter, defaultdict, deque
from functools import lru_cache
from itertools import combinations, permutations, product
from math import comb, factorial, log2


def pair_index(n: int) -> dict[tuple[int, int], int]:
    return {pair: idx for idx, pair in enumerate(combinations(range(n), 2))}


PAIR_INDEX: dict[int, dict[tuple[int, int], int]] = {}
PAIR_LIST: dict[int, list[tuple[int, int]]] = {}
PERMS: dict[int, list[tuple[int, ...]]] = {}


def pairs(n: int) -> list[tuple[int, int]]:
    if n not in PAIR_LIST:
        PAIR_LIST[n] = list(combinations(range(n), 2))
    return PAIR_LIST[n]


def idx(n: int, i: int, j: int) -> int:
    if n not in PAIR_INDEX:
        PAIR_INDEX[n] = pair_index(n)
    if i > j:
        i, j = j, i
    return PAIR_INDEX[n][(i, j)]


def perms(n: int) -> list[tuple[int, ...]]:
    if n not in PERMS:
        PERMS[n] = list(permutations(range(n)))
    return PERMS[n]


def has_arc(mask: int, n: int, i: int, j: int) -> bool:
    if i == j:
        return False
    if i < j:
        return ((mask >> idx(n, i, j)) & 1) == 1
    return ((mask >> idx(n, j, i)) & 1) == 0


def set_arc_matrix_to_mask(adj: list[list[int]]) -> int:
    n = len(adj)
    out = 0
    for i, j in pairs(n):
        if adj[i][j]:
            out |= 1 << idx(n, i, j)
    return out


def transitive_mask(n: int) -> int:
    return (1 << comb(n, 2)) - 1


def cycle3_count(mask: int, n: int) -> int:
    total = 0
    for a, b, c in combinations(range(n), 3):
        ab = has_arc(mask, n, a, b)
        bc = has_arc(mask, n, b, c)
        ca = has_arc(mask, n, c, a)
        if (ab and bc and ca) or ((not ab) and (not bc) and (not ca)):
            total += 1
    return total


def hamiltonian_count(mask: int, n: int) -> int:
    @lru_cache(maxsize=None)
    def dp(vset: int, last: int) -> int:
        if vset == (1 << last):
            return 1
        prev_set = vset ^ (1 << last)
        total = 0
        for prev in range(n):
            if ((prev_set >> prev) & 1) and has_arc(mask, n, prev, last):
                total += dp(prev_set, prev)
        return total

    full = (1 << n) - 1
    return sum(dp(full, last) for last in range(n))


def score_hist(mask: int, n: int) -> dict[int, int]:
    scores = [0] * n
    for i, j in pairs(n):
        if has_arc(mask, n, i, j):
            scores[i] += 1
        else:
            scores[j] += 1
    return dict(sorted(Counter(scores).items()))


def scc_sizes(mask: int, n: int) -> list[int]:
    graph = [[] for _ in range(n)]
    rev = [[] for _ in range(n)]
    for i, j in pairs(n):
        if has_arc(mask, n, i, j):
            graph[i].append(j)
            rev[j].append(i)
        else:
            graph[j].append(i)
            rev[i].append(j)

    seen = [False] * n
    order: list[int] = []

    def dfs(v: int) -> None:
        seen[v] = True
        for w in graph[v]:
            if not seen[w]:
                dfs(w)
        order.append(v)

    for v in range(n):
        if not seen[v]:
            dfs(v)

    seen = [False] * n
    sizes = []

    def rdfs(v: int) -> int:
        seen[v] = True
        total = 1
        for w in rev[v]:
            if not seen[w]:
                total += rdfs(w)
        return total

    for v in reversed(order):
        if not seen[v]:
            sizes.append(rdfs(v))
    return sorted(sizes, reverse=True)


def permute_mask(mask: int, n: int, perm: tuple[int, ...]) -> int:
    """New label a is old vertex perm[a]."""
    out = 0
    for a, b in pairs(n):
        if has_arc(mask, n, perm[a], perm[b]):
            out |= 1 << idx(n, a, b)
    return out


def canonical_mask(mask: int, n: int) -> int:
    return min(permute_mask(mask, n, p) for p in perms(n))


@lru_cache(maxsize=None)
def tournament_classes(n: int) -> tuple[int, ...]:
    seen = set()
    reps = []
    for mask in range(1 << comb(n, 2)):
        can = canonical_mask(mask, n)
        if can not in seen:
            seen.add(can)
            reps.append(can)
    return tuple(sorted(reps))


def automorphisms(mask: int, n: int) -> list[tuple[int, ...]]:
    return [p for p in perms(n) if permute_mask(mask, n, p) == mask]


def cycle_count_perm(p: tuple[int, ...]) -> int:
    seen = [False] * len(p)
    cycles = 0
    for i in range(len(p)):
        if seen[i]:
            continue
        cycles += 1
        j = i
        while not seen[j]:
            seen[j] = True
            j = p[j]
    return cycles


def matching_weight_distribution(n: int, cap: int | None = None) -> Counter[int]:
    @lru_cache(maxsize=None)
    def rec(mask: int) -> tuple[tuple[int, int], ...]:
        if mask == 0:
            return ((0, 1),)
        i = (mask & -mask).bit_length() - 1
        rest = mask ^ (1 << i)
        out = Counter(dict(rec(rest)))  # leave i unmatched
        jmask = rest
        while jmask:
            j = (jmask & -jmask).bit_length() - 1
            w = j - i - 1
            sub = rest ^ (1 << j)
            for sub_w, count in rec(sub):
                new_w = sub_w + w
                if cap is None or new_w <= cap:
                    out[new_w] += count
            jmask ^= 1 << j
        if cap is not None:
            out = Counter({k: v for k, v in out.items() if k <= cap})
        return tuple(sorted(out.items()))

    return Counter(dict(rec((1 << n) - 1)))


def enumerate_matchings(n: int) -> list[list[tuple[int, int]]]:
    out: list[list[tuple[int, int]]] = []

    def rec(mask: int, edges: list[tuple[int, int]]) -> None:
        if mask == 0:
            out.append(edges[:])
            return
        i = (mask & -mask).bit_length() - 1
        rest = mask ^ (1 << i)
        rec(rest, edges)
        jmask = rest
        while jmask:
            j = (jmask & -jmask).bit_length() - 1
            edges.append((i, j))
            rec(rest ^ (1 << j), edges)
            edges.pop()
            jmask ^= 1 << j

    rec((1 << n) - 1, [])
    return out


def mask_from_reversed_matching(n: int, matching: list[tuple[int, int]]) -> int:
    mask = transitive_mask(n)
    for a, b in matching:
        if a > b:
            a, b = b, a
        mask &= ~(1 << idx(n, a, b))
    return mask


def validate_thm410(max_n: int = 8) -> list[tuple[int, int]]:
    failures = []
    for n in range(1, max_n + 1):
        for matching in enumerate_matchings(n):
            actual = cycle3_count(mask_from_reversed_matching(n, matching), n)
            formula = sum(b - a - 1 for a, b in matching)
            if actual != formula:
                failures.append((n, len(failures)))
                break
    return failures


def ordered_path_cover_counts(mask: int, n: int) -> list[int]:
    """counts[r] = ordered covers of all vertices by r directed paths."""

    @lru_cache(maxsize=None)
    def hp_subset(vset: int, last: int) -> int:
        if vset == (1 << last):
            return 1
        prev_set = vset ^ (1 << last)
        total = 0
        for prev in range(n):
            if ((prev_set >> prev) & 1) and has_arc(mask, n, prev, last):
                total += hp_subset(prev_set, prev)
        return total

    @lru_cache(maxsize=None)
    def path_count(vset: int) -> int:
        return sum(hp_subset(vset, last) for last in range(n) if (vset >> last) & 1)

    full = (1 << n) - 1
    cover: list[Counter[int]] = [Counter() for _ in range(1 << n)]
    cover[0][0] = 1
    for vset in range(1, 1 << n):
        sub = vset
        while sub:
            pc = path_count(sub)
            rest = vset ^ sub
            for r, val in cover[rest].items():
                cover[vset][r + 1] += pc * val
            sub = (sub - 1) & vset
    return [cover[full].get(r, 0) for r in range(n + 1)]


def macro_word_count(base: int, k: int, run_counts: tuple[int, ...]) -> int:
    @lru_cache(maxsize=None)
    def rec(rem: tuple[int, ...], last: int) -> int:
        if sum(rem) == 0:
            return 1
        total = 0
        for j in range(k):
            if rem[j] == 0:
                continue
            if last != k and not has_arc(base, k, last, j):
                continue
            nxt = list(rem)
            nxt[j] -= 1
            total += rec(tuple(nxt), j)
        return total

    return rec(run_counts, k)


def block_h_formula(base: int, k: int, fibers: list[tuple[int, int]]) -> int:
    covers = [ordered_path_cover_counts(mask, size) for mask, size in fibers]
    ranges = [range(1, size + 1) for _, size in fibers]
    total = 0
    for run_counts in product(*ranges):
        macro = macro_word_count(base, k, tuple(run_counts))
        if macro == 0:
            continue
        ways = macro
        for i, r in enumerate(run_counts):
            ways *= covers[i][r]
        total += ways
    return total


def lex_product_mask(base: int, k: int, fibers: list[tuple[int, int]]) -> int:
    sizes = [size for _, size in fibers]
    offsets = [0]
    for size in sizes:
        offsets.append(offsets[-1] + size)
    n = offsets[-1]
    adj = [[0] * n for _ in range(n)]

    for bi, (fiber_mask, size) in enumerate(fibers):
        off = offsets[bi]
        for a, b in pairs(size):
            if has_arc(fiber_mask, size, a, b):
                adj[off + a][off + b] = 1
            else:
                adj[off + b][off + a] = 1

    for bi in range(k):
        for bj in range(k):
            if bi == bj:
                continue
            if has_arc(base, k, bi, bj):
                for a in range(sizes[bi]):
                    for b in range(sizes[bj]):
                        adj[offsets[bi] + a][offsets[bj] + b] = 1
    return set_arc_matrix_to_mask(adj)


def square_template_table(max_n: int = 5) -> list[dict[str, object]]:
    rows = []
    for n in range(1, max_n + 1):
        reps = tournament_classes(n)
        q = len(reps)
        aut_hist = Counter()
        cycle_index_shadow = Counter()
        template_count = 0
        for rep in reps:
            auts = automorphisms(rep, n)
            aut_hist[len(auts)] += 1
            fixed_sum = 0
            for a in auts:
                c = cycle_count_perm(a)
                cycle_index_shadow[c] += 1
                fixed_sum += q**c
            template_count += fixed_sum // len(auts)
        rows.append(
            {
                "n": n,
                "A": q,
                "uniform_pairs": q * q,
                "raw_vertex_assignments": q ** (n + 1),
                "square_templates": template_count,
                "aut_hist": dict(sorted(aut_hist.items())),
                "cycle_shadow": dict(sorted(cycle_index_shadow.items())),
            }
        )
    return rows


def validate_block_formula() -> list[dict[str, object]]:
    t3 = transitive_mask(3)
    c3 = 0
    c3 |= 1 << idx(3, 0, 1)
    c3 |= 1 << idx(3, 1, 2)
    # bit (0,2)=0, so 2 -> 0; this is the directed 3-cycle.
    samples = [
        ("chain3[chain3,chain3,chain3]", t3, [t3, t3, t3]),
        ("cycle3[chain3,chain3,chain3]", c3, [t3, t3, t3]),
        ("cycle3[cycle3,chain3,cycle3]", c3, [c3, t3, c3]),
        ("chain3[cycle3,chain3,cycle3]", t3, [c3, t3, c3]),
    ]
    rows = []
    for name, base, fiber_masks in samples:
        fibers = [(m, 3) for m in fiber_masks]
        product_mask = lex_product_mask(base, 3, fibers)
        direct = hamiltonian_count(product_mask, 9)
        formula = block_h_formula(base, 3, fibers)
        naive = hamiltonian_count(base, 3)
        for m in fiber_masks:
            naive *= hamiltonian_count(m, 3)
        rows.append(
            {
                "name": name,
                "direct_H": direct,
                "block_formula_H": formula,
                "naive_product_H": naive,
                "ok": direct == formula,
            }
        )
    return rows


def method_tournament() -> tuple[list[str], int]:
    methods = [
        ("thm410_interval_matching_kernel", (10, 10, 9, 9)),
        ("block_path_cover_macro_dp", (10, 9, 10, 8)),
        ("substitution_square_cycle_index", (9, 9, 9, 9)),
        ("tiling_succession_dp", (9, 8, 9, 8)),
        ("crt_burnside_partition", (8, 8, 8, 7)),
        ("generic_canonical_augmentation", (7, 7, 8, 6)),
        ("raw_labeled_held_karp", (4, 3, 5, 3)),
    ]
    n = len(methods)
    adj = [[0] * n for _ in range(n)]
    for i, j in pairs(n):
        si = methods[i][1]
        sj = methods[j][1]
        if si >= sj:
            adj[i][j] = 1
        else:
            adj[j][i] = 1
    return [name for name, _ in methods], set_arc_matrix_to_mask(adj)


def print_header(title: str) -> None:
    print("=" * 78)
    print(title)
    print("=" * 78)


def main() -> None:
    print_header("S652 THM-410 / square-node tournament enumeration speedups")

    print("\nTHM-410 validation")
    failures = validate_thm410(8)
    print(f"  exhaustive matching validation through n=8: failures={len(failures)}")
    print("  theorem kernel: c3(T_M) = sum_{(a,b) in M} (b-a-1)")

    print("\nInterval-matching low-c3 carrier")
    print("  c3<=10 is the H=21-relevant cap from H >= 1+2*c3.")
    print(f"  {'n':>2} {'all_matchings':>14} {'c3<=10':>12} {'dist[0..10]':>42}")
    for n in range(3, 15):
        dist = matching_weight_distribution(n, cap=10)
        total = sum(matching_weight_distribution(n).values())
        low = sum(dist.values())
        small = [dist.get(i, 0) for i in range(11)]
        print(f"  {n:>2} {total:>14,} {low:>12,} {str(small):>42}")

    print("\nSquare-substitution Burnside templates")
    print("  A=A000568(n); template = sum_[base B] (1/|Aut(B)|) sum_g A^cycles(g).")
    print("  This counts fiber-class assignments to base vertices, modulo base symmetry.")
    print(f"  {'n':>2} {'A':>5} {'A^2 uniform':>12} {'raw A^(n+1)':>16} {'templates':>16}  aut_order_hist")
    for row in square_template_table(5):
        print(
            f"  {row['n']:>2} {row['A']:>5} {row['uniform_pairs']:>12,} "
            f"{row['raw_vertex_assignments']:>16,} {row['square_templates']:>16,}  "
            f"{row['aut_hist']}"
        )

    print("\nBlock path-cover formula for H on lexicographic products")
    print("  H(base[fibers]) = sum_run_counts macro_words(base,r) * prod_i PathCovers(F_i,r_i).")
    print("  The naive product H(base)*prod H(fiber_i) fails once the outer base has cycles.")
    for row in validate_block_formula():
        print(
            f"  {row['name']:<34} direct={row['direct_H']:<8} "
            f"formula={row['block_formula_H']:<8} naive={row['naive_product_H']:<6} ok={row['ok']}"
        )

    print("\nScale reading")
    for n in range(3, 6):
        q = len(tournament_classes(n))
        templates = square_template_table(n)[-1]["square_templates"]
        raw_labeled_bits = comb(n * n, 2)
        print(
            f"  n={n}: n^2 vertices={n*n}, log2(raw labeled tournaments)={raw_labeled_bits}, "
            f"log2(square templates)={log2(templates):.2f}, A(n)^2={q*q}"
        )

    print("\nTournament Analysis over speedup methods")
    labels, method_mask = method_tournament()
    for i, label in enumerate(labels):
        beaten = [labels[j] for j in range(len(labels)) if i != j and has_arc(method_mask, len(labels), i, j)]
        print(f"  {label:<34} beats={len(beaten)}")
    print(f"  score_hist={score_hist(method_mask, len(labels))}")
    print(f"  directed_3cycles={cycle3_count(method_mask, len(labels))}")
    print(f"  scc_sizes={scc_sizes(method_mask, len(labels))}")
    print(f"  hamiltonian_paths={hamiltonian_count(method_mask, len(labels))}")

    print("\nAssumption challenge")
    print("  Vertices need not be the original tournament vertices.")
    print("  Useful vertex sets here include reversed intervals, block modules, macro-runs,")
    print("  ordered path-cover segments, automorphism cycles, H-values, and proof obligations.")
    print("  THM-410 preserves exact cyclic-triangle witnesses but destroys nonmatching flips.")
    print("  Square substitution preserves modular decomposition and block H-data but destroys")
    print("  prime/indecomposable tournaments outside the imprimitive slice.")


if __name__ == "__main__":
    main()
