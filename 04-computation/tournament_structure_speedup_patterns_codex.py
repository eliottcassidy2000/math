#!/usr/bin/env python3
"""Efficient tournament structure computations and the first trace correction.

This scout has two goals:

1. Turn the recent "efficiency becomes proof" theme into reusable formulas.
   For a tournament adjacency matrix A, directed k-cycles satisfy

       c_k = tr(A^k) / k,  for k = 3,4,5.

   The reason is structural, not accidental: tournaments have no loops or
   directed 2-cycles, so no non-simple closed walk of length < 6 exists.
   Length 6 is the first failure, and the non-simple closed walks are exactly
   rotated figure-eights: either one directed triangle twice, or two directed
   triangles with nonempty vertex intersection.  Therefore

       c_6 = (tr(A^6) - 3*c_3 - 6*p_33^meet) / 6,

   where p_33^meet is the number of unordered pairs of distinct directed
   triangles with nonempty intersection.  This is the first place where trace
   computation must retain cycle-placement data.

2. Ask which fast scalar invariants explain each other.  Exhaustively over all
   labelled n=6 tournaments, we build an information tournament on invariants.

Tournament Analysis
-------------------
Vertices: invariants `score`, `c3`, `c4`, `c5`, `c6`, `H`.
Pairwise observable: normalized mutual explanation
    U(Y|X) = 1 - entropy(Y|X)/entropy(Y).
Gauge: orient X -> Y when X explains a larger fraction of Y than Y explains X.
Tie Hamiltonian path: score -> c3 -> c4 -> c5 -> c6 -> H.
Reported fingerprints: scores, score histogram, directed 3-cycles, SCC sizes,
Hamiltonian-path count, and edge list.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from itertools import combinations, permutations
from math import comb, log2
from random import Random
from time import perf_counter


INFO_ATTRS = ("score", "c3", "c4", "c5", "c6", "H")


def out_bits_from_mask(n: int, mask: int) -> list[int]:
    out = [0] * n
    bit = 0
    for i in range(n):
        for j in range(i + 1, n):
            if (mask >> bit) & 1:
                out[i] |= 1 << j
            else:
                out[j] |= 1 << i
            bit += 1
    return out


def random_out_bits(n: int, rng: Random) -> list[int]:
    out = [0] * n
    for i in range(n):
        for j in range(i + 1, n):
            if rng.getrandbits(1):
                out[i] |= 1 << j
            else:
                out[j] |= 1 << i
    return out


def adjacency_matrix(out: list[int]) -> list[list[int]]:
    n = len(out)
    return [[1 if (out[i] >> j) & 1 else 0 for j in range(n)] for i in range(n)]


def matmul(a: list[list[int]], b: list[list[int]]) -> list[list[int]]:
    n = len(a)
    out = [[0] * n for _ in range(n)]
    for i in range(n):
        row = out[i]
        for k in range(n):
            aik = a[i][k]
            if not aik:
                continue
            bk = b[k]
            for j in range(n):
                row[j] += aik * bk[j]
    return out


def trace_powers(out: list[int], max_power: int = 6) -> dict[int, int]:
    a = adjacency_matrix(out)
    power = a
    traces: dict[int, int] = {1: sum(power[i][i] for i in range(len(a)))}
    for k in range(2, max_power + 1):
        power = matmul(power, a)
        traces[k] = sum(power[i][i] for i in range(len(a)))
    return traces


def triangle_counts_through_vertices(out: list[int]) -> list[int]:
    n = len(out)
    all_mask = (1 << n) - 1
    counts: list[int] = []
    for v in range(n):
        in_v = all_mask ^ out[v] ^ (1 << v)
        total = 0
        x = out[v]
        while x:
            low = x & -x
            a = low.bit_length() - 1
            total += (out[a] & in_v).bit_count()
            x ^= low
        counts.append(total)
    return counts


def directed_triangle_masks(out: list[int]) -> list[int]:
    n = len(out)
    masks: list[int] = []
    for i in range(n):
        for j in range(i + 1, n):
            for k in range(j + 1, n):
                if ((out[i] >> j) & 1) and ((out[j] >> k) & 1) and ((out[k] >> i) & 1):
                    masks.append((1 << i) | (1 << j) | (1 << k))
                elif ((out[i] >> k) & 1) and ((out[k] >> j) & 1) and ((out[j] >> i) & 1):
                    masks.append((1 << i) | (1 << j) | (1 << k))
    return masks


def disjoint_triangle_pair_count(tri_masks: list[int]) -> int:
    total = 0
    for i, left in enumerate(tri_masks):
        for right in tri_masks[i + 1 :]:
            if left & right == 0:
                total += 1
    return total


def score_sequence(out: list[int]) -> tuple[int, ...]:
    return tuple(sorted(bits.bit_count() for bits in out))


def c3_from_score(out: list[int]) -> int:
    n = len(out)
    return comb(n, 3) - sum(comb(d, 2) for d in (bits.bit_count() for bits in out))


def fast_cycle_counts(out: list[int]) -> dict[str, int]:
    traces = trace_powers(out, 6)
    tau = triangle_counts_through_vertices(out)
    tri_masks = directed_triangle_masks(out)
    c3 = len(tri_masks)
    p33_disjoint = disjoint_triangle_pair_count(tri_masks)
    p33_meet = comb(c3, 2) - p33_disjoint
    return {
        "c3": traces[3] // 3,
        "c3_score": c3_from_score(out),
        "c4": traces[4] // 4,
        "c5": traces[5] // 5,
        "c6": (traces[6] - 3 * c3 - 6 * p33_meet) // 6,
        "trace6": traces[6],
        "half_return_pairs": sum(t * t for t in tau),
        "p33_disjoint": p33_disjoint,
        "p33_meet": p33_meet,
        "tau": tuple(tau),
    }


def brute_simple_cycle_count(out: list[int], k: int) -> int:
    n = len(out)
    total = 0
    for subset in combinations(range(n), k):
        root = subset[0]
        rest = subset[1:]
        for perm in permutations(rest):
            path = (root,) + perm
            ok = True
            for i in range(k - 1):
                if not ((out[path[i]] >> path[i + 1]) & 1):
                    ok = False
                    break
            if ok and ((out[path[-1]] >> root) & 1):
                total += 1
    return total


def hamiltonian_paths_dp(out: list[int]) -> int:
    n = len(out)
    size = 1 << n
    dp = [[0] * n for _ in range(size)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(size):
        for last in range(n):
            count = dp[mask][last]
            if not count:
                continue
            choices = out[last] & ~mask
            while choices:
                low = choices & -choices
                nxt = low.bit_length() - 1
                dp[mask | low][nxt] += count
                choices ^= low
    return sum(dp[-1])


def hamiltonian_paths_brute(out: list[int]) -> int:
    n = len(out)
    total = 0
    for perm in permutations(range(n)):
        if all((out[perm[i]] >> perm[i + 1]) & 1 for i in range(n - 1)):
            total += 1
    return total


def entropy(values: list[object]) -> float:
    n = len(values)
    counts = Counter(values)
    return -sum((c / n) * log2(c / n) for c in counts.values())


def conditional_entropy(records: list[dict[str, object]], x: str, y: str) -> float:
    n = len(records)
    buckets: dict[object, list[object]] = defaultdict(list)
    for rec in records:
        buckets[rec[x]].append(rec[y])
    total = 0.0
    for vals in buckets.values():
        total += (len(vals) / n) * entropy(vals)
    return total


def uncertainty(records: list[dict[str, object]], x: str, y: str) -> float:
    hy = entropy([rec[y] for rec in records])
    if hy == 0:
        return 1.0
    return 1.0 - conditional_entropy(records, x, y) / hy


def information_tournament(records: list[dict[str, object]]) -> dict[str, object]:
    attrs = list(INFO_ATTRS)
    tie_rank = {attr: i for i, attr in enumerate(attrs)}
    adj = {a: set() for a in attrs}
    strengths = {}
    for a, b in combinations(attrs, 2):
        ab = uncertainty(records, a, b)
        ba = uncertainty(records, b, a)
        strengths[(a, b)] = ab
        strengths[(b, a)] = ba
        if ab > ba:
            adj[a].add(b)
        elif ba > ab:
            adj[b].add(a)
        elif tie_rank[a] < tie_rank[b]:
            adj[a].add(b)
        else:
            adj[b].add(a)
    return {
        "attrs": attrs,
        "adj": adj,
        "strengths": strengths,
    }


def directed_3_cycles(vertices: list[str], adj: dict[str, set[str]]) -> int:
    total = 0
    for a, b, c in combinations(vertices, 3):
        if b in adj[a] and c in adj[b] and a in adj[c]:
            total += 1
        elif c in adj[a] and b in adj[c] and a in adj[b]:
            total += 1
    return total


def scc_sizes(vertices: list[str], adj: dict[str, set[str]]) -> list[int]:
    radj = {v: set() for v in vertices}
    for v in vertices:
        for w in adj[v]:
            radj[w].add(v)

    seen: set[str] = set()
    order: list[str] = []

    def dfs(v: str) -> None:
        seen.add(v)
        for w in adj[v]:
            if w not in seen:
                dfs(w)
        order.append(v)

    for v in vertices:
        if v not in seen:
            dfs(v)

    seen.clear()
    sizes: list[int] = []

    def rdfs(v: str) -> int:
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


def hamiltonian_path_count_tournament(vertices: list[str], adj: dict[str, set[str]]) -> int:
    idx = {v: i for i, v in enumerate(vertices)}
    n = len(vertices)
    dp: dict[tuple[int, str], int] = {}
    for v in vertices:
        dp[(1 << idx[v], v)] = 1
    for mask in range(1 << n):
        for last in vertices:
            count = dp.get((mask, last), 0)
            if not count:
                continue
            for nxt in adj[last]:
                bit = 1 << idx[nxt]
                if mask & bit:
                    continue
                dp[(mask | bit, nxt)] = dp.get((mask | bit, nxt), 0) + count
    full = (1 << n) - 1
    return sum(dp.get((full, v), 0) for v in vertices)


def exhaustive_n6_records() -> list[dict[str, object]]:
    n = 6
    records: list[dict[str, object]] = []
    for mask in range(1 << comb(n, 2)):
        out = out_bits_from_mask(n, mask)
        counts = fast_cycle_counts(out)
        records.append(
            {
                "score": score_sequence(out),
                "c3": counts["c3"],
                "c4": counts["c4"],
                "c5": counts["c5"],
                "c6": counts["c6"],
                "H": hamiltonian_paths_dp(out),
            }
        )
    return records


def mixed_bucket_report(records: list[dict[str, object]], keys: tuple[str, ...], target: str) -> dict[str, int]:
    buckets: dict[tuple[object, ...], set[object]] = defaultdict(set)
    sizes: Counter[tuple[object, ...]] = Counter()
    for rec in records:
        key = tuple(rec[k] for k in keys)
        buckets[key].add(rec[target])
        sizes[key] += 1
    mixed = {k: vals for k, vals in buckets.items() if len(vals) > 1}
    return {
        "buckets": len(buckets),
        "mixed_buckets": len(mixed),
        "records_in_mixed_buckets": sum(sizes[k] for k in mixed),
        "max_target_values_in_bucket": max((len(vals) for vals in mixed.values()), default=1),
    }


def benchmark() -> None:
    rng = Random(2495)
    out14 = random_out_bits(14, rng)
    out8 = random_out_bits(8, rng)

    t0 = perf_counter()
    fast5 = fast_cycle_counts(out14)["c5"]
    t1 = perf_counter()
    brute5 = brute_simple_cycle_count(out14, 5)
    t2 = perf_counter()

    fast_counts = fast_cycle_counts(out14)
    fast6 = fast_counts["c6"]
    t3 = perf_counter()
    brute6 = brute_simple_cycle_count(out14, 6)
    t4 = perf_counter()

    h0 = perf_counter()
    h_dp = hamiltonian_paths_dp(out8)
    h1 = perf_counter()
    h_brute = hamiltonian_paths_brute(out8)
    h2 = perf_counter()

    print("\n=== Benchmarks on fixed random tournaments ===")
    print(f"n=14 c5 trace={fast5}, brute={brute5}, match={fast5 == brute5}")
    print(f"  trace_time={t1 - t0:.6f}s, brute_time={t2 - t1:.6f}s, speedup={(t2 - t1) / max(t1 - t0, 1e-12):.1f}x")
    print(f"n=14 c6 corrected={fast6}, brute={brute6}, match={fast6 == brute6}")
    print(f"  corrected_time={t3 - t2:.6f}s, brute_time={t4 - t3:.6f}s, speedup={(t4 - t3) / max(t3 - t2, 1e-12):.1f}x")
    print(
        "  trace6_decomposition: "
        f"trace6={fast_counts['trace6']}, 3*c3={3 * fast_counts['c3']}, "
        f"6*p33_meet={6 * fast_counts['p33_meet']}, "
        f"6*c6={6 * fast6}"
    )
    print(f"  diagnostic: half-return rooted triangle-pairs sum tau_v^2={fast_counts['half_return_pairs']} (insufficient alone)")
    print(f"n=8 H dp={h_dp}, brute={h_brute}, match={h_dp == h_brute}")
    print(f"  dp_time={h1 - h0:.6f}s, brute_time={h2 - h1:.6f}s, speedup={(h2 - h1) / max(h1 - h0, 1e-12):.1f}x")


def validate_trace_formulas() -> None:
    print("=== Formula validation ===")
    mismatches = Counter()
    checked = 0
    for n in range(3, 7):
        for mask in range(1 << comb(n, 2)):
            out = out_bits_from_mask(n, mask)
            fast = fast_cycle_counts(out)
            if fast["c3"] != fast["c3_score"]:
                mismatches["c3_score"] += 1
            for k in (3, 4, 5, 6):
                brute = brute_simple_cycle_count(out, k) if n >= k else 0
                if fast[f"c{k}"] != brute:
                    mismatches[f"c{k}"] += 1
            checked += 1
    print(f"exhaustive n=3..6 checked={checked}, mismatches={dict(mismatches)}")

    rng = Random(917)
    random_mismatches = Counter()
    for n in (7, 8, 9):
        for _ in range(80):
            out = random_out_bits(n, rng)
            fast = fast_cycle_counts(out)
            for k in (3, 4, 5, 6):
                brute = brute_simple_cycle_count(out, k)
                if fast[f"c{k}"] != brute:
                    random_mismatches[(n, k)] += 1
    print(f"random n=7..9 samples per n=80, mismatches={dict(random_mismatches)}")


def print_information_tournament(records: list[dict[str, object]]) -> None:
    info = information_tournament(records)
    attrs = info["attrs"]
    adj = info["adj"]
    strengths = info["strengths"]
    outdegrees = {a: len(adj[a]) for a in attrs}
    score_hist = Counter(outdegrees.values())
    edges = []
    for a, b in combinations(attrs, 2):
        if b in adj[a]:
            edges.append(f"{a}>{b}")
        else:
            edges.append(f"{b}>{a}")
    print("\n=== Tournament Analysis: invariant mutual-explanation on all n=6 tournaments ===")
    print("vertices:", attrs)
    print("observable: U(Y|X)=1-H(Y|X)/H(Y), orient X->Y when X explains Y better")
    print("tie_path: score -> c3 -> c4 -> c5 -> c6 -> H")
    print("outdegrees:", outdegrees)
    print("score_hist:", dict(sorted(score_hist.items())))
    print("directed_3_cycles:", directed_3_cycles(attrs, adj))
    print("scc_sizes:", scc_sizes(attrs, adj))
    print("hamiltonian_paths:", hamiltonian_path_count_tournament(attrs, adj))
    print("edges:", ", ".join(edges))
    print("selected explanation strengths:")
    for a, b in (("score", "c3"), ("score", "c5"), ("c5", "H"), ("c6", "H"), ("H", "c6")):
        print(f"  U({b}|{a})={strengths[(a, b)]:.4f}")


def print_exhaustive_patterns(records: list[dict[str, object]]) -> None:
    print("\n=== Exhaustive n=6 structural buckets ===")
    print(f"records={len(records)}")
    for attr in INFO_ATTRS:
        vals = [rec[attr] for rec in records]
        print(f"{attr}: distinct={len(set(vals))}, entropy={entropy(vals):.4f}")
    for keys in (
        ("score",),
        ("score", "c5"),
        ("score", "c5", "c6"),
        ("c3", "c4", "c5", "c6"),
        ("score", "c3", "c4", "c5", "c6"),
    ):
        print(f"keys={keys} -> H mixing {mixed_bucket_report(records, keys, 'H')}")
    for keys in (("score",), ("c3", "c4", "c5"), ("c3", "c4", "c5", "c6")):
        print(f"keys={keys} -> c6 mixing {mixed_bucket_report(records, keys, 'c6')}")


def main() -> None:
    validate_trace_formulas()
    benchmark()
    records = exhaustive_n6_records()
    print_exhaustive_patterns(records)
    print_information_tournament(records)
    print("\n=== Interpretation ===")
    print("1. Trace formulas are exact through length 5 because non-simple closed walks need two cycles.")
    print("2. Length 6 is the first correction: double triangles and intersecting triangle pairs.")
    print("3. c6 is still fast, but it is the first trace invariant needing cycle-placement data p33_meet.")
    print("4. At n=6, the low cycle vector (c3,c4,c5,c6) determines H; score+c5+c6 still mixes.")
    print("5. The next speed frontier is to test whether corrected trace vectors keep compressing H at n>=7.")


if __name__ == "__main__":
    main()
