#!/usr/bin/env python3
"""Merged tiling bucket constraints.

kind-pasteur-2026-05-29-S5

This script formalizes the "merged tiling bucket" viewpoint.

Let Q_m be the fixed-base-path tiling cube with m = C(n-1, 2), and let

    pi : Q_m -> G_n/Z_2

send a tiling to its merged tournament isomorphism class.  The fibers of pi
are the merged buckets.  For every Hamming layer d, flipping any d-set of
tiles gives a symmetric integer transport matrix W_d between buckets.

The all-n theorem checked here in small cases:

1. Bucket size parity:
   - self-complementary bucket size is odd;
   - non-self-complementary bucket size is 2 mod 4.

2. Hamming-layer parity:
   - row_sum_i(W_d) = bucket_i * C(m, d);
   - W_d is symmetric;
   - diagonal entries are even;
   - cross-outflow parity from bucket i is
         bucket_i * C(m, d) mod 2,
     hence it is odd exactly when bucket i is SC and C(m, d) is odd.

The proof is independent of enumeration; enumeration here catches mistakes
in the old S202 narration and gives exact n=3..6 tables.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from itertools import combinations, permutations
from math import comb


def tile_pairs(n: int) -> list[tuple[int, int]]:
    """Explorer order: lo increases, hi decreases."""
    return [(hi, lo) for lo in range(n - 2) for hi in range(n - 1, lo + 1, -1)]


def tiling_to_tournament(bits: int, n: int, tiles: list[tuple[int, int]]) -> list[list[int]]:
    A = [[0] * n for _ in range(n)]

    # Base Hamiltonian path: n-1 -> n-2 -> ... -> 0.
    for lo in range(n - 1):
        A[lo + 1][lo] = 1

    for k, (hi, lo) in enumerate(tiles):
        if bits & (1 << k):
            A[lo][hi] = 1
        else:
            A[hi][lo] = 1
    return A


def complement_tournament(A: list[list[int]]) -> list[list[int]]:
    n = len(A)
    return [[0 if i == j else 1 - A[i][j] for j in range(n)] for i in range(n)]


def tournament_canon(A: list[list[int]], perms: list[tuple[int, ...]]) -> tuple[int, ...]:
    n = len(A)
    best = None
    for p in perms:
        key = tuple(A[p[i]][p[j]] for i in range(n) for j in range(i + 1, n))
        if best is None or key < best:
            best = key
    assert best is not None
    return best


def merged_canon_and_sc(
    A: list[list[int]], perms: list[tuple[int, ...]]
) -> tuple[tuple[int, ...], bool]:
    c = tournament_canon(A, perms)
    cc = tournament_canon(complement_tournament(A), perms)
    return min(c, cc), c == cc


def masks_by_weight(m: int) -> dict[int, list[int]]:
    by_d: dict[int, list[int]] = defaultdict(list)
    for d in range(1, m + 1):
        for combo in combinations(range(m), d):
            mask = 0
            for k in combo:
                mask |= 1 << k
            by_d[d].append(mask)
    return by_d


def build_bucket_data(n: int):
    tiles = tile_pairs(n)
    m = len(tiles)
    total = 1 << m
    perms = list(permutations(range(n)))

    bucket_id: dict[tuple[int, ...], int] = {}
    bucket_sc: dict[int, bool] = {}
    tiling_bucket: dict[int, int] = {}

    for bits in range(total):
        A = tiling_to_tournament(bits, n, tiles)
        canon, is_sc = merged_canon_and_sc(A, perms)
        if canon not in bucket_id:
            bid = len(bucket_id)
            bucket_id[canon] = bid
            bucket_sc[bid] = is_sc
        bid = bucket_id[canon]
        if bucket_sc[bid] != is_sc:
            raise AssertionError(f"mixed SC status inside merged bucket {bid}")
        tiling_bucket[bits] = bid

    bucket_sizes = Counter(tiling_bucket.values())
    return m, total, tiling_bucket, bucket_sizes, bucket_sc


def analyze_layer(
    m: int,
    total: int,
    tiling_bucket: dict[int, int],
    bucket_sizes: Counter[int],
    bucket_sc: dict[int, bool],
    d: int,
    masks: list[int],
):
    nb = len(bucket_sizes)
    W = [[0] * nb for _ in range(nb)]
    for bits in range(total):
        i = tiling_bucket[bits]
        for mask in masks:
            j = tiling_bucket[bits ^ mask]
            W[i][j] += 1

    binom_md = comb(m, d)
    row_viol = 0
    sym_viol = 0
    diag_odd = 0
    parity_viol = 0
    odd_cross_rows = 0
    directed_cross = 0
    directed_self = 0
    type_cross = Counter()

    for i in range(nb):
        row_sum = sum(W[i])
        expected_row = bucket_sizes[i] * binom_md
        if row_sum != expected_row:
            row_viol += 1

        if W[i][i] % 2:
            diag_odd += 1

        cross_out = row_sum - W[i][i]
        directed_cross += cross_out
        directed_self += W[i][i]
        if cross_out % 2:
            odd_cross_rows += 1

        expected_parity = (bucket_sizes[i] * binom_md) % 2
        if cross_out % 2 != expected_parity:
            parity_viol += 1

        for j in range(i + 1, nb):
            if W[i][j] != W[j][i]:
                sym_viol += 1
            if W[i][j]:
                if bucket_sc[i] and bucket_sc[j]:
                    etype = "SC-SC"
                elif bucket_sc[i] or bucket_sc[j]:
                    etype = "SC-NS"
                else:
                    etype = "NS-NS"
                type_cross[etype] += W[i][j]

    expected_odd_rows = sum(
        1 for i in range(nb) if bucket_sc[i] and (binom_md % 2 == 1)
    )

    return {
        "d": d,
        "binom": binom_md,
        "binom_odd": binom_md % 2,
        "row_viol": row_viol,
        "sym_viol": sym_viol,
        "diag_odd": diag_odd,
        "parity_viol": parity_viol,
        "odd_cross_rows": odd_cross_rows,
        "expected_odd_rows": expected_odd_rows,
        "cross_lines": directed_cross // 2,
        "self_lines": directed_self // 2,
        "type_cross": dict(type_cross),
        "total_lines": total * binom_md // 2,
        "parity_lower_bound_cross_lines": expected_odd_rows // 2,
    }


def analyze_n(n: int):
    m, total, tiling_bucket, bucket_sizes, bucket_sc = build_bucket_data(n)
    nb = len(bucket_sizes)
    sc_count = sum(1 for v in bucket_sc.values() if v)
    ns_count = nb - sc_count

    bucket_parity_violations = []
    bucket_hist = Counter()
    for bid in range(nb):
        size = bucket_sizes[bid]
        is_sc = bucket_sc[bid]
        bucket_hist[("SC" if is_sc else "NS", size)] += 1
        if is_sc and size % 2 != 1:
            bucket_parity_violations.append((bid, "SC", size))
        if not is_sc and size % 4 != 2:
            bucket_parity_violations.append((bid, "NS", size))

    layer_masks = masks_by_weight(m)
    layers = [
        analyze_layer(m, total, tiling_bucket, bucket_sizes, bucket_sc, d, layer_masks[d])
        for d in range(1, m + 1)
    ]

    return {
        "n": n,
        "m": m,
        "total": total,
        "nb": nb,
        "sc_count": sc_count,
        "ns_count": ns_count,
        "bucket_sum": sum(bucket_sizes.values()),
        "bucket_hist": bucket_hist,
        "bucket_parity_violations": bucket_parity_violations,
        "layers": layers,
    }


def print_report(results):
    print("Merged tiling bucket constraints (kind-pasteur-2026-05-29-S5)")
    print()
    print("Theorem schema:")
    print("  bucket_size(M) is odd for SC buckets and 2 mod 4 for NS buckets.")
    print("  For layer d, W_d has row sums bucket_size(M)*C(m,d),")
    print("  symmetric entries, even diagonal, and cross-row parity")
    print("  bucket_size(M)*C(m,d) mod 2.")

    for r in results:
        print()
        print("=" * 88)
        print(
            f"n={r['n']} m={r['m']} tilings={r['total']} "
            f"merged_buckets={r['nb']} SC={r['sc_count']} NS={r['ns_count']}"
        )
        print(f"bucket sum = {r['bucket_sum']} (expected {r['total']})")
        print(f"bucket parity violations = {len(r['bucket_parity_violations'])}")
        print("bucket histogram (type,size -> count):")
        for key, count in sorted(r["bucket_hist"].items(), key=lambda kv: (kv[0][0], kv[0][1])):
            print(f"  {key}: {count}")

        print()
        print(
            "layer d | C(m,d)%2 | odd rows actual/expected | "
            "cross lines | parity LB | violations(row,sym,diag,parity)"
        )
        for layer in r["layers"]:
            print(
                f"{layer['d']:7d} | {layer['binom_odd']:8d} | "
                f"{layer['odd_cross_rows']:4d}/{layer['expected_odd_rows']:<4d} | "
                f"{layer['cross_lines']:11d} | "
                f"{layer['parity_lower_bound_cross_lines']:9d} | "
                f"({layer['row_viol']},{layer['sym_viol']},"
                f"{layer['diag_odd']},{layer['parity_viol']})"
            )
            tc = layer["type_cross"]
            print(
                "          type cross-lines: "
                f"SC-SC={tc.get('SC-SC', 0)}, "
                f"SC-NS={tc.get('SC-NS', 0)}, "
                f"NS-NS={tc.get('NS-NS', 0)}"
            )

        all_clean = (
            not r["bucket_parity_violations"]
            and all(
                layer["row_viol"] == 0
                and layer["sym_viol"] == 0
                and layer["diag_odd"] == 0
                and layer["parity_viol"] == 0
                and layer["odd_cross_rows"] == layer["expected_odd_rows"]
                for layer in r["layers"]
            )
        )
        print(f"all constraints clean: {all_clean}")


def main() -> None:
    results = [analyze_n(n) for n in range(3, 7)]
    print_report(results)


if __name__ == "__main__":
    main()
