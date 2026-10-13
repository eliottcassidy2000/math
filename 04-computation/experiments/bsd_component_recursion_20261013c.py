"""Exact SCC factorization of the tournament Monsky obstruction.

Run: python 04-computation/experiments/bsd_component_recursion_20261013c.py
No arithmetic claim beyond the cited Monsky family in the paired note.
"""
from collections import Counter
from itertools import product, permutations
from math import comb, prod, isqrt
import json
from bsd_selmer_growth_20261013a import tournament_a, phi3_rows, rank, monsky, factors

CHECKS = 0


def check(value, message):
    global CHECKS
    CHECKS += 1
    if not value:
        raise RuntimeError(message)


def components(a):
    """Independent reachability closure; source-to-sink SCC order."""
    k = len(a)
    reach = [(1 << i) | sum(a[i][j] << j for j in range(k) if i != j)
             for i in range(k)]
    for t in range(k):
        for i in range(k):
            if (reach[i] >> t) & 1:
                reach[i] |= reach[t]
    remaining = set(range(k))
    blocks = []
    while remaining:
        i = min(remaining)
        block = sorted(j for j in remaining if (reach[i] >> j) & (reach[j] >> i) & 1)
        blocks.append(block)
        remaining.difference_update(block)
    return sorted(blocks, key=lambda b: -sum(a[b[0]][v] for other in blocks for v in other if v not in b))


def induced(a, vertices):
    b = [[a[i][j] if i != j else 0 for j in vertices] for i in vertices]
    for i in range(len(b)):
        b[i][i] = sum(b[i]) % 2
    return b


def join(a, b):
    n, m = len(a), len(b)
    c = [[0] * (n + m) for _ in range(n + m)]
    for i in range(n):
        for j in range(n):
            c[i][j] = a[i][j]
        c[i][i] ^= m % 2
        for j in range(n, n + m):
            c[i][j] = 1
    for i in range(m):
        for j in range(m):
            c[n + i][n + j] = b[i][j]
    return c


def dim(a):
    return len(a) - rank(phi3_rows(a), len(a))


def main():
    all_counts = {0: Counter({0: 1})}
    strong_counts = {0: Counter()}
    strong_zero = None
    for n in range(1, 7):
        all_counts[n], strong_counts[n] = Counter(), Counter()
        for code in range(1 << comb(n, 2)):
            a = tournament_a(n, code)
            blocks = components(a)
            check(all(a[u][v] == 1 for i, left in enumerate(blocks)
                      for right in blocks[i + 1:] for u in left for v in right),
                  "literal source-to-sink component order")
            ph = phi3_rows(a)
            expected = [0] * n
            subtotal = 0
            for block in blocks:
                small = induced(a, block)
                sp = phi3_rows(small)
                subtotal += dim(small)
                for i, v in enumerate(block):
                    expected[v] = sum(((sp[i] >> j) & 1) << w for j, w in enumerate(block))
            check(ph == expected, "entire Phi3 matrix must split, not just nullity")
            d = dim(a)
            check(d == subtotal, "SCC additivity")
            all_counts[n][d] += 1
            if len(blocks) == 1:
                strong_counts[n][d] += 1
                if n > 1 and d == 0 and strong_zero is None:
                    strong_zero = {"order": n, "code": code, "matrix": a}
        # Unique source SCC gives a labelled convolution identity.
        reconstructed = Counter()
        for k in range(1, n + 1):
            for d, u in strong_counts[k].items():
                for e, v in all_counts[n - k].items():
                    reconstructed[d + e] += comb(n, k) * u * v
        check(reconstructed == all_counts[n], "labelled SCC convolution")

    # Exhaustive independent two-block matrix identity, through sizes 4+3.
    joins = 0
    for n in range(1, 5):
        for m in range(1, 4):
            for a_code, b_code in product(range(1 << comb(n, 2)), range(1 << comb(m, 2))):
                a, b = tournament_a(n, a_code), tournament_a(m, b_code)
                expected = phi3_rows(a) + [v << n for v in phi3_rows(b)]
                check(phi3_rows(join(a, b)) == expected, "ordered join matrix identity")
                joins += 1

    # Arithmetic labels are authentic, distinct, all 3 mod 8.
    witness = []
    for ps in ([3, 11, 19], [3, 11, 59]):
        n = prod(ps)
        got, a, rows = monsky(n)
        check(got == ps and all(p % 8 == 3 for p in got), "arithmetic sidecar")
        d = dim(a)
        check(d == 2 * len(ps) - rank(rows, 2 * len(ps)), "full Monsky kernel agrees")
        witness.append({"n": n, "primes": ps, "components": components(a), "pure_selmer_dimension": d})
    # Find an actual fourth prime destroying the cyclic triple's obstruction.
    extension = None
    for p in range(27, 2000, 8):
        if any(p % d == 0 for d in range(2, isqrt(p) + 1)) or p in (3, 11, 19):
            continue
        n = 627 * p
        ps, a, rows = monsky(n)
        if dim(a) == 0:
            extension = {"new_prime": p, "new_curve_n": n, "primes": ps,
                         "pure_selmer_dimension": 2 * len(ps) - rank(rows, 2 * len(ps))}
            break
    check(extension is not None and extension["pure_selmer_dimension"] == 0,
          "mixed extension defeats monotone obstruction inference")
    check(strong_zero is not None, "strong cycle does not force Phi3 obstruction")
    check(components(tournament_a(2, 0)) == [[1], [0]], "reverse-label component-order hostile")
    strong_a = strong_zero["matrix"]
    paths = [p for p in permutations(range(4)) if all(strong_a[p[i]][p[i + 1]] for i in range(3))]
    signed = sum((-1) ** sum(p[i] > p[j] for i in range(4) for j in range(i + 1, 4)) for p in paths)
    check(len(paths) == 5 and abs(signed) == 1, "strong four Hamiltonian correction")
    strong_zero["hamiltonian_paths"] = [list(p) for p in paths]
    strong_zero["signed_path_sum"] = signed
    check(sum(sum(all_counts[n].values()) for n in range(1, 6)) == 1099,
          "all tournaments on one through five labelled vertices")
    report = {"status": "PROVED binary matrix identity; FINITE-EXACT controls; CITED Monsky interface",
              "universe": "all labelled tournaments on 1..6 vertices; independent joins 1..4 by 1..3",
              "all_histograms": {n: dict(sorted(c.items())) for n, c in all_counts.items()},
              "strong_histograms": {n: dict(sorted(c.items())) for n, c in strong_counts.items()},
              "first_strong_invertible": strong_zero,
              "arithmetic_controls": witness, "different_curve_hostile": extension,
              "independent_joins": joins, "checks": CHECKS}
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
