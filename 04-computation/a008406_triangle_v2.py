#!/usr/bin/env python3
"""
A008406: Triangle T(n,k) = number of simple graphs on n nodes with k edges.

Formula:
  T(n,k) = [x^k] (1/n!) * sum_{sigma in S_n} GF_sigma(x)

where GF_sigma(x) = prod_{orbit O of sigma on 2-subsets} (1 + x^{|O|})

For cycle type lambda = [(r_a, m_a)]:

1. Within a single cycle of length r:
   - r odd:  (1+x^r)^{(r-1)/2}
   - r even: (1+x^r)^{(r-2)/2} * (1+x^{r/2})
   For m cycles: raise to m-th power.

2. Cross pairs between two cycles of lengths r, s:
   gcd(r,s) orbits, each of size lcm(r,s).
   GF: (1+x^{lcm(r,s)})^{gcd(r,s)}

Combined for compressed partition:
- Within-cycle for (r, m):
  r odd:  (1+x^r)^{m*(r-1)/2}
  r even: (1+x^r)^{m*(r-2)/2} * (1+x^{r/2})^m
- Same-part cross for (r, m):
  (1+x^r)^{C(m,2)*r}
- Cross-part for (r_a, m_a) and (r_b, m_b):
  (1+x^{lcm(r_a,r_b)})^{m_a*m_b*gcd(r_a,r_b)}

Author: opus-2026-03-08-S48
"""

from math import gcd, factorial
from fractions import Fraction
from time import perf_counter
import sys


def lcm(a, b):
    return a * b // gcd(a, b)


def binomial(n, k):
    if k < 0 or k > n:
        return 0
    if k == 0 or k == n:
        return 1
    k = min(k, n - k)
    result = 1
    for i in range(k):
        result = result * (n - i) // (i + 1)
    return result


def poly_mul_1_plus_xd_pow_e(poly, d, e, max_deg):
    """Multiply polynomial by (1+x^d)^e, truncated to max_deg.

    (1+x^d)^e = sum_{j=0}^{e} C(e,j) x^{j*d}
    """
    if e == 0 or d > max_deg:
        return poly

    # Build (1+x^d)^e as a polynomial
    factor = [0] * (max_deg + 1)
    for j in range(e + 1):
        deg = j * d
        if deg > max_deg:
            break
        factor[deg] = binomial(e, j)

    # Multiply
    result = [0] * (max_deg + 1)
    for i in range(max_deg + 1):
        if poly[i] == 0:
            continue
        for j in range(max_deg + 1 - i):
            if factor[j] == 0:
                continue
            result[i + j] += poly[i] * factor[j]

    return result


def compute_gf_for_partition(compressed, max_deg):
    """Compute GF(x) = prod of (1+x^{orbit_size}) terms for a partition."""
    poly = [0] * (max_deg + 1)
    poly[0] = 1

    d = len(compressed)

    for idx in range(d):
        r, m = compressed[idx]

        # Within-cycle terms: for each of m cycles of length r
        if r % 2 == 1:
            # r odd: (1+x^r)^{m*(r-1)/2}
            exp = m * (r - 1) // 2
            if exp > 0:
                poly = poly_mul_1_plus_xd_pow_e(poly, r, exp, max_deg)
        else:
            # r even: (1+x^r)^{m*(r-2)/2} * (1+x^{r/2})^m
            exp1 = m * (r - 2) // 2
            if exp1 > 0:
                poly = poly_mul_1_plus_xd_pow_e(poly, r, exp1, max_deg)
            poly = poly_mul_1_plus_xd_pow_e(poly, r // 2, m, max_deg)

        # Same-part cross-cycle terms: C(m,2) pairs, each pair gives r orbits of size r
        # (1+x^r)^{C(m,2)*r}
        cross_same = m * (m - 1) // 2 * r
        if cross_same > 0:
            poly = poly_mul_1_plus_xd_pow_e(poly, r, cross_same, max_deg)

    # Cross-part terms
    for i in range(d):
        ri, mi = compressed[i]
        for j in range(i + 1, d):
            rj, mj = compressed[j]
            g = gcd(ri, rj)
            l = lcm(ri, rj)
            exp = mi * mj * g
            if exp > 0 and l <= max_deg:
                poly = poly_mul_1_plus_xd_pow_e(poly, l, exp, max_deg)

    return poly


def a008406_row(n):
    """Compute row n of A008406."""
    max_edges = n * (n - 1) // 2
    fact_n = factorial(n)

    # Accumulate: (1/n!) * sum permcount * GF
    # = sum (1/z_lambda) * GF_lambda
    total_poly = [Fraction(0)] * (max_edges + 1)
    count = 0

    def enumerate_partitions(remaining, max_part, parts_so_far):
        nonlocal count
        if remaining == 0:
            count += 1
            gf = compute_gf_for_partition(parts_so_far, max_edges)
            z = 1
            for r, m in parts_so_far:
                z *= r ** m * factorial(m)
            coeff = Fraction(1, z)

            for k in range(max_edges + 1):
                if gf[k] != 0:
                    total_poly[k] += coeff * gf[k]
            return

        for p in range(min(remaining, max_part), 0, -1):
            max_m = remaining // p
            for m in range(1, max_m + 1):
                parts_so_far.append((p, m))
                enumerate_partitions(remaining - m * p, p - 1, parts_so_far)
                parts_so_far.pop()

    enumerate_partitions(n, n, [])

    row = []
    for k in range(max_edges + 1):
        val = total_poly[k]
        assert val.denominator == 1, f"n={n}, k={k}: non-integer {val}"
        row.append(int(val))

    return row, count


# Known A008406 rows
known = {
    1: [1],
    2: [1, 1],
    3: [1, 1, 1, 1],
    4: [1, 1, 2, 3, 2, 1, 1],
    5: [1, 1, 2, 4, 6, 6, 6, 4, 2, 1, 1],
    6: [1, 1, 2, 5, 9, 15, 21, 24, 24, 21, 15, 9, 5, 2, 1, 1],
}

if __name__ == "__main__":
    max_n = int(sys.argv[1]) if len(sys.argv) > 1 else 12

    print("A008406: T(n,k) = graphs on n nodes with k edges")
    print("=" * 60)

    for n in range(1, max_n + 1):
        t0 = perf_counter()
        row, nparts = a008406_row(n)
        dt = perf_counter() - t0

        expected = known.get(n)
        if expected:
            match = "OK" if row == expected else f"FAIL"
        else:
            match = ""

        # Print compact if row is long
        if len(row) > 20:
            print(f"  n={n}: [{', '.join(str(x) for x in row[:10])}, ..., {', '.join(str(x) for x in row[-5:])}] ({len(row)} terms) [{dt:.3f}s]  {match}")
        else:
            print(f"  n={n}: {row}  [{dt:.3f}s, {nparts} parts]  {match}")

        if dt > 120:
            print("  Stopping: too slow")
            break

    # B-file output
    if max_n >= 20:
        print("\n=== Generating b-file ===")
        with open("b008406.txt", "w") as f:
            idx = 0
            for n in range(1, max_n + 1):
                row, _ = a008406_row(n)
                for k in range(len(row)):
                    f.write(f"{idx} {row[k]}\n")
                    idx += 1
        print(f"Written b008406.txt with {idx} entries (rows 1..{max_n})")

    print("\nDone.")
