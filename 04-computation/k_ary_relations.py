#!/usr/bin/env python3
"""
Compute k-ary relations on n unlabeled nodes.

a(n) = (1/n!) * sum_{sigma in S_n} 2^{orbits of sigma on [n]^k}

For a permutation sigma with cycle type [(r_i, m_i)]:
  - L = lcm of all r_i
  - f(t) = number of fixed points of sigma^t = sum of r_i where r_i | t
  - orbits on [n]^k = (1/L) * sum_{t=0}^{L-1} f(t)^k

Known OEIS sequences:
  k=2: A000595 (binary relations)
  k=3: A000662
  k=4: A001377
  k=5: A051241
  k>=6: not yet in OEIS

Author: opus-2026-03-09-S50
"""

from math import gcd, factorial
from fractions import Fraction
from functools import reduce
import sys
import time


def lcm(a, b):
    return a * b // gcd(a, b)


def k_ary_relations(n, k):
    """Number of k-ary relations on n unlabeled nodes."""
    if n == 0:
        return 1

    result = Fraction(0)

    def enum_parts(remaining, max_part, parts):
        nonlocal result
        if remaining == 0:
            # Expand cycle lengths
            cycles = []
            for r, m in parts:
                cycles.extend([r] * m)

            # z_lambda = prod r^m * m!
            z = 1
            for r, m in parts:
                z *= r**m * factorial(m)

            # L = lcm of cycle lengths
            L = reduce(lcm, [r for r, m in parts])

            # Orbits on [n]^k via Burnside
            total_fixed = 0
            for t in range(L):
                if t == 0:
                    ft = n  # identity fixes all n elements
                else:
                    ft = sum(c for c in cycles if t % c == 0)
                total_fixed += ft ** k

            num_orbits = Fraction(total_fixed, L)
            assert num_orbits.denominator == 1

            result += Fraction(1, z) * (2 ** int(num_orbits))
            return

        for p in range(min(remaining, max_part), 0, -1):
            for m in range(1, remaining // p + 1):
                parts.append((p, m))
                enum_parts(remaining - m * p, p - 1, parts)
                parts.pop()

    enum_parts(n, n, [])
    assert result.denominator == 1
    return int(result)


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python3 k_ary_relations.py <k> [max_n] [output_file]")
        print("  k: arity (2=binary, 3=ternary, etc.)")
        print("  max_n: compute a(0)..a(max_n) (default 30)")
        print("  output_file: b-file name (default: b_kary_k<k>.txt)")
        sys.exit(1)

    k = int(sys.argv[1])
    max_n = int(sys.argv[2]) if len(sys.argv) > 2 else 30
    outfile = sys.argv[3] if len(sys.argv) > 3 else f"b_kary_k{k}.txt"

    # Map to OEIS sequence IDs
    oeis_map = {2: "A000595", 3: "A000662", 4: "A001377", 5: "A051241"}
    seq_name = oeis_map.get(k, f"k={k}-ary relations")

    print(f"Computing {seq_name}: {k}-ary relations on n unlabeled nodes")
    print(f"Output: {outfile}, n=0..{max_n}")
    print()

    with open(outfile, 'w') as f:
        for n in range(max_n + 1):
            t0 = time.time()
            val = k_ary_relations(n, k)
            dt = time.time() - t0
            f.write(f"{n} {val}\n")
            f.flush()

            s = str(val)
            display = s[:40] + '...' if len(s) > 40 else s
            print(f"a({n}) = {display} ({len(s)} digits)  [{dt:.3f}s]")

            if dt > 300:
                print("Stopping: too slow")
                break

    print(f"\nWritten {outfile}")
