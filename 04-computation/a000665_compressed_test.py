#!/usr/bin/env python3
"""
Test compressed partition form of edges formula for A000665 (3-uniform hypergraphs).
Verify that the compressed form matches the expanded form.
"""

from math import gcd, factorial, ceil
from fractions import Fraction


def lcm(a, b):
    return a * b // gcd(a, b)

def lcm3(a, b, c):
    return lcm(lcm(a, b), c)


def edges_expanded(p):
    """Edges formula on expanded partition vector p (with repetitions)."""
    n = len(p)
    result = sum(ceil((p[i] - 1) * (p[i] - 2) / 6) for i in range(n))
    for i in range(1, n):
        for j in range(i):
            c, d = p[i], p[j]
            g = gcd(c, d)
            result += g * (c + d - 2 + ((c - d) // g) % 2) // 2
            for k in range(j):
                result += c * d * p[k] // lcm3(c, d, p[k])
    return result


def edges_compressed(parts_mults):
    """Edges formula on compressed form [(r, m), ...] sorted by r descending."""
    d = len(parts_mults)

    # Single-cycle terms
    single = 0
    for r, m in parts_mults:
        single += m * ceil((r - 1) * (r - 2) / 6)

    # Pair terms
    pair = 0
    for i in range(d):
        r, mr = parts_mults[i]
        # Same part pairs: C(m, 2) * r * (r - 1)
        pair += mr * (mr - 1) // 2 * r * (r - 1)
        # Cross-part pairs
        for j in range(i + 1, d):
            s, ms = parts_mults[j]
            g = gcd(r, s)
            pair_val = g * (r + s - 2 + ((r - s) // g) % 2) // 2
            pair += mr * ms * pair_val

    # Triple terms
    triple = 0
    for i in range(d):
        r, mr = parts_mults[i]
        # All three from same part: C(m, 3) * r^2
        triple += mr * (mr - 1) * (mr - 2) // 6 * r * r
        # Two from part i, one from part j (i != j)
        for j in range(d):
            if j == i:
                continue
            s, ms = parts_mults[j]
            # C(mr, 2) * ms * r * gcd(r, s)
            triple += mr * (mr - 1) // 2 * ms * r * gcd(r, s)
        # One from each of three distinct parts
        for j in range(i + 1, d):
            s, ms = parts_mults[j]
            for k in range(j + 1, d):
                t, mt = parts_mults[k]
                triple += mr * ms * mt * r * s * t // lcm3(r, s, t)

    return single + pair + triple


def compress(p):
    """Convert expanded partition to compressed form."""
    from collections import Counter
    c = Counter(p)
    return sorted(c.items(), reverse=True)


def test_edges():
    """Test compressed vs expanded edges for all partitions of n."""
    from itertools import combinations_with_replacement

    for n in range(3, 16):
        # Generate all partitions of n
        parts = []
        def gen(remaining, max_part, current):
            if remaining == 0:
                parts.append(tuple(current))
                return
            for p in range(min(remaining, max_part), 0, -1):
                current.append(p)
                gen(remaining - p, p, current)
                current.pop()
        gen(n, n, [])

        all_ok = True
        for p in parts:
            e_exp = edges_expanded(list(p))
            e_comp = edges_compressed(compress(list(p)))
            if e_exp != e_comp:
                print(f"MISMATCH at n={n}, p={p}: expanded={e_exp}, compressed={e_comp}")
                all_ok = False

        status = "OK" if all_ok else "FAIL"
        print(f"n={n:2d}: {len(parts):6d} partitions  {status}")


def test_full():
    """Verify A000665 using compressed edges."""
    known = {0: 1, 1: 1, 2: 1, 3: 2, 4: 5, 5: 34, 6: 2136, 7: 7013320, 8: 1788782616656}

    for n in range(3, 13):
        total = Fraction(0)
        fact_n = factorial(n)

        def enum(remaining, max_part, parts_so_far):
            nonlocal total
            if remaining == 0:
                e = edges_compressed(parts_so_far)
                pc = 1
                # Compute permcount from compressed form
                s = sum(r * m for r, m in parts_so_far)
                z = 1
                for r, m in parts_so_far:
                    z *= r ** m * factorial(m)
                pc = factorial(s) // z
                total += Fraction(pc * (1 << e), fact_n)
                return
            for pi in range(max_part, 0, -1):
                if pi > remaining:
                    continue
                max_m = remaining // pi
                for m in range(1, max_m + 1):
                    parts_so_far.append((pi, m))
                    enum(remaining - m * pi, pi - 1, parts_so_far)
                    parts_so_far.pop()

        enum(n, n, [])
        result = int(total)
        exp = known.get(n, "?")
        match = "OK" if result == exp else f"FAIL (expected {exp})"
        print(f"a({n}) = {result}  {match}")


if __name__ == "__main__":
    print("Testing compressed vs expanded edges formula...")
    test_edges()
    print()
    print("Testing A000665 with compressed form...")
    test_full()
