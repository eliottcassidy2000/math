#!/usr/bin/env python3
"""
Compute derived OEIS sequences via Euler/inverse Euler transforms.

Given b-files for base sequences, compute:
  A001349 = InvEuler(A000088)  (connected graphs from all graphs)
  A003085 = InvEuler(A000273)  (weakly connected digraphs from all digraphs)
  A003084 from A000273 (via logarithmic derivative)
  A051337 = InvEuler-like(A000568) (strongly connected tournaments)

The Euler transform:
  If b(n) = connected count, then B(x) = sum b(n) x^n
  and A(x) = prod_{k>=1} 1/(1-x^k)^{b(k)} = sum a(n) x^n
  gives a(n) = total count (connected + disconnected)

Inverse: given a(n), find b(n).
  Let A(x) = sum a(n) x^n (a(0)=1).
  Then log A(x) = sum_{k>=1} b(k) * sum_{m>=1} x^{km}/m
  So n * [x^n] log A(x) = sum_{d|n} d * b(d) / (n/d)
  ... this is related to the Euler transform inversion.

Actually, the standard inversion:
  a(0) = 1 (or the base count for 0 nodes)
  c(n) = sum_{d|n} d * b(d)  (Dirichlet series of b weighted by d)
  n * a(n) = sum_{k=1}^{n} c(k) * a(n-k)

Inverse:
  c(n) = n * a(n) - sum_{k=1}^{n-1} c(k) * a(n-k)
  b(n) = (1/n) * sum_{d|n} mu(n/d) * c(d)

Author: opus-2026-03-08-S48
"""

import os
import sys
from collections import defaultdict

def load_bfile(path):
    vals = {}
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) >= 2:
                n, v = int(parts[0]), int(parts[1])
                vals[n] = v
    return vals

def mobius(n):
    """Möbius function μ(n)."""
    if n == 1:
        return 1
    factors = {}
    m = n
    for p in range(2, int(m**0.5) + 1):
        while m % p == 0:
            factors[p] = factors.get(p, 0) + 1
            m //= p
    if m > 1:
        factors[m] = factors.get(m, 0) + 1
    for e in factors.values():
        if e > 1:
            return 0
    return (-1) ** len(factors)

def inverse_euler_transform(a, max_n):
    """Given a(0), a(1), ..., a(max_n), compute b(1), b(2), ..., b(max_n).

    a(n) = total count, b(n) = connected count.
    a(0) must be 1.

    Uses the standard inversion:
      c(n) = n*a(n) - sum_{k=1}^{n-1} c(k)*a(n-k)
      b(n) = (1/n) * sum_{d|n} mu(n/d) * c(d)
    """
    assert a.get(0, 1) == 1, f"a(0) must be 1, got {a.get(0)}"

    c = {}
    b = {}

    for n in range(1, max_n + 1):
        if n not in a:
            break
        c[n] = n * a[n] - sum(c.get(k, 0) * a.get(n - k, 0) for k in range(1, n))

        # b(n) = (1/n) * sum_{d|n} mu(n/d) * c(d)
        s = sum(mobius(n // d) * c.get(d, 0) for d in range(1, n + 1) if n % d == 0)
        assert s % n == 0, f"b({n}) not integer: {s}/{n}"
        b[n] = s // n

    return b

def ogf_inverse(T, max_n):
    """Compute strongly connected count from tournament count via OGF inversion.

    If B(x) = sum T(n) x^n, then SC(n) comes from 1 - 1/B(x).
    This is: b[0] = 1, b[n] = -sum_{k=1}^n T[k]*b[n-k] for n >= 1
    SC(n) = -b[n] for n >= 1.
    """
    from fractions import Fraction
    b = [Fraction(0)] * (max_n + 1)
    b[0] = Fraction(1)
    for n in range(1, max_n + 1):
        s = Fraction(0)
        for k in range(1, n + 1):
            if k in T:
                s += T[k] * b[n - k]
        b[n] = -s

    sc = {}
    for n in range(1, max_n + 1):
        sc[n] = int(-b[n])
    return sc


def save_bfile(vals, path, start_n=None):
    if start_n is None:
        start_n = min(vals.keys())
    with open(path, 'w') as f:
        for n in sorted(vals.keys()):
            if n >= start_n:
                f.write(f"{n} {vals[n]}\n")
    print(f"  Saved {len([n for n in vals if n >= start_n])} values to {path}")


if __name__ == "__main__":
    base_dir = os.path.dirname(os.path.abspath(__file__))

    # =============================================
    # A001349 = InvEuler(A000088) — connected graphs
    # =============================================
    bfile = os.path.join(base_dir, 'b000088.txt')
    if os.path.exists(bfile):
        print("Computing A001349 (connected graphs) from A000088...")
        a088 = load_bfile(bfile)
        max_n = max(a088.keys())
        b = inverse_euler_transform(a088, max_n)

        # Validate against known values
        known = {1: 1, 2: 1, 3: 2, 4: 6, 5: 21, 6: 112, 7: 853, 8: 11117,
                 9: 261080, 10: 11716571}
        print("  Validation:")
        for n, exp in sorted(known.items()):
            got = b.get(n, "?")
            ok = "OK" if got == exp else f"FAIL (got {got})"
            print(f"    b({n}) = {exp}  {ok}")

        outpath = os.path.join(base_dir, 'b001349.txt')
        save_bfile({0: 1, **b}, outpath, start_n=0)
        print(f"  A001349: {max_n} terms (OEIS has 76)")

    # =============================================
    # A003085 = InvEuler(A000273) — weakly connected digraphs
    # =============================================
    bfile = os.path.join(base_dir, 'b000273.txt')
    if os.path.exists(bfile):
        print("\nComputing A003085 (weakly connected digraphs) from A000273...")
        a273 = load_bfile(bfile)
        max_n = max(a273.keys())
        b = inverse_euler_transform(a273, max_n)

        known = {1: 1, 2: 2, 3: 13, 4: 199, 5: 9364, 6: 1530843,
                 7: 880471142}
        print("  Validation:")
        for n, exp in sorted(known.items()):
            got = b.get(n, "?")
            ok = "OK" if got == exp else f"FAIL (got {got})"
            print(f"    b({n}) = {exp}  {ok}")

        outpath = os.path.join(base_dir, 'b003085.txt')
        save_bfile(b, outpath, start_n=1)
        print(f"  A003085: {max_n} terms (OEIS has 64)")

    # =============================================
    # A051337 from A000568 — strongly connected tournaments
    # =============================================
    bfile = os.path.join(base_dir, 'b000568.txt')
    if os.path.exists(bfile):
        print("\nComputing A051337 (strongly connected tournaments) from A000568...")
        a568 = load_bfile(bfile)
        max_n = max(a568.keys())
        sc = ogf_inverse(a568, max_n)

        known = {1: 1, 2: 0, 3: 1, 4: 1, 5: 6, 6: 35, 7: 353,
                 8: 6008, 9: 178133, 10: 9355949}
        print("  Validation:")
        for n, exp in sorted(known.items()):
            got = sc.get(n, "?")
            ok = "OK" if got == exp else f"FAIL (got {got})"
            print(f"    sc({n}) = {exp}  {ok}")

        outpath = os.path.join(base_dir, 'b051337.txt')
        save_bfile(sc, outpath, start_n=1)
        print(f"  A051337: {max_n} terms (OEIS has 50)")

    # =============================================
    # Connected binary relations from A000595
    # Actually: A000595 = labeled binary relations on n points
    # Its inverse Euler transform gives connected binary relations
    # Let's check if that's a known sequence
    # =============================================
    bfile = os.path.join(base_dir, 'b000595.txt')
    if os.path.exists(bfile):
        print("\nComputing InvEuler(A000595) — connected binary relations...")
        a595 = load_bfile(bfile)
        max_n = max(a595.keys())
        b = inverse_euler_transform(a595, max_n)

        print("  First values:")
        for n in range(1, min(max_n + 1, 11)):
            print(f"    b({n}) = {b.get(n, '?')}")

        # This might be A001173 or similar
        outpath = os.path.join(base_dir, 'b_inv_euler_a000595.txt')
        save_bfile(b, outpath, start_n=1)

    # =============================================
    # Connected symmetric relations from A000666
    # =============================================
    bfile = os.path.join(base_dir, 'b000666.txt')
    if os.path.exists(bfile):
        print("\nComputing InvEuler(A000666) — connected symmetric relations...")
        a666 = load_bfile(bfile)
        max_n = max(a666.keys())
        b = inverse_euler_transform(a666, max_n)

        print("  First values:")
        for n in range(1, min(max_n + 1, 11)):
            print(f"    b({n}) = {b.get(n, '?')}")

        outpath = os.path.join(base_dir, 'b_inv_euler_a000666.txt')
        save_bfile(b, outpath, start_n=1)

    print("\nAll transforms complete.")
