"""
paley_betti_pattern.py — Analyze patterns in Paley tournament Betti numbers

Known: T_3: (1, 1, 0), T_7: (1, 0, 0, 0, 6, 0, 0)
Pattern: nonzero Betti at degree 0 and degree p-3?

T_3: b₀=1, b₁=1 (deg p-2 = 1)
T_7: b₀=1, b₄=6 (deg p-3 = 4)

Conjectured: for Paley prime p = 3 mod 4, b_k = 0 for 0 < k < p-3,
and b_{p-3} = some function of p.

Also check: is there a duality b_k = b_{n-2-k} for VT tournaments?

Author: kind-pasteur-2026-03-10-S50
"""
import sys
import time
import numpy as np
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

from tournament_utils import (
    full_chain_complex_modp, RANK_PRIME
)


def paley_adj(p):
    """Construct Paley tournament for prime p = 3 mod 4."""
    qr = set()
    for i in range(1, p):
        qr.add((i * i) % p)

    # Verify p = 3 mod 4 (so -1 is NOT a QR)
    if p % 4 != 3:
        print(f"  WARNING: p={p} = {p%4} mod 4, NOT a valid Paley prime!")
        return None

    A = [[0] * p for _ in range(p)]
    for i in range(p):
        for j in range(p):
            if i != j:
                if (j - i) % p in qr:
                    A[i][j] = 1
    return A


def betti_duality_check(bettis, n):
    """Check if b_k = b_{n-2-k} (Poincaré-type duality)."""
    max_deg = max(bettis.keys()) if bettis else 0
    dual_holds = True
    for k in bettis:
        dual_k = n - 2 - k
        if dual_k >= 0 and dual_k <= max_deg:
            if bettis.get(k, 0) != bettis.get(dual_k, 0):
                dual_holds = False
    return dual_holds


def main():
    print("="*70)
    print("PALEY TOURNAMENT BETTI PATTERNS")
    print("="*70)

    # Valid Paley primes (= 3 mod 4): 3, 7, 11, 19, 23, 31, 43, 47, ...
    paley_primes = [3, 7, 11, 19, 23]

    results = {}
    for p in paley_primes:
        print(f"\n--- T_{p} (Paley, p = 3 mod 4) ---")
        A = paley_adj(p)
        if A is None:
            continue

        max_deg = min(p - 1, 6)  # limit computation
        if p >= 19:
            max_deg = 4  # even more limited for large p
        if p >= 23:
            max_deg = 3

        t0 = time.time()
        cc = full_chain_complex_modp(A, p, max_deg)
        elapsed = time.time() - t0

        bettis = cc['bettis']
        chi = sum((-1)**k * v for k, v in bettis.items())
        print(f"  Betti (up to deg {max_deg}): {dict(sorted(bettis.items()))}")
        print(f"  Euler char: {chi}")
        print(f"  Time: {elapsed:.2f}s")

        # Duality check
        dual = betti_duality_check(bettis, p)
        print(f"  Duality b_k = b_{{p-2-k}}: {dual}")

        # Pattern: where is the first nonzero b_k for k > 0?
        nonzero = [(k, v) for k, v in sorted(bettis.items()) if k > 0 and v > 0]
        if nonzero:
            print(f"  Nonzero b_k (k>0): {nonzero}")
            print(f"  First nonzero at degree {nonzero[0][0]} (= p-{p-nonzero[0][0]})")
        else:
            print(f"  ALL b_k = 0 for k > 0 (within computed range)")

        results[p] = bettis

    # Summary table
    print(f"\n{'='*70}")
    print("SUMMARY:")
    for p, bettis in sorted(results.items()):
        nonzero = [(k, v) for k, v in sorted(bettis.items()) if k > 0 and v > 0]
        nz_str = ", ".join(f"b_{k}={v}" for k, v in nonzero) if nonzero else "all zero"
        print(f"  T_{p}: b₀=1, {nz_str}  chi={sum((-1)**k * v for k, v in bettis.items())}")


if __name__ == '__main__':
    main()
    print("\nDONE.")
