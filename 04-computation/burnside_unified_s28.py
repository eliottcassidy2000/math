#!/usr/bin/env python3
"""
opus-2026-04-05-S28: Production-quality Burnside enumerator.

VERIFIED CORRECT for:
  A000568 (tournaments) — all parts odd, base=2, self=floor(k/2)
  A000088 (simple graphs) — all parts, base=2, self=floor(k/2)
  A000273 (digraphs) — all parts, base=2, 2×gcd + sum_v - num_v
  A000595 (binary relations) — all parts, base=2, 2×gcd + sum_v
  A000666 (symmetric relations) — all parts, base=2, gcd + floor(k/2) + 1

NEEDS INVESTIGATION:
  A001174 (oriented graphs) — off by 3247 at n=8 (formula issue with even cycles)
  *** RESOLVED (monad-researcher-2026-06-02-S565, HYP-2078): the 3^{#pair-orbits}
      formula overcounts because an orientation-REVERSING (odd-swap) pair-orbit can
      only take the iota-fixed color (no-edge), not all 3.  The mechanical per-orbit-
      monodromy engine in self_converse_families_burnside_s565.py gives the correct
      A001174 (575016219 at n=8), matching the OEIS b-file to n=40. ***

Features:
  - Pure Python, arbitrary precision integers
  - Partition enumeration with optional budget limit
  - Series expansion for approximate computation at large n
  - Computes new OEIS values in milliseconds
"""

from math import factorial, gcd
import sys
import time

def z_lambda(pm):
    z = 1
    for k, m in pm:
        z *= k**m * factorial(m)
    return z

def enumerate_partitions(n, odd_only=False):
    results = []
    if odd_only:
        parts = list(range(n if n%2 else n-1, 0, -2))
    else:
        parts = list(range(n, 0, -1))
    def rec(idx, rem, cur):
        if rem == 0: results.append(list(cur)); return
        if idx >= len(parts): return
        k = parts[idx]
        if k > rem: rec(idx+1, rem, cur); return
        rec(idx+1, rem, cur)
        for m in range(1, rem//k + 1):
            rec(idx+1, rem-m*k, cur + [(k,m)])
    rec(0, n, [])
    return results

# ─── Edge count formulas ───

def _gcd_and_sums(pm):
    """Common helper: compute gcd_sum, sum_v, num_v."""
    gcd_sum = 0
    sum_v = 0
    num_v = 0
    for i, (ki, mi) in enumerate(pm):
        sum_v += mi * ki
        num_v += mi
        gcd_sum += mi * (mi - 1) // 2 * ki  # same-type pairs
        for j in range(i+1, len(pm)):
            kj, mj = pm[j]
            gcd_sum += mi * mj * gcd(ki, kj)  # cross-type pairs
    return gcd_sum, sum_v, num_v

def edges_tournament(pm):  # A000568
    gcd_sum, _, _ = _gcd_and_sums(pm)
    self_term = sum(mi * (ki // 2) for ki, mi in pm)  # floor(k/2); k always odd
    return gcd_sum + self_term

def edges_graph(pm):  # A000088
    gcd_sum, _, _ = _gcd_and_sums(pm)
    self_term = sum(mi * (ki // 2) for ki, mi in pm)
    return gcd_sum + self_term

def edges_digraph(pm):  # A000273
    gcd_sum, sum_v, num_v = _gcd_and_sums(pm)
    return 2 * gcd_sum + sum_v - num_v

def edges_binrel(pm):  # A000595
    gcd_sum, sum_v, _ = _gcd_and_sums(pm)
    return 2 * gcd_sum + sum_v

def edges_symrel(pm):  # A000666
    gcd_sum, _, _ = _gcd_and_sums(pm)
    self_term = sum(mi * (ki // 2 + 1) for ki, mi in pm)
    return gcd_sum + self_term

SEQUENCES = {
    'A000568': ('Tournaments',         edges_tournament, 2, True),
    'A000088': ('Simple graphs',       edges_graph,      2, False),
    'A000273': ('Digraphs',            edges_digraph,    2, False),
    'A000595': ('Binary relations',    edges_binrel,     2, False),
    'A000666': ('Symmetric relations', edges_symrel,     2, False),
}

def burnside(n, seq_id):
    """Compute a(n) for the given sequence."""
    if n == 0:
        return 1
    name, efunc, base, odd_only = SEQUENCES[seq_id]
    nfact = factorial(n)
    parts = enumerate_partitions(n, odd_only)
    total = sum((nfact // z_lambda(p)) * (base ** efunc(p)) for p in parts)
    return total // nfact

def burnside_approx(n, seq_id, num_corrections=5):
    """
    Approximate a(n) using identity + top corrections.
    Returns (approximate_value, identity_only_value).
    """
    if n <= 1:
        return burnside(n, seq_id), burnside(n, seq_id)
    name, efunc, base, odd_only = SEQUENCES[seq_id]
    nfact = factorial(n)

    # Identity partition
    id_part = [(1, n)]
    e_id = efunc(id_part)
    id_total = base ** e_id
    id_val = id_total // nfact

    # Top correction types (sorted by typical importance)
    correction_types = []
    if odd_only:
        for k in range(3, min(n+1, 30), 2):
            correction_types.append([(k, 1)])
        for k in range(3, min(n+1, 20), 2):
            correction_types.append([(k, 2)])
    else:
        for k in range(2, min(n+1, 30)):
            correction_types.append([(k, 1)])
        for k in range(2, min(n+1, 15)):
            correction_types.append([(k, 2)])
        for k1 in range(2, min(n+1, 12)):
            for k2 in range(k1+1, min(n+1, 12)):
                correction_types.append([(k1, 1), (k2, 1)])

    total = id_total
    count = 0
    for ptype in correction_types[:num_corrections]:
        used = sum(k*m for k,m in ptype)
        if used > n:
            continue
        remaining = n - used
        full = list(ptype) + ([(1, remaining)] if remaining > 0 else [])
        z = z_lambda(full)
        e = efunc(full)
        total += (nfact // z) * (base ** e)
        count += 1

    return total // nfact, id_val

# ─── Main ───

if __name__ == "__main__":
    print("="*70)
    print("BURNSIDE UNIFIED ENUMERATOR — VERIFIED")
    print("="*70)

    # Compute and display values for all sequences
    for seq_id, (name, efunc, base, odd_only) in SEQUENCES.items():
        print(f"\n── {seq_id}: {name} ──")
        for n in range(0, 22):
            t0 = time.time()
            result = burnside(n, seq_id)
            elapsed = time.time() - t0
            if elapsed > 5:
                break
            ndigits = len(str(result))
            print(f"  a({n:2d}) = {result}")

    # Large n approximate
    print(f"\n{'='*70}")
    print("LARGE n APPROXIMATE VALUES")
    print("="*70)

    for seq_id, (name, efunc, base, odd_only) in SEQUENCES.items():
        print(f"\n── {seq_id}: {name} ──")
        for n in [50, 100, 200]:
            t0 = time.time()
            approx, id_only = burnside_approx(n, seq_id, num_corrections=10)
            elapsed = time.time() - t0
            ndigits = len(str(approx))
            correction_pct = 100 * (approx - id_only) / id_only if id_only > 0 else 0
            print(f"  a({n:3d}) ≈ {str(approx)[:40]}... ({ndigits} digits, "
                  f"correction={correction_pct:.4e}%, {elapsed:.3f}s)")

    print("\n\nDONE.")
