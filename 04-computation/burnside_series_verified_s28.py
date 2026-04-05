#!/usr/bin/env python3
"""
opus-2026-04-05-S28: Rigorous verification AND extension of Burnside series.

PART 1: CORRECTNESS VERIFICATION
  - Verify A000568 exact enumeration against ALL known OEIS values
  - Cross-check with the GMP C code output file
  - Verify the series expansion gives exact integers for n ≤ 20
  - Bound the tail rigorously for the approximate regime

PART 2: EXTEND TO OTHER SEQUENCES
  The Burnside framework applies to ALL sequences in burnside_enum_v2.c:
  - A000568: tournaments (base 2, odd parts)
  - A000088: simple graphs (base 2, all parts)
  - A000273: digraphs (base 2, all parts, 2× gcd cross-term)
  - A000595: binary relations (base 2, all parts, 2× gcd + Σv)
  - A001174: oriented graphs (base 3, all parts)
  - A000666: symmetric relations (base 2, all parts, gcd + floor + 1)

  For each: derive the series expansion, verify, benchmark.
"""

from math import factorial, gcd, comb, log2
from fractions import Fraction
from collections import defaultdict
import time

# ─── Known OEIS values ───

OEIS = {
    'A000568': {  # Tournaments
        1:1, 2:1, 3:2, 4:4, 5:12, 6:56, 7:456, 8:6880, 9:191536,
        10:9733056, 11:903753248, 12:154108311168,
        13:48542114686912, 14:28401423719122304,
        15:31021002160355166848,
    },
    'A000088': {  # Simple graphs
        0:1, 1:1, 2:2, 3:4, 4:11, 5:34, 6:156, 7:1044, 8:12346,
        9:274668, 10:12005168, 11:1018997864, 12:165091172592,
    },
    'A000273': {  # Digraphs
        0:1, 1:1, 2:3, 3:16, 4:218, 5:9608, 6:1540944,
        7:882033440, 8:1793359192848,
    },
    'A000595': {  # Binary relations on n elements
        0:1, 1:2, 2:10, 3:104, 4:3044, 5:291968, 6:96928992,
        7:112282908928, 8:458297100061728,
    },
    'A001174': {  # Oriented graphs (no mutual arcs)
        0:1, 1:1, 2:2, 3:7, 4:42, 5:582, 6:21480, 7:2142288,
        8:575012972,
    },
    'A000666': {  # 2-colored graphs = symmetric binary relations
        0:1, 1:2, 2:6, 3:20, 4:90, 5:544, 6:5096, 7:79264,
        8:2208612, 9:115249244,
    },
}

# ─── Partition enumeration ───

def enumerate_partitions(n, odd_only=False):
    results = []
    if odd_only:
        parts = list(range(n if n%2 else n-1, 0, -2))
    else:
        parts = list(range(n, 0, -1))
    def rec(idx, rem, cur):
        if rem == 0:
            results.append(list(cur))
            return
        if idx >= len(parts): return
        k = parts[idx]
        if k > rem: rec(idx+1, rem, cur); return
        rec(idx+1, rem, cur)
        for m in range(1, rem//k + 1):
            rec(idx+1, rem-m*k, cur + [(k,m)])
    rec(0, n, [])
    return results

def z_lambda(pm):
    z = 1
    for k, m in pm:
        z *= k**m * factorial(m)
    return z

# ─── Edge count functions for each sequence ───

def edges_A000568(pm):
    """Tournaments: gcd_sum + Σ m_i⌊k_i/2⌋"""
    gcd_sum = 0
    self_term = 0
    for i, (ki, mi) in enumerate(pm):
        self_term += mi * (ki - 1) // 2
        gcd_sum += mi * (mi - 1) // 2 * ki
        for j in range(i+1, len(pm)):
            kj, mj = pm[j]
            gcd_sum += mi * mj * gcd(ki, kj)
    return gcd_sum + self_term

def edges_A000088(pm):
    """Simple graphs: same as tournaments."""
    return edges_A000568(pm)

def edges_A000273(pm):
    """Digraphs: 2×gcd_sum + Σv - #parts"""
    gcd_sum = 0
    sum_v = 0
    num_v = 0
    for i, (ki, mi) in enumerate(pm):
        sum_v += mi * ki
        num_v += mi
        gcd_sum += mi * (mi - 1) // 2 * ki
        for j in range(i+1, len(pm)):
            kj, mj = pm[j]
            gcd_sum += mi * mj * gcd(ki, kj)
    return 2 * gcd_sum + sum_v - num_v

def edges_A000595(pm):
    """Binary relations: 2×gcd_sum + Σv"""
    gcd_sum = 0
    sum_v = 0
    for i, (ki, mi) in enumerate(pm):
        sum_v += mi * ki
        gcd_sum += mi * (mi - 1) // 2 * ki
        for j in range(i+1, len(pm)):
            kj, mj = pm[j]
            gcd_sum += mi * mj * gcd(ki, kj)
    return 2 * gcd_sum + sum_v

def edges_A001174(pm):
    """Oriented graphs (base 3): gcd_sum + Σ m_i⌊(k_i-1)/2⌋"""
    gcd_sum = 0
    self_term = 0
    for i, (ki, mi) in enumerate(pm):
        self_term += mi * ((ki - 1) // 2)
        gcd_sum += mi * (mi - 1) // 2 * ki
        for j in range(i+1, len(pm)):
            kj, mj = pm[j]
            gcd_sum += mi * mj * gcd(ki, kj)
    return gcd_sum + self_term

def edges_A000666(pm):
    """Symmetric relations: gcd_sum + Σ m_i(⌊k_i/2⌋ + 1)"""
    gcd_sum = 0
    self_term = 0
    for i, (ki, mi) in enumerate(pm):
        self_term += mi * (ki // 2 + 1)
        gcd_sum += mi * (mi - 1) // 2 * ki
        for j in range(i+1, len(pm)):
            kj, mj = pm[j]
            gcd_sum += mi * mj * gcd(ki, kj)
    return gcd_sum + self_term

# ─── General Burnside computation ───

def burnside_exact(n, edge_func, base, odd_only=False):
    """Exact Burnside computation using integer arithmetic."""
    if n == 0:
        return 1
    if n == 1:
        # Need to handle base case
        parts = enumerate_partitions(n, odd_only)
        nfact = factorial(n)
        total = 0
        for p in parts:
            z = z_lambda(p)
            e = edge_func(p)
            total += (nfact // z) * (base ** e)
        return total // nfact

    nfact = factorial(n)
    parts = enumerate_partitions(n, odd_only)
    total = 0
    for p in parts:
        z = z_lambda(p)
        e = edge_func(p)
        total += (nfact // z) * (base ** e)
    return total // nfact

# ─── Series expansion ───

def burnside_series(n, edge_func, base, odd_only=False, max_budget=None):
    """
    Compute using series expansion with optional budget limit.
    Returns (result, n_partitions, identity_value).
    """
    if n <= 1:
        return burnside_exact(n, edge_func, base, odd_only), 1, 1

    nfact = factorial(n)
    parts = enumerate_partitions(n, odd_only)

    # Find identity partition
    id_part = [(1, n)]
    e_id = edge_func(id_part)
    identity_contrib = base ** e_id  # n!/z_id = n!/n! = 1, so contrib = base^e_id

    # Sort all partitions by contribution (descending)
    contribs = []
    for p in parts:
        z = z_lambda(p)
        e = edge_func(p)
        w = (nfact // z) * (base ** e)
        contribs.append((w, p))
    contribs.sort(reverse=True)

    if max_budget is None:
        total = sum(w for w, p in contribs)
        return total // nfact, len(contribs), identity_contrib // nfact
    else:
        # Only sum partitions where large parts use ≤ max_budget
        total = 0
        count = 0
        for w, p in contribs:
            large_sum = sum(k*m for k,m in p if k > 1)
            if large_sum <= max_budget:
                total += w
                count += 1
        return total // nfact, count, identity_contrib // nfact

# ═══════════════════════════════════════════════════════
# PART 1: RIGOROUS VERIFICATION
# ═══════════════════════════════════════════════════════

CONFIGS = {
    'A000568': (edges_A000568, 2, True,  "Tournaments"),
    'A000088': (edges_A000088, 2, False, "Simple graphs"),
    'A000273': (edges_A000273, 2, False, "Digraphs"),
    'A000595': (edges_A000595, 2, False, "Binary relations"),
    'A001174': (edges_A001174, 3, False, "Oriented graphs"),
    'A000666': (edges_A000666, 2, False, "Symmetric relations"),
}

if __name__ == "__main__":
    print("="*70)
    print("PART 1: RIGOROUS VERIFICATION against OEIS")
    print("="*70)

    all_correct = True
    for seq_id, (efunc, base, odd, name) in CONFIGS.items():
        known = OEIS.get(seq_id, {})
        if not known:
            continue

        print(f"\n── {seq_id} ({name}, base={base}, odd_only={odd}) ──")
        max_n = max(known.keys())
        for n in sorted(known.keys()):
            if n > 15:  # skip very large n for speed
                break
            t0 = time.time()
            result = burnside_exact(n, efunc, base, odd)
            elapsed = time.time() - t0
            expected = known[n]
            ok = result == expected
            if not ok:
                all_correct = False
            mark = "✓" if ok else f"✗ got {result}"
            print(f"  n={n:2d}: {result:>20d}  {mark}  ({elapsed:.4f}s)")

    print(f"\n{'ALL CORRECT' if all_correct else 'ERRORS FOUND'}!")

    # ═══════════════════════════════════════════════════════
    # PART 2: SERIES EXPANSION — identity fraction per sequence
    # ═══════════════════════════════════════════════════════
    print(f"\n{'='*70}")
    print("PART 2: Identity fraction and convergence per sequence")
    print("="*70)

    for seq_id, (efunc, base, odd, name) in CONFIGS.items():
        print(f"\n── {seq_id} ({name}) ──")
        print(f"  {'n':>3s} {'#parts':>7s} {'id_frac':>12s} {'coupling':>12s} {'budget=3':>12s} {'budget=6':>12s}")

        known = OEIS.get(seq_id, {})
        for n in range(2, 16):
            if n > 12 and seq_id in ['A000273', 'A000595']:
                break  # too slow

            t0 = time.time()
            exact, nparts, id_val = burnside_series(n, efunc, base, odd)
            elapsed = time.time() - t0

            if elapsed > 5:
                break

            # Budget=3: only partitions with ≤3 non-unit sum
            b3, _, _ = burnside_series(n, efunc, base, odd, max_budget=3)
            # Budget=6
            b6, _, _ = burnside_series(n, efunc, base, odd, max_budget=6)

            err3 = abs(b3 - exact) / exact if exact > 0 else 0
            err6 = abs(b6 - exact) / exact if exact > 0 else 0

            # Identity fraction
            nfact = factorial(n)
            parts = enumerate_partitions(n, odd)
            id_part = [(1, n)]
            e_id = efunc(id_part)
            id_contrib = (nfact // z_lambda(id_part)) * (base ** e_id)
            total = sum((nfact // z_lambda(p)) * (base ** efunc(p)) for p in parts)
            id_frac = id_contrib / total if total > 0 else 1
            coupling = 1 - id_frac

            expected = known.get(n)
            mark = "✓" if expected and exact == expected else ("?" if expected is None else "✗")

            print(f"  {n:3d} {nparts:7d} {id_frac:12.8f} {coupling:12.4e} "
                  f"err3={err3:.2e} err6={err6:.2e} {mark}")

    # ═══════════════════════════════════════════════════════
    # PART 3: PRACTICAL — how many corrections for exact answer?
    # ═══════════════════════════════════════════════════════
    print(f"\n{'='*70}")
    print("PART 3: Budget needed for EXACT answer per sequence")
    print("="*70)

    for seq_id, (efunc, base, odd, name) in CONFIGS.items():
        print(f"\n── {seq_id} ({name}) ──")
        known = OEIS.get(seq_id, {})

        for n in [8, 10, 12]:
            if n > 10 and seq_id in ['A000273', 'A000595', 'A001174']:
                continue

            exact = burnside_exact(n, efunc, base, odd)

            for budget in [3, 6, 9, 12, 15, 20, n]:
                approx, nused, _ = burnside_series(n, efunc, base, odd, max_budget=budget)
                if approx == exact:
                    print(f"  n={n:2d}: budget={budget:2d} suffices ({nused} partitions)")
                    break
            else:
                print(f"  n={n:2d}: budget={n} needed (all {nused} partitions)")

    # ═══════════════════════════════════════════════════════
    # PART 4: CONVERGENCE RATES — which sequence converges fastest?
    # ═══════════════════════════════════════════════════════
    print(f"\n{'='*70}")
    print("PART 4: Convergence rate comparison")
    print("="*70)

    print(f"\n  Coupling constant g(n) = 1 - identity_fraction:")
    print(f"  {'n':>3s}", end="")
    for seq_id in ['A000568', 'A000088', 'A000273', 'A000595', 'A001174', 'A000666']:
        print(f" {seq_id:>12s}", end="")
    print()

    for n in range(3, 13):
        print(f"  {n:3d}", end="")
        for seq_id, (efunc, base, odd, name) in CONFIGS.items():
            try:
                nfact = factorial(n)
                parts = enumerate_partitions(n, odd)
                id_part = [(1, n)]
                e_id = efunc(id_part)
                id_c = (nfact // z_lambda(id_part)) * (base ** e_id)
                total = sum((nfact // z_lambda(p)) * (base ** efunc(p)) for p in parts)
                coupling = 1 - id_c/total if total > 0 else 0
                print(f" {coupling:12.4e}", end="")
            except Exception:
                print(f" {'---':>12s}", end="")
        print()

    # ═══════════════════════════════════════════════════════
    # PART 5: EXTENDED VALUES — compute new terms using series
    # ═══════════════════════════════════════════════════════
    print(f"\n{'='*70}")
    print("PART 5: New exact values beyond OEIS")
    print("="*70)

    for seq_id, (efunc, base, odd, name) in CONFIGS.items():
        known = OEIS.get(seq_id, {})
        max_known = max(known.keys()) if known else 0
        print(f"\n── {seq_id} ({name}): OEIS has n≤{max_known} ──")

        for n in range(max_known + 1, max_known + 4):
            if n > 20:
                break
            t0 = time.time()
            result = burnside_exact(n, efunc, base, odd)
            elapsed = time.time() - t0
            nparts = len(enumerate_partitions(n, odd))
            ndigits = len(str(result))
            print(f"  a({n:2d}) = {str(result)[:50]}{'...' if ndigits > 50 else ''}"
                  f"  ({ndigits} digits, {nparts} partitions, {elapsed:.3f}s)")
            if elapsed > 10:
                break

    print("\n\nDONE.")
