#!/usr/bin/env python3
"""
opus-2026-04-05-S27: Prime vs composite n — systematic differences.

Multiple levels where prime/composite n matters:

1. BURNSIDE FORMULA: a(n) = (1/n!) Σ_{λ⊢n, odd parts} (n!/z_λ) 2^{e(λ)}
   For PRIME n: the only partition with a part of size n is (n) itself.
   The cycle type (n) contributes iff n is odd. Its edge count e((n)) = (n-1)/2.
   For COMPOSITE n: partitions like (n/p)^p exist, creating additional symmetry.

2. CYCLE STRUCTURE: Cyclic tournaments C_n exist only when n is odd.
   For prime p: C_p = quadratic residue tournament QR_p (unique up to iso).
   QR_p has special spectral properties (circulant, Paley graph connection).

3. METAGRAPH G_n: The isomorphism class count A000568(n) has different
   divisibility properties at prime vs composite n.

4. PATH HOMOLOGY: Tang-Yau (2602.04140) showed that GLMY boundary maps
   for circulant digraphs decompose differently for prime vs composite n.
   At prime p: each eigenspace is 1-dimensional → clean decomposition.
   At composite n: eigenspaces can be higher-dimensional → more complex.

5. SELF-COMPLEMENTARY: SC tournament count has different parity behavior
   at prime vs composite n.

This script explores all five levels computationally.
"""

from math import factorial, gcd, comb, log2
from itertools import combinations, permutations
from collections import Counter, defaultdict
from sympy import isprime, factorint, totient, divisors
import time

# ─── Known values ───
A000568 = {
    1:1, 2:1, 3:2, 4:4, 5:12, 6:56, 7:456, 8:6880,
    9:191536, 10:9733056, 11:903753248, 12:154108311168,
    13:48542114686912, 14:28401423719122304,
    15:31021002160355166848,
    16:63530415842308265100288,
    17:244912778438520759443245824,
    18:1783398846284777975419600287232,
    19:24605641171260376770598003978281472,
    20:645022068557873570931850526424042500096,
}

# ─── Utilities ───

def all_tournaments(n):
    pairs = [(i,j) for i in range(n) for j in range(i+1, n)]
    m = len(pairs)
    for bits in range(2**m):
        T = [[0]*n for _ in range(n)]
        for k, (i,j) in enumerate(pairs):
            if (bits >> k) & 1: T[i][j] = 1
            else: T[j][i] = 1
        yield T

def hamiltonian_path_count(T):
    n = len(T)
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n): dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)) or dp[mask][v] == 0: continue
            for u in range(n):
                if mask & (1 << u): continue
                if T[v][u]: dp[mask | (1 << u)][u] += dp[mask][v]
    return sum(dp[(1 << n) - 1])

def count_3cycles(T):
    n = len(T)
    scores = [sum(T[i]) for i in range(n)]
    return comb(n, 3) - sum(s*(s-1)//2 for s in scores)

def score_sequence(T):
    n = len(T)
    return tuple(sorted(sum(T[i]) for i in range(n)))

def is_self_complementary(T):
    """Check if T is isomorphic to its complement."""
    n = len(T)
    # T^op reverses all arcs
    Top = [[T[j][i] for j in range(n)] for i in range(n)]
    # Check all permutations
    for perm in permutations(range(n)):
        match = True
        for i in range(n):
            for j in range(n):
                if i != j and T[perm[i]][perm[j]] != Top[i][j]:
                    match = False
                    break
            if not match:
                break
        if match:
            return True
    return False

def tournament_to_canon(T):
    n = len(T)
    best = None
    for perm in permutations(range(n)):
        form = tuple(T[perm[i]][perm[j]] for i in range(n) for j in range(n))
        if best is None or form < best:
            best = form
    return best

# ═══════════════════════════════════════════════════════
# LEVEL 1: Burnside formula — partition structure
# ═══════════════════════════════════════════════════════

print("="*70)
print("LEVEL 1: Burnside partition structure — prime vs composite")
print("="*70)

def enumerate_odd_partitions(n):
    results = []
    parts = list(range(n if n % 2 else n-1, 0, -2))
    def recurse(idx, rem, cur):
        if rem == 0:
            results.append(list(cur))
            return
        if idx >= len(parts): return
        k = parts[idx]
        if k > rem:
            recurse(idx+1, rem, cur)
            return
        recurse(idx+1, rem, cur)
        for m in range(1, rem//k + 1):
            recurse(idx+1, rem-m*k, cur + [(k,m)])
    recurse(0, n, [])
    return results

def z_lambda(pm):
    z = 1
    for k, m in pm:
        z *= k**m * factorial(m)
    return z

def edge_count(pm):
    total = 0
    for i, (ki,mi) in enumerate(pm):
        total += mi*(ki-1)//2
        total += mi*(mi-1)//2*ki
        for j in range(i+1, len(pm)):
            kj, mj = pm[j]
            total += mi*mj*gcd(ki,kj)
    return total

for n in range(3, 21):
    parts = enumerate_odd_partitions(n)
    nparts = len(parts)

    # Key: does a partition with a single part (n) exist? Only if n is odd.
    has_full_cycle = (n % 2 == 1)

    # Find the partition with the LARGEST part
    max_part = max(max(k for k, m in p) for p in parts)

    # Count partitions by largest part
    by_max = Counter(max(k for k, m in p) for p in parts)

    prime = isprime(n)
    label = "PRIME" if prime else "comp"

    # The "full cycle" partition (n) contributes:
    if has_full_cycle:
        e_full = edge_count([(n, 1)])
        cn2 = n*(n-1)//2
        # Relative contribution: (n!/z) * 2^e / (2^{C(n,2)})
        # z = n, n!/z = (n-1)!, 2^e = 2^{(n-1)/2}
        # ratio = (n-1)! * 2^{(n-1)/2} / 2^{n(n-1)/2}
        # = (n-1)! * 2^{(n-1)/2 - n(n-1)/2}
        # = (n-1)! * 2^{-(n-1)^2/2}
        delta_e = e_full - cn2
        ratio_log2 = sum(log2(k) for k in range(1, n)) + delta_e  # log2((n-1)!) + delta_e
        full_cycle_info = f"full-cycle: e={e_full}, Δe={delta_e}, log₂(ratio)={ratio_log2:.1f}"
    else:
        full_cycle_info = "no full-cycle (n even)"

    print(f"  n={n:2d} [{label}]: {nparts:4d} partitions, max_part={max_part}, "
          f"by_max_part={dict(sorted(by_max.items(), reverse=True)[:4])}, "
          f"{full_cycle_info}")

# ═══════════════════════════════════════════════════════
# LEVEL 2: H distribution — prime vs composite
# ═══════════════════════════════════════════════════════

print(f"\n{'='*70}")
print("LEVEL 2: H distribution shape — prime vs composite")
print("="*70)

for n in range(3, 9):
    H_vals = []
    for T in all_tournaments(n):
        H_vals.append(hamiltonian_path_count(T))

    H_counter = Counter(H_vals)
    achievable = sorted(H_counter.keys())
    max_H = max(achievable)
    mean_H = sum(H_vals) / len(H_vals)

    # Missing odd values
    all_odd = set(range(1, max_H+1, 2))
    missing = sorted(all_odd - set(achievable))

    # H mod n distribution
    h_mod_n = Counter(h % n for h in H_vals)

    prime = isprime(n)
    label = "PRIME" if prime else "comp"

    print(f"\n  n={n} [{label}]: {len(set(H_vals))} distinct H values, "
          f"max={max_H}, E[H]={mean_H:.2f}")
    print(f"    H values: {achievable}")
    print(f"    Missing odd: {missing}")
    print(f"    H mod {n}: {dict(sorted(h_mod_n.items()))}")

# ═══════════════════════════════════════════════════════
# LEVEL 3: Self-complementary tournaments
# ═══════════════════════════════════════════════════════

print(f"\n{'='*70}")
print("LEVEL 3: Self-complementary count — prime vs composite")
print("="*70)

# SC tournaments exist only for odd n
# A051337: 1, 1, 1, 3, 3, 15, 27, ...
for n in range(3, 10, 2):
    sc_count = 0
    iso_classes = {}
    for T in all_tournaments(n):
        canon = tournament_to_canon(T)
        if canon not in iso_classes:
            iso_classes[canon] = T

    sc_classes = 0
    for canon, T in iso_classes.items():
        if is_self_complementary(T):
            sc_classes += 1

    total_classes = len(iso_classes)
    prime = isprime(n)
    label = "PRIME" if prime else "comp "
    print(f"  n={n} [{label}]: {sc_classes}/{total_classes} SC classes "
          f"({100*sc_classes/total_classes:.1f}%)")

# ═══════════════════════════════════════════════════════
# LEVEL 4: Cycle structure specifics
# ═══════════════════════════════════════════════════════

print(f"\n{'='*70}")
print("LEVEL 4: 3-cycle distribution — prime vs composite")
print("="*70)

for n in range(3, 9):
    c3_vals = []
    for T in all_tournaments(n):
        c3_vals.append(count_3cycles(T))

    c3_counter = Counter(c3_vals)
    max_c3 = max(c3_vals)
    mean_c3 = sum(c3_vals) / len(c3_vals)

    # Theoretical: E[c3] = C(n,3)/4
    expected_c3 = comb(n, 3) / 4

    prime = isprime(n)
    label = "PRIME" if prime else "comp"

    print(f"  n={n} [{label}]: E[c₃]={mean_c3:.2f} (theory={expected_c3:.2f}), "
          f"max={max_c3}, range={sorted(c3_counter.keys())}")

# ═══════════════════════════════════════════════════════
# LEVEL 5: A000568(n) modular properties
# ═══════════════════════════════════════════════════════

print(f"\n{'='*70}")
print("LEVEL 5: a(n) modular properties — prime vs composite")
print("="*70)

for n in range(3, 21):
    a = A000568.get(n)
    if a is None: continue
    prime = isprime(n)
    label = "PRIME" if prime else "comp"

    # a(n) mod small primes
    mods = {p: a % p for p in [2, 3, 5, 7, 11, 13]}

    # a(n) mod n
    a_mod_n = a % n

    # v_2(a(n)) = 2-adic valuation
    v2 = 0
    temp = a
    while temp % 2 == 0 and temp > 0:
        v2 += 1
        temp //= 2

    # Ratio a(n)/a(n-1)
    a_prev = A000568.get(n-1)
    ratio = a / a_prev if a_prev else 0

    print(f"  n={n:2d} [{label}]: a mod n={a_mod_n:4d}, v₂={v2:2d}, "
          f"mod2={mods[2]} mod3={mods[3]} mod5={mods[5]} mod7={mods[7]}, "
          f"ratio={ratio:.4f}")

# ═══════════════════════════════════════════════════════
# LEVEL 6: The "n-cycle" partition — prime n special
# ═══════════════════════════════════════════════════════

print(f"\n{'='*70}")
print("LEVEL 6: The n-cycle partition contribution — prime vs composite")
print("="*70)

print("\nFor ODD n, the partition (n) [single n-cycle] contributes:")
print("  (n-1)! × 2^{(n-1)/2} to the Burnside sum.")
print("  This is the CYCLIC TOURNAMENT contribution (QR_p for prime p).\n")

for n in range(3, 22, 2):
    # Single n-cycle partition
    e_n = (n-1) // 2
    z_n = n
    nfact = factorial(n)
    cn2 = n*(n-1)//2

    # Contribution to Burnside sum
    contrib = (nfact // z_n) * (2 ** e_n)

    # As fraction of identity
    identity = 2 ** cn2
    ratio = contrib / identity

    # As fraction of total (approximate)
    a_n = A000568.get(n)
    if a_n:
        frac_of_total = contrib / (a_n * nfact)
    else:
        frac_of_total = ratio  # approximate

    prime = isprime(n)
    label = "PRIME" if prime else "comp"

    # For PRIME p: (p) is the ONLY partition with largest part p
    # For COMPOSITE n: there are also partitions like (n/d)^d

    # Count divisor-related partitions
    divs = [d for d in range(3, n+1, 2) if n % d == 0 and d < n]
    div_parts = []
    for d in divs:
        m = n // d
        if d % 2 == 1:  # d must be odd for odd-part partition
            div_parts.append(f"({d}^{m})")

    print(f"  n={n:2d} [{label}]: e={(n-1)//2:3d}, ratio_to_id={ratio:.4e}, "
          f"frac_of_a(n)={frac_of_total:.4e}"
          f"{', divisor_parts: ' + ', '.join(div_parts) if div_parts else ''}")

# ═══════════════════════════════════════════════════════
# LEVEL 7: The KEY difference — divisor sum structure
# ═══════════════════════════════════════════════════════

print(f"\n{'='*70}")
print("LEVEL 7: Divisor sum in edge formula — prime vs composite")
print("="*70)

print("""
The edge count uses gcd(k_i, k_j) between cycle lengths.
For PRIME p: if k_i = p, then gcd(p, k_j) = p if k_j = p, else 1.
  → The p-cycle interacts MINIMALLY with other parts (only gcd=1 cross-terms).
  → Clean "orthogonal" contribution.

For COMPOSITE n = ab: gcd(a,b) > 1 when both parts share factors.
  → More complex interactions between cycle types.
  → Cross-terms create "entanglement" in the Burnside sum.
""")

# Compute: for each n, the average gcd between pairs of odd divisors
for n in range(3, 25):
    odd_divs = [d for d in range(1, n+1, 2) if n % d == 0]
    if len(odd_divs) < 2:
        avg_gcd = 0
    else:
        gcd_sum = sum(gcd(a, b) for a in odd_divs for b in odd_divs if a < b)
        n_pairs = len(odd_divs) * (len(odd_divs) - 1) // 2
        avg_gcd = gcd_sum / n_pairs if n_pairs > 0 else 0

    prime = isprime(n)
    label = "P" if prime else "C"
    divs_str = ",".join(map(str, odd_divs))
    print(f"  n={n:2d} [{label}]: odd divisors={{{divs_str}:16s}}, "
          f"avg_gcd={avg_gcd:.2f}, #divs={len(odd_divs)}")

# ═══════════════════════════════════════════════════════
# LEVEL 8: Quantitative — correction term ratios
# ═══════════════════════════════════════════════════════

print(f"\n{'='*70}")
print("LEVEL 8: 3-cycle correction relative size — prime vs composite")
print("="*70)

print("  The 3-cycle correction ratio R₃(n) = 16n(n-1)(n-2)/(3×4^n)")
print("  This is the SAME formula for prime and composite n.")
print("  But the HIGHER corrections differ!\n")

for n in range(3, 25):
    R3 = 16 * n * (n-1) * (n-2) / (3 * 4**n)

    # For composite odd n with factor 3: the partition (3^{n/3}) exists
    # and contributes more than at prime n (where no such regular partition exists)
    has_3_power = (n % 3 == 0 and n > 3)

    # The (3^{n/3}) partition:
    if has_3_power and n % 2 == 1:
        m3 = n // 3
        full = [(3, m3)]
        e3 = edge_count(full)
        cn2 = n*(n-1)//2
        z3 = z_lambda(full)
        nfact = factorial(n)
        ratio_3power = (nfact // z3) * 2**(e3 - cn2) if e3 >= cn2 else (nfact // z3) / 2**(cn2 - e3)
        ratio_3power_float = (nfact / z3) * 2**(e3 - cn2)
    else:
        ratio_3power_float = 0

    prime = isprime(n)
    label = "P" if prime else "C"
    extra = f", R(3^{n//3})={ratio_3power_float:.4e}" if has_3_power and n%2==1 else ""
    print(f"  n={n:2d} [{label}]: R₃={R3:.6e}{extra}")

print("\n\nDONE.")
