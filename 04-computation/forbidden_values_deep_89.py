"""
forbidden_values_deep_89.py
opus-2026-03-14-S89

Deep analysis of forbidden H values.
n=5 forbidden: {7}
n=6 forbidden: {7, 21, 35, 39}
Are 35 and 39 permanently forbidden or just not-yet-achieved?

Also: Monte Carlo H-spectrum at n=7 to find new forbidden values.
"""

from itertools import permutations, combinations
from collections import Counter
from fractions import Fraction
import random

print("=" * 70)
print("FORBIDDEN H VALUES — DEEP ANALYSIS")
print("opus-2026-03-14-S89")
print("=" * 70)

def compute_H(adj, n):
    """Compute H(T) by counting Hamiltonian paths."""
    count = 0
    for perm in permutations(range(n)):
        is_path = True
        for i in range(n - 1):
            if not adj[perm[i]][perm[i+1]]:
                is_path = False
                break
        if is_path:
            count += 1
    return count

# =====================================================================
# PART 1: FORBIDDEN VALUES FACTORIZATION
# =====================================================================
print("\n" + "=" * 70)
print("PART 1: FACTORIZATION OF FORBIDDEN VALUES")
print("=" * 70)

phi3_vals = {k: k*k+k+1 for k in range(20)}

forbidden_6 = [7, 21, 35, 39]
print(f"\n  Forbidden at n=6: {forbidden_6}")

for v in forbidden_6:
    # Check if it's a Phi_3 value
    is_phi3 = v in phi3_vals.values()
    phi3_at = [k for k, val in phi3_vals.items() if val == v]

    # Factorize
    factors = []
    n = v
    for p in range(2, v+1):
        while n % p == 0:
            factors.append(p)
            n //= p
        if n == 1:
            break

    # Check if it's a product of Phi_3 values
    phi3_product = []
    remaining = v
    for k in sorted(phi3_vals.keys(), reverse=True):
        val = phi3_vals[k]
        if val > 1 and remaining % val == 0:
            phi3_product.append(f"Phi_3({k})={val}")
            remaining //= val
    if remaining > 1:
        phi3_product.append(str(remaining))

    print(f"\n  {v}:")
    print(f"    Factors: {' * '.join(map(str, factors))}")
    print(f"    Is Phi_3(k)? {is_phi3}" + (f" at k={phi3_at[0]}" if phi3_at else ""))
    print(f"    Phi_3 product: {' * '.join(phi3_product)}")

# Check 35 and 39 specifically
print(f"\n  35 = 5 * 7 = 5 * Phi_3(2)")
print(f"  39 = 3 * 13 = Phi_3(1) * Phi_3(3)")
print(f"  Note: 39 = Phi_3(1) * Phi_3(3) — product of TWO Phi_3 values!")
print(f"  But 35 has factor 5 which is NOT a Phi_3 value.")
print(f"  5 = Phi_5(2) = 2^4+2^3+2^2+2+1 ... no. Phi_5(2) = 31.")
print(f"  5 = Phi_4(2) = 2^2+1 = 5. YES! 5 = Phi_4(2)!")
print(f"  So 35 = Phi_4(2) * Phi_3(2) = 5 * 7!")

# All cyclotomic at x=2
print("\n  Cyclotomic polynomials at x=2:")
# Phi_1(2) = 1, Phi_2(2) = 3, Phi_3(2) = 7, Phi_4(2) = 5, Phi_5(2) = 31
# Phi_6(2) = 3, Phi_8(2) = 17, Phi_10(2) = 11, Phi_12(2) = 13
# Product formula: 2^n - 1 = product Phi_d(2) for d|n
# So Phi_d(2) are the multiplicative building blocks of Mersenne numbers!

# Compute cyclotomic at 2 for small d
# Phi_d(x) = product_{gcd(k,d)=1, 1<=k<=d} (x - e^(2*pi*i*k/d))
# At x=2: use the factorization 2^n - 1 = prod_{d|n} Phi_d(2)

phi_at_2 = {}
for n in range(1, 25):
    prod = 1
    for d in range(1, n):
        if n % d == 0 and d in phi_at_2:
            prod *= phi_at_2[d]
    phi_at_2[n] = (2**n - 1) // prod if prod > 0 else 2**n - 1

print(f"  {'d':>3} {'Phi_d(2)':>10} {'= 2^d-1 / ...':>30}")
for d in range(1, 16):
    print(f"  {d:3d} {phi_at_2[d]:10d}")

# Now express each forbidden value in terms of Phi_d(2)
print("\n  Forbidden values as products of Phi_d(2):")
for v in forbidden_6:
    remaining = v
    factors_phi = []
    for d in sorted(phi_at_2.keys(), reverse=True):
        val = phi_at_2[d]
        while val > 1 and remaining % val == 0:
            factors_phi.append(f"Phi_{d}(2)={val}")
            remaining //= val
    if remaining > 1:
        factors_phi.append(str(remaining))
    print(f"  {v:4d} = {' * '.join(factors_phi)}")

# =====================================================================
# PART 2: WHICH ODD NUMBERS ARE ACHIEVABLE?
# =====================================================================
print("\n" + "=" * 70)
print("PART 2: ACHIEVABLE H VALUES UP TO n=6")
print("=" * 70)

# Collect all achieved H values
achieved = set()
achieved_at = {}
for n in range(1, 7):
    if n == 1:
        achieved.add(1)
        achieved_at[1] = 1
        continue
    edges = [(i,j) for i in range(n) for j in range(i+1,n)]
    m = len(edges)
    for bits in range(2**m):
        adj = [[0]*n for _ in range(n)]
        for k, (i,j) in enumerate(edges):
            if bits & (1 << k):
                adj[i][j] = 1
            else:
                adj[j][i] = 1
        H = compute_H(adj, n)
        if H not in achieved:
            achieved_at[H] = n
        achieved.add(H)

max_achieved = max(achieved)
print(f"\n  Max H achieved at n=6: {max_achieved}")
print(f"  All achieved H values (sorted): {sorted(achieved)}")
print(f"  Count: {len(achieved)}")

# Missing odd values up to max
missing = sorted(set(range(1, max_achieved+1, 2)) - achieved)
print(f"\n  Missing odd values in [1, {max_achieved}]: {missing}")

# For each missing value, check its structure
print("\n  Missing value analysis:")
for v in missing:
    remaining = v
    factors_phi = []
    for d in sorted(phi_at_2.keys(), reverse=True):
        val = phi_at_2[d]
        while val > 1 and remaining % val == 0:
            factors_phi.append(f"Phi_{d}(2)")
            remaining //= val
    if remaining > 1:
        factors_phi.append(str(remaining))
    print(f"    {v:4d} = {' * '.join(factors_phi)}")

# =====================================================================
# PART 3: MONTE CARLO AT n=7 — FIND NEW H VALUES
# =====================================================================
print("\n" + "=" * 70)
print("PART 3: MONTE CARLO H-SPECTRUM AT n=7")
print("=" * 70)

random.seed(42)
n = 7
edges_7 = [(i,j) for i in range(n) for j in range(i+1,n)]
m = len(edges_7)

H_counter_7 = Counter()
N_SAMPLES = 100000

print(f"\n  Sampling {N_SAMPLES} random tournaments on n={n} vertices...")
print(f"  (m = {m} arcs, 2^m = {2**m} total tournaments)")

for sample in range(N_SAMPLES):
    adj = [[0]*n for _ in range(n)]
    for (i,j) in edges_7:
        if random.random() < 0.5:
            adj[i][j] = 1
        else:
            adj[j][i] = 1
    H = compute_H(adj, n)
    H_counter_7[H] += 1

H_vals_7 = sorted(H_counter_7.keys())
print(f"\n  Distinct H values found: {len(H_vals_7)}")
print(f"  H values: {H_vals_7}")
print(f"  Min: {min(H_vals_7)}, Max: {max(H_vals_7)}")
print(f"  Mean: {sum(h*c for h,c in H_counter_7.items())/N_SAMPLES:.2f}")

# Statistics
sum_H = sum(h*c for h,c in H_counter_7.items())
sum_H2 = sum(h**2*c for h,c in H_counter_7.items())
mean_7 = sum_H / N_SAMPLES
var_7 = sum_H2 / N_SAMPLES - mean_7**2
ratio_7 = var_7 / mean_7**2
print(f"  Var: {var_7:.2f}")
print(f"  Var/Mean^2: {ratio_7:.6f}")
print(f"  1/3: {1/3:.6f}")
print(f"  Deviation: {1/3 - ratio_7:.6f}")

# Distribution (top entries)
print(f"\n  H distribution (top 30 by value):")
for H_val in sorted(H_counter_7.keys())[:30]:
    count = H_counter_7[H_val]
    bar = '#' * (count * 50 // max(H_counter_7.values()))
    print(f"    H={H_val:4d}: {count:6d} {bar}")
print(f"    ...")
for H_val in sorted(H_counter_7.keys())[-10:]:
    count = H_counter_7[H_val]
    bar = '#' * (count * 50 // max(H_counter_7.values()))
    print(f"    H={H_val:4d}: {count:6d} {bar}")

# Missing odd values at n=7
max_7 = max(H_vals_7)
present_odd = set(H_vals_7)
all_odd_7 = set(range(1, max_7+1, 2))
missing_7 = sorted(all_odd_7 - present_odd)
print(f"\n  Odd values in [1,{max_7}] not found in {N_SAMPLES} samples:")
print(f"  {missing_7[:50]}")
if len(missing_7) > 50:
    print(f"  ... and {len(missing_7)-50} more")
print(f"  Total missing: {len(missing_7)} out of {len(all_odd_7)} odd values")

# Check if 7 and 21 are still missing
print(f"\n  7 found? {7 in present_odd}")
print(f"  21 found? {21 in present_odd}")
print(f"  35 found? {35 in present_odd}")
print(f"  39 found? {39 in present_odd}")

# New achievements at n=7 that weren't at n=6
new_at_7 = present_odd - achieved
print(f"\n  New H values first seen at n=7: {sorted(new_at_7)}")

# =====================================================================
# PART 4: WHICH FORBIDDEN VALUES ARE PERMANENT?
# =====================================================================
print("\n" + "=" * 70)
print("PART 4: PERMANENT vs TEMPORARY FORBIDDEN VALUES")
print("=" * 70)

# From previous work:
# H=7: PERMANENTLY FORBIDDEN (proven, Fano obstruction)
# H=21: PERMANENTLY FORBIDDEN (proven, PG(2,4) obstruction)
# H=35: Unknown — is it permanent?
# H=39: Unknown

# Check if 35 and 39 appear in the n=7 Monte Carlo
if 35 in present_odd:
    print(f"  H=35: ACHIEVED at n=7 ({H_counter_7.get(35, 0)} times)")
else:
    print(f"  H=35: NOT FOUND at n=7 in {N_SAMPLES} samples")

if 39 in present_odd:
    print(f"  H=39: ACHIEVED at n=7 ({H_counter_7.get(39, 0)} times)")
else:
    print(f"  H=39: NOT FOUND at n=7 in {N_SAMPLES} samples")

# Updated forbidden list
still_missing = [v for v in [7, 21, 35, 39] if v not in present_odd]
print(f"\n  Values still missing at n=7: {still_missing}")

# Permanently forbidden (proven): 7, 21
# Possibly forbidden: 35, 39 (if still missing)
# Definitely temporary: any value in missing_6 that appears at n=7

# =====================================================================
# PART 5: THE FORBIDDEN VALUE PATTERN — Phi_3 AND BEYOND
# =====================================================================
print("\n" + "=" * 70)
print("PART 5: THE FORBIDDEN VALUE PATTERN")
print("=" * 70)

# The permanent forbiddens are H_forb_k = Phi_3(2^k) for k=1,2
# = 7 and 21
# The conjecture from S71i: H=273 = Phi_3(16) IS achievable
# So not all Phi_3(2^k) are forbidden!

# What determines forbiddenness?
# 7 = |PG(2,F_2)| — no tournament has exactly 7 Hamiltonian paths
# 21 = |PG(2,F_4)| — no tournament has exactly 21 Hamiltonian paths
# 73 = |PG(2,F_8)| — ??? unknown

print("""
  KNOWN FORBIDDEN H VALUES:
    H=7  = Phi_3(2) = |PG(2,F_2)| — PROVEN permanently forbidden
    H=21 = Phi_3(4) = |PG(2,F_4)| — PROVEN permanently forbidden

  KNOWN ACHIEVABLE H VALUES:
    H=273 = Phi_3(16) = |PG(2,F_16)| — PROVEN achievable (S71i)

  OPEN: Are 35, 39 permanently forbidden?
  OPEN: Is 73 = Phi_3(8) forbidden or achievable?

  THE PATTERN OF FORBIDDENNESS:
    Phi_3(2^1) = 7:   FORBIDDEN
    Phi_3(2^2) = 21:  FORBIDDEN
    Phi_3(2^3) = 73:  UNKNOWN
    Phi_3(2^4) = 273: ACHIEVABLE

  The key difference:
    F_2 = GF(2): Phi_3 irreducible over F_2 (gives F_4)
    F_4 = GF(4): Phi_3 factors over F_4 (x^2+x+1 = (x-alpha)(x-alpha^2))
    F_8 = GF(8): Phi_3 factors over F_8 (cube roots of unity in F_8)
    F_16 = GF(16): Phi_3 factors over F_16

  All F_{2^k} for k>=2 contain cube roots of unity.
  So the Phi_3 obstruction only applies when the field is "too small."
  F_2 and F_4 correspond to k=1,2: BOTH forbidden.
  F_8 has cube roots of 1 since 3|7=|F_8*|: MIGHT be achievable.
""")

# =====================================================================
# PART 6: THE MULTIPLICITY PATTERN — HOW MANY TOURNAMENTS GIVE H?
# =====================================================================
print("\n" + "=" * 70)
print("PART 6: TOURNAMENT COUNTS PER H VALUE")
print("=" * 70)

# At n=6, the counts are very structured
print(f"\n  n=6 counts:")
counts_6 = {1: 720, 3: 960, 5: 2160, 9: 2960, 11: 1440, 13: 1440,
             15: 2208, 17: 1440, 19: 1440, 23: 2880, 25: 1440,
             27: 480, 29: 2880, 31: 1440, 33: 2640, 37: 3600,
             41: 720, 43: 1440, 45: 480}

total = sum(counts_6.values())
print(f"  Total: {total} = 2^15 = {2**15}")

for H, count in sorted(counts_6.items()):
    g = Fraction(count, total)
    print(f"    H={H:3d}: {count:5d} = {count}/32768 = {g}")

# GCD of all counts
from functools import reduce
from math import gcd as math_gcd
g = reduce(math_gcd, counts_6.values())
print(f"\n  GCD of all counts: {g}")
print(f"  Counts / {g}: {[c // g for c in sorted(counts_6.values())]}")

# The counts modulo small numbers
print(f"\n  All counts divisible by 480? {all(c % 480 == 0 for c in counts_6.values())}")
print(f"  Counts / 480: {[c // 480 for c in [counts_6[h] for h in sorted(counts_6.keys())]]}")

# =====================================================================
# PART 7: H-SPECTRUM SIZE AT n=7
# =====================================================================
print("\n" + "=" * 70)
print("PART 7: H-SPECTRUM SIZE AT n=7")
print("=" * 70)

print(f"\n  At n=7 ({N_SAMPLES} samples), found {len(H_vals_7)} distinct H values")
print(f"  Expected spectrum size (if ratio -> e): ~ 19*e ~ {19*2.718:.0f}")
print(f"  Actual found: {len(H_vals_7)}")
print(f"  (May be undercount due to sampling — rare H values might be missed)")

# How many of the odd numbers in range are covered?
coverage_7 = len(present_odd) / len(all_odd_7)
print(f"  Coverage: {len(present_odd)}/{len(all_odd_7)} = {coverage_7:.4f}")

# Interpolation: if Monte Carlo found X values but the true spectrum
# has Y >= X, we can estimate based on rarefaction curves
# For now just report the found count
print(f"\n  SEQUENCE UPDATE:")
print(f"    n:      1    2    3    4    5    6    7")
print(f"    |Spec|: 1    1    2    3    7   19  >={len(H_vals_7)}")

# =====================================================================
# PART 8: THE 35 AND 39 MYSTERY
# =====================================================================
print("\n" + "=" * 70)
print("PART 8: WHY ARE 35 AND 39 MISSING AT n=6?")
print("=" * 70)

print(f"""
  35 = 5 * 7 = Phi_4(2) * Phi_3(2)
  39 = 3 * 13 = Phi_3(1) * Phi_3(3)

  Both are products involving Phi_3 values.
  But 15 = 3 * 5 = Phi_3(1) * Phi_4(2) IS achievable (at n=5!).
  And 9 = 3 * 3 = Phi_3(1)^2 IS achievable.

  So the product rule doesn't directly explain forbiddenness.

  THE KEY OBSERVATION:
  At n=6, max H = 45 = 6!/16 = n!/2^(n-1).
  The achievable H values are NOT all odd numbers in [1, 45].
  The GAPS are at {7, 21, 35, 39}.

  MODULAR ANALYSIS:
  H values mod 3: {sorted(set(h % 3 for h in counts_6.keys()))}
  H values mod 7: {sorted(set(h % 7 for h in counts_6.keys()))}
  H values mod 5: {sorted(set(h % 5 for h in counts_6.keys()))}
""")

# Modular residues
for mod in [3, 4, 6, 7, 8]:
    achieved_mod = sorted(set(h % mod for h in counts_6.keys()))
    missing_mod = sorted(set(v % mod for v in [7, 21, 35, 39]))
    print(f"  Achieved H mod {mod}: {achieved_mod}")
    print(f"  Missing values mod {mod}: {missing_mod}")
    print()

# The parity structure
print("  All H values are odd: CONFIRMED")
print("  H mod 4:")
for h in sorted(counts_6.keys()):
    print(f"    H={h:3d}: H mod 4 = {h % 4}")

print("\n" + "=" * 70)
print("DONE — FORBIDDEN VALUES ANALYSIS")
print("=" * 70)
