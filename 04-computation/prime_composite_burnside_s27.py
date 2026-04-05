#!/usr/bin/env python3
"""
opus-2026-04-05-S27: Prime vs composite in the Burnside perturbation expansion.

THE KEY THEORETICAL QUESTION:
In the series a(n) = (2^{C(n,2)}/n!) × [1 + Σ_k R_k(n)],
how do the correction terms R_k differ for prime vs composite n?

For PRIME p:
  - The partition (p) [single p-cycle] exists (p odd prime)
  - Its edge count = (p-1)/2
  - Ratio: (p-1)! × 2^{(p-1)/2 - p(p-1)/2} = (p-1)! × 2^{-(p-1)^2/2}
  - This is NEGLIGIBLE (super-exponentially small in p)
  - BUT: no partition like (p/3)^3 exists (p is prime → no non-trivial divisors)

For COMPOSITE n = ab:
  - Partition (a^b) exists if a is odd
  - Edge count for (a^b): C(b,2)×a + b×(a-1)/2 = b(b-1)a/2 + b(a-1)/2
  - This creates a CLUSTER of correction terms not present at prime n

PREDICTION: The correction terms at composite n should be LARGER than at
nearby primes, because composite n has more "resonant" partitions.

This would manifest as: a(n)/leading_term having bigger corrections at
composite n than at prime n.
"""

from math import factorial, gcd, comb
from sympy import isprime, factorint
import time

def z_lambda(pm):
    z = 1
    for k, m in pm:
        z *= k**m * factorial(m)
    return z

def edge_count(pm):
    total = 0
    for i, (ki, mi) in enumerate(pm):
        total += mi*(ki-1)//2
        total += mi*(mi-1)//2*ki
        for j in range(i+1, len(pm)):
            kj, mj = pm[j]
            total += mi*mj*gcd(ki,kj)
    return total

def enumerate_odd_partitions(n):
    results = []
    parts = list(range(n if n%2 else n-1, 0, -2))
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

# ═══════════════════════════════════════════════════════
# MAIN ANALYSIS: correction term structure
# ═══════════════════════════════════════════════════════

print("="*70)
print("PRIME vs COMPOSITE in BURNSIDE PERTURBATION EXPANSION")
print("="*70)

print("\n── Partition counts and identity fractions ──\n")
print(f"{'n':>3s} {'P/C':>4s} {'#parts':>7s} {'#odd_div':>8s} "
      f"{'id_frac':>12s} {'coupling':>12s} {'ratio_to_prev':>14s}")
print("-"*72)

prev_coupling = None
for n in range(3, 40):
    parts = enumerate_odd_partitions(n)
    nparts = len(parts)
    nfact = factorial(n)
    cn2 = n*(n-1)//2

    # Identity contribution
    identity = 2**cn2

    # Total Burnside sum
    total = 0
    for p in parts:
        z = z_lambda(p)
        e = edge_count(p)
        total += (nfact // z) * (2**e)

    id_frac = identity / total
    coupling = 1 - id_frac

    label = "P" if isprime(n) else "C"

    # Number of odd divisors of n (relevant for "resonant" partitions)
    odd_divs = len([d for d in range(1, n+1, 2) if n % d == 0])

    ratio = coupling / prev_coupling if prev_coupling and prev_coupling > 0 else 0
    prev_coupling = coupling

    print(f"{n:3d} [{label}] {nparts:7d} {odd_divs:8d} "
          f"{id_frac:12.8f} {coupling:12.6e} {ratio:14.4f}")

# ═══════════════════════════════════════════════════════
# DETAILED: Which partition types contribute MOST at each n?
# ═══════════════════════════════════════════════════════

print(f"\n{'='*70}")
print("TOP NON-IDENTITY CORRECTION by n (prime vs composite)")
print("="*70)

for n in range(5, 30):
    parts = enumerate_odd_partitions(n)
    nfact = factorial(n)
    cn2 = n*(n-1)//2
    identity = 2**cn2

    # Compute each partition's contribution relative to identity
    contribs = []
    for p in parts:
        z = z_lambda(p)
        e = edge_count(p)
        w = (nfact // z) * (2**e)
        if p == [(1, n)]:
            continue  # skip identity
        ratio = w / identity
        contribs.append((ratio, p))

    contribs.sort(reverse=True)

    label = "PRIME" if isprime(n) else "comp "

    # Show top 3 corrections
    top3 = contribs[:3]
    top3_str = "; ".join(f"{'+'.join(f'{k}^{m}' for k,m in p):12s} r={r:.4e}" for r, p in top3)
    print(f"  n={n:2d} [{label}]: {top3_str}")

# ═══════════════════════════════════════════════════════
# THE KEY TEST: Composite "resonance" partitions
# ═══════════════════════════════════════════════════════

print(f"\n{'='*70}")
print("RESONANCE PARTITIONS — only exist at composite n")
print("="*70)

print("""
A "resonance" partition is (d^{n/d}) where d|n, d>1, d odd.
These exist ONLY at composite n.
At prime p: the only single-part partition with part size dividing p is (1^p) and (p).
No "resonance" partitions exist.

This means composite n has EXTRA correction terms not available at prime n.
""")

for n in range(3, 30):
    # Find resonance partitions: (d^{n/d}) with d odd, d|n, d>1
    resonances = []
    for d in range(3, n+1, 2):
        if n % d == 0:
            m = n // d
            p = [(d, m)]
            z = z_lambda(p)
            e = edge_count(p)
            nfact = factorial(n)
            cn2 = n*(n-1)//2
            w = (nfact // z) * (2**e)
            ratio = w / (2**cn2)
            resonances.append((d, m, ratio))

    label = "P" if isprime(n) else "C"
    if resonances:
        res_str = ", ".join(f"({d}^{m}): {r:.4e}" for d, m, r in resonances)
        print(f"  n={n:2d} [{label}]: {res_str}")
    else:
        print(f"  n={n:2d} [{label}]: (none)")

# ═══════════════════════════════════════════════════════
# THE PUNCHLINE: excess coupling at composite n
# ═══════════════════════════════════════════════════════

print(f"\n{'='*70}")
print("EXCESS COUPLING at composite n vs nearby primes")
print("="*70)

print("""
Define: excess(n) = coupling(n) / interpolated_coupling(prev_prime, next_prime)
If excess > 1: composite n has MORE correction than expected from prime trend.
""")

couplings = {}
for n in range(3, 40):
    parts = enumerate_odd_partitions(n)
    nfact = factorial(n)
    cn2 = n*(n-1)//2
    identity = 2**cn2
    total = sum((nfact // z_lambda(p)) * (2**edge_count(p)) for p in parts)
    couplings[n] = 1 - identity/total

# For each composite n, interpolate between neighboring primes
primes = [n for n in range(3, 40) if isprime(n)]
for n in range(4, 38):
    if isprime(n):
        continue
    # Find nearest prime before and after
    prev_p = max(p for p in primes if p < n)
    next_p = min(p for p in primes if p > n)

    # Linear interpolation of log(coupling)
    import math
    lc_prev = math.log(couplings[prev_p]) if couplings[prev_p] > 0 else -100
    lc_next = math.log(couplings[next_p]) if couplings[next_p] > 0 else -100
    t = (n - prev_p) / (next_p - prev_p)
    lc_interp = lc_prev + t * (lc_next - lc_prev)
    c_interp = math.exp(lc_interp)

    excess = couplings[n] / c_interp if c_interp > 0 else 0
    factors = dict(factorint(n))
    fstr = "×".join(f"{p}^{e}" if e > 1 else str(p) for p, e in factors.items())
    print(f"  n={n:2d} ({fstr:8s}): coupling={couplings[n]:.6e}, "
          f"interp={c_interp:.6e}, excess={excess:.4f}")

print("\n\nDONE.")
