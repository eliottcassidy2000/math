#!/usr/bin/env python3
"""
opus-2026-04-05-S29: Abstract Burnside Perturbation Theory.

DEEPER INVESTIGATIONS:

1. THE SPECTRAL GAP: Why does the identity dominate? Because the orbit count
   e(λ) has a UNIQUE MAXIMUM at λ=(1^n), and the gap to the next partition
   grows with n. This gap determines the convergence rate.

2. THE EULER PRODUCT: The Burnside sum factors as a product over cycle lengths.
   The perturbation expansion = Taylor expansion of log(product).

3. CONTINUOUS BASE: What happens as we vary the base from 1 to ∞?
   At base=1: all partitions contribute equally → a(n) = p(n) (partition count!)
   At base=∞: only identity matters → a(n) → 1
   Is there a phase transition?

4. THE DEFECT ENERGY SPECTRUM: Each partition type has a "defect energy"
   Δe = e_max - e(λ). The spectrum of defect energies determines everything.

5. UNIVERSAL FORMULAS: The leading correction for ANY base-b Burnside problem
   on C(n,k) structures has the form poly(n)/b^{α(k)·n}.

6. MÖBIUS INVERSION: The perturbation expansion is related to Möbius inversion
   on the partition lattice. What is the Möbius function telling us?
"""

from math import factorial, gcd, comb, log, log2, exp
from collections import defaultdict
import time

def z_lambda(pm):
    z = 1
    for k, m in pm:
        z *= k**m * factorial(m)
    return z

def enumerate_odd_partitions(n):
    results = []
    parts = list(range(n if n%2 else n-1, 0, -2))
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

def enumerate_partitions(n):
    results = []
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

def edge_count_general(pm):
    """General edge count for C(n,2) structures: gcd_sum + floor(k/2) self-term."""
    total = 0
    for i, (ki, mi) in enumerate(pm):
        total += mi * (ki // 2)
        total += mi * (mi - 1) // 2 * ki
        for j in range(i+1, len(pm)):
            kj, mj = pm[j]
            total += mi * mj * gcd(ki, kj)
    return total

# ═══════════════════════════════════════════════════════
# 1. THE SPECTRAL GAP
# ═══════════════════════════════════════════════════════

print("="*70)
print("1. SPECTRAL GAP — why identity dominates")
print("="*70)

print("""
The orbit count e(λ) counts independent edge orbits under permutation type λ.
The identity (1^n) gives e = C(n,2) (every edge is its own orbit).
The SPECTRAL GAP Δ = C(n,2) - e(next best) determines convergence.
""")

for n in range(3, 25):
    # Compute e for all partitions
    parts = enumerate_partitions(n)
    e_vals = []
    for p in parts:
        e = edge_count_general(p)
        e_vals.append((e, p))
    e_vals.sort(reverse=True)

    e_max = e_vals[0][0]
    e_second = e_vals[1][0] if len(e_vals) > 1 else e_max

    gap = e_max - e_second
    second_part = e_vals[1][1]

    # What's the second-best partition type?
    sp_str = "+".join(f"{k}^{m}" for k, m in second_part)

    print(f"  n={n:2d}: e_max={e_max:5d}, gap={gap:3d}, next={sp_str:12s}  "
          f"(gap/n = {gap/n:.2f})")

# ═══════════════════════════════════════════════════════
# 2. CONTINUOUS BASE — the phase diagram
# ═══════════════════════════════════════════════════════

print(f"\n{'='*70}")
print("2. PHASE DIAGRAM — varying the base continuously")
print("="*70)

print("""
a(n; b) = (1/n!) Σ_λ (n!/z_λ) × b^{e(λ)}

At b=1: a(n;1) = (1/n!) Σ_λ n!/z_λ = number of conjugacy classes? No...
Actually at b=1: b^e = 1 for all e, so a(n;1) = (1/n!) Σ_λ n!/z_λ = p(n)? No...
Σ_λ n!/z_λ = Σ_{σ ∈ S_n} 1 = n!
So a(n;1) = n!/n! = 1. Every structure is "the same" when coloring = trivial.

At b=0: b^e = 0 for e>0, = 1 for e=0. Only partitions with e=0 contribute.
e=0 iff... the only partition with e=0 is (n) for n=1. For n≥2, all partitions
have e≥1. So a(n;0) = 0 for n≥2. Trivially.

So the range is b ∈ [0, ∞) with a(n;b) interpolating from 0 to ∞.
""")

n = 10
parts = enumerate_partitions(n)
nfact = factorial(n)

print(f"Phase diagram at n={n}:")
for b_float in [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 10.0]:
    # Use Fraction for exact computation where possible, float otherwise
    total = 0.0
    for p in parts:
        z = z_lambda(p)
        e = edge_count_general(p)
        total += (nfact / z) * (b_float ** e)
    a_val = total / nfact
    # Identity fraction
    id_e = edge_count_general([(1, n)])
    id_contrib = b_float ** id_e
    id_frac = id_contrib / total * nfact if total > 0 else 1

    print(f"  b={b_float:5.1f}: a({n};b) = {a_val:15.2f}, id_frac = {id_frac:.6f}, "
          f"coupling = {1-id_frac:.4e}")

# Critical base: where does the coupling cross 50%?
print(f"\nCritical base b_c where coupling = 50% (n={n}):")
lo, hi = 1.0, 5.0
for _ in range(50):
    mid = (lo + hi) / 2
    total = sum((nfact / z_lambda(p)) * (mid ** edge_count_general(p)) for p in parts)
    id_e = edge_count_general([(1, n)])
    id_frac = (mid ** id_e) / total * nfact
    if id_frac > 0.5:
        hi = mid
    else:
        lo = mid
b_c = (lo + hi) / 2
print(f"  b_c(n={n}) ≈ {b_c:.6f}")

# Find b_c for several n
print(f"\nCritical base b_c vs n:")
for n in range(3, 16):
    parts = enumerate_partitions(n)
    nfact = factorial(n)

    lo, hi = 1.0, 10.0
    for _ in range(50):
        mid = (lo + hi) / 2
        total = sum((nfact / z_lambda(p)) * (mid ** edge_count_general(p)) for p in parts)
        id_e = edge_count_general([(1, n)])
        id_frac = (mid ** id_e) / total * nfact
        if id_frac > 0.5:
            hi = mid
        else:
            lo = mid
    b_c = (lo + hi) / 2
    print(f"  n={n:2d}: b_c = {b_c:.6f}")

# ═══════════════════════════════════════════════════════
# 3. THE DEFECT ENERGY SPECTRUM
# ═══════════════════════════════════════════════════════

print(f"\n{'='*70}")
print("3. DEFECT ENERGY SPECTRUM — structure of corrections")
print("="*70)

n = 15
parts = enumerate_partitions(n)
e_max = edge_count_general([(1, n)])

# Compute defect energy for each partition
defects = []
for p in parts:
    e = edge_count_general(p)
    delta = e_max - e
    z = z_lambda(p)
    nfact = factorial(n)
    weight = nfact // z  # multiplicity
    defects.append((delta, weight, p))

defects.sort()

print(f"Defect energy spectrum at n={n} (e_max = {e_max}):\n")
print(f"  {'Δe':>4s} {'weight':>12s} {'partition':>25s} {'density_of_states':>20s}")
print(f"  {'─'*4} {'─'*12} {'─'*25} {'─'*20}")

# Group by defect energy
by_delta = defaultdict(list)
for d, w, p in defects:
    by_delta[d].append((w, p))

for delta in sorted(by_delta.keys())[:20]:
    items = by_delta[delta]
    total_weight = sum(w for w, p in items)
    n_types = len(items)
    # Show first partition type
    _, first_p = items[0]
    p_str = "+".join(f"{k}^{m}" for k, m in first_p)
    print(f"  {delta:4d} {total_weight:12d} {p_str:>25s} ({n_types} type{'s' if n_types>1 else ''})")

# ═══════════════════════════════════════════════════════
# 4. THE EULER PRODUCT STRUCTURE
# ═══════════════════════════════════════════════════════

print(f"\n{'='*70}")
print("4. EULER PRODUCT — factored perturbation expansion")
print("="*70)

print("""
The Burnside sum factors as:
  Σ_λ (n!/z_λ) b^{e(λ)} = n! × coeff of x^n in F(x)
where F(x) = Π_{k≥1} G_k(x)
and G_k(x) = Σ_{m≥0} b^{e_self(k,m)} x^{km} / (k^m m!)

The INTERACTION between different k values contributes cross-terms.
But the SELF contribution of each k is independent.

This gives a MULTIPLICATIVE structure to the perturbation expansion:
  S(n) = Π_k S_k(n) × cross_corrections

where S_k(n) captures the contribution of cycle length k.
""")

# Compute the self-contribution of each cycle length
n = 12
print(f"Self-contributions at n={n}:")
parts_full = enumerate_partitions(n)
e_max = edge_count_general([(1, n)])

# The contribution of partitions using ONLY parts of size k
for k in range(1, n+1):
    # Partitions of n using only parts of size k: (k^{n/k}) if k|n
    if n % k != 0:
        continue
    m = n // k
    p = [(k, m)]
    e = edge_count_general(p)
    z = z_lambda(p)
    nfact = factorial(n)
    ratio = (nfact / z) * (2**(e - e_max))
    print(f"  k={k:2d}: λ=({k}^{m}), e={e:4d}, Δe={e_max-e:4d}, ratio={ratio:.6e}")

# ═══════════════════════════════════════════════════════
# 5. UNIVERSAL LEADING CORRECTION FORMULA
# ═══════════════════════════════════════════════════════

print(f"\n{'='*70}")
print("5. UNIVERSAL LEADING CORRECTION — depends on smallest non-trivial cycle")
print("="*70)

print("""
For base b and structures on C(n,2) edges:
  Leading correction from transpositions (k=2):
    R₂(n) = C(n,2) × b^{-1} × ...
  Leading correction from 3-cycles (k=3):
    R₃(n) = C(n,3) × (2/3) × b^{-(2n-4)}

General: the smallest allowed cycle determines the convergence rate.
- All parts allowed (graphs, digraphs): k=2 → rate ~ b^{-1} per n → SLOW
- Odd parts only (tournaments): k=3 → rate ~ b^{-2} per n → FASTER
""")

# Compute: for each allowed-part restriction, what's the convergence rate?
print("Convergence rate (ratio of coupling constants at successive n):\n")
print(f"  {'n':>3s}  {'all parts':>12s}  {'odd only':>12s}  {'k≥3':>12s}  {'k≥5':>12s}")

prev = {}
for n in range(4, 20):
    ratios = {}
    for label, odd_only, min_k in [('all', False, 1), ('odd', True, 1),
                                     ('k≥3', False, 3), ('k≥5', False, 5)]:
        if odd_only:
            parts = enumerate_odd_partitions(n)
        else:
            parts = [p for p in enumerate_partitions(n)
                     if all(k >= min_k or k == 1 for k, m in p)]
            # Actually need to enumerate only partitions with parts ≥ min_k or = 1
            # Simpler: just use all partitions and filter
            parts_raw = enumerate_partitions(n)
            parts = [p for p in parts_raw if all(k >= min_k or k == 1 for k, m in p)]

        nfact = factorial(n)
        e_max = edge_count_general([(1, n)])
        total = sum((nfact // z_lambda(p)) * (2**edge_count_general(p)) for p in parts)
        coupling = 1 - (2**e_max) / total

        if label in prev and prev[label] > 0:
            ratios[label] = coupling / prev[label]
        else:
            ratios[label] = 0
        prev[label] = coupling

    if all(v > 0 for v in ratios.values()):
        print(f"  {n:3d}  {ratios['all']:12.6f}  {ratios['odd']:12.6f}  "
              f"{ratios.get('k≥3',0):12.6f}  {ratios.get('k≥5',0):12.6f}")

# ═══════════════════════════════════════════════════════
# 6. THE DEEP ABSTRACTION: DEFECT = REPRESENTATION
# ═══════════════════════════════════════════════════════

print(f"\n{'='*70}")
print("6. DEEP ABSTRACTION — each correction = an irrep of S_n")
print("="*70)

print("""
In the Burnside sum, each conjugacy class λ of S_n contributes:
  (n!/z_λ) × b^{e(λ)}

The conjugacy classes of S_n are indexed by partitions of n.
The representation theory of S_n associates to each partition λ:
  - An irreducible representation ρ_λ of dimension d_λ = n!/z_λ × (something)
  - The character χ_λ evaluated at various conjugacy classes

The Burnside sum is NOT the character sum — it's a DIFFERENT weighting.
But the PARTITION INDEXING is the same!

Key insight: the "defect energy" Δe(λ) = C(n,2) - e(λ) can be expressed as
a function of the partition λ. And the convergence rate is determined by the
MINIMAL defect, which corresponds to the SMALLEST non-trivial cycle.

For S_n representation theory:
- The trivial representation (1^n) → identity partition → Δe = 0
- The sign representation (n) → full cycle → Δe = (n-1)(n-2)/2 (LARGE)
- The standard representation → transposition (2,1^{n-2}) → Δe = 1

So the perturbation expansion is ORDERED by representation-theoretic complexity,
from trivial (identity) to most complex (full cycle).
""")

# Verify: is Δe related to any standard representation-theoretic quantity?
print("Defect energy vs representation-theoretic quantities:\n")
print(f"  {'partition':>15s} {'Δe':>5s} {'dim(ρ_λ)':>12s} {'Δe/dim':>8s}")

n = 8
parts = enumerate_partitions(n)
e_max = edge_count_general([(1, n)])

# Hook length formula for dimension of irrep
def hook_lengths(partition):
    """Compute dimension of irrep from hook length formula.
    partition = list of parts in decreasing order."""
    # Young diagram
    p = sorted([k for k, m in partition for _ in range(m)], reverse=True)
    n = sum(p)
    # Conjugate partition
    conj = []
    for j in range(p[0]):
        conj.append(sum(1 for pi in p if pi > j))
    # Hook lengths
    hooks = []
    for i in range(len(p)):
        for j in range(p[i]):
            hook = p[i] - j + conj[j] - i - 1
            hooks.append(hook)
    dim = factorial(n)
    for h in hooks:
        dim //= h
    return dim

for p in parts[:12]:
    e = edge_count_general(p)
    delta = e_max - e
    try:
        dim = hook_lengths(p)
    except Exception:
        dim = -1
    ratio = delta / dim if dim > 0 else 0
    p_str = "+".join(f"{k}^{m}" for k, m in p)
    print(f"  {p_str:>15s} {delta:5d} {dim:12d} {ratio:8.4f}")

# ═══════════════════════════════════════════════════════
# 7. THE GENERATING FUNCTION
# ═══════════════════════════════════════════════════════

print(f"\n{'='*70}")
print("7. GENERATING FUNCTION — the partition function as formal series")
print("="*70)

print("""
The exponential generating function of a(n; b):
  F(x; b) = Σ_{n≥0} a(n; b) × x^n

For simple graphs (A000088 at b=2):
  F(x; 2) = 1 + x + 2x² + 4x³/3! + 11x⁴/4! + ...

The Burnside structure gives:
  F(x; b) = Σ_{n≥0} (1/n!) Σ_{λ⊢n} (n!/z_λ) b^{e(λ)} × x^n
           = Σ_λ (b^{e(λ)}/z_λ) × x^{|λ|}

This factors as:
  F(x; b) = Π_{k≥1} [Σ_{m≥0} b^{s(k,m)} × (x^k)^m / (k^m × m!)]

where s(k,m) = m(k-1)/2 + C(m,2)k = self-edge contribution of (k^m).

The CROSS-EDGE contributions couple different factors, so the product
is not exact. But the perturbation expansion captures the corrections.
""")

# Compute the "bare" product (no cross-terms) vs exact
for n in [5, 8, 10]:
    parts = enumerate_partitions(n)
    nfact = factorial(n)
    e_max = edge_count_general([(1, n)])

    # Exact
    total_exact = sum((nfact // z_lambda(p)) * (2**edge_count_general(p)) for p in parts)
    a_exact = total_exact // nfact

    # "Bare" product: only self-contributions (ignore cross-terms)
    # This means: for partition (k1^m1, k2^m2, ...), use
    # e_bare = Σ_i [m_i(k_i-1)/2 + C(m_i,2)k_i + m_i(m_i)k_i/2... no]
    # Actually, e_bare = e_full - Σ_{i<j} m_i m_j gcd(k_i, k_j)
    total_bare = 0
    for p in parts:
        e_full = edge_count_general(p)
        # Subtract cross-terms
        cross = 0
        for i in range(len(p)):
            for j in range(i+1, len(p)):
                cross += p[i][1] * p[j][1] * gcd(p[i][0], p[j][0])
        e_bare = e_full - cross
        total_bare += (nfact // z_lambda(p)) * (2**e_bare)

    a_bare = total_bare // nfact
    ratio = a_bare / a_exact if a_exact > 0 else 0
    print(f"  n={n:2d}: a_exact={a_exact}, a_bare(no cross)={a_bare}, ratio={ratio:.6f}")

print("\n\nDONE.")
