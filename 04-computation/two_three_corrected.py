"""
two_three_corrected.py -- kind-pasteur-2026-03-13-S62

Now that we've confirmed:
  - tr(A^k)/k is EXACT for k=3,5 but OVERCOUNTS for k>=7
  - (H-1)/2 mod 2 = alpha_1 mod 2 with EXACT counts (100/100 at n=7)

Explore the full 2-adic and 3-adic structure of H with CORRECT cycle counts.

Key identity chain:
  H = 1 + 2*alpha_1 + 4*alpha_2 + 8*alpha_3 + ...    (OCF at x=2)
  H mod 2 = 1                                           (Redei)
  (H-1)/2 mod 2 = alpha_1 mod 2                         (just proved)
  H mod 3 = I(Omega, -1) mod 3                          (since 2 = -1 mod 3)

This script explores:
1. Full 2-adic expansion: what determines (H-1)/2^k mod 2 for k=1,2,3?
2. alpha_1 mod 3: does it determine H mod 3?
3. The interplay: (alpha_1 mod 2, alpha_1 mod 3) vs (H mod 2, H mod 3)
4. alpha_2 structure: vertex-disjoint pairs of directed cycles
5. The characteristic polynomial of Omega as a graph
"""

import numpy as np
from itertools import combinations, permutations
from collections import Counter, defaultdict

def bits_to_adj(bits, n):
    A = np.zeros((n, n), dtype=int)
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx):
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    return A

def count_ham_paths(A, n):
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for mask_size in range(2, n+1):
        for mask in range(1 << n):
            if bin(mask).count('1') != mask_size:
                continue
            for v in range(n):
                if not (mask & (1 << v)):
                    continue
                prev_mask = mask ^ (1 << v)
                total = 0
                for u in range(n):
                    if (prev_mask & (1 << u)) and A[u][v]:
                        total += dp.get((prev_mask, u), 0)
                if total:
                    dp[(mask, v)] = total
    full = (1 << n) - 1
    return sum(dp.get((full, v), 0) for v in range(n))

def count_ham_cycles_exact(A, n):
    """Count directed Hamiltonian cycles using Held-Karp DP.
    Fix starting vertex 0, count paths 0->...->v->0."""
    if n < 3:
        return 0
    dp = {}
    dp[(1, 0)] = 1
    for mask_size in range(2, n+1):
        for mask in range(1 << n):
            if bin(mask).count('1') != mask_size:
                continue
            if not (mask & 1):
                continue
            for v in range(n):
                if not (mask & (1 << v)):
                    continue
                if v == 0 and mask_size < n:
                    continue
                if v == 0:
                    continue
                prev_mask = mask ^ (1 << v)
                total = 0
                for u in range(n):
                    if (prev_mask & (1 << u)) and A[u][v] and dp.get((prev_mask, u), 0):
                        total += dp.get((prev_mask, u), 0)
                if total:
                    dp[(mask, v)] = total
    full = (1 << n) - 1
    total_cycles = 0
    for v in range(1, n):
        if A[v][0] and dp.get((full, v), 0):
            total_cycles += dp[(full, v)]
    return total_cycles

def build_omega_graph(A, n):
    """Build Omega(T) with directed cycles as vertices (with multiplicity).
    Returns (cycle_list, adjacency_matrix_of_Omega) where each cycle is
    (vertex_set, direction_index)."""
    cycles = []  # list of (frozenset of vertices, cycle_index)

    for size in range(3, n+1, 2):
        for combo in combinations(range(n), size):
            verts = list(combo)
            sub = np.zeros((size, size), dtype=int)
            for a in range(size):
                for b in range(size):
                    sub[a][b] = A[verts[a]][verts[b]]
            if size <= 5:
                # trace is exact for k=3,5
                c = int(np.trace(np.linalg.matrix_power(sub, size))) // size
            else:
                # use exact Hamiltonian cycle counting
                c = count_ham_cycles_exact(sub, size)
            for idx in range(c):
                cycles.append(frozenset(combo))

    # Build conflict graph: two cycles conflict iff they share a vertex
    m = len(cycles)
    adj = np.zeros((m, m), dtype=int)
    for i in range(m):
        for j in range(i+1, m):
            if cycles[i] & cycles[j]:
                adj[i][j] = adj[j][i] = 1

    return cycles, adj

def independence_poly(adj, m):
    """Compute independence polynomial coefficients of graph with m vertices."""
    # Use inclusion-exclusion / brute force for small m
    alpha = [0] * (m + 1)
    for mask in range(1 << m):
        # Check if mask is an independent set
        verts = [i for i in range(m) if mask & (1 << i)]
        is_indep = True
        for i in range(len(verts)):
            for j in range(i+1, len(verts)):
                if adj[verts[i]][verts[j]]:
                    is_indep = False
                    break
            if not is_indep:
                break
        if is_indep:
            alpha[len(verts)] += 1
    return alpha

# ============================================================
# PART 1: Full 2-adic expansion of H at n=5 (exhaustive)
# ============================================================
print("=" * 70)
print("PART 1: Full 2-adic expansion of H at n=5 (exhaustive)")
print("=" * 70)

n = 5
total_bits = n*(n-1)//2

# Track: H, alpha_1, alpha_2, and the 2-adic digits of H
h_data = []

for bits in range(2**total_bits):
    A = bits_to_adj(bits, n)
    H = count_ham_paths(A, n)
    cycles, adj = build_omega_graph(A, n)
    alpha = independence_poly(adj, len(cycles))

    # Verify OCF
    I_at_2 = sum(alpha[k] * (2**k) for k in range(len(alpha)))
    assert I_at_2 == H, f"OCF failed: I={I_at_2}, H={H}"

    h_data.append({
        'bits': bits, 'H': H,
        'alpha': alpha,
        'alpha_1': alpha[1] if len(alpha) > 1 else 0,
        'alpha_2': alpha[2] if len(alpha) > 2 else 0,
        'alpha_3': alpha[3] if len(alpha) > 3 else 0,
    })

print(f"  OCF verified for all {2**total_bits} tournaments at n=5")

# 2-adic digits
print(f"\n  2-adic structure of H:")
for k in range(5):
    bit_k_dist = Counter()
    for d in h_data:
        H = d['H']
        digit = (H >> k) & 1
        bit_k_dist[digit] += 1
    print(f"    bit {k} of H: {dict(sorted(bit_k_dist.items()))}")

# Check: bit 0 always 1 (Redei), bit 1 = alpha_1 mod 2
print(f"\n  Verification of 2-adic digit identities:")
bit1_match = sum(1 for d in h_data if ((d['H'] >> 1) & 1) == (d['alpha_1'] % 2))
print(f"    bit 1 = alpha_1 mod 2: {bit1_match}/{len(h_data)}")

# What determines bit 2?
# (H-1)/2 = alpha_1 + 2*alpha_2 + 4*alpha_3 + ...
# ((H-1)/2 - alpha_1) / 2 = alpha_2 + 2*alpha_3 + ...
# So bit 2 of H = bit 1 of (H-1)/2 = (alpha_1 + 2*alpha_2) bit 1
#  = ((alpha_1 // 2) + alpha_2) mod 2 ... wait

# Let's be more careful.
# H = 1 + 2*alpha_1 + 4*alpha_2 + 8*alpha_3
# H in binary: bit 0 = 1 (always)
#              bit 1 = alpha_1 mod 2 (carry from bit 0 is 0 since bit 0 = 1)
# Wait, this isn't right because there can be carries!
# H = 1 + 2*a1 + 4*a2 + 8*a3
# This is NOT a binary decomposition because a1, a2, a3 can be > 1

print(f"\n  Carry analysis:")
print(f"  H = 1 + 2*a1 + 4*a2 + 8*a3 (a_k can be > 1)")
print(f"  This is NOT a binary decomposition -- carries propagate!")

# Let's check directly
carry_data = Counter()
for d in h_data:
    a1 = d['alpha_1']
    a2 = d['alpha_2']
    a3 = d.get('alpha_3', 0)
    # H = 1 + 2*a1 + 4*a2 + 8*a3
    reconstructed = 1 + 2*a1 + 4*a2 + 8*a3
    assert reconstructed == d['H'], f"Reconstruct failed"

    # Binary representation
    h_bin = format(d['H'], '08b')

    # "Ideal" bits without carries
    b0 = 1
    b1 = a1 % 2
    b2 = a2 % 2
    b3 = a3 % 2

    # Actual bits
    ab0 = d['H'] & 1
    ab1 = (d['H'] >> 1) & 1
    ab2 = (d['H'] >> 2) & 1
    ab3 = (d['H'] >> 3) & 1

    carry_data[(ab1 == b1, ab2 == b2, ab3 == b3)] += 1

print(f"  Carry propagation (bit1_ok, bit2_ok, bit3_ok): {dict(carry_data)}")

# Hmm, let's think about this differently.
# H = 1 + 2*a1 + 4*a2 + 8*a3
# Bit 0 of H: 1 (always, since 2*a1 + 4*a2 + 8*a3 is even)
# Bit 1 of H: a1 mod 2 (since (H-1)/2 = a1 + 2*a2 + 4*a3, bit 0 of that = a1 mod 2)
# Bit 2 of H: bit 1 of (H-1)/2 = bit 1 of (a1 + 2*a2 + 4*a3)
#           = ((a1 + 2*a2 + 4*a3) >> 1) & 1
#           = ((a1 >> 1) + a2 + 2*a3 + carry) & 1
# The carry from bit 0 of (a1 + 2*a2 + ...) is already accounted for by a1 mod 2

# Let's use the recursive extraction
print(f"\n  Recursive 2-adic extraction:")
for d in h_data[:10]:
    H = d['H']
    a1 = d['alpha_1']
    a2 = d['alpha_2']
    a3 = d.get('alpha_3', 0)

    r0 = (H - 1) // 2  # = a1 + 2*a2 + 4*a3
    r1 = (r0 - a1) // 2  # = a2 + 2*a3
    r2 = (r1 - a2) // 2  # = a3

    assert r0 == a1 + 2*a2 + 4*a3
    assert r1 == a2 + 2*a3
    assert r2 == a3

# This means (H-1)/2 = a1 + 2*a2 + 4*a3
# So: (H-1)/2 mod 2 = a1 mod 2  (PROVED)
#     ((H-1)/2 - a1)/2 mod 2 = a2 mod 2
# But we need a1 EXACTLY to extract a2 mod 2...
# Unless there's a mod-2 identity for bit 2

print("  Confirmed: (H-1)/2 = a1 + 2*a2 + 4*a3 exactly at n=5")
print("  So extracting a_k requires knowing all lower a_j exactly")
print("  No purely H-based mod formula for bit 2 without knowing a1")

# ============================================================
# PART 2: alpha_1, alpha_2 joint distribution at n=5
# ============================================================
print("\n" + "=" * 70)
print("PART 2: Joint distribution of (alpha_1, alpha_2) at n=5")
print("=" * 70)

joint = Counter()
for d in h_data:
    joint[(d['alpha_1'], d['alpha_2'])] += 1

print(f"  (alpha_1, alpha_2) distribution:")
for key in sorted(joint.keys()):
    a1, a2 = key
    h_val = 1 + 2*a1 + 4*a2
    print(f"    a1={a1}, a2={a2}: count={joint[key]:>5}, H={h_val}")

# ============================================================
# PART 3: H mod 3 = I(Omega, -1) mod 3 exhaustive check at n=5
# ============================================================
print("\n" + "=" * 70)
print("PART 3: H mod 3 = I(Omega, -1) mod 3 at n=5")
print("=" * 70)

matches = 0
for d in h_data:
    a = d['alpha']
    I_neg1 = sum((-1)**k * a[k] for k in range(len(a)))
    if d['H'] % 3 == I_neg1 % 3:
        matches += 1

print(f"  Matches: {matches}/{len(h_data)}")

# What's the actual distribution of I(Omega, -1)?
i_neg1_dist = Counter()
for d in h_data:
    a = d['alpha']
    I_neg1 = sum((-1)**k * a[k] for k in range(len(a)))
    i_neg1_dist[I_neg1] += 1

print(f"  I(Omega, -1) distribution: {dict(sorted(i_neg1_dist.items()))}")

# ============================================================
# PART 4: H mod 6 = CRT(1, I(Omega,-1) mod 3)
# ============================================================
print("\n" + "=" * 70)
print("PART 4: H mod 6 structure")
print("=" * 70)

h_mod6 = Counter()
for d in h_data:
    h_mod6[d['H'] % 6] += 1

print(f"  H mod 6 distribution: {dict(sorted(h_mod6.items()))}")
print(f"  Possible H mod 6: {sorted(h_mod6.keys())}")
print(f"  Note: H odd => H mod 6 in {{1, 3, 5}}")

# CRT: H = 1 mod 2, H mod 3 = I(Omega,-1) mod 3
# So H mod 6 is determined by I(Omega,-1) mod 3:
#   I(Omega,-1) = 0 mod 3 => H = 3 mod 6
#   I(Omega,-1) = 1 mod 3 => H = 1 mod 6
#   I(Omega,-1) = 2 mod 3 => H = 5 mod 6
crt_check = 0
for d in h_data:
    a = d['alpha']
    I_neg1 = sum((-1)**k * a[k] for k in range(len(a)))
    r3 = I_neg1 % 3
    expected_mod6 = {0: 3, 1: 1, 2: 5}[r3]
    if d['H'] % 6 == expected_mod6:
        crt_check += 1

print(f"  CRT verification: {crt_check}/{len(h_data)}")

# ============================================================
# PART 5: n=7 -- The full 2-adic and 3-adic structure with exact counts
# ============================================================
print("\n" + "=" * 70)
print("PART 5: n=7 -- full structure with EXACT cycle counts")
print("=" * 70)

n = 7
total_bits = n*(n-1)//2
np.random.seed(42)

n_samples = 200
results = []

for trial in range(n_samples):
    bits = np.random.randint(0, 1 << total_bits)
    A = bits_to_adj(bits, n)
    H = count_ham_paths(A, n)
    cycles, adj = build_omega_graph(A, n)
    alpha = independence_poly(adj, len(cycles))

    # Verify OCF
    I_at_2 = sum(alpha[k] * (2**k) for k in range(len(alpha)))
    if I_at_2 != H:
        print(f"  OCF MISMATCH at trial {trial}: I={I_at_2}, H={H}")
        continue

    I_neg1 = sum((-1)**k * alpha[k] for k in range(len(alpha)))

    results.append({
        'H': H, 'alpha': alpha,
        'alpha_1': alpha[1] if len(alpha) > 1 else 0,
        'alpha_2': alpha[2] if len(alpha) > 2 else 0,
        'alpha_3': alpha[3] if len(alpha) > 3 else 0,
        'I_neg1': I_neg1,
    })

    if trial < 5:
        print(f"  Trial {trial}: H={H}, alpha={alpha[:6]}, I(-1)={I_neg1}")

print(f"\n  OCF verified: {len(results)}/{n_samples}")

# 2-adic check
bit1_ok = sum(1 for r in results if ((r['H'] >> 1) & 1) == (r['alpha_1'] % 2))
print(f"  bit 1 of H = alpha_1 mod 2: {bit1_ok}/{len(results)}")

# H mod 3 = I(-1) mod 3
mod3_ok = sum(1 for r in results if r['H'] % 3 == r['I_neg1'] % 3)
print(f"  H mod 3 = I(-1) mod 3: {mod3_ok}/{len(results)}")

# Joint distribution of (alpha_1 mod 2, alpha_1 mod 3)
a1_mod6 = Counter()
for r in results:
    a1_mod6[(r['alpha_1'] % 2, r['alpha_1'] % 3)] += 1
print(f"\n  alpha_1 (mod 2, mod 3) distribution: {dict(sorted(a1_mod6.items()))}")

# H mod 6 distribution
h_mod6 = Counter()
for r in results:
    h_mod6[r['H'] % 6] += 1
print(f"  H mod 6 distribution: {dict(sorted(h_mod6.items()))}")

# alpha_1 range
a1_vals = [r['alpha_1'] for r in results]
a2_vals = [r['alpha_2'] for r in results]
a3_vals = [r['alpha_3'] for r in results]
print(f"\n  alpha_1 range: [{min(a1_vals)}, {max(a1_vals)}], mean={np.mean(a1_vals):.1f}")
print(f"  alpha_2 range: [{min(a2_vals)}, {max(a2_vals)}], mean={np.mean(a2_vals):.1f}")
print(f"  alpha_3 range: [{min(a3_vals)}, {max(a3_vals)}], mean={np.mean(a3_vals):.1f}")

# ============================================================
# PART 6: The "3 evaluations" theorem
# ============================================================
print("\n" + "=" * 70)
print("PART 6: The three key evaluations of I(Omega, x)")
print("=" * 70)

print("""
I(Omega, 0) = 1                (empty set contribution)
I(Omega, 1) = #independent sets in Omega
I(Omega, 2) = H                (OCF -- Hamiltonian paths)
I(Omega, -1) = H mod 3 equiv   (reduced Euler char + 1)
I(Omega, 3) = ?                (what does this count?)

The evaluations at x = 0, 1, 2 are all "counting" things.
At x = -1: topological (Euler characteristic of Ind(Omega)).

KEY: x=2 simultaneously satisfies:
  - x = 0 mod 2 => Redei: H always odd
  - x = -1 mod 3 => topological: H mod 3 = chi_tilde(Ind(Omega)) + 1
  - x = 2 mod 5 => ???
  - x = 2 mod 7 => ???

What about x = 2 mod p for other primes p?
""")

# Check I(Omega, 2) mod 5 at n=5
for p in [5, 7]:
    # H mod p = I(Omega, 2) mod p = I(Omega, 2 mod p) mod p
    # Since polynomial evaluation: I(Omega, x) mod p = I(Omega, x mod p) mod p
    # At x=2: I(Omega, 2) mod p
    # Since 2 mod 5 = 2, 2 mod 7 = 2: nothing special

    # But what about I(Omega, x) mod p for the CHARACTERISTIC of 2?
    # 2^1 = 2, 2^2 = 4 = -1 mod 5, 2^3 = 3 mod 5, 2^4 = 1 mod 5
    # So 2 has order 4 mod 5

    # More interesting: since I(Omega, x) = sum alpha_k x^k,
    # I(Omega, 2) mod p = sum alpha_k * 2^k mod p
    # The periodic structure of 2^k mod p affects which alpha_k contribute
    pass

# Compute I(Omega, x) at several values for n=5
print("  I(Omega, x) at several values for n=5:")
print(f"  {'bits':>5} {'H':>5} {'I(0)':>5} {'I(1)':>5} {'I(2)':>5} {'I(3)':>5} {'I(-1)':>5} {'I(4)':>5}")
for bits in [0, 10, 40, 100, 200, 500, 1023]:
    A = bits_to_adj(bits, 5)
    cycles, adj = build_omega_graph(A, 5)
    alpha = independence_poly(adj, len(cycles))

    I_vals = {}
    for x in [0, 1, 2, 3, -1, 4]:
        I_vals[x] = sum(alpha[k] * (x**k) for k in range(len(alpha)))

    print(f"  {bits:>5} {I_vals[2]:>5} {I_vals[0]:>5} {I_vals[1]:>5} "
          f"{I_vals[2]:>5} {I_vals[3]:>5} {I_vals[-1]:>5} {I_vals[4]:>5}")

# ============================================================
# PART 7: The 3-cycle contribution to H mod 3
# ============================================================
print("\n" + "=" * 70)
print("PART 7: Direct formula for H mod 3 in terms of cycle counts")
print("=" * 70)

n = 5
total_bits = n*(n-1)//2

print(f"  At n=5:")
print(f"  H = 1 + 2*a1 + 4*a2 (a3=0 at n=5)")
print(f"  H mod 3 = (1 - a1 + a2) mod 3 = I(Omega,-1) mod 3")
print(f"")
print(f"  So: H = 0 mod 3 iff a1 - a2 = 1 mod 3")
print(f"      H = 1 mod 3 iff a1 - a2 = 0 mod 3")
print(f"      H = 2 mod 3 iff a1 - a2 = 2 mod 3")

# Verify
tab = Counter()
for d in h_data:
    a1 = d['alpha_1']
    a2 = d['alpha_2']
    diff = (a1 - a2) % 3
    h3 = d['H'] % 3
    tab[(h3, diff)] += 1

print(f"\n  Cross-check (H mod 3, (a1-a2) mod 3):")
for key in sorted(tab.keys()):
    print(f"    {key}: {tab[key]}")

# ============================================================
# PART 8: The "alpha ratio" a2/a1 and connection to graph density
# ============================================================
print("\n" + "=" * 70)
print("PART 8: Omega graph statistics at n=7")
print("=" * 70)

if results:
    omega_sizes = []
    for r in results:
        alpha = r['alpha']
        omega_size = sum(alpha[k] for k in range(1, len(alpha)))  # total cycles
        # Actually alpha_1 = number of vertices in Omega
        omega_sizes.append(r['alpha_1'])

    print(f"  |Omega| (number of directed odd cycles):")
    print(f"    min={min(omega_sizes)}, max={max(omega_sizes)}, mean={np.mean(omega_sizes):.1f}")

    # Max independent set size
    max_indep = []
    for r in results:
        alpha = r['alpha']
        max_k = max(k for k in range(len(alpha)) if alpha[k] > 0)
        max_indep.append(max_k)

    print(f"  Max independent set size in Omega:")
    print(f"    min={min(max_indep)}, max={max(max_indep)}, mean={np.mean(max_indep):.1f}")
    print(f"    distribution: {dict(Counter(max_indep))}")

    # alpha(Omega) = independence number = max_k with alpha_k > 0
    # This is the max number of pairwise vertex-disjoint directed odd cycles
    # At n=7: max possible = 2 (two 3-cycles use 6 vertices, leaving 1)
    #         or 1 seven-cycle
    # So alpha(Omega) <= 2 at n=7

# ============================================================
# PART 9: Does alpha_2 have a mod-2 or mod-3 pattern?
# ============================================================
print("\n" + "=" * 70)
print("PART 9: alpha_2 modular patterns")
print("=" * 70)

a2_mod2 = Counter()
a2_mod3 = Counter()
for r in results:
    a2_mod2[r['alpha_2'] % 2] += 1
    a2_mod3[r['alpha_2'] % 3] += 1

print(f"  alpha_2 mod 2: {dict(sorted(a2_mod2.items()))}")
print(f"  alpha_2 mod 3: {dict(sorted(a2_mod3.items()))}")

# alpha_2 vs H mod 12
h_mod12_by_a2 = defaultdict(Counter)
for r in results:
    h_mod12_by_a2[r['alpha_2'] % 2][r['H'] % 12] += 1

print(f"\n  H mod 12 conditioned on alpha_2 parity:")
for a2p in sorted(h_mod12_by_a2.keys()):
    print(f"    alpha_2={'even' if a2p==0 else 'odd'}: {dict(sorted(h_mod12_by_a2[a2p].items()))}")

# ============================================================
# PART 10: The number 2*3 = 6 in the structure
# ============================================================
print("\n" + "=" * 70)
print("PART 10: The role of 6 = 2*3")
print("=" * 70)

print("""
Summary of 2-3 appearances:
  - H always odd (factor 2: Redei)
  - H mod 3 = I(Omega, -1) mod 3 (factor 3: topology)
  - tr(Sigma^3) always divisible by 6 (proved: symmetric trace identity)
  - tr(Sigma^3) gaps always divisible by 12 = 2*6 = 4*3
  - Worpitzky: k! | W_k (so 2|W_2, 6|W_3)
  - sigma:lambda:delta = 2:1:1 (factor 2 in sigma)
  - CRT: H mod 6 determined by (H mod 2, H mod 3) = (1, I(Omega,-1) mod 3)

The number 6 = 2*3 is the "modular period" of the 2-3 structure.
H mod 6 captures ALL information from both the 2-adic (Redei) and
3-adic (topology) structures at the first level.

The NEXT level would be:
  H mod 12 = CRT(H mod 4, H mod 3)
  H mod 4 determined by: bit 0 = 1, bit 1 = alpha_1 mod 2
  So H mod 4 = 1 + 2*(alpha_1 mod 2) = 1 or 3

  H mod 12 = CRT(H mod 4, H mod 3)
    = CRT(1 + 2*(a1 mod 2), I(Omega,-1) mod 3)
""")

# Verify H mod 4 at n=7
h_mod4 = Counter()
for r in results:
    h_mod4[r['H'] % 4] += 1
print(f"  H mod 4 distribution: {dict(sorted(h_mod4.items()))}")
print(f"  (Should be only 1 and 3)")

# H mod 12
h_mod12 = Counter()
for r in results:
    h_mod12[r['H'] % 12] += 1
print(f"  H mod 12 distribution: {dict(sorted(h_mod12.items()))}")

# Predict H mod 12 from (alpha_1 mod 2, I(-1) mod 3)
pred_ok = 0
for r in results:
    a1p = r['alpha_1'] % 2
    i_neg1_mod3 = r['I_neg1'] % 3

    # H mod 4 = 1 + 2*a1p
    h4 = 1 + 2 * a1p

    # CRT: find x s.t. x = h4 mod 4, x = i_neg1_mod3 mod 3, 0 <= x < 12
    for x in range(12):
        if x % 4 == h4 and x % 3 == i_neg1_mod3:
            if r['H'] % 12 == x:
                pred_ok += 1
            break

print(f"  H mod 12 predicted from (a1 mod 2, I(-1) mod 3): {pred_ok}/{len(results)}")

print("\n\nDone.")
