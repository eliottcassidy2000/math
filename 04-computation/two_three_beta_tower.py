"""
two_three_beta_tower.py -- kind-pasteur-2026-03-13-S62

The beta basis for the 2-3 bridge.

DISCOVERED: The change of basis x = y - 1 transforms I(Omega, x) into
a 3-adic expansion:

  I(Omega, x) = 1 + a1*x + a2*x^2     (alpha basis, x=2 gives H)

  I(Omega, y-1) = b0 + b1*y + b2*y^2   (beta basis, y=3 gives H)

where:
  b0 = I(Omega, -1) = 1 - a1 + a2     (TOPOLOGICAL: Euler char + 1)
  b1 = a1 - 2*a2                       (DERIVATIVE: I'(Omega, -1))
  b2 = a2                              (CURVATURE: I''(Omega,-1)/2)

The 3-adic tower:
  H mod 3  = b0 mod 3
  H mod 9  = (b0 + 3*b1) mod 9
  H mod 27 = (b0 + 3*b1 + 9*b2) mod 27 = H mod 27 (exact at n<=7)

The 2-adic tower (in alpha basis):
  H mod 2  = 1 (Redei)
  H mod 4  = 1 + 2*(a1 mod 2)
  H mod 8  = 1 + 2*a1 + 4*(a2 mod 2) mod 8

Combined via CRT:
  H mod 6  = CRT(H mod 2, H mod 3)
  H mod 12 = CRT(H mod 4, H mod 3)
  H mod 36 = CRT(H mod 4, H mod 9)

This script explores these towers systematically at n=7.
"""

import numpy as np
from itertools import combinations
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

def count_ham_cycles_exact(A_sub, k):
    if k < 3:
        return 0
    dp = {}
    dp[(1, 0)] = 1
    for mask_size in range(2, k+1):
        for mask in range(1 << k):
            if bin(mask).count('1') != mask_size:
                continue
            if not (mask & 1):
                continue
            for v in range(k):
                if not (mask & (1 << v)):
                    continue
                if v == 0:
                    continue
                prev_mask = mask ^ (1 << v)
                total = 0
                for u in range(k):
                    if (prev_mask & (1 << u)) and A_sub[u][v]:
                        total += dp.get((prev_mask, u), 0)
                if total:
                    dp[(mask, v)] = total
    full = (1 << k) - 1
    total_cycles = 0
    for v in range(1, k):
        if A_sub[v][0] and dp.get((full, v), 0):
            total_cycles += dp[(full, v)]
    return total_cycles

def get_alpha_1_2(A, n):
    """Compute alpha_1 and alpha_2 with EXACT cycle counting."""
    cycles = []
    for size in range(3, n+1, 2):
        for combo in combinations(range(n), size):
            verts = list(combo)
            sub = np.zeros((size, size), dtype=int)
            for a in range(size):
                for b in range(size):
                    sub[a][b] = A[verts[a]][verts[b]]
            if size <= 5:
                c = int(np.trace(np.linalg.matrix_power(sub, size))) // size
            else:
                c = count_ham_cycles_exact(sub, size)
            for _ in range(c):
                cycles.append(frozenset(combo))

    alpha_1 = len(cycles)
    alpha_2 = 0
    for i in range(len(cycles)):
        for j in range(i+1, len(cycles)):
            if not (cycles[i] & cycles[j]):
                alpha_2 += 1

    return alpha_1, alpha_2

# ============================================================
# PART 1: The beta basis at n=7
# ============================================================
print("=" * 70)
print("PART 1: Beta basis at n=7")
print("=" * 70)

n = 7
total_bits = n*(n-1)//2
np.random.seed(42)

results = []
for trial in range(400):
    bits = np.random.randint(0, 1 << total_bits)
    A = bits_to_adj(bits, n)
    H = count_ham_paths(A, n)
    a1, a2 = get_alpha_1_2(A, n)

    b0 = 1 - a1 + a2
    b1 = a1 - 2*a2
    b2 = a2

    results.append({
        'H': H, 'a1': a1, 'a2': a2,
        'b0': b0, 'b1': b1, 'b2': b2, 'bits': bits
    })

# Verify H = b0 + 3*b1 + 9*b2
ok = sum(1 for r in results if r['H'] == r['b0'] + 3*r['b1'] + 9*r['b2'])
print(f"  H = b0 + 3*b1 + 9*b2: {ok}/{len(results)}")

# Beta statistics
print(f"\n  Beta coefficient ranges:")
print(f"    b0 = I(Omega,-1): [{min(r['b0'] for r in results)}, {max(r['b0'] for r in results)}]")
print(f"    b1 = a1 - 2*a2:  [{min(r['b1'] for r in results)}, {max(r['b1'] for r in results)}]")
print(f"    b2 = a2:         [{min(r['b2'] for r in results)}, {max(r['b2'] for r in results)}]")

# ============================================================
# PART 2: The 3-adic tower
# ============================================================
print("\n" + "=" * 70)
print("PART 2: The 3-adic tower")
print("=" * 70)

# Level 1: H mod 3 = b0 mod 3
ok3 = sum(1 for r in results if r['H'] % 3 == r['b0'] % 3)
print(f"  H mod 3 = b0 mod 3: {ok3}/{len(results)}")

# Level 2: H mod 9 = (b0 + 3*b1) mod 9
ok9 = sum(1 for r in results if r['H'] % 9 == (r['b0'] + 3*r['b1']) % 9)
print(f"  H mod 9 = (b0 + 3*b1) mod 9: {ok9}/{len(results)}")

# Level 3: H mod 27 = (b0 + 3*b1 + 9*b2) mod 27
ok27 = sum(1 for r in results if r['H'] % 27 == (r['b0'] + 3*r['b1'] + 9*r['b2']) % 27)
print(f"  H mod 27 = H mod 27: {ok27}/{len(results)} (tautological at n<=7)")

# Distribution of b0 mod 3
b0_mod3 = Counter(r['b0'] % 3 for r in results)
print(f"\n  b0 mod 3: {dict(sorted(b0_mod3.items()))}")

# Distribution of b1 mod 3
b1_mod3 = Counter(r['b1'] % 3 for r in results)
print(f"  b1 mod 3: {dict(sorted(b1_mod3.items()))}")

# H mod 9 distribution
h_mod9 = Counter(r['H'] % 9 for r in results)
print(f"  H mod 9: {dict(sorted(h_mod9.items()))}")

# ============================================================
# PART 3: The 2-adic tower
# ============================================================
print("\n" + "=" * 70)
print("PART 3: The 2-adic tower")
print("=" * 70)

# Level 1: H mod 2 = 1
ok2 = sum(1 for r in results if r['H'] % 2 == 1)
print(f"  H mod 2 = 1: {ok2}/{len(results)}")

# Level 2: H mod 4 = 1 + 2*(a1 mod 2)
ok4 = sum(1 for r in results if r['H'] % 4 == 1 + 2*(r['a1'] % 2))
print(f"  H mod 4 = 1 + 2*(a1 mod 2): {ok4}/{len(results)}")

# Level 3: H mod 8
# H = 1 + 2*a1 + 4*a2
# H mod 8 = (1 + 2*a1 + 4*a2) mod 8
# = (1 + 2*(a1 mod 4) + 4*(a2 mod 2)) mod 8
ok8 = sum(1 for r in results
          if r['H'] % 8 == (1 + 2*r['a1'] + 4*r['a2']) % 8)
print(f"  H mod 8 = (1+2*a1+4*a2) mod 8: {ok8}/{len(results)} (tautological)")

# More interesting: what determines H mod 8?
# Need a1 mod 4 and a2 mod 2
h_mod8 = Counter(r['H'] % 8 for r in results)
print(f"\n  H mod 8 distribution: {dict(sorted(h_mod8.items()))}")
print(f"  H mod 8 possible values: {sorted(h_mod8.keys())}")
print(f"  All odd: {all(v % 2 == 1 for v in h_mod8.keys())}")

# ============================================================
# PART 4: The CRT tower
# ============================================================
print("\n" + "=" * 70)
print("PART 4: CRT combinations")
print("=" * 70)

# H mod 6 = CRT(H mod 2, H mod 3)
h_mod6 = Counter(r['H'] % 6 for r in results)
print(f"  H mod 6: {dict(sorted(h_mod6.items()))}")

# H mod 12 = CRT(H mod 4, H mod 3)
h_mod12 = Counter(r['H'] % 12 for r in results)
print(f"  H mod 12: {dict(sorted(h_mod12.items()))}")

# H mod 36 = CRT(H mod 4, H mod 9)
h_mod36 = Counter(r['H'] % 36 for r in results)
print(f"  H mod 36: {dict(sorted(h_mod36.items()))}")

# H mod 24 = CRT(H mod 8, H mod 3)
h_mod24 = Counter(r['H'] % 24 for r in results)
print(f"  H mod 24: {dict(sorted(h_mod24.items()))}")

# H mod 72 = CRT(H mod 8, H mod 9)
h_mod72 = Counter(r['H'] % 72 for r in results)
print(f"  H mod 72: {dict(sorted(h_mod72.items()))}")

# ============================================================
# PART 5: What determines b0 (the topological content)?
# ============================================================
print("\n" + "=" * 70)
print("PART 5: The topological content b0 = I(Omega, -1)")
print("=" * 70)

# b0 = 1 - a1 + a2
# At n=7: a1 ranges ~2-61, a2 ranges ~0-12
# b0 = 1 - a1 + a2 can be VERY NEGATIVE

b0_dist = Counter(r['b0'] for r in results)
print(f"  b0 range: [{min(r['b0'] for r in results)}, {max(r['b0'] for r in results)}]")
print(f"  b0 mean: {np.mean([r['b0'] for r in results]):.1f}")

# b0 is MOSTLY NEGATIVE (since a1 >> a2 in general)
neg = sum(1 for r in results if r['b0'] < 0)
zero = sum(1 for r in results if r['b0'] == 0)
pos = sum(1 for r in results if r['b0'] > 0)
print(f"  b0 < 0: {neg}, b0 = 0: {zero}, b0 > 0: {pos}")

# b0 mod 3 = H mod 3
# b0 mod 2 = (1 - a1 + a2) mod 2 = (1 + a1 + a2) mod 2
# (since -1 = 1 mod 2)
# = 1 + a1 + a2 mod 2
b0_mod2 = Counter(r['b0'] % 2 for r in results)
print(f"\n  b0 mod 2: {dict(sorted(b0_mod2.items()))}")

# Compare with a1 + a2 parity
a1_a2_parity = Counter((r['a1'] + r['a2']) % 2 for r in results)
print(f"  (a1+a2) mod 2: {dict(sorted(a1_a2_parity.items()))}")
ok_b0_mod2 = sum(1 for r in results if r['b0'] % 2 == (1 + r['a1'] + r['a2']) % 2)
print(f"  b0 mod 2 = 1 + a1 + a2 mod 2: {ok_b0_mod2}/{len(results)}")

# ============================================================
# PART 6: The derivative b1 and its sign
# ============================================================
print("\n" + "=" * 70)
print("PART 6: The derivative b1 = a1 - 2*a2")
print("=" * 70)

# b1 represents I'(Omega, -1) = derivative at x=-1
# If Omega has independent structure, b1 tells us about the "slope" of I

b1_vals = [r['b1'] for r in results]
print(f"  b1 range: [{min(b1_vals)}, {max(b1_vals)}]")
print(f"  b1 mean: {np.mean(b1_vals):.1f}")

b1_neg = sum(1 for r in results if r['b1'] < 0)
b1_zero = sum(1 for r in results if r['b1'] == 0)
b1_pos = sum(1 for r in results if r['b1'] > 0)
print(f"  b1 < 0: {b1_neg}, b1 = 0: {b1_zero}, b1 > 0: {b1_pos}")

# b1 = a1 - 2*a2
# b1 > 0 iff a1 > 2*a2 (i.e., more cycles than twice the disjoint pairs)
# This is ALWAYS true at n=5 (since a2=0, b1=a1>=0)
# At n=7: b1 can be negative? Let's check
print(f"\n  b1 < 0 cases (a1 < 2*a2):")
neg_cases = [r for r in results if r['b1'] < 0]
for r in neg_cases[:5]:
    print(f"    a1={r['a1']}, a2={r['a2']}, b1={r['b1']}, H={r['H']}")

# ============================================================
# PART 7: The ratio b1/b0 and convexity
# ============================================================
print("\n" + "=" * 70)
print("PART 7: Convexity structure")
print("=" * 70)

# I(Omega, x) is a polynomial with positive leading coefficient
# Its graph passes through:
#   (x=-1, b0), (x=0, 1), (x=1, 1+a1), (x=2, H)
# The "convexity" at x=-1 is b2 = a2 >= 0 (ALWAYS)

# Newton's inequality: alpha_k^2 >= alpha_{k-1} * alpha_{k+1}
# At n=7: alpha_0 * alpha_2 <= alpha_1^2 => a2 <= a1^2
newton_ok = sum(1 for r in results if r['a2'] <= r['a1']**2)
print(f"  Newton inequality a2 <= a1^2: {newton_ok}/{len(results)}")

# Also: alpha_1^2 >= alpha_0 * alpha_2 => a1^2 >= a2
# This is equivalent to a2 <= a1^2, same thing
# And: 4*alpha_0*alpha_2 <= alpha_1^2 (sufficient for real roots)
unimodal_ok = sum(1 for r in results if 4*r['a2'] <= r['a1']**2)
print(f"  4*a2 <= a1^2 (real roots sufficient): {unimodal_ok}/{len(results)}")

# The ratio a2/a1^2 (Newton ratio)
newton_ratios = [r['a2'] / r['a1']**2 if r['a1'] > 0 else 0 for r in results]
print(f"  Newton ratio a2/a1^2: [{min(newton_ratios):.4f}, {max(newton_ratios):.4f}], mean={np.mean(newton_ratios):.4f}")

# ============================================================
# PART 8: Connection to sigma matrix
# ============================================================
print("\n" + "=" * 70)
print("PART 8: b0, b1, b2 in terms of sigma matrix")
print("=" * 70)

# sigma(u,v) counts common successors + common predecessors
# Let's compute c3, c5, c7 and see if b0, b1, b2 have sigma expressions

np.random.seed(42)
sigma_data = []

for trial in range(100):
    bits = np.random.randint(0, 1 << total_bits)
    A = bits_to_adj(bits, n)

    # sigma matrix
    sigma = np.zeros((n, n), dtype=int)
    for u in range(n):
        for v in range(u+1, n):
            s = 0
            for w in range(n):
                if w == u or w == v:
                    continue
                if (A[u][w] and A[v][w]) or (A[w][u] and A[w][v]):
                    s += 1
            sigma[u][v] = sigma[v][u] = s

    tr2 = int(np.trace(sigma @ sigma))
    tr3 = int(np.trace(sigma @ sigma @ sigma))

    # Score sequence
    scores = sorted([sum(A[i]) for i in range(n)])

    # c3, c5 (exact)
    c3 = int(np.trace(np.linalg.matrix_power(A, 3))) // 3
    c5 = 0
    for combo in combinations(range(n), 5):
        sub = np.zeros((5, 5), dtype=int)
        verts = list(combo)
        for a in range(5):
            for b in range(5):
                sub[a][b] = A[verts[a]][verts[b]]
        c5 += int(np.trace(np.linalg.matrix_power(sub, 5))) // 5

    # c7 (exact)
    c7 = count_ham_cycles_exact(A, n)

    sigma_data.append({
        'c3': c3, 'c5': c5, 'c7': c7,
        'tr2': tr2, 'tr3': tr3,
        'scores': tuple(scores)
    })

# Check: is tr(Sigma^2) determined by score sequence?
score_to_tr2 = defaultdict(set)
for d in sigma_data:
    score_to_tr2[d['scores']].add(d['tr2'])

all_determined = all(len(v) == 1 for v in score_to_tr2.values())
print(f"  tr(Sigma^2) determined by score sequence: {all_determined}")

# Relationship between c3 and tr(Sigma^2)?
# Known: sum sigma(u,v) = c3-related
# sigma(u,v) = |{w: both u,v succeed w}| + |{w: both u,v are succeeded by w}|
# = number of common out-neighbors + common in-neighbors

# Known identity: sum_{u<v} sigma(u,v) = 3*c3 + some correction?
# Actually sum_{u<v} sigma(u,v) = sum_w (C(d_w^+, 2) + C(d_w^-, 2))
# where d_w^+ = out-degree of w, d_w^- = in-degree of w
# = sum_w C(d_w,2) + C(n-1-d_w, 2) = sum_w d_w(d_w-1)/2 + (n-1-d_w)(n-2-d_w)/2

# And tr(Sigma^2) = 2 * sum_{u<v} sigma(u,v)^2

# Let me check: c3 = (sum_i d_i^2 - sum d_i * (sum d_i / (n-1))) / ... no
# Standard: c3 = (n*(n-1)*(2n-1)/6 - sum s_i^2) / ... hmm
# Actually c3 = sum_{u<v} sigma(u,v) - ... no.

# Let me just check correlations
c3_vals = [d['c3'] for d in sigma_data]
c5_vals = [d['c5'] for d in sigma_data]
c7_vals = [d['c7'] for d in sigma_data]
tr2_vals = [d['tr2'] for d in sigma_data]
tr3_vals = [d['tr3'] for d in sigma_data]

corr_c3_tr2 = np.corrcoef(c3_vals, tr2_vals)[0, 1]
corr_c5_tr2 = np.corrcoef(c5_vals, tr2_vals)[0, 1]
corr_c7_tr2 = np.corrcoef(c7_vals, tr2_vals)[0, 1]
corr_c3_tr3 = np.corrcoef(c3_vals, tr3_vals)[0, 1]

print(f"\n  Correlations:")
print(f"    c3 vs tr(Sigma^2): r = {corr_c3_tr2:.4f}")
print(f"    c5 vs tr(Sigma^2): r = {corr_c5_tr2:.4f}")
print(f"    c7 vs tr(Sigma^2): r = {corr_c7_tr2:.4f}")
print(f"    c3 vs tr(Sigma^3): r = {corr_c3_tr3:.4f}")

# Linear model: c3 from tr2?
# c3 = (1/12) * (sum s_i * (s_i - 1) ... ) = function of score sequence
# tr2 = 2 * sum sigma^2 = also a function of score sequence + more

# ============================================================
# PART 9: The key theorem — mod structure summary
# ============================================================
print("\n" + "=" * 70)
print("PART 9: COMPLETE MOD STRUCTURE SUMMARY")
print("=" * 70)

print("""
AT n=7 (and any n where alpha_k = 0 for k >= 3):

H(T) = I(Omega(T), 2) = 1 + 2*a1 + 4*a2

ALPHA BASIS (centered at x=0):
  I(x) = 1 + a1*x + a2*x^2
  a1 = total # directed odd cycles in T (with multiplicity)
  a2 = # vertex-disjoint pairs of directed odd cycles

BETA BASIS (centered at x=-1, i.e., y = x+1):
  I(y-1) = b0 + b1*y + b2*y^2
  b0 = 1 - a1 + a2 = I(Omega, -1) = reduced Euler char + 1
  b1 = a1 - 2*a2 = I'(Omega, -1) = derivative at topology point
  b2 = a2 = I''(Omega, -1)/2 = curvature at topology point

2-ADIC TOWER (from alpha basis):
  H mod 2 = 1                              [REDEI]
  H mod 4 = 1 + 2*(a1 mod 2)               [cycle count parity]
  H mod 8 = (1 + 2*a1 + 4*a2) mod 8        [requires a1 mod 4, a2 mod 2]
  v_2(H-1) = 1 + v_2(a1 + 2*a2)

3-ADIC TOWER (from beta basis):
  H mod 3 = b0 mod 3 = I(Omega,-1) mod 3   [TOPOLOGY]
  H mod 9 = (b0 + 3*b1) mod 9              [topology + derivative]
  H mod 27 = H mod 27                       [exact at n<=7]

CRT COMBINATIONS:
  H mod 6  = CRT(1, b0 mod 3)              [Redei + topology]
  H mod 12 = CRT(1+2*(a1%2), b0%3)          [cycle parity + topology]
  H mod 36 = CRT(H%4, H%9)                 [requires b1]
  H mod 72 = CRT(H%8, H%9)                 [requires a1%4, a2%2, b1%3]

THE 2-3 BRIDGE:
  x = 2 = -1 + 3
  H = I(Omega, -1 + 3) = b0 + 3*b1 + 9*b2
  The 3 in x = -1 + 3 IS the modular gap between
  the counting point (x=2) and the topological point (x=-1).

INTERPRETATION:
  H is decomposed into three layers:
  1. b0 (TOPOLOGICAL): the Euler characteristic of Ind(Omega)
     - encodes the global structure of cycle conflicts
     - determines H mod 3
  2. b1 (DIFFERENTIAL): the sensitivity of I at the topology point
     - b1 = a1 - 2*a2 measures "excess cycles beyond disjoint pairs"
     - determines H mod 9 (together with b0)
  3. b2 (GEOMETRIC): the curvature / number of disjoint pairs
     - determines the exact value of H (given b0, b1)
     - b2 = a2 >= 0 always (I is convex at x=-1)
""")

# Final verification table
print("  Verification table (first 20 samples):")
print(f"  {'H':>5} {'a1':>4} {'a2':>3} {'b0':>4} {'b1':>4} {'b2':>3} "
      f"{'H%3':>3} {'b0%3':>4} {'H%9':>3} {'(b0+3b1)%9':>11}")
for r in results[:20]:
    b0_3b1_mod9 = (r['b0'] + 3*r['b1']) % 9
    print(f"  {r['H']:>5} {r['a1']:>4} {r['a2']:>3} {r['b0']:>4} {r['b1']:>4} {r['b2']:>3} "
          f"{r['H']%3:>3} {r['b0']%3:>4} {r['H']%9:>3} {b0_3b1_mod9:>11}")

print("\n\nDone.")
