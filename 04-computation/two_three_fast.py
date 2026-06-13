"""
two_three_fast.py -- kind-pasteur-2026-03-13-S62

FAST version of the 2-3 exploration.
At n=7, |Omega| can be 30-50 cycles, making brute-force IP(Omega) infeasible.
But we only need alpha_1 (total cycles) and alpha_2 (disjoint pairs) since
the max independent set size at n=7 is at most 2 (two disjoint 3-cycles
use 6 vertices, leaving 1).

Wait -- can we have alpha_3 > 0 at n=7?
  - Three mutually disjoint 3-cycles need 9 vertices. n=7 < 9. NO.
  - A 3-cycle and a disjoint 5-cycle need 8 > 7. NO.
  - A 7-cycle is disjoint from nothing (uses all vertices). NO.
So alpha_k = 0 for k >= 3 at n=7.

Therefore: H = 1 + 2*alpha_1 + 4*alpha_2 EXACTLY at n=7.
And: I(Omega, -1) = 1 - alpha_1 + alpha_2 EXACTLY at n=7.

This makes everything fast: just count cycles and disjoint pairs.
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
    """Held-Karp for directed Ham cycles on k-vertex subgraph."""
    if k < 3:
        return 0
    dp = {}
    dp[(1, 0)] = 1
    for mask in range(1, 1 << k):
        if not (mask & 1):
            continue
        for v in range(k):
            if not (mask & (1 << v)):
                continue
            if v == 0:
                continue
            key = (mask, v)
            if key not in dp:
                continue
            cnt = dp[key]
            for w in range(k):
                if mask & (1 << w):
                    continue
                if A_sub[v][w]:
                    nkey = (mask | (1 << w), w)
                    dp[nkey] = dp.get(nkey, 0) + cnt
    # Also build forward from start
    dp2 = {}
    dp2[(1, 0)] = 1
    for mask_size in range(2, k+1):
        for mask in range(1 << k):
            if bin(mask).count('1') != mask_size:
                continue
            if not (mask & 1):
                continue
            for v in range(k):
                if not (mask & (1 << v)):
                    continue
                if v == 0 and mask_size < k:
                    continue
                if v == 0:
                    continue
                prev_mask = mask ^ (1 << v)
                total = 0
                for u in range(k):
                    if (prev_mask & (1 << u)) and A_sub[u][v]:
                        total += dp2.get((prev_mask, u), 0)
                if total:
                    dp2[(mask, v)] = total
    full = (1 << k) - 1
    total_cycles = 0
    for v in range(1, k):
        if A_sub[v][0] and dp2.get((full, v), 0):
            total_cycles += dp2[(full, v)]
    return total_cycles

def get_directed_cycles(A, n):
    """Get all directed odd cycles as list of (frozenset_of_vertices, multiplicity).
    Returns flat list where each directed cycle is its own entry."""
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
    return cycles

def compute_alpha_1_2(cycles):
    """Compute alpha_1 and alpha_2 from cycle list."""
    alpha_1 = len(cycles)
    alpha_2 = 0
    for i in range(len(cycles)):
        for j in range(i+1, len(cycles)):
            if not (cycles[i] & cycles[j]):
                alpha_2 += 1
    return alpha_1, alpha_2

# ============================================================
# PART 1: Verify H = 1 + 2*a1 + 4*a2 at n=5 and n=7
# ============================================================
print("=" * 70)
print("PART 1: H = 1 + 2*a1 + 4*a2 verification")
print("=" * 70)

for n in [5, 7]:
    total_bits = n*(n-1)//2
    np.random.seed(42)

    n_samples = min(2**total_bits, 500) if n <= 5 else 300

    ok = 0
    for trial in range(n_samples):
        if n <= 5 and trial < 2**total_bits:
            bits = trial
        else:
            bits = np.random.randint(0, 1 << total_bits)

        A = bits_to_adj(bits, n)
        H = count_ham_paths(A, n)
        cycles = get_directed_cycles(A, n)
        a1, a2 = compute_alpha_1_2(cycles)

        reconstructed = 1 + 2*a1 + 4*a2
        if reconstructed == H:
            ok += 1
        elif trial < 5:
            print(f"  MISMATCH n={n} trial={trial}: H={H}, 1+2*{a1}+4*{a2}={reconstructed}")

    print(f"  n={n}: H = 1 + 2*a1 + 4*a2 verified: {ok}/{n_samples}")

# ============================================================
# PART 2: Full 2-adic and 3-adic structure at n=5 (exhaustive)
# ============================================================
print("\n" + "=" * 70)
print("PART 2: Full structure at n=5 (exhaustive)")
print("=" * 70)

n = 5
total_bits = n*(n-1)//2

h_data = []
for bits in range(2**total_bits):
    A = bits_to_adj(bits, n)
    H = count_ham_paths(A, n)
    cycles = get_directed_cycles(A, n)
    a1, a2 = compute_alpha_1_2(cycles)
    h_data.append({'bits': bits, 'H': H, 'a1': a1, 'a2': a2})

# Joint distribution
joint = Counter()
for d in h_data:
    joint[(d['a1'], d['a2'])] += 1

print(f"  (alpha_1, alpha_2) distribution:")
print(f"  {'a1':>4} {'a2':>4} {'count':>6} {'H':>5} {'H%3':>4} {'H%6':>4} {'H%12':>5}")
for key in sorted(joint.keys()):
    a1, a2 = key
    h_val = 1 + 2*a1 + 4*a2
    print(f"  {a1:>4} {a2:>4} {joint[key]:>6} {h_val:>5} {h_val%3:>4} {h_val%6:>4} {h_val%12:>5}")

# 2-adic digits
print(f"\n  2-adic structure:")
print(f"  bit 0 of H: always 1 (Redei)")
bit1 = Counter((d['H'] >> 1) & 1 for d in h_data)
print(f"  bit 1 of H: {dict(bit1)} (= alpha_1 mod 2)")
bit2 = Counter((d['H'] >> 2) & 1 for d in h_data)
print(f"  bit 2 of H: {dict(bit2)}")
bit3 = Counter((d['H'] >> 3) & 1 for d in h_data)
print(f"  bit 3 of H: {dict(bit3)}")

# H mod 3, 6, 12
print(f"\n  H mod 3: {dict(Counter(d['H'] % 3 for d in h_data))}")
print(f"  H mod 6: {dict(Counter(d['H'] % 6 for d in h_data))}")

# I(Omega, -1) = 1 - a1 + a2
print(f"\n  I(Omega, -1) = 1 - a1 + a2:")
i_neg1_dist = Counter()
for d in h_data:
    i_neg1 = 1 - d['a1'] + d['a2']
    i_neg1_dist[i_neg1] += 1
print(f"  distribution: {dict(sorted(i_neg1_dist.items()))}")

# Verify H mod 3 = I(Omega, -1) mod 3
ok = sum(1 for d in h_data if d['H'] % 3 == (1 - d['a1'] + d['a2']) % 3)
print(f"  H mod 3 = I(Omega, -1) mod 3: {ok}/{len(h_data)}")

# ============================================================
# PART 3: n=7 — larger sample with 2-adic and 3-adic analysis
# ============================================================
print("\n" + "=" * 70)
print("PART 3: n=7 — 2-adic and 3-adic structure (300 samples)")
print("=" * 70)

n = 7
total_bits = n*(n-1)//2
np.random.seed(42)

results = []
for trial in range(300):
    bits = np.random.randint(0, 1 << total_bits)
    A = bits_to_adj(bits, n)
    H = count_ham_paths(A, n)
    cycles = get_directed_cycles(A, n)
    a1, a2 = compute_alpha_1_2(cycles)
    results.append({'H': H, 'a1': a1, 'a2': a2, 'bits': bits})

# Verify OCF
ocf_ok = sum(1 for r in results if r['H'] == 1 + 2*r['a1'] + 4*r['a2'])
print(f"  OCF (H = 1+2*a1+4*a2) verified: {ocf_ok}/{len(results)}")

# 2-adic
bit1_ok = sum(1 for r in results if ((r['H']>>1)&1) == (r['a1']%2))
print(f"  bit 1 = a1 mod 2: {bit1_ok}/{len(results)}")

# H mod 3
mod3_ok = sum(1 for r in results if r['H']%3 == (1-r['a1']+r['a2'])%3)
print(f"  H mod 3 = (1-a1+a2) mod 3: {mod3_ok}/{len(results)}")

# Statistics
a1_vals = [r['a1'] for r in results]
a2_vals = [r['a2'] for r in results]
print(f"\n  alpha_1: [{min(a1_vals)}, {max(a1_vals)}], mean={np.mean(a1_vals):.1f}")
print(f"  alpha_2: [{min(a2_vals)}, {max(a2_vals)}], mean={np.mean(a2_vals):.1f}")

# H distribution
h_dist = Counter(r['H'] for r in results)
print(f"  H range: [{min(h_dist)}, {max(h_dist)}]")
print(f"  H mod 2: {dict(Counter(r['H']%2 for r in results))}")
print(f"  H mod 3: {dict(Counter(r['H']%3 for r in results))}")
print(f"  H mod 4: {dict(Counter(r['H']%4 for r in results))}")
print(f"  H mod 6: {dict(Counter(r['H']%6 for r in results))}")
print(f"  H mod 12: {dict(Counter(r['H']%12 for r in results))}")

# ============================================================
# PART 4: The alpha_1 - alpha_2 relationship
# ============================================================
print("\n" + "=" * 70)
print("PART 4: alpha_1 vs alpha_2 at n=7")
print("=" * 70)

# Is there a functional relationship between a1 and a2?
a1_to_a2 = defaultdict(list)
for r in results:
    a1_to_a2[r['a1']].append(r['a2'])

print(f"  a1 -> range of a2:")
for a1 in sorted(a1_to_a2.keys()):
    a2_list = a1_to_a2[a1]
    a2_set = sorted(set(a2_list))
    print(f"    a1={a1:>3}: a2 in {a2_set} (H in {[1+2*a1+4*a for a in a2_set]})")

# ============================================================
# PART 5: The "2-3 CRT tower"
# ============================================================
print("\n" + "=" * 70)
print("PART 5: The 2-3 CRT tower at n=7")
print("=" * 70)

print("""
Level 0: H mod 1 = 0 (trivial)
Level 1: H mod 2 = 1 (Redei)
Level 2: H mod 3 = (1-a1+a2) mod 3 = I(Omega,-1) mod 3
Level 3: H mod 4 = 1 + 2*(a1 mod 2) (from bit 0=1, bit 1=a1 mod 2)
Level 4: H mod 6 = CRT(H mod 2, H mod 3) = CRT(1, I(Omega,-1) mod 3)
Level 5: H mod 12 = CRT(H mod 4, H mod 3)

The number 2 controls the PARITY structure (Redei, bit extraction).
The number 3 controls the TOPOLOGICAL structure (independence complex).
""")

# CRT verification
crt_ok = 0
for r in results:
    a1p = r['a1'] % 2
    i_neg1_mod3 = (1 - r['a1'] + r['a2']) % 3

    h4 = 1 + 2 * a1p  # H mod 4

    # CRT: find x s.t. x = h4 mod 4, x = i_neg1_mod3 mod 3, 0 <= x < 12
    for x in range(12):
        if x % 4 == h4 and x % 3 == i_neg1_mod3:
            if r['H'] % 12 == x:
                crt_ok += 1
            break

print(f"  H mod 12 = CRT(H mod 4, H mod 3): {crt_ok}/{len(results)}")

# What H mod 12 values are possible?
print(f"  Possible H mod 12 values: {sorted(set(r['H']%12 for r in results))}")

# ============================================================
# PART 6: The 2-valuation of H-1 and its meaning
# ============================================================
print("\n" + "=" * 70)
print("PART 6: 2-adic valuation of H-1")
print("=" * 70)

v2_dist = Counter()
for r in results:
    h_minus_1 = r['H'] - 1
    if h_minus_1 == 0:
        v2 = float('inf')
    else:
        v2 = 0
        temp = h_minus_1
        while temp % 2 == 0:
            v2 += 1
            temp //= 2
    v2_dist[v2] += 1

print(f"  v_2(H-1) distribution: {dict(sorted(v2_dist.items()))}")
print(f"  v_2(H-1) = 1: H = 3 mod 4 => a1 odd")
print(f"  v_2(H-1) >= 2: H = 1 mod 4 => a1 even")

# When a1 is even, what determines v_2(H-1)?
# H-1 = 2*a1 + 4*a2 = 2*(a1 + 2*a2)
# v_2(H-1) = 1 + v_2(a1 + 2*a2)
# If a1 odd: v_2(a1+2*a2) = 0, so v_2(H-1) = 1
# If a1 even: v_2(a1+2*a2) = v_2(a1/2 + a2) ... depends on a1/2 + a2

for r in results[:20]:
    h_m1 = r['H'] - 1
    a1, a2 = r['a1'], r['a2']
    v2 = 0
    temp = h_m1
    while temp > 0 and temp % 2 == 0:
        v2 += 1
        temp //= 2

# ============================================================
# PART 7: H mod 9 and higher 3-adic structure
# ============================================================
print("\n" + "=" * 70)
print("PART 7: 3-adic structure at n=7")
print("=" * 70)

h_mod9 = Counter()
for r in results:
    h_mod9[r['H'] % 9] += 1
print(f"  H mod 9: {dict(sorted(h_mod9.items()))}")

# I(Omega, -1) mod 9
# I(Omega, -1) = 1 - a1 + a2
# H mod 9 = I(Omega, 2) mod 9
# Since 2^1=2, 2^2=4, 2^3=8=-1, 2^4=16=7=-2, 2^5=32=5, 2^6=64=1 (mod 9)
# So 2 has order 6 mod 9
# H mod 9 = 1 + 2*a1 + 4*a2 mod 9

# But I(Omega, -1) mod 9 = 1 - a1 + a2 mod 9
# Is there a relationship?
# H mod 9 = 1 + 2*a1 + 4*a2 mod 9
# I(-1) mod 9 = 1 - a1 + a2 mod 9
# Difference: H - I(-1) = 3*a1 + 3*a2 = 3*(a1+a2) mod 9
# So: H = I(-1) + 3*(a1+a2) mod 9
# OR: H = I(-1) mod 3 (which we know), and H mod 9 requires knowing a1+a2 mod 3

print(f"\n  H = I(-1) + 3*(a1+a2) mod 9:")
ok = sum(1 for r in results
         if r['H'] % 9 == ((1-r['a1']+r['a2']) + 3*(r['a1']+r['a2'])) % 9)
print(f"  Verified: {ok}/{len(results)}")

# Simplify: 1 - a1 + a2 + 3*a1 + 3*a2 = 1 + 2*a1 + 4*a2. Tautology!
print(f"  (This is tautological: 1-a1+a2+3(a1+a2) = 1+2a1+4a2)")

# What's non-trivial: does (a1+a2) mod 3 have any structure?
a1_plus_a2_mod3 = Counter()
for r in results:
    a1_plus_a2_mod3[(r['a1']+r['a2']) % 3] += 1
print(f"\n  (a1+a2) mod 3 distribution: {dict(sorted(a1_plus_a2_mod3.items()))}")

# ============================================================
# PART 8: The DEEP identity: what is I(Omega, 3)?
# ============================================================
print("\n" + "=" * 70)
print("PART 8: I(Omega, 3) at n=5 and n=7")
print("=" * 70)

# At n=5: I(Omega, 3) = 1 + 3*a1 + 9*a2
# At n=7: I(Omega, 3) = 1 + 3*a1 + 9*a2 (same since a_k=0 for k>=3)

for n_test in [5, 7]:
    total_bits_test = n_test*(n_test-1)//2
    np.random.seed(42)

    i3_dist = Counter()
    n_samp = min(2**total_bits_test, 500) if n_test <= 5 else 300

    for trial in range(n_samp):
        if n_test <= 5 and trial < 2**total_bits_test:
            bits = trial
        else:
            bits = np.random.randint(0, 1 << total_bits_test)

        A = bits_to_adj(bits, n_test)
        cycles = get_directed_cycles(A, n_test)
        a1, a2 = compute_alpha_1_2(cycles)
        i3 = 1 + 3*a1 + 9*a2

        i3_dist[i3] += 1

    print(f"\n  n={n_test}: I(Omega, 3) distribution:")
    for val in sorted(i3_dist.keys()):
        if i3_dist[val] >= 3 or n_test >= 7:
            # Extract a1, a2 from I(3) = 1 + 3*a1 + 9*a2
            i3_m1 = val - 1
            a2_from_i3 = i3_m1 // 9
            a1_from_i3 = (i3_m1 - 9*a2_from_i3) // 3
            if i3_dist[val] >= 3:
                print(f"    I(3)={val:>5} (a1={a1_from_i3}, a2={a2_from_i3}): count={i3_dist[val]}")

    # I(Omega, 3) mod 3 = 1 always (since 3*a1 + 9*a2 = 0 mod 3)
    i3_mod3 = Counter(v % 3 for v in i3_dist.keys() for _ in range(i3_dist[v]))
    print(f"  I(Omega, 3) mod 3: {dict(sorted(i3_mod3.items()))}")
    print(f"  Always 1 mod 3: {all(v%3 == 1 for v in i3_dist.keys())}")

    # I(Omega, 3) / 3 rounded
    # I(Omega, 3) = 1 + 3*(a1 + 3*a2) so (I(3)-1)/3 = a1 + 3*a2
    # This is a DIFFERENT weighting than (I(2)-1)/2 = a1 + 2*a2
    # The ratio I(3)/I(2) = (1+3*a1+9*a2)/(1+2*a1+4*a2)

    if n_test == 5:
        print(f"\n  Ratio I(3)/I(2) at n=5:")
        for d in h_data[:10]:
            a1, a2 = d['a1'], d['a2']
            i2 = 1 + 2*a1 + 4*a2
            i3 = 1 + 3*a1 + 9*a2
            print(f"    a1={a1}, a2={a2}: I(2)={i2}, I(3)={i3}, ratio={i3/i2:.4f}")

# ============================================================
# PART 9: The "2 vs 3" weight in independent sets
# ============================================================
print("\n" + "=" * 70)
print("PART 9: Why x=2? The minimum positive integer for OCF")
print("=" * 70)

print("""
At x=1: I(Omega, 1) = sum alpha_k = total # independent sets.
  This counts ALL independent sets equally.

At x=2: I(Omega, 2) = H = Hamiltonian paths.
  Weight 2^k gives independent sets of size k a weight 2^k.
  The "2" gives EXPONENTIAL weight to larger independent sets.

At x=3: I(Omega, 3) = 1 + 3*a1 + 9*a2.
  Even MORE weight on larger independent sets.

KEY QUESTION: Why does x=2 count Hamiltonian paths?

ANSWER (from Grinberg-Stanley):
  H(T) = sum over sigma in S_n with all nontrivial cycles being directed
  odd cycles in T, of 2^{number of nontrivial cycles of sigma}

  The factor 2 comes from: each directed odd cycle has 2 orientations,
  and the GS formula picks one for each cycle independently.

  If we used 3 orientations (impossible in tournaments), we'd get I(Omega, 3).
  If we used 1 orientation (fixed), we'd get I(Omega, 1).

  x=2 is the UNIQUE value that counts Hamiltonian paths because
  tournaments have EXACTLY 2 directed Hamiltonian cycles per odd vertex set
  (forward and backward, when a Hamiltonian cycle exists).

Wait: at k=3, a 3-vertex tournament has AT MOST 1 directed 3-cycle
(if it's a cyclic tournament) with 2 directions: a->b->c->a and a->c->b->a.
But those are DIFFERENT directed cycles. So a 3-vertex cyclic tournament
contributes 2 to alpha_1 (2 directed cycles on the same vertex set).

No wait -- I need to check: does a 3-vertex cyclic tournament
have 1 or 2 directed Hamiltonian cycles?

If a->b->c->a: this IS a directed Ham cycle.
The reverse a->c->b->a would need c->b, b->a, a->c, which is NOT in T
(since T is a tournament, either a->b or b->a but not both).

So a 3-vertex tournament has AT MOST 1 directed 3-cycle (Hamiltonian cycle).
NOT 2! In a tournament, the "reverse" of a 3-cycle is NOT a valid cycle.

Then why factor 2? The GS formula says:
  For each sigma with all nontrivial cycles being directed odd cycles,
  contribute 2^(# nontrivial cycles).

  A permutation sigma = (abc)(identity on rest) contributes 2^1 = 2
  IF abc is a directed cycle of T. The factor 2 is NOT "two orientations"
  but rather the fixed point of the formula.

  More precisely: sigma = (1)(2)...(n) (identity) contributes 2^0 = 1.
  This is the "1" in H = 1 + 2*a1 + ...

  sigma = (abc)(identity on rest) contributes 2^1 = 2 for each directed cycle.
  This gives 2*a1.

  sigma = (abc)(def)(identity on rest) contributes 2^2 = 4 for each pair
  of directed cycles on disjoint vertex sets. This gives 4*a2.

So the factor 2 in the GS formula is NOT about "two orientations."
It's a combinatorial factor from the signed version of the permanent.

ACTUALLY: Let me re-read GS. Their Corollary 20 says:
  ham(D) = sum_{sigma: all cycles odd D-cycles} 2^{psi(sigma)}
  where psi(sigma) = # cycles of sigma that are NOT fixed points.

  For identity: psi = 0, contribution = 1.
  For (abc)(rest fixed): psi = 1, contribution = 2.
  For (abc)(def)(rest fixed): psi = 2, contribution = 4.

  This IS 2^(# nontrivial cycles) = x^(# vertices in I.S.) at x=2
  via the correspondence: nontrivial cycles <-> independent set in Omega.

The 2 ultimately comes from the relationship between the permanent
and the independence polynomial, evaluated at x=2.
""")

# ============================================================
# PART 10: Summary — the complete 2-3 picture
# ============================================================
print("\n" + "=" * 70)
print("PART 10: SUMMARY — The Complete 2-3 Picture")
print("=" * 70)

print("""
THE NUMBER 2:
  1. H = I(Omega, 2): the OCF evaluation point
  2. 2^k weight from GS formula: each nontrivial cycle contributes factor 2
  3. H is always ODD: I(Omega, 2) = I(Omega, 0) = 1 (mod 2)
  4. bit 1 of H = alpha_1 mod 2: parity of total directed odd cycles
  5. H mod 4 = 1 or 3, determined by alpha_1 parity
  6. v_2(H-1) encodes the 2-adic structure of the cycle decomposition

THE NUMBER 3:
  1. 3-cycles are the "atoms" of tournament structure
  2. H mod 3 = I(Omega, -1) mod 3 (since 2 = -1 mod 3)
  3. I(Omega, -1) = reduced Euler char of Ind(Omega) + 1 = topological
  4. tr(Sigma^3) always divisible by 6 (factor 3 from symmetric trace)
  5. 9 | F(T, omega) universally for n >= 6 (THM-085)
  6. THM-086: 3 | c_j for j < val(n)
  7. sigma:lambda:delta = 2:1:1 — the "3" in common neighborhood

THE INTERPLAY:
  1. 6 = 2*3 = lcm(2,3): H mod 6 is the CRT combination
  2. 12 = 4*3: H mod 12 = CRT(H mod 4, H mod 3) captures bit 1 + topology
  3. tr(Sigma^3) gaps always divisible by 12 = 4*3
  4. Worpitzky: 6 | W_3 (from 3! = 6)
  5. At n=7: H = 1 + 2*a1 + 4*a2, and H mod 3 = 1 - a1 + a2 mod 3
     So 3 "mixes" a1 and a2 while 2 "separates" them by bit position
""")

print("\nDone.")
