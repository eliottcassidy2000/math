#!/usr/bin/env python3
"""
nine_squared_deep.py — opus-2026-03-14-S71e

DEEP DIVE: 9 = 3² and the structure of α₃

Questions:
1. Why is α₂^{37} = 0 at n=9? (3+7=10 > 9, so impossible!)
2. What's the maximum α₃ achievable?
3. Does the "3-cycle packing" interpretation of 9=3² give insight?
4. How does I(-1) = 1 - α₁ + α₂ - α₃ behave at n=9?
5. The k-nacci / Vandermonde connection

Also explores: Can we predict α₃ from score sequence?
"""

import sys
import numpy as np
from itertools import combinations
from collections import Counter, defaultdict
import time
sys.stdout.reconfigure(line_buffering=True)

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
    for ms in range(2, n+1):
        for mask in range(1 << n):
            if bin(mask).count('1') != ms: continue
            for v in range(n):
                if not (mask & (1 << v)): continue
                pm = mask ^ (1 << v)
                t = 0
                for u in range(n):
                    if (pm & (1 << u)) and A[u][v]:
                        t += dp.get((pm, u), 0)
                if t:
                    dp[(mask, v)] = t
    return sum(dp.get(((1 << n) - 1, v), 0) for v in range(n))

def count_ham_cycles(A, n):
    if n < 3: return 0
    full_mask = (1 << n) - 1
    dp = {(1 << 0, 0): 1}
    for ms in range(2, n+1):
        for mask in range(1 << n):
            if bin(mask).count('1') != ms: continue
            if not (mask & 1): continue
            for v in range(n):
                if not (mask & (1 << v)): continue
                if v == 0 and ms < n: continue
                pm = mask ^ (1 << v)
                if not (pm & 1): continue
                t = 0
                for u in range(n):
                    if (pm & (1 << u)) and A[u][v]:
                        t += dp.get((pm, u), 0)
                if t:
                    dp[(mask, v)] = t
    total = 0
    for v in range(1, n):
        if A[v][0] and (full_mask, v) in dp:
            total += dp[(full_mask, v)]
    return total

def count_directed_k_cycles(A, n, k):
    if k > n: return 0
    total = 0
    for combo in combinations(range(n), k):
        verts = list(combo)
        sub = np.zeros((k, k), dtype=int)
        for i in range(k):
            for j in range(k):
                sub[i][j] = A[verts[i]][verts[j]]
        total += count_ham_cycles(sub, k)
    return total

def count_disjoint_pairs(A, n, k1, k2):
    """Count pairs of vertex-disjoint directed cycles of lengths k1, k2."""
    if k1 + k2 > n: return 0
    total = 0
    for c1 in combinations(range(n), k1):
        s1 = set(c1)
        sub1 = np.zeros((k1, k1), dtype=int)
        for i in range(k1):
            for j in range(k1):
                sub1[i][j] = A[c1[i]][c1[j]]
        hc1 = count_ham_cycles(sub1, k1)
        if hc1 == 0: continue
        remaining = [v for v in range(n) if v not in s1]
        for c2 in combinations(remaining, k2):
            sub2 = np.zeros((k2, k2), dtype=int)
            for i in range(k2):
                for j in range(k2):
                    sub2[i][j] = A[c2[i]][c2[j]]
            hc2 = count_ham_cycles(sub2, k2)
            total += hc1 * hc2
    # Divide by 2 if k1==k2 (each unordered pair counted twice)
    if k1 == k2:
        total //= 2
    return total

def count_disjoint_triples_333(A, n):
    """Count triples of vertex-disjoint 3-cycles."""
    if n < 9: return 0
    total = 0
    for c1 in combinations(range(n), 3):
        s1 = set(c1)
        sub1 = np.zeros((3, 3), dtype=int)
        for i in range(3):
            for j in range(3):
                sub1[i][j] = A[c1[i]][c1[j]]
        hc1 = count_ham_cycles(sub1, 3)
        if hc1 == 0: continue
        rem1 = [v for v in range(n) if v not in s1]
        for c2 in combinations(rem1, 3):
            s2 = set(c2)
            sub2 = np.zeros((3, 3), dtype=int)
            for i in range(3):
                for j in range(3):
                    sub2[i][j] = A[c2[i]][c2[j]]
            hc2 = count_ham_cycles(sub2, 3)
            if hc2 == 0: continue
            rem2 = [v for v in rem1 if v not in s2]
            for c3 in combinations(rem2, 3):
                sub3 = np.zeros((3, 3), dtype=int)
                for i in range(3):
                    for j in range(3):
                        sub3[i][j] = A[c3[i]][c3[j]]
                hc3 = count_ham_cycles(sub3, 3)
                total += hc1 * hc2 * hc3
    # Each unordered triple counted 3! = 6 times
    return total // 6

# ======================================================================
# PART 1: WHY α₂^{37} = 0 at n=9
# ======================================================================
print("=" * 70)
print("PART 1: CYCLE PAIR BUDGET AT EACH n")
print("=" * 70)

print("""
  Cycle pair type (k1,k2) requires k1+k2 vertices.
  At n=9, possible pairs:
    (3,3): 6 vertices — YES (exists since n=6)
    (3,5): 8 vertices — YES (exists since n=8)
    (3,7): 10 vertices — NO! 3+7=10 > 9
    (5,5): 10 vertices — NO! 5+5=10 > 9
    (5,7): 12 vertices — NO!
    (7,7): 14 vertices — NO!

  So at n=9: α₂ = α₂^{33} + α₂^{35} only.
  This explains why α₂^{37} = 0.

  At n=10 = 5·2:
    (3,3): YES     (3,5): YES     (3,7): YES (first time!)
    (5,5): YES (first time! — this is the 5·2 connection)
    So α₂ = α₂^{33} + α₂^{35} + α₂^{37} + α₂^{55}

  At n=11:
    (3,3): YES     (3,5): YES     (3,7): YES
    (5,5): YES     (5,7): 12>11 NO
    So α₂ same as n=10.

  FIRST APPEARANCE OF EACH PAIR TYPE:
""")

for k1 in range(3, 12, 2):
    for k2 in range(k1, 12, 2):
        n_min = k1 + k2
        if n_min <= 20:
            print(f"    ({k1},{k2}): first at n={n_min}")

# ======================================================================
# PART 2: 3-CYCLE PACKING AND 9 = 3²
# ======================================================================
print("\n" + "=" * 70)
print("PART 2: 3-CYCLE PACKING — THE 9=3² STRUCTURE")
print("=" * 70)

n = 9
tb = n*(n-1)//2
np.random.seed(42)

# Find tournaments that MAXIMIZE α₃ (perfect 3-cycle packing)
print("\n  Searching for max α₃ at n=9 (100 random tournaments)...")

data = []
t0 = time.time()
for trial in range(100):
    bits = np.random.randint(0, 1 << tb)
    A = bits_to_adj(bits, n)

    # Count all cycle types
    dc3 = count_directed_k_cycles(A, n, 3)
    dc5 = count_directed_k_cycles(A, n, 5)
    dc7 = count_directed_k_cycles(A, n, 7)
    dc9 = count_ham_cycles(A, n)
    a1 = dc3 + dc5 + dc7 + dc9

    dp33 = count_disjoint_pairs(A, n, 3, 3)
    dp35 = count_disjoint_pairs(A, n, 3, 5)
    a2 = dp33 + dp35

    a3 = count_disjoint_triples_333(A, n)

    scores = tuple(sorted([sum(A[i]) for i in range(n)]))
    data.append({
        'dc3': dc3, 'dc5': dc5, 'dc7': dc7, 'dc9': dc9,
        'a1': a1, 'dp33': dp33, 'dp35': dp35, 'a2': a2, 'a3': a3,
        'scores': scores
    })

print(f"  Done: {time.time()-t0:.1f}s")

# Sort by α₃
data.sort(key=lambda d: d['a3'], reverse=True)

print("\n  Top 5 by α₃:")
for i in range(5):
    d = data[i]
    print(f"    α₃={d['a3']:3d}, dc3={d['dc3']:3d}, α₂={d['a2']:3d} (33:{d['dp33']},35:{d['dp35']}), scores={d['scores']}")

print("\n  Bottom 5 by α₃ (α₃=0):")
zeros = [d for d in data if d['a3'] == 0]
for d in zeros[:5]:
    print(f"    α₃={d['a3']}, dc3={d['dc3']:3d}, α₂={d['a2']:3d}, scores={d['scores']}")

# ======================================================================
# PART 3: 9 = 3² and the TRANSITIVE-REGULAR SPECTRUM
# ======================================================================
print("\n" + "=" * 70)
print("PART 3: SCORE SEQUENCE → α₃ PREDICTION")
print("=" * 70)

# Group by score sequence
score_groups = defaultdict(list)
for d in data:
    score_groups[d['scores']].append(d)

print(f"\n  {len(score_groups)} distinct score sequences in 100 samples")
print(f"\n  Score → α₃ ranges:")
for scores in sorted(score_groups.keys()):
    vals = [d['a3'] for d in score_groups[scores]]
    if len(vals) >= 2:
        print(f"    {scores}: α₃ in [{min(vals)}, {max(vals)}] ({len(vals)} samples)")

# ======================================================================
# PART 4: I(-1) AT n=9
# ======================================================================
print("\n" + "=" * 70)
print("PART 4: I(Omega, -1) AT n=9")
print("=" * 70)

print("\n  I(-1) = 1 - α₁ + α₂ - α₃")
print("  At n≤8: α₃=0, so I(-1)=1-α₁+α₂")
print("  At n=9: the -α₃ term makes I(-1) even MORE negative")

I_neg1_vals = [1 - d['a1'] + d['a2'] - d['a3'] for d in data]
print(f"\n  I(-1) statistics at n=9 (100 samples):")
print(f"    mean: {np.mean(I_neg1_vals):.1f}")
print(f"    range: [{min(I_neg1_vals)}, {max(I_neg1_vals)}]")
print(f"    I(-1) > 0: {sum(1 for v in I_neg1_vals if v > 0)}")
print(f"    I(-1) = 0: {sum(1 for v in I_neg1_vals if v == 0)}")
print(f"    I(-1) < 0: {sum(1 for v in I_neg1_vals if v < 0)}")

# For transitive tournament: dc3=0 → α₁=0, α₂=0, α₃=0
# I(-1) = 1. So transitive still has I(-1)=1>0.
# But at n=8 already I(-1)≤-6. What about n=9 transitive?
print("\n  Transitive tournament at n=9: bits=0")
A_trans = bits_to_adj(0, n)
dc3_t = count_directed_k_cycles(A_trans, n, 3)
H_trans = count_ham_paths(A_trans, n)
print(f"    dc3={dc3_t}, H={H_trans}")
print(f"    (Transitive: dc3=0, α₁=0, α₂=0, α₃=0, I(-1)=1)")
print(f"    But n=8 showed I(-1)∈[-158,-6] for ALL tournaments.")
print(f"    At n=9: is there any T with I(-1)≥0?")

# Check a few near-transitive
print("\n  Near-transitive tournaments (flip 1 edge from transitive):")
for flip_idx in range(5):
    bits = 1 << flip_idx  # flip one edge
    A = bits_to_adj(bits, n)
    dc3 = count_directed_k_cycles(A, n, 3)
    dc5 = count_directed_k_cycles(A, n, 5)
    dc7 = count_directed_k_cycles(A, n, 7)
    dc9 = count_ham_cycles(A, n)
    a1 = dc3 + dc5 + dc7 + dc9
    dp33 = count_disjoint_pairs(A, n, 3, 3)
    dp35 = count_disjoint_pairs(A, n, 3, 5)
    a2 = dp33 + dp35
    a3 = count_disjoint_triples_333(A, n)
    I_neg1 = 1 - a1 + a2 - a3
    print(f"    flip {flip_idx}: dc3={dc3}, α₁={a1}, α₂={a2}, α₃={a3}, I(-1)={I_neg1}")

# ======================================================================
# PART 5: THE k-NACCI / VANDERMONDE CONNECTION
# ======================================================================
print("\n" + "=" * 70)
print("PART 5: k-NACCI CONVERGENCE AND EVALUATION POINTS")
print("=" * 70)

print("""
  The k-nacci golden ratio phi_k approaches 2 as k -> infinity:
    phi_k = 2 - 1/2^k + O(k/4^k)

  The weighted k-nacci ratio psi_k approaches 3:
    psi_k = 3 - (2/3)^k + O(k*(4/9)^k)

  KEY INSIGHT: The evaluation points for I(Omega, x) are 2, 3, 4, 5, ...
  The k-nacci roots APPROACH these evaluation points from below.

  Precisely:
    k-nacci root phi_k: 2 - epsilon_2(k)  where epsilon_2 = 1/2^k
    weighted root psi_k: 3 - epsilon_3(k)  where epsilon_3 = (2/3)^k

  At tournament size n:
    - We need evaluation points 2, 3, ..., floor(n/3)+1
    - The k-nacci roots for k ~ n are VERY close to these integers
    - The "error" epsilon_m(n) = ((m-1)/m)^n

  For n=9 = 3^2:
    eps_2(9) = 1/512 ~ 0.002  (effectively 2)
    eps_3(9) = (2/3)^9 ~ 0.026 (close to 3)
    eps_4(9) = (3/4)^9 ~ 0.075 (somewhat close to 4)

  The PRECISION HIERARCHY:
    2 is known to precision 1/2^n
    3 is known to precision (2/3)^n
    4 is known to precision (3/4)^n

  The ratio of consecutive precisions: (err_{m+1}/err_m) = (m/(m-1))^n · ((m-1)/m)
  For large n, this ratio ~ (m/(m-1))^n, growing EXPONENTIALLY.

  So "knowing 2" is EXPONENTIALLY more precise than "knowing 3",
  which is exponentially more precise than "knowing 4".

  THIS IS WHY 2 AND 3 ARE SPECIAL:
    - 2 is the dominant eigenvalue (most information)
    - 3 gives the first correction (needed only when α₂ > 0, i.e., n ≥ 6)
    - 4 gives the second correction (needed only when α₃ > 0, i.e., n ≥ 9)

  The "keys to the universe" hierarchy:
    Level 0: know 2 → know H (sufficient for n ≤ 5)
    Level 1: know 2,3 → know α₁,α₂ (sufficient for n ≤ 8)
    Level 2: know 2,3,4 → know α₁,α₂,α₃ (sufficient for n ≤ 11)
    Level k: know 2,...,k+2 → know α₁,...,α_{k+1} (sufficient for n ≤ 3k+2)

  But 4 = 2², so "knowing 4" IS knowing 2 at second order!
  And 8 = 2³, and 9 = 3² — these are COMBINATIONS of knowing 2 and 3.

  The "truly new" information enters at primes:
    2 (Level 0), 3 (Level 1), 5 (Level 3, n ≥ 15), 7 (Level 5, n ≥ 21)
""")

# Compute the "information content" at each level
print("  Information hierarchy — Vandermonde determinants and prime content:")
for k in range(1, 8):
    points = list(range(2, k+3))
    # Vandermonde det = product of (p_j - p_i) for i < j
    det = 1
    for i in range(len(points)):
        for j in range(i+1, len(points)):
            det *= (points[j] - points[i])
    # Factor
    d = abs(det)
    primes_in = []
    for p in [2, 3, 5, 7, 11, 13, 17, 19]:
        if d % p == 0:
            exp = 0
            while d % p == 0:
                exp += 1
                d //= p
            primes_in.append(f"{p}^{exp}" if exp > 1 else str(p))
    n_sufficient = 3*k + 2
    print(f"    Level {k}: points {points}, det={abs(det)}, "
          f"primes={'*'.join(primes_in)}, sufficient n<={n_sufficient}")

# ======================================================================
# PART 6: THE 3^2 STRUCTURE — PERFECT PACKING
# ======================================================================
print("\n" + "=" * 70)
print("PART 6: PERFECT 3-CYCLE PACKING AT n=9")
print("=" * 70)

# A "perfect 3-cycle packing" means partitioning all 9 vertices into
# three disjoint 3-cycles. This is a covering of the vertex set.
#
# How many ways can this be done? Depends on the tournament.
# The MAXIMUM α₃ is when the most triples of disjoint 3-cycles exist.

# For the regular tournament C_9^{1,2,3,4} (circulant):
# Let's construct it manually
print("\n  Circulant tournament C_9^{1,2,3,4}:")
A_circ = np.zeros((9, 9), dtype=int)
for i in range(9):
    for d in [1, 2, 3, 4]:
        A_circ[i][(i+d) % 9] = 1
scores_circ = tuple(sorted([sum(A_circ[i]) for i in range(9)]))
print(f"    Scores: {scores_circ}")

dc3_c = count_directed_k_cycles(A_circ, 9, 3)
dc5_c = count_directed_k_cycles(A_circ, 9, 5)
dc7_c = count_directed_k_cycles(A_circ, 9, 7)
dc9_c = count_ham_cycles(A_circ, 9)
a1_c = dc3_c + dc5_c + dc7_c + dc9_c
dp33_c = count_disjoint_pairs(A_circ, 9, 3, 3)
dp35_c = count_disjoint_pairs(A_circ, 9, 3, 5)
a2_c = dp33_c + dp35_c
a3_c = count_disjoint_triples_333(A_circ, 9)
H_c = count_ham_paths(A_circ, 9)

print(f"    dc3={dc3_c}, dc5={dc5_c}, dc7={dc7_c}, dc9={dc9_c}")
print(f"    α₁={a1_c}, α₂={a2_c} (33:{dp33_c}, 35:{dp35_c}), α₃={a3_c}")
print(f"    H={H_c}")
print(f"    I(-1) = {1 - a1_c + a2_c - a3_c}")

# How does α₃ relate to dc3?
print("\n  Correlation: dc3 vs α₃")
dc3_vals = [d['dc3'] for d in data]
a3_vals = [d['a3'] for d in data]
corr = np.corrcoef(dc3_vals, a3_vals)[0,1]
print(f"    corr(dc3, α₃) = {corr:.4f}")

# α₃ vs α₂^{33}
dp33_vals = [d['dp33'] for d in data]
corr2 = np.corrcoef(dp33_vals, a3_vals)[0,1]
print(f"    corr(α₂^{{33}}, α₃) = {corr2:.4f}")

# α₃ / C(dc3, 3) — what fraction of 3-cycle triples are disjoint?
print("\n  Fraction of 3-cycle triples that are disjoint:")
for d in data[:10]:
    from math import comb
    if d['dc3'] >= 3:
        total_triples = comb(d['dc3'], 3)
        frac = d['a3'] / total_triples if total_triples > 0 else 0
        print(f"    dc3={d['dc3']:3d}, C(dc3,3)={total_triples:6d}, α₃={d['a3']:3d}, ratio={frac:.6f}")

# ======================================================================
# PART 7: THE VANDERMONDE AS GENERATING FUNCTION
# ======================================================================
print("\n" + "=" * 70)
print("PART 7: VANDERMONDE GENERATING FUNCTION INTERPRETATION")
print("=" * 70)

print("""
  I(Omega, x) = sum_k alpha_k x^k is a GENERATING FUNCTION.

  The "keys" 2 and 3 are the evaluation points that EXTRACT the coefficients.

  Consider the Z-transform / probability generating function viewpoint:
    I(x) evaluated at x=2 weights each independent k-set by 2^k.
    I(x) evaluated at x=3 weights each independent k-set by 3^k.

  So H = I(2) = sum_k alpha_k 2^k is a BINARY-WEIGHTED count.
  And I_3 = I(3) = sum_k alpha_k 3^k is a TERNARY-WEIGHTED count.

  The Hamiltonian paths are the BINARY counting of cycle independence!

  Now: what is the NATURAL basis for tournament information?
  If we define x_k = alpha_k / alpha_max_k where alpha_max_k is the max
  possible alpha_k, then x_k in [0,1] and
    I(x) = 1 + alpha_max_1 x_1 x + alpha_max_2 x_2 x^2 + ...

  This normalizes the information content of each level.

  At n=9:
    alpha_max_1 ~ 819 (from data)
    alpha_max_2 ~ 323
    alpha_max_3 ~ 16

  The ratio alpha_max_k / alpha_max_{k-1} shrinks rapidly.
  Higher-order information is increasingly RARE.

  THIS is why knowing 2 and 3 "gives the keys":
  The alpha_k values are SPARSE at high k, so the higher evaluation
  points contribute very little to the total information.
""")

# Compute information concentration
print("  Information concentration (how much of I(2) comes from each level):")
for d in data[:5]:
    total = 1
    for k in range(1, 4):
        if k == 1: val = d['a1']
        elif k == 2: val = d['a2']
        elif k == 3: val = d['a3']
        total += val * (2**k)
    H_val = total
    pct1 = 2 * d['a1'] / H_val * 100 if H_val > 0 else 0
    pct2 = 4 * d['a2'] / H_val * 100 if H_val > 0 else 0
    pct3 = 8 * d['a3'] / H_val * 100 if H_val > 0 else 0
    print(f"    H={H_val}: level1={pct1:.1f}%, level2={pct2:.1f}%, level3={pct3:.1f}%")

# ======================================================================
# PART 8: THE DEEP 2-3 DUALITY
# ======================================================================
print("\n" + "=" * 70)
print("PART 8: WHY 2 AND 3 — THE DEEP DUALITY")
print("=" * 70)

print("""
  THESIS: 2 and 3 are the keys because they are the SIMPLEST PRIMES
  and the FIRST TWO POINTS of the Vandermonde interpolation scheme.

  But there's a DEEPER reason connected to tournament structure:

  1. "2" = the binary nature of tournaments (each edge has 2 orientations)
     → H = I(2) counts Hamiltonian paths (binary evaluation)
     → Redei's theorem: H ≡ 1 (mod 2) — PARITY is the first invariant

  2. "3" = the smallest odd cycle (triangle)
     → I(3) uses ternary weights — EXACTLY the weight needed to "see" triangles
     → dc3 = C(n,3) - sum C(s_i,2) — triangle count is the FIRST non-trivial invariant
     → I_3 ≡ 1 (mod 3) — TERNARY parity is the second invariant

  3. "6 = 2*3" = the Vandermonde determinant = LCM of the first two primes
     → 6 = 3! — also the first non-trivial factorial
     → The integrality constraint mod 6 combines BOTH parities

  4. The k-nacci approach rates:
     → phi_k → 2 at rate 1/2^k (exponential in k)
     → psi_k → 3 at rate (2/3)^k (exponential, but slower)
     → The RATIO of precision losses = (4/3)^k
     → At n=9: 3 is 13x less precisely known than 2
     → At n=20: 3 is 1046x less precisely known than 2

  5. The 9 = 3^2 connection:
     → At n=9, α₃ appears → need evaluation at x=4=2^2
     → "4 is the square of 2" — second-order knowledge of 2
     → The Vandermonde (2,3,4) has det 48 = 2^4 * 3 — STILL only primes 2,3
     → So the "3^2 threshold" is handled by "2^2 knowledge"

  6. The 10 = 5*2 connection:
     → At n=10, (5,5) pairs appear — the PRIME 5 gets its own pairing structure
     → But we STILL don't need evaluation at x=5 (since α₄ doesn't appear until n=12)
     → 5 enters the Vandermonde only at k=3 (points 2,3,4,5, det=1440=2^5*3^2*5)
     → This happens at n=9+, but the ACTUAL need for x=5 is n=15+ (α₅ first at n=15)

  SUMMARY: "2 and 3 are the keys to the universe" means:
  - At small n (≤8): LITERALLY TRUE. Two evaluation points suffice.
  - At n=9-11: Still true, because 4=2^2 carries no new prime information.
  - At n≥15: 5 enters as a genuinely new prime. The universe expands.
  - At ALL n: 2 and 3 are the DOMINANT evaluation points, carrying
    exponentially more information than higher points.
    The "keys" are always the most important, even when other keys exist.
""")

print("Done.")
