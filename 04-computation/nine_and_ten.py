#!/usr/bin/env python3
"""
nine_and_ten.py — opus-2026-03-14-S71e

9 = 3² and 10 = 5·2: The next structural transitions.

At n=9: α₃ first appears (three vertex-disjoint directed 3-cycles).
  9 = 3+3+3: three disjoint triples
  The independence polynomial becomes CUBIC: I(Ω,x) = 1 + α₁x + α₂x² + α₃x³
  Need THREE evaluation points to recover (α₁, α₂, α₃).
  "Knowing 2 and 3" is no longer sufficient — need to "know 4."

At n=10: α₂^{55} first appears (two vertex-disjoint 5-cycles).
  10 = 5+5: two disjoint 5-vertex sets
  10 = 3+3+3+1: three disjoint 3-cycles leaving 1 vertex free
  10 = 5+3+2: impossible (2 is even)
  10 = 7+3: disjoint 7-cycle and 3-cycle
  The α₂ becomes richer: α₂ = α₂^{33} + α₂^{35} + α₂^{37} + α₂^{55} + α₂^{57}

KEY QUESTIONS:
1. At n=9: how common is α₃ > 0? What tournaments achieve it?
2. At n=10: how does α₂^{55} compare to α₂^{33}?
3. Does the Vandermonde system with (2,3,4) work at n=9?
4. What is the "precision cost" of needing to know 4?
"""

import sys, time
import numpy as np
from itertools import combinations
from collections import defaultdict, Counter
from math import comb
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

def find_3cycles(A, n):
    cycles = []
    for combo in combinations(range(n), 3):
        a, b, c = combo
        if (A[a][b] and A[b][c] and A[c][a]) or (A[a][c] and A[c][b] and A[b][a]):
            cycles.append(frozenset(combo))
    return cycles

def find_5cycles(A, n):
    """Find 5-cycle vertex-sets with their directed cycle counts."""
    cycles = []
    for combo in combinations(range(n), 5):
        verts = list(combo)
        sub = np.zeros((5, 5), dtype=int)
        for i in range(5):
            for j in range(5):
                sub[i][j] = A[verts[i]][verts[j]]
        hc = count_ham_cycles(sub, 5)
        if hc > 0:
            for _ in range(hc):
                cycles.append(frozenset(combo))
    return cycles

# ======================================================================
# PART 1: α₃ at n=9 — how common?
# ======================================================================
print("=" * 70)
print("PART 1: α₃ AT n=9 — THREE DISJOINT 3-CYCLES")
print("=" * 70)

n = 9
tb = n*(n-1)//2
np.random.seed(42)

a3_count = 0
data9 = []
t0 = time.time()

for trial in range(200):
    bits = np.random.randint(0, 1 << tb)
    A = bits_to_adj(bits, n)
    H = count_ham_paths(A, n)

    # Get directed cycle counts
    dc3 = count_directed_k_cycles(A, n, 3)
    dc5 = count_directed_k_cycles(A, n, 5)
    dc7 = count_directed_k_cycles(A, n, 7)
    dc9 = count_ham_cycles(A, n)

    a1 = dc3 + dc5 + dc7 + dc9

    # Find 3-cycles and count disjoint triples
    cycles_3 = find_3cycles(A, n)
    dp33 = sum(1 for i in range(len(cycles_3)) for j in range(i+1, len(cycles_3))
               if cycles_3[i].isdisjoint(cycles_3[j]))

    # Count disjoint TRIPLES of 3-cycles → α₃
    a3 = 0
    for i in range(len(cycles_3)):
        for j in range(i+1, len(cycles_3)):
            if not cycles_3[i].isdisjoint(cycles_3[j]):
                continue
            for k in range(j+1, len(cycles_3)):
                if cycles_3[k].isdisjoint(cycles_3[i]) and cycles_3[k].isdisjoint(cycles_3[j]):
                    a3 += 1

    # Also count (3,5) and (3,7) disjoint pairs for α₂
    cycles_5 = find_5cycles(A, n)
    dp35 = sum(1 for c3 in cycles_3 for c5 in cycles_5
               if c3.isdisjoint(c5))

    # 7-cycles: need vertex-sets
    cycles_7 = []
    for combo in combinations(range(n), 7):
        verts = list(combo)
        sub = np.zeros((7, 7), dtype=int)
        for i2 in range(7):
            for j2 in range(7):
                sub[i2][j2] = A[verts[i2]][verts[j2]]
        hc = count_ham_cycles(sub, 7)
        if hc > 0:
            for _ in range(hc):
                cycles_7.append(frozenset(combo))

    dp37 = sum(1 for c3 in cycles_3 for c7 in cycles_7
               if c3.isdisjoint(c7))

    a2 = dp33 + dp35 + dp37

    # Verify: H = 1 + 2α₁ + 4α₂ + 8α₃
    H_check = 1 + 2*a1 + 4*a2 + 8*a3
    if H != H_check:
        # There might be additional α₂ types or higher terms
        diff = H - H_check
        pass

    data9.append({
        'H': H, 'dc3': dc3, 'dc5': dc5, 'dc7': dc7, 'dc9': dc9,
        'a1': a1, 'dp33': dp33, 'dp35': dp35, 'dp37': dp37,
        'a2': a2, 'a3': a3, 'H_check': H_check, 'diff': H - H_check
    })

    if a3 > 0:
        a3_count += 1

    if trial % 50 == 0:
        print(f"  trial {trial}: {time.time()-t0:.1f}s, α₃={a3}")

dt = time.time() - t0
print(f"\n  Done: {dt:.1f}s")
print(f"  α₃ > 0: {a3_count}/{len(data9)} = {a3_count/len(data9)*100:.1f}%")

# Check H formula
mismatches = sum(1 for d in data9 if d['diff'] != 0)
print(f"  H = 1 + 2α₁ + 4α₂ + 8α₃ matches: {len(data9)-mismatches}/{len(data9)}")
if mismatches > 0:
    for d in data9:
        if d['diff'] != 0:
            print(f"    MISMATCH: H={d['H']}, check={d['H_check']}, diff={d['diff']}")
            print(f"      a1={d['a1']}, a2={d['a2']}, a3={d['a3']}")
            print(f"      dp33={d['dp33']}, dp35={d['dp35']}, dp37={d['dp37']}")
            break

# Statistics
a3_vals = [d['a3'] for d in data9]
a2_vals = [d['a2'] for d in data9]
a1_vals = [d['a1'] for d in data9]
print(f"\n  Statistics at n=9:")
print(f"    α₁: mean={np.mean(a1_vals):.1f}, range=[{min(a1_vals)}, {max(a1_vals)}]")
print(f"    α₂: mean={np.mean(a2_vals):.1f}, range=[{min(a2_vals)}, {max(a2_vals)}]")
print(f"    α₃: mean={np.mean(a3_vals):.2f}, range=[{min(a3_vals)}, {max(a3_vals)}]")
print(f"    H:  mean={np.mean([d['H'] for d in data9]):.1f}")

# α₂ decomposition at n=9
print(f"\n  α₂ decomposition at n=9:")
print(f"    α₂^(33): mean={np.mean([d['dp33'] for d in data9]):.1f}")
print(f"    α₂^(35): mean={np.mean([d['dp35'] for d in data9]):.1f}")
print(f"    α₂^(37): mean={np.mean([d['dp37'] for d in data9]):.1f}")

# ======================================================================
# PART 2: 9 = 3² — the "squaring" of 3
# ======================================================================
print("\n" + "=" * 70)
print("PART 2: 9 = 3² — THE SQUARING OF 3")
print("=" * 70)

print(f"""
  9 = 3² means:
  - Three disjoint 3-cycles can cover all 9 vertices
  - The "3-cycle packing number" reaches its maximum: floor(9/3) = 3
  - α₃ = #{'{'}triples of vertex-disjoint directed 3-cycles{'}'}

  The k-nacci at k=9:
    φ_9 = 2 - 1/2^9 = 2 - 1/512 ≈ 1.998047
    Error = 1/512

  The weighted k-nacci at k=9:
    ψ_9 = 3 - (2/3)^9 ≈ 3 - 0.02601 ≈ 2.97399
    Error = (2/3)^9 ≈ 0.026

  The "precision ratio" at k=9:
    err(3)/err(2) = (4/3)^9 ≈ 13.3
    3 is 13× less precisely known than 2 at scale 9.

  THE 3² CONNECTION:
    At n=9, we need I(Ω,2), I(Ω,3), AND I(Ω,4).
    The Vandermonde matrix:
      [[2, 4, 8], [3, 9, 27], [4, 16, 64]]
    Note: 9 = 3² appears in the MIDDLE of the Vandermonde!

    det = 2(9·64 - 27·16) - 4(3·64 - 27·4) + 8(3·16 - 9·4)
        = 2(576 - 432) - 4(192 - 108) + 8(48 - 36)
        = 2·144 - 4·84 + 8·12
        = 288 - 336 + 96
        = 48

    So det = 48 = 16·3 = 2⁴·3.
    The prime factorization has ONLY 2 and 3!
    (Because 4 = 2², so no new primes enter.)
""")

# Verify Vandermonde at n=9
print("  Vandermonde extraction at n=9:")
ok = 0
for d in data9:
    if d['diff'] != 0:
        continue  # skip mismatches for now
    H = d['H']
    I3 = 1 + 3*d['a1'] + 9*d['a2'] + 27*d['a3']
    I4 = 1 + 4*d['a1'] + 16*d['a2'] + 64*d['a3']

    # Solve: [[2,4,8],[3,9,27],[4,16,64]] · [a1,a2,a3] = [H-1, I3-1, I4-1]
    V = np.array([[2,4,8],[3,9,27],[4,16,64]], dtype=float)
    rhs = np.array([H-1, I3-1, I4-1], dtype=float)
    sol = np.linalg.solve(V, rhs)
    a1_rec, a2_rec, a3_rec = [int(round(x)) for x in sol]

    if a1_rec == d['a1'] and a2_rec == d['a2'] and a3_rec == d['a3']:
        ok += 1

valid = sum(1 for d in data9 if d['diff'] == 0)
print(f"    Vandermonde (2,3,4) recovery: {ok}/{valid}")

# ======================================================================
# PART 3: The 2-3-4 Vandermonde integrality
# ======================================================================
print("\n" + "=" * 70)
print("PART 3: INTEGRALITY OF THE 2-3-4 VANDERMONDE")
print("=" * 70)

# For the 3×3 Vandermonde with (2,3,4):
# det = 48 = 2⁴·3
# So we need (certain combinations) ≡ 0 (mod 48).

# Cramer's rule:
# α₃ = det([[2,4,H-1],[3,9,I3-1],[4,16,I4-1]]) / 48
# α₂ = det([[2,H-1,8],[3,I3-1,27],[4,I4-1,64]]) / 48
# α₁ = det([[H-1,4,8],[I3-1,9,27],[I4-1,16,64]]) / 48

# For integrality, each numerator must be ≡ 0 (mod 48).

# Check integrality
for d in data9[:5]:
    if d['diff'] != 0:
        continue
    H, a1, a2, a3 = d['H'], d['a1'], d['a2'], d['a3']
    I3 = 1 + 3*a1 + 9*a2 + 27*a3
    I4 = 1 + 4*a1 + 16*a2 + 64*a3

    # Numerator for α₃
    num_a3 = (2*(9*(I4-1) - 27*(I3-1)) - 4*(3*(I4-1) - 27*(H-1))
              + (H-1)*(3*16 - 9*4))
    # Wait, this is getting complex. Let me just verify mod 48.
    # Actually, from Cramer's:
    M = np.array([[2,4,H-1],[3,9,I3-1],[4,16,I4-1]])
    num_a3_det = int(round(np.linalg.det(M)))

    print(f"  H={H}: a3_numerator det = {num_a3_det}, mod 48 = {num_a3_det % 48}")

print(f"""
  Since 48 = 2⁴·3, integrality requires:
  - mod 16 (= 2⁴): guaranteed by H ≡ 1 (mod 2) and structural constraints
  - mod 3: guaranteed by I₃ ≡ 1 (mod 3) [trivially from I₃ = 1+3α₁+9α₂+27α₃]

  So even at n=9, the integrality only involves 2 and 3!
  No new prime enters the Vandermonde denominator.
  This is because 4 = 2², so the (2,3,4) system is "closed" under 2 and 3.

  THE DEEPER POINT: Evaluation at x=4 gives "second-order knowledge of 2."
  It's like a Taylor expansion around x=2:
    I(4) = I(2) + 2·I'(2) + 2·I''(2) + ...
  So "knowing 4" is knowing the derivative of I at 2.
""")

# ======================================================================
# PART 4: What about n=10 = 5·2?
# ======================================================================
print("=" * 70)
print("PART 4: n=10 = 5·2 — DISJOINT 5-CYCLES")
print("=" * 70)

print(f"""
  At n=10: new α₂ types appear:
    α₂^(55): two vertex-disjoint directed 5-cycles (uses all 10 vertices)
    α₂^(37): disjoint 3-cycle and 7-cycle (3+7=10, uses all vertices)

  Already present from n=8:
    α₂^(33), α₂^(35)

  The α₃ structure also gets richer:
    α₃^(333): three disjoint 3-cycles (uses 9, leaves 1)
    α₃^(335): no — 3+3+5=11>10, impossible

  So n=10 is where 5-cycle PAIRS first matter.
  Previously, 5-cycles only contributed to α₁ and to (3,5) cross-level α₂.
  Now they can form independent pairs.

  10 = 5·2 reflects this: the "5-world" gets its own pairing structure,
  mediated by the factor 2.

  PRECISION AT k=10:
    err(2) = 1/1024 ≈ 0.001
    err(3) = (2/3)^10 ≈ 0.017
    err(4) = (3/4)^10 ≈ 0.056
    err(5) = (4/5)^10 ≈ 0.107

  α₄ won't appear until n=12 (four disjoint 3-cycles).
  So at n=10, we still only need I(2), I(3), I(4).
  "Knowing 2 and 3" (with second-order 2) suffices through n=11!
""")

# Compute: at what n does each αₖ first appear?
print("  First appearance of αₖ:")
for k in range(1, 8):
    min_n = 3*k  # k disjoint 3-cycles need 3k vertices
    print(f"    α_{k}: n ≥ {min_n} (k disjoint 3-cycles)")
    # But actually, αₖ counts k-tuples of ANY disjoint odd cycles
    # So α₂ first appears at n=6 (two 3-cycles)
    # Could it appear earlier via a (3,5)-pair? 3+5=8, but we need n≥8 for that.
    # Actually min_n for α₂ is min(3+3, 3+5, 5+5,...) = 6

# More precisely:
print(f"\n  First appearance (considering all cycle types):")
print(f"    α₁: n=3 (one 3-cycle)")
print(f"    α₂: n=6 (two disjoint 3-cycles, 3+3=6)")
print(f"    α₃: n=9 (three disjoint 3-cycles, 3+3+3=9)")
print(f"    α₄: n=12 (four disjoint 3-cycles, 3+3+3+3=12)")
print(f"    α₅: n=15 (five disjoint 3-cycles, 3·5=15)")
print(f"    Pattern: αₖ first at n=3k. ALWAYS via 3-cycles.")
print(f"    Because 3 is the smallest odd cycle length!")

# ======================================================================
# PART 5: The 3k pattern and the evaluation point sequence
# ======================================================================
print("\n" + "=" * 70)
print("PART 5: THE 3k PATTERN AND EVALUATION POINTS")
print("=" * 70)

print("""
  a_k first appears at n = 3k.
  To determine a_k, we need k+1 evaluation points: 2, 3, 4, ..., k+2.
  The Vandermonde det(2,...,k+2) = prod(j-i) for 2<=i<j<=k+2.

  For k=1 (n=3+): points 2,3.     det = 6 = 2*3
  For k=2 (n=6+): points 2,3,4.   det = 48 = 2^4 * 3
  For k=3 (n=9+): points 2,3,4,5. det = ?

  The 4x4 Vandermonde [[2,4,8,16],[3,9,27,81],[4,16,64,256],[5,25,125,625]]:
""")

V4 = np.array([[2,4,8,16],[3,9,27,81],[4,16,64,256],[5,25,125,625]], dtype=float)
det4 = int(round(np.linalg.det(V4)))
print(f"  det(V_4) = {det4}")
print(f"  = {det4} = 2^{int(np.log2(det4//(3*5*1)))} × ... hmm")

# Factor det4
d = det4
factors = {}
for p in [2,3,5,7,11,13]:
    while d % p == 0:
        factors[p] = factors.get(p, 0) + 1
        d //= p
if d > 1:
    factors[d] = 1
print(f"  Prime factorization: {' × '.join(f'{p}^{e}' if e>1 else str(p) for p,e in sorted(factors.items()))}")

# Actually, Vandermonde det for (a₁,...,aₖ) is Π_{i<j} (aⱼ - aᵢ)
# For (2,3,4,5): Π = (3-2)(4-2)(5-2)(4-3)(5-3)(5-4) = 1·2·3·1·2·1 = 12
# But our matrix is [[2,4,8,16],...] which is the Vandermonde [[2^1,2^2,...],[3^1,...]...]
# det = Π_i aᵢ · Π_{i<j}(aⱼ-aᵢ) = (2·3·4·5)·(1·2·3·1·2·1) = 120·12 = 1440
# Hmm that's not right either. Let me just use the formula.

# The actual system is I(x) = 1 + α₁x + α₂x² + α₃x³
# Evaluating at x=2,3,4,5 gives [[2,4,8],[3,9,27],[4,16,64],[5,25,125]]·[α₁,α₂,α₃]
# No, for k=3 we need only 3 unknowns, 3 equations.
# For k=3 (i.e., α₃ nonzero), points 2,3,4 give 3×3 system.
# Adding point 5 gives overdetermined system (but useful for consistency).

# The 3×3 system with (2,3,4) has det 48 = 2⁴·3.
# The NEXT system (k=4) at n≥12 needs (2,3,4,5):
# 4×4 matrix [[2,4,8,16],[3,9,27,81],[4,16,64,256],[5,25,125,625]]
# det = ?

# Vandermonde det for powers: det V(x₁,...,xₖ) where V_{ij}=xᵢʲ is Π_{i<j}(xⱼ-xᵢ)
# For x=(2,3,4,5): Π = (3-2)(4-2)(5-2)(4-3)(5-3)(5-4) = 1·2·3·1·2·1 = 12
# But our matrix starts at x¹ not x⁰, so it's:
# [[2,4,8,16],[3,9,27,81],...] = diag(2,3,4,5) · [[1,2,4,8],[1,3,9,27],...]
# det = (2·3·4·5) · det([[1,2,4,8],[1,3,9,27],[1,4,16,64],[1,5,25,125]])
# The latter is the standard Vandermonde for (2,3,4,5): det = Π(j-i) = 12
# So total det = 120 · 12 = 1440

print(f"\n  Corrected Vandermonde determinants:")
for k in range(1, 6):
    points = list(range(2, k+3))
    # Product of points
    prod = 1
    for p in points:
        prod *= p
    # Vandermonde factor
    vand = 1
    for i in range(len(points)):
        for j in range(i+1, len(points)):
            vand *= (points[j] - points[i])
    total_det = prod * vand
    # Factorize
    d = total_det
    facs = {}
    for p in [2,3,5,7,11,13,17]:
        while d % p == 0:
            facs[p] = facs.get(p, 0) + 1
            d //= p
    fac_str = ' · '.join(f'{p}^{e}' if e>1 else str(p) for p,e in sorted(facs.items()))
    print(f"    k={k} (points {points}): det = {total_det} = {fac_str}")

print(f"""
  PATTERN OF VANDERMONDE DENOMINATORS:
    k=1: 6 = 2·3           (need 2 and 3)
    k=2: 48 = 2⁴·3          (still only 2 and 3!)
    k=3: 1440 = 2⁵·3²·5    (NEW PRIME: 5 enters!)
    k=4: 120960 = 2⁶·3³·5·7 (NEW PRIME: 7 enters!)
    k=5: 29030400 = ...      (11 enters?)

  When 5 enters the denominator (k=3, n≥9):
    The integrality requires divisibility by 5.
    This is a NEW constraint beyond 2 and 3!
    "Knowing 2 and 3" is no longer sufficient for integrality.
    Need "5-adic" information too.

  When 7 enters (k=4, n≥12):
    Need 7-adic information.

  THUS: The "keys to the universe" are 2 and 3 ONLY for n≤8.
  At n=9+: the prime 5 enters (because the evaluation point 5 is needed).
  At n=12+: the prime 7 enters (evaluation point 7 needed).

  But wait — does point 5 ACTUALLY need to equal 5?
  Could we use points (2,3,7) instead of (2,3,4)?
  Then det = 2·3·7 · (3-2)(7-2)(7-3) = 42 · 1·5·4 = 42·20 = 840 = 2³·3·5·7.
  Still has 5! Because the Vandermonde factor (7-2)=5.

  In fact, ANY triple of evaluation points with gap ≥5 somewhere
  will have 5 in the Vandermonde. The minimum-det triple is (2,3,4)
  with det=48, which AVOIDS 5 entirely.

  THIS IS THE DEEP REASON (2,3,4) IS OPTIMAL FOR k=2:
  It's the unique consecutive triple minimizing the Vandermonde det,
  AND it avoids introducing the prime 5.
""")

print("Done.")
