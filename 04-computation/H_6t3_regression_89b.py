#!/usr/bin/env python3
"""
H ≈ 6·t3 + INTERCEPT — THE LINEAR H-TRIANGLE FORMULA
opus-2026-03-14-S89b

DISCOVERY: H(T) ≈ a(n) + 6·t3(T) where t3 = number of 3-cycles.

The slope is EXACTLY 6 = |S_3| for n=6. Is this exact for all n?
The intercept a(n) varies with n.

Also: the transitive tournament flip graph is (n-1)-regular.
This means each transitive tournament has exactly n-1 neighbors
that are also transitive. WHY?

A transitive tournament T_σ corresponds to a permutation σ ∈ S_n
(the score ordering). Two transitive tournaments T_σ and T_τ
are one-flip apart iff σ and τ differ by an ADJACENT TRANSPOSITION.
The number of adjacent transpositions from any permutation is n-1
(swap positions i and i+1 for i = 1, ..., n-1).

So the transitive flip graph IS the CAYLEY GRAPH of S_n
with generators = adjacent transpositions!
This is the PERMUTOHEDRON!
"""

from math import factorial, comb
from collections import Counter, defaultdict
from fractions import Fraction

def compute_H(n, adj):
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            if (mask, v) not in dp:
                continue
            for u in range(n):
                if mask & (1 << u):
                    continue
                if adj.get((v, u), 0) == 1:
                    new_mask = mask | (1 << u)
                    dp[(new_mask, u)] = dp.get((new_mask, u), 0) + dp[(mask, v)]
    full_mask = (1 << n) - 1
    return sum(dp.get((full_mask, v), 0) for v in range(n))

def tournament_from_bits(n, bits):
    adj = {}
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if (bits >> idx) & 1:
                adj[(i,j)] = 1
                adj[(j,i)] = 0
            else:
                adj[(i,j)] = 0
                adj[(j,i)] = 1
            idx += 1
    return adj

def count_3cycles(n, adj):
    count = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if (adj.get((i,j),0) and adj.get((j,k),0) and adj.get((k,i),0)) or \
                   (adj.get((i,k),0) and adj.get((k,j),0) and adj.get((j,i),0)):
                    count += 1
    return count

print("=" * 70)
print("H ≈ 6·t3 + INTERCEPT — THE LINEAR H-TRIANGLE FORMULA")
print("opus-2026-03-14-S89b")
print("=" * 70)

# ======================================================================
# PART 1: EXACT LINEAR REGRESSION H vs t3
# ======================================================================
print("\n" + "=" * 70)
print("PART 1: EXACT REGRESSION H = a + b·t3")
print("=" * 70)

for n in range(3, 7):
    m = n * (n-1) // 2

    # Collect (t3, H) pairs for all tournaments
    t3_vals = []
    h_vals = []
    for bits in range(2**m):
        adj = tournament_from_bits(n, bits)
        h = compute_H(n, adj)
        t3 = count_3cycles(n, adj)
        t3_vals.append(t3)
        h_vals.append(h)

    # Exact regression using Fraction
    N = len(h_vals)
    sum_t3 = sum(t3_vals)
    sum_h = sum(h_vals)
    sum_t3h = sum(t3_vals[i] * h_vals[i] for i in range(N))
    sum_t3_sq = sum(t**2 for t in t3_vals)

    mean_t3 = Fraction(sum_t3, N)
    mean_h = Fraction(sum_h, N)

    # b = Cov(t3, H) / Var(t3)
    cov = Fraction(sum_t3h, N) - mean_t3 * mean_h
    var_t3 = Fraction(sum_t3_sq, N) - mean_t3**2

    if var_t3 > 0:
        b = cov / var_t3
        a = mean_h - b * mean_t3
    else:
        b = Fraction(0)
        a = mean_h

    # R² = Cov²/(Var_t3 * Var_H)
    var_h = Fraction(sum(h**2 for h in h_vals), N) - mean_h**2
    r_sq = (cov * cov) / (var_t3 * var_h) if var_t3 > 0 and var_h > 0 else Fraction(0)

    print(f"\n  n={n}:")
    print(f"    Slope b = {b} = {float(b):.6f}")
    print(f"    Intercept a = {a} = {float(a):.6f}")
    print(f"    R² = {r_sq} = {float(r_sq):.6f}")
    print(f"    Mean(t3) = {mean_t3} = {float(mean_t3):.4f}")
    print(f"    Mean(H) = {mean_h} = {float(mean_h):.4f}")

    # Check: is the slope EXACTLY 6?
    if b == 6:
        print(f"    *** SLOPE IS EXACTLY 6 ***")
    elif b == 2:
        print(f"    *** SLOPE IS EXACTLY 2 ***")
    else:
        print(f"    Slope in simplest form: {b}")

    # Check residuals: is H - a - b*t3 always the same for fixed t3?
    residuals_by_t3 = defaultdict(list)
    for i in range(N):
        residual = Fraction(h_vals[i]) - a - b * t3_vals[i]
        residuals_by_t3[t3_vals[i]].append(residual)

    non_zero_residuals = False
    for t3 in sorted(residuals_by_t3.keys()):
        resids = residuals_by_t3[t3]
        unique_resids = sorted(set(resids))
        if len(unique_resids) > 1 or (len(unique_resids) == 1 and unique_resids[0] != 0):
            if not non_zero_residuals and unique_resids != [Fraction(0)]:
                non_zero_residuals = True
            if len(unique_resids) > 1:
                print(f"    t3={t3}: residuals = {[float(r) for r in unique_resids[:5]]}... ({len(unique_resids)} distinct)")

    if not non_zero_residuals:
        print(f"    *** H = a + b·t3 EXACTLY (zero residuals) ***")

# ======================================================================
# PART 2: WHY IS THE SLOPE 2 FOR n=3,4 AND 6 FOR n≥5?
# ======================================================================
print("\n" + "=" * 70)
print("PART 2: THE SLOPE PATTERN")
print("=" * 70)

print("""
  Observed slopes:
    n=3: b = 2
    n=4: b = 2
    n=5: b = ?
    n=6: b = 6

  The slope b = Cov(H, t3) / Var(t3).

  For n=3: t3 ∈ {0, 1}, H ∈ {1, 3}
    Perfectly correlated: H = 1 + 2·t3 exactly.

  For n=4: t3 ∈ {0, 1, 2}, H ∈ {1, 3, 5}
    Still perfect: H = 1 + 2·t3 exactly.

  For n=5: t3 determines H for t3 ≤ 3 but NOT for t3=4.
    The regression is no longer exact.

  Mean(H) = n!/2^{n-1}. Mean(t3) = C(n,3)/4.
  The ratio: Mean(H)/Mean(t3) = n!/(2^{n-1} · C(n,3)/4)
                                = 4n!/(2^{n-1} · C(n,3))
                                = 4·6·(n-3)!/(2^{n-1})
                                = 24·(n-3)!/2^{n-1}
""")

for n in range(3, 10):
    mean_h = Fraction(factorial(n), 2**(n-1))
    mean_t3 = Fraction(comb(n, 3), 4)
    if mean_t3 > 0:
        ratio = mean_h / mean_t3
    else:
        ratio = Fraction(0)
    print(f"  n={n}: Mean(H)/Mean(t3) = {ratio} = {float(ratio):.4f}")

# ======================================================================
# PART 3: PERMUTOHEDRON = TRANSITIVE FLIP GRAPH
# ======================================================================
print("\n" + "=" * 70)
print("PART 3: TRANSITIVE FLIP GRAPH = PERMUTOHEDRON")
print("=" * 70)

print("""
  THEOREM: The subgraph of Q_m induced on transitive tournaments
  is isomorphic to the permutohedron Π_n.

  PROOF SKETCH:
  1. Each transitive tournament T corresponds to a unique
     total ordering σ of the vertices: σ(1) beats σ(2) beats ... beats σ(n).
     This is a permutation σ ∈ S_n.

  2. Flipping arc (i,j) in T_σ changes the order of i and j.
     If i and j are ADJACENT in σ (σ^{-1}(i) = σ^{-1}(j) ± 1),
     then the result is another transitive tournament T_τ
     where τ is obtained by swapping i and j in σ.

  3. If i and j are NOT adjacent in σ, flipping (i,j) creates
     a 3-CYCLE (with the vertex between them in the ordering).
     So the result is NOT transitive.

  4. Therefore: T_σ and T_τ are one-flip-apart transitive iff
     σ and τ differ by an adjacent transposition.

  5. This is EXACTLY the definition of the permutohedron!
     The Cayley graph of S_n with generators {s_1, ..., s_{n-1}}
     where s_i = (i, i+1).

  COROLLARY: The transitive flip graph is (n-1)-regular.
  (Each permutation has n-1 adjacent transpositions.)

  The permutohedron has:
  - n! vertices
  - n!·(n-1)/2 edges
  - Diameter = C(n,2) = m (the maximum: reverse the ordering)
""")

# Verify the diameter claim
for n in range(3, 7):
    m = n * (n-1) // 2

    # Find transitive tournaments
    trans_bits = []
    for bits in range(2**m):
        adj = tournament_from_bits(n, bits)
        scores = sorted([sum(adj.get((i,j),0) for j in range(n) if j != i) for i in range(n)])
        if scores == list(range(n)):
            trans_bits.append(bits)

    # Find the two antipodal transitive tournaments
    # (the one with 0→1→2→...→n-1 and the reverse n-1→...→1→0)
    max_dist = 0
    for i in range(len(trans_bits)):
        for j in range(i+1, len(trans_bits)):
            d = bin(trans_bits[i] ^ trans_bits[j]).count('1')
            max_dist = max(max_dist, d)

    print(f"  n={n}: #trans={len(trans_bits)}, #edges={len(trans_bits)*(n-1)//2}, "
          f"diameter={max_dist} (expected C({n},2)={m})")

# ======================================================================
# PART 4: H ON THE PERMUTOHEDRON — ALWAYS 1
# ======================================================================
print("\n" + "=" * 70)
print("PART 4: H IS CONSTANT (=1) ON THE PERMUTOHEDRON")
print("=" * 70)

print("""
  Since all transitive tournaments have H = 1 (exactly one
  Hamiltonian path: the score ordering itself), H is CONSTANT
  on the permutohedron.

  This makes H a "distance function" FROM the permutohedron:
  H(T) = 1 iff T is transitive iff T is ON the permutohedron.

  The LEVEL SETS of H:
  - Level 1 = the permutohedron (n! vertices)
  - Level 3 = tournaments one "step above" in H
  - Level max = the "center" of the hypercube (most 3-cycles)

  This is a STRATIFICATION of Q_m by concentric "shells" around
  the permutohedron, where each shell has higher H.

  The permutohedron is the "equator" of a higher-dimensional
  analog of the sphere, with H playing the role of latitude.
""")

# ======================================================================
# PART 5: H = 1 + 2·t3 FOR n ≤ 4 (EXACT)
# ======================================================================
print("\n" + "=" * 70)
print("PART 5: WHY H = 1 + 2·t3 IS EXACT FOR n ≤ 4")
print("=" * 70)

print("""
  For n ≤ 4, the 3-cycle count t3 COMPLETELY determines H.

  WHY? At n=3: there is only one triangle {0,1,2}.
  Either it's a 3-cycle (t3=1, H=3) or not (t3=0, H=1).
  H = 1 + 2·t3 exactly.

  At n=4: there are C(4,3) = 4 triangles.
  The possible t3 values are 0, 1, 2.
  (t3=3 is impossible at n=4! Because 3 overlapping 3-cycles
   on 4 vertices force the fourth triangle to also be a 3-cycle.)
  Wait, let me check: can t3=3 at n=4?
""")

# Check: at n=4, what t3 values appear?
n = 4
m = 6
t3_h = defaultdict(set)
for bits in range(2**m):
    adj = tournament_from_bits(n, bits)
    h = compute_H(n, adj)
    t3 = count_3cycles(n, adj)
    t3_h[t3].add(h)

print(f"  n=4: t3 → H mapping:")
for t3 in sorted(t3_h.keys()):
    print(f"    t3={t3}: H = {sorted(t3_h[t3])}")

print("""
  So at n=4: t3 ∈ {0, 1, 2} and H = 1 + 2·t3 exactly.
  t3 cannot be 3 at n=4 (one can verify: having 3 out of 4
  triangles as 3-cycles forces the fourth to also be a 3-cycle,
  giving t3=4, which is the regular tournament — but wait,
  n=4 has C(4,3)=4 triangles).
""")

# Check: at n=4, can we have t3=3?
print(f"  Checking: t3 values that occur at n=4:")
all_t3 = Counter()
for bits in range(2**m):
    adj = tournament_from_bits(n, bits)
    t3 = count_3cycles(n, adj)
    all_t3[t3] += 1

for t3 in sorted(all_t3.keys()):
    print(f"    t3={t3}: {all_t3[t3]} tournaments")

# ======================================================================
# PART 6: THE FORMULA H = f(t3, ...) FOR n ≥ 5
# ======================================================================
print("\n" + "=" * 70)
print("PART 6: WHAT ADDITIONAL INVARIANT IS NEEDED FOR n ≥ 5?")
print("=" * 70)

# At n=5, t3=4 has H ∈ {11, 13, 15}. What distinguishes these?
n = 5
m = 10

tournaments_t3_4 = []
for bits in range(2**m):
    adj = tournament_from_bits(n, bits)
    t3 = count_3cycles(n, adj)
    if t3 == 4:
        h = compute_H(n, adj)
        scores = tuple(sorted(sum(adj.get((i,j),0) for j in range(n) if j != i) for i in range(n)))
        tournaments_t3_4.append((bits, h, scores))

print(f"\n  n=5, t3=4: {len(tournaments_t3_4)} tournaments")
h_by_scores = defaultdict(list)
for bits, h, scores in tournaments_t3_4:
    h_by_scores[(scores, h)] = h_by_scores.get((scores, h), 0) + 1

# Count by (scores, H)
score_h_dist = defaultdict(lambda: Counter())
for bits, h, scores in tournaments_t3_4:
    score_h_dist[scores][h] += 1

for scores in sorted(score_h_dist.keys()):
    h_dist = score_h_dist[scores]
    print(f"    Scores {scores}: H distribution = {dict(sorted(h_dist.items()))}")

# What about 5-cycles? At n=5, there are directed 5-cycles.
# Maybe t5 (number of 5-cycles) discriminates.
print(f"\n  Checking 5-cycle counts for t3=4 tournaments:")

def count_5cycles(n, adj):
    """Count directed 5-cycles."""
    count = 0
    from itertools import permutations
    # A 5-cycle on all 5 vertices: (v0→v1→v2→v3→v4→v0)
    # There are (5-1)!/2 = 12 directed 5-cycles on 5 labeled vertices
    for perm in permutations(range(n)):
        is_cycle = all(adj.get((perm[i], perm[(i+1) % n]), 0) for i in range(n))
        if is_cycle:
            count += 1
    # Each 5-cycle is counted n times (rotations) and not counted for reverse
    return count // n

t5_h = defaultdict(list)
for bits, h, scores in tournaments_t3_4:
    adj = tournament_from_bits(n, bits)
    t5 = count_5cycles(n, adj)
    t5_h[t5].append(h)

for t5 in sorted(t5_h.keys()):
    h_vals = t5_h[t5]
    h_counter = Counter(h_vals)
    print(f"    t5={t5}: {dict(sorted(h_counter.items()))}")

# ======================================================================
# PART 7: THE EXACT FORMULA H = 1 + 2t3 + ε(T) WHERE ε DEPENDS ON...
# ======================================================================
print("\n" + "=" * 70)
print("PART 7: THE RESIDUAL ε(T) = H - 1 - 2t3")
print("=" * 70)

for n in range(3, 7):
    m = n * (n-1) // 2

    residuals = Counter()
    for bits in range(2**m):
        adj = tournament_from_bits(n, bits)
        h = compute_H(n, adj)
        t3 = count_3cycles(n, adj)
        eps = h - 1 - 2 * t3
        residuals[eps] += 1

    print(f"\n  n={n}: Residual ε = H - 1 - 2t3:")
    for eps in sorted(residuals.keys()):
        print(f"    ε={eps:3d}: {residuals[eps]:5d} tournaments ({residuals[eps]/2**m:.4f})")

    # Is ε always ≥ 0?
    min_eps = min(residuals.keys())
    max_eps = max(residuals.keys())
    mean_eps = sum(eps * count for eps, count in residuals.items()) / 2**m
    print(f"    Range: [{min_eps}, {max_eps}], Mean(ε) = {mean_eps:.4f}")

    # ε = 0 count:
    zero_count = residuals.get(0, 0)
    print(f"    ε = 0: {zero_count}/{2**m} ({zero_count/2**m:.4f})")

# ======================================================================
# PART 8: WHAT IS ε?
# ======================================================================
print("\n" + "=" * 70)
print("PART 8: ε AS A FUNCTION OF HIGHER CYCLE COUNTS")
print("=" * 70)

# For n=6, compute (H, t3, t5, score_seq) and see what determines H
n = 6
m = 15

# We can't easily count 5-cycles at n=6 (too expensive). But we can
# use the score sequence as a proxy.

# At n=6: the formula H = 1 + 2t3 gives residual ε.
# Is ε related to the score sequence? Or to some other invariant?

# Let's group by (t3, score_seq) and check if H is determined
from itertools import permutations as perms

print(f"\n  n={n}: Grouping by (t3, score_sequence) → H:")
group = defaultdict(lambda: Counter())

for bits in range(2**m):
    adj = tournament_from_bits(n, bits)
    h = compute_H(n, adj)
    t3 = count_3cycles(n, adj)
    scores = tuple(sorted(sum(adj.get((i,j),0) for j in range(n) if j != i) for i in range(n)))
    group[(t3, scores)][h] += 1

non_determined = 0
total_groups = 0
for key in sorted(group.keys()):
    h_dist = group[key]
    total_groups += 1
    if len(h_dist) > 1:
        non_determined += 1
        if non_determined <= 10:
            print(f"    t3={key[0]}, scores={key[1]}: H={dict(sorted(h_dist.items()))}")

print(f"\n  Total groups: {total_groups}")
print(f"  Non-determined (multiple H values): {non_determined}")
print(f"  Fraction determined: {1 - non_determined/total_groups:.4f}")

print("\n" + "=" * 70)
print("DONE — H ≈ 6·t3 REGRESSION AND PERMUTOHEDRON")
print("=" * 70)
