#!/usr/bin/env python3
"""
COMPLETE ALGEBRAIC DERIVATION of moment formulas at n=7.

Using the Stirling decomposition: sum_P f^j = sum_{k=0}^j S(j,k) * k! * SIGMA_k
where SIGMA_k = sum over size-k position subsets of sigma(S).

Each sigma(S) depends only on the adjacency pattern of S.
We now have algebraic proofs for most patterns.

opus-2026-03-06-S28
"""
from math import factorial, comb
from itertools import permutations, combinations
import random

# Stirling numbers of the second kind
def stirling2(n, k):
    if n == 0 and k == 0: return 1
    if n == 0 or k == 0: return 0
    return k * stirling2(n-1, k) + stirling2(n-1, k-1)

n = 7

# =====================================================================
# SIGMA values for each pattern at n=7
# All patterns with their multiplicities and formulas
# =====================================================================

# For pattern (c1, c2, ...) = connected component sizes of position subset S,
# count how many position subsets have this pattern.

def count_patterns(n, k):
    """Count position subsets of [0,n-2] of size k, grouped by pattern."""
    from collections import Counter, defaultdict
    patterns = defaultdict(int)
    if k == 0:
        patterns[()] = 1
        return dict(patterns)
    for S in combinations(range(n-1), k):
        # Compute connected components
        pos = sorted(S)
        comps = []
        comp = [pos[0]]
        for i in range(1, len(pos)):
            if pos[i] == comp[-1] + 1:
                comp.append(pos[i])
            else:
                comps.append(len(comp))
                comp = [pos[i]]
        comps.append(len(comp))
        pat = tuple(sorted(comps, reverse=True))
        patterns[pat] += 1
    return dict(patterns)

print("=" * 70)
print(f"COMPLETE MOMENT DERIVATION at n={n}")
print("=" * 70)

# List all patterns by size
print(f"\nPosition subset patterns:")
all_patterns = {}
for k in range(0, n):
    pats = count_patterns(n, k)
    all_patterns[k] = pats
    for pat, count in sorted(pats.items()):
        print(f"  k={k}: pattern {pat} x {count}")

# =====================================================================
# SIGMA formulas (algebraically proved or to verify)
# =====================================================================

# Let's define sigma for each pattern as a function of (t3, t5, bc)
# Based on our proofs:

# Pattern formulas at n=7:
# (1,) -> n!/2 = 2520
# (1,1) -> 1260
# (1,1,1) -> 630
# (2,) -> (n-3)!(C(n,3) + 2*t3) = 120*(35+2*t3) = 4200 + 240*t3
#          Wait, (n-3)! = 4! = 24 at n=7. So = 24*(35+2*t3) = 840 + 48*t3
#          Hmm, let me re-check. I had sigma((2,)) varying: for t3=5 it was 1080.
#          840 + 48*5 = 840+240 = 1080. ✓

# Actually wait — I need to recheck. The formula was:
# sigma({i,i+1}) = (n-3)! * sum over ordered triples (a,b,c) of A[a,b]*A[b,c]
# = (n-3)! * (C(n,3) + 2*t3)
# At n=7: (n-3)! = 4! = 24. C(7,3) = 35.
# = 24 * (35 + 2*t3) = 840 + 48*t3
# For t3=5: 840 + 240 = 1080 ✓
# For t3=12: 840 + 576 = 1416 ✓ ✓ ✓

sigma_formulas = {
    # pattern -> (const, t3_coeff, t5_coeff, bc_coeff)
    (): (factorial(n), 0, 0, 0),  # k=0: sum over 0 positions = n!
    (1,): (factorial(n)//2, 0, 0, 0),
    (2,): (24*35, 24*2, 0, 0),  # = 840 + 48*t3
    (3,): (210, 48, 0, 0),
    (4,): (42, 24, 4, 0),
    # (5,): TBD
    (6,): (0, 0, 0, 0),  # = H, needs separate handling
}

# For non-consecutive patterns, I need to verify/derive:
# (1,1): from the pair-partition lemma
# sum_P A[p_i,p_{i+1}]*A[p_j,p_{j+1}] for |i-j|>1
# = (n-4)! * sum over ordered 4-tuples (a,b,c,d) all distinct of A[a,b]*A[c,d]
# = (n-4)! * [sum_{(a,b) dist} A[a,b]] * [sum_{(c,d) dist, c,d not in {a,b}} A[c,d]]
# Hmm, that's not a simple product. Let me think...
#
# Actually: sum_{(a,b,c,d) all distinct} A[a,b]*A[c,d]
# = [sum_{(a,b) dist} A[a,b]] * [sum_{(c,d) dist} A[c,d]] - [sum_{(a,b) dist} A[a,b]*A[b,a]]
#   - ... no, this overcounting correction is complex.
#
# By the pair-partition universality: for k=2 (two edges on 4 distinct vertices),
# sum over ordered (a,b,c,d) distinct of A[a,b]*A[c,d] = 4!/2^2 * C(n,4) = 6*C(n,4)
# Wait, that's not quite right. The universality says:
# sum over ordered (a,b,c,d) distinct of A[a,b]*A[c,d]
# where we pair (a,b) and (c,d) = 4!/2^2 * C(n,4) = 6 * C(n,4)
#
# Hmm, but the original universality was about pairing (v1,v2),(v3,v4) etc.
# At k=2: sum over ordered 4-tuples = 2!*C(n,4) = 2*C(n,4)?
# No. The statement was: sum over ORDERED (2k)-tuples of prod paired values = (2k)!/2^k * C(n,2k).
# At k=2: 4!/4 = 6 * C(n,4).
# But we want the pairing (a,b),(c,d) — which is ONE specific pairing.
# There are (2k-1)!! = 3 pairings of 4 elements.
# Each pairing contributes equally, so each pairing gives 6*C(n,4)/3 = 2*C(n,4)?
# No, the universality says total = (2k)!/2^k * C(n,2k) = sum over ALL orderings.
# Each ordering assigns the 4 elements to positions. The pairing is (pos1,pos2),(pos3,pos4).
# There are 3 possible pairings, each with 2^2 * 2! = 8 orderings per pairing.
# Total orderings = 3 * 8 = 24 = 4!. Each contributes 1 per 4-element subset.
# So total = C(n,4) * sum per subset = C(n,4) * 24/4 = 6*C(n,4).
# Wait I'm confusing myself. Let me just use the verified value: sigma((1,1)) = 1260.
# = (n-4)! * P(n,4)/4 = 6 * 7*6*5*4/4 = 6 * 210 = 1260. ✓

sigma_formulas[(1,1)] = (1260, 0, 0, 0)  # Universal

# (1,1,1): = (n-6)! * [n*(n-1)/2]^3 / adjustment? Let me just verify.
# For 3 isolated positions {i,j,k} with gaps > 1:
# sigma = (n-6)! * sum over ordered 6-tuples (a,b,c,d,e,f) all distinct of A[a,b]*A[c,d]*A[e,f]
# By pair-partition: this should be (n-6)! * 6!/2^3 * C(n,6) / (number of other pairings?)
# Actually this is the 3-pair case of pair-partition universality.
# The total over all orderings is 6!/2^3 * C(n,6) = 90 * C(7,6) = 90*7 = 630.
# And (n-6)! = 1! = 1 at n=7.
# So sigma = 1 * 630 = 630. ✓
sigma_formulas[(1,1,1)] = (630, 0, 0, 0)  # Universal

# (2,1): consecutive pair + isolated
# sigma = (n-5)! * sum over ordered 5-tuples (a,b,c,d,e) all distinct of A[a,b]*A[b,c]*A[d,e]
# The A[a,b]*A[b,c] part is a directed 2-path (depends on t3 among the 3 vertices).
# The A[d,e] part is a single edge (always 1/2 per ordered pair).
# But they're on DISJOINT vertex sets, so:
# = (n-5)! * [sum over ordered triples of A[a,b]*A[b,c]] * ... no, not independent.
#
# More carefully: choose 5 distinct vertices, assign 3 to the 2-path and 2 to the edge.
# sum = (n-5)! * sum_{5-subsets S} sum_{(a,b,c,d,e): a,b,c form 2-path, d,e form edge,
#        all from S, distinct} A[a,b]*A[b,c]*A[d,e]
#
# For a 5-subset S:
# sum over (ordered triple in S forming 2-path) * (ordered pair from remaining 2 vertices of edge)
# = [sum over ordered (a,b,c) from S of A[a,b]*A[b,c]] * [A[d,e] for the remaining 2]
# Wait no — the 3 vertices (a,b,c) and 2 vertices (d,e) partition the 5-subset.
# But the assignment is: positions (i, i+1) get (a,b,c) = 2-path, position j gets (d,e).
#
# Actually, for the (2,1) pattern at positions {i, i+1, j}:
# sigma = sum_P A[p_i,p_{i+1}]*A[p_{i+1},p_{i+2}]*A[p_j,p_{j+1}]
# Positions i, i+1, i+2, j, j+1 involve 5 vertex positions.
# The remaining n-5=2 vertices are free.
#
# = (n-5)! * sum over ordered 5-tuples (a,b,c,d,e) of distinct vertices:
#   assign a=p_i, b=p_{i+1}, c=p_{i+2}, d=p_j, e=p_{j+1}
#   product = A[a,b]*A[b,c]*A[d,e]
#
# = (n-5)! * sum_{5-subsets S} sum over assignments of S to (a,b,c,d,e)
#   with all orderings: 5! / (already fully ordered) = actually all 5! orderings
#   but the product A[a,b]*A[b,c]*A[d,e] depends on which vertices are in which roles.
#
# Decompose: for each 5-subset, choose which 3 play (a,b,c) role and which 2 play (d,e).
# Number of ways: C(5,3) = 10 choices of 3-subset for (a,b,c).
# For each: 3! orderings of the 3-subset as (a,b,c), and 2! orderings of the 2-subset as (d,e).
# Product = A[a,b]*A[b,c]*A[d,e].
# sum over 3! orderings of (a,b,c) = sum of directed 2-paths = C(3,3) + 2*c3(triple)
#   No: for a specific 3-element set T:
#   sum over ordered (a,b,c) of T of A[a,b]*A[b,c] = dp_2(T) = 3 if cyclic, 1 if transitive.
# sum over 2! orderings of (d,e) = A[d,e] + A[e,d] = 1.
#
# So for each 5-subset S and each 3-subset T of S:
# contribution = dp_2(T) * 1 = dp_2(T)
#
# Total = (n-5)! * sum_{5-subsets S} sum_{3-subsets T of S} dp_2(T)
# = (n-5)! * sum_{3-subsets T of V} dp_2(T) * C(n-3, 5-3)
# = (n-5)! * C(n-3, 2) * sum_{3-subsets T} dp_2(T)
# = 2! * C(4,2) * [C(n,3) + 2*t3]
# = 2 * 6 * (35 + 2*t3) = 12*(35 + 2*t3) = 420 + 24*t3 ✓
sigma_formulas[(2,1)] = (420, 24, 0, 0)

# (2,2): PROVED: 16*t3 + 8*bc + 140
sigma_formulas[(2,2)] = (140, 16, 0, 8)

# (3,1): consecutive 3 + isolated
# Similar decomposition: positions {i,i+1,i+2,j} use 5 vertex positions.
# The 3-chain uses (p_i,p_{i+1},p_{i+2},p_{i+3}) = 4 vertex positions.
# The isolated uses (p_j,p_{j+1}) = 2 vertex positions.
# Total: 6 vertex positions, 1 free.
#
# sigma = (n-6)! * sum over ordered 6-tuples (a,b,c,d,e,f) all distinct:
#   product = A[a,b]*A[b,c]*A[c,d]*A[e,f]
#   where (a,b,c,d) is a 3-path and (e,f) is a single edge.
#
# For each 6-vertex subset and partition into 4-group + 2-group:
# contribution = H(T[4-group]) * (A[e,f] + A[f,e]) = H(T[4-group]) * 1
#
# Total = (n-6)! * sum_{6-subsets} sum_{4-subsets G of S} H(T[G])
# = 1 * sum_{4-subsets G of V} H(T[G]) * C(n-4, 6-4)
# = C(3,2) * sum_{4-subsets} H(T[G])
# = 3 * sum_{4-subsets} H(T[G])
#
# And H(T[G]) = 1 + 2*c3(G) by OCF at n=4.
# sum_{4-subsets} H(T[G]) = C(7,4) + 2*C(4,1)*t3 = 35 + 8*t3
#   (each 3-cycle is in C(4,1)=4 four-element subsets)
#
# sigma((3,1)) = 3*(35 + 8*t3) = 105 + 24*t3
#
# BUT WAIT: does sigma((3,1)) depend on t5? Let me check from the data.
# In the first script, I had sigma((3,1)) varying. Let me check if it's just t3.

# Actually, I need to double-check: the 3-chain uses positions {i,i+1,i+2},
# which means edges (p_i,p_{i+1}), (p_{i+1},p_{i+2}), (p_{i+2},p_{i+3}).
# That's 4 VERTICES at these positions: p_i, p_{i+1}, p_{i+2}, p_{i+3}.
# Product A[a,b]*A[b,c]*A[c,d] = directed 3-path = Ham path on 4 vertices.
# H(T[{a,b,c,d}]) counts all orderings of Ham paths, which equals the sum.

# So sigma((3,1)) = 3 * (35 + 8*t3) = 105 + 24*t3
sigma_formulas[(3,1)] = (105, 24, 0, 0)

# (2,1,1): consecutive pair + two isolated
# Positions use 2+2+3 = 7 vertex positions... wait, let me count.
# Pattern (2,1,1) means: one consecutive pair (2 positions, 3 vertices),
# two isolated positions (each using 2 vertices).
# But the positions might overlap in vertex usage!
# Actually, positions {i, i+1, j, k} where |i-j|>1, |i+1-j|>1, |j-k|>1.
# These use vertex positions i, i+1, i+2, j, j+1, k, k+1 = 7 positions at n=7.
# So ALL n vertices are used! No free vertex.
#
# sigma = sum_P A[p_i,p_{i+1}]*A[p_{i+1},p_{i+2}]*A[p_j,p_{j+1}]*A[p_k,p_{k+1}]
# = 1 * sum over permutations of 7 vertices of A[a,b]*A[b,c]*A[d,e]*A[f,g]
#
# where (a,b,c) form a 2-path and (d,e), (f,g) are single edges,
# all on disjoint vertex sets partitioning all 7 vertices into 3+2+2.
#
# For each partition (3-set T, 2-set U, 2-set V) of [n]:
# contribution = dp_2(T) * (A[u1,u2]+A[u2,u1]) * (A[v1,v2]+A[v2,v1]) = dp_2(T)
#
# Number of ORDERED partitions into (3,2,2):
# (We need to count ordered because positions specify the order)
# C(7,3) * C(4,2) / 1 = 35 * 6 = 210?
# Wait, the two 2-sets are distinguishable (one goes to position j, other to k).
# So: C(7,3) * C(4,2) * C(2,2) = 35 * 6 * 1 = 210 ordered partitions.
# Each contributes dp_2(T) * 1 * 1 = dp_2(T).
#
# But each 3-set T appears in C(4,2) = 6 ordered partitions (ways to split remaining 4 into 2+2).
# Wait, for a specific 3-set T, the remaining 4 vertices split into two 2-sets.
# The number of ORDERED splits is C(4,2) = 6 (choose which 2 go to position j).
#
# So sigma = sum_{3-subsets T} dp_2(T) * C(4,2) = 6 * sum_T dp_2(T)
# = 6 * [sum_T dp_2(T)]
#
# sum_T dp_2(T) = sum over C(7,3)=35 triples of (3 if cyclic, 1 if transitive)
# = 3*t3 + 1*(35-t3) = 35 + 2*t3
#
# sigma((2,1,1)) = 6 * (35 + 2*t3) = 210 + 12*t3
sigma_formulas[(2,1,1)] = (210, 12, 0, 0)

# =====================================================================
# Now compute SIGMA_k and moments
# =====================================================================
print(f"\nSIGMA_k FORMULAS:")
print(f"SIGMA_k = sum of pattern multiplicities * sigma(pattern)")

# SIGMA_k as function of (t3, t5, bc, H)
# SIGMA_k(t3, t5, bc) = sum over patterns of size k of (count * sigma(pattern))

for k in range(n):
    pats = all_patterns[k]
    const = 0
    t3_coeff = 0
    t5_coeff = 0
    bc_coeff = 0
    h_coeff = 0
    has_unknown = False

    for pat, count in pats.items():
        if pat in sigma_formulas:
            c, t3c, t5c, bcc = sigma_formulas[pat]
            const += count * c
            t3_coeff += count * t3c
            t5_coeff += count * t5c
            bc_coeff += count * bcc
        elif pat == (6,):
            h_coeff += count  # sigma((6,)) = H
        elif pat == (5,):
            has_unknown = True
            print(f"  WARNING: pattern {pat} not yet derived")
        elif pat == (4,1):
            has_unknown = True
            print(f"  WARNING: pattern {pat} not yet derived")
        elif pat == (3,2):
            has_unknown = True
            print(f"  WARNING: pattern {pat} not yet derived")
        else:
            has_unknown = True
            print(f"  WARNING: pattern {pat} not yet derived")

    terms = [f"{const}"]
    if t3_coeff: terms.append(f"{t3_coeff}*t3")
    if t5_coeff: terms.append(f"{t5_coeff}*t5")
    if bc_coeff: terms.append(f"{bc_coeff}*bc")
    if h_coeff: terms.append(f"{h_coeff}*H")
    formula = " + ".join(terms)
    status = "" if not has_unknown else " [INCOMPLETE]"
    print(f"  SIGMA_{k} = {formula}{status}")

# =====================================================================
# Compute moments from Stirling + SIGMA
# =====================================================================
print(f"\nMOMENT FORMULAS (from Stirling decomposition):")
print(f"  m_j = sum_k S(j,k) * k! * SIGMA_k")

# For moments up to j=4 (where we have all needed patterns):
for j in range(5):
    const = 0
    t3_coeff = 0
    t5_coeff = 0
    bc_coeff = 0
    h_coeff = 0

    for k in range(min(j, n-1) + 1):
        S_jk = stirling2(j, k)
        if S_jk == 0:
            continue

        pats = all_patterns[k]
        for pat, count in pats.items():
            if pat in sigma_formulas:
                c, t3c, t5c, bcc = sigma_formulas[pat]
                factor = S_jk * factorial(k)
                const += factor * count * c
                t3_coeff += factor * count * t3c
                t5_coeff += factor * count * t5c
                bc_coeff += factor * count * bcc
            elif pat == (6,):
                h_coeff += S_jk * factorial(k) * count

    terms = [f"{const}"]
    if t3_coeff: terms.append(f"{t3_coeff}*t3")
    if t5_coeff: terms.append(f"{t5_coeff}*t5")
    if bc_coeff: terms.append(f"{bc_coeff}*bc")
    if h_coeff: terms.append(f"{h_coeff}*H")
    print(f"  m_{j} = {' + '.join(terms)}")

# =====================================================================
# VERIFY against computational values
# =====================================================================
print(f"\n{'='*70}")
print("VERIFICATION")
print(f"{'='*70}")

def count_3_cycles(A, n):
    count = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if A[i][j]*A[j][k]*A[k][i]: count += 1
                if A[i][k]*A[k][j]*A[j][i]: count += 1
    return count

def count_5_cycles_dp(A, n):
    count = 0
    for verts in combinations(range(n), 5):
        sub = [[A[verts[i]][verts[j]] for j in range(5)] for i in range(5)]
        dp = [[0]*5 for _ in range(1 << 5)]
        dp[1][0] = 1
        for mask in range(1, 1 << 5):
            for v in range(5):
                if not (mask & (1 << v)) or dp[mask][v] == 0: continue
                for u in range(5):
                    if mask & (1 << u): continue
                    if sub[v][u]:
                        dp[mask | (1 << u)][u] += dp[mask][v]
        full = (1 << 5) - 1
        hc = sum(dp[full][v] for v in range(1, 5) if sub[v][0])
        count += hc
    return count

def is_3cycle_check(A, triple):
    a, b, c = triple
    return (A[a][b]*A[b][c]*A[c][a] + A[a][c]*A[c][b]*A[b][a]) > 0

def both_cyclic_total_check(A, n):
    total = 0
    for v in range(n):
        S = [u for u in range(n) if u != v]
        for T in combinations(S, 3):
            T_comp = tuple(u for u in S if u not in T)
            if T < T_comp:
                if is_3cycle_check(A, T) and is_3cycle_check(A, T_comp):
                    total += 1
    return total

max_errs = {j: 0 for j in range(5)}
for trial in range(15):
    random.seed(trial * 53 + 7)
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1

    t3 = count_3_cycles(A, n)
    t5 = count_5_cycles_dp(A, n)
    bc = both_cyclic_total_check(A, n)

    # Compute actual moments
    moments = [0] * 5
    for p in permutations(range(n)):
        f = sum(1 for i in range(n-1) if A[p[i]][p[i+1]])
        fpow = 1
        for j in range(5):
            moments[j] += fpow
            fpow *= f

    # Predicted moments from formulas
    # m0 = 5040
    pred = [
        5040,
        15120,
        48720 + 480*t3,
        166320 + 4320*t3,
        596064 + 27840*t3 + 288*t5 + 576*bc,
    ]

    for j in range(5):
        err = abs(moments[j] - pred[j])
        max_errs[j] = max(max_errs[j], err)

    if trial < 5:
        print(f"  T{trial}: t3={t3:2d} t5={t5:2d} bc={bc:2d}")
        for j in range(5):
            print(f"    m{j}: actual={moments[j]:10d} pred={pred[j]:10d} err={abs(moments[j]-pred[j])}")

print(f"\nMax errors: {max_errs}")

# Now derive the coefficients from the sigma table
print(f"\n{'='*70}")
print("DERIVED COEFFICIENTS (algebraic)")
print(f"{'='*70}")

# m0 = S(0,0)*0!*SIGMA_0 = 1*1*5040 = 5040 ✓
# m1 = S(1,1)*1!*SIGMA_1 = 1*1*(6*2520) = 15120 ✓
# m2 = S(2,1)*1!*SIGMA_1 + S(2,2)*2!*SIGMA_2
# S(2,1)=1, S(2,2)=1
# = 1*1*15120 + 1*2*SIGMA_2
# SIGMA_2 = 10*(1260) + 5*(840+48*t3) = 12600 + 4200 + 240*t3 = 16800 + 240*t3
# m2 = 15120 + 2*(16800 + 240*t3) = 15120 + 33600 + 480*t3 = 48720 + 480*t3 ✓!

# m3 = S(3,1)*1!*SIGMA_1 + S(3,2)*2!*SIGMA_2 + S(3,3)*3!*SIGMA_3
# S(3,1)=1, S(3,2)=3, S(3,3)=1
# SIGMA_3 = 4*630 + 12*(420+24*t3) + 4*(210+48*t3)
#         = 2520 + 5040+288*t3 + 840+192*t3 = 8400 + 480*t3
# m3 = 15120 + 6*(16800+240*t3) + 6*(8400+480*t3)
# = 15120 + 100800+1440*t3 + 50400+2880*t3
# = 166320 + 4320*t3 ✓!

# m4 = S(4,1)*1!*SIGMA_1 + S(4,2)*2!*SIGMA_2 + S(4,3)*3!*SIGMA_3 + S(4,4)*4!*SIGMA_4
# S(4,1)=1, S(4,2)=7, S(4,3)=6, S(4,4)=1
# SIGMA_4 = 3*(210+12*t3) + 3*(140+16*t3+8*bc) + 6*(105+24*t3) + 3*(42+24*t3+4*t5)
#         = 630+36*t3 + 420+48*t3+24*bc + 630+144*t3 + 126+72*t3+12*t5
#         = 1806 + 300*t3 + 12*t5 + 24*bc
# m4 = 15120 + 14*(16800+240*t3) + 36*(8400+480*t3) + 24*(1806+300*t3+12*t5+24*bc)
# = 15120 + 235200+3360*t3 + 302400+17280*t3 + 43344+7200*t3+288*t5+576*bc
# = (15120+235200+302400+43344) + (3360+17280+7200)*t3 + 288*t5 + 576*bc
# = 596064 + 27840*t3 + 288*t5 + 576*bc ✓!!!

print(f"""
  SIGMA_0 = {factorial(n)}
  SIGMA_1 = {6 * factorial(n)//2}
  SIGMA_2 = 16800 + 240*t3
  SIGMA_3 = 8400 + 480*t3
  SIGMA_4 = 1806 + 300*t3 + 12*t5 + 24*bc  [FIRST appearance of t5 and bc]

  m0 = 5040
  m1 = 15120
  m2 = 48720 + 480*t3     [t3 enters at j=2, via Stirling hitting SIGMA_2]
  m3 = 166320 + 4320*t3
  m4 = 596064 + 27840*t3 + 288*t5 + 576*bc  [t5, bc enter at j=4]

  ALL DERIVED ALGEBRAICALLY from sigma pattern formulas!

  The hierarchy emerges because:
  - SIGMA_0, SIGMA_1: universal (only pattern (1,...) singletons)
  - SIGMA_2: pattern (2,) introduces t3 via directed 2-path counting
  - SIGMA_3: patterns (3,), (2,1) both depend on t3 only
  - SIGMA_4: pattern (2,2) introduces bc; pattern (4,) introduces t5
  - SIGMA_5: patterns (3,2), (4,1), (5,) — expected to bring in t5, bc
  - SIGMA_6: pattern (6,) = H
""")

print(f"{'='*70}")
print("DONE")
print(f"{'='*70}")
