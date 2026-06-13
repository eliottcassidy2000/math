import os; os.environ['PYTHONIOENCODING'] = 'utf-8'
"""
c5_formula_regular.py
kind-pasteur-2026-03-07-S39b

Investigate closed-form formula for c_5 in regular tournaments.

At n=7, c_5 takes 3 values {28, 36, 42} for the 3 iso-classes of regular
tournaments. Can c_5 be expressed in terms of:
1. tr(A^5) where A is the adjacency matrix?
2. Spectral invariants of A?
3. Number of certain substructures?

Savchenko (arXiv:2403.07629): c_k for DRT satisfies closed forms.
The question is whether a UNIVERSAL formula exists for all regular tournaments.
"""

import sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '03-artifacts', 'code'))
from tournament_lib import tournament_from_bits
from itertools import combinations
from collections import defaultdict
from math import comb


def count_directed_k_cycles_dp(T, k):
    """Count directed k-cycles using subset DP."""
    n = len(T)
    if k > n or k < 3:
        return 0
    count = 0
    for verts in combinations(range(n), k):
        v = list(verts)
        dp = [[0] * k for _ in range(1 << k)]
        dp[1][0] = 1
        for mask in range(1, 1 << k):
            for last in range(k):
                if dp[mask][last] == 0 or not (mask & (1 << last)):
                    continue
                for nxt in range(1, k):
                    if mask & (1 << nxt):
                        continue
                    if T[v[last]][v[nxt]]:
                        dp[mask | (1 << nxt)][nxt] += dp[mask][last]
        full = (1 << k) - 1
        for last in range(1, k):
            if T[v[last]][v[0]]:
                count += dp[full][last]
    return count


def matrix_power_trace(T, k):
    """Compute tr(A^k) where A is the adjacency matrix of T."""
    n = len(T)
    # Compute A^k via repeated multiplication
    A = [row[:] for row in T]
    Ak = [[int(i == j) for j in range(n)] for i in range(n)]
    for _ in range(k):
        new = [[0]*n for _ in range(n)]
        for i in range(n):
            for j in range(n):
                for l in range(n):
                    new[i][j] += Ak[i][l] * A[l][j]
        Ak = new
    return sum(Ak[i][i] for i in range(n))


def matrix_product_trace(T, powers):
    """Compute tr(A^p1 * A^p2 * ...) = tr(A^(p1+p2+...))."""
    total_power = sum(powers)
    return matrix_power_trace(T, total_power)


def count_homomorphic_cycles(T, k):
    """Count closed walks of length k (= tr(A^k)).
    This counts homomorphic cycles, not just simple ones."""
    return matrix_power_trace(T, k)


# ============================================================
# Analysis at n=7
# ============================================================
print("=" * 70)
print("c_5 FORMULA FOR REGULAR TOURNAMENTS")
print("=" * 70)

n = 7
m = n * (n - 1) // 2

# Collect data for all regular n=7 tournaments
data = []
for bits in range(1 << m):
    T = tournament_from_bits(n, bits)
    scores = [sum(T[i]) for i in range(n)]
    if not all(s == 3 for s in scores):
        continue

    c3 = count_directed_k_cycles_dp(T, 3)
    c5 = count_directed_k_cycles_dp(T, 5)
    c7 = count_directed_k_cycles_dp(T, 7)

    tr3 = matrix_power_trace(T, 3)
    tr5 = matrix_power_trace(T, 5)

    # Skew adjacency: S[i][j] = T[i][j] - T[j][i]
    S = [[T[i][j] - T[j][i] for j in range(n)] for i in range(n)]
    # tr(S^k) for odd k
    Sk = [[int(i == j) for j in range(n)] for i in range(n)]
    for _ in range(5):
        new = [[0]*n for _ in range(n)]
        for i in range(n):
            for j in range(n):
                for l in range(n):
                    new[i][j] += Sk[i][l] * S[l][j]
        Sk = new
    trS5 = sum(Sk[i][i] for i in range(n))

    data.append({
        'bits': bits, 'c3': c3, 'c5': c5, 'c7': c7,
        'tr3': tr3, 'tr5': tr5, 'trS5': trS5,
    })

# Group by c5
c5_groups = defaultdict(list)
for d in data:
    c5_groups[d['c5']].append(d)

print(f"\nn=7: {len(data)} regular tournaments, 3 classes by c5:")
for c5val in sorted(c5_groups.keys()):
    group = c5_groups[c5val]
    tr3_vals = set(d['tr3'] for d in group)
    tr5_vals = set(d['tr5'] for d in group)
    trS5_vals = set(d['trS5'] for d in group)
    c7_vals = set(d['c7'] for d in group)
    print(f"\n  c_5 = {c5val}: {len(group)} tournaments")
    print(f"    c_3 = {group[0]['c3']}, c_7 in {sorted(c7_vals)}")
    print(f"    tr(A^3) in {sorted(tr3_vals)}")
    print(f"    tr(A^5) in {sorted(tr5_vals)}")
    print(f"    tr(S^5) in {sorted(trS5_vals)} (S = skew adjacency)")

# Check if tr(A^5) or tr(S^5) uniquely determines c_5
print(f"\n--- Relationship: tr(A^5) vs c_5 ---")
tr5_to_c5 = defaultdict(set)
for d in data:
    tr5_to_c5[d['tr5']].add(d['c5'])
for tr5 in sorted(tr5_to_c5.keys()):
    c5s = tr5_to_c5[tr5]
    print(f"  tr(A^5) = {tr5} -> c_5 in {sorted(c5s)}")

print(f"\n--- Relationship: tr(S^5) vs c_5 ---")
trS5_to_c5 = defaultdict(set)
for d in data:
    trS5_to_c5[d['trS5']].add(d['c5'])
for trS5 in sorted(trS5_to_c5.keys()):
    c5s = trS5_to_c5[trS5]
    print(f"  tr(S^5) = {trS5} -> c_5 in {sorted(c5s)}")

# tr(A^5) counts closed walks of length 5.
# Decompose: tr(A^5) = 5*c_5 + (walks through repeated vertices)
# For regular tournament: tr(A^3) = 3*c_3 (no repeated vertices possible in 3-walk)
# Check: tr(A^3) = 3*c_3?
print(f"\n--- Check: tr(A^3) = 3*c_3? ---")
for d in data[:3]:
    print(f"  tr(A^3) = {d['tr3']}, 3*c_3 = {3*d['c3']}, match: {d['tr3'] == 3*d['c3']}")

# For k=5: tr(A^5) = 5*c_5 + correction terms from non-simple walks
# Non-simple 5-walks that close: include walks with repeated vertices
# The correction depends on local structure (number of common neighbors, etc.)
print(f"\n--- Correction: tr(A^5) - 5*c_5 ---")
corrections = set()
for d in data:
    corr = d['tr5'] - 5 * d['c5']
    corrections.add((d['c5'], corr))
for c5, corr in sorted(corrections):
    print(f"  c_5 = {c5}: tr(A^5) - 5*c_5 = {corr}")

# Is the correction constant for regular tournaments?
corr_vals = set(d['tr5'] - 5*d['c5'] for d in data)
print(f"\n  All correction values: {sorted(corr_vals)}")
if len(corr_vals) == 1:
    print(f"  CONSTANT correction = {corr_vals.pop()}")
    print(f"  => c_5 = (tr(A^5) - correction) / 5")
else:
    print(f"  NOT constant -- correction varies with tournament structure")

# Deeper: the correction for 5-walks comes from walks that visit some vertex twice.
# For regular tournament with degree d=(n-1)/2:
# Walk i->j->k->l->m->i with some repeated vertex.
# The simplest: i->j->i->k->l->i (3 visits to i)
# Or: i->j->k->j->l->i (2 visits to j)

# A cleaner formula might use tr(S^5) where S is the skew-adjacency matrix
# S[i][j] = 1 if i->j, S[i][j] = -1 if j->i
# For tournaments: S = 2A - J + I (where J is all-ones, I identity)
# No wait: S[i][j] = T[i][j] - T[j][i] for i != j, S[i][i] = 0
# So S = A - A^T (the antisymmetric part)

# Savchenko uses the "signed permanental" formula which involves S
print(f"\n--- Formula attempt: c_5 = (tr(S^5) + correction) / 10 ---")
for d in data[:3]:
    print(f"  tr(S^5) = {d['trS5']}, c_5 = {d['c5']}, tr(S^5)/10 = {d['trS5']/10}")

# Check if tr(S^5) = 10*c_5 - (something involving c_3)
for_check = set()
for d in data:
    # tr(S^5) should be 10*c_5 - 10*n*c_3 + ... (from expansion)
    # Actually: for the skew matrix, tr(S^k) = (-1)^{(k-1)/2} * 2 * k * c_k
    # for odd k, if all cycles are simple. But this is only for k-cycles.
    # The correction involves walks through repeated vertices.
    # For k=3: tr(S^3) = 6*c_3 (since each directed 3-cycle contributes +1 from
    # S[a][b]*S[b][c]*S[c][a] = 1*1*1 = 1, and the reverse contributes
    # S[a][c]*S[c][b]*S[b][a] = (-1)(-1)(-1) = -1. Net per triple: 1-(-1)=2? No...
    # Actually S[a][b] = 1 if a->b (tournament edge), -1 if b->a.
    # For directed cycle a->b->c->a: S[a][b]*S[b][c]*S[c][a] = 1*1*1 = 1
    # For reverse a->c->b->a: S[a][c]*S[c][b]*S[b][a] = (-1)(-1)(-1) = -1
    # Hmm, wait: if a->b, b->c, c->a then a->c is FALSE (c->a), so
    # S[a][c] = -1. S[c][b] = -1 (b->c so S[c][b] = -1). S[b][a] = -1 (a->b).
    # Product = (-1)^3 = -1.
    pass

# Direct test
print(f"\n--- Direct: tr(S^3) vs c_3 ---")
for d in data[:1]:
    # Compute tr(S^3) for this tournament
    T = tournament_from_bits(n, d['bits'])
    S = [[T[i][j] - T[j][i] for j in range(n)] for i in range(n)]
    Sk = [[int(i == j) for j in range(n)] for i in range(n)]
    for _ in range(3):
        new = [[0]*n for _ in range(n)]
        for i in range(n):
            for j in range(n):
                for l in range(n):
                    new[i][j] += Sk[i][l] * S[l][j]
        Sk = new
    trS3 = sum(Sk[i][i] for i in range(n))
    print(f"  tr(S^3) = {trS3}, c_3 = {d['c3']}")
    print(f"  Ratio: {trS3}/{d['c3']} = {trS3/d['c3'] if d['c3'] != 0 else 'N/A'}")

# Compute tr(S^3) for all regular
trS3_vals = set()
for d in data:
    T = tournament_from_bits(n, d['bits'])
    S = [[T[i][j] - T[j][i] for j in range(n)] for i in range(n)]
    S3 = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            for k in range(n):
                for l in range(n):
                    pass  # too slow, let me just check tr(S^3)
    break

# Actually: tr(S^3) = sum_{i,j,k} S[i][j]*S[j][k]*S[k][i]
# For the directed cycle a->b->c->a: S[a][b]*S[b][c]*S[c][a] = 1*1*1 = 1
# For the NON-cycle triple (transitive) a->b, b->c, a->c:
#   S[a][b]=1, S[b][c]=1, S[c][a]=1 (c->a so S[c][a]=1 if c->a, but a->c => S[c][a]=-1)
#   Wait: a->c means T[a][c]=1, T[c][a]=0. So S[c][a] = T[c][a]-T[a][c] = 0-1 = -1.
#   Product = 1*1*(-1) = -1.
# So for ANY triple {a,b,c}:
#   If cyclic (directed 3-cycle exists): contributes +1 to tr(S^3)
#   If transitive: contributes -1 to tr(S^3)
# Wait but we sum over ALL ordered triples (i,j,k), not just vertex SETS.
# Each 3-element set appears in 3! = 6 orderings.
# For a cyclic triple: 3 orderings give +1 (the 3 rotations of the cycle),
#   3 give -1 (the 3 rotations of the reverse).
# Hmm, let me think more carefully...

# For directed triple a->b->c->a:
#   (a,b,c): S[a][b]*S[b][c]*S[c][a] = 1*1*1 = 1 (this IS the cycle)
#   (b,c,a): S[b][c]*S[c][a]*S[a][b] = 1*1*1 = 1 (same cycle rotated)
#   (c,a,b): S[c][a]*S[a][b]*S[b][c] = 1*1*1 = 1 (same cycle rotated)
#   (a,c,b): S[a][c]*S[c][b]*S[b][a] = (-1)*(-1)*(-1) = -1 (reverse direction)
#   (c,b,a): S[c][b]*S[b][a]*S[a][c] = (-1)*(-1)*(-1) = -1
#   (b,a,c): S[b][a]*S[a][c]*S[c][b] = (-1)*(-1)*(-1) = -1
# Total: 3*(+1) + 3*(-1) = 0 per cyclic triple?!

# Hmm wait, that can't be right. Let me reconsider:
# Actually S[i][i] = 0, and we sum over DISTINCT i,j,k (including i=j etc).
# Oh wait, the diagonal is 0 so terms with repeated indices are 0.
# So tr(S^3) = sum over DISTINCT ordered triples... but also including
# permutations where i=k, etc. No: tr(S^3) = sum_i (S^3)[i][i]
# = sum_i sum_j sum_k S[i][j]*S[j][k]*S[k][i]
# This includes i=j, j=k, etc., but S[i][i]=0 kills those.
# So it's sum over ordered triples with i,j,k all distinct.
#
# For cyclic triple a->b->c->a:
# Orderings (i,j,k) where all distinct and from {a,b,c}: 6 orderings.
# I computed above: total = 3 - 3 = 0.
# That means tr(S^3) = 0 for ALL tournaments? That's wrong.

# Actually let me reconsider. For (a,c,b):
# S[a][c]*S[c][b]*S[b][a]
# a->c is FALSE (c->a in the cycle). So S[a][c] = T[a][c]-T[c][a].
# In cycle a->b->c->a: a->b, b->c, c->a. So T[a][c]=0, T[c][a]=1.
# S[a][c] = 0-1 = -1. OK that's what I had.
# S[c][b]: T[c][b]=0 (b->c), T[b][c]=1. S[c][b] = 0-1 = -1. OK.
# S[b][a]: T[b][a]=0 (a->b), T[a][b]=1. S[b][a] = 0-1 = -1. OK.
# Product = (-1)^3 = -1. OK.

# For transitive a->b, b->c, a->c:
# (a,b,c): S[a][b]*S[b][c]*S[c][a] = 1*1*(-1) = -1
# (b,c,a): S[b][c]*S[c][a]*S[a][b] = 1*(-1)*1 = -1
# (c,a,b): S[c][a]*S[a][b]*S[b][c] = (-1)*1*1 = -1
# (a,c,b): S[a][c]*S[c][b]*S[b][a] = 1*(-1)*(-1) = 1
# (c,b,a): S[c][b]*S[b][a]*S[a][c] = (-1)*(-1)*1 = 1
# (b,a,c): S[b][a]*S[a][c]*S[c][b] = (-1)*1*(-1) = 1
# Total: 3*(-1) + 3*(+1) = 0.

# So EVERY triple contributes 0 to tr(S^3)?!
# That means tr(S^3) = 0 for ALL tournaments.

# Let me verify:
T_check = tournament_from_bits(7, 0)  # first tournament
S_check = [[T_check[i][j] - T_check[j][i] for j in range(7)] for i in range(7)]
# Compute S^3
S2 = [[sum(S_check[i][k]*S_check[k][j] for k in range(7)) for j in range(7)] for i in range(7)]
S3c = [[sum(S2[i][k]*S_check[k][j] for k in range(7)) for j in range(7)] for i in range(7)]
trS3_check = sum(S3c[i][i] for i in range(7))
print(f"\n  Verification: tr(S^3) for bits=0: {trS3_check}")

# So tr(S^k) for odd k: the k=3 case is always 0.
# For k=5, the above data shows tr(S^5) varies: {-280, -360, -420}.
# So tr(S^5) DOES distinguish the 3 regular n=7 classes.

# Relation: is c_5 = -tr(S^5) / 10?
print(f"\n--- Testing: c_5 = -tr(S^5) / 10 ---")
for c5val in sorted(c5_groups.keys()):
    trS5s = set(d['trS5'] for d in c5_groups[c5val])
    for s5 in trS5s:
        ratio = -s5 / (10 * c5val) if c5val else 'N/A'
        print(f"  c_5 = {c5val}, tr(S^5) = {s5}, -trS5/(10*c5) = {ratio}")


# General formula: c_k = (-1)^{(k+1)/2} * tr(S^k) / (2k) for simple cycles?
# Actually the contribution of a directed k-cycle to tr(S^k) is:
# k rotations * (+1)^k = k (since all edges follow cycle direction)
# plus k rotations in reverse * (-1)^k = k*(-1)^k
# Total: k*(1 + (-1)^k) = 0 for odd k, 2k for even k.
# But that gives 0 for all odd k, which contradicts tr(S^5) != 0.

# The issue: tr(S^5) includes non-simple walks too!
# A non-simple walk of length 5 can visit some vertex twice.

# So: tr(S^5) = (simple 5-cycle contribution) + (non-simple walk contribution)
# = 0 + (non-simple walks)
# The non-simple walks of length 5 involve visiting 3 or 4 distinct vertices.

# For 3 distinct vertices: walk i->j->k->j->i->j type... but S has 0 diagonal
# so i->i is impossible. Walk must have all consecutive distinct.
# Length 5, all consecutive distinct, but visiting 3 vertices:
# e.g., i->j->k->i->j->k->i... but length 5 closing back to start.
# Like i->j->i->k->i->j closing? No, length 5 means 5 edges, returning to start.

# Possible with 3 vertices: must visit each vertex multiple times.
# i->j->k->i->j->i: this has length 5 and returns to i.
# But intermediate: i,j,k,i,j,i uses i 3 times. Non-simple.
# S[i][j]*S[j][k]*S[k][i]*S[i][j]*S[j][i] = 1*1*1*1*(-1) = -1
# (assuming cycle i->j->k->i)

# This is getting complex. Let me just verify the formula numerically.
print(f"\n--- Numerical formula search ---")
print("Trying: c_5 = a * tr(A^5) + b * n * c_3 + c * n")
print("We have 3 data points (one per class):")

classes = []
for c5val in sorted(c5_groups.keys()):
    tr5 = list(set(d['tr5'] for d in c5_groups[c5val]))[0]
    c3 = list(set(d['c3'] for d in c5_groups[c5val]))[0]
    classes.append((c5val, tr5, c3))
    print(f"  c_5={c5val}, tr(A^5)={tr5}, c_3={c3}")

# c_5 = a * tr(A^5) + b * n * c_3 + c * n
# 28 = a * 308 + b * 98 + c * 7
# 36 = a * 348 + b * 98 + c * 7
# 42 = a * 378 + b * 98 + c * 7

# From equations 1 and 2: 8 = 40*a => a = 1/5
# From equations 2 and 3: 6 = 30*a => a = 1/5. Consistent!
# From equation 1: 28 = 308/5 + b*98 + 7c = 61.6 + 98b + 7c
# 28 - 61.6 = 98b + 7c => -33.6 = 98b + 7c

# Hmm, non-integer. Let me try: c_5 = (tr(A^5) - correction) / 5
# 28 = (308 - corr) / 5 => corr = 308 - 140 = 168
# 36 = (348 - corr) / 5 => corr = 348 - 180 = 168
# 42 = (378 - corr) / 5 => corr = 378 - 210 = 168
# CONSTANT correction = 168!

print(f"\n  DISCOVERY: c_5 = (tr(A^5) - 168) / 5 for n=7 regular")
print(f"  Checking: 168 = ?")
print(f"  168 = 24 * 7 = n * 24 = n * (n-1)*(n-2) * ... nah")
print(f"  168 = 8 * 21 = 8 * C(7,2)")
print(f"  168 = n * (n-1) * (n-2) / (n-3)? = 7*6*5/4 = 52.5 no")

# The correction 168 for n=7:
# Non-simple 5-walks that close: walks of length 5 returning to start,
# visiting at most 4 distinct vertices.
# For regular n=7 with degree 3:
# Number of such walks = 168 (constant for all regular tournaments)

# Check: is correction = n * d * (d-1) * (something)?
# d = 3 for n=7. 7 * 3 * 2 = 42. 168/42 = 4. So 168 = 4*42 = 4*n*d*(d-1)/... no.
# 168 = 7 * 24. d*(d+1)*(2d+1)/... 3*4*7/... = 84/... hmm.
# 168 = 7 * 3 * 8 = n * d * 2^3
# 168 = C(8,3) = 56 * 3 = ... no, C(8,3) = 56.
# 7*24 = 7 * 4! = 7! / (7-4)! / ... no.
# Actually 168 = 7! / (7-4)! / ... = 7*6*5*4 / ... 840/5 = 168!
# 168 = P(7,3) - something? P(7,3) = 210. No.
# 168 = 7 * 4! = 168. Yes!
# Or: 168 = n * (n-3)!  = 7 * 24 = 168. YES!

n_val = 7
d_val = (n_val - 1) // 2
from math import factorial
print(f"  168 = n * (n-3)! = {n_val} * {factorial(n_val-3)} = {n_val * factorial(n_val-3)}")
print(f"  Alternative: 168 = n * 4! = 7 * 24")

# So the formula might be:
# c_5 = (tr(A^5) - n * (n-3)!) / 5 for regular tournaments
# This would need verification at other n.

# Let's check at n=5
print(f"\n--- Check at n=5 ---")
n5 = 5
m5 = n5 * (n5 - 1) // 2
for bits in range(1 << m5):
    T5 = tournament_from_bits(n5, bits)
    if not all(sum(T5[i]) == 2 for i in range(n5)):
        continue
    c5_val = count_directed_k_cycles_dp(T5, 5)
    tr5_val = matrix_power_trace(T5, 5)
    predicted_corr = n5 * factorial(n5 - 3)
    predicted_c5 = (tr5_val - predicted_corr) / 5
    print(f"  bits={bits}: c_5={c5_val}, tr(A^5)={tr5_val}, predicted c_5={(tr5_val - predicted_corr)/5}")
    break  # Just check one

# More careful: for n=5, d=2.
# n*(n-3)! = 5*2! = 10
# c_5 = (tr(A^5) - 10) / 5
# For the cyclic tournament T_5: c_5 = 2 (two directed 5-cycles: clockwise and counter... wait)
# Actually in T_5 = 0->1->2->3->4->0 (and 0->2, 1->3, 2->4, 3->0, 4->1):
# Directed 5-cycle: 0->1->2->3->4->0 is one. Others?
# 0->2->4->1->3->0: check edges. 0->2 yes, 2->4 yes, 4->1 yes, 1->3 yes, 3->0 yes. YES!
# So there are at least 2. How many total directed 5-cycles on 5 vertices?
# For cyclic T_5: the group Z_5 acts, so cycles come in orbits of size 5/gcd.
# A directed 5-cycle visits all 5 vertices. Total = (5-1)! / 5 * (valid fraction) = ...
# Actually (5-1)! = 24 permutations of the cycle. Each directed cycle = 5 rotations.
# So at most 24/5 = 4.8 -> at most 4 directed cycles (excluding rotations).
# We found 2. The other two would be 0->1->3->4->2->0 (check: 1->3 yes, 3->4 yes, 4->2 yes, 2->0 yes) -> YES!
# And 0->2->1->4->3->0 (check: 0->2 yes, 2->1? T[2][1]: 2 beats 1? In T_5, out-neighbors of 2 are 3,4. So T[2][1]=0. NO.)
# So 3 directed 5-cycles? Let me just compute.

for bits in range(1 << m5):
    T5 = tournament_from_bits(n5, bits)
    if not all(sum(T5[i]) == 2 for i in range(n5)):
        continue
    c5_val = count_directed_k_cycles_dp(T5, 5)
    tr5_val = matrix_power_trace(T5, 5)
    print(f"  Regular n=5: bits={bits}, c_5={c5_val}, tr(A^5)={tr5_val}")


print("\nDone.")
