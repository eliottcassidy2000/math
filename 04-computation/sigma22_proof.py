#!/usr/bin/env python3
"""
ALGEBRAIC PROOF: sigma((2,2)) = 16*t3 + 8*bc + 140 at n=7.

DERIVATION:
  sigma({0,1,3,4}) = sum_P A[p_0,p_1]*A[p_1,p_2]*A[p_3,p_4]*A[p_4,p_5]

  This involves 6 vertex positions (0-5), with one free vertex p_6.
  The two groups {p_0,p_1,p_2} and {p_3,p_4,p_5} form two directed 2-paths.

  For a 3-vertex set S, define dp_2(S) = #{directed 2-paths in T[S]}:
    dp_2(S) = 3 if S is a cyclic triple
    dp_2(S) = 1 if S is a transitive triple

  sigma = 2 * sum_v sum_{unord part (S,S') of V\{v}} dp_2(S) * dp_2(S')
        (factor 2 from ordered vs unordered partitions)

  For each v:
    sum_{partitions} dp_2(S)*dp_2(S') = 9*bc(v) + 3*mixed(v) + 1*trans_pair(v)
    where bc(v) + mixed(v) + trans_pair(v) = C(6,3)/2 = 10

  Using mixed(v) = c_3(T-v) - 2*bc(v):
    = 10 + 4*bc(v) + 2*c_3(T-v)

  Summing over v:
    sum_v c_3(T-v) = (n-3)*t3 = 4*t3
    sum_v bc(v) = bc (total complementary-cyclic count)

  sigma = 2*(7*10 + 4*bc + 2*4*t3) = 140 + 8*bc + 16*t3

opus-2026-03-06-S28
"""
from itertools import permutations, combinations
from math import factorial
import random

def count_3_cycles(A, n):
    count = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if A[i][j]*A[j][k]*A[k][i]: count += 1
                if A[i][k]*A[k][j]*A[j][i]: count += 1
    return count

def count_3_cycles_sub(A, verts):
    count = 0
    for i in range(len(verts)):
        for j in range(i+1, len(verts)):
            for k in range(j+1, len(verts)):
                a, b, c = verts[i], verts[j], verts[k]
                if A[a][b]*A[b][c]*A[c][a]: count += 1
                if A[a][c]*A[c][b]*A[b][a]: count += 1
    return count

def is_3cycle(A, triple):
    a, b, c = triple
    return (A[a][b]*A[b][c]*A[c][a] + A[a][c]*A[c][b]*A[b][a]) > 0

def both_cyclic_total(A, n):
    total = 0
    for v in range(n):
        S = [u for u in range(n) if u != v]
        for T in combinations(S, 3):
            T_comp = tuple(u for u in S if u not in T)
            if T < T_comp:
                if is_3cycle(A, T) and is_3cycle(A, T_comp):
                    total += 1
    return total

def compute_sigma(A, n, positions):
    total = 0
    for p in permutations(range(n)):
        prod = 1
        for i in positions:
            prod *= A[p[i]][p[i+1]]
        total += prod
    return total

# =====================================================================
n = 7
print("=" * 70)
print(f"PROOF VERIFICATION: sigma((2,2)) = 16*t3 + 8*bc + 140")
print("=" * 70)

all_pass = True
for trial in range(15):
    random.seed(trial * 37 + 11)
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1

    t3 = count_3_cycles(A, n)
    bc = both_cyclic_total(A, n)

    # All (2,2) position subsets
    pos_22 = [(0,1,3,4), (0,1,4,5), (1,2,4,5)]
    for pos in pos_22:
        actual = compute_sigma(A, n, pos)
        predicted = 16*t3 + 8*bc + 140
        match = actual == predicted
        if not match:
            all_pass = False
        if trial < 5 or not match:
            print(f"  T{trial:2d}: t3={t3:2d} bc={bc:2d} pos={pos} "
                  f"actual={actual} pred={predicted} {'OK' if match else 'FAIL'}")

if all_pass:
    print(f"\nALL {15*3} TESTS PASS!")

# =====================================================================
# Now verify sigma((3,)): 3 consecutive positions = directed 3-path
# =====================================================================
print(f"\n{'='*70}")
print("sigma((3,)) = ?")
print(f"{'='*70}")

# sigma({0,1,2}) = sum_P A[p_0,p_1]*A[p_1,p_2]*A[p_2,p_3]
# = (n-4)! * sum_{ordered 4-tuples} A[a,b]*A[b,c]*A[c,d]
# = (n-4)! * #{directed 3-paths on n vertices}

# How does #{directed 3-paths} relate to t3?
# A directed 3-path a->b->c->d involves 4 vertices.
# For each 4-vertex subtournament, count directed 3-paths.

# For a 4-vertex tournament with score sequence (s0,s1,s2,s3):
# Each 4-vertex set has C(4,2)=6 edges. The tournament can have
# 0, 1, or 4 directed 3-cycles. The number of directed 3-paths
# depends on the score sequence.

# Actually, for 4 vertices: #directed 3-paths = #Ham paths in the 4-vertex sub.
# And at n=4, H(T) = c3 + 1 for any 4-vertex tournament (by OCF: H = 1 + 2*alpha_1,
# alpha_1 = #odd cycles = c3, and alpha_2=0 since two disjoint 3-cycles need 6 > 4).

# Wait, H(T) counts UNDIRECTED Ham paths, i.e., both directions. No, H(T) counts
# DIRECTED Hamiltonian paths. And a directed 3-path on 4 vertices IS a directed Ham path.
# So #{directed 3-paths on 4 vertices} = H(T[4-vertex sub]).

# And by OCF at n=4 (even): wait, OCF only works at odd n? No, OCF works at all n.
# H(T) = I(Omega(T), 2) for ALL tournaments.

# At n=4: Omega has only 3-cycles (no 5-cycles possible). Disjoint pair impossible.
# So alpha_1 = c3(4-vertex sub), alpha_2 = 0.
# H = 1 + 2*c3.

# But wait, for a 4-vertex tournament:
# Transitive: c3=0, H=1. Direct check: only 1 Ham path (along score order).
# One 3-cycle: c3=1, H=3.
# Four 3-cycles: c3=4, H=9? No, max c3 at n=4 is 4 (all C(4,3)=4 triples cyclic)?
# No — a 4-vertex regular tournament has score (1.5,...), not integer. So regular is impossible.
# Actually at n=4, possible score sequences are (0,1,2,3), (1,1,1,3), (0,2,2,2), (1,1,2,2).
# c3 values: (0,1,2,3)→0, (1,1,1,3)→1, (0,2,2,2)→1, (1,1,2,2)→2.
# So H = 1 + 2*c3 = 1,3,3,5.
# Total # directed Ham paths at n=4: max is 5 (two 3-cycles).

# Now, sum over 4-vertex subsets of H(T[S]):
# sum_S H(T[S]) = sum_S (1 + 2*c3(S)) = C(n,4) + 2*sum_S c3(S)
# And sum_S c3(S) = C(n-3, 4-3) * t3 = (n-3)*t3 = 4*t3.
# So sum_S H(T[S]) = C(7,4) + 2*(n-3)*t3 = 35 + 8*t3.

# sigma((3,)) = (n-4)! * sum of ordered 4-tuples of directed 3-paths
# A directed 3-path on 4 vertices: there are H(T[S]) for each 4-set S,
# but each path uses a SPECIFIC ordering of the 4 vertices.
# Actually, sum over ordered 4-tuples (a,b,c,d) all distinct of A[a,b]*A[b,c]*A[c,d]
# = sum over 4-element subsets S of #{directed Ham paths in T[S]}.

# Wait, that's the count of ORDERED directed Ham paths = H(T[S]) since a directed Ham path
# is an ordered sequence.

# So sigma((3,)) = (n-4)! * sum_S H(T[S]) = 6 * (35 + 8*t3) = 210 + 48*t3

for trial in range(10):
    random.seed(trial * 37 + 11)
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1

    t3 = count_3_cycles(A, n)
    actual = compute_sigma(A, n, (0, 1, 2))
    predicted = 210 + 48*t3
    print(f"  T{trial}: t3={t3:2d}, sigma(0,1,2)={actual}, pred={predicted}, "
          f"{'OK' if actual == predicted else 'FAIL'}")

# =====================================================================
# sigma((4,)): 4 consecutive positions = directed 4-path on 5 vertices
# =====================================================================
print(f"\n{'='*70}")
print("sigma((4,)) = ?")
print(f"{'='*70}")

# sigma((4,)) = (n-5)! * sum_S H(T[S]) for 5-vertex subsets S
# At n=7: (n-5)! = 2! = 2
# H(T[5-vertex sub]) = 1 + 2*(c3(S) + c5(S)) by OCF at n=5

# sum_S H(T[S]) = sum_S [1 + 2*c3(S) + 2*c5(S)]
# = C(7,5) + 2*C(4,2)*t3 + 2*sum_S c5(S)
# Wait, sum_S c3(S) where S is a 5-vertex subset: each 3-cycle is in C(n-3,5-3)=C(4,2)=6 subsets.
# So sum_S c3(S) = 6*t3.
# sum_S c5(S) = C(n-5,5-5)*t5 = 1*t5 = t5. (Each 5-cycle is in exactly 1 5-vertex subset — wait,
# no, it's in C(n-5, 0) = 1 subset at n=7 since 5-vertex subset = exactly the 5 vertices of the cycle.)
# Wait, at n=7, a 5-vertex subset contains C(7-5, 5-5)=C(2,0)=1... no. A 5-cycle uses 5 vertices.
# The only 5-vertex subset containing those 5 vertices IS those 5 vertices. So each 5-cycle
# appears in exactly C(7-5, 5-5) = C(2,0) = 1 subset. So sum_S c5(S) = t5.

# sum_S H(T[S]) = C(7,5) + 2*6*t3 + 2*t5 = 21 + 12*t3 + 2*t5
# sigma((4,)) = 2 * (21 + 12*t3 + 2*t5) = 42 + 24*t3 + 4*t5

from itertools import combinations as comb_it

def count_5_cycles_dp(A, n):
    count = 0
    for verts in comb_it(range(n), 5):
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

for trial in range(10):
    random.seed(trial * 37 + 11)
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1

    t3 = count_3_cycles(A, n)
    t5 = count_5_cycles_dp(A, n)
    actual = compute_sigma(A, n, (0, 1, 2, 3))
    predicted = 42 + 24*t3 + 4*t5
    print(f"  T{trial}: t3={t3:2d} t5={t5:2d}, sigma(0,1,2,3)={actual}, "
          f"pred={predicted}, {'OK' if actual == predicted else 'FAIL'}")

# =====================================================================
# sigma((6,)): all 6 positions = H (Hamiltonian path count)
# =====================================================================
print(f"\n{'='*70}")
print("sigma((6,)) = H")
print(f"{'='*70}")

def ham_count_dp(A, n):
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)) or dp[mask][v] == 0: continue
            for u in range(n):
                if (mask & (1 << u)) or not A[v][u]: continue
                dp[mask | (1 << u)][u] += dp[mask][v]
    return sum(dp[(1<<n)-1][v] for v in range(n))

for trial in range(5):
    random.seed(trial * 37 + 11)
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1

    H = ham_count_dp(A, n)
    sigma_6 = compute_sigma(A, n, (0, 1, 2, 3, 4, 5))
    print(f"  T{trial}: H={H}, sigma(0,...,5)={sigma_6}, match={H==sigma_6}")

# =====================================================================
# COMPLETE TABLE: all sigma formulas at n=7
# =====================================================================
print(f"\n{'='*70}")
print("COMPLETE SIGMA TABLE AT n=7")
print(f"{'='*70}")
print(f"""
  Pattern   | Formula                        | Depends on
  ----------|--------------------------------|----------
  (1,)      | n!/2 = 2520                    | universal
  (1,1)     | (n-4)! * n(n-1)(n-2)(n-3)/4    | universal
            | = 6 * 210 = 1260               |
  (2,)      | (n-3)! * (C(n,3) + 2*t3)      | t3
  (1,1,1)   | 630                            | universal
  (2,1)     | 24*t3 + 420                    | t3
  (3,)      | 48*t3 + 210                    | t3
  (2,1,1)   | 12*t3 + ... (TBD)             | t3
  (2,2)     | 16*t3 + 8*bc + 140            | t3, bc
  (3,1)     | 24*t3 + 2*t5 + ... (TBD)      | t3, t5
  (4,)      | 24*t3 + 4*t5 + 42             | t3, t5
  (3,2)     | TBD                            | t3, t5, bc?
  (4,1)     | TBD                            | t3, t5
  (5,)      | TBD                            | t3, t5, bc?
  (6,)      | H                              | H
""")

print(f"{'='*70}")
print("DONE")
print(f"{'='*70}")
