#!/usr/bin/env python3
"""
pfaffian_structure.py - Deep analysis of the Pfaffian-Betti connection.

WHY does |Pf(S)| determine the "topological phase" of a tournament?

Exhaustive n=6 data:
  β₁>0: |Pf| ∈ {1,3} only (4800 tournaments)
  β₃>0: |Pf| ∈ {7,9} only (320 tournaments)
  β=0: |Pf| ∈ {1,3,5,7,9} (27648 tournaments)

The Pfaffian Pf(S) counts SIGNED PERFECT MATCHINGS of K_6 weighted by tournament arcs.
There are 15 perfect matchings of K_6.

KEY ALGEBRAIC IDENTITY: For a tournament T with skew-adjacency S,
  Pf(S) = Σ_{M ∈ PM(K_n)} ε(M) · ∏_{(i,j)∈M} T[i,j]^{sign}

where T[i,j]^{sign} = +1 if i→j, -1 if j→i.

QUESTION 1: What is the combinatorial meaning of |Pf(S)|?
QUESTION 2: Is there a formula |Pf(S)|² = f(t₃, score, ...)?
QUESTION 3: Does |Pf(S)| relate to the independence polynomial I(Ω,x)?

Author: opus-2026-03-07-S46f
"""
import sys
sys.path.insert(0, '/Users/e/Documents/GitHub/math/04-computation')
from path_homology_v2 import path_betti_numbers
from itertools import combinations, permutations
from collections import Counter, defaultdict
import numpy as np
from math import factorial

def count_t3(A, n):
    return sum(1 for i,j,k in combinations(range(n), 3)
               if (A[i][j] and A[j][k] and A[k][i]) or
                  (A[i][k] and A[k][j] and A[j][i]))

def score_seq(A, n):
    return tuple(sorted(sum(A[i][j] for j in range(n)) for i in range(n)))

def skew_adj(A, n):
    S = np.zeros((n, n), dtype=int)
    for i in range(n):
        for j in range(n):
            if i != j:
                S[i][j] = 1 if A[i][j] else -1
    return S

def pfaffian_6(S):
    """Compute Pfaffian of 6x6 skew-symmetric matrix via expansion."""
    # Pf(S) for 6x6 = sum over 15 perfect matchings of K_6
    n = 6
    # Perfect matchings of {0,1,2,3,4,5}:
    # Enumerate by first choosing partner of 0
    pf = 0
    remaining = list(range(1, n))
    for partner_idx, p0 in enumerate(remaining):
        rest = [x for x in remaining if x != p0]
        sign_0 = (-1) ** partner_idx
        # Now match rest = [a,b,c,d] (4 elements)
        # 3 matchings: (a,b)(c,d), (a,c)(b,d), (a,d)(b,c)
        a, b, c, d = rest
        matchings = [
            ((a,b), (c,d), 1),
            ((a,c), (b,d), -1),
            ((a,d), (b,c), 1),
        ]
        for (p1, q1), (p2, q2), sign_inner in matchings:
            sign = sign_0 * sign_inner
            val = S[0, p0] * S[p1, q1] * S[p2, q2]
            pf += sign * val
    return pf

def count_t5(A, n):
    t = 0
    for combo in combinations(range(n), 5):
        for perm in permutations(combo):
            if all(A[perm[i]][perm[(i+1)%5]] for i in range(5)):
                t += 1
    return t // 5  # Each 5-cycle counted 5 times

def alpha2(A, n):
    """Count vertex-disjoint pairs of directed odd cycles."""
    cycles_3 = []
    for combo in combinations(range(n), 3):
        i, j, k = combo
        if (A[i][j] and A[j][k] and A[k][i]) or (A[i][k] and A[k][j] and A[j][i]):
            cycles_3.append(frozenset(combo))

    count = 0
    for i in range(len(cycles_3)):
        for j in range(i+1, len(cycles_3)):
            if len(cycles_3[i] & cycles_3[j]) == 0:
                count += 1
    return count

# === EXHAUSTIVE n=6 ===
print("=" * 70)
print("n=6: EXHAUSTIVE Pfaffian analysis")
print("=" * 70)

n = 6
m = n*(n-1)//2
total = 1 << m

# Collect (Pf, t3, score, beta_type)
data = []
pf_by_class = defaultdict(list)

for bits in range(total):
    A = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if (bits >> idx) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1

    S = skew_adj(A, n)
    pf = pfaffian_6(S)
    t3 = count_t3(A, n)
    score = score_seq(A, n)

    beta = path_betti_numbers(A, n)
    b1 = int(beta[1]) if len(beta) > 1 else 0
    b3 = int(beta[3]) if len(beta) > 3 else 0

    btype = 'b3' if b3 > 0 else ('b1' if b1 > 0 else 'b0')

    data.append({
        'pf': pf, 't3': t3, 'score': score, 'btype': btype
    })

    key = (t3, score)
    pf_by_class[key].append((pf, btype))

    if (bits + 1) % 10000 == 0:
        print(f"  ... {bits+1}/{total}")

# Analyze Pf² vs (t3, score)
print(f"\nTotal: {len(data)}")

# Is Pf determined by (t3, score)?
print("\n" + "=" * 70)
print("IS Pf DETERMINED BY (t₃, score)?")
print("=" * 70)

for key in sorted(pf_by_class.keys()):
    pf_vals = [x[0] for x in pf_by_class[key]]
    pf_set = sorted(set(pf_vals))
    btypes = Counter(x[1] for x in pf_by_class[key])
    if len(pf_set) > 1 or btypes.get('b3', 0) > 0 or btypes.get('b1', 0) > 0:
        print(f"  t₃={key[0]}, score={key[1]}: Pf ∈ {pf_set}, β: {dict(btypes)}")

# Is |Pf|² a function of t3 alone?
print("\n" + "=" * 70)
print("|Pf|² vs t₃")
print("=" * 70)

pf2_by_t3 = defaultdict(Counter)
for d in data:
    pf2_by_t3[d['t3']][d['pf']**2] += 1

for t3 in sorted(pf2_by_t3.keys()):
    dist = dict(sorted(pf2_by_t3[t3].items()))
    print(f"  t₃={t3}: |Pf|² dist = {dist}")

# Now the KEY: is there a FORMULA for Pf(S) in terms of cycle counts?
print("\n" + "=" * 70)
print("PFAFFIAN FORMULA")
print("=" * 70)

# For n=6: Pf(S) is a degree-3 polynomial in the S[i,j] entries.
# Since S[i,j] = 2*A[i,j] - 1 for i<j, we have S[i,j] = 1 - 2*(1-A[i,j])
# The Pfaffian involves all perfect matchings of K_6.

# Key identity: det(S) = Pf(S)²
# For a tournament: det(S) can be expressed in terms of cycle invariants.

# Let me compute det(S) analytically for tournaments.
# S = 2A - J + I for a tournament (where J is all-ones, I is identity... wait)
# Actually S[i,j] = A[i,j] - A[j,i] = A[i,j] - (1-A[i,j]) = 2A[i,j] - 1 for i≠j
# S = 2A - J + I (where J is all-ones with zero diagonal... no)
# S[i,i] = 0, S[i,j] = 2A[i,j] - 1 for i≠j
# So S = 2A - (J - I) where J is all-1 matrix with 0 diagonal

# For a tournament with adjacency A and complement A^c = J-I-A:
# S = A - A^T = A - (J-I-A) = 2A - J + I
# Actually: A+A^T = J-I, so A^T = J-I-A, and S = A-A^T = A-(J-I-A) = 2A-J+I

# det(S) = det(2A - J + I) where J has J[i,j] = 1 for all i≠j, J[i,i] = 0

# This is hard to simplify directly. Let me check numerically:
# Is Pf related to t₃ and t₅?

# At n=6: |Pf| ∈ {1,3,5,7,9}
# t₃ range: 0 to 10 (max for n=6)
# Hypothesis: |Pf| = ??? function of (t₃, t₅, α₂)

# Let me compute t₅ and α₂ for each tournament and check
print("\nComputing t₅, α₂ for representative tournaments...")

# Sample from each (|Pf|, btype) class
import random
random.seed(42)

for target_pf in [1, 3, 5, 7, 9]:
    for target_btype in ['b0', 'b1', 'b3']:
        matches = [d for d in data if abs(d['pf']) == target_pf and d['btype'] == target_btype]
        if not matches:
            continue
        # Take first match
        d = matches[0]
        # Reconstruct tournament
        bits = data.index(d)
        A = [[0]*n for _ in range(n)]
        idx = 0
        for i in range(n):
            for j in range(i+1, n):
                if (bits >> idx) & 1:
                    A[i][j] = 1
                else:
                    A[j][i] = 1
                idx += 1
        t5 = count_t5(A, n)
        a2 = alpha2(A, n)
        H = sum(1 for perm in permutations(range(n))
                if all(A[perm[i]][perm[i+1]] for i in range(n-1)))
        print(f"  |Pf|={target_pf}, {target_btype}: t₃={d['t3']}, t₅={t5}, α₂={a2}, H={H}, score={d['score']}")

# Check: is there a formula Pf² = f(t₃)?
# From the t₃ data above, check if Pf² is determined by t₃ + score
print("\n" + "=" * 70)
print("FORMULA SEARCH: Pf vs (t₃, score)")
print("=" * 70)

pf_by_t3_score = defaultdict(set)
for d in data:
    pf_by_t3_score[(d['t3'], d['score'])].add(abs(d['pf']))

# Check if (t₃, score) determines |Pf|
determined = True
for key, vals in pf_by_t3_score.items():
    if len(vals) > 1:
        determined = False
        print(f"  MULTIPLE |Pf|: t₃={key[0]}, score={key[1]}: |Pf| ∈ {sorted(vals)}")

if determined:
    print("  YES: (t₃, score) completely determines |Pf| at n=6!")
    # Print the formula
    print("\n  FORMULA (t₃, score) → |Pf|:")
    for key in sorted(pf_by_t3_score.keys()):
        pf_val = list(pf_by_t3_score[key])[0]
        print(f"    t₃={key[0]}, score={key[1]}: |Pf| = {pf_val}")
else:
    print("  NO: (t₃, score) does NOT determine |Pf|")

# Check simpler: is |Pf| determined by score alone?
pf_by_score = defaultdict(set)
for d in data:
    pf_by_score[d['score']].add(abs(d['pf']))

print("\n|Pf| vs score sequence:")
for score in sorted(pf_by_score.keys()):
    vals = sorted(pf_by_score[score])
    print(f"  score={score}: |Pf| ∈ {vals}")

print("\n" + "=" * 70)
print("INTERPRETATION")
print("=" * 70)
print("""
The Pfaffian Pf(S) of the skew-adjacency matrix of a tournament encodes
the tournament's "matching character" — a signed count of how well the
tournament's arcs align with perfect matchings of the complete graph.

CONJECTURE: |Pf(S)| is determined by the score sequence alone at n=6.
If true, this gives a very simple formula: score → |Pf| → constraints on β.

The β₁/β₃ separation by |Pf| then reduces to:
- LOW |Pf| (score close to irregular): favorable for 1-dimensional topology (β₁)
- HIGH |Pf| (score close to regular): favorable for 3-dimensional topology (β₃)

This connects to the known result that regular tournaments maximize H (= I(Ω,2)),
and β₃>0 correlates with high H. The Pfaffian provides the ALGEBRAIC BRIDGE.
""")
