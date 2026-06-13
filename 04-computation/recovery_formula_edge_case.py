#!/usr/bin/env python3
"""
Check the edge case: when d_i = d_j, can the formula A_{ij}=1 iff d_i-d_j=n-1-σ-2λ
still work?

When d_i = d_j, the formula becomes: A_{ij}=1 iff 0 = n-1-σ-2λ, i.e., σ+2λ = n-1.

But for A_{ij}=0 (A_{ji}=1): the formula also gives 0 = d_j-d_i = n-1-σ-2λ.

So when d_i = d_j AND σ+2λ = n-1: both A_{ij}=1 and A_{ij}=0 satisfy the formula!
And when d_i = d_j AND σ+2λ ≠ n-1: neither direction satisfies the formula!

Wait, that can't be right since the formula was verified computationally.
Let me re-check.

opus-2026-03-13-S71c
"""
import numpy as np
from collections import defaultdict

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

n = 7
tb = n*(n-1)//2
np.random.seed(42)

equal_degree_count = 0
formula_ok = 0
formula_bad = 0

for trial in range(10000):
    bits = np.random.randint(0, 1 << tb)
    A = bits_to_adj(bits, n)
    scores = [int(sum(A[i])) for i in range(n)]

    L = np.zeros((n, n), dtype=int)
    S = np.zeros((n, n), dtype=int)
    for u in range(n):
        for v in range(u+1, n):
            for w in range(n):
                if w == u or w == v: continue
                if (A[u][v] and A[v][w] and A[w][u]) or (A[v][u] and A[u][w] and A[w][v]):
                    L[u][v] += 1; L[v][u] += 1
                if A[u][w] and A[v][w]: S[u][v] += 1; S[v][u] += 1
                if A[w][u] and A[w][v]: S[u][v] += 1; S[v][u] += 1

    for i in range(n):
        for j in range(i+1, n):
            if scores[i] == scores[j]:
                equal_degree_count += 1
                lam = L[i][j]
                sig = S[i][j]
                rhs = n - 1 - sig - 2*lam
                lhs = scores[i] - scores[j]  # = 0

                aij = A[i][j]
                if aij == 1 and lhs == rhs:
                    formula_ok += 1
                elif aij == 0 and lhs != rhs:
                    formula_ok += 1
                elif aij == 1 and lhs != rhs:
                    formula_bad += 1
                    if formula_bad <= 3:
                        print(f"  FAIL A=1: ({i},{j}), d=({scores[i]},{scores[j]}), λ={lam}, σ={sig}, rhs={rhs}")
                elif aij == 0 and lhs == rhs:
                    formula_bad += 1
                    if formula_bad <= 3:
                        print(f"  FAIL A=0: ({i},{j}), d=({scores[i]},{scores[j]}), λ={lam}, σ={sig}, rhs={rhs}")

print(f"Equal degree pairs checked: {equal_degree_count}")
print(f"Formula OK: {formula_ok}")
print(f"Formula BAD: {formula_bad}")

# Now let's understand: when d_i = d_j, what determines the direction?
# From the algebra:
# If A_{ij}=1: λ = d_j - co, σ = co + ci = co + (n-2-d_i+A_{ij}-d_j+A_{ji}+2*co)/...
# Hmm, let me just check what sigma+2*lambda equals when d_i=d_j.

# If A_{ij}=1 and d_i=d_j: σ+2λ = σ + 2(d_j-co) = σ + 2d_j - 2co = (co+ci) + 2d_j - 2co
#   = ci - co + 2d_j
# And co - ci = d_i + d_j - n + 1 = 2d_i - n + 1 (since d_i=d_j)
# So ci = co - (2d_i - n + 1)
# σ+2λ = (co-(2d_i-n+1)) - co + 2d_i = -(2d_i-n+1) + 2d_i = n-1
# So σ+2λ = n-1 ALWAYS when d_i=d_j and A_{ij}=1!

# Similarly, if A_{ij}=0 (A_{ji}=1) and d_i=d_j:
# λ = P_{ij} = d_i - 0 - co = d_i - co
# σ+2λ = (co+ci) + 2(d_i-co) = ci - co + 2d_i = -(2d_i-n+1) + 2d_i = n-1
# SAME! So σ+2λ = n-1 ALWAYS when d_i=d_j, regardless of A_{ij}!

# This means: when d_i=d_j, rhs=0=lhs ALWAYS, so the formula says "A_{ij}=1" for BOTH directions!
# But we verified 0 failures... so what's going on?

# Wait — my verification code: "aij == 1 and lhs == rhs" → formula_ok.
# When d_i=d_j: lhs=0 and rhs=0 always. So if A_{ij}=1, formula says "yes" → OK.
# If A_{ij}=0, lhs=0 and rhs=0, formula says "A_{ij}=1" but actually A_{ij}=0 → FAIL!

# But I checked "aij == 0 and lhs != rhs" → formula_ok. When lhs=rhs=0, this is FALSE,
# so we DON'T count it as OK. And "aij == 0 and lhs == rhs" → formula_bad. YES!
# So there SHOULD be failures!

# Let me re-examine: my code checks the IFF both ways.
# For i<j only. A_{ij}=0 means j→i.
# formula_bad increments when (aij==0 and lhs==rhs) — this SHOULD happen when d_i=d_j.

print(f"\nExpected: all d_i=d_j pairs with A_ij=0 should be formula_bad.")
print(f"But formula_bad = {formula_bad}")

# Hmm, maybe I miscounted. Let me recheck with explicit tracking.
count_aij1_deq = 0
count_aij0_deq = 0
for trial in range(100):
    bits = np.random.randint(0, 1 << tb)
    A = bits_to_adj(bits, n)
    scores = [int(sum(A[i])) for i in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if scores[i] == scores[j]:
                if A[i][j] == 1: count_aij1_deq += 1
                else: count_aij0_deq += 1

print(f"\n100 samples: d_i=d_j pairs: A_ij=1: {count_aij1_deq}, A_ij=0: {count_aij0_deq}")
print(f"(About half should be each direction by symmetry)")

print(f"\nDone.")
