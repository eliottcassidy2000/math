#!/usr/bin/env python3
"""
Investigate Alekseyev-type formulas for NUD perms with unit ascent marking.

Alekseyev's formula for Hertzsprung (no adjacencies at all):
A002464(n) = n! + Σ_{k=1}^{n} (-1)^k Σ_{t=1}^{k} C(k-1,t-1) C(n-k,t) 2^t (n-k)!

For NUD (no unit descents), the cluster method gives:
NUD(n) = Σ_{k=0}^{n-1} (-1)^k * #{ways to force k unit descents} * (remaining perms)

Key insight: in Alekseyev's formula, the "2^t" counts how each cluster can be
ascending or descending. For NUD, clusters are ONLY descending (we forbid descents).
So the NUD analog should replace 2^t with 1^t = 1.

Let me verify: NUD analog of Alekseyev:
NUD(n) = n! + Σ_{k=1}^{n-1} (-1)^k Σ_{t=1}^{k} C(k-1,t-1) C(n-k,t) 1^t (n-k)!
       = n! + Σ_{k=1}^{n-1} (-1)^k Σ_{t=1}^{k} C(k-1,t-1) C(n-k,t) (n-k)!
       = Σ_{k=0}^{n-1} (-1)^k Σ_{t=1}^{k} C(k-1,t-1) C(n-k,t) (n-k)!   [with k=0 giving n!]

Then for W_u(n) = Σ_j u^j N(n,j), we mark unit ASCENTS with u.
When we do inclusion-exclusion on unit descents, the remaining permutation
(with descent clusters merged into single objects) can still have unit ascents.
A unit ascent at position i requires σ(i+1) = σ(i) + 1. If position i is inside
a descent cluster, it can't have a unit ascent (the cluster forces descent).
If position i is between a cluster endpoint and a free element (or between two
free elements), a unit ascent is possible.

This is getting complex. Let me just verify the NUD formula first.
"""
from math import factorial, comb
from fractions import Fraction as F

# NUD(n) formula verification
def nud_alekseyev(n):
    """NUD(n) via Alekseyev-type formula (only descending clusters)."""
    total = factorial(n)  # k=0 term
    for k in range(1, n):
        for t in range(1, k+1):
            term = (-1)**k * comb(k-1, t-1) * comb(n-k, t) * factorial(n-k)
            total += term
    return total

# Known NUD values
nud_known = {1: 1, 2: 1, 3: 3, 4: 11, 5: 53, 6: 309, 7: 2119, 8: 16687,
             9: 148329, 10: 1468457, 11: 16019531, 12: 190899411}

print("=== Verify NUD Alekseyev formula ===")
for n in range(1, 13):
    computed = nud_alekseyev(n)
    known = nud_known[n]
    print(f"  n={n}: formula={computed}, known={known}, {'OK' if computed == known else 'FAIL'}")

# Good! Now for A002464, the formula uses 2^t instead of 1^t:
def hertzsprung_alekseyev(n):
    """A002464(n) via Alekseyev formula."""
    total = factorial(n)
    for k in range(1, n):
        for t in range(1, k+1):
            term = (-1)**k * comb(k-1, t-1) * comb(n-k, t) * 2**t * factorial(n-k)
            total += term
    return total

hertz_known = {1: 1, 2: 0, 3: 0, 4: 2, 5: 14, 6: 90, 7: 646, 8: 5242,
               9: 47622, 10: 479306, 11: 5296790, 12: 63779034}

print("\n=== Verify Hertzsprung Alekseyev formula ===")
for n in range(1, 13):
    computed = hertzsprung_alekseyev(n)
    known = hertz_known[n]
    print(f"  n={n}: formula={computed}, known={known}, {'OK' if computed == known else 'FAIL'}")

# Now the key question: what happens if we replace 2^t with (1+u)^t?
# With u=0: we get NUD (forbid descents, don't mark ascents) — should give A000255
# With u=-1: we get... perms with no descents AND marked ascents removed?
# Actually, let me think about this differently.
#
# In the Alekseyev formula:
# - We IE over "adjacency positions" (positions where |σ(i+1)-σ(i)|=1)
# - Each cluster of t consecutive adjacency positions contributes 2^t because
#   each cluster can be ascending or descending
#
# For NUD: we IE over "descent adjacency positions" only
# - Each cluster contributes 1^t = 1 (only descending direction)
#
# For W_u: we want to mark ascending adjacencies with u, while IE-ing over descending
# This is NOT the same as replacing 2^t with (1+u)^t, because:
# - The (1+u)^t replacement would mean: for each cluster, choose ascending (weight u)
#   or descending (weight 1). But we're NOT including-excluding over ascending ones.

# Let me try a different approach. What if we compute:
# F(n, u) = Σ_σ u^{asc_adj(σ)} * [no desc_adj]
# Using IE on desc_adj:
# F(n, u) = Σ_{D} (-1)^|D| Σ_σ [desc at all pos in D] * u^{asc_adj(σ)}

# For a fixed descent-forcing set D, we need to compute:
# Σ_σ [desc at all D] * u^{asc_adj(σ outside D)}
# The ascent adjacencies at positions NOT in D's clusters contribute.

# This is tractable for small D sets. Let me compute it for small n.

# Actually, let me try the formula (1+u)^t instead of 2^t and see what happens.
def wu_formula(n, u):
    """Hypothetical: replace 2^t with (1+u)^t in Alekseyev."""
    total = F(0)
    for k in range(0, n):
        for t in range(max(1, k > 0), k+1):
            if k == 0:
                total += F(factorial(n))  # base term
                break
            term = F((-1)**k) * comb(k-1, t-1) * comb(n-k, t) * F(1+u)**t * factorial(n-k)
            total += term
    return total

print("\n=== Test: (1+u)^t formula ===")
for u in range(3):
    print(f"\n  u={u}:")
    for n in range(1, 13):
        val = wu_formula(n, u)
        print(f"    n={n}: {val}")

# Check u=0 against NUD
print("\n=== u=0 should give NUD ===")
for n in range(1, 10):
    val = wu_formula(n, 0)
    print(f"  n={n}: formula={val}, NUD={nud_known[n]}, {'OK' if val == nud_known[n] else 'FAIL'}")

# Check u=-1 against Hertzsprung... wait, that would give 0^t = 0 for t>=1.
# So u=-1 gives: just n! (all clusters vanish). That's wrong for Hertzsprung.

# Actually, the Hertzsprung formula ADDS another factor. Let me reconsider.
# Hertzsprung = perms with no ascending adj AND no descending adj
# = IE over BOTH types of adjacencies
# The Alekseyev formula for Hertzsprung uses 2^t because each of the t clusters
# in the IE can represent either type.

# For our problem: NUD = IE only over descending adj. Correct, verified above.
# W_u = NUD weighted by u^{ascending adj} = need to additionally track ascending adj.

# The (1+u)^t hypothesis says: replace each "descending cluster contributes 1"
# with "each position in the cluster could alternatively be ascending (contributing u)"
# But this conflates two different things.

# Let me check numerically whether (1+u)^t gives W(n).
# Need N(n,j) data
Nnj = {
    1: [1], 2: [0, 1], 3: [0, 2, 1], 4: [2, 5, 3, 1],
    5: [14, 20, 14, 4, 1], 6: [90, 115, 72, 26, 5, 1],
    7: [646, 790, 467, 168, 41, 6, 1],
    8: [5242, 6217, 3557, 1285, 319, 59, 7, 1],
    9: [47622, 55160, 30968, 11120, 2834, 536, 80, 8, 1],
    10: [479306, 545135, 301970, 107918, 27752, 5432, 830, 104, 9, 1],
}

def Wu_exact(n, u):
    d = Nnj[n]
    return sum(F(u)**j * d[j] for j in range(len(d)))

print("\n=== Compare (1+u)^t formula with exact W_u ===")
for u in [0, 1, 2, 3]:
    print(f"\n  u={u}:")
    for n in range(1, 11):
        formula = wu_formula(n, u)
        exact = Wu_exact(n, u)
        match = formula == exact
        if not match:
            print(f"    n={n}: formula={formula}, exact={exact}, DIFF={formula-exact}")
        else:
            print(f"    n={n}: {exact} OK")
