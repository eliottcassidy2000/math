#!/usr/bin/env python3
"""
CONNECTION: Hook Schur coefficients of U_T and transfer matrix M.

Irving-Omar (Prop 26): [s_{(i,1^{n-i})}] U_T = #{sigma : Des_T(sigma) = {i,...,n-1}}

For tournament T, Des_T(sigma) at position k means T has edge sigma_k -> sigma_{k+1}.

So h_i = #{permutations where first i-1 consecutive pairs are T^op-arcs
           and last n-i consecutive pairs are T-arcs}

This is: #{orderings where the "break point" from backward to forward is at position i}

Our M[a,b] splits the vertex set into left (ending at a) and right (starting at b).
The hook coefficient splits the permutation at a POSITION.

QUESTION: Is there a formula connecting sum of h_i to traces of M or M^k?

kind-pasteur-2026-03-06-S25 (continuation)
"""

from itertools import permutations
import numpy as np

def make_circulant(n, S):
    T = {}
    for i in range(n):
        for j in range(n):
            if i != j:
                T[(i,j)] = 1 if ((j - i) % n) in S else 0
    return T


def compute_hook_schur(T, n):
    """Compute h_i = [s_{(i,1^{n-i})}] U_T for i = 1,...,n."""
    h = [0] * (n + 1)  # h[1],...,h[n]
    for sigma in permutations(range(n)):
        # Compute Des_T(sigma) = {k : T has sigma[k] -> sigma[k+1]}
        des = set()
        for k in range(n - 1):
            if T.get((sigma[k], sigma[k + 1]), 0) == 1:
                des.add(k + 1)  # 1-indexed positions

        # Check if des = {i, i+1, ..., n-1} for some i
        if len(des) == 0:
            h[n] += 1  # Des = empty, corresponds to s_{(n)}
        else:
            min_des = min(des)
            max_des = max(des)
            if max_des == n - 1 and len(des) == n - 1 - min_des + 1:
                # des = {min_des, min_des+1, ..., n-1}
                h[min_des] += 1
    return h


def compute_M_entry(T, n, a, b):
    if a == b:
        val = 0
        for perm in permutations(range(n)):
            prod = 1
            for k in range(n - 1):
                prod *= T.get((perm[k], perm[k + 1]), 0)
            if prod != 0:
                pos = list(perm).index(a)
                val += (-1) ** pos * prod
        return val
    else:
        U = [v for v in range(n) if v != a and v != b]
        val = 0
        for mask in range(1 << len(U)):
            S = [U[k] for k in range(len(U)) if mask & (1 << k)]
            R = [U[k] for k in range(len(U)) if not (mask & (1 << k))]
            sign = (-1) ** len(S)
            S_set = set(S) | {a}
            R_set = set(R) | {b}
            ea = 0
            for p in permutations(sorted(S_set)):
                if p[-1] != a: continue
                prod = 1
                for k in range(len(p) - 1):
                    prod *= T.get((p[k], p[k + 1]), 0)
                ea += prod
            if len(S_set) == 1: ea = 1
            bb2 = 0
            for p in permutations(sorted(R_set)):
                if p[0] != b: continue
                prod = 1
                for k in range(len(p) - 1):
                    prod *= T.get((p[k], p[k + 1]), 0)
                bb2 += prod
            if len(R_set) == 1: bb2 = 1
            val += sign * ea * bb2
        return val


# ============================================================
# Compute hook Schur coefficients for various tournaments
# ============================================================
print("=" * 70)
print("Hook Schur coefficients [s_{(i,1^{n-i})}] U_T")
print("=" * 70)

for n, name, S_set in [(3, "T_3 cyclic", {1}),
                        (5, "Paley T_5", {1, 2}),
                        (5, "circ {1,3}", {1, 3}),
                        (7, "Paley T_7", {1, 2, 4})]:
    T = make_circulant(n, S_set)
    h = compute_hook_schur(T, n)
    H = h[1]  # = ham(T)

    print(f"\n  {name} (n={n}): H = {H}")
    print(f"    Hook coefficients: h = {h[1:]}")
    print(f"    Palindromic: {all(h[i] == h[n + 1 - i] for i in range(1, n + 1))}")
    print(f"    Sum h_i = {sum(h[1:])}")
    print(f"    n! = {int(np.prod(np.arange(1, n+1)))}")

    # Sum h_i should equal... what? Not n!, because most permutations
    # don't have descent set of the form {i,...,n-1}
    # But sum over ALL Schur coefficients = sum over all descent sets = n!

    # For VT tournament: h_i should be constant (independent of i)?
    # No, palindromic but not constant.

    # Check: what's the alternating sum sum (-1)^i h_i?
    alt_sum = sum((-1) ** (i - 1) * h[i] for i in range(1, n + 1))
    print(f"    Alternating sum: {alt_sum}")

    # What about sum i * h_i?
    weighted = sum(i * h[i] for i in range(1, n + 1))
    print(f"    sum i*h_i = {weighted}")


# ============================================================
# Transitive tournament: hook Schur vs transfer matrix
# ============================================================
print("\n" + "=" * 70)
print("Transitive tournament: hook Schur + transfer matrix")
print("=" * 70)

for n in [3, 5]:
    T = {}
    for i in range(n):
        for j in range(n):
            if i != j:
                T[(i, j)] = 1 if i < j else 0

    h = compute_hook_schur(T, n)
    H = h[1]

    M = np.zeros((n, n), dtype=int)
    for a in range(n):
        for b in range(n):
            M[a, b] = compute_M_entry(T, n, a, b)

    print(f"\n  Transitive T_{n}: H = {H}")
    print(f"    Hook coefficients: h = {h[1:]}")
    print(f"    M =")
    for row in M:
        print(f"      {list(row)}")
    print(f"    tr(M) = {int(np.trace(M))}")
    print(f"    Eigenvalues of M: {sorted(np.linalg.eigvalsh(M))[::-1]}")


# ============================================================
# KEY TEST: Does M^2 or M eigenvalues relate to hook coefficients?
# ============================================================
print("\n" + "=" * 70)
print("Transfer matrix eigenvalues vs hook coefficients")
print("=" * 70)

for n, name, S_set in [(5, "Paley T_5", {1, 2}),
                        (5, "circ {1,3}", {1, 3})]:
    T = make_circulant(n, S_set)
    h = compute_hook_schur(T, n)

    M = np.zeros((n, n), dtype=int)
    for a in range(n):
        for b in range(n):
            M[a, b] = compute_M_entry(T, n, a, b)

    evals = sorted(np.linalg.eigvalsh(M))[::-1]

    print(f"\n  {name} (n={n}):")
    print(f"    Hook: h = {h[1:]}")
    print(f"    M eigenvalues: {[round(e, 4) for e in evals]}")
    print(f"    M = {int(M[0,0])} * I (scalar)")

    # For scalar M, all eigenvalues = H/n
    # The hook coefficients don't seem to directly give eigenvalues

    # But: the hook GENERATING function
    # sum h_i x^i should relate to the transfer matrix somehow

    # For Paley T_5 with H=15:
    # h = [15, 5, 5, 5, 15] (if palindromic)
    # Actually let's check

    print(f"    h generating polynomial: ", end="")
    print(" + ".join(f"{h[i]}*x^{i}" for i in range(1, n + 1) if h[i] != 0))


# ============================================================
# Test non-VT tournament
# ============================================================
print("\n" + "=" * 70)
print("Non-VT regular tournament: hook Schur + M")
print("=" * 70)

import random
random.seed(42)
n = 5
pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]

for trial in range(100000):
    T = {}
    for (i, j) in pairs:
        if random.random() < 0.5:
            T[(i, j)] = 1
            T[(j, i)] = 0
        else:
            T[(i, j)] = 0
            T[(j, i)] = 1

    if not all(sum(T.get((i, j), 0) for j in range(n) if j != i) == 2 for i in range(n)):
        continue

    # Check not circulant
    scores = sorted([sum(T.get((i, j), 0) for j in range(n) if j != i) for i in range(n)])
    if scores != [2, 2, 2, 2, 2]:
        continue

    h = compute_hook_schur(T, n)
    M = np.zeros((n, n), dtype=int)
    for a in range(n):
        for b in range(n):
            M[a, b] = compute_M_entry(T, n, a, b)

    is_scalar = np.allclose(M, M[0, 0] * np.eye(n))

    print(f"  H={h[1]}, hook={h[1:]}, M[0,0]={M[0,0]}, scalar={is_scalar}")
    break


# ============================================================
# DEEPER: Break point distribution and M
# ============================================================
print("\n" + "=" * 70)
print("Break point vs endpoints: joint distribution")
print("=" * 70)

# For each Ham path in T, record (start, end, break_points)
# A "break point" at position k means the arc at position k is a T^op-arc
# (i.e., sigma[k-1] -> sigma[k] is NOT in T, rather sigma[k] -> sigma[k-1] is)

n = 5
T = make_circulant(n, {1, 2})
print(f"\n  Paley T_5:")

# Compute for each permutation that is a Ham path:
# the set of "ascent positions" (where T has reverse arc)
endpoint_break = {}
for sigma in permutations(range(n)):
    prod = 1
    for k in range(n - 1):
        prod *= T.get((sigma[k], sigma[k + 1]), 0)
    if prod == 0:
        continue

    start = sigma[0]
    end = sigma[-1]
    # All steps are T-arcs (since prod=1), so Des = {1,...,n-1}
    # The hook coefficient h_1 = H

# Hmm, for a Ham path in T, ALL consecutive pairs are T-arcs by definition
# So Des_T(sigma) = {1,2,...,n-1} for every Ham path sigma
# Therefore h_1 = H, and h_i = 0 for i > 1... but that contradicts palindromicity!

# Wait, the hook coefficients count ALL permutations, not just Ham paths!
# Most permutations are NOT Ham paths (they use some reverse arcs)

print("  Clarification: h_i counts ALL permutations, not just Ham paths.")
print("  h_1 = #{sigma with Des_T = {1,...,n-1}} = #{Ham paths in T} = H")
print("  h_n = #{sigma with Des_T = empty} = #{Ham paths in T^op} = H")
print("  h_i (1 < i < n) counts permutations with mixed descent structure")

# So the hook coefficients refine n! (total permutations), not H
# h_1 = h_n = H, and sum h_i < n! typically (most perms don't have
# descent set of the form {i,...,n-1})

# The connection to M is more subtle.
# M[a,a] = sum over Ham paths of (-1)^{pos(a)}
# h_1 = sum over Ham paths of 1

# For the transfer matrix, we need to work with WEIGHTED permutations
# (weight = product of T[sigma_k, sigma_{k+1}])

# In the c-tournament generalization: t_{ij} = r + s_{ij}
# The hook coefficients would become polynomials in r and s

print("\n  INSIGHT: In the c-tournament generalization,")
print("  h_i(r, s) = sum_{sigma with Des_T = {i,...,n-1}} prod t(sigma_k, sigma_{k+1})")
print("  This is a polynomial in r and the s_{ij} variables.")
print("  h_1 = H(T) (total weighted Ham path count)")
print("  The EVEN R-POWERS property should apply to h_i as well!")
print("  Because h_i = sum of products, each using all n-1 edges")


print("\n" + "=" * 70)
print("DONE")
print("=" * 70)
