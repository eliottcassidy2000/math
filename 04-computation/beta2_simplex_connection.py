#!/usr/bin/env python3
"""
beta2_simplex_connection.py - Connection between transitive tournament and simplex

OBSERVATION: For the transitive tournament T_n (i->j iff i<j):
  dim(Om_p) = C(n, p+1) for all p

This means the path chain complex is ISOMORPHIC to the simplicial
chain complex of the (n-1)-simplex Delta^{n-1}!

The simplex is contractible, so ALL betti numbers are 0 (except beta_0=1).

For general tournaments, arc flips deform this chain complex.
The question: which deformations preserve exactness at dim 2?

Key idea: The simplex chain complex has the property that
  rk(d_{p+1}) + rk(d_p) = dim(C_p) for all p >= 1
  (exactness everywhere except dim 0)

For tournaments, we've shown exactness holds at dim 2 always,
but can fail at dim 1 (beta1 in {0,1}) and dim 3 (beta3 >= 0 at n>=6).

QUESTION: Is there a "homotopy" argument that preserves
exactness at dim 2 while allowing non-exactness elsewhere?

Author: kind-pasteur-2026-03-08-S41
"""
import sys, os, time
import numpy as np
from collections import Counter, defaultdict
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

_saved = sys.stdout
sys.stdout = open(os.devnull, 'w', encoding='utf-8')
from path_homology_v2 import (
    enumerate_allowed_paths, compute_omega_basis,
    build_full_boundary_matrix, path_betti_numbers
)
sys.stdout = _saved


def build_adj(n, bits):
    A = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx):
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    return A


def transitive_adj(n):
    """Build the transitive tournament: i->j iff i<j."""
    return [[1 if j > i else 0 for j in range(n)] for i in range(n)]


# ============================================================
# Verify: transitive tournament = simplex chain complex
# ============================================================
print("=" * 70)
print("TRANSITIVE TOURNAMENT = SIMPLEX?")
print("=" * 70)

from math import comb

for n in range(3, 9):
    A = transitive_adj(n)
    betti = path_betti_numbers(A, n, max_dim=min(n-1, 6))

    # Compute Om dims
    a_prev = [(i,) for i in range(n)]
    om_dims = [n]
    for p in range(1, min(n, 7)):
        ap = enumerate_allowed_paths(A, n, p)
        om = compute_omega_basis(A, n, p, ap, a_prev)
        d = om.shape[1] if om.ndim == 2 else 0
        om_dims.append(d)
        a_prev = ap

    expected = [comb(n, p+1) for p in range(len(om_dims))]

    match = om_dims == expected
    print(f"n={n}: Om = {om_dims}, C(n,p+1) = {expected}, match = {match}")
    print(f"       betti = {betti}")

print(f"\nThe transitive tournament on n vertices has path chain complex")
print(f"ISOMORPHIC to the (n-1)-simplex chain complex!")
print(f"All betti numbers = 0 except beta_0 = 1.")


# ============================================================
# How does ONE arc flip change the dimensions?
# ============================================================
print(f"\n{'='*70}")
print("ARC FLIP FROM TRANSITIVE: DIMENSION CHANGES")
print("=" * 70)

n = 5
A_trans = transitive_adj(n)

print(f"Transitive T_5: Om = ", end="")
a_prev = [(i,) for i in range(n)]
om_dims_trans = [n]
for p in range(1, 5):
    ap = enumerate_allowed_paths(A_trans, n, p)
    om = compute_omega_basis(A_trans, n, p, ap, a_prev)
    d = om.shape[1] if om.ndim == 2 else 0
    om_dims_trans.append(d)
    a_prev = ap
print(om_dims_trans)

# Flip each possible arc and see what happens
print(f"\nAfter one flip:")
for u in range(n):
    for v in range(u+1, n):
        if not A_trans[u][v]:
            continue
        A_flip = [row[:] for row in A_trans]
        A_flip[u][v] = 0
        A_flip[v][u] = 1
        scores = tuple(sorted([sum(row) for row in A_flip]))

        a_prev = [(i,) for i in range(n)]
        om_dims = [n]
        for p in range(1, 5):
            ap = enumerate_allowed_paths(A_flip, n, p)
            om = compute_omega_basis(A_flip, n, p, ap, a_prev)
            d = om.shape[1] if om.ndim == 2 else 0
            om_dims.append(d)
            a_prev = ap

        betti = path_betti_numbers(A_flip, n, max_dim=4)
        delta = [om_dims[p] - om_dims_trans[p] for p in range(len(om_dims))]
        print(f"  Flip {u}->{v}: Om={om_dims}, delta={delta}, betti={betti}")


# ============================================================
# Pattern: how do Om dimensions change under repeated flips?
# ============================================================
print(f"\n{'='*70}")
print("OMEGA DIMENSION LANDSCAPE AT n=5")
print("=" * 70)

n = 5
total = 1 << (n*(n-1)//2)

om_patterns = Counter()
for bits in range(total):
    A = build_adj(n, bits)
    a_prev = [(i,) for i in range(n)]
    om_dims = [n]
    for p in range(1, 5):
        ap = enumerate_allowed_paths(A, n, p)
        om = compute_omega_basis(A, n, p, ap, a_prev)
        d = om.shape[1] if om.ndim == 2 else 0
        om_dims.append(d)
        a_prev = ap
    om_patterns[tuple(om_dims)] += 1

print(f"Distinct Om dimension patterns at n=5: {len(om_patterns)}")
for pat in sorted(om_patterns.keys()):
    cnt = om_patterns[pat]
    chi = sum((-1)**p * pat[p] for p in range(len(pat)))
    print(f"  {list(pat)}: {cnt} tournaments, chi={chi}")


# ============================================================
# KEY IDENTITY: Does dim(Om_2) = C(n-1,2) + c3 - <something>?
# We showed Om_2 = |A_2| - n_bp
# What is |A_2|? And n_bp?
# ============================================================
print(f"\n{'='*70}")
print("DECOMPOSITION OF Om_2")
print("=" * 70)

n = 5
for bits in [0, 341, 10]:
    A = build_adj(n, bits)
    scores = tuple(sorted([sum(row) for row in A]))
    c3 = sum(1 for i in range(n) for j in range(i+1,n) for k in range(j+1,n)
             if (A[i][j] and A[j][k] and A[k][i]) or (A[j][i] and A[k][j] and A[i][k]))

    a2 = enumerate_allowed_paths(A, n, 2)

    # Count by type
    tt = sum(1 for p in a2 if A[p[0]][p[2]])  # transitive triple (a->c exists)
    jk = sum(1 for p in a2 if A[p[2]][p[0]])  # junk (c->a)
    print(f"bits={bits}, scores={scores}, c3={c3}")
    print(f"  |A2|={len(a2)}, TT={tt}, Junk={jk}")
    print(f"  |A2| = sum d_b^out * d_b^in = sum d_b * (n-1-d_b)")

    # The number of backward pairs with intermediaries
    bp = 0
    for a in range(n):
        for c in range(n):
            if a != c and A[c][a]:
                if any(A[a][b] and A[b][c] for b in range(n) if b != a and b != c):
                    bp += 1
    print(f"  n_bp={bp}")
    print(f"  Om2 = |A2| - n_bp = {len(a2)} - {bp} = {len(a2) - bp}")


# ============================================================
# EULER CHARACTERISTIC IDENTITY
# ============================================================
print(f"\n{'='*70}")
print("EULER CHARACTERISTIC ANALYSIS")
print("=" * 70)

# From the data above:
# chi = 1 when beta1 = 0
# chi = 0 when beta1 = 1
# Since beta2 = beta3 = beta4 = 0 at n=5:
# chi = 1 - beta1

# So: sum(-1)^p dim(Om_p) = 1 - beta1
# dim(Om_0) - dim(Om_1) + dim(Om_2) - dim(Om_3) + dim(Om_4) = 1 - beta1

# For transitive: 5 - 10 + 10 - 5 + 1 = 1 (beta1=0)
# For regular:    5 - 10 + 10 - 10 + 5 = 0 (beta1=1)

print(f"chi = 1 - beta1 at n=5 (since beta2=beta3=beta4=0)")
print(f"Transitive: 5 - 10 + 10 - 5 + 1 = 1")
print(f"Regular:    5 - 10 + 10 - 10 + 5 = 0")

# At n=6: beta3 can be 1, so chi = 1 - beta1 - beta3
# Let's verify at n=6
import random
random.seed(42)
n = 6

print(f"\nVerifying chi = 1 - beta1 - beta3 at n=6:")
checks = 0
failures = 0
for _ in range(200):
    bits = random.randint(0, (1 << (n*(n-1)//2)) - 1)
    A = build_adj(n, bits)

    betti = path_betti_numbers(A, n, max_dim=5)

    a_prev = [(i,) for i in range(n)]
    om_dims = [n]
    for p in range(1, 6):
        ap = enumerate_allowed_paths(A, n, p)
        om = compute_omega_basis(A, n, p, ap, a_prev)
        d = om.shape[1] if om.ndim == 2 else 0
        om_dims.append(d)
        a_prev = ap

    chi_dims = sum((-1)**p * om_dims[p] for p in range(len(om_dims)))
    chi_betti = sum((-1)**p * betti[p] for p in range(len(betti)))
    expected = 1 - betti[1] - betti[3] + (betti[4] if len(betti) > 4 else 0)

    if chi_dims != chi_betti:
        failures += 1
    checks += 1

print(f"  {checks} checks, {failures} chi mismatches")

# What IS beta1 in terms of tournament structure?
print(f"\n\nbeta1 analysis at n=5:")
n = 5
beta1_data = defaultdict(list)
for bits in range(1 << (n*(n-1)//2)):
    A = build_adj(n, bits)
    scores = tuple(sorted([sum(row) for row in A]))
    c3 = sum(1 for i in range(n) for j in range(i+1,n) for k in range(j+1,n)
             if (A[i][j] and A[j][k] and A[k][i]) or (A[j][i] and A[k][j] and A[i][k]))
    betti = path_betti_numbers(A, n, max_dim=1)
    beta1_data[(scores, c3)].append(betti[1])

for (sc, c3) in sorted(beta1_data.keys()):
    vals = Counter(beta1_data[(sc, c3)])
    print(f"  {sc}, c3={c3}: beta1 = {dict(sorted(vals.items()))}")


print("\n\nDone.")
