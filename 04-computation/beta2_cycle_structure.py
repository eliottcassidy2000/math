#!/usr/bin/env python3
"""
beta2_cycle_structure.py - Analyze ker(d_2) structure and 3-cycle coupling

KEY INSIGHT from computation:
- beta_2 = 0 is dimension-specific (beta_3 can be nonzero at n>=7)
- At dimension 2, Omega constraints involve NON-ALLOWED 1-PATHS
- In a tournament, these are EXACTLY the "reverse edges"
- Each constraint couples 2-paths sharing a non-allowed face

This script explores:
1. What do ker(d_2) generators look like?
2. Do they decompose by 3-cycle triples?
3. How does the constraint coupling graph work?

Author: kind-pasteur-2026-03-08-S42
"""
import sys, os
import numpy as np
from collections import defaultdict
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

_saved = sys.stdout
sys.stdout = open(os.devnull, 'w', encoding='utf-8')
from path_homology_v2 import (
    enumerate_allowed_paths, build_full_boundary_matrix,
    compute_omega_basis
)
sys.stdout = _saved


def build_adj(n, bits):
    A = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx): A[i][j] = 1
            else: A[j][i] = 1
            idx += 1
    return A


def get_ker_generators(A, n):
    """Get explicit ker(d_2|Omega_2) generators in A_2 coordinates."""
    paths2 = enumerate_allowed_paths(A, n, 2)
    paths1 = enumerate_allowed_paths(A, n, 1)
    omega2 = compute_omega_basis(A, n, 2, paths2, paths1)
    dim_O2 = omega2.shape[1] if omega2.ndim == 2 else 0

    if dim_O2 == 0:
        return [], paths2

    D2 = build_full_boundary_matrix(paths2, paths1)
    D2_omega = D2 @ omega2
    U, S, Vt = np.linalg.svd(D2_omega, full_matrices=True)
    rank = sum(s > 1e-8 for s in S)

    if rank == dim_O2:
        return [], paths2

    ker_basis_omega = Vt[rank:].T  # columns in Omega_2 coords
    ker_in_A2 = omega2 @ ker_basis_omega  # columns in A_2 coords

    generators = []
    for col in range(ker_in_A2.shape[1]):
        v = ker_in_A2[:, col]
        gen = {}
        for j, path in enumerate(paths2):
            if abs(v[j]) > 1e-10:
                gen[path] = v[j]
        generators.append(gen)

    return generators, paths2


# ============================================================
# PART 1: Explicit generators at n=4
# ============================================================
print("=" * 70)
print("ker(d_2|Omega_2) GENERATORS AT n=4")
print("=" * 70)

n = 4
for bits in range(1 << (n*(n-1)//2)):
    A = build_adj(n, bits)
    scores = tuple(sorted([sum(row) for row in A]))

    gens, paths2 = get_ker_generators(A, n)
    if len(gens) == 0:
        continue

    print(f"\nT#{bits} scores={scores}, ker(d_2) dim={len(gens)}")

    # Adjacency
    edges = []
    for i in range(n):
        for j in range(n):
            if A[i][j]:
                edges.append(f"{i}->{j}")
    print(f"  Edges: {', '.join(edges)}")

    for gi, gen in enumerate(gens):
        print(f"  Generator #{gi+1}:")
        for (a, b, c), coeff in sorted(gen.items()):
            ac_type = "TT" if A[a][c] else "NT"
            triple = tuple(sorted([a, b, c]))
            print(f"    ({a},{b},{c}): {coeff:+.4f} [{ac_type}] triple={{{triple[0]},{triple[1]},{triple[2]}}}")

    # Only show first 3
    if bits > 20:
        break


# ============================================================
# PART 2: Do generators only use 3-cycle triples? (n=5 check)
# ============================================================
print(f"\n{'='*70}")
print("DO GENERATORS ONLY USE 3-CYCLE TRIPLES?")
print("=" * 70)

n = 5
uses_transitive = 0
only_cycles = 0

for bits in range(1 << (n*(n-1)//2)):
    A = build_adj(n, bits)
    gens, paths2 = get_ker_generators(A, n)

    if len(gens) == 0:
        continue

    this_only_cycles = True
    for gen in gens:
        for (a, b, c) in gen:
            triple = tuple(sorted([a, b, c]))
            i, j, k = triple
            is_cycle = ((A[i][j] and A[j][k] and A[k][i]) or
                       (A[j][i] and A[i][k] and A[k][j]))
            if not is_cycle:
                this_only_cycles = False
                break
        if not this_only_cycles:
            break

    if this_only_cycles:
        only_cycles += 1
    else:
        uses_transitive += 1

print(f"\nn=5: Tournaments with ker(d_2)>0:")
print(f"  Generators use ONLY 3-cycle triples: {only_cycles}")
print(f"  Generators use some TRANSITIVE triples: {uses_transitive}")


# ============================================================
# PART 3: 3-cycle chain analysis
# ============================================================
print(f"\n{'='*70}")
print("3-CYCLE CHAIN ANALYSIS")
print("=" * 70)

# For a 3-cycle triple {a,b,c} with a->b->c->a, define:
# The 3 allowed 2-paths are: (a,b,c), (b,c,a), (c,a,b)
# In Omega_2, these are constrained:
#   Non-allowed pair (c,a) [a->c is wrong]: constraint on mediators b with c->b->a? No...
# Wait, non-allowed 1-path = ordered pair (x,y) where y->x (edge goes "wrong way").
# For 3-cycle a->b->c->a:
# Non-allowed: (b,a) since a->b, (c,b) since b->c, (a,c) since c->a.
# For (b,a): mediators d with b->d and d->a. In {a,b,c}: d=c works? b->c yes, c->a yes. YES.
# But d can also be OUTSIDE {a,b,c}! This is key.

# For a larger tournament, the constraint for (b,a) involves ALL d in V\{a,b} with b->d->a.
# Some d are in the 3-cycle {a,b,c} (d=c), some are external.

# Let's count: for each non-allowed pair, how many mediators are in-cycle vs external?
print("\nn=5 regular tournament: mediator analysis")
for bits in [76]:  # first regular
    A = build_adj(n, bits)
    scores = tuple(sorted([sum(row) for row in A]))
    print(f"  T#{bits} scores={scores}")

    # Find 3-cycles
    cycles = []
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if A[i][j] and A[j][k] and A[k][i]:
                    cycles.append((i, j, k))
                elif A[j][i] and A[i][k] and A[k][j]:
                    cycles.append((j, i, k))  # directed: j->i->k->j

    print(f"  3-cycles (directed): {cycles}")

    # For each non-allowed pair, show mediators
    for a in range(n):
        for c in range(n):
            if a != c and A[c][a] == 1:
                meds = [b for b in range(n) if b != a and b != c
                        and A[a][b] == 1 and A[b][c] == 1]
                if len(meds) > 0:
                    # Which mediators are in a shared 3-cycle?
                    in_cycle_meds = []
                    ext_meds = []
                    for b in meds:
                        triple = tuple(sorted([a, b, c]))
                        is_3cycle = any(set(cy) == set(triple) for cy in cycles)
                        if is_3cycle:
                            in_cycle_meds.append(b)
                        else:
                            ext_meds.append(b)
                    print(f"    ({a},{c}): meds={meds}, in-cycle={in_cycle_meds}, external={ext_meds}")


# ============================================================
# PART 4: The filling map — explicit at n=4
# ============================================================
print(f"\n{'='*70}")
print("EXPLICIT FILLING MAP AT n=4")
print("=" * 70)

# For n=4 with score (1,1,2,2), ker(d_2) has dim 1.
# We need to find w in Omega_3 with d_3(w) = z.
# Let's do this explicitly.

n = 4
for bits in range(1 << (n*(n-1)//2)):
    A = build_adj(n, bits)
    scores = tuple(sorted([sum(row) for row in A]))
    if scores != (1, 1, 2, 2):
        continue

    gens, paths2 = get_ker_generators(A, n)
    if len(gens) == 0:
        continue

    paths3 = enumerate_allowed_paths(A, n, 3)
    omega3 = compute_omega_basis(A, n, 3, paths3, paths2)
    dim_O3 = omega3.shape[1] if omega3.ndim == 2 else 0

    print(f"\nT#{bits} scores={scores}")
    print(f"  dim(Omega_3) = {dim_O3}")
    print(f"  Allowed 3-paths: {paths3}")

    # The cycle generator
    z = gens[0]
    print(f"  z = {z}")

    # Build d_3 restricted to Omega_3
    D3 = build_full_boundary_matrix(paths3, paths2)
    if dim_O3 > 0:
        D3_omega = D3 @ omega3

        # z as a vector in A_2 coords
        path2_idx = {p: i for i, p in enumerate(paths2)}
        z_vec = np.zeros(len(paths2))
        for path, coeff in z.items():
            z_vec[path2_idx[path]] = coeff

        # Solve D3_omega @ w = z_vec (for w in Omega_3 coords)
        w, residual, _, _ = np.linalg.lstsq(D3_omega, z_vec, rcond=None)
        error = np.max(np.abs(D3_omega @ w - z_vec))

        if error < 1e-8:
            print(f"  FILLING FOUND! w in Omega_3:")
            w_in_A3 = omega3 @ w
            for j, path in enumerate(paths3):
                if abs(w_in_A3[j]) > 1e-10:
                    print(f"    w[{path}] = {w_in_A3[j]:.4f}")
        else:
            print(f"  NO filling in Omega_3! error={error:.6f}")

    break


print("\n\nDone.")
