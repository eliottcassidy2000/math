#!/usr/bin/env python3
"""
Connection between sigma-algebra hierarchy and path homology Betti numbers.

Questions:
1. Is β_1 determined by lambda? (β_1 depends on c3 structure, which IS lambda-determined)
2. Is β_3 determined by lambda? (β_3 depends on 7-cycle structure — maybe not?)
3. Is β_3 determined by (lambda, sigma)?

We know:
- β_1 ∈ {0, 1} (HYP-320)
- β_2 = 0 always (HYP-322)
- β_3 ∈ {0, 1} for n ≤ 8 (observed)
- β_1 * β_3 = 0 (seesaw, HYP-316)

opus-2026-03-13-S71c
"""
import sys, time
import numpy as np
from itertools import combinations, permutations
from collections import defaultdict
sys.stdout.reconfigure(line_buffering=True)

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

def lambda_sigma(A, n):
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
    return L, S

def compute_betti(A, n, max_d=3):
    """Compute Betti numbers β_0, β_1, ..., β_{max_d} for tournament with adjacency A."""
    # Enumerate allowed d-paths
    # A d-path is a sequence (v_0,...,v_d) where v_i→v_{i+1} for all i
    # and all vertices are distinct.
    # "Allowed" = no two consecutive arrows define a total ordering on any 3 vertices.
    # Actually for GLMY path homology: allowed = elementary (all vertices distinct).

    # For tournaments, every path on distinct vertices is allowed.
    # The chain spaces: C_d = span of d-paths (v_0,...,v_d) with all v_i distinct and v_i→v_{i+1}
    # Omega_d = ker(constraint: faces must also be allowed or cancel)

    # For simplicity, compute using the boundary matrix approach.
    # ∂_d: C_d → C_{d-1} with ∂(v_0,...,v_d) = Σ_i (-1)^i d_i(v_0,...,v_d)
    # where d_i removes or merges vertex i.

    # For GLMY:
    # d_0(v_0,...,v_d) = (v_1,...,v_d)
    # d_d(v_0,...,v_d) = (v_0,...,v_{d-1})
    # d_i for 0<i<d: (v_0,...,v_{i-1}·v_i,...,v_d) where v_{i-1}·v_i means
    #   the path skipping v_i: v_{i-1} → v_{i+1} directly.
    #   This is only valid if v_{i-1} → v_{i+1} exists.

    # Actually for TOURNAMENT path homology, the simplicial version:
    # d_i just removes vertex v_i from the sequence.
    # The resulting path (v_0,...,v̂_i,...,v_d) is a (d-1)-path.
    # It's "allowed" iff v_{i-1}→v_{i+1} (for interior faces) or it just truncates (face 0 or d).

    # Let me use the standard approach from path_homology_v2.py.
    # Enumerate all d-paths with distinct vertices where each consecutive pair has an edge.

    paths = {0: [(v,) for v in range(n)]}
    for d in range(1, max_d + 2):
        pd = []
        for p in paths[d-1]:
            last = p[-1]
            used = set(p)
            for v in range(n):
                if v not in used and A[last][v] == 1:
                    pd.append(p + (v,))
        paths[d] = pd
        if not pd:
            break

    actual_max = max(d for d in paths if paths[d])

    # Build boundary matrices
    # For degree d, boundary ∂_d: C_d → C_{d-1}
    # Faces: d_i(v_0,...,v_d) = (v_0,...,v̂_i,...,v_d) with sign (-1)^i
    # A face is "allowed" if the resulting sequence is a valid (d-1)-path.

    # For the constraint/Omega approach:
    # Ω_d = {chains in C_d such that ∂c has only allowed faces}
    # i.e., the "junk" faces (non-allowed faces) must cancel.

    path_idx = {}
    for d in range(actual_max + 1):
        path_idx[d] = {p: i for i, p in enumerate(paths[d])}

    allowed_set = {d: set(paths[d]) for d in range(actual_max + 1)}

    betti = []
    omega_dims = [n]  # Ω_0 = C_0 = n vertices
    bd_ranks = [0]    # rank(∂_0) = 0

    for d in range(1, min(max_d + 1, actual_max + 1)):
        n_d = len(paths[d])
        n_dm1 = len(paths[d-1])
        if n_d == 0:
            omega_dims.append(0)
            bd_ranks.append(0)
            continue

        # Compute faces
        junk_set = set()
        face_junk = []  # per path in paths[d]
        face_allowed = []

        for p in paths[d]:
            jf = []
            af = []
            for fi in range(d + 1):
                face = p[:fi] + p[fi+1:]
                sign = 1 if fi % 2 == 0 else -1
                # Check if face is an allowed (d-1)-path
                if face in allowed_set[d-1]:
                    af.append((face, sign))
                else:
                    junk_set.add(face)
                    jf.append((face, sign))
            face_junk.append(jf)
            face_allowed.append(af)

        junk_list = sorted(junk_set)
        n_junk = len(junk_list)
        junk_idx = {j: i for i, j in enumerate(junk_list)}

        # Constraint matrix (junk rows × path cols)
        C = np.zeros((n_junk, n_d), dtype=float)
        for j, jf in enumerate(face_junk):
            for face, sign in jf:
                C[junk_idx[face], j] += sign

        rank_c = int(np.linalg.matrix_rank(C)) if n_junk > 0 else 0
        omega_d = n_d - rank_c
        omega_dims.append(omega_d)

        # Combined matrix (junk + allowed rows × path cols)
        CB = np.zeros((n_junk + n_dm1, n_d), dtype=float)
        CB[:n_junk, :] = C
        for j, af in enumerate(face_allowed):
            for face, sign in af:
                row = n_junk + path_idx[d-1][face]
                CB[row, j] += sign

        rank_cb = int(np.linalg.matrix_rank(CB))
        bd_rank = rank_cb - rank_c
        bd_ranks.append(bd_rank)

    # Compute Betti
    for d in range(min(max_d + 1, len(omega_dims))):
        od = omega_dims[d]
        rd = bd_ranks[d] if d < len(bd_ranks) else 0
        rd1 = bd_ranks[d+1] if d+1 < len(bd_ranks) else 0
        betti.append(od - rd - rd1)

    return betti

# Test at n=7
n = 7
tb = n*(n-1)//2
np.random.seed(42)

print(f"n={n}: checking if Betti numbers are determined by lambda/sigma")
print(f"Computing Betti through β_3 for sampled tournaments...")

lam_betti = defaultdict(set)
sig_betti = defaultdict(set)

t0 = time.time()
for trial in range(2000):
    bits = np.random.randint(0, 1 << tb)
    A = bits_to_adj(bits, n)
    L, S = lambda_sigma(A, n)

    betti = compute_betti(A, n, max_d=3)
    betti_key = tuple(betti)

    lam_key = tuple(L[i][j] for i in range(n) for j in range(i+1, n))
    sig_key = lam_key + tuple(S[i][j] for i in range(n) for j in range(i+1, n))

    lam_betti[lam_key].add(betti_key)
    sig_betti[sig_key].add(betti_key)

    if trial % 500 == 0:
        dt = time.time() - t0
        print(f"  trial {trial}: {dt:.1f}s")

dt = time.time() - t0
print(f"\nDone: {len(lam_betti)} lambda groups, {len(sig_betti)} sigma groups, {dt:.1f}s")

lam_ambig = sum(1 for v in lam_betti.values() if len(v) > 1)
sig_ambig = sum(1 for v in sig_betti.values() if len(v) > 1)

print(f"\nLambda determines Betti? Ambiguous: {lam_ambig}")
print(f"(Lambda,sigma) determines Betti? Ambiguous: {sig_ambig}")

if lam_ambig > 0:
    count = 0
    for key, vals in sorted(lam_betti.items()):
        if len(vals) > 1:
            print(f"  Lambda ambiguity: betti values = {sorted(vals)}")
            count += 1
            if count >= 5: break

# Check which Betti numbers individually are lambda-determined
for bd in range(4):
    lam_bd = defaultdict(set)
    for key, bvals in lam_betti.items():
        for bv in bvals:
            if bd < len(bv):
                lam_bd[key].add(bv[bd])
    ambig_bd = sum(1 for v in lam_bd.values() if len(v) > 1)
    print(f"  β_{bd} lambda-determined? Ambiguous: {ambig_bd}")

print(f"\nDone.")
