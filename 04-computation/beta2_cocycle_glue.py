#!/usr/bin/env python3
"""
COCYCLE GLUE ANALYSIS — Algebraic mechanism of flip obstruction in tournament path homology

KEY FACT (corrected): At n=5, beta_1=0 tournaments with exactly 3 bad vertices have those
3 bad vertices forming a single 3-cycle. The flip that creates beta_1=1 involves an arc
between a GOOD and BAD vertex, creating 3 NEW 3-cycles and making all 5 vertices bad.

This script performs:
1. COCYCLE ANALYSIS (cohomological): Z^1, B^1, restriction maps
2. AFTER THE FLIP: what changes when a single arc flip creates beta_1=1
3. GLUING QUESTION: net constraint changes in allowed paths and Omega
4. BOUNDARY ANALYSIS: is the original 3-cycle a boundary? what changes?
5. DIMENSION TRACKING: full chain complex dimensions before/after

Uses path_homology_v2.py core functions.
"""

import numpy as np
from itertools import permutations, combinations
from math import comb
from collections import defaultdict
import sys

# ========== CORE FUNCTIONS (from path_homology_v2.py) ==========

def enumerate_allowed_paths(A, n, p):
    if p < 0: return []
    if p == 0: return [(v,) for v in range(n)]
    adj = [[] for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if A[i][j] == 1: adj[i].append(j)
    paths = []
    stack = []
    for start in range(n):
        stack.append(([start], 1 << start))
        while stack:
            path, visited = stack.pop()
            if len(path) == p + 1:
                paths.append(tuple(path))
                continue
            v = path[-1]
            for u in adj[v]:
                if not (visited & (1 << u)):
                    stack.append((path + [u], visited | (1 << u)))
    return paths

def boundary_coeffs(path):
    return [((-1)**i, path[:i] + path[i+1:]) for i in range(len(path))]

def build_full_boundary_matrix(allowed_p, allowed_pm1):
    if not allowed_p or not allowed_pm1:
        return np.zeros((max(len(allowed_pm1), 0), max(len(allowed_p), 0)))
    idx_pm1 = {path: i for i, path in enumerate(allowed_pm1)}
    M = np.zeros((len(allowed_pm1), len(allowed_p)))
    for j, path in enumerate(allowed_p):
        for sign, face in boundary_coeffs(path):
            if face in idx_pm1: M[idx_pm1[face], j] += sign
    return M

def compute_omega_basis(A, n, p, allowed_p, allowed_pm1):
    dim_Ap = len(allowed_p)
    if dim_Ap == 0: return np.zeros((0, 0))
    if p == 0: return np.eye(dim_Ap)
    allowed_pm1_set = set(allowed_pm1)
    non_allowed_faces = {}
    na_count = 0
    for j, path in enumerate(allowed_p):
        for sign, face in boundary_coeffs(path):
            if len(set(face)) == len(face) and face not in allowed_pm1_set:
                if face not in non_allowed_faces:
                    non_allowed_faces[face] = na_count
                    na_count += 1
    if na_count == 0: return np.eye(dim_Ap)
    P = np.zeros((na_count, dim_Ap))
    for j, path in enumerate(allowed_p):
        for sign, face in boundary_coeffs(path):
            if face in non_allowed_faces: P[non_allowed_faces[face], j] += sign
    U, S, Vt = np.linalg.svd(P, full_matrices=True)
    rank = sum(s > 1e-10 for s in S)
    null_space = Vt[rank:].T
    if null_space.shape[1] == 0: return np.zeros((dim_Ap, 0))
    return null_space

def compute_full_chain_complex(A, n, max_dim=None):
    if max_dim is None: max_dim = n - 1
    allowed = {}
    for p in range(-1, max_dim + 2):
        allowed[p] = [] if p < 0 else enumerate_allowed_paths(A, n, p)
    omega = {}
    for p in range(max_dim + 2):
        omega[p] = compute_omega_basis(A, n, p, allowed[p], allowed[p-1])
    return {'allowed': allowed, 'omega': omega}

def compute_chain_dims(A, n, max_dim=3):
    """Compute dim(Omega_k), dim(Z_k), dim(B_k), beta_k for k=0..max_dim."""
    cc = compute_full_chain_complex(A, n, max_dim=max_dim)
    allowed = cc['allowed']
    omega = cc['omega']
    results = []
    for k in range(max_dim + 1):
        dim_Ok = omega[k].shape[1] if omega[k].ndim == 2 else 0
        # Z_k = ker(d_k on Omega_k)
        if dim_Ok > 0 and k > 0:
            bd_k = build_full_boundary_matrix(allowed[k], allowed[k-1])
            bd_k_o = bd_k @ omega[k]
            sv = np.linalg.svd(bd_k_o, compute_uv=False)
            rank_dk = sum(s > 1e-8 for s in sv)
            dim_Zk = dim_Ok - rank_dk
        elif k == 0:
            dim_Zk = dim_Ok
        else:
            dim_Zk = 0
        # B_k = im(d_{k+1} on Omega_{k+1})
        dim_Okp1 = omega[k+1].shape[1] if k+1 in omega and omega[k+1].ndim == 2 else 0
        if dim_Okp1 > 0:
            bd_kp1 = build_full_boundary_matrix(allowed[k+1], allowed[k])
            bd_kp1_o = bd_kp1 @ omega[k+1]
            sv2 = np.linalg.svd(bd_kp1_o, compute_uv=False)
            dim_Bk = sum(s > 1e-8 for s in sv2)
        else:
            dim_Bk = 0
        beta_k = max(0, dim_Zk - dim_Bk)
        results.append({'dim_A': len(allowed[k]), 'dim_O': dim_Ok,
                       'dim_Z': dim_Zk, 'dim_B': dim_Bk, 'beta': beta_k})
    return results, cc

def compute_cocycle_spaces(A, n, p):
    """Compute Z^p (cocycles) and B^p (coboundaries) in Omega coordinates.
    Returns: allowed_p, Z_p_basis, B_p_basis, H_p_dim, omega_p, allowed_dict
    """
    allowed = {}
    for q in range(max(p-2, -1), p+3):
        allowed[q] = [] if q < 0 else enumerate_allowed_paths(A, n, q)

    omega = {}
    for q in range(max(p-1, 0), p+2):
        omega[q] = compute_omega_basis(A, n, q, allowed[q], allowed[q-1])

    dim_omega_p = omega[p].shape[1] if omega[p].ndim == 2 and omega[p].shape[1] > 0 else 0
    if dim_omega_p == 0:
        return allowed[p], np.zeros((0,0)), np.zeros((0,0)), 0, np.zeros((0,0)), allowed

    # d_{p+1}^Omega: Omega_{p+1} -> A_p, then project to Omega_p coords
    dim_omega_pp1 = omega[p+1].shape[1] if p+1 in omega and omega[p+1].ndim == 2 else 0
    if dim_omega_pp1 > 0:
        bd_pp1 = build_full_boundary_matrix(allowed[p+1], allowed[p])
        bd_pp1_o = bd_pp1 @ omega[p+1]  # |A_p| x dim_O_{p+1}
        # Project into Omega_p coords
        d_pp1_in_Op = omega[p].T @ bd_pp1_o  # dim_O_p x dim_O_{p+1}
    else:
        d_pp1_in_Op = np.zeros((dim_omega_p, 0))

    # Z^p = ker(delta^p) where delta^p(f) = f . d_{p+1}^Omega
    # delta^p: (Omega_p)* -> (Omega_{p+1})*, f |-> d_{p+1}^{Omega,T} f
    # ker(delta^p) = nullspace of d_pp1_in_Op^T
    if d_pp1_in_Op.shape[1] > 0:
        U, S_vals, Vt = np.linalg.svd(d_pp1_in_Op.T, full_matrices=True)
        rank_delta = sum(s > 1e-8 for s in S_vals)
        Z_p = Vt[rank_delta:].T  # columns are basis of Z^p
    else:
        Z_p = np.eye(dim_omega_p)

    # B^p = im(delta^{p-1}) where delta^{p-1}(g) = g . d_p^Omega
    # d_p^Omega: Omega_p -> Omega_{p-1} in Omega coords
    dim_omega_pm1 = omega[p-1].shape[1] if p-1 in omega and omega[p-1].ndim == 2 else 0
    if dim_omega_pm1 > 0:
        bd_p = build_full_boundary_matrix(allowed[p], allowed[p-1])
        bd_p_o = bd_p @ omega[p]  # |A_{p-1}| x dim_O_p
        d_p_in_O = omega[p-1].T @ bd_p_o  # dim_O_{p-1} x dim_O_p
        # B^p = column space of d_p_in_O^T
        B_mat = d_p_in_O.T  # dim_O_p x dim_O_{p-1}
        U_B, S_B, Vt_B = np.linalg.svd(B_mat, full_matrices=False)
        rank_B = sum(s > 1e-8 for s in S_B)
        B_p = U_B[:, :rank_B]
    else:
        B_p = np.zeros((dim_omega_p, 0))

    dim_Zp = Z_p.shape[1] if Z_p.ndim == 2 else 0
    dim_Bp = B_p.shape[1] if B_p.ndim == 2 else 0

    return allowed[p], Z_p, B_p, dim_Zp - dim_Bp, omega[p], allowed

# ========== TOURNAMENT TOOLS ==========

def all_tournaments(n):
    edges = [(i,j) for i in range(n) for j in range(i+1,n)]
    m = len(edges)
    for mask in range(1 << m):
        A = [[0]*n for _ in range(n)]
        for idx, (i,j) in enumerate(edges):
            if (mask >> idx) & 1: A[i][j] = 1
            else: A[j][i] = 1
        yield [row[:] for row in A]

def find_3cycles(A, n):
    cycles = set()
    for i in range(n):
        for j in range(n):
            if j == i: continue
            for k in range(n):
                if k == i or k == j: continue
                if A[i][j] and A[j][k] and A[k][i]:
                    cycles.add(frozenset([i,j,k]))
    return [tuple(sorted(c)) for c in cycles]

def find_3cycles_directed(A, n):
    """Return all directed 3-cycles as (i,j,k) with i->j->k->i, i<j<k canonical."""
    result = []
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if A[i][j] and A[j][k] and A[k][i]:
                    result.append(((i,j,k), 'forward'))  # i->j->k->i
                elif A[j][i] and A[i][k] and A[k][j]:
                    result.append(((i,j,k), 'backward'))  # j->i->k->j = i<-j, i->k, k<-j? No
                    # Actually: A[j][i]=j->i, A[i][k]=i->k, A[k][j]=k->j => j->i->k->j
                    pass
    return result

def bad_vertices(A, n):
    bad = set()
    for c in find_3cycles(A, n):
        bad.update(c)
    return bad

def score_seq(A, n):
    return tuple(sorted([sum(A[i]) for i in range(n)]))

def count_3cycles(A, n):
    return len(find_3cycles(A, n))

def flip_arc(A, n, u, v):
    """Flip arc u->v to v->u. Assumes A[u][v]=1."""
    A2 = [row[:] for row in A]
    A2[u][v] = 0
    A2[v][u] = 1
    return A2

def subtournament(A, n, vertices):
    vlist = sorted(vertices)
    m = len(vlist)
    idx = {v: i for i, v in enumerate(vlist)}
    B = [[0]*m for _ in range(m)]
    for i in vlist:
        for j in vlist:
            if i != j and A[i][j]:
                B[idx[i]][idx[j]] = 1
    return B, vlist

# ========== FIND NON-COBOUNDARY COCYCLE ==========

def find_non_coboundary(Z_basis, B_basis, omega, allowed_p, label=""):
    """Find a representative non-coboundary cocycle and return it in A_p coordinates."""
    dim_Z = Z_basis.shape[1] if Z_basis.ndim == 2 else 0
    dim_B = B_basis.shape[1] if B_basis.ndim == 2 else 0
    if dim_Z <= dim_B:
        return None  # H^p = 0
    # Find vector in Z not in B
    for col in range(dim_Z):
        z = Z_basis[:, col:col+1]
        if dim_B > 0:
            proj = B_basis @ (B_basis.T @ z)
            w = z - proj
        else:
            w = z.copy()
        if np.linalg.norm(w) > 1e-8:
            w = w / np.linalg.norm(w)
            # Convert to A_p coordinates
            w_Ap = omega @ w
            return w, w_Ap
    return None


# ========== MAIN ANALYSIS ==========

def pr(s=""): print(s)
def sep(title, ch='=', w=78):
    pr(f"\n{ch*w}")
    pr(f"  {title}")
    pr(f"{ch*w}")

sep("COCYCLE GLUE ANALYSIS: Flip Obstruction in Tournament Path Homology")
pr("n=5, beta_1=0 tournaments with exactly 3 bad vertices")
pr("All 3 bad vertices form a SINGLE 3-cycle (always true at n=5)")
pr("Flip: one arc between good/bad vertex -> creates beta_1=1")
pr()

n = 5

# ========== PHASE 1: COLLECT ALL RELEVANT TOURNAMENTS ==========
pr("Phase 1: Collecting tournaments...")

tournaments_b0_3bad = []  # (A, bad_set)
for A in all_tournaments(n):
    dims, cc = compute_chain_dims(A, n, max_dim=1)
    b1 = dims[1]['beta']
    if b1 != 0:
        continue
    bad = bad_vertices(A, n)
    if len(bad) != 3:
        continue
    tournaments_b0_3bad.append(A)

pr(f"Found {len(tournaments_b0_3bad)} tournaments with beta_1=0 and 3 bad vertices")

# ========== PHASE 2: FIND ALL FLIPS THAT CREATE BETA_1=1 ==========
pr("\nPhase 2: Finding single-arc flips that create beta_1=1...")

flip_pairs = []  # (A_before, A_after, flipped_arc_old, flipped_arc_new)
seen_after = set()

for A in tournaments_b0_3bad:
    for i in range(n):
        for j in range(n):
            if i == j or not A[i][j]:
                continue
            A2 = flip_arc(A, n, i, j)
            dims2, _ = compute_chain_dims(A2, n, max_dim=1)
            if dims2[1]['beta'] == 1:
                # Canonicalize to avoid duplicates
                key = tuple(tuple(r) for r in A2)
                if key not in seen_after:
                    seen_after.add(key)
                    flip_pairs.append((A, A2, (i,j), (j,i)))

pr(f"Found {len(flip_pairs)} distinct (before, after) flip pairs")

# Group by isomorphism class (use score sequence as rough proxy)
by_score = defaultdict(list)
for item in flip_pairs:
    A, A2, old_arc, new_arc = item
    s1 = score_seq(A, n)
    s2 = score_seq(A2, n)
    by_score[(s1, s2)].append(item)

pr(f"Score sequence classes: {len(by_score)}")
for (s1, s2), items in sorted(by_score.items()):
    pr(f"  {s1} -> {s2}: {len(items)} pairs")

# ========== PHASE 3: DETAILED ANALYSIS ON REPRESENTATIVE EXAMPLES ==========
# Take first from each score class, and also a few extras

examples = []
for (s1, s2), items in sorted(by_score.items()):
    examples.append(items[0])

pr(f"\nAnalyzing {len(examples)} representative examples in detail...")

for ex_idx, (A, A2, old_arc, new_arc) in enumerate(examples):
    sep(f"EXAMPLE {ex_idx+1}: score {score_seq(A,n)} -> {score_seq(A2,n)}", ch='-', w=70)
    pr(f"  Flipped arc: {old_arc[0]}->{old_arc[1]}  became  {new_arc[0]}->{new_arc[1]}")

    bad_old = bad_vertices(A, n)
    bad_new = bad_vertices(A2, n)
    cycles_old = find_3cycles(A, n)
    cycles_new = find_3cycles(A2, n)
    pr(f"  3-cycles before: {cycles_old}  (bad vertices: {sorted(bad_old)})")
    pr(f"  3-cycles after:  {cycles_new}  (bad vertices: {sorted(bad_new)})")
    new_cycles = [c for c in cycles_new if c not in cycles_old]
    pr(f"  NEW 3-cycles created by flip: {new_cycles}")

    # ---------- PART 1: COCYCLE ANALYSIS ----------
    pr(f"\n  === PART 1: Cohomology (degree 1) ===")

    allowed_1, Z1, B1, H1, omega_1, all_allowed = compute_cocycle_spaces(A, n, 1)
    dim_Z1 = Z1.shape[1] if Z1.ndim == 2 else 0
    dim_B1 = B1.shape[1] if B1.ndim == 2 else 0
    pr(f"  T (before):  dim Z^1 = {dim_Z1},  dim B^1 = {dim_B1},  H^1 = {H1}")

    allowed_1p, Z1p, B1p, H1p, omega_1p, all_allowed_p = compute_cocycle_spaces(A2, n, 1)
    dim_Z1p = Z1p.shape[1] if Z1p.ndim == 2 else 0
    dim_B1p = B1p.shape[1] if B1p.ndim == 2 else 0
    pr(f"  T' (after):  dim Z^1 = {dim_Z1p},  dim B^1 = {dim_B1p},  H^1 = {H1p}")
    pr(f"  Change: dZ^1 = {dim_Z1p - dim_Z1:+d},  dB^1 = {dim_B1p - dim_B1:+d},  dH^1 = {H1p - H1:+d}")

    # Find the non-coboundary cocycle in T'
    ncb = find_non_coboundary(Z1p, B1p, omega_1p, allowed_1p)
    if ncb is not None:
        w_Omega, w_Ap = ncb
        pr(f"\n  New non-coboundary cocycle w' in T':")
        nz_edges = [(allowed_1p[i], w_Ap[i,0]) for i in range(len(allowed_1p)) if abs(w_Ap[i,0]) > 1e-8]
        involves_new = [(e, c) for e, c in nz_edges if e == new_arc]
        not_new = [(e, c) for e, c in nz_edges if e != new_arc]
        if involves_new:
            pr(f"    Coefficient on NEW arc {new_arc}: {involves_new[0][1]:.6f}")
        else:
            pr(f"    Coefficient on NEW arc {new_arc}: 0 (not involved!)")
        pr(f"    Other nonzero edges ({len(not_new)}):")
        for e, c in sorted(not_new, key=lambda x: -abs(x[1])):
            in_bad = "bad" if (e[0] in bad_old and e[1] in bad_old) else "mixed" if (e[0] in bad_old or e[1] in bad_old) else "good"
            pr(f"      {e}: {c:+.6f}  [{in_bad}]")

    # ---------- PART 1b: RESTRICTION TO SUBTOURNAMENTS ----------
    pr(f"\n  === PART 1b: Restriction to T\\v ===")

    for v in range(n):
        remaining = [u for u in range(n) if u != v]
        B, vmap = subtournament(A, n, remaining)
        m = len(remaining)
        _, Z1_sub, B1_sub, H1_sub, omega_sub, _ = compute_cocycle_spaces(B, m, 1)
        dim_Zs = Z1_sub.shape[1] if Z1_sub.ndim == 2 else 0
        dim_Bs = B1_sub.shape[1] if B1_sub.ndim == 2 else 0
        status = "BAD" if v in bad_old else "GOOD"

        B2, _ = subtournament(A2, n, remaining)
        _, Z1_sub2, B1_sub2, H1_sub2, _, _ = compute_cocycle_spaces(B2, m, 1)
        dim_Zs2 = Z1_sub2.shape[1] if Z1_sub2.ndim == 2 else 0
        dim_Bs2 = B1_sub2.shape[1] if B1_sub2.ndim == 2 else 0

        # T\v same or different? (only differs if flipped arc involves v)
        same = "SAME" if (v not in old_arc) else "DIFF"

        pr(f"    T\\{v} [{status}] ({same}): Z^1={dim_Zs} B^1={dim_Bs} H^1={H1_sub}"
           f"  |  T'\\{v}: Z^1={dim_Zs2} B^1={dim_Bs2} H^1={H1_sub2}")

    # ---------- PART 2: ALLOWED PATH CHANGES ----------
    pr(f"\n  === PART 2: Allowed path changes ===")

    cc_old = compute_full_chain_complex(A, n, max_dim=4)
    cc_new = compute_full_chain_complex(A2, n, max_dim=4)

    for p in range(5):
        old_set = set(cc_old['allowed'][p])
        new_set = set(cc_new['allowed'][p])
        gained = sorted(new_set - old_set)
        lost = sorted(old_set - new_set)
        if gained or lost:
            pr(f"    A_{p}: +{len(gained)} -{len(lost)} (was {len(old_set)}, now {len(new_set)})")
            if p <= 2:
                for path in gained[:6]:
                    pr(f"      + {path}")
                if len(gained) > 6: pr(f"      ... ({len(gained)} total gained)")
                for path in lost[:6]:
                    pr(f"      - {path}")
                if len(lost) > 6: pr(f"      ... ({len(lost)} total lost)")
        else:
            pr(f"    A_{p}: unchanged ({len(old_set)})")

    # ---------- PART 3: OMEGA DIMENSION CHANGES ----------
    pr(f"\n  === PART 3: Omega and chain complex dimensions ===")

    dims_old, _ = compute_chain_dims(A, n, max_dim=4)
    dims_new, _ = compute_chain_dims(A2, n, max_dim=4)

    pr(f"    k | dim_A  dim_O  dim_Z  dim_B  beta | dim_A' dim_O' dim_Z' dim_B' beta'| dO   dZ   dB   db")
    pr(f"   ---+--------------------------------------+---------------------------------------+-----------------")
    for k in range(min(5, n)):
        d = dims_old[k]
        dp = dims_new[k]
        dO = dp['dim_O'] - d['dim_O']
        dZ = dp['dim_Z'] - d['dim_Z']
        dB = dp['dim_B'] - d['dim_B']
        db = dp['beta'] - d['beta']
        pr(f"    {k} | {d['dim_A']:>5} {d['dim_O']:>5} {d['dim_Z']:>5} {d['dim_B']:>5} {d['beta']:>5}"
           f" | {dp['dim_A']:>5} {dp['dim_O']:>5} {dp['dim_Z']:>5} {dp['dim_B']:>5} {dp['beta']:>5}"
           f" | {dO:>+4} {dZ:>+4} {dB:>+4} {db:>+4}")

    # ---------- PART 4: BOUNDARY ANALYSIS ----------
    pr(f"\n  === PART 4: Boundary analysis (homological) ===")

    # In T' (after flip), check each 3-cycle
    for cyc in cycles_new:
        a, b, c = cyc
        # Find the directed orientation
        if A2[a][b] and A2[b][c] and A2[c][a]:
            directed = (a, b, c)  # a->b->c->a
        elif A2[b][a] and A2[a][c] and A2[c][b]:
            directed = (b, a, c)  # b->a->c->b
        elif A2[a][c] and A2[c][b] and A2[b][a]:
            directed = (a, c, b)  # a->c->b->a
        else:
            pr(f"    WARNING: could not find orientation for {cyc}")
            continue

        # The 3-cycle as a 1-chain: sum of its 3 edges
        cycle_edges = [(directed[0], directed[1]),
                       (directed[1], directed[2]),
                       (directed[2], directed[0])]

        # Express in A_1(T')
        idx_1p = {e: i for i, e in enumerate(allowed_1p)}
        cycle_vec = np.zeros((len(allowed_1p), 1))
        valid = True
        for e in cycle_edges:
            if e in idx_1p:
                cycle_vec[idx_1p[e], 0] = 1.0
            else:
                valid = False

        if not valid:
            continue

        # Check if in Omega_1
        if omega_1p.shape[1] > 0:
            cycle_in_O = omega_1p.T @ cycle_vec
            recon = omega_1p @ cycle_in_O
            residual = np.linalg.norm(cycle_vec - recon)
            in_omega = residual < 1e-8
        else:
            in_omega = False

        # Check if boundary
        omega_2p = cc_new['omega'][2]
        dim_O2p = omega_2p.shape[1] if omega_2p.ndim == 2 else 0
        if in_omega and dim_O2p > 0:
            bd2 = build_full_boundary_matrix(cc_new['allowed'][2], cc_new['allowed'][1])
            im_bd2 = bd2 @ omega_2p
            combined = np.column_stack([im_bd2, cycle_vec])
            rank_im = np.linalg.matrix_rank(im_bd2, tol=1e-8)
            rank_comb = np.linalg.matrix_rank(combined, tol=1e-8)
            is_bdy = (rank_comb == rank_im)
        elif in_omega:
            is_bdy = False
        else:
            is_bdy = None

        is_new = cyc in new_cycles
        tag = "NEW" if is_new else "old"
        pr(f"    3-cycle {directed[0]}->{directed[1]}->{directed[2]}->{directed[0]} [{tag}]:"
           f"  in_Omega_1={in_omega}, is_boundary={is_bdy}")

    # Now check the ORIGINAL 3-cycle in the ORIGINAL tournament
    pr(f"\n    In original T (beta_1=0):")
    for cyc in cycles_old:
        a, b, c = cyc
        if A[a][b] and A[b][c] and A[c][a]:
            directed = (a, b, c)
        elif A[b][a] and A[a][c] and A[c][b]:
            directed = (b, a, c)
        elif A[a][c] and A[c][b] and A[b][a]:
            directed = (a, c, b)
        else:
            continue

        cycle_edges = [(directed[0], directed[1]),
                       (directed[1], directed[2]),
                       (directed[2], directed[0])]
        idx_1 = {e: i for i, e in enumerate(allowed_1)}
        cycle_vec_orig = np.zeros((len(allowed_1), 1))
        for e in cycle_edges:
            if e in idx_1:
                cycle_vec_orig[idx_1[e], 0] = 1.0

        if omega_1.shape[1] > 0:
            cyc_in_O = omega_1.T @ cycle_vec_orig
            recon = omega_1 @ cyc_in_O
            residual = np.linalg.norm(cycle_vec_orig - recon)
            in_omega = residual < 1e-8
        else:
            in_omega = False

        omega_2_old = cc_old['omega'][2]
        dim_O2_old = omega_2_old.shape[1] if omega_2_old.ndim == 2 else 0
        if in_omega and dim_O2_old > 0:
            bd2_old = build_full_boundary_matrix(cc_old['allowed'][2], cc_old['allowed'][1])
            im_bd2_old = bd2_old @ omega_2_old
            combined = np.column_stack([im_bd2_old, cycle_vec_orig])
            rank_im_old = np.linalg.matrix_rank(im_bd2_old, tol=1e-8)
            rank_comb_old = np.linalg.matrix_rank(combined, tol=1e-8)
            is_bdy_old = (rank_comb_old == rank_im_old)
            pr(f"    3-cycle {directed[0]}->{directed[1]}->{directed[2]}->{directed[0]}: "
               f"in_Omega_1={in_omega}, is_boundary={is_bdy_old}")

            # Find the 2-chain that kills it
            if is_bdy_old:
                # Solve: im_bd2_old @ c = cycle_vec_orig
                # First project cycle_vec_orig into A_1 coords and find coefficients
                coeffs, _, _, _ = np.linalg.lstsq(im_bd2_old, cycle_vec_orig, rcond=None)
                # Express back in Omega_2 coords
                chain_2 = omega_2_old @ coeffs
                nz_2paths = [(cc_old['allowed'][2][i], chain_2[i,0])
                             for i in range(len(cc_old['allowed'][2]))
                             if abs(chain_2[i,0]) > 1e-8]
                pr(f"      Killing 2-chain ({len(nz_2paths)} paths):")
                for path, coeff in sorted(nz_2paths, key=lambda x: -abs(x[1])):
                    pr(f"        {path}: {coeff:+.6f}")
        elif in_omega:
            pr(f"    3-cycle {directed[0]}->{directed[1]}->{directed[2]}->{directed[0]}: "
               f"in_Omega_1={in_omega}, is_boundary=False (no Omega_2)")

    # Check which 2-paths in the killing chain are LOST after the flip
    pr(f"\n    2-paths lost by flip:")
    old_2 = set(cc_old['allowed'][2])
    new_2 = set(cc_new['allowed'][2])
    lost_2 = sorted(old_2 - new_2)
    gained_2 = sorted(new_2 - old_2)
    pr(f"      Lost: {lost_2[:8]}")
    if len(lost_2) > 8: pr(f"        ... ({len(lost_2)} total)")
    pr(f"      Gained: {gained_2[:8]}")
    if len(gained_2) > 8: pr(f"        ... ({len(gained_2)} total)")

    # ---------- PART 5: THE GLUING QUESTION ----------
    pr(f"\n  === PART 5: The gluing question ===")

    # The flip removes arc old_arc and adds new_arc.
    # This changes which paths are allowed and which Omega constraints exist.

    # Key question: what SPECIFIC constraint changes?
    # Non-allowed faces determine Omega. Compare.

    for p in [1, 2, 3]:
        old_paths = cc_old['allowed'][p]
        new_paths = cc_new['allowed'][p]
        old_faces_set = set(cc_old['allowed'][p-1])
        new_faces_set = set(cc_new['allowed'][p-1])

        # Count non-allowed faces for old vs new
        na_old = set()
        for path in old_paths:
            for sign, face in boundary_coeffs(path):
                if len(set(face)) == len(face) and face not in old_faces_set:
                    na_old.add(face)

        na_new = set()
        for path in new_paths:
            for sign, face in boundary_coeffs(path):
                if len(set(face)) == len(face) and face not in new_faces_set:
                    na_new.add(face)

        pr(f"    p={p}: non-allowed (p-1)-faces constraining Omega_{p}:")
        pr(f"      Before: {len(na_old)}, After: {len(na_new)}, Delta: {len(na_new)-len(na_old):+d}")
        lost_na = na_old - na_new
        gained_na = na_new - na_old
        if lost_na:
            pr(f"      Lost constraints: {sorted(lost_na)[:5]}")
        if gained_na:
            pr(f"      New constraints:  {sorted(gained_na)[:5]}")


# ========== AGGREGATE SUMMARY ==========
sep("AGGREGATE SUMMARY ACROSS ALL FLIP PAIRS")

# Collect dimension changes across ALL pairs
dO_counts = defaultdict(lambda: defaultdict(int))  # k -> delta -> count
dZ_counts = defaultdict(lambda: defaultdict(int))
dB_counts = defaultdict(lambda: defaultdict(int))
db_counts = defaultdict(lambda: defaultdict(int))

# Also track: does the original 3-cycle become a non-boundary after flip?
original_cycle_status = {'boundary_before': 0, 'nonboundary_after': 0, 'boundary_after': 0}

pr("\nProcessing all flip pairs for aggregate statistics...")
for A, A2, old_arc, new_arc in flip_pairs:
    dims_old, _ = compute_chain_dims(A, n, max_dim=3)
    dims_new, _ = compute_chain_dims(A2, n, max_dim=3)
    for k in range(4):
        dO = dims_new[k]['dim_O'] - dims_old[k]['dim_O']
        dZ = dims_new[k]['dim_Z'] - dims_old[k]['dim_Z']
        dB = dims_new[k]['dim_B'] - dims_old[k]['dim_B']
        db = dims_new[k]['beta'] - dims_old[k]['beta']
        dO_counts[k][dO] += 1
        dZ_counts[k][dZ] += 1
        dB_counts[k][dB] += 1
        db_counts[k][db] += 1

pr(f"\n{len(flip_pairs)} total flip pairs analyzed")
pr(f"\nDimension change distributions:")
for k in range(4):
    pr(f"\n  k={k}:")
    pr(f"    dOmega: {dict(sorted(dO_counts[k].items()))}")
    pr(f"    dZ:     {dict(sorted(dZ_counts[k].items()))}")
    pr(f"    dB:     {dict(sorted(dB_counts[k].items()))}")
    pr(f"    dbeta:  {dict(sorted(db_counts[k].items()))}")

# Key insight analysis
pr(f"\n{'='*78}")
pr(f"  KEY MECHANISM ANALYSIS")
pr(f"{'='*78}")

pr(f"""
The algebraic mechanism of the flip obstruction:

1. COCYCLE DIMENSION: When we flip one arc, Z^1 changes by dZ and B^1 by dB.
   The net H^1 = dim(Z^1) - dim(B^1) goes from 0 to 1.

2. What determines the mechanism:
   - If dZ > dB (Z^1 grows faster): new cocycle CREATED by the flip
   - If dZ = dB but dB < 0 (both shrink, B more): old coboundary DESTROYED
   - If dZ > 0 and dB = 0: pure cocycle creation
   - If dZ = 0 and dB < 0: pure coboundary destruction

Observed patterns:""")

# Classify the mechanism for each pair
mechanisms = defaultdict(int)
for A, A2, old_arc, new_arc in flip_pairs:
    _, Z1_old, B1_old, _, _, _ = compute_cocycle_spaces(A, n, 1)
    _, Z1_new, B1_new, _, _, _ = compute_cocycle_spaces(A2, n, 1)
    dZ = (Z1_new.shape[1] if Z1_new.ndim == 2 else 0) - (Z1_old.shape[1] if Z1_old.ndim == 2 else 0)
    dB = (B1_new.shape[1] if B1_new.ndim == 2 else 0) - (B1_old.shape[1] if B1_old.ndim == 2 else 0)
    if dZ > dB:
        if dB == 0:
            mech = "pure cocycle creation (dZ>0, dB=0)"
        elif dB > 0:
            mech = f"cocycle outpaces coboundary (dZ={dZ:+d}, dB={dB:+d})"
        else:
            mech = f"mixed: cocycle up, coboundary down (dZ={dZ:+d}, dB={dB:+d})"
    elif dZ == dB:
        mech = f"impossible: dZ=dB={dZ} but H^1 changed"
    else:
        if dZ == 0:
            mech = f"pure coboundary destruction (dZ=0, dB={dB:+d})"
        else:
            mech = f"coboundary drops faster (dZ={dZ:+d}, dB={dB:+d})"
    mechanisms[mech] += 1

for mech, count in sorted(mechanisms.items(), key=lambda x: -x[1]):
    pr(f"  {mech}: {count} cases ({100*count/len(flip_pairs):.1f}%)")

# Check: is the new non-boundary cocycle always supported on the new arc?
pr(f"\n  Does the new non-coboundary cocycle always use the NEW arc?")
uses_new_arc = 0
not_uses = 0
for A, A2, old_arc, new_arc in flip_pairs[:50]:  # sample
    _, Z1p, B1p, H1p, omega_1p, _ = compute_cocycle_spaces(A2, n, 1)
    allowed_1p = enumerate_allowed_paths(A2, n, 1)
    ncb = find_non_coboundary(Z1p, B1p, omega_1p, allowed_1p)
    if ncb is not None:
        _, w_Ap = ncb
        idx_1p = {e: i for i, e in enumerate(allowed_1p)}
        if new_arc in idx_1p:
            coeff = abs(w_Ap[idx_1p[new_arc], 0])
            if coeff > 1e-8:
                uses_new_arc += 1
            else:
                not_uses += 1
        else:
            not_uses += 1
    else:
        not_uses += 1

pr(f"    Uses new arc: {uses_new_arc}/{uses_new_arc + not_uses}")
pr(f"    Does NOT use new arc: {not_uses}/{uses_new_arc + not_uses}")
pr(f"    (Note: cocycle basis is not unique; this tests one representative)")

pr(f"\nDone.")
