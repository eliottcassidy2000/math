#!/usr/bin/env python3
"""
beta2_arcflip_invariance.py - Verify beta_2 arc-flip invariance at n=7,8

KEY DISCOVERY: Under any arc flip u->v to v->u in a tournament,
  delta(dim Z_2) = delta(rk d_3)  EXACTLY.

This means beta_2 = dim(Z_2) - rk(d_3) is an ARC-FLIP INVARIANT.
Since beta_2(transitive tournament) = 0, this proves beta_2 = 0 for ALL tournaments.

This script verifies the invariance at n=7 (sampling) and n=8 (sampling).
Also computes the exact delta values to find a formula.

Author: kind-pasteur-2026-03-08-S41
"""
import sys, time, os
import numpy as np
import random
from collections import Counter, defaultdict
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

_saved = sys.stdout
sys.stdout = open(os.devnull, 'w', encoding='utf-8')
from path_homology_v2 import (
    enumerate_allowed_paths, compute_omega_basis,
    build_full_boundary_matrix
)
sys.stdout = _saved

random.seed(42)

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

def compute_exactness_data(A, n):
    """Compute dim(Omega_2), dim(Z_2), rk(d_3), beta_2."""
    paths = {}
    omega = {}
    for p in range(5):
        paths[p] = enumerate_allowed_paths(A, n, p)
        if p == 0:
            omega[p] = np.eye(n)
        elif len(paths[p]) > 0 and len(paths[p-1]) > 0:
            omega[p] = compute_omega_basis(A, n, p, paths[p], paths[p-1])
        else:
            omega[p] = np.zeros((max(1, len(paths[p])), 0))

    dim_om2 = omega[2].shape[1] if omega[2].ndim == 2 else 0
    dim_om3 = omega[3].shape[1] if omega[3].ndim == 2 else 0

    if dim_om2 == 0:
        return {'dim_om2': 0, 'dim_om3': 0, 'dim_z2': 0, 'rk_d3': 0, 'beta2': 0,
                'rk_d2': 0, 'dim_om1': 0}

    bd2 = build_full_boundary_matrix(paths[2], paths[1])
    bd2_om = bd2 @ omega[2]
    sv2 = np.linalg.svd(bd2_om, compute_uv=False)
    rk_d2 = int(np.sum(np.abs(sv2) > 1e-8))
    dim_z2 = dim_om2 - rk_d2

    if dim_om3 > 0:
        bd3 = build_full_boundary_matrix(paths[3], paths[2])
        bd3_om = bd3 @ omega[3]
        im3_coords, _, _, _ = np.linalg.lstsq(omega[2], bd3_om, rcond=None)
        rk_d3 = np.linalg.matrix_rank(im3_coords, tol=1e-8)
    else:
        rk_d3 = 0

    beta2 = dim_z2 - rk_d3

    return {
        'dim_om2': dim_om2, 'dim_om3': dim_om3,
        'dim_z2': dim_z2, 'rk_d3': rk_d3, 'beta2': beta2,
        'rk_d2': rk_d2, 'dim_om1': len(paths[1]),
    }


def flip_arc(bits, i, j, n):
    idx = 0
    for a in range(n):
        for b in range(a+1, n):
            if a == i and b == j:
                return bits ^ (1 << idx)
            idx += 1
    return bits


def count_common_out(A, n, u, v):
    """Count vertices w with both u->w and v->w (w != u,v)."""
    return sum(1 for w in range(n) if w != u and w != v and A[u][w] and A[v][w])

def count_common_in(A, n, u, v):
    """Count vertices w with both w->u and w->v (w != u,v)."""
    return sum(1 for w in range(n) if w != u and w != v and A[w][u] and A[w][v])


print("=" * 70)
print("ARC-FLIP INVARIANCE OF beta_2")
print("beta_2 = dim(Z_2) - rk(d_3) is constant under arc flips!")
print("=" * 70)

for n in [7, 8]:
    print(f"\n--- n={n} sampled ---")
    n_arcs = n*(n-1)//2
    n_samples = 500 if n <= 7 else 100

    t0 = time.time()
    mismatches = 0
    total_flips = 0
    delta_dist = Counter()

    # Also track: delta values by local topology
    # (common_out, common_in, score_src, score_tgt) -> deltas
    local_deltas = defaultdict(list)

    for trial in range(n_samples):
        if trial % 100 == 0 and trial > 0:
            dt = time.time() - t0
            print(f"  ... {trial}/{n_samples} ({dt:.0f}s), mismatches={mismatches}")

        bits = random.randint(0, (1 << n_arcs) - 1)
        A = build_adj(n, bits)
        d = compute_exactness_data(A, n)
        scores = [sum(row) for row in A]

        # Pick a random arc to flip (instead of all arcs, to save time)
        for _ in range(5):  # 5 random flips per tournament
            i = random.randint(0, n-2)
            j = random.randint(i+1, n-1)

            bits2 = flip_arc(bits, i, j, n)
            A2 = build_adj(n, bits2)
            d2 = compute_exactness_data(A2, n)

            delta_z2 = d2['dim_z2'] - d['dim_z2']
            delta_im3 = d2['rk_d3'] - d['rk_d3']
            delta_om2 = d2['dim_om2'] - d['dim_om2']
            delta_om3 = d2['dim_om3'] - d['dim_om3']
            delta_rk2 = d2['rk_d2'] - d['rk_d2']

            delta_dist[(delta_z2, delta_im3)] += 1
            total_flips += 1

            if delta_z2 != delta_im3:
                mismatches += 1
                if mismatches <= 3:
                    print(f"    MISMATCH! n={n}, bits={bits}, flip ({i},{j})")
                    print(f"      T: Z2={d['dim_z2']}, im3={d['rk_d3']}, beta2={d['beta2']}")
                    print(f"      T': Z2={d2['dim_z2']}, im3={d2['rk_d3']}, beta2={d2['beta2']}")

            # Local topology
            if A[i][j] == 1:
                src, tgt = i, j
            else:
                src, tgt = j, i
            cout = count_common_out(A, n, src, tgt)
            cin = count_common_in(A, n, src, tgt)
            local_deltas[(scores[src], scores[tgt], cout, cin)].append(
                (delta_om2, delta_om3, delta_z2, delta_im3, delta_rk2))

    dt = time.time() - t0
    print(f"\n  Processed {n_samples} tournaments x 5 flips = {total_flips} flips in {dt:.0f}s")
    print(f"  Mismatches (dZ2 != dim3): {mismatches}")

    print(f"\n  Delta distribution:")
    for (dz, di), count in sorted(delta_dist.items()):
        match = "OK" if dz == di else "FAIL"
        pct = 100*count/total_flips
        print(f"    (dZ2={dz:+d}, dim3={di:+d}): {count} ({pct:.1f}%) {match}")

    # Check if delta_om2 determines delta_z2
    print(f"\n  Local topology analysis (top entries by frequency):")
    top_keys = sorted(local_deltas.keys(), key=lambda k: -len(local_deltas[k]))[:15]
    for key in top_keys:
        entries = local_deltas[key]
        ds, dt_score, cout, cin = key
        dz2_vals = [e[2] for e in entries]
        dz2_counter = Counter(dz2_vals)
        dom2_vals = [e[0] for e in entries]
        dom2_counter = Counter(dom2_vals)
        print(f"    scores=({ds},{dt_score}), c_out={cout}, c_in={cin}: "
              f"N={len(entries)}, dOm2={dict(dom2_counter)}, dZ2={dict(dz2_counter)}")

# Also verify: can we compute delta_dim(Omega_2) analytically?
# Under flip u->v to v->u:
# A TT triple (a,b,c) requires a->b, b->c, a->c.
# Gained TT: triples that are TT in T' but not in T
# Lost TT: triples that are TT in T but not in T'
print(f"\n{'='*70}")
print("ANALYTIC DELTA(Omega_2) CHECK at n=5")
print("=" * 70)

n = 5
n_arcs = n*(n-1)//2
total = 1 << n_arcs
errors = 0

for bits in range(min(total, 200)):
    A = build_adj(n, bits)
    scores = [sum(row) for row in A]

    for i in range(n):
        for j in range(i+1, n):
            bits2 = flip_arc(bits, i, j, n)
            A2 = build_adj(n, bits2)

            # Count TT triples in T and T'
            tt_T = sum(1 for a in range(n) for b in range(n) for c in range(n)
                       if a!=b and b!=c and a!=c and A[a][b] and A[b][c] and A[a][c])
            tt_T2 = sum(1 for a in range(n) for b in range(n) for c in range(n)
                        if a!=b and b!=c and a!=c and A2[a][b] and A2[b][c] and A2[a][c])

            # Analytic formula for delta_TT under flip u->v to v->u
            if A[i][j] == 1:
                u, v = i, j
            else:
                u, v = j, i
            # Lost: TT triples containing arc u->v
            # (a,u,v): a->u, u->v, a->v. Lost if a->v.
            # (u,b,v): u->b, b->v, u->v. Lost always (u->v used).
            #   Wait, (u,b,v) uses arc u->b and b->v and shortcut u->v
            # (u,v,c): u->v, v->c, u->c. Lost if u->c.
            # Also: (a,b,u) where b=v is not possible since u->v is the flipped arc
            # Actually need to be more careful...
            # A TT triple (a,b,c) uses arcs a->b, b->c, a->c.
            # It involves the arc u->v if (a,b)=(u,v) or (b,c)=(u,v) or (a,c)=(u,v).

            # Case 1: a=u, b=v. Triple (u,v,c). Requires u->v (flipped), v->c, u->c.
            lost_case1 = sum(1 for c in range(n) if c!=u and c!=v and A[v][c] and A[u][c])

            # Case 2: b=u, c=v. Triple (a,u,v). Requires a->u, u->v (flipped), a->v.
            lost_case2 = sum(1 for a in range(n) if a!=u and a!=v and A[a][u] and A[a][v])

            # Case 3: a=u, c=v. Triple (u,b,v). Requires u->b, b->v, u->v (shortcut, flipped).
            lost_case3 = sum(1 for b in range(n) if b!=u and b!=v and A[u][b] and A[b][v])

            total_lost = lost_case1 + lost_case2 + lost_case3

            # Gained: TT triples that use arc v->u (new arc)
            # Case 1: a=v, b=u. Triple (v,u,c). Requires v->u (new), u->c, v->c.
            gained_case1 = sum(1 for c in range(n) if c!=u and c!=v and A[u][c] and A[v][c])

            # Case 2: b=v, c=u. Triple (a,v,u). Requires a->v, v->u (new), a->u.
            gained_case2 = sum(1 for a in range(n) if a!=u and a!=v and A[a][v] and A[a][u])

            # Case 3: a=v, c=u. Triple (v,b,u). Requires v->b, b->u, v->u (shortcut, new).
            gained_case3 = sum(1 for b in range(n) if b!=u and b!=v and A[v][b] and A[b][u])

            total_gained = gained_case1 + gained_case2 + gained_case3

            predicted_delta = total_gained - total_lost
            actual_delta = tt_T2 - tt_T

            if predicted_delta != actual_delta:
                errors += 1
                if errors <= 3:
                    print(f"  ERROR at bits={bits}, flip u={u}->v={v}: "
                          f"predicted={predicted_delta}, actual={actual_delta}")

print(f"  Errors in delta(Omega_2) formula: {errors}/10000")
if errors == 0:
    print("  Formula VERIFIED: delta(Omega_2) = gained - lost TT triples involving flipped arc")

# Now check: is there a closed-form for gained/lost?
# lost_case1 = |{c: v->c AND u->c}| = common_out(u,v) in T
# lost_case2 = |{a: a->u AND a->v}| = common_in(u,v) in T
# lost_case3 = |{b: u->b AND b->v}| = directed paths u->b->v in T
# gained_case1 = |{c: u->c AND v->c}| = common_out(u,v) in T = lost_case1!
# gained_case2 = |{a: a->v AND a->u}| = common_in(u,v) in T = lost_case2!
# gained_case3 = |{b: v->b AND b->u}| = directed paths v->b->u in T
print("\n  Symmetry observation:")
print("  lost_case1 = gained_case1 (common out-neighbors)")
print("  lost_case2 = gained_case2 (common in-neighbors)")
print("  delta(Omega_2) = gained_case3 - lost_case3")
print("                 = |{b: v->b->u}| - |{b: u->b->v}|")
print("  = (# 2-paths v->b->u) - (# 2-paths u->b->v)")

print("\nDone.")
