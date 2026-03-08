#!/usr/bin/env python3
"""
beta2_relative_structure.py - Analyze the structure of the relative complex

For the pair (T, T\v), the relative chain complex is:
  Omega_p^rel = Omega_p(T) / Omega_p(T\{v})

where Omega_p(T\{v}) is identified with paths in Omega_p(T) not using v.

Key question: WHY is H_2^rel = 0 for tournaments?

Approach:
1. Compute dim(Omega_p^rel) for p=1,2,3
2. Compute rk(d_p^rel)
3. See if there's a dimensional reason (exact sequence becomes trivially exact)

For a vertex v with out-neighborhood N+(v) and in-neighborhood N-(v):
- 2-paths using v: (a,v,c) with a in N-, c in N+;
  plus (v,b,c) with b in N+, b->c; plus (a,b,v) with a->b, b in N-
- The relative complex "sees" only the v-dependent part

Author: kind-pasteur-2026-03-08-S41
"""
import sys, time, os, random
import numpy as np
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

_saved = sys.stdout
sys.stdout = open(os.devnull, 'w', encoding='utf-8')
from path_homology_v2 import (
    enumerate_allowed_paths, compute_omega_basis,
    build_full_boundary_matrix
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

def random_tournament(n):
    return random.getrandbits(n*(n-1)//2)

def compute_relative_dims(A, n, v):
    """Compute dimensions of relative chain complex at vertex v."""
    results = {}

    # Compute Omega_p and paths for p=0,1,2,3,4
    a = {}
    om = {}
    for p in range(5):
        a[p] = enumerate_allowed_paths(A, n, p)
        if p == 0:
            om[p] = np.eye(n)
        elif p >= 1 and len(a[p]) > 0 and len(a[p-1]) > 0:
            om[p] = compute_omega_basis(A, n, p, a[p], a[p-1])
        else:
            om[p] = np.zeros((max(1, len(a[p])), 0))

    for p in range(5):
        dim_om = om[p].shape[1] if om[p].ndim == 2 and om[p].shape[1] > 0 else 0

        if dim_om == 0:
            results[p] = {'dim_om': 0, 'dim_rel': 0, 'dim_nv': 0}
            continue

        # v-path indices
        v_idx = [i for i in range(len(a[p])) if v in a[p][i]]

        if not v_idx:
            results[p] = {'dim_om': dim_om, 'dim_rel': 0, 'dim_nv': dim_om}
            continue

        # Restriction to v-paths
        R_v = om[p][v_idx, :]
        sv = np.linalg.svd(R_v, compute_uv=False)
        rank_rv = int(np.sum(np.abs(sv) > 1e-8))

        dim_nv = dim_om - rank_rv  # dim(Omega_p ∩ no-v)
        dim_rel = rank_rv

        results[p] = {
            'dim_om': dim_om,
            'dim_rel': dim_rel,
            'dim_nv': dim_nv,
            'n_vpaths': len(v_idx),
            'n_paths': len(a[p])
        }

    # Compute boundary ranks in relative complex
    for p in [2, 3]:
        if results[p]['dim_rel'] == 0 or results[p-1].get('dim_rel', 0) == 0:
            results[p]['rk_bd_rel'] = 0
            continue

        dim_p = om[p].shape[1] if om[p].ndim == 2 else 0
        dim_pm1 = om[p-1].shape[1] if om[p-1].ndim == 2 else 0

        if dim_p == 0 or dim_pm1 == 0:
            results[p]['rk_bd_rel'] = 0
            continue

        # Boundary d_p: A_p -> A_{p-1}
        bd = build_full_boundary_matrix(a[p], a[p-1])
        bd_om = bd @ om[p]  # |A_{p-1}| x dim_p

        # Express in Omega_{p-1} coordinates
        im_p, _, _, _ = np.linalg.lstsq(om[p-1], bd_om, rcond=None)
        # im_p: dim_{p-1} x dim_p

        # V_{p-1} = Omega_{p-1} ∩ no-v
        v_idx_pm1 = [i for i in range(len(a[p-1])) if v in a[p-1][i]]
        R_v_pm1 = om[p-1][v_idx_pm1, :]
        sv_pm1 = np.linalg.svd(R_v_pm1, compute_uv=False)
        rk_pm1 = int(np.sum(np.abs(sv_pm1) > 1e-8))

        if rk_pm1 < dim_pm1:
            Uv, Sv, Vtv = np.linalg.svd(R_v_pm1, full_matrices=True)
            V_pm1 = Vtv[rk_pm1:].T  # basis of V_{p-1}

            # V_p = Omega_p ∩ no-v
            v_idx_p = [i for i in range(len(a[p])) if v in a[p][i]]
            R_v_p = om[p][v_idx_p, :]
            sv_p = np.linalg.svd(R_v_p, compute_uv=False)
            rk_p = int(np.sum(np.abs(sv_p) > 1e-8))

            if rk_p < dim_p:
                Uvp, Svp, Vtvp = np.linalg.svd(R_v_p, full_matrices=True)
                V_p = Vtvp[rk_p:].T
            else:
                V_p = np.zeros((dim_p, 0))

            # rk(d_p^rel) = rk of im_p modulo V_{p-1}, restricted to quotient Omega_p/V_p
            # = rk([im_p | V_{p-1}]) - dim(V_{p-1}) minus rk(im_p restricted to V_p modulo V_{p-1})

            # Simpler: compute rank of the induced map on quotients
            # Map: Omega_p/V_p -> Omega_{p-1}/V_{p-1}
            # Represented by im_p in coordinates, modulo V_{p-1} in target and V_p in source

            if V_p.shape[1] > 0:
                # Image of V_p under im_p (lives in V_{p-1} by d(no-v) = no-v)
                im_Vp = im_p @ V_p
                # Combine im_p, V_{p-1}, im_Vp
                all_cols = np.hstack([im_p, V_pm1])
                rk_all = np.linalg.matrix_rank(all_cols, tol=1e-8)

                Vp_cols = np.hstack([im_Vp, V_pm1])
                rk_Vp = np.linalg.matrix_rank(Vp_cols, tol=1e-8)

                results[p]['rk_bd_rel'] = rk_all - rk_Vp
            else:
                all_cols = np.hstack([im_p, V_pm1])
                rk_all = np.linalg.matrix_rank(all_cols, tol=1e-8)
                rk_V = np.linalg.matrix_rank(V_pm1, tol=1e-8)
                results[p]['rk_bd_rel'] = rk_all - rk_V
        else:
            # V_{p-1} = 0, so relative = absolute
            results[p]['rk_bd_rel'] = np.linalg.matrix_rank(im_p, tol=1e-8)

    return results


# ===== Exhaustive n=5 =====
print("=" * 70)
print("RELATIVE COMPLEX STRUCTURE")
print("=" * 70)

for n in [5]:
    print(f"\n--- n={n} exhaustive ---")
    n_arcs = n*(n-1)//2
    total = 1 << n_arcs

    # Collect statistics
    rel_dim_stats = {}  # (d_v, dim_rel_1, dim_rel_2, dim_rel_3) -> count

    for bits in range(total):
        A = build_adj(n, bits)

        for v in range(n):
            d_v = sum(A[v])  # out-degree
            d_v_in = n - 1 - d_v

            res = compute_relative_dims(A, n, v)

            key = (d_v,
                   res.get(1, {}).get('dim_rel', 0),
                   res.get(2, {}).get('dim_rel', 0),
                   res.get(3, {}).get('dim_rel', 0),
                   res.get(2, {}).get('rk_bd_rel', 0),
                   res.get(3, {}).get('rk_bd_rel', 0))
            rel_dim_stats[key] = rel_dim_stats.get(key, 0) + 1

    print(f"\n  (d_v, rel_1, rel_2, rel_3, rk_d2_rel, rk_d3_rel): count")
    for key in sorted(rel_dim_stats.keys()):
        d_v, r1, r2, r3, rk2, rk3 = key
        h2_rel = r2 - rk2 - rk3
        print(f"    d_v={d_v}: rel=({r1},{r2},{r3}), rk_d=(_, {rk2},{rk3}), "
              f"H_2^rel={h2_rel}, count={rel_dim_stats[key]}")

# ===== n=6 sample =====
print(f"\n--- n=6 sample (5000 tournaments) ---")
n = 6
random.seed(42)
rel_dim_stats_6 = {}

for trial in range(5000):
    bits = random_tournament(n)
    A = build_adj(n, bits)

    for v in range(n):
        d_v = sum(A[v])
        res = compute_relative_dims(A, n, v)

        key = (d_v,
               res.get(1, {}).get('dim_rel', 0),
               res.get(2, {}).get('dim_rel', 0),
               res.get(3, {}).get('dim_rel', 0),
               res.get(2, {}).get('rk_bd_rel', 0),
               res.get(3, {}).get('rk_bd_rel', 0))
        rel_dim_stats_6[key] = rel_dim_stats_6.get(key, 0) + 1

print(f"\n  (d_v, rel_1, rel_2, rel_3, rk_d2_rel, rk_d3_rel): count")
for key in sorted(rel_dim_stats_6.keys()):
    d_v, r1, r2, r3, rk2, rk3 = key
    h2_rel = r2 - rk2 - rk3
    print(f"    d_v={d_v}: rel=({r1},{r2},{r3}), rk_d=(_, {rk2},{rk3}), "
          f"H_2^rel={h2_rel}, count={rel_dim_stats_6[key]}")

# ===== Key question: is dim_rel_2 = rk_d2_rel + rk_d3_rel always? =====
print(f"\n{'='*70}")
print("CHECK: Is Omega_2^rel = ker(d_2^rel) filled by im(d_3^rel)?")
print("=" * 70)

for label, stats in [("n=5 exhaustive", rel_dim_stats), ("n=6 sample", rel_dim_stats_6)]:
    print(f"\n  {label}:")
    for key in sorted(stats.keys()):
        d_v, r1, r2, r3, rk2, rk3 = key
        ker2_rel = r2 - rk2  # ker(d_2^rel)
        h2_rel = ker2_rel - rk3
        if h2_rel != 0:
            print(f"    *** H_2^rel != 0: d_v={d_v}, rel=({r1},{r2},{r3}), "
                  f"ker={ker2_rel}, im={rk3}")
        else:
            # Check if ker(d_2^rel) = 0 (no cycles) or im(d_3^rel) fills it
            if ker2_rel == 0:
                status = "ker=0 (no cycles)"
            else:
                status = f"ker={ker2_rel}, im fills exactly"
            print(f"    d_v={d_v}: {status}, count={stats[key]}")

# ===== Formula investigation =====
print(f"\n{'='*70}")
print("FORMULA: dim(Omega_p^rel) as function of d_v?")
print("=" * 70)

for label, stats in [("n=5", rel_dim_stats), ("n=6", rel_dim_stats_6)]:
    print(f"\n  {label}:")
    by_dv = {}
    for key in stats:
        d_v = key[0]
        r1, r2, r3 = key[1], key[2], key[3]
        if d_v not in by_dv:
            by_dv[d_v] = set()
        by_dv[d_v].add((r1, r2, r3))

    for d_v in sorted(by_dv):
        vals = sorted(by_dv[d_v])
        if len(vals) == 1:
            print(f"    d_v={d_v}: CONSTANT rel dims = {vals[0]}")
        else:
            print(f"    d_v={d_v}: VARIES, {len(vals)} patterns: {vals[:5]}...")

print("\nDone.")
