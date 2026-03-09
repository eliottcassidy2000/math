#!/usr/bin/env python3
"""beta2_vertex_analysis.py - Analyze which vertices give injective i_*

For the inductive proof of beta_2=0, we need: for every tournament T,
there exists a vertex v such that i_*: H_1(T\v) -> H_1(T) is injective.

KEY FINDINGS FROM beta2_relative_induction.py:
  - n=5: 100% of tournaments have a good vertex
  - n=6: 100% of tournaments have a good vertex
  - n=7,8: 100% (sampled)
  - ker(i_*) = 3 at n=5, = 6 at n=6 when nonzero
  - Bad vertices: 840/5120 at n=5, 35328/196608 at n=6

This script characterizes good vs bad vertices.

Author: kind-pasteur-2026-03-08-S43
"""
import sys, os, random, time
import numpy as np
from collections import Counter
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

_saved = sys.stdout
sys.stdout = open(os.devnull, 'w', encoding='utf-8')
from path_homology_v2 import (
    enumerate_allowed_paths, build_full_boundary_matrix,
    compute_omega_basis
)
sys.stdout = _saved

random.seed(42)


def all_tournaments_gen(n):
    edges = [(i, j) for i in range(n) for j in range(i + 1, n)]
    m = len(edges)
    for mask in range(1 << m):
        A = [[0] * n for _ in range(n)]
        for idx, (i, j) in enumerate(edges):
            if (mask >> idx) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
        yield A


def get_induced(A, n, vertices):
    vlist = sorted(vertices)
    m = len(vlist)
    B = [[0] * m for _ in range(m)]
    for i in range(m):
        for j in range(m):
            B[i][j] = A[vlist[i]][vlist[j]]
    return B, vlist


def check_vertex(A, n, v):
    """Check i_* injectivity and return detailed info."""
    d_out = sum(A[v])
    d_in = n - 1 - d_out

    others = [x for x in range(n) if x != v]
    B_sub, vlist = get_induced(A, n, others)
    np2 = len(others)

    # b1(T\v)
    p1tv = enumerate_allowed_paths(B_sub, np2, 1)
    p0tv = [(i,) for i in range(np2)]
    p2tv = enumerate_allowed_paths(B_sub, np2, 2)

    if not p1tv:
        return {'d_out': d_out, 'b1_Tv': 0, 'b1_T': 0, 'ker': 0, 'good': True}

    o1tv = compute_omega_basis(B_sub, np2, 1, p1tv, p0tv)
    do1tv = o1tv.shape[1] if o1tv.ndim == 2 else 0
    if do1tv == 0:
        return {'d_out': d_out, 'b1_Tv': 0, 'b1_T': 0, 'ker': 0, 'good': True}

    D1tv = build_full_boundary_matrix([tuple(p) for p in p1tv], p0tv)
    D1tv_o = D1tv @ o1tv
    sv = np.linalg.svd(D1tv_o, compute_uv=False)
    rk1tv = int(sum(s > 1e-8 for s in sv))

    if p2tv:
        o2tv = compute_omega_basis(B_sub, np2, 2, p2tv, p1tv)
        do2tv = o2tv.shape[1] if o2tv.ndim == 2 else 0
        if do2tv > 0:
            D2tv = build_full_boundary_matrix([tuple(p) for p in p2tv],
                                              [tuple(p) for p in p1tv])
            D2tv_o = D2tv @ o2tv
            sv2 = np.linalg.svd(D2tv_o, compute_uv=False)
            rk2tv = int(sum(s > 1e-8 for s in sv2))
        else:
            rk2tv = 0
    else:
        rk2tv = 0

    b1tv = do1tv - rk1tv - rk2tv
    if b1tv == 0:
        return {'d_out': d_out, 'b1_Tv': 0, 'b1_T': 0, 'ker': 0, 'good': True}

    # Z_1(T\v) cycles
    _, _, Vt = np.linalg.svd(D1tv_o, full_matrices=True)
    z1tv = o1tv @ Vt[rk1tv:].T

    # Embed to T
    p1T = enumerate_allowed_paths(A, n, 1)
    p0T = [(i,) for i in range(n)]
    p2T = enumerate_allowed_paths(A, n, 2)

    p1T_idx = {tuple(p): i for i, p in enumerate(p1T)}
    emb = np.zeros((len(p1T), len(p1tv)))
    for j, ptv in enumerate(p1tv):
        pt = tuple(vlist[x] for x in ptv)
        if pt in p1T_idx:
            emb[p1T_idx[pt], j] = 1

    z1_in_T = emb @ z1tv

    # im(d_2|Omega_2(T))
    o2T = compute_omega_basis(A, n, 2, p2T, p1T)
    do2T = o2T.shape[1] if o2T.ndim == 2 else 0
    if do2T > 0:
        D2T = build_full_boundary_matrix([tuple(p) for p in p2T],
                                         [tuple(p) for p in p1T])
        imd2T = D2T @ o2T
    else:
        imd2T = np.zeros((len(p1T), 0))

    # b1(T)
    o1T = compute_omega_basis(A, n, 1, p1T, p0T)
    do1T = o1T.shape[1] if o1T.ndim == 2 else 0
    if do1T > 0:
        D1T = build_full_boundary_matrix([tuple(p) for p in p1T], p0T)
        D1T_o = D1T @ o1T
        sv1T = np.linalg.svd(D1T_o, compute_uv=False)
        rk1T = int(sum(s > 1e-8 for s in sv1T))
        rk2T = np.linalg.matrix_rank(imd2T, tol=1e-8) if imd2T.shape[1] > 0 else 0
        b1T = do1T - rk1T - rk2T
    else:
        b1T = 0

    # ker(i_*)
    ker = 0
    for k in range(z1_in_T.shape[1]):
        z = z1_in_T[:, k]
        if np.max(np.abs(z)) < 1e-12:
            ker += 1
            continue
        if imd2T.shape[1] > 0:
            c, _, _, _ = np.linalg.lstsq(imd2T, z, rcond=None)
            err = np.max(np.abs(imd2T @ c - z))
            if err < 1e-6:
                ker += 1

    return {
        'd_out': d_out, 'b1_Tv': b1tv, 'b1_T': b1T,
        'ker': ker, 'good': (ker == 0),
        'delta_b1': b1tv - b1T,
    }


# ============================================================
# Part 1: n=5 characterization
# ============================================================
print("=" * 70)
print("VERTEX CHARACTERIZATION at n=5")
print("=" * 70)

n = 5
good_data = []
bad_data = []

for A in all_tournaments_gen(n):
    scores = tuple(sorted([sum(row) for row in A]))
    for v in range(n):
        r = check_vertex(A, n, v)
        r['scores'] = scores
        if r['good']:
            good_data.append(r)
        else:
            bad_data.append(r)

print(f"Good: {len(good_data)}, Bad: {len(bad_data)}")

# Degree distribution
print(f"\nBy out-degree:")
for d in range(n):
    g = sum(1 for r in good_data if r['d_out'] == d)
    b = sum(1 for r in bad_data if r['d_out'] == d)
    t = g + b
    print(f"  d_out={d}: good={g}/{t} ({100*g/t:.0f}%)" if t > 0 else f"  d_out={d}: 0")

# b1 patterns
print(f"\nFor BAD vertices:")
print(f"  b1(T\\v) values: {Counter(r['b1_Tv'] for r in bad_data)}")
print(f"  b1(T) values: {Counter(r['b1_T'] for r in bad_data)}")
print(f"  delta_b1 = b1(T\\v)-b1(T): {Counter(r['delta_b1'] for r in bad_data)}")
print(f"  ker(i_*) values: {Counter(r['ker'] for r in bad_data)}")

# Score patterns for bad
print(f"\nScore sequences with bad vertices:")
sc_bad = Counter(r['scores'] for r in bad_data)
for sc, cnt in sc_bad.most_common():
    total_this_sc = sum(1 for r in good_data + bad_data if r['scores'] == sc)
    print(f"  {sc}: {cnt}/{total_this_sc} bad")

# Key question: b1(T\v) <= b1(T) iff good?
print(f"\nb1(T\\v) <= b1(T) vs goodness:")
le_good = sum(1 for r in good_data if r.get('delta_b1', 0) <= 0)
le_bad = sum(1 for r in bad_data if r.get('delta_b1', 0) <= 0)
gt_good = sum(1 for r in good_data if r.get('delta_b1', 0) > 0)
gt_bad = sum(1 for r in bad_data if r.get('delta_b1', 0) > 0)
print(f"  b1(T\\v) <= b1(T): good={le_good}, bad={le_bad}")
print(f"  b1(T\\v) > b1(T):  good={gt_good}, bad={gt_bad}")

# Is good <==> b1(T\v) <= b1(T)?
print(f"\n  => b1(T\\v) <= b1(T) implies goodness? {'YES' if le_bad == 0 else 'NO'}")
print(f"  => goodness implies b1(T\\v) <= b1(T)? {'YES' if gt_good == 0 else 'NO'}")


# ============================================================
# Part 2: Does EVERY tournament have a vertex with b1(T\v) <= b1(T)?
# ============================================================
print(f"\n{'=' * 70}")
print("DOES EVERY TOURNAMENT HAVE v WITH b1(T\\v) <= b1(T)?")
print("=" * 70)

for n in [5, 6]:
    t0 = time.time()
    gen = list(all_tournaments_gen(n))
    has_monotone_v = 0
    no_monotone_v = 0

    for A in gen:
        found = False
        for v in range(n):
            r = check_vertex(A, n, v)
            if r.get('delta_b1', 0) <= 0:
                found = True
                break
        if found:
            has_monotone_v += 1
        else:
            no_monotone_v += 1

    elapsed = time.time() - t0
    print(f"n={n} ({elapsed:.0f}s): has_monotone={has_monotone_v}, "
          f"no_monotone={no_monotone_v}/{len(gen)}")


# ============================================================
# Part 3: Compute ker(i_*) = dim, not just binary
# And check: ker(i_*) = C(d_out, 2) or similar?
# ============================================================
print(f"\n{'=' * 70}")
print("ker(i_*) FORMULA at n=5")
print("=" * 70)

n = 5
for d_out in range(n):
    subset = [r for r in bad_data if r['d_out'] == d_out]
    if not subset:
        continue
    kers = [r['ker'] for r in subset]
    b1tvs = [r['b1_Tv'] for r in subset]
    print(f"d_out={d_out}: ker(i_*)={Counter(kers)}, b1_Tv={Counter(b1tvs)}")

# At n=5, ker(i_*) = 3 always when bad.
# C(d_out-1, 2) at d_out=2: C(1,2) = 0. No.
# (n-1)(n-2)/2 = 6 at n=5. ker = 3 = 6/2.
# Actually let me check: which cycles die?
print(f"\nDetailed analysis of ONE bad (T,v) pair:")
for A in all_tournaments_gen(n):
    scores = tuple(sorted([sum(row) for row in A]))
    for v in range(n):
        r = check_vertex(A, n, v)
        if not r['good']:
            d_out = r['d_out']
            d_in = n - 1 - d_out
            P = [a for a in range(n) if a != v and A[v][a] == 1]
            Q = [b for b in range(n) if b != v and A[b][v] == 1]
            print(f"T scores={scores}, v={v}, d_out={d_out}")
            print(f"P (out-nbrs of v) = {P}")
            print(f"Q (in-nbrs of v) = {Q}")
            print(f"b1(T\\v) = {r['b1_Tv']}, b1(T) = {r['b1_T']}")
            print(f"ker(i_*) = {r['ker']}")
            print(f"delta_b1 = {r['delta_b1']}")

            # The cycles that die are 1-cycles in T\v that become
            # boundaries in T. These must be "filled" by 2-paths
            # involving v in T.

            # The key: removing v from T "opens up" some 1-cycles
            # that were boundaries (filled by 2-chains through v).
            # When we put v back, those 2-chains become available again.
            break
    else:
        continue
    break


print("\n\nDone.")
