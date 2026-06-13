#!/usr/bin/env python3
"""beta2_rk_d2_formula.py - Find formula for rk(d_2) in tournaments

KEY EQUATIONS:
- dim(Z_1) = (n-1)(n-2)/2  [constant for tournaments]
- b1(T) = dim(Z_1) - rk(d_2) = (n-1)(n-2)/2 - rk(d_2)
- dim(Omega_1) = n(n-1)/2 - (n-1) = (n-1)(n-2)/2  [arcs - (n-1)]
  Wait, actually: Omega_1 = ker(d_0) restricted to allowed 1-paths.
  For tournaments: every edge is an allowed 1-path.
  Omega_1 basis: need d_0 on 1-paths. d_0(i,j) = j - i.
  So Omega_1 = ker(d_0|A_1) has dim = #edges - rk(d_0).
  #edges = n(n-1)/2. rk(d_0) = n-1 (connected graph). So dim(Omega_1) = n(n-1)/2 - (n-1) = (n-1)(n-2)/2.

  And Z_1 = ker(d_1|Omega_1). Since Omega_1 is already defined as ker(d_0),
  actually I'm confused. Let me re-derive.

  Path homology chain complex:
    A_p = span of allowed p-paths
    Omega_p = {x in A_p : d_p(x) in Omega_{p-1}}
    d_p: Omega_p -> Omega_{p-1} is the restriction

  For p=1: A_1 = allowed 1-paths (= directed edges of T)
  Omega_0 = A_0 = R^n (vertices)
  d_1(i,j) = j - i. This is always in Omega_0 = A_0. So Omega_1 = A_1.
  dim(Omega_1) = n(n-1)/2 (number of directed edges in tournament).

  Wait, n(n-1)/2 edges for a tournament on n vertices? Each pair has exactly one edge, so yes, C(n,2) = n(n-1)/2.

  rk(d_1) = n - 1 (tournament is connected).
  Z_1 = ker(d_1|Omega_1): dim = n(n-1)/2 - (n-1) = (n-1)(n-2)/2. Correct.

  b1 = dim(Z_1) - rk(d_2) = (n-1)(n-2)/2 - rk(d_2).

Now: what determines rk(d_2)?

d_2: Omega_2 -> Omega_1 maps 2-chains to 1-chains.
rk(d_2) = dim(im(d_2)).

im(d_2) is a subspace of Z_1 (since d_1 o d_2 = 0 in the chain complex? Actually d_1 o d_2 = 0 only in Omega_*).

Wait: in path homology, d_{p-1} o d_p = 0 on Omega_p. So im(d_2) is in ker(d_1) = Z_1.

So: rk(d_2) <= dim(Z_1) = (n-1)(n-2)/2.
And: b1 = dim(Z_1) - rk(d_2) >= 0.

For tournaments: b1 is typically 0 or small.
b1 = 0 iff rk(d_2) = (n-1)(n-2)/2.

Is rk(d_2) always close to (n-1)(n-2)/2?

Author: kind-pasteur-2026-03-08-S43
"""
import sys, os, random, time
import numpy as np
from collections import Counter, defaultdict
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


def random_tournament(n):
    A = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            if random.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1
    return A


def get_induced(A, n, vertices):
    vlist = sorted(vertices)
    m = len(vlist)
    B = [[0] * m for _ in range(m)]
    for i in range(m):
        for j in range(m):
            B[i][j] = A[vlist[i]][vlist[j]]
    return B, vlist


def compute_dimensions(A, n):
    """Compute key dimensions of the chain complex."""
    paths1 = enumerate_allowed_paths(A, n, 1)
    paths0 = [(i,) for i in range(n)]
    paths2 = enumerate_allowed_paths(A, n, 2)
    paths3 = enumerate_allowed_paths(A, n, 3)

    dim_A1 = len(paths1)

    omega1 = compute_omega_basis(A, n, 1, paths1, paths0) if paths1 else np.zeros((0, 0))
    dim_O1 = omega1.shape[1] if omega1.ndim == 2 and omega1.shape[1] > 0 else 0

    if paths2 and paths1:
        omega2 = compute_omega_basis(A, n, 2, paths2, paths1)
        dim_O2 = omega2.shape[1] if omega2.ndim == 2 and omega2.shape[1] > 0 else 0
    else:
        omega2 = np.zeros((0, 0))
        dim_O2 = 0

    if paths3 and paths2:
        omega3 = compute_omega_basis(A, n, 3, paths3, paths2)
        dim_O3 = omega3.shape[1] if omega3.ndim == 2 and omega3.shape[1] > 0 else 0
    else:
        dim_O3 = 0

    # Ranks of boundary maps
    rk_d1 = 0
    if dim_O1 > 0:
        D1 = build_full_boundary_matrix([tuple(p) for p in paths1], paths0)
        D1_om = D1 @ omega1
        sv = np.linalg.svd(D1_om, compute_uv=False)
        rk_d1 = int(sum(s > 1e-8 for s in sv))

    rk_d2 = 0
    if dim_O2 > 0:
        D2 = build_full_boundary_matrix([tuple(p) for p in paths2],
                                        [tuple(p) for p in paths1])
        D2_om = D2 @ omega2
        sv2 = np.linalg.svd(D2_om, compute_uv=False)
        rk_d2 = int(sum(s > 1e-8 for s in sv2))

    rk_d3 = 0
    if dim_O3 > 0 and paths3:
        D3 = build_full_boundary_matrix([tuple(p) for p in paths3],
                                        [tuple(p) for p in paths2])
        D3_om = D3 @ omega3
        sv3 = np.linalg.svd(D3_om, compute_uv=False)
        rk_d3 = int(sum(s > 1e-8 for s in sv3))

    z1_dim = dim_O1 - rk_d1
    z2_dim = dim_O2 - rk_d2

    b0 = n - rk_d1  # actually beta_0 = dim(ker d_0|Omega_0) - ... but for connected graph = 1
    b1 = z1_dim - rk_d2
    b2 = z2_dim - rk_d3

    return {
        'dim_A1': dim_A1, 'dim_A2': len(paths2), 'dim_A3': len(paths3) if paths3 else 0,
        'dim_O1': dim_O1, 'dim_O2': dim_O2, 'dim_O3': dim_O3,
        'rk_d1': rk_d1, 'rk_d2': rk_d2, 'rk_d3': rk_d3,
        'z1_dim': z1_dim, 'z2_dim': z2_dim,
        'b1': b1, 'b2': b2,
    }


# ============================================================
# Part 1: rk(d_2) analysis
# ============================================================
print("=" * 70)
print("rk(d_2) FORMULA SEARCH")
print("=" * 70)

for n in [4, 5, 6]:
    gen = list(all_tournaments_gen(n))
    z1_formula = (n - 1) * (n - 2) // 2

    rk_d2_vals = []
    b1_vals = []
    dim_O2_vals = []

    for A in gen:
        d = compute_dimensions(A, n)
        rk_d2_vals.append(d['rk_d2'])
        b1_vals.append(d['b1'])
        dim_O2_vals.append(d['dim_O2'])

    print(f"\nn={n}: {len(gen)} tournaments")
    print(f"  dim(Z_1) = {z1_formula} (constant)")
    print(f"  rk(d_2): min={min(rk_d2_vals)}, max={max(rk_d2_vals)}, "
          f"dist={Counter(rk_d2_vals)}")
    print(f"  b1: min={min(b1_vals)}, max={max(b1_vals)}, "
          f"dist={Counter(b1_vals)}")
    print(f"  dim(Omega_2): min={min(dim_O2_vals)}, max={max(dim_O2_vals)}")

    # Formula candidates for rk(d_2)
    # rk(d_2) = z1_dim iff b1 = 0
    # rk(d_2) = z1_dim - b1 always (definition)
    # But what determines b1?

    # Check: is b1 = number of directed 3-cycles? Or related to cycle count?
    for A in gen[:1]:
        d = compute_dimensions(A, n)
        c3 = 0
        for i in range(n):
            for j in range(n):
                for k in range(n):
                    if i != j and j != k and i != k:
                        if A[i][j] == 1 and A[j][k] == 1 and A[k][i] == 1:
                            c3 += 1
        c3 //= 3  # each cycle counted 3 times
        print(f"  Example: b1={d['b1']}, c3={c3}")


# ============================================================
# Part 2: rk(d_2) vs rk(d_2|T\v) drop analysis
# ============================================================
print(f"\n{'=' * 70}")
print("RANK DROP: rk(d_2(T)) - rk(d_2(T\\v))")
print("=" * 70)

n = 5
drops_by_dout = defaultdict(list)
z1_n = (n - 1) * (n - 2) // 2
z1_nm1 = (n - 2) * (n - 3) // 2

for A in all_tournaments_gen(n):
    d_T = compute_dimensions(A, n)
    rk_T = d_T['rk_d2']

    for v in range(n):
        others = [x for x in range(n) if x != v]
        B_sub, vlist = get_induced(A, n, others)
        d_Tv = compute_dimensions(B_sub, len(others))
        rk_Tv = d_Tv['rk_d2']

        drop = rk_T - rk_Tv
        d_out = sum(A[v])
        drops_by_dout[d_out].append(drop)

print(f"n={n}: dim(Z_1(T))={z1_n}, dim(Z_1(T\\v))={z1_nm1}")
print(f"  Threshold for b1 monotonicity: drop <= {n-2}")
print(f"\nDrop distribution by d_out:")
for d in sorted(drops_by_dout.keys()):
    vals = drops_by_dout[d]
    dist = Counter(vals)
    total = len(vals)
    good = sum(v for k, v in dist.items() if k <= n - 2)
    print(f"  d_out={d}: {dict(sorted(dist.items()))}, good={good}/{total}")

# Theoretical lower bound for drop:
# rk(d_2(T)) counts the number of "linearly independent boundary equations"
# from Omega_2(T).
# rk(d_2(T\v)) counts the same for T\v.
# The difference is the "new" boundaries contributed by 2-chains involving v.
#
# Omega_2(T) includes:
# 1. All Omega_2(T\v) elements (not involving v) -> contribute rk(d_2(T\v))
# 2. New elements involving v -> contribute additional rank
#
# The additional rank from v-involving elements is at most:
# dim(Omega_2(T)) - dim(Omega_2(T\v))
# And at least: rk(d_2(T)) - rk(d_2(T\v)) = the drop.

print(f"\n{'=' * 70}")
print("dim(Omega_2) ANALYSIS")
print("=" * 70)

n = 5
for A in all_tournaments_gen(n):
    d_T = compute_dimensions(A, n)
    for v in range(n):
        others = [x for x in range(n) if x != v]
        B_sub, vlist = get_induced(A, n, others)
        d_Tv = compute_dimensions(B_sub, len(others))

        delta_O2 = d_T['dim_O2'] - d_Tv['dim_O2']
        drop = d_T['rk_d2'] - d_Tv['rk_d2']

        if drop > n - 2:
            d_out = sum(A[v])
            scores = sorted([sum(row) for row in A])
            print(f"  BAD: v={v}, d_out={d_out}, scores={scores}")
            print(f"    dim_O2(T)={d_T['dim_O2']}, dim_O2(T\\v)={d_Tv['dim_O2']}, "
                  f"delta_O2={delta_O2}")
            print(f"    rk_d2(T)={d_T['rk_d2']}, rk_d2(T\\v)={d_Tv['rk_d2']}, "
                  f"drop={drop}")
            print(f"    b1(T)={d_T['b1']}, b1(T\\v)={d_Tv['b1']}")
            break
    else:
        continue
    break


# ============================================================
# Part 3: b1 formula for tournaments
# ============================================================
print(f"\n{'=' * 70}")
print("b1 AS FUNCTION OF TOURNAMENT STRUCTURE")
print("=" * 70)

# b1 is related to the number of "independent" directed cycles
# that are not fillable by 2-chains.

# Hypothesis: b1 = number of directed 3-cycles - dim(Omega_2)?
# Or b1 = some function of the cycle counts.

n = 5
b1_c3_data = []
for A in all_tournaments_gen(n):
    d = compute_dimensions(A, n)

    # Count directed 3-cycles
    c3 = 0
    for i in range(n):
        for j in range(i + 1, n):
            for k in range(j + 1, n):
                fwd = A[i][j] + A[j][k] + A[k][i]
                rev = A[j][i] + A[k][j] + A[i][k]
                if fwd == 3:
                    c3 += 1
                if rev == 3:
                    c3 += 1

    # Count TT (transitive) triples
    tt = 0
    for i in range(n):
        for j in range(n):
            for k in range(n):
                if len({i, j, k}) == 3:
                    if A[i][j] == 1 and A[j][k] == 1 and A[i][k] == 1:
                        tt += 1
    tt //= 1  # each TT triple counted once as (i,j,k) with i->j->k, i->k

    b1_c3_data.append((d['b1'], c3, d['dim_O2'], d['rk_d2']))

# Check: is b1 a function of c3?
b1_by_c3 = defaultdict(set)
for b1, c3, dO2, rk in b1_c3_data:
    b1_by_c3[c3].add(b1)

print(f"n={n}: b1 as function of c3 (number of directed 3-cycles):")
for c3 in sorted(b1_by_c3.keys()):
    print(f"  c3={c3}: b1 in {sorted(b1_by_c3[c3])}")

# Check if b1 = 0 iff c3 >= some threshold
b1_zero = [(c3, dO2, rk) for b1, c3, dO2, rk in b1_c3_data if b1 == 0]
b1_pos = [(c3, dO2, rk) for b1, c3, dO2, rk in b1_c3_data if b1 > 0]
if b1_zero:
    print(f"\n  b1=0: c3 range [{min(c for c,_,_ in b1_zero)}, {max(c for c,_,_ in b1_zero)}]")
if b1_pos:
    print(f"  b1>0: c3 range [{min(c for c,_,_ in b1_pos)}, {max(c for c,_,_ in b1_pos)}]")


print("\n\nDone.")
