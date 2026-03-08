#!/usr/bin/env python3
"""
beta2_h2rel_characterize.py - Characterize when H2(T,T\v) > 0.

From exhaustive data:
- H2(T,T\v) > 0 iff beta1(T\v) = 1 (at n=5,6)
- H2(T,T\v) is always dim 1 when positive

QUESTION 1: Is H2(T,T\v) > 0 <=> beta1(T\v) = 1? (for interior v)
QUESTION 2: Why is dim H2(T,T\v) always exactly 1?
QUESTION 3: What is the relative 2-cycle explicitly?

For beta_2=0 proof: if H2(T,T\v) = 0 for some v, then LES gives
H2(T) = 0 trivially (no need for delta-injectivity at that v).
So we only need delta-injectivity when H2(T,T\v) > 0, which requires
beta1(T\v) = 1.

Author: kind-pasteur-2026-03-08-S42
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
            if bits & (1 << idx): A[i][j] = 1
            else: A[j][i] = 1
            idx += 1
    return A


def dim_om(om):
    return om.shape[1] if om.ndim == 2 and om.shape[0] > 0 else 0


def compute_h2_rel_dim(A, n, v):
    """Compute dim H2(T, T\v) using the quotient chain complex."""
    others = [i for i in range(n) if i != v]
    n1 = n - 1
    A_sub = [[A[others[i]][others[j]] for j in range(n1)] for i in range(n1)]

    ap0 = [(i,) for i in range(n)]
    ap1 = enumerate_allowed_paths(A, n, 1)
    ap2 = enumerate_allowed_paths(A, n, 2)
    ap1_list = [tuple(p) for p in ap1]

    om2 = compute_omega_basis(A, n, 2, ap2, ap1) if ap2 else np.zeros((0,0))
    d2 = dim_om(om2)
    if d2 == 0:
        return 0

    # v-arc projection of boundary
    v_arc_mask = np.array([1.0 if v in p else 0.0 for p in ap1_list])
    bd2 = build_full_boundary_matrix(ap2, ap1)
    om1 = compute_omega_basis(A, n, 1, ap1, ap0)
    bd2_om = bd2 @ om2
    v_proj = np.diag(v_arc_mask) @ bd2_om
    rk_v = np.linalg.matrix_rank(v_proj, tol=1e-8)

    # Relative 2-cycles = kernel of v-arc projection in Omega_2
    z2_rel_dim = d2 - rk_v

    # Subtract Omega_2(T\v) dimension
    ap2_sub = enumerate_allowed_paths(A_sub, n1, 2)
    ap1_sub = enumerate_allowed_paths(A_sub, n1, 1)
    ap0_sub = [(i,) for i in range(n1)]
    om2_sub = compute_omega_basis(A_sub, n1, 2, ap2_sub, ap1_sub) if ap2_sub else np.zeros((0,0))
    d2_sub = dim_om(om2_sub)

    # Relative dim = z2_rel_dim - d2_sub (approx, need to subtract intersection)
    # Actually need to compute properly using quotient

    # Better: compute via LES
    # H2(T,T\v) = ker(delta) where delta: H2(T,T\v) -> H1(T\v)
    # From LES: dim H2(T,T\v) = dim H2(T) + dim(ker(H1(T\v)->H1(T)))
    #                          + dim(coker(H2(T\v)->H2(T)))
    # Actually just use H2(T,T\v) = Z2_rel / B2_rel
    # where Z2_rel = {x in Om2 : del(x) in Om1(T\v)} / Om2(T\v)
    # and B2_rel = image of del_3 in the quotient

    # Simple formula from LES:
    # ... -> H2(T\v) -> H2(T) -> H2(T,T\v) -> H1(T\v) -> H1(T) -> ...
    # dim H2(T,T\v) = dim H2(T) - dim(im H2(T\v)->H2(T))
    #               + dim(ker delta)
    # Since we assume beta2=0 by induction:
    # dim H2(T,T\v) = 0 + dim(ker(H1(T\v)->H1(T)))
    # Wait that's using beta2(T)=0 which we're trying to prove...

    # Let me just compute it directly
    betti_T = path_betti_numbers(A, n, max_dim=2)
    betti_sub = path_betti_numbers(A_sub, n1, max_dim=2)

    return betti_T[2], betti_sub[1]


# ============================================================
# Q1: Correlation between H2(T,T\v) > 0 and beta1(T\v)
# ============================================================
print("=" * 70)
print("Q1: WHEN IS H2(T,T\\v) > 0?")
print("=" * 70)

for n in [4, 5, 6]:
    m = n*(n-1)//2
    total = 1 << m if n <= 6 else 1000

    joint = Counter()  # (beta2_T, beta1_sub, is_interior) -> count

    t0 = time.time()
    for bits in (range(total) if n <= 6 else []):
        A = build_adj(n, bits)
        scores = [sum(A[i][j] for j in range(n) if j != i) for i in range(n)]
        betti_T = path_betti_numbers(A, n, max_dim=2)
        b2 = betti_T[2]

        for v in range(n):
            dv = scores[v]
            is_int = 1 <= dv <= n-2
            others = [i for i in range(n) if i != v]
            A_sub = [[A[others[i]][others[j]] for j in range(n-1)] for i in range(n-1)]
            b_sub = path_betti_numbers(A_sub, n-1, max_dim=1)
            b1_sub = b_sub[1]

            joint[(b2, b1_sub, is_int)] += 1

    elapsed = time.time() - t0
    print(f"\nn={n} ({total} tournaments, {elapsed:.0f}s):")
    print(f"  (beta2_T, beta1(T\\v), is_interior) -> count:")
    for key in sorted(joint.keys()):
        b2, b1s, is_int = key
        vtype = "INT" if is_int else "EXT"
        print(f"    beta2={b2}, beta1(T\\v)={b1s}, {vtype}: {joint[key]}")


# ============================================================
# Q2: LES dimension analysis
# ============================================================
print(f"\n{'='*70}")
print("Q2: LES DIMENSION ANALYSIS")
print("=" * 70)

# From the LES:
# H2(T\v) -> H2(T) -> H2(T,T\v) -> H1(T\v) ->^j H1(T) -> H1(T,T\v) -> H0(T\v) -> H0(T)
#
# beta0(T\v) = beta0(T) = 1 (connected)
# H1(T,T\v): from LES, dim = dim(ker(H0(T\v)->H0(T))) + dim(coker j) + dim H1(T)
#
# Key: dim H2(T,T\v) = dim(ker delta)
#   where delta: H2(T,T\v) -> H1(T\v)
#   By exactness: ker(delta) = im(H2(T) -> H2(T,T\v))
#
# If beta2(T) = 0: dim H2(T,T\v) = dim(ker delta) + h2_rel_free
#   where h2_rel_free = dim(things killed by delta going to H1)
#
# The LES gives:
#   0 -> H2(T) -> H2(T,T\v) -> H1(T\v) -> H1(T) -> ...
#   (when beta2(T\v)=0)
#
# So: rank of alpha = dim H2(T)
#     rank of delta = dim H2(T,T\v) - dim H2(T)
#     kernel of j = im delta = dim H2(T,T\v) - dim H2(T)
#     rank of j = beta1(T\v) - (dim H2(T,T\v) - dim H2(T))

# Let's compute all the LES terms explicitly at n=5
n = 5
print(f"\nn={n}: Full LES analysis")
print(f"{'bits':>6} {'scores':<16} {'b2T':>3} {'b1T':>3} {'b1Tv':>4} {'h2rel':>5} {'rk_delta':>8} {'rk_j':>4}")

count = 0
for bits in range(1 << (n*(n-1)//2)):
    A = build_adj(n, bits)
    scores = tuple(sorted([sum(A[i][j] for j in range(n) if j != i) for i in range(n)]))
    betti_T = path_betti_numbers(A, n, max_dim=2)
    b2T = betti_T[2]
    b1T = betti_T[1]

    for v in range(n):
        dv = sum(A[v])
        if dv == 0 or dv == n-1:
            continue

        others = [i for i in range(n) if i != v]
        A_sub = [[A[others[i]][others[j]] for j in range(4)] for i in range(4)]
        b_sub = path_betti_numbers(A_sub, 4, max_dim=2)
        b1_sub = b_sub[1]
        b2_sub = b_sub[2]

        if b1_sub == 0:
            continue

        # Compute rank of j: H1(T\v) -> H1(T)
        # This is the map induced by inclusion T\v -> T
        ap1_T = enumerate_allowed_paths(A, n, 1)
        ap0_T = [(i,) for i in range(n)]
        om1_T = compute_omega_basis(A, n, 1, ap1_T, ap0_T)
        d1_T = dim_om(om1_T)

        ap1_sub = enumerate_allowed_paths(A_sub, 4, 1)
        ap0_sub = [(i,) for i in range(4)]
        om1_sub = compute_omega_basis(A_sub, 4, 1, ap1_sub, ap0_sub)
        d1_sub = dim_om(om1_sub)

        # Embed Z1(T\v) into Z1(T), then quotient by B1
        remap = {i: others[i] for i in range(4)}
        bd1_T = build_full_boundary_matrix(ap1_T, ap0_T)
        rk_d1_T = np.linalg.matrix_rank(bd1_T @ om1_T, tol=1e-8)
        z1_T_dim = d1_T - rk_d1_T

        bd1_sub = build_full_boundary_matrix(ap1_sub, ap0_sub)
        rk_d1_sub = np.linalg.matrix_rank(bd1_sub @ om1_sub, tol=1e-8)
        z1_sub_dim = d1_sub - rk_d1_sub

        # H1 = Z1/B1. Map j: H1(T\v) -> H1(T).
        # j maps a 1-cycle in T\v to the same 1-cycle viewed in T.
        # rk(j) = dim(im j).
        # By LES: rk(delta) + rk(j) = beta1(T\v)
        # and rk(alpha) + rk(delta) = dim H2(T,T\v)
        # and rk(alpha) = beta2(T) (= 0 if our conjecture holds)

        # So: dim H2(T,T\v) = beta2(T) + rk(delta) = 0 + (beta1_sub - rk(j))
        # => dim H2(T,T\v) = beta1_sub - rk(j)

        # Compute rk(j) = rank of map H1(T\v) -> H1(T)
        # If beta1(T) = 0, then H1(T) = 0, so rk(j) = 0
        # If beta1(T) = 1 and beta1(T\v) = 1, rk(j) in {0,1}

        rk_j = "?"
        if b1T == 0:
            rk_j = 0
            h2_rel_predicted = b1_sub - 0
        else:
            # Need to compute explicitly
            h2_rel_predicted = "?"

        if count < 20:
            print(f"{bits:>6} {str(scores):<16} {b2T:>3} {b1T:>3} {b1_sub:>4} {h2_rel_predicted:>5} {'':>8} {rk_j:>4}")
        count += 1

print(f"\n... ({count} total interior cases with beta1(T\\v)=1)")

# KEY FORMULA:
# If beta2(T) = 0 (conjecture), then:
#   dim H2(T,T\v) = beta1(T\v) - rk(j)
#   where j: H1(T\v) -> H1(T) is the inclusion map
#
# delta injective means: rk(delta) = dim H2(T,T\v)
# By LES: rk(delta) = beta1(T\v) - rk(j)
# So delta injective iff rk(delta) = beta1(T\v) - rk(j)
# which is ALWAYS true by LES!
#
# Wait - does this mean delta-injectivity is EQUIVALENT to beta2(T)=0
# (given beta2(T\v)=0)?

print(f"\n{'='*70}")
print("KEY INSIGHT: LES IMPLIES delta-injectivity <=> beta2(T)=0")
print("=" * 70)
print("""
From the LES with beta2(T\\v) = 0:
  0 -> H2(T) ->^alpha H2(T,T\\v) ->^delta H1(T\\v) ->^j H1(T) -> ...

Exactness gives:
  im(alpha) = ker(delta)
  rk(alpha) = beta2(T)
  rk(delta) = dim H2(T,T\\v) - beta2(T)

delta injective <=> ker(delta) = 0 <=> im(alpha) = 0 <=> beta2(T) = 0.

So proving delta-injectivity for interior v IS proving beta2(T)=0.
It's NOT a weaker statement — it's EQUIVALENT.

The question is: can we prove it WITHOUT assuming beta2=0?
""")

# ============================================================
# Q3: What makes the H1(T\v)->H1(T) map special?
# ============================================================
print(f"\n{'='*70}")
print("Q3: RANK OF j: H1(T\\v) -> H1(T)")
print("=" * 70)

n = 5
rk_j_dist = Counter()

for bits in range(1 << (n*(n-1)//2)):
    A = build_adj(n, bits)
    scores = [sum(A[i][j] for j in range(n) if j != i) for i in range(n)]
    betti_T = path_betti_numbers(A, n, max_dim=1)
    b1T = betti_T[1]

    for v in range(n):
        dv = scores[v]
        if dv == 0 or dv == n-1:
            continue

        others = [i for i in range(n) if i != v]
        A_sub = [[A[others[i]][others[j]] for j in range(4)] for i in range(4)]
        b_sub = path_betti_numbers(A_sub, 4, max_dim=1)
        b1_sub = b_sub[1]

        if b1_sub == 0:
            continue

        # Compute rk(j)
        if b1T == 0:
            rk_j_dist[(b1T, b1_sub, 0, dv)] += 1
        elif b1T == 1 and b1_sub == 1:
            # Need to check if the H1 generator of T\v maps nontrivially to H1(T)
            # The 3-cycle in T\v: is it homologous to zero in T?
            # A 3-cycle (a,b,c) in T\v is a boundary in T iff it's the boundary
            # of some 2-chain in Omega_2(T).

            # Actually, compute directly
            ap1_T = enumerate_allowed_paths(A, n, 1)
            ap0_T = [(i,) for i in range(n)]
            ap2_T = enumerate_allowed_paths(A, n, 2)
            om1_T = compute_omega_basis(A, n, 1, ap1_T, ap0_T)
            om2_T = compute_omega_basis(A, n, 2, ap2_T, ap1_T) if ap2_T else np.zeros((0,0))

            ap1_sub = enumerate_allowed_paths(A_sub, 4, 1)
            ap0_sub = [(i,) for i in range(4)]
            om1_sub = compute_omega_basis(A_sub, 4, 1, ap1_sub, ap0_sub)

            bd1_T = build_full_boundary_matrix(ap1_T, ap0_T) @ om1_T
            rk_d1 = np.linalg.matrix_rank(bd1_T, tol=1e-8)
            _, _, Vt1 = np.linalg.svd(bd1_T, full_matrices=True)
            z1_T = Vt1[rk_d1:]  # Z1(T) basis in Om1(T) coords

            # B1(T) = im(d2)
            d2_T = dim_om(om2_T)
            if d2_T > 0:
                bd2_T = build_full_boundary_matrix(ap2_T, ap1_T)
                b1_coords = np.linalg.lstsq(om1_T, bd2_T @ om2_T, rcond=None)[0]
                b1_in_z1 = z1_T @ b1_coords
                rk_b1 = np.linalg.matrix_rank(b1_in_z1, tol=1e-8)
            else:
                rk_b1 = 0

            # Embed H1 generator of T\v into T
            # Z1(T\v) generator
            bd1_sub = build_full_boundary_matrix(ap1_sub, ap0_sub) @ om1_sub
            rk_d1_sub = np.linalg.matrix_rank(bd1_sub, tol=1e-8)
            _, _, Vt1_sub = np.linalg.svd(bd1_sub, full_matrices=True)
            z1_sub = Vt1_sub[rk_d1_sub:]

            # Map z1_sub generator to T coords
            ap1_T_list = [tuple(p) for p in ap1_T]
            remap = {i: others[i] for i in range(4)}
            d1_sub = dim_om(om1_sub)
            embed1 = np.zeros((len(ap1_T_list), d1_sub))
            for j in range(d1_sub):
                for k, p in enumerate(ap1_sub):
                    pT = tuple(remap[x] for x in p)
                    if pT in ap1_T_list:
                        embed1[ap1_T_list.index(pT), j] = om1_sub[k, j]
            psi = np.linalg.lstsq(om1_T, embed1, rcond=None)[0]

            # z1_sub_gen in T coords
            z1_sub_gen_T = psi @ z1_sub.T  # (d1_T, z1_sub_dim)

            # Project to Z1(T)/B1(T) = H1(T)
            z1_sub_in_z1T = z1_T @ z1_sub_gen_T  # (z1_T_dim, z1_sub_dim)
            # Check if it's in B1(T)
            if d2_T > 0 and rk_b1 > 0:
                combined_check = np.hstack([b1_in_z1.T, z1_sub_in_z1T])
                rk_comb = np.linalg.matrix_rank(combined_check, tol=1e-8)
                rk_j_val = rk_comb - rk_b1
            else:
                rk_j_val = np.linalg.matrix_rank(z1_sub_in_z1T, tol=1e-8)

            rk_j_dist[(b1T, b1_sub, rk_j_val, dv)] += 1

print("(beta1_T, beta1_sub, rk_j, d+) -> count:")
for key in sorted(rk_j_dist.keys()):
    b1T, b1s, rj, dv = key
    print(f"  beta1(T)={b1T}, beta1(T\\v)={b1s}, rk(j)={rj}, d+={dv}: {rk_j_dist[key]}")

# If rk(j) = 0 when beta1(T)=0, then H2(T,T\v) = 1 (= beta1_sub).
# If rk(j) = 1 when beta1(T)=1 and beta1(T\v)=1, then H2(T,T\v) = 0.
# If rk(j) = 0 when beta1(T)=1 and beta1(T\v)=1, then H2(T,T\v) = 1.


print("\n\nDone.")
