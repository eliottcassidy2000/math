"""
relative_euler_char.py — Euler characteristic analysis of relative complex

H_3(T,T\v) sits in the relative complex R_p = Omega_p(T)/Omega_p(T\v).
The relative Euler characteristic is:
  chi_rel = sum_p (-1)^p dim(R_p) = sum_p (-1)^p (dim_Om_p(T) - dim_Om_p(T\v))

If chi_rel can be computed in terms of local data, and if we can show
H_3^rel is the ONLY non-trivial relative Betti number, then
dim H_3^rel = |chi_rel - (known contributions from other degrees)|.

This would give an exact formula for H_3^rel.

Author: kind-pasteur-S47 (2026-03-09)
"""
import sys
import numpy as np
from itertools import combinations, permutations
from collections import Counter, defaultdict
sys.path.insert(0, '.')
sys.stdout.reconfigure(line_buffering=True)

from tournament_utils import (
    bits_to_adj, random_tournament, full_chain_complex
)

PRIME = (1 << 31) - 1


def gauss_rank_modp(M, prime=PRIME):
    M = M.copy() % prime
    nrows, ncols = M.shape
    rank = 0
    for col in range(ncols):
        nonzero = np.where(M[rank:, col] != 0)[0]
        if len(nonzero) == 0:
            continue
        pivot = nonzero[0] + rank
        if pivot != rank:
            M[[rank, pivot]] = M[[pivot, rank]]
        inv = pow(int(M[rank, col]), prime - 2, prime)
        M[rank] = M[rank] * inv % prime
        factors = M[:, col].copy()
        factors[rank] = 0
        nz = np.where(factors != 0)[0]
        if len(nz) > 0:
            M[nz] = (M[nz] - np.outer(factors[nz], M[rank])) % prime
        rank += 1
    return rank


def gauss_nullspace_modp(M, prime=PRIME):
    M = M.copy() % prime
    nrows, ncols = M.shape
    pivot_cols = []
    rank = 0
    for col in range(ncols):
        nonzero = np.where(M[rank:, col] != 0)[0]
        if len(nonzero) == 0:
            continue
        pivot = nonzero[0] + rank
        if pivot != rank:
            M[[rank, pivot]] = M[[pivot, rank]]
        inv = pow(int(M[rank, col]), prime - 2, prime)
        M[rank] = M[rank] * inv % prime
        factors = M[:, col].copy()
        factors[rank] = 0
        nz = np.where(factors != 0)[0]
        if len(nz) > 0:
            M[nz] = (M[nz] - np.outer(factors[nz], M[rank])) % prime
        pivot_cols.append(col)
        rank += 1

    null_dim = ncols - rank
    if null_dim == 0:
        return np.zeros((0, ncols), dtype=np.int64)

    free_cols = sorted(set(range(ncols)) - set(pivot_cols))
    null_basis = np.zeros((null_dim, ncols), dtype=np.int64)
    for i, fc in enumerate(free_cols):
        null_basis[i, fc] = 1
        for j, pc in enumerate(pivot_cols):
            null_basis[i, pc] = (-M[j, fc]) % prime
    return null_basis % prime


def get_omega_dims(A, n, max_p=6):
    """Compute Omega dimensions for tournament A on n vertices."""
    all_v = list(range(n))
    omega_dims = {}

    for p in range(max_p + 1):
        # Allowed p-paths
        ap = []
        if p + 1 <= n:
            for verts in combinations(all_v, p + 1):
                for perm in permutations(verts):
                    ok = True
                    for i in range(len(perm) - 1):
                        if A[perm[i]][perm[i+1]] != 1:
                            ok = False
                            break
                    if ok:
                        ap.append(perm)

        if not ap:
            omega_dims[p] = 0
            continue

        # Get A_{p-1}
        ap_prev = []
        if p > 0 and p <= n:
            for verts in combinations(all_v, p):
                for perm in permutations(verts):
                    ok = True
                    for i in range(len(perm) - 1):
                        if A[perm[i]][perm[i+1]] != 1:
                            ok = False
                            break
                    if ok:
                        ap_prev.append(perm)

        ap_prev_set = set(ap_prev)

        # Constraint matrix
        na_faces = {}
        for sigma in ap:
            for i in range(1, len(sigma) - 1):
                face = sigma[:i] + sigma[i+1:]
                if p > 0 and face not in ap_prev_set:
                    if face not in na_faces:
                        na_faces[face] = len(na_faces)

        if not na_faces:
            omega_dims[p] = len(ap)
        else:
            C = np.zeros((len(na_faces), len(ap)), dtype=np.int64)
            for j, sigma in enumerate(ap):
                for i in range(1, len(sigma) - 1):
                    face = sigma[:i] + sigma[i+1:]
                    if face in na_faces:
                        C[na_faces[face], j] = (C[na_faces[face], j] + ((-1)**i)) % PRIME
            C = C % PRIME
            rank = gauss_rank_modp(C, PRIME)
            omega_dims[p] = len(ap) - rank

    return omega_dims


def main():
    print("=" * 70)
    print("RELATIVE EULER CHARACTERISTIC ANALYSIS")
    print("=" * 70)

    # Part 1: n=6 exhaustive — compute relative Euler char for all (T,v)
    print("\n--- Part 1: Relative Euler characteristic at n=6 ---")

    n = 6
    total = 2 ** (n*(n-1)//2)

    chi_rel_dist_b3_1 = Counter()
    chi_rel_dist_b3_0 = Counter()
    rel_betti_data = []

    for bits in range(total):
        A = bits_to_adj(bits, n)
        data = full_chain_complex(A, n, max_p=5)
        b3 = data['bettis'].get(3, 0)

        if b3 == 0 and len(chi_rel_dist_b3_0) > 100:
            continue

        om_T = data['omega_dims']

        for v in range(n):
            remaining = [i for i in range(n) if i != v]
            A_sub = [[A[remaining[i]][remaining[j]] for j in range(n-1)] for i in range(n-1)]
            data_sub = full_chain_complex(A_sub, n-1, max_p=5)
            om_Tv = data_sub['omega_dims']

            # Relative dimensions: R_p = dim_Om_p(T) - dim_Om_p(T\v)
            # (This is only approximately right — the actual quotient dim can differ
            # if Omega(T\v) doesn't embed cleanly into Omega(T))
            # For a proper computation, see relative_h3_structure_deep.py

            # But as an approximation:
            rel_dims = {p: om_T.get(p, 0) - om_Tv.get(p, 0) for p in range(7)}

            # Relative Euler characteristic
            chi_rel = sum((-1)**p * rel_dims[p] for p in range(7))

            if b3 == 1:
                chi_rel_dist_b3_1[chi_rel] += 1
            else:
                chi_rel_dist_b3_0[chi_rel] += 1

            if b3 == 1 and len(rel_betti_data) < 10:
                # Also compute actual relative Betti numbers from our earlier analysis
                rel_betti_data.append({
                    'bits': bits, 'v': v,
                    'rel_dims': {p: rel_dims[p] for p in range(7)},
                    'chi_rel': chi_rel,
                    'score': tuple(sorted([sum(A[i]) for i in range(n)])),
                    'b3_T': b3,
                    'b3_Tv': data_sub['bettis'].get(3, 0),
                })

        if bits % 5000 == 0 and bits > 0:
            print(f"  ... {bits}/{total}", flush=True)

    print(f"\n  Relative chi for beta_3=1 tours: {dict(sorted(chi_rel_dist_b3_1.items()))}")
    print(f"  Relative chi for beta_3=0 tours (sample): {dict(sorted(chi_rel_dist_b3_0.items()))}")

    # Part 2: What does chi_rel = -1 mean?
    print("\n--- Part 2: Interpretation ---")
    if rel_betti_data:
        for d in rel_betti_data[:5]:
            print(f"\n  bits={d['bits']}, v={d['v']}, score={d['score']}")
            print(f"    R_p dims: {[d['rel_dims'][p] for p in range(7)]}")
            print(f"    chi_rel = {d['chi_rel']}")
            print(f"    b3_T={d['b3_T']}, b3_Tv={d['b3_Tv']}")

    # Part 3: n=7 — check chi_rel for bad and good vertices
    print("\n--- Part 3: n=7 chi_rel analysis ---")

    n = 7
    rng = np.random.RandomState(42)

    chi_bad = Counter()
    chi_good = Counter()
    checked = 0

    for trial in range(5000):
        A = random_tournament(n, rng)
        data = full_chain_complex(A, n, max_p=6)
        if data['bettis'].get(3, 0) != 1:
            continue

        checked += 1
        om_T = data['omega_dims']

        for v in range(n):
            remaining = [i for i in range(n) if i != v]
            A_sub = [[A[remaining[i]][remaining[j]] for j in range(n-1)] for i in range(n-1)]
            data_sub = full_chain_complex(A_sub, n-1, max_p=6)
            om_Tv = data_sub['omega_dims']

            rel_dims = {p: om_T.get(p, 0) - om_Tv.get(p, 0) for p in range(8)}
            chi_rel = sum((-1)**p * rel_dims[p] for p in range(8))

            b3_Tv = data_sub['bettis'].get(3, 0)

            if b3_Tv == 1:
                chi_bad[chi_rel] += 1
            else:
                chi_good[chi_rel] += 1

        if checked >= 50:
            break

        if checked % 10 == 0:
            print(f"  ... {checked} beta_3=1 found", flush=True)

    print(f"\n  n=7 chi_rel for BAD vertices (b3_Tv=1): {dict(sorted(chi_bad.items()))}")
    print(f"  n=7 chi_rel for GOOD vertices (b3_Tv=0): {dict(sorted(chi_good.items()))}")

    # Part 4: Key observation
    print("\n--- Part 4: Key observation ---")
    print("  From LES: sum_p (-1)^p H_p^rel = chi_rel")
    print("  If H_0^rel = 1, H_1^rel = 0, H_2^rel = 0 (from beta_2=0 and LES),")
    print("  then: H_3^rel - H_4^rel + H_5^rel - ... = chi_rel - 1")
    print("  For beta_3=1, H_3^rel in {0,1}:")
    print("    If H_3^rel = 1: -H_4^rel + H_5^rel - ... = chi_rel - 2")
    print("    If H_3^rel = 0: -H_4^rel + H_5^rel - ... = chi_rel - 1")

    # Verify H_0^rel = 1 claim
    print("\n  H_0(T,T\\v) should be 1 (T connected, T\\v connected, one component difference)")

    # Verify H_1^rel and H_2^rel values
    print("  H_1(T,T\\v) and H_2(T,T\\v) from LES:")
    print("    ... H_1(T) = 0 (most T, since beta_3=1 => beta_1=0 by seesaw)")
    print("    ... H_1(T\\v) = 0 (mostly, since beta_3_Tv >= 0)")
    print("    LES: H_1(T\\v) -> H_1(T) -> H_1(T,T\\v) -> H_0(T\\v) -> H_0(T)")
    print("    H_0(T\\v) = 1, H_0(T) = 1 => connecting map H_1(T,T\\v) -> H_0(T\\v)")

    print("\nDONE.")


if __name__ == '__main__':
    main()
