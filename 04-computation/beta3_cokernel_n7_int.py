"""
beta3_cokernel_n7_int.py — Check integer structure of cokernel at n=7

At n=6, 3*cokernel is integer. Does the same hold at n=7, or is there a
different multiplier? Also check Type A vs Type B structure at n=7.

Author: kind-pasteur-S47 (2026-03-09)
"""
import sys
import numpy as np
from fractions import Fraction
from collections import defaultdict
sys.path.insert(0, '.')
sys.stdout.reconfigure(line_buffering=True)

from tournament_utils import (
    random_tournament, enumerate_all_allowed, compute_omega_basis_numpy,
    boundary_faces, full_chain_complex
)


def get_cokernel_raw(A, n):
    """Get raw cokernel generator (unnormalized). Returns (vec, ap3, data) or None."""
    data = full_chain_complex(A, n, max_p=5)
    if data['bettis'].get(3, 0) != 1:
        return None

    ap = data['ap']
    omega_bases = data['omega_bases']
    omega_dims = data['omega_dims']

    bd3 = np.zeros((len(ap[2]), len(ap[3])))
    idx2 = {path: i for i, path in enumerate(ap[2])}
    for j, path in enumerate(ap[3]):
        for sign, face in boundary_faces(path):
            if face in idx2:
                bd3[idx2[face], j] += sign

    O3 = omega_bases[3]
    d3_om = bd3 @ O3
    U3, S3, Vt3 = np.linalg.svd(d3_om, full_matrices=True)
    rank_d3 = int(sum(s > 1e-8 for s in S3))
    ker_d3_basis = Vt3[rank_d3:]

    if omega_dims.get(4, 0) > 0:
        bd4 = np.zeros((len(ap[3]), len(ap.get(4, []))))
        idx3 = {path: i for i, path in enumerate(ap[3])}
        for j, path in enumerate(ap.get(4, [])):
            for sign, face in boundary_faces(path):
                if face in idx3:
                    bd4[idx3[face], j] += sign
        O4 = omega_bases[4]
        d4_om = bd4 @ O4
        O3pinv = np.linalg.pinv(O3)
        d4_omega3 = O3pinv @ d4_om
        im_d4_in_ker = ker_d3_basis @ d4_omega3
        U4, S4, Vt4 = np.linalg.svd(im_d4_in_ker, full_matrices=True)
        rank_d4_in_ker = int(sum(s > 1e-8 for s in S4))
        coker_basis = U4[:, rank_d4_in_ker:]
    else:
        coker_basis = np.eye(ker_d3_basis.shape[0])

    coker_omega3 = coker_basis.T @ ker_d3_basis
    coker_A3 = (coker_omega3 @ O3.T).flatten()
    return coker_A3, data


def find_integer_multiplier(vec, max_mult=100):
    """Find smallest integer k such that k*vec is approximately integer."""
    nonzero = vec[np.abs(vec) > 1e-8]
    if len(nonzero) == 0:
        return 1

    # Normalize by max |coeff|
    norm = np.max(np.abs(nonzero))
    normed = nonzero / norm

    for k in range(1, max_mult + 1):
        scaled = k * normed
        if np.allclose(scaled, np.round(scaled), atol=0.01):
            return k
    return None


def main():
    print("=" * 70)
    print("BETA_3 COKERNEL INTEGER STRUCTURE AT n=7")
    print("=" * 70)

    n = 7
    rng = np.random.RandomState(42)
    found = 0
    multipliers = []
    omega4_dims = []

    for trial in range(3000):
        A = random_tournament(n, rng)
        result = get_cokernel_raw(A, n)
        if result is None:
            continue

        coker_A3, data = result
        found += 1

        # Find integer multiplier
        k = find_integer_multiplier(coker_A3)
        multipliers.append(k)
        omega4_dims.append(data['omega_dims'].get(4, 0))

        if found <= 10:
            norm = np.max(np.abs(coker_A3))
            normed = coker_A3 / norm if norm > 1e-10 else coker_A3

            # Check various multipliers
            for m in [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 15, 20]:
                scaled = m * normed
                rounded = np.round(scaled)
                residual = np.max(np.abs(scaled - rounded))
                if residual < 0.01:
                    vals = sorted(set(int(round(x)) for x in scaled if abs(x) > 0.01))
                    print(f"  #{found}: mult={m} works! omega4={data['omega_dims'].get(4,0)}, "
                          f"ker_d3={data['kers'][3]}, values={vals[:10]}")
                    break
            else:
                print(f"  #{found}: no small integer multiplier found (omega4={data['omega_dims'].get(4,0)}, "
                      f"ker_d3={data['kers'][3]})")
                # Show the nonzero ratios
                nonzero = normed[np.abs(normed) > 1e-8]
                ratios = sorted(set(round(x, 4) for x in nonzero))
                print(f"    Ratios: {ratios[:15]}")

        if found >= 50:
            break

        if (trial + 1) % 500 == 0:
            print(f"  {trial+1} checked, {found} with beta_3=1", flush=True)

    print(f"\n  Total found: {found} beta_3=1 from {trial+1} random n=7 tournaments")
    print(f"\n  Integer multiplier distribution:")
    from collections import Counter
    mult_dist = Counter(multipliers)
    for m, c in sorted(mult_dist.items()):
        print(f"    k={m}: {c} tournaments ({100*c/found:.1f}%)")

    print(f"\n  Omega_4 dimension distribution:")
    o4_dist = Counter(omega4_dims)
    for d, c in sorted(o4_dist.items()):
        print(f"    dim(Omega_4)={d}: {c} tournaments")

    print("\nDONE.")


if __name__ == '__main__':
    main()
