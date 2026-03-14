"""
euler_char_glmy.py -- kind-pasteur-2026-03-14-S69
Correct GLMY path homology Euler characteristic for tournaments.

Uses the correct boundary map with junk-face projection (from betti_alpha_v2.py).
Key question: what does chi = sum (-1)^k beta_k equal?

If chi = 1 always, path homology categorifies "1" (trivially).
If chi depends on the tournament, it's a new invariant.
If chi relates to H somehow, that's the Khovanov analogy!
"""

import numpy as np
from itertools import permutations
from collections import Counter, defaultdict
import sys

sys.stdout.reconfigure(encoding='utf-8')

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

def compute_H_dp(A, n):
    """Compute H(T) using dynamic programming (Held-Karp)."""
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for ms in range(2, n+1):
        for mask in range(1 << n):
            if bin(mask).count('1') != ms:
                continue
            for v in range(n):
                if not (mask & (1 << v)):
                    continue
                pm = mask ^ (1 << v)
                t = 0
                for u in range(n):
                    if (pm & (1 << u)) and A[u][v]:
                        t += dp.get((pm, u), 0)
                if t:
                    dp[(mask, v)] = t
    return sum(dp.get(((1 << n) - 1, v), 0) for v in range(n))

def compute_betti(A, n, max_d=None):
    """
    Compute GLMY path homology Betti numbers for tournament with adjacency A.
    Correct implementation with junk-face projection.
    (Adapted from betti_alpha_v2.py / betti_sigma_hierarchy.py)
    """
    if max_d is None:
        max_d = n - 1

    # Generate allowed paths at each degree
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

    actual_max = max(d for d in paths if paths[d]) if any(paths[d] for d in paths) else 0

    path_idx = {}
    for d in range(actual_max + 1):
        path_idx[d] = {p: i for i, p in enumerate(paths[d])}

    allowed_set = {d: set(paths[d]) for d in range(actual_max + 1)}

    omega_dims = [n]
    bd_ranks = [0]

    for d in range(1, min(max_d + 1, actual_max + 1)):
        n_d = len(paths[d])
        n_dm1 = len(paths[d-1])
        if n_d == 0:
            omega_dims.append(0)
            bd_ranks.append(0)
            continue

        junk_set = set()
        face_junk = []
        face_allowed = []

        for p in paths[d]:
            jf = []; af = []
            for fi in range(d + 1):
                face = p[:fi] + p[fi+1:]
                sign = 1 if fi % 2 == 0 else -1
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

        # C = constraint matrix (junk faces)
        C = np.zeros((n_junk, n_d), dtype=float)
        for j, jf in enumerate(face_junk):
            for face, sign in jf:
                C[junk_idx[face], j] += sign

        rank_c = int(np.linalg.matrix_rank(C)) if n_junk > 0 else 0
        omega_d = n_d - rank_c
        omega_dims.append(omega_d)

        # CB = combined constraint + boundary matrix
        CB = np.zeros((n_junk + n_dm1, n_d), dtype=float)
        CB[:n_junk, :] = C
        for j, af in enumerate(face_allowed):
            for face, sign in af:
                row = n_junk + path_idx[d-1][face]
                CB[row, j] += sign

        rank_cb = int(np.linalg.matrix_rank(CB))
        bd_rank = rank_cb - rank_c
        bd_ranks.append(bd_rank)

    betti = []
    for d in range(min(max_d + 1, len(omega_dims))):
        od = omega_dims[d]
        rd = bd_ranks[d] if d < len(bd_ranks) else 0
        rd1 = bd_ranks[d+1] if d+1 < len(bd_ranks) else 0
        betti.append(od - rd - rd1)

    return betti, omega_dims

def main():
    print("=" * 70)
    print("GLMY PATH HOMOLOGY EULER CHARACTERISTIC")
    print("kind-pasteur-2026-03-14-S69")
    print("=" * 70)

    for n in [3, 4, 5, 6]:
        print(f"\n{'='*70}")
        print(f"n = {n}")
        print(f"{'='*70}")

        chi_by_H = defaultdict(list)
        betti_by_H = defaultdict(list)
        omega_by_H = defaultdict(list)

        total_bits = n * (n - 1) // 2
        total = 2 ** total_bits
        max_check = total if n <= 5 else 2000

        count = 0
        for bits in range(total):
            count += 1
            if count > max_check:
                break

            A = bits_to_adj(bits, n)
            H = compute_H_dp(A, n)
            betti, omega_dims = compute_betti(A, n)
            chi = sum((-1)**k * b for k, b in enumerate(betti))

            chi_by_H[H].append(chi)
            betti_by_H[H].append(tuple(betti))
            omega_by_H[H].append(tuple(omega_dims))

        # Summary by H
        print(f"\n  H -> Betti numbers -> Euler characteristic:")
        for H in sorted(chi_by_H.keys()):
            chis = sorted(set(chi_by_H[H]))
            bettis = set(betti_by_H[H])
            print(f"\n  H={H:3d} ({len(chi_by_H[H]):4d} tournaments):")
            print(f"    chi values: {chis}")
            for b in sorted(bettis)[:4]:
                chi_b = sum((-1)**k * v for k, v in enumerate(b))
                print(f"    beta = {list(b)}, chi = {chi_b}")

        # Global chi analysis
        all_chis = []
        for chis in chi_by_H.values():
            all_chis.extend(chis)

        print(f"\n  GLOBAL CHI ANALYSIS:")
        print(f"    Distinct chi values: {sorted(set(all_chis))}")
        print(f"    Chi distribution: {dict(Counter(all_chis).most_common())}")

        chi_is_1 = all(c == 1 for c in all_chis)
        print(f"    Chi always = 1? {chi_is_1}")

        # Check: chi = f(H) for some function?
        chi_determines = True
        for H in chi_by_H:
            if len(set(chi_by_H[H])) > 1:
                chi_determines = False
                break
        print(f"    Chi determined by H? {chi_determines}")

        # Check: chi mod 2
        print(f"    Chi always odd? {all(c % 2 == 1 for c in all_chis)}")

        # Omega dimensions (vector space dims before quotienting)
        print(f"\n  OMEGA DIMENSIONS (path space sizes):")
        for H in sorted(omega_by_H.keys()):
            omegas = set(omega_by_H[H])
            if len(omegas) <= 3:
                for o in sorted(omegas):
                    print(f"    H={H:3d}: Omega = {list(o)}")
            else:
                print(f"    H={H:3d}: {len(omegas)} distinct Omega sequences")

        # Chi vs H correlation
        if len(set(all_chis)) > 1:
            H_vals = []
            chi_vals = []
            for H in chi_by_H:
                for c in chi_by_H[H]:
                    H_vals.append(H)
                    chi_vals.append(c)
            corr = np.corrcoef(H_vals, chi_vals)[0,1]
            print(f"\n  Correlation(H, chi) = {corr:.6f}")

    # Check at n=7 (sample)
    print(f"\n{'='*70}")
    print(f"n = 7 (sample 500)")
    print(f"{'='*70}")

    n = 7
    total_bits = n * (n - 1) // 2
    np.random.seed(42)

    chi_by_H = defaultdict(list)
    betti_by_H = defaultdict(list)

    for trial in range(500):
        bits = np.random.randint(0, 1 << total_bits)
        A = bits_to_adj(bits, n)
        H = compute_H_dp(A, n)
        betti, _ = compute_betti(A, n, max_d=6)
        chi = sum((-1)**k * b for k, b in enumerate(betti))

        chi_by_H[H].append(chi)
        betti_by_H[H].append(tuple(betti))

    print(f"\n  H -> chi (sample):")
    for H in sorted(chi_by_H.keys()):
        chis = sorted(set(chi_by_H[H]))
        print(f"    H={H:3d}: chi in {chis}, count={len(chi_by_H[H])}")

    all_chis = []
    for chis in chi_by_H.values():
        all_chis.extend(chis)

    print(f"\n  GLOBAL (n=7):")
    print(f"    Distinct chi: {sorted(set(all_chis))}")
    print(f"    Chi distribution: {dict(Counter(all_chis).most_common(10))}")

    H_vals = []
    chi_vals = []
    for H in chi_by_H:
        for c in chi_by_H[H]:
            H_vals.append(H)
            chi_vals.append(c)
    corr = np.corrcoef(H_vals, chi_vals)[0,1]
    print(f"    Correlation(H, chi) = {corr:.6f}")

    print(f"\n{'='*70}")
    print("DONE")
    print(f"{'='*70}")

if __name__ == '__main__':
    main()
