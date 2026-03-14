"""
chi_paley_pattern.py -- kind-pasteur-2026-03-14-S69
Test: does chi(Paley T_p) = p for all Paley primes p?

Known: chi(T_3) = 0, chi(T_7) = 7.
If chi(T_p) = p, this would be a remarkable connection:
  chi = p = 1 + beta_4 (since beta_1=beta_2=beta_3=0 for Paley at p>=7)

Also explore: chi for non-Paley tournaments at the same n.
What about regular non-Paley tournaments?

And: verify Vassiliev type at n=6 (expected: type 4).
"""

import numpy as np
from collections import Counter, defaultdict
import sys

sys.stdout.reconfigure(encoding='utf-8')

def paley_tournament(p):
    qr = set()
    for x in range(1, p):
        qr.add((x * x) % p)
    A = np.zeros((p, p), dtype=int)
    for i in range(p):
        for j in range(p):
            if i != j and ((j - i) % p) in qr:
                A[i][j] = 1
    return A

def compute_betti_full(A, n, max_d=None):
    if max_d is None:
        max_d = n - 1

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

        C = np.zeros((n_junk, n_d), dtype=float)
        for j, jf in enumerate(face_junk):
            for face, sign in jf:
                C[junk_idx[face], j] += sign

        rank_c = int(np.linalg.matrix_rank(C)) if n_junk > 0 else 0
        omega_d = n_d - rank_c
        omega_dims.append(omega_d)

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
    print("CHI(PALEY) PATTERN TEST")
    print("kind-pasteur-2026-03-14-S69")
    print("=" * 70)

    # Paley primes p = 3 mod 4
    for p in [3, 7]:
        print(f"\n--- Paley T_{p} ---")
        A = paley_tournament(p)
        betti, omega = compute_betti_full(A, p)
        chi = sum((-1)**k * b for k, b in enumerate(betti))
        print(f"  Betti = {betti}")
        print(f"  Omega = {omega}")
        print(f"  Chi = {chi}")
        print(f"  Chi == p? {chi == p}")
        palin = (omega == omega[::-1])
        print(f"  Omega palindromic? {palin}")

    # Check transitive tournament at various n for comparison
    print("\n--- TRANSITIVE vs PALEY chi ---")
    for n in [3, 5, 7]:
        # Transitive
        A_trans = np.zeros((n, n), dtype=int)
        for i in range(n):
            for j in range(i+1, n):
                A_trans[i][j] = 1
        b_trans, o_trans = compute_betti_full(A_trans, n)
        chi_trans = sum((-1)**k * b for k, b in enumerate(b_trans))
        print(f"  n={n} transitive: chi={chi_trans}, betti={b_trans}")

    # Random regular tournaments at n=7 for comparison
    print("\n--- RANDOM n=7 tournaments: chi distribution ---")
    n = 7
    np.random.seed(42)
    chi_dist = Counter()

    for trial in range(100):
        bits = np.random.randint(0, 1 << (n*(n-1)//2))
        A = np.zeros((n, n), dtype=int)
        idx = 0
        for i in range(n):
            for j in range(i+1, n):
                if bits & (1 << idx):
                    A[i][j] = 1
                else:
                    A[j][i] = 1
                idx += 1

        betti, omega = compute_betti_full(A, n)
        chi = sum((-1)**k * b for k, b in enumerate(betti))
        chi_dist[chi] += 1

    print(f"  Chi distribution (100 random n=7): {dict(sorted(chi_dist.items()))}")

    # Cyclic interval tournament (should also be interesting)
    print("\n--- CYCLIC INTERVAL at various n ---")
    for n in [5, 7, 9]:
        A = np.zeros((n, n), dtype=int)
        half = (n - 1) // 2
        for i in range(n):
            for d in range(1, half + 1):
                j = (i + d) % n
                A[i][j] = 1

        betti, omega = compute_betti_full(A, n, max_d=min(n-1, 6))
        chi = sum((-1)**k * b for k, b in enumerate(betti))
        print(f"  n={n} cyclic interval: chi={chi}, betti={betti}")
        print(f"    Omega = {omega}")
        palin = (omega == omega[::-1])
        print(f"    Palindromic? {palin}")

    print("\n" + "=" * 70)
    print("DONE")
    print("=" * 70)

if __name__ == '__main__':
    main()
