"""
chi_special_tournaments.py -- kind-pasteur-2026-03-14-S69
Compute GLMY Euler characteristic for special tournaments:
- Paley T_7 (beta_4 = 6, expect chi = 7)
- Regular n=7 tournaments (H-maximizers)
- H-maximizer at n=6 (H=45)
- 3-cycle reversal invariance at n=4

Also: deeper analysis of chi = 1 - beta_1 - beta_3 + beta_4 formula.
"""

import numpy as np
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

def paley_tournament(p):
    """Construct Paley tournament on p vertices (p prime, p = 3 mod 4)."""
    # QR = quadratic residues mod p
    qr = set()
    for x in range(1, p):
        qr.add((x * x) % p)
    # T[i][j] = 1 iff (j-i) mod p is a QR
    A = np.zeros((p, p), dtype=int)
    for i in range(p):
        for j in range(p):
            if i != j and ((j - i) % p) in qr:
                A[i][j] = 1
    return A

def compute_H_dp(A, n):
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

def compute_betti_full(A, n, max_d=None):
    """Full GLMY Betti computation with junk-face projection."""
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
    print("SPECIAL TOURNAMENT EULER CHARACTERISTICS")
    print("kind-pasteur-2026-03-14-S69")
    print("=" * 70)

    # ====================================
    # PALEY T_7
    # ====================================
    print("\n--- PALEY T_7 (7 vertices, H=189, regular) ---")
    A7 = paley_tournament(7)
    H7 = compute_H_dp(A7, 7)
    print(f"  H(T_7) = {H7}")
    betti7, omega7 = compute_betti_full(A7, 7, max_d=6)
    chi7 = sum((-1)**k * b for k, b in enumerate(betti7))
    print(f"  Betti = {betti7}")
    print(f"  Omega dims = {omega7}")
    print(f"  Chi = {chi7}")
    print(f"  Expected: beta_4 = 6, chi = 1 + 0 + 0 + 6 = 7")

    # ====================================
    # PALEY T_3
    # ====================================
    print("\n--- PALEY T_3 (3 vertices, H=3, 3-cycle) ---")
    A3 = paley_tournament(3)
    H3 = compute_H_dp(A3, 3)
    betti3, omega3 = compute_betti_full(A3, 3)
    chi3 = sum((-1)**k * b for k, b in enumerate(betti3))
    print(f"  H(T_3) = {H3}")
    print(f"  Betti = {betti3}")
    print(f"  Chi = {chi3}")

    # ====================================
    # TRANSITIVE T_n
    # ====================================
    for n in [3, 4, 5, 6, 7]:
        print(f"\n--- TRANSITIVE T_{n} (H=1) ---")
        A = np.zeros((n, n), dtype=int)
        for i in range(n):
            for j in range(i+1, n):
                A[i][j] = 1
        betti, omega = compute_betti_full(A, n)
        chi = sum((-1)**k * b for k, b in enumerate(betti))
        print(f"  Betti = {betti}")
        print(f"  Omega = {omega}")
        print(f"  Chi = {chi}")

    # ====================================
    # n=4 EXHAUSTIVE: verify 3-cycle-reversal invariance for chi
    # ====================================
    print("\n--- n=4: 3-CYCLE REVERSAL and CHI ---")
    n = 4
    tb = n * (n - 1) // 2
    chi_change_count = 0
    total_reversals = 0

    from itertools import combinations, permutations

    for bits in range(1 << tb):
        A = bits_to_adj(bits, n)
        betti_orig, _ = compute_betti_full(A, n)
        chi_orig = sum((-1)**k * b for k, b in enumerate(betti_orig))

        # Find all 3-cycles
        for a, b, c in combinations(range(n), 3):
            for perm in permutations([a, b, c]):
                x, y, z = perm
                if A[x][y] == 1 and A[y][z] == 1 and A[z][x] == 1:
                    # Reverse the 3-cycle
                    A_rev = A.copy()
                    A_rev[x][y] = 0; A_rev[y][x] = 1
                    A_rev[y][z] = 0; A_rev[z][y] = 1
                    A_rev[z][x] = 0; A_rev[x][z] = 1

                    betti_rev, _ = compute_betti_full(A_rev, n)
                    chi_rev = sum((-1)**k * b for k, b in enumerate(betti_rev))

                    total_reversals += 1
                    if chi_orig != chi_rev:
                        chi_change_count += 1
                    break  # one cycle per triple

    print(f"  Total 3-cycle reversals: {total_reversals}")
    print(f"  Chi changes: {chi_change_count}")
    print(f"  Chi invariant under 3-cycle reversal: {chi_change_count == 0}")

    # ====================================
    # n=5: check chi vs beta_1 relationship
    # ====================================
    print("\n--- n=5: CHI vs BETA_1, BETA_3 ---")
    n = 5
    tb = n * (n - 1) // 2

    chi_beta_relation = Counter()
    c3_chi = defaultdict(list)

    for bits in range(1 << tb):
        A = bits_to_adj(bits, n)
        betti, _ = compute_betti_full(A, n)
        chi = sum((-1)**k * b for k, b in enumerate(betti))

        b1 = betti[1] if len(betti) > 1 else 0
        b3 = betti[3] if len(betti) > 3 else 0

        chi_beta_relation[(b1, b3, chi)] += 1

        # Count 3-cycles
        c3 = int(np.trace(A @ A @ A)) // 3
        c3_chi[c3].append(chi)

    print(f"  (beta_1, beta_3, chi) distribution:")
    for key, count in sorted(chi_beta_relation.items()):
        print(f"    beta_1={key[0]}, beta_3={key[1]}, chi={key[2]}: {count} tournaments")

    print(f"\n  Formula check: chi = 1 - beta_1 + beta_3?")
    for key, count in sorted(chi_beta_relation.items()):
        b1, b3, chi = key
        expected = 1 - b1 + b3
        match = (expected == chi)
        print(f"    beta_1={b1}, beta_3={b3}: chi={chi}, 1-b1+b3={expected}, match={match}")

    print(f"\n  c3 vs chi:")
    for c3 in sorted(c3_chi.keys()):
        vals = c3_chi[c3]
        dist = Counter(vals)
        print(f"    c3={c3}: chi distribution = {dict(dist)}")

    # ====================================
    # n=6 exhaustive: chi vs beta detailed
    # ====================================
    print("\n--- n=6: CHI vs BETA (sample 5000) ---")
    n = 6
    tb = n * (n - 1) // 2
    np.random.seed(123)

    chi_beta6 = Counter()
    for _ in range(5000):
        bits = np.random.randint(0, 1 << tb)
        A = bits_to_adj(bits, n)
        betti, _ = compute_betti_full(A, n, max_d=5)
        chi = sum((-1)**k * b for k, b in enumerate(betti))

        b1 = betti[1] if len(betti) > 1 else 0
        b3 = betti[3] if len(betti) > 3 else 0
        b4 = betti[4] if len(betti) > 4 else 0

        chi_beta6[(b1, b3, b4, chi)] += 1

    print(f"  (beta_1, beta_3, beta_4, chi) distribution:")
    for key, count in sorted(chi_beta6.items()):
        b1, b3, b4, chi = key
        expected = 1 - b1 + b3 - b4
        match = (expected == chi)
        print(f"    b1={b1}, b3={b3}, b4={b4}: chi={chi}, 1-b1+b3-b4={expected}, match={match}, count={count}")

    # ====================================
    # KEY QUESTION: Does chi relate to number of 3-cycles?
    # ====================================
    print("\n--- n=5: EXACT relationship chi vs c3 ---")
    print("  chi = 0 iff beta_1 = 1 (from above)")
    print("  beta_1 = 1 iff tournament has a 3-cycle, i.e., c3 >= 1")

    n = 5
    tb = n * (n - 1) // 2
    for bits in range(1 << tb):
        A = bits_to_adj(bits, n)
        c3 = int(np.trace(A @ A @ A)) // 3
        betti, _ = compute_betti_full(A, n)
        b1 = betti[1] if len(betti) > 1 else 0

        if (c3 == 0 and b1 != 0) or (c3 > 0 and b1 == 0 and betti[3] == 0):
            # Check if this is an exception
            H = compute_H_dp(A, n)
            print(f"  EXCEPTION: bits={bits}, c3={c3}, beta_1={b1}, H={H}, betti={betti}")

    # Also check at n=4
    print("\n--- n=4: beta_1 vs c3 ---")
    n = 4
    tb = n * (n - 1) // 2
    b1_c3 = Counter()
    for bits in range(1 << tb):
        A = bits_to_adj(bits, n)
        c3 = int(np.trace(A @ A @ A)) // 3
        betti, _ = compute_betti_full(A, n)
        b1 = betti[1] if len(betti) > 1 else 0
        b1_c3[(c3, b1)] += 1

    for key, count in sorted(b1_c3.items()):
        print(f"    c3={key[0]}, beta_1={key[1]}: {count} tournaments")

    print(f"\n  At n=4: beta_1 = 1 iff c3 >= 2? or c3 >= 1?")

    print("\n" + "=" * 70)
    print("ANALYSIS COMPLETE")
    print("=" * 70)

if __name__ == '__main__':
    main()
