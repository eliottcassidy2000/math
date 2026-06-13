"""
poincare_duality_tournament.py -- kind-pasteur-2026-03-14-S69
Explore Poincare duality for tournament path homology.

Key observation: Paley T_7 has Omega = [7,21,42,63,63,42,21] which is
ALMOST palindromic (missing the last term). The cyclic interval T_5 has
Omega = [5,10,10,10,5] which IS palindromic.

In topology, Poincare duality for a closed oriented n-manifold:
  H_k(M) ≅ H_{n-k}(M)  =>  beta_k = beta_{n-k}

For tournament path homology:
  Does Omega_k = Omega_{n-1-k}? (path space duality)
  Does beta_k = beta_{n-1-k}? (homology duality)

Explore: which tournaments have palindromic Omega or Betti?
What structural property characterizes palindromy?
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

def compute_omega_and_betti(A, n, max_d=None):
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

    # Raw path counts (before junk projection)
    raw_counts = [len(paths.get(d, [])) for d in range(n)]

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

    return omega_dims, betti, raw_counts

def main():
    print("=" * 70)
    print("POINCARE DUALITY FOR TOURNAMENT PATH HOMOLOGY")
    print("kind-pasteur-2026-03-14-S69")
    print("=" * 70)

    # ============================================
    # PART 1: Raw path counts vs Omega dims
    # ============================================
    print("\n" + "=" * 70)
    print("PART 1: RAW PATH COUNTS vs OMEGA DIMS")
    print("  Raw = #{d-paths}, Omega = #{regular d-paths}")
    print("  Transitive has raw = C(n, d+1) * (d+1)! = P(n, d+1)")
    print("=" * 70)

    for n in [3, 4, 5, 6, 7]:
        print(f"\n--- TRANSITIVE n={n} ---")
        A = np.zeros((n, n), dtype=int)
        for i in range(n):
            for j in range(i+1, n):
                A[i][j] = 1
        omega, betti, raw = compute_omega_and_betti(A, n)
        print(f"  Raw    = {raw}")
        print(f"  Omega  = {omega}")
        print(f"  Betti  = {betti}")
        print(f"  Raw palindromic: {raw == raw[::-1]}")
        print(f"  Omega palindromic: {omega == omega[::-1]}")

    # ============================================
    # PART 2: Palindrome classification at n=5
    # ============================================
    print("\n" + "=" * 70)
    print("PART 2: PALINDROME CLASSIFICATION at n=5")
    print("=" * 70)

    n = 5
    palin_raw = 0
    palin_omega = 0
    total = 2 ** (n*(n-1)//2)

    palin_by_score = defaultdict(list)
    nonpalin_by_score = defaultdict(list)

    for bits in range(total):
        A = bits_to_adj(bits, n)
        H = compute_H_dp(A, n)
        omega, betti, raw = compute_omega_and_betti(A, n)
        scores = tuple(sorted([sum(A[i]) for i in range(n)]))

        raw_p = (raw == raw[::-1])
        omega_p = (omega == omega[::-1])

        if raw_p:
            palin_raw += 1
        if omega_p:
            palin_omega += 1
            palin_by_score[scores].append((H, omega, betti))
        else:
            nonpalin_by_score[scores].append((H, omega, betti))

    print(f"  Total: {total}")
    print(f"  Raw palindromic: {palin_raw}")
    print(f"  Omega palindromic: {palin_omega}")

    print(f"\n  Palindromic tournaments by score sequence:")
    for scores, items in sorted(palin_by_score.items()):
        H_vals = set(it[0] for it in items)
        omega_vals = set(tuple(it[1]) for it in items)
        print(f"    score={scores}: {len(items)} tournaments, H in {sorted(H_vals)}")
        for o in sorted(omega_vals):
            print(f"      Omega = {list(o)}")

    print(f"\n  Which score sequences have palindromic tournaments:")
    for scores in sorted(set(list(palin_by_score.keys()) + list(nonpalin_by_score.keys()))):
        n_p = len(palin_by_score.get(scores, []))
        n_np = len(nonpalin_by_score.get(scores, []))
        print(f"    {scores}: {n_p} palindromic, {n_np} non-palindromic")

    # ============================================
    # PART 3: Check palindrome for circulant tournaments
    # ============================================
    print("\n" + "=" * 70)
    print("PART 3: CIRCULANT TOURNAMENT PALINDROMY")
    print("  Circulant: A[i][j] = 1 iff (j-i) mod n in S")
    print("=" * 70)

    for n in [5, 7, 9]:
        print(f"\n--- n = {n} ---")
        half = (n - 1) // 2

        # Generate some circulant connection sets
        from itertools import combinations
        all_S = []
        for S in combinations(range(1, n), half):
            # Check: S and n-S partition {1,...,n-1}
            S_set = set(S)
            complement = set((n - s) % n for s in S)
            if S_set == complement or len(S_set & complement) > 0:
                continue  # not a valid tournament connection set
            # Valid: S and {n-s : s in S} are disjoint and union = {1,...,n-1}
            if S_set.isdisjoint(complement) and len(S_set) + len(complement) == n - 1:
                all_S.append(S)

        print(f"  {len(all_S)} valid circulant connection sets")

        for S in all_S[:10]:  # limit
            A = np.zeros((n, n), dtype=int)
            for i in range(n):
                for s in S:
                    j = (i + s) % n
                    A[i][j] = 1

            H = compute_H_dp(A, n)
            max_d = min(n-1, 6)
            omega, betti, raw = compute_omega_and_betti(A, n, max_d=max_d)
            omega_p = (omega == omega[::-1])
            raw_p = (raw[:len(omega)] == raw[:len(omega)][::-1])

            print(f"    S={S}: H={H}, Omega={omega}, palindromic={omega_p}")

    # ============================================
    # PART 4: Poincare duality polynomial
    # ============================================
    print("\n" + "=" * 70)
    print("PART 4: POINCARE POLYNOMIAL P(T,t) = sum Omega_k * t^k")
    print("  Duality: P(T,t) = t^{n-1} * P(T, 1/t)?")
    print("  Equivalently: Omega_k = Omega_{n-1-k}")
    print("=" * 70)

    # Compute P(T, t) for special tournaments
    for label, A_fn, n_val in [
        ("Transitive", lambda n: np.triu(np.ones((n,n), dtype=int), 1), 7),
        ("Paley T_7", lambda n: None, 7),  # special
        ("Cyclic interval n=5", lambda n: None, 5),
    ]:
        n = n_val
        if label == "Paley T_7":
            A = np.zeros((7, 7), dtype=int)
            qr = set()
            for x in range(1, 7):
                qr.add((x*x) % 7)
            for i in range(7):
                for j in range(7):
                    if i != j and ((j-i) % 7) in qr:
                        A[i][j] = 1
        elif label == "Cyclic interval n=5":
            A = np.zeros((5, 5), dtype=int)
            for i in range(5):
                for d in [1, 2]:
                    A[i][(i+d) % 5] = 1
        else:
            A = A_fn(n)

        omega, betti, raw = compute_omega_and_betti(A, n)
        chi = sum((-1)**k * b for k, b in enumerate(betti))
        H = compute_H_dp(A, n)

        print(f"\n  {label}:")
        print(f"    n={n}, H={H}")
        print(f"    Omega = {omega}")
        print(f"    Betti = {betti}")
        print(f"    Chi = {chi}")
        print(f"    sum(Omega) = {sum(omega)}")
        print(f"    Alt sum((-1)^k * Omega_k) = {sum((-1)**k * o for k, o in enumerate(omega))}")
        print(f"    P(1) = {sum(omega)} (total dim of path space)")

        # Check palindrome
        if omega == omega[::-1]:
            print(f"    *** Omega is PALINDROMIC (Poincare duality) ***")
        else:
            # Check "almost palindromic" (shifted by last element)
            if len(omega) >= 2 and omega[:-1] == omega[-2::-1]:
                print(f"    Omega is 'almost palindromic' (Omega[0:-1] reversed = Omega[-2::-1])")

    # ============================================
    # PART 5: What determines the Poincare polynomial?
    # ============================================
    print("\n" + "=" * 70)
    print("PART 5: WHAT DETERMINES OMEGA (POINCARE POLYNOMIAL)?")
    print("  Is Omega determined by H? By score? By cycle counts?")
    print("=" * 70)

    for n in [5, 6]:
        print(f"\n--- n = {n} ---")
        omega_by_H = defaultdict(set)
        omega_by_score = defaultdict(set)
        c3_by_omega = defaultdict(set)

        count = 0
        for bits in range(2**(n*(n-1)//2)):
            count += 1
            if n == 6 and count > 3000:
                break

            A = bits_to_adj(bits, n)
            H = compute_H_dp(A, n)
            scores = tuple(sorted([sum(A[i]) for i in range(n)]))
            c3 = int(np.trace(A @ A @ A)) // 3

            omega, _, _ = compute_omega_and_betti(A, n)
            omega_t = tuple(omega)

            omega_by_H[H].add(omega_t)
            omega_by_score[scores].add(omega_t)
            c3_by_omega[omega_t].add(c3)

        print(f"  H -> #distinct Omega:")
        for H in sorted(omega_by_H.keys()):
            print(f"    H={H:3d}: {len(omega_by_H[H])} distinct Omega")

        print(f"\n  score -> #distinct Omega:")
        for s in sorted(omega_by_score.keys()):
            print(f"    {s}: {len(omega_by_score[s])}")

    print("\n" + "=" * 70)
    print("DONE")
    print("=" * 70)

if __name__ == '__main__':
    main()
