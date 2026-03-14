"""
tropical_quantum_tournament.py -- kind-pasteur-2026-03-14-S71
Creative exploration: tropical geometry and quantum mechanics of tournaments.

TROPICAL TOURNAMENTS:
- Replace (×, +) with (min, +) in the path weight algebra
- The "tropical permanent" = minimum weight Hamiltonian path
- The "tropical H" = number of minimum weight paths (or the weight itself)
- Tournament as a weighted digraph with weights 0/1 on arcs
- The tropical semiring perspective on H(T)

QUANTUM TOURNAMENTS:
- Replace T[i][j] ∈ {0,1} with amplitudes z[i][j] ∈ C
- The "quantum H" = permanent of the amplitude matrix
- Connection to BosonSampling (Aaronson-Arkhipov)
- The tournament permanent vs. determinant

GAME THEORY:
- Tournament = preference ordering (i beats j means i is preferred over j)
- H(T) = number of consistent total orderings
- Connection to Kemeny distance and social choice

INFORMATION THEORY:
- H as a measure of "how close to a total order" the tournament is
- Entropy of the tournament
- Mutual information between arc variables
"""

import numpy as np
from itertools import permutations, combinations
from collections import Counter, defaultdict
import sys, math

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

def permanent(M):
    """Compute permanent of matrix M (Ryser's formula for small n)."""
    n = len(M)
    if n == 0:
        return 1
    if n == 1:
        return M[0][0]
    # Ryser's formula
    total = 0
    for S in range(1, 1 << n):
        col_sums = np.zeros(n, dtype=complex)
        bits_set = bin(S).count('1')
        for j in range(n):
            if S & (1 << j):
                for i in range(n):
                    col_sums[i] += M[i][j]
        prod = 1
        for i in range(n):
            prod *= col_sums[i]
        if (n - bits_set) % 2 == 1:
            total -= prod
        else:
            total += prod
    return total * ((-1)**n)

def main():
    print("=" * 70)
    print("TROPICAL & QUANTUM TOURNAMENTS")
    print("kind-pasteur-2026-03-14-S71")
    print("=" * 70)

    # ========================================
    # PART 1: PERMANENT OF TOURNAMENT ADJACENCY
    # ========================================
    print("\n" + "=" * 70)
    print("PART 1: PERMANENT OF TOURNAMENT ADJACENCY MATRIX")
    print("  perm(A) = sum over permutations sigma: prod A[i][sigma(i)]")
    print("  This counts CYCLE COVERS, not Hamiltonian paths!")
    print("  perm(A) = sum over all cycle covers of T")
    print("  Compare: H(T) = sum over Ham PATHS (not cycle covers)")
    print("=" * 70)

    for n in [3, 4, 5]:
        print(f"\n--- n = {n} ---")
        perm_by_H = defaultdict(list)

        for bits in range(2**(n*(n-1)//2)):
            A = bits_to_adj(bits, n)
            H = compute_H_dp(A, n)
            p = round(permanent(A.astype(float)).real)
            perm_by_H[H].append(p)

        for H in sorted(perm_by_H.keys()):
            vals = sorted(set(perm_by_H[H]))
            print(f"  H={H:3d}: perm(A) in {vals}")

        # KEY: Does perm(A) relate to H?
        # perm(A) counts cycle covers. H counts ham paths.
        # A cycle cover of a tournament = collection of disjoint directed cycles
        # covering all vertices. This is exactly a permutation where every cycle
        # is a directed cycle of T.

        # Connection to det(A):
        # det(A) = sum sgn(sigma) prod A[i][sigma(i)]
        # perm(A) = sum prod A[i][sigma(i)] (no sign)
        # det(A) = perm_signed(A)

    # ========================================
    # PART 2: PERMANENT vs DETERMINANT vs H
    # ========================================
    print("\n" + "=" * 70)
    print("PART 2: PERMANENT, DETERMINANT, and H")
    print("  perm(A) = unsigned cycle cover count")
    print("  det(A) = signed cycle cover count")
    print("  H(T) = Hamiltonian path count")
    print("  Question: Is there a FORMULA relating them?")
    print("=" * 70)

    for n in [3, 4, 5]:
        print(f"\n--- n = {n} ---")
        for bits in range(2**(n*(n-1)//2)):
            A = bits_to_adj(bits, n)
            H = compute_H_dp(A, n)
            det_A = round(np.linalg.det(A.astype(float)))
            perm_A = round(permanent(A.astype(float)).real)

            # Also: permanent of J-I-A = permanent of A^op (complement)
            A_op = 1 - A - np.eye(n, dtype=int)
            perm_Aop = round(permanent(A_op.astype(float)).real)

            # H and perm of the "path matrix"?
            # Actually, H(T) = perm of the transfer matrix M? No, M is n×n.
            # H = sum of all M[i][j] = trace of J*M... no.

            if bits <= 10 or H == max(compute_H_dp(bits_to_adj(b, n), n) for b in [bits]):
                # Print a few examples
                pass  # too noisy

        # Correlation
        data = []
        for bits in range(2**(n*(n-1)//2)):
            A = bits_to_adj(bits, n)
            H = compute_H_dp(A, n)
            perm_A = round(permanent(A.astype(float)).real)
            det_A = round(np.linalg.det(A.astype(float)))
            data.append((H, perm_A, det_A))

        H_vals = [d[0] for d in data]
        perm_vals = [d[1] for d in data]
        det_vals = [d[2] for d in data]

        corr_hp = np.corrcoef(H_vals, perm_vals)[0,1]
        corr_hd = np.corrcoef(H_vals, det_vals)[0,1]
        print(f"  Correlation(H, perm(A)) = {corr_hp:.6f}")
        print(f"  Correlation(H, det(A)) = {corr_hd:.6f}")

        # CHECK: Is perm(A) + perm(A^op) constant?
        # For tournament: A + A^op = J - I. So perm(A) + perm(A^op) relates to
        # the permanent of J-I (which is the number of derangements times something)
        perm_sum = set()
        for bits in range(2**(n*(n-1)//2)):
            A = bits_to_adj(bits, n)
            A_op = 1 - A - np.eye(n, dtype=int)
            s = round(permanent(A.astype(float)).real) + round(permanent(A_op.astype(float)).real)
            perm_sum.add(s)
        print(f"  perm(A) + perm(A^op) values: {sorted(perm_sum)}")

        # And H(T) + H(T^op)? = 2*H(T) since H(T^op) = H(T)
        # What about perm(A) - perm(A^op)?
        perm_diff_by_H = defaultdict(set)
        for bits in range(2**(n*(n-1)//2)):
            A = bits_to_adj(bits, n)
            H = compute_H_dp(A, n)
            A_op = 1 - A - np.eye(n, dtype=int)
            diff = round(permanent(A.astype(float)).real) - round(permanent(A_op.astype(float)).real)
            perm_diff_by_H[H].add(diff)
        print(f"  H -> perm(A)-perm(A^op):")
        for H in sorted(perm_diff_by_H.keys()):
            print(f"    H={H:3d}: diff in {sorted(perm_diff_by_H[H])}")

    # ========================================
    # PART 3: TROPICAL TOURNAMENT — Shortest Ham Path
    # ========================================
    print("\n" + "=" * 70)
    print("PART 3: TROPICAL TOURNAMENT — SHORTEST HAMILTONIAN PATH")
    print("  Assign weight w(i->j) = 0 if i->j, infinity if j->i")
    print("  Then min-weight Ham path = min # 'wrong-direction' arcs = 0")
    print("  More interesting: w(i->j) = |i-j| (displacement cost)")
    print("  Or: w(i->j) = 1 - T[i][j] (number of backward arcs)")
    print("=" * 70)

    for n in [4, 5, 6]:
        print(f"\n--- n = {n} ---")
        backward_by_H = defaultdict(list)
        min_backward_by_H = defaultdict(list)

        count = 0
        for bits in range(2**(n*(n-1)//2)):
            count += 1
            if n >= 6 and count > 3000:
                break

            A = bits_to_adj(bits, n)
            H = compute_H_dp(A, n)

            # For each Ham path, count backward arcs (j->i where i<j)
            min_back = n
            total_back = 0
            num_paths = 0
            for perm in permutations(range(n)):
                valid = True
                for i in range(n-1):
                    if A[perm[i]][perm[i+1]] != 1:
                        valid = False
                        break
                if valid:
                    back = sum(1 for i in range(n-1) if perm[i] > perm[i+1])
                    min_back = min(min_back, back)
                    total_back += back
                    num_paths += 1

            if num_paths > 0:
                backward_by_H[H].append(total_back / num_paths)
                min_backward_by_H[H].append(min_back)

        print(f"  H -> min backward arcs in any Ham path:")
        for H in sorted(min_backward_by_H.keys())[:15]:
            vals = sorted(set(min_backward_by_H[H]))
            mean_back = np.mean(backward_by_H[H])
            print(f"    H={H:3d}: min_back in {vals}, mean_back={mean_back:.2f}")

    # ========================================
    # PART 4: KEMENY DISTANCE AND SOCIAL CHOICE
    # ========================================
    print("\n" + "=" * 70)
    print("PART 4: KEMENY DISTANCE — TOURNAMENT AS VOTING SYSTEM")
    print("  Kemeny distance = min # arc flips to make transitive")
    print("  = min # pairwise comparisons to change for a total order")
    print("  Connection to H: more transitive (low Kemeny) = low H?")
    print("=" * 70)

    for n in [4, 5, 6]:
        print(f"\n--- n = {n} ---")
        kemeny_by_H = defaultdict(list)

        count = 0
        for bits in range(2**(n*(n-1)//2)):
            count += 1
            if n >= 6 and count > 3000:
                break

            A = bits_to_adj(bits, n)
            H = compute_H_dp(A, n)

            # Kemeny distance = min over all total orders sigma:
            # #{(i,j) : sigma(i) < sigma(j) but A[j][i] = 1}
            # = min over permutations: # disagreements
            min_kemeny = n * (n-1) // 2
            for perm in permutations(range(n)):
                disagree = 0
                for i in range(n):
                    for j in range(i+1, n):
                        # In total order perm: perm^{-1}(i) < perm^{-1}(j)?
                        # Position of vertex i in perm
                        pass
                # Actually: Kemeny distance = C(n,2) - max # agreements
                # Agreement with perm: arc i->j agrees if i appears before j in perm
                agreements = 0
                for pos_a in range(n):
                    for pos_b in range(pos_a+1, n):
                        if A[perm[pos_a]][perm[pos_b]] == 1:
                            agreements += 1
                kemeny = n*(n-1)//2 - agreements
                min_kemeny = min(min_kemeny, kemeny)

            kemeny_by_H[H].append(min_kemeny)

        print(f"  H -> Kemeny distance:")
        for H in sorted(kemeny_by_H.keys())[:15]:
            vals = sorted(set(kemeny_by_H[H]))
            print(f"    H={H:3d}: Kemeny in {vals}")

        # Correlation
        H_vals = []
        K_vals = []
        for H, kems in kemeny_by_H.items():
            for k in kems:
                H_vals.append(H)
                K_vals.append(k)
        if len(set(K_vals)) > 1:
            corr = np.corrcoef(H_vals, K_vals)[0,1]
            print(f"\n  Correlation(H, Kemeny) = {corr:.6f}")

    # ========================================
    # PART 5: TOURNAMENT ENTROPY
    # ========================================
    print("\n" + "=" * 70)
    print("PART 5: TOURNAMENT ENTROPY AND INFORMATION CONTENT")
    print("  Score entropy: S = -sum (s_i/total) log(s_i/total)")
    print("  H-entropy: log2(H) measures 'disorder' of tournament")
    print("  Question: Is log2(H) ~ n * score_entropy?")
    print("=" * 70)

    for n in [5, 6, 7]:
        print(f"\n--- n = {n} ---")
        data = []
        count = 0

        for bits in range(2**(n*(n-1)//2)):
            count += 1
            if count > 5000:
                break
            if n >= 7 and count > 1000:
                break

            A = bits_to_adj(bits, n)
            H = compute_H_dp(A, n)
            scores = [sum(A[i]) for i in range(n)]
            total = sum(scores)

            # Score entropy
            if total > 0:
                probs = [s/total for s in scores if s > 0]
                entropy = -sum(p * math.log2(p) for p in probs if p > 0)
            else:
                entropy = 0

            # Score variance
            score_var = np.var(scores)

            # c3 count
            c3 = int(np.trace(A @ A @ A)) // 3

            log_H = math.log2(H) if H > 0 else 0

            data.append({
                'H': H, 'logH': log_H, 'entropy': entropy,
                'score_var': score_var, 'c3': c3,
            })

        H_vals = np.array([d['H'] for d in data])
        logH_vals = np.array([d['logH'] for d in data])
        ent_vals = np.array([d['entropy'] for d in data])
        var_vals = np.array([d['score_var'] for d in data])
        c3_vals = np.array([d['c3'] for d in data])

        print(f"  Correlation(logH, score_entropy) = {np.corrcoef(logH_vals, ent_vals)[0,1]:.6f}")
        print(f"  Correlation(logH, score_variance) = {np.corrcoef(logH_vals, var_vals)[0,1]:.6f}")
        print(f"  Correlation(logH, c3) = {np.corrcoef(logH_vals, c3_vals)[0,1]:.6f}")
        print(f"  Correlation(H, score_entropy) = {np.corrcoef(H_vals, ent_vals)[0,1]:.6f}")
        print(f"  Correlation(H, score_variance) = {np.corrcoef(H_vals, var_vals)[0,1]:.6f}")

        # For H-maximizers
        max_H = max(H_vals)
        max_data = [d for d in data if d['H'] == max_H]
        if max_data:
            print(f"\n  H-maximizer (H={max_H}):")
            print(f"    score_entropy = {max_data[0]['entropy']:.4f}")
            print(f"    score_var = {max_data[0]['score_var']:.4f}")
            print(f"    c3 = {max_data[0]['c3']}")
            print(f"    log2(H) = {max_data[0]['logH']:.4f}")

        # Max entropy tournament?
        max_ent = max(ent_vals)
        max_ent_data = [d for d in data if abs(d['entropy'] - max_ent) < 0.001]
        if max_ent_data:
            H_at_max_ent = [d['H'] for d in max_ent_data]
            print(f"  Max entropy = {max_ent:.4f}: H in {sorted(set(H_at_max_ent))[:5]}")

    # ========================================
    # PART 6: QUANTUM PERMANENT — COMPLEX AMPLITUDES
    # ========================================
    print("\n" + "=" * 70)
    print("PART 6: QUANTUM TOURNAMENT — COMPLEX PERMANENT")
    print("  Replace T[i][j] with z[i][j] = T[i][j] * e^{i*theta}")
    print("  Quantum H = permanent of amplitude matrix")
    print("  Connection to BosonSampling")
    print("=" * 70)

    n = 4
    print(f"\n--- n = {n}, phase theta ---")
    # For each tournament, compute perm(A * e^{i*theta}) as function of theta
    for bits in [0, 15, 63]:  # a few tournaments
        A = bits_to_adj(bits, n)
        H = compute_H_dp(A, n)
        scores = tuple(sorted([sum(A[i]) for i in range(n)]))

        print(f"\n  bits={bits}, H={H}, scores={scores}")

        for theta in [0, np.pi/6, np.pi/4, np.pi/3, np.pi/2, np.pi]:
            # Amplitude matrix: z[i][j] = A[i][j] * e^{i*theta} + (1-A[i][j]) * 0
            Z = A.astype(complex) * np.exp(1j * theta)
            p = permanent(Z)
            print(f"    theta={theta:.4f}: perm = {p.real:.4f} + {p.imag:.4f}i, |perm|={abs(p):.4f}")

        # Also: the "unitary tournament" — normalize rows to unit norm
        # Each row of A has out-degree d_i arcs. Normalize to d_i^{-1/2}
        # to get a sub-unitary matrix

    # ========================================
    # PART 7: THE H-LANDSCAPE ON THE HYPERCUBE
    # ========================================
    print("\n" + "=" * 70)
    print("PART 7: H-LANDSCAPE ON THE TOURNAMENT HYPERCUBE")
    print("  Tournaments = vertices of C(n,2)-dim hypercube")
    print("  H(T) is a function on this hypercube")
    print("  Question: How many local maxima? Saddle points?")
    print("=" * 70)

    for n in [4, 5]:
        print(f"\n--- n = {n} ---")
        m = n * (n-1) // 2

        local_max = 0
        local_min = 0
        saddle = 0
        total = 2**m

        for bits in range(total):
            A = bits_to_adj(bits, n)
            H = compute_H_dp(A, n)

            # Neighbors: all tournaments differing in exactly one arc
            is_max = True
            is_min = True
            for arc_bit in range(m):
                nbr_bits = bits ^ (1 << arc_bit)
                A_nbr = bits_to_adj(nbr_bits, n)
                H_nbr = compute_H_dp(A_nbr, n)
                if H_nbr > H:
                    is_max = False
                if H_nbr < H:
                    is_min = False

            if is_max:
                local_max += 1
            if is_min:
                local_min += 1
            if not is_max and not is_min:
                saddle += 1

        print(f"  Total tournaments: {total}")
        print(f"  Local maxima: {local_max} ({100*local_max/total:.1f}%)")
        print(f"  Local minima: {local_min} ({100*local_min/total:.1f}%)")
        print(f"  Saddle points: {saddle} ({100*saddle/total:.1f}%)")
        print(f"  Global max H: known from earlier")

        # How many DISTINCT H values are local maxima?
        max_H_vals = set()
        for bits in range(total):
            A = bits_to_adj(bits, n)
            H = compute_H_dp(A, n)
            is_max = True
            for arc_bit in range(m):
                nbr_bits = bits ^ (1 << arc_bit)
                A_nbr = bits_to_adj(nbr_bits, n)
                H_nbr = compute_H_dp(A_nbr, n)
                if H_nbr > H:
                    is_max = False
                    break
            if is_max:
                max_H_vals.add(H)

        print(f"  H values at local maxima: {sorted(max_H_vals)}")

    print("\n" + "=" * 70)
    print("DONE")
    print("=" * 70)

if __name__ == '__main__':
    main()
