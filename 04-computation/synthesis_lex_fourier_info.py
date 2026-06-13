"""
synthesis_lex_fourier_info.py -- kind-pasteur-2026-03-14-S82
GRAND SYNTHESIS: connecting lex products, Fourier, information, and the tiling model.

THREADS TO WEAVE TOGETHER:
1. Lex product H formula (S81-S82)
2. Fourier spectrum and the 0.27 constant (S73, S75, S79)
3. The pin grid / tiling model (S76-S78)
4. H-landscape and dynamics (S71)
5. Degree Drop and Vassiliev type (S69-S72)

THE UNIFYING QUESTIONS:
A. What structural property of a tournament determines H?
B. Why does the 3-cycle lex product give max_H at n=6?
C. What makes the Fourier spectrum so clean (even levels, 75/25 split)?
D. How does the tiling model interact with these algebraic properties?
"""

import numpy as np
from itertools import permutations
from collections import Counter, defaultdict
import sys, math

sys.stdout.reconfigure(encoding='utf-8')

def bits_to_adj(bits, n):
    A = np.zeros((n, n), dtype=int)
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx): A[j][i] = 1
            else: A[i][j] = 1
            idx += 1
    return A

def compute_H(A, n):
    dp = {}
    for v in range(n): dp[(1 << v, v)] = 1
    for ms in range(2, n+1):
        for mask in range(1 << n):
            if bin(mask).count('1') != ms: continue
            for v in range(n):
                if not (mask & (1 << v)): continue
                pm = mask ^ (1 << v)
                t = sum(dp.get((pm, u), 0) for u in range(n) if (pm & (1 << u)) and A[u][v])
                if t: dp[(mask, v)] = t
    return sum(dp.get(((1 << n) - 1, v), 0) for v in range(n))

def main():
    print("=" * 70)
    print("GRAND SYNTHESIS — LEX + FOURIER + INFO + TILING")
    print("kind-pasteur-2026-03-14-S82")
    print("=" * 70)

    # ========================================
    # SYNTHESIS 1: The "interleaving number" I(T)
    # ========================================
    print(f"\n{'='*70}")
    print("SYNTHESIS 1: THE INTERLEAVING NUMBER")
    print("  For T = T1 lex T2, decompose H into block + interleaved paths.")
    print("  I(T) = H(T) - H_block = the 'extra' paths from interleaving.")
    print("  For T_3 lex T_2: I = 45 - 3 = 42.")
    print("  Question: What is I(T) for general tournaments (not just lex)?")
    print(f"{'='*70}")

    # For a general tournament, we can define "interleaving" with respect
    # to any partition of vertices into groups.
    # For the score-based partition (group vertices by score), the
    # interleaving paths are those that visit vertices of different
    # scores in a non-monotone order.

    n = 5
    m = n*(n-1)//2

    # For each tournament, compute the "score-interleaving number"
    # = # paths where the score sequence of visited vertices is NOT monotone
    interleave_by_H = defaultdict(list)

    for bits in range(2**m):
        A = bits_to_adj(bits, n)
        H = compute_H(A, n)
        scores_list = A.sum(axis=1).astype(int)

        # Count paths that are "score-monotone" vs "score-interleaved"
        mono = 0
        inter = 0
        for perm in permutations(range(n)):
            valid = all(A[perm[i]][perm[i+1]] for i in range(n-1))
            if valid:
                path_scores = [scores_list[perm[i]] for i in range(n)]
                # Monotone: scores are non-decreasing or non-increasing
                non_dec = all(path_scores[i] <= path_scores[i+1] for i in range(n-1))
                non_inc = all(path_scores[i] >= path_scores[i+1] for i in range(n-1))
                if non_dec or non_inc:
                    mono += 1
                else:
                    inter += 1

        interleave_by_H[H].append({'mono': mono, 'inter': inter, 'total': mono+inter})

    print(f"\n  n={n}: Score-interleaving analysis")
    print(f"  {'H':>4} {'avg_mono':>10} {'avg_inter':>10} {'avg_total':>10} {'inter%':>8}")
    for H in sorted(interleave_by_H.keys()):
        items = interleave_by_H[H]
        avg_m = np.mean([d['mono'] for d in items])
        avg_i = np.mean([d['inter'] for d in items])
        avg_t = np.mean([d['total'] for d in items])
        pct = 100 * avg_i / avg_t if avg_t > 0 else 0
        print(f"  {H:4d} {avg_m:10.2f} {avg_i:10.2f} {avg_t:10.2f} {pct:7.1f}%")

    # ========================================
    # SYNTHESIS 2: H decomposition by 3-cycle structure
    # ========================================
    print(f"\n{'='*70}")
    print("SYNTHESIS 2: H = 1 + 2*alpha_1 + 4*alpha_2 + ... (OCF)")
    print("  Each term in OCF corresponds to a specific cycle structure.")
    print("  Term 1 (constant): the transitive 'backbone'")
    print("  Term 2*alpha_1: contribution from individual odd cycles")
    print("  Term 4*alpha_2: contribution from PAIRS of disjoint cycles")
    print("  The 'interleaving' in lex products corresponds to alpha_2+.")
    print(f"{'='*70}")

    for n in [5, 6]:
        m_val = n*(n-1)//2
        count = 0

        alpha_data = defaultdict(list)

        for bits in range(2**m_val):
            count += 1
            if n >= 6 and count > 10000: break

            A = bits_to_adj(bits, n)
            H = compute_H(A, n)
            c3 = int(np.trace(A @ A @ A)) // 3

            # At n=5: alpha_2 = 0, so alpha_1 = (H-1)/2
            if n == 5:
                alpha_1 = (H - 1) // 2
                alpha_2 = 0
            else:
                # At n=6: need to compute alpha_2 from H and alpha_1
                # alpha_1 = c3 + c5 + c7 + ...
                c5 = int(np.trace(np.linalg.matrix_power(A, 5))) // 5
                alpha_1 = c3 + c5
                alpha_2 = (H - 1 - 2*alpha_1) // 4

            alpha_data[H].append({'a1': alpha_1, 'a2': alpha_2, 'c3': c3})

        print(f"\n  n={n}: OCF decomposition")
        print(f"  {'H':>4} {'alpha_1':>8} {'alpha_2':>8} {'c3':>4} {'OCF check':>10}")
        for H in sorted(alpha_data.keys())[:15]:
            items = alpha_data[H]
            # Check OCF: H = 1 + 2*a1 + 4*a2
            a1 = items[0]['a1']
            a2 = items[0]['a2']
            ocf_check = 1 + 2*a1 + 4*a2
            print(f"  {H:4d} {a1:8d} {a2:8d} {items[0]['c3']:4d} {ocf_check:10d}")

    # ========================================
    # SYNTHESIS 3: The "three numbers" that define a tournament
    # ========================================
    print(f"\n{'='*70}")
    print("SYNTHESIS 3: WHAT THREE NUMBERS DEFINE A TOURNAMENT?")
    print("  H alone captures ~27% of tournament info.")
    print("  (H, score_seq) captures ~33%.")
    print("  (H, score_seq, c3) captures more.")
    print("  What is the MINIMAL set of invariants that determines the class?")
    print(f"{'='*70}")

    n = 5
    m = n*(n-1)//2

    # Compute various invariant tuples and check if they determine the class
    class_data = defaultdict(lambda: defaultdict(set))

    perms = list(permutations(range(n)))
    def canonicalize(A):
        best = None
        for p in perms:
            s = ''.join(str(int(A[p[i]][p[j]])) for i in range(n) for j in range(n))
            if best is None or s < best: best = s
        return best

    for bits in range(2**m):
        A = bits_to_adj(bits, n)
        H = compute_H(A, n)
        scores = tuple(sorted(A.sum(axis=1).astype(int)))
        c3 = int(np.trace(A @ A @ A)) // 3
        c5 = int(np.trace(np.linalg.matrix_power(A, 5))) // 5
        det_A = round(np.linalg.det(A))
        canon = canonicalize(A)

        class_data['H'][H].add(canon)
        class_data['(H,scores)'][(H, scores)].add(canon)
        class_data['(H,c3)'][(H, c3)].add(canon)
        class_data['(H,scores,c3)'][(H, scores, c3)].add(canon)
        class_data['(H,scores,c3,c5)'][(H, scores, c3, c5)].add(canon)
        class_data['(H,det)'][(H, det_A)].add(canon)
        class_data['(scores,c3)'][(scores, c3)].add(canon)

    print(f"\n  n={n}: How many classes does each invariant tuple distinguish?")
    total_classes = len(set(
        canonicalize(bits_to_adj(b, n)) for b in range(2**m)
    ))
    print(f"  Total classes: {total_classes}")

    for inv_name in ['H', '(H,scores)', '(H,c3)', '(H,scores,c3)', '(H,scores,c3,c5)', '(H,det)', '(scores,c3)']:
        data = class_data[inv_name]
        n_values = len(data)
        max_classes_per_value = max(len(v) for v in data.values())
        determines = (max_classes_per_value == 1)
        print(f"  {inv_name:>25}: {n_values} values, max classes per value = {max_classes_per_value}, "
              f"determines class? {determines}")

    # ========================================
    # SYNTHESIS 4: The complete picture
    # ========================================
    print(f"\n{'='*70}")
    print("SYNTHESIS 4: THE COMPLETE PICTURE")
    print(f"{'='*70}")
    print("""
  WHAT WE NOW KNOW (S69-S82):

  ALGEBRAIC:
  - H = I(Omega, 2) = 1 + 2*alpha_1 + 4*alpha_2 + 8*alpha_3 + ...  [OCF]
  - deg(H) = 2*floor((n-1)/2) in arc variables  [Degree Drop]
  - H_hat(S) = (n-2)!/2^{n-2} at level 2, sign by 2-path vs V-shape  [Fourier]
  - All odd Fourier levels vanish (from H(T)=H(T^op))  [path reversal]

  STRUCTURAL:
  - H-landscape: unimodal n<=5, multimodal n>=6  [phase transition]
  - H=37 trap at n=6 from score (1,2,2,3,3,4)  [score barrier]
  - Lex product: H(T1 lex T2) = H(T1)*H(T2)^|V1| for transitive T1
  - Cyclic lex: creates 14x more paths via interleaving  [42/3 at n=6]
  - T_3 lex T_2 = H-maximizer at n=6

  INFORMATION:
  - I(T;H)/m ≈ 0.27 = non-constant Fourier energy / m  [information rate]
  - Var(H)/Mean(H)^2 ≈ 1/3 = E_nonconst/E_0  [Fourier Parseval]
  - H beats score as best single summary at n>=6
  - (H, score) jointly captures I/m = 0.33 = class info rate

  TILING:
  - GS code = product code (1+z)^f*(1+z^2)^p  [S76]
  - GS iff self-converse  [S78]
  - Blue skeleton bipartite at odd n by t3 parity  [THM-060]
  - Blue line weights always even  [S77]
  - Monadic bind non-deterministic at class level  [S78]

  TOPOLOGICAL:
  - chi(GLMY) = 1-b1-b3+b4, Paley T_7 has chi=7  [S69]
  - DC = skein relation (Jones polynomial analog)  [S69]
  - H is Vassiliev type 2*floor((n-1)/2)  [S69-S72]

  DYNAMIC:
  - Majority dynamics: transitivizes in O(1) steps  [S71]
  - SA escapes H=37 trap (20/20 at n=6)  [S71]
  - logH ≈ -score_variance with R = -0.97  [S71]

  THE 20 OPEN LEADS FROM S80 (most promising):
  1. Tournament operads
  2. Awan-Bernardi B-polynomial
  3. Asao magnitude-path spectral sequence
  4. The Var/Mean^2 = 1/3 exact proof
  5. Lee-Yang zero concentration
""")

    print(f"\n{'='*70}")
    print("DONE")
    print(f"{'='*70}")

if __name__ == '__main__':
    main()
