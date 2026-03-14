"""
signed_redei_proof.py -- kind-pasteur-2026-03-14-S70
Prove F(T, -1) is always odd for any tournament.

F(T, x) = sum over Ham paths P: x^{fwd(P)}
F(T, -1) = sum over Ham paths P: (-1)^{fwd(P)}

Redei: F(T, 1) = H(T) is always odd (classical).
"Signed Redei": F(T, -1) is also always odd.

Approach 1: Direct from F(T, x) = sum over paths.
  F(T, 1) = H(T) and F(T, -1) are both evaluations of a polynomial.
  F(T, 1) - F(T, -1) = 2 * sum_{fwd odd} count = even.
  So F(T, -1) = F(T, 1) - 2 * (#paths with odd fwd) = H - 2k.
  H is odd, 2k is even, so F(-1) is odd. QED!

Wait, that proves it immediately! F(-1) = H - 2*(#paths with odd fwd) is odd
because H is odd (Redei) and 2*(anything) is even.

Actually: F(1) = H = total paths (odd by Redei).
F(-1) = (even_fwd_paths) - (odd_fwd_paths).
F(1) = (even_fwd_paths) + (odd_fwd_paths).
F(1) - F(-1) = 2 * odd_fwd_paths.
So F(-1) = F(1) - 2 * odd_fwd_paths = odd - even = odd. QED!

This is trivial! The "signed Redei" follows directly from classical Redei.
But let's explore deeper: what about F(T, omega) for roots of unity?

Also investigate:
- F(T, i) where i = sqrt(-1): what is |F(T, i)|^2?
- F(T, omega_3) where omega_3 = e^{2pi*i/3}: does 9 | |F(omega_3)|^2?
- F(T, -1) mod 4: what determines the residue?
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

def fwd_polynomial_full(A, n):
    """Compute F(T,x) as dict fwd_count -> number_of_paths."""
    poly = Counter()
    for perm in permutations(range(n)):
        valid = True
        for i in range(n-1):
            if A[perm[i]][perm[i+1]] != 1:
                valid = False
                break
        if valid:
            fwd = sum(1 for i in range(n-1) if perm[i] < perm[i+1])
            poly[fwd] += 1
    return poly

def evaluate_F(poly, x):
    """Evaluate F(T,x) at a specific x value."""
    return sum(x**k * c for k, c in poly.items())

def main():
    print("=" * 70)
    print("SIGNED REDEI THEOREM & F(T, roots of unity)")
    print("kind-pasteur-2026-03-14-S70")
    print("=" * 70)

    # ====================================
    # PART 1: Proof that F(-1) is always odd
    # ====================================
    print("\n" + "=" * 70)
    print("PART 1: PROOF THAT F(T, -1) IS ALWAYS ODD")
    print("=" * 70)
    print("""
  THEOREM: F(T, -1) is odd for any tournament T on n >= 1 vertices.

  PROOF:
  Let e = #{paths with even fwd}, o = #{paths with odd fwd}.
  Then F(1) = e + o = H(T) and F(-1) = e - o.
  By Redei's theorem, H(T) = e + o is ODD.
  Therefore F(-1) = (e + o) - 2o = H - 2o is odd - even = ODD. QED.

  Note: this works because the parity of F(-1) equals the parity of H,
  which is always odd. The SIGNED Redei theorem is a TRIVIAL consequence
  of the UNSIGNED Redei theorem. No new proof needed!
""")

    # Verify computationally
    for n in [3, 4, 5, 6]:
        print(f"  Verification at n={n}:", end=" ")
        all_odd = True
        count = 0
        for bits in range(2**(n*(n-1)//2)):
            count += 1
            if n >= 6 and count > 5000:
                break
            A = bits_to_adj(bits, n)
            poly = fwd_polynomial_full(A, n)
            F_neg1 = evaluate_F(poly, -1)
            if F_neg1 % 2 == 0:
                all_odd = False
                break
        print(f"{'ALL ODD' if all_odd else 'FAILURE!'} ({count} checked)")

    # ====================================
    # PART 2: F(-1) mod 4 structure
    # ====================================
    print("\n" + "=" * 70)
    print("PART 2: F(-1) MOD 4 STRUCTURE")
    print("  F(-1) = H - 2o where o = #odd-fwd paths")
    print("  F(-1) mod 4 = H mod 4 - 2(o mod 2)")
    print("  Since H mod 4 = 1 + 2*alpha_1 mod 4:")
    print("  F(-1) mod 4 depends on alpha_1 parity AND o parity")
    print("=" * 70)

    for n in [4, 5, 6]:
        print(f"\n--- n = {n} ---")
        mod4_dist = Counter()
        mod4_by_H = defaultdict(Counter)

        count = 0
        for bits in range(2**(n*(n-1)//2)):
            count += 1
            if n >= 6 and count > 5000:
                break

            A = bits_to_adj(bits, n)
            poly = fwd_polynomial_full(A, n)
            H = sum(poly.values())
            F_neg1 = int(evaluate_F(poly, -1))

            mod4_dist[F_neg1 % 4] += 1
            mod4_by_H[H].update([F_neg1 % 4])

        print(f"  F(-1) mod 4 distribution: {dict(sorted(mod4_dist.items()))}")

        print(f"  F(-1) mod 4 by H:")
        for H in sorted(mod4_by_H.keys())[:10]:
            print(f"    H={H:3d}: {dict(sorted(mod4_by_H[H].items()))}")

    # ====================================
    # PART 3: F(T, i) where i = sqrt(-1)
    # ====================================
    print("\n" + "=" * 70)
    print("PART 3: F(T, i) — GAUSSIAN EVALUATION")
    print("  F(T, i) = sum paths P: i^{fwd(P)}")
    print("  |F(T, i)|^2 is a tournament invariant")
    print("=" * 70)

    for n in [4, 5]:
        print(f"\n--- n = {n} ---")
        norm_sq_dist = Counter()
        norm_sq_by_H = defaultdict(list)

        for bits in range(2**(n*(n-1)//2)):
            A = bits_to_adj(bits, n)
            poly = fwd_polynomial_full(A, n)
            H = sum(poly.values())

            F_i = evaluate_F(poly, 1j)
            norm_sq = round(abs(F_i)**2)
            norm_sq_dist[norm_sq] += 1
            norm_sq_by_H[H].append(norm_sq)

        print(f"  |F(i)|^2 distribution: {dict(sorted(norm_sq_dist.items()))}")
        print(f"\n  |F(i)|^2 by H:")
        for H in sorted(norm_sq_by_H.keys()):
            vals = sorted(set(norm_sq_by_H[H]))
            print(f"    H={H:3d}: |F(i)|^2 in {vals}")

    # ====================================
    # PART 4: F(T, omega_3) — cube root evaluation
    # ====================================
    print("\n" + "=" * 70)
    print("PART 4: F(T, omega_3) — CUBE ROOT OF UNITY")
    print("  omega = e^{2*pi*i/3}")
    print("  F(T, omega) relates to THM-085 (9 | F(T, omega) for n>=6)")
    print("=" * 70)

    omega = np.exp(2j * np.pi / 3)

    for n in [4, 5, 6]:
        print(f"\n--- n = {n} ---")
        norm_dist = Counter()

        count = 0
        for bits in range(2**(n*(n-1)//2)):
            count += 1
            if n >= 6 and count > 3000:
                break
            A = bits_to_adj(bits, n)
            poly = fwd_polynomial_full(A, n)

            F_omega = evaluate_F(poly, omega)
            norm_sq = round(abs(F_omega)**2)
            norm_dist[norm_sq] += 1

        print(f"  |F(omega_3)|^2 values: {sorted(norm_dist.keys())[:20]}")
        print(f"  All divisible by 9: {all(v % 9 == 0 for v in norm_dist.keys() if v > 0)}")
        print(f"  Min nonzero: {min(v for v in norm_dist.keys() if v > 0) if any(v > 0 for v in norm_dist.keys()) else 'N/A'}")

    # ====================================
    # PART 5: DEEP — F(-1) vs cycle structure
    # ====================================
    print("\n" + "=" * 70)
    print("PART 5: F(-1) vs CYCLE STRUCTURE")
    print("  F(-1) = H - 2o. What determines o (odd-fwd path count)?")
    print("  Can we express o in terms of cycle counts?")
    print("=" * 70)

    for n in [4, 5]:
        print(f"\n--- n = {n} ---")
        # Collect (H, F(-1), c3, o_fwd)
        data = []

        for bits in range(2**(n*(n-1)//2)):
            A = bits_to_adj(bits, n)
            poly = fwd_polynomial_full(A, n)
            H = sum(poly.values())
            F_neg1 = int(evaluate_F(poly, -1))
            o_fwd = sum(c for k, c in poly.items() if k % 2 == 1)
            c3 = int(np.trace(A @ A @ A)) // 3

            data.append({
                'H': H, 'F(-1)': F_neg1, 'o': o_fwd,
                'c3': c3, 'bits': bits
            })

        # Check: is o determined by (H, c3)?
        o_by_H_c3 = defaultdict(set)
        for d in data:
            o_by_H_c3[(d['H'], d['c3'])].add(d['o'])

        all_determined = all(len(v) == 1 for v in o_by_H_c3.values())
        print(f"  o determined by (H, c3)? {all_determined}")

        if not all_determined:
            for key, vals in sorted(o_by_H_c3.items()):
                if len(vals) > 1:
                    print(f"    (H={key[0]}, c3={key[1]}): o in {sorted(vals)}")

        # Formula attempt: F(-1) = H - 2o, and H = 1 + 2*alpha_1 + 4*alpha_2
        # Is there a pattern: o = f(alpha_1, alpha_2)?
        # At n=4: alpha_2 = 0 always, so o depends only on c3 (= alpha_1)
        # At n=5: alpha_2 = 0 (trivially, need >=6 vertices for disjoint cycles)
        #   So o should depend on alpha_1 = c3 + c5

        if n == 5:
            # Check: o vs writhe
            writhe_o = defaultdict(set)
            for d in data:
                A = bits_to_adj(d['bits'], n)
                asc = sum(A[i][j] for i in range(n) for j in range(i+1, n))
                writhe = 2*asc - n*(n-1)//2
                writhe_o[(d['H'], abs(writhe))].add(d['o'])

            print(f"\n  o determined by (H, |writhe|)? {all(len(v)==1 for v in writhe_o.values())}")

    print("\n" + "=" * 70)
    print("DONE")
    print("=" * 70)

if __name__ == '__main__':
    main()
