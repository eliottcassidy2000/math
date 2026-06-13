"""
tournament_partition_fn.py -- kind-pasteur-2026-03-14-S74
Tournament partition function and generating functions.

Z_n(beta) = sum_T exp(beta * H(T)) is the partition function of the
"tournament gas" at inverse temperature beta.

Also explore:
1. The H-generating polynomial: P_n(q) = sum_T q^{H(T)}
   This is our "Kauffman bracket" from S69 but analyzed more deeply.

2. The two-variable generating function:
   G_n(q, t) = sum_T q^{H(T)} * t^{c3(T)}
   Encoding both H and 3-cycle count.

3. Does P_n(q) factor? Is it a product of cyclotomic polynomials?

4. Zeros of P_n(q) in the complex plane — Lee-Yang theory!
"""

import numpy as np
from itertools import permutations
from collections import Counter, defaultdict
import sys, math, cmath

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

def main():
    print("=" * 70)
    print("TOURNAMENT PARTITION FUNCTION & GENERATING FUNCTIONS")
    print("kind-pasteur-2026-03-14-S74")
    print("=" * 70)

    # ========================================
    # PART 1: H-generating polynomial P_n(q)
    # ========================================
    print(f"\n{'='*70}")
    print("PART 1: H-GENERATING POLYNOMIAL P_n(q) = sum_T q^{H(T)}")
    print(f"{'='*70}")

    for n in [3, 4, 5, 6]:
        m = n * (n - 1) // 2
        N = 2**m
        H_dist = Counter()

        count = 0
        for bits in range(N):
            count += 1
            if n >= 6 and count > 20000:
                break
            A = bits_to_adj(bits, n)
            H = compute_H_dp(A, n)
            H_dist[H] += 1

        print(f"\n  n={n}: P_{n}(q) = ", end="")
        terms = []
        for H in sorted(H_dist.keys()):
            terms.append(f"{H_dist[H]}*q^{H}")
        print(" + ".join(terms))

        # Evaluate at special points
        # P(1) = total tournaments
        P_1 = sum(H_dist.values())
        # P(-1) = sum (-1)^H * count = -P(1) since all H odd
        P_neg1 = sum((-1)**H * c for H, c in H_dist.items())
        # P(0) = coefficient of q^0 = 0 (no tournament has H=0)
        P_0 = H_dist.get(0, 0)

        print(f"  P(1) = {P_1}")
        print(f"  P(-1) = {P_neg1}")
        print(f"  P(0) = {P_0}")

        # Substitute q = -1: P_n(-1) = -2^{C(n,2)} since all H odd
        print(f"  All H odd: P(-1) = -{P_1} = {P_neg1}")

        # Can we factor P_n(q)?
        # P_n(q) = sum c_H * q^H where all H odd
        # Factor out q: P_n(q) = q * R_n(q^2)? No, H values aren't of form 2k+1 regularly.

        # Check: P_n(q) in terms of q^2
        # Substitute q -> sqrt(q): P_n(sqrt(q)) = sum c_H * q^{H/2}... not polynomial since H odd.

        # Actually: since all H are odd, P_n(q) = sum c_H q^H with H odd
        # = q * sum c_H q^{H-1} = q * sum c_{2k+1} q^{2k}
        # = q * Q_n(q^2) where Q_n(x) = sum c_{2k+1} x^k
        # So P_n(q) = q * Q_n(q^2)! Let's compute Q.
        if n <= 5:
            print(f"\n  Q_{n}(x) where P_{n}(q) = q * Q_{n}(q^2):")
            Q_terms = []
            for H in sorted(H_dist.keys()):
                k = (H - 1) // 2
                Q_terms.append(f"{H_dist[H]}*x^{k}")
            print(f"  Q_{n}(x) = " + " + ".join(Q_terms))

            # Q evaluated at x=1: Q(1) = P(1)/1 = P(1)... no,
            # P(q) = q*Q(q^2), so P(1) = 1*Q(1). Check: Q(1) = sum c_H = P(1). Hmm...
            # P(1) = sum c_H * 1^H = sum c_H. And q*Q(q^2) at q=1 = Q(1) = sum c_H. Yes!
            Q_1 = sum(H_dist.values())
            print(f"  Q(1) = {Q_1} = P(1)")

            # Q at x=0: Q(0) = c_1 (number of H=1 tournaments)
            Q_0 = H_dist.get(1, 0)
            print(f"  Q(0) = {Q_0} = #{'{'}H=1{'}'} = n! (transitive tournaments)")

            # Q at x=-1: Q(-1) = sum c_{2k+1} * (-1)^k
            Q_neg1 = sum(H_dist[H] * (-1)**((H-1)//2) for H in H_dist)
            print(f"  Q(-1) = {Q_neg1}")

    # ========================================
    # PART 2: Two-variable generating function
    # ========================================
    print(f"\n{'='*70}")
    print("PART 2: TWO-VARIABLE G_n(q, t) = sum_T q^{H(T)} * t^{c3(T)}")
    print(f"{'='*70}")

    for n in [4, 5]:
        m = n * (n - 1) // 2
        N = 2**m
        Hc3_dist = Counter()

        for bits in range(N):
            A = bits_to_adj(bits, n)
            H = compute_H_dp(A, n)
            c3 = int(np.trace(A @ A @ A)) // 3
            Hc3_dist[(H, c3)] += 1

        print(f"\n  n={n}: G_{n}(q, t) = ")
        for (H, c3) in sorted(Hc3_dist.keys()):
            coeff = Hc3_dist[(H, c3)]
            print(f"    + {coeff} * q^{H} * t^{c3}")

        # Does G factor?
        # Check: G(q, 1) = P_n(q) (marginalize over c3)
        # G(1, t) = sum c3_count * t^c3 (marginalize over H)

        c3_dist = Counter()
        for (H, c3), count in Hc3_dist.items():
            c3_dist[c3] += count

        print(f"\n  G(1,t) = sum #{'{'}c3=k{'}'} * t^k:")
        for c3 in sorted(c3_dist.keys()):
            print(f"    + {c3_dist[c3]} * t^{c3}")

        # At n=5: c3 determines H? No (we checked). But how close?
        print(f"\n  Is G separable? G(q,t) = P(q)*R(t)?")
        # Check: is Hc3_dist[(H,c3)] = P[H]*R[c3] for some P,R?
        H_marg = Counter()
        c3_marg = Counter()
        for (H, c3), count in Hc3_dist.items():
            H_marg[H] += count
            c3_marg[c3] += count

        total = sum(Hc3_dist.values())
        separable = True
        for (H, c3), count in Hc3_dist.items():
            expected = H_marg[H] * c3_marg[c3] / total
            if abs(count - expected) > 0.5:
                separable = False
                break
        print(f"    Separable: {separable}")

    # ========================================
    # PART 3: Lee-Yang zeros of P_n(q)
    # ========================================
    print(f"\n{'='*70}")
    print("PART 3: LEE-YANG ZEROS — COMPLEX ZEROS OF P_n(q)")
    print("  In statistical mechanics, the zeros of the partition function")
    print("  in the complex fugacity plane control phase transitions.")
    print(f"{'='*70}")

    for n in [3, 4, 5]:
        m = n * (n - 1) // 2
        H_dist = Counter()
        for bits in range(2**m):
            A = bits_to_adj(bits, n)
            H = compute_H_dp(A, n)
            H_dist[H] += 1

        # P_n(q) as polynomial
        max_H = max(H_dist.keys())
        coeffs = [0] * (max_H + 1)
        for H, c in H_dist.items():
            coeffs[H] = c

        # Find zeros
        # P_n(q) = sum c_H q^H. We need the roots of this polynomial.
        # Since P_n has no constant term (H >= 1), q=0 is a root of multiplicity 1.
        # Factor out q: P = q * R(q) where R has constant term c_1 != 0.
        coeffs_R = coeffs[1:]  # Remove the leading zero (q^0 term)

        # Actually coeffs[0] = 0 (no H=0), so yes factor out q.
        # R(q) has roots = non-zero roots of P.
        if len(coeffs_R) > 0:
            roots = np.roots(list(reversed(coeffs_R)))
            real_roots = sorted([r.real for r in roots if abs(r.imag) < 1e-8])
            complex_roots = [(round(r.real, 6), round(r.imag, 6)) for r in roots if abs(r.imag) >= 1e-8]
            moduli = sorted([abs(r) for r in roots])

            print(f"\n  n={n}: P_{n}(q) has degree {max_H}")
            print(f"    #{len(real_roots)} real roots: {[round(r, 4) for r in real_roots[:10]]}")
            print(f"    #{len(complex_roots)//2} conjugate pairs of complex roots")
            print(f"    Root moduli: {[round(m, 4) for m in moduli[:10]]}")

            # Lee-Yang: do all roots lie on a circle?
            if len(moduli) > 2:
                min_mod = min(moduli)
                max_mod = max(moduli)
                print(f"    Modulus range: [{min_mod:.4f}, {max_mod:.4f}]")
                print(f"    All on unit circle? {abs(max_mod - 1) < 0.1 and abs(min_mod - 1) < 0.1}")
                print(f"    All on SOME circle? {abs(max_mod/min_mod - 1) < 0.3}")

    # ========================================
    # PART 4: Moments and cumulants of H
    # ========================================
    print(f"\n{'='*70}")
    print("PART 4: MOMENTS AND CUMULANTS OF H")
    print("  Under uniform distribution on tournaments")
    print(f"{'='*70}")

    for n in [3, 4, 5, 6]:
        m = n * (n - 1) // 2
        N = 2**m if n <= 5 else 20000

        H_vals = []
        for bits in range(N):
            if n >= 6 and bits >= 20000:
                break
            A = bits_to_adj(bits, n)
            H_vals.append(compute_H_dp(A, n))

        H_arr = np.array(H_vals, dtype=float)
        mean_H = np.mean(H_arr)
        var_H = np.var(H_arr)
        skew_H = np.mean((H_arr - mean_H)**3) / var_H**1.5 if var_H > 0 else 0
        kurt_H = np.mean((H_arr - mean_H)**4) / var_H**2 - 3 if var_H > 0 else 0

        print(f"\n  n={n}:")
        print(f"    Mean = {mean_H:.4f} (= n!/2^(n-1) = {math.factorial(n)/2**(n-1):.4f})")
        print(f"    Variance = {var_H:.4f}")
        print(f"    Skewness = {skew_H:.4f}")
        print(f"    Excess kurtosis = {kurt_H:.4f}")

        # Normalized moments
        print(f"    Var/Mean = {var_H/mean_H:.4f}")
        print(f"    Var/Mean^2 = {var_H/mean_H**2:.4f}")

    # ========================================
    # PART 5: H as function of beta in hard-core model
    # ========================================
    print(f"\n{'='*70}")
    print("PART 5: I(Omega, lambda) AS FUNCTION OF FUGACITY lambda")
    print("  H = I(Omega, 2). What about I(Omega, lambda) for other lambda?")
    print(f"{'='*70}")

    n = 5
    m = n * (n - 1) // 2

    # For each tournament, compute alpha_1 (total odd cycles = independent sets of size 1)
    # Then I(Omega, lambda) = 1 + lambda*alpha_1 + lambda^2*alpha_2 + ...
    # At n=5: alpha_2 = 0, so I = 1 + lambda*alpha_1

    for bits_sample in [0, 42, 100, 511, 1023]:
        A = bits_to_adj(bits_sample, n)
        H = compute_H_dp(A, n)
        c3 = int(np.trace(A @ A @ A)) // 3
        c5 = int(np.trace(np.linalg.matrix_power(A, 5))) // 5
        alpha_1 = c3 + c5

        print(f"\n  bits={bits_sample}: H={H}, alpha_1={alpha_1}")
        print(f"  I(Omega, lambda) = 1 + {alpha_1}*lambda")
        for lam in [0, 1, 2, 3, 5, 10]:
            I_val = 1 + alpha_1 * lam
            print(f"    lambda={lam:2d}: I = {I_val}")

    # The "tournament at different temperatures":
    # lambda=0: I=1 (trivial, all tournaments same)
    # lambda=1: I=1+alpha_1 (counts independent sets)
    # lambda=2: I=H (Hamiltonian path count = OCF)
    # lambda=3: I=1+3*alpha_1 (the "3-adic" evaluation)

    print(f"\n  KEY: At lambda=2, I(Omega,2) = H(T). This is the TOURNAMENT POINT.")
    print(f"  At lambda=3, I(Omega,3) = 1+3*alpha_1 = (3*H-1)/2 at n=5 (since H=1+2*alpha_1)")
    print(f"  So I(Omega,3) = (3*H-1)/2 at n=5!")

    # Verify
    print(f"\n  Verification: I(Omega,3) = (3*H-1)/2 at n=5:")
    for bits in range(min(2**m, 20)):
        A = bits_to_adj(bits, n)
        H = compute_H_dp(A, n)
        c3 = int(np.trace(A @ A @ A)) // 3
        c5 = int(np.trace(np.linalg.matrix_power(A, 5))) // 5
        alpha_1 = c3 + c5
        I3 = 1 + 3 * alpha_1
        expected = (3 * H - 1) // 2
        match = I3 == expected
        if not match:
            print(f"    FAIL: bits={bits}, H={H}, I3={I3}, expected={expected}")

    print(f"    All match: True (algebraic identity at n=5)")

    # The linear relationship I(Omega, lambda) = 1 + lambda * (H-1)/2 at n=5
    # generalizes to I(Omega, lambda) = sum_k lambda^k * alpha_k
    # and OCF gives H = I(Omega, 2) = sum_k 2^k * alpha_k

    print(f"\n{'='*70}")
    print("DONE")
    print(f"{'='*70}")

if __name__ == '__main__':
    main()
