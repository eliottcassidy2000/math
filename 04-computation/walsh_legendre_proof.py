#!/usr/bin/env python3
"""
walsh_legendre_proof.py -- Towards proving the Walsh-Legendre Sign Law

THEOREM (empirical, p=7,11): For circulant tournaments on Z_p with p=3 mod 4:
  sign(h_hat[S]) = chi(prod_{i in S}(i+1))
where chi is the Legendre symbol mod p, and h_hat[S] is the Walsh coefficient
of H(sigma) on the orientation cube {+1,-1}^m (m = (p-1)/2).

PROOF STRATEGY:
The Walsh coefficient is:
  h_hat[S] = (1/2^m) sum_{sigma} H(sigma) * prod_{i in S} sigma_i

By circulant symmetry, H(sigma) depends on sigma through the connection set S.
The connection set is:
  S(sigma) = {(i+1) : sigma_i = +1} union {(p-i-1) : sigma_i = -1}

Note: sigma_i = +1 iff gap (i+1) is in S. So sigma_i = chi(g_i) where
g_i = actual gap used from chord pair i (either i+1 or p-i-1).

The H-path count for a circulant tournament on Z_p with connection set S is:
  H = (sum over Hamilton paths) 1

By the transfer matrix method (Tang-Yau, THM-125), for circulant tournaments
with p prime, H decomposes into eigenspace contributions:
  H/p = sum_k h_k
where h_k is the contribution from eigenspace k of Z_p.

The eigenspace k is the 1-dimensional subspace of the DFT at frequency k.
For k=0: h_0 depends only on the CONNECTION SET'S PROPERTIES (size, etc.)
For k!=0: h_k depends on the Fourier transform of the indicator of S.

KEY INSIGHT: S_hat(k) = sum_{s in S} omega^{ks} where omega = exp(2pi*i/p).
For the Interval S = {1,...,m}: S_hat(k) is a geometric sum.
For the Paley S = QR: S_hat(k) = Gauss sum (related to chi(k)).

The Walsh-Legendre sign law says that when we Fourier-analyze H over the
orientation cube (sigma space), the sign structure is determined by chi.

This script investigates:
1. The relationship between H and the Fourier spectrum of S
2. How changing sigma affects S_hat(k)
3. Whether the sign law follows from Gauss sum identities

Author: kind-pasteur-2026-03-12-S59c
"""

import cmath
from math import comb
from collections import defaultdict


def build_adj(p, S):
    S_set = set(S)
    A = [[0]*p for _ in range(p)]
    for i in range(p):
        for s in S_set:
            A[i][(i + s) % p] = 1
    return A


def count_ham_paths(A, n):
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            if dp[mask][v] == 0:
                continue
            for w in range(n):
                if mask & (1 << w):
                    continue
                if A[v][w]:
                    dp[mask | (1 << w)][w] += dp[mask][v]
    full = (1 << n) - 1
    return sum(dp[full][v] for v in range(n))


def fourier_spectrum(S, p):
    """Compute the Fourier transform of the indicator of S."""
    omega = cmath.exp(2j * cmath.pi / p)
    S_hat = []
    for k in range(p):
        val = sum(omega ** (k * s) for s in S)
        S_hat.append(val)
    return S_hat


def sigma_to_S(sigma, p):
    m = (p - 1) // 2
    S = []
    for i in range(m):
        if sigma[i] == 1:
            S.append(i + 1)
        else:
            S.append(p - i - 1)
    return sorted(S)


def main():
    print("=" * 70)
    print("TOWARDS PROVING THE WALSH-LEGENDRE SIGN LAW")
    print("=" * 70)

    for p in [7, 11]:
        m = (p - 1) // 2
        omega = cmath.exp(2j * cmath.pi / p)
        QR = set(j for j in range(1, p) if pow(j, (p-1)//2, p) == 1)

        print(f"\n{'='*70}")
        print(f"p={p}, m={m}")
        print(f"{'='*70}")

        # 1. Fourier spectrum for each sigma
        print(f"\n  1. FOURIER SPECTRUM OF S(sigma)")

        pairs = [(s, p - s) for s in range(1, m + 1)]
        H_dict = {}
        S_hat_dict = {}

        for bits in range(1 << m):
            sigma = tuple((1 if bits & (1 << i) else -1) for i in range(m))
            S = sigma_to_S(sigma, p)
            A = build_adj(p, S)
            H = count_ham_paths(A, p)
            H_dict[sigma] = H
            S_hat = fourier_spectrum(S, p)
            S_hat_dict[sigma] = S_hat

        # 2. The key quantity: |S_hat(k)|^2
        # For a circulant tournament, the transfer matrix eigenvalue at frequency k
        # depends on S_hat(k).
        print(f"\n  2. |S_hat(k)|^2 FOR EACH sigma")
        print(f"    (each row: sigma -> |S_hat(1)|^2, ..., |S_hat(m)|^2)")

        for sigma in sorted(H_dict, key=lambda s: -H_dict[s]):
            S_hat = S_hat_dict[sigma]
            H = H_dict[sigma]
            powers = [abs(S_hat[k])**2 for k in range(1, m+1)]
            powers_str = " ".join(f"{x:>6.1f}" for x in powers)
            print(f"    sigma={sigma}: H={H:>8}, |S_hat|^2 = [{powers_str}]")

        # 3. Connection to Walsh: how does sigma_i affect S_hat?
        print(f"\n  3. EFFECT OF FLIPPING sigma_i ON S_hat(k)")
        # When sigma_i flips (+1 to -1), gap (i+1) is replaced by gap (p-i-1).
        # Effect on S_hat(k):
        #   Delta S_hat(k) = omega^{k*(p-i-1)} - omega^{k*(i+1)}
        #                  = omega^{-k*(i+1)} - omega^{k*(i+1)}  [since omega^{kp}=1]
        #                  = -2i*sin(2*pi*k*(i+1)/p)

        print(f"\n    For chord i (gap g=i+1), flipping sigma_i changes S_hat(k) by:")
        print(f"    Delta(k) = omega^{{-kg}} - omega^{{kg}} = -2i*sin(2*pi*k*g/p)")
        print()
        for i in range(m):
            g = i + 1
            print(f"    chord {i} (g={g}):", end="")
            for k in range(1, m+1):
                delta = omega**(-k*g) - omega**(k*g)
                print(f"  k={k}: {delta:.3f}", end="")
            print()

        # 4. H as a function of S_hat
        # The transfer matrix T for a circulant on Z_p with eigenvalue decomposition:
        #   T = sum_k lambda_k P_k
        # where P_k is the projection onto the k-th DFT eigenspace.
        # For an m x m transfer matrix at frequency k, the eigenvalue is S_hat(k).
        # But it's more complex: the transfer matrix has SIZE (p-1) x (p-1),
        # with eigenvalues related to S_hat.

        # Actually, by THM-125, for circulant tournaments the "symbol matrix"
        # is constant. The eigenvalue decomposition gives:
        #   H/p = sum_k det(I + M_k)  (or similar)
        # where M_k depends on S_hat(k).

        # For NOW, let's just check: is H a symmetric function of |S_hat(k)|^2?
        print(f"\n  4. IS H DETERMINED BY THE |S_hat(k)|^2 SPECTRUM?")

        # Group sigmas by their |S_hat|^2 spectrum
        spectrum_groups = defaultdict(list)
        for sigma, H in H_dict.items():
            S_hat = S_hat_dict[sigma]
            spectrum = tuple(round(abs(S_hat[k])**2, 4) for k in range(1, p))
            spectrum_groups[spectrum].append((sigma, H))

        for spec, items in sorted(spectrum_groups.items(), key=lambda x: -x[1][0][1]):
            H_vals = sorted(set(H for _, H in items))
            sigmas = [s for s, _ in items]
            print(f"    spectrum={spec[:m]}: H={H_vals}, {len(sigmas)} orientations")

        # 5. Does |S_hat(k)|^2 have a formula in terms of sigma?
        print(f"\n  5. |S_hat(k)|^2 FORMULA")
        # |S_hat(k)|^2 = sum_{s1 in S} sum_{s2 in S} omega^{k(s1-s2)}
        #              = m + sum_{s1 != s2} omega^{k(s1-s2)}
        #              = m + sum_{d != 0} N(d) * omega^{kd}
        # where N(d) = |{(s1,s2) in S^2 : s1 - s2 = d}| = auto-correlation of S at shift d.

        # For connection set S = {g_1,...,g_m}:
        # N(d) = |{(i,j) : g_i - g_j = d mod p}|

        # When sigma changes, S changes, and so does N(d).

        # The auto-correlation N(d) at gap d:
        # For Interval: N(d) depends on how many pairs (g_i, g_j) with g_i - g_j = d
        # For Paley: N(d) is constant (difference set property)

        print(f"\n    Auto-correlation N(d) for a few orientations:")
        for sigma in [tuple([1]*m), tuple([-1]*m)]:
            S = sigma_to_S(sigma, p)
            N = defaultdict(int)
            for s1 in S:
                for s2 in S:
                    d = (s1 - s2) % p
                    N[d] += 1
            print(f"    sigma={sigma}, S={S}")
            print(f"      N(d) = {dict(sorted(N.items()))}")

        # 6. Walsh coefficient in terms of N(d)
        print(f"\n  6. CONNECTING WALSH TO AUTO-CORRELATION")
        # H(sigma) = f(S(sigma)) = f(auto-corr of S(sigma)) approximately
        # The Walsh coefficient h_hat[{i}] measures the linear effect of sigma_i on H.
        # It's zero by path reversal symmetry.
        # h_hat[{i,j}] measures the effect of the PAIR interaction between chords i and j.

        # When we flip chord i: gap (i+1) changes to gap (p-i-1).
        # The auto-correlation N changes. The change in N(d) for each d:
        #   For all other gaps g_j in S (j != i):
        #     Old contributions: omega^{k*((i+1)-g_j)} + omega^{k*(g_j-(i+1))}
        #     New contributions: omega^{k*((p-i-1)-g_j)} + omega^{k*(g_j-(p-i-1))}
        #     Difference: replace (i+1) with (p-i-1) = -(i+1) mod p
        #     So all differences involving chord i change sign!

        # This is the key: N(d) -> N(d) with sign flips on differences involving chord i.
        # The effect on |S_hat(k)|^2 depends on how much chord i contributes to the auto-correlation.

        # 7. The Gauss sum connection
        print(f"\n  7. GAUSS SUM CONNECTION")
        # For Paley: S = QR_p. The Gauss sum is:
        #   tau(chi) = sum_{a=1}^{p-1} chi(a) omega^a
        # Key identity: |tau(chi)|^2 = p.
        # And: sum_{a in QR} omega^{ka} = (tau(chi) * chi(k) - 1) / 2
        #   (for k != 0)

        # So S_hat_Paley(k) = (tau * chi(k) - 1) / 2

        # |S_hat_Paley(k)|^2 = |tau * chi(k) - 1|^2 / 4
        # = (|tau|^2 - tau*chi(k) - tau_bar*chi(k) + 1) / 4
        # = (p + 1 - 2*Re(tau*chi(k))) / 4

        # For chi(k) = +1: (p + 1 - 2*Re(tau)) / 4
        # For chi(k) = -1: (p + 1 + 2*Re(tau)) / 4

        # So |S_hat_Paley|^2 takes EXACTLY TWO values depending on chi(k)!

        tau = sum((1 if j in QR else -1) * omega**j for j in range(1, p))
        print(f"    Gauss sum tau = {tau:.6f}")
        print(f"    |tau|^2 = {abs(tau)**2:.2f} (should be {p})")
        print(f"    Re(tau) = {tau.real:.6f}")

        if p % 4 == 3:
            # At p=3 mod 4: tau = i*sqrt(p) (purely imaginary)
            print(f"    p=3 mod 4: tau should be i*sqrt(p) = {1j * p**0.5:.6f}")
            print(f"    Re(tau) = {tau.real:.6f} (should be ~0)")

            # Therefore |S_hat_Paley(k)|^2 = (p + 1) / 4 for ALL k != 0
            # (since Re(tau) = 0)
            paley_sigma = tuple(1 if (i+1) in QR else -1 for i in range(m))
            S_paley = sigma_to_S(paley_sigma, p)
            S_hat_paley = fourier_spectrum(S_paley, p)

            print(f"\n    Paley |S_hat(k)|^2:")
            for k in range(1, p):
                val = abs(S_hat_paley[k])**2
                predicted = (p + 1) / 4
                print(f"      k={k}: {val:.4f} (predicted: {predicted:.4f})")

        # 8. The Walsh sign from the perspective of the transfer matrix
        print(f"\n  8. TRANSFER MATRIX EIGENVALUE ANALYSIS")
        # For a circulant tournament, the transfer matrix has eigenvalues
        # related to S_hat(k). The H count is:
        #   H = p * prod_{k=1}^{m} (something involving S_hat(k))
        # or more precisely,
        #   H/p = sum_k |contribution from eigenspace k|
        # where the contribution depends on S_hat(k).

        # By THM-125, the symbol matrix is constant for circulant tournaments.
        # This means the eigenspace contribution h_k depends ONLY on S_hat(k).

        # Since |S_hat(k)|^2 is the same for all nonzero k at the Paley point
        # (it's (p+1)/4), the eigenspace contributions are all equal.
        # H_Paley / p = p * h_single  (each of p-1 nonzero eigenspaces contributes h_single)
        # Actually wait, the eigenspaces pair up: k and p-k give conjugate contributions.

        # Let me just compute H/p for each orientation and see the structure
        print(f"\n    H/p for each orientation:")
        for sigma in sorted(H_dict, key=lambda s: -H_dict[s])[:5]:
            H = H_dict[sigma]
            S = sigma_to_S(sigma, p)
            S_hat = fourier_spectrum(S, p)
            shat_sq = tuple(round(abs(S_hat[k])**2, 2) for k in range(1, m+1))
            label = ""
            paley_sigma = tuple(1 if (i+1) in QR else -1 for i in range(m))
            if sigma == paley_sigma or sigma == tuple(-s for s in paley_sigma):
                label = " [PALEY]"
            elif sigma == tuple([1]*m):
                label = " [INTERVAL]"
            print(f"    sigma={sigma}: H/p={H/p:.2f}, |S_hat|^2={shat_sq}{label}")


if __name__ == '__main__':
    main()
