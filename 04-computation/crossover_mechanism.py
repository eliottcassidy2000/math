"""
crossover_mechanism.py — Why Paley loses to Interval at large p

The trace alternation theorem (kind-pasteur-S56c) shows:
  - C_k(Paley) > C_k(Interval) when k ≡ 1 mod 4
  - C_k(Paley) < C_k(Interval) when k ≡ 3 mod 4, k > 3
  - C_3 is always a tie (universal for all circulants on Z_p)

The OCF says H(T) = I(Ω(T), 2) where Ω(T) is the odd-cycle complex.
The cycle counts C_k feed into the independence polynomial through
the graph structure of Ω.

Key question: How does the OCF's 2^|S| weighting interact with
the trace alternation to create the crossover?

Hypothesis: The crossover happens when the CUMULATIVE advantage of
interval at k≡3(4) lengths, amplified by 2^k, exceeds Paley's
advantage at k≡1(4) lengths.

Also investigates: dihedral group interpretation of the interval tournament.

Author: opus-2026-03-12-S60
"""
import sys
import time
import numpy as np
from collections import defaultdict
from fractions import Fraction
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)


def hamiltonian_paths_dp(A, n):
    dp = defaultdict(lambda: defaultdict(int))
    for v in range(n):
        dp[1 << v][v] = 1
    full = (1 << n) - 1
    for mask in range(1, full + 1):
        if not dp[mask]:
            continue
        for v in dp[mask]:
            if dp[mask][v] == 0:
                continue
            for w in range(n):
                if mask & (1 << w):
                    continue
                if A[v][w]:
                    dp[mask | (1 << w)][w] += dp[mask][v]
    return sum(dp[full][v] for v in range(n))


def circulant_adj(n, S):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for s in S:
            A[i][(i+s)%n] = 1
    return A


def paley_set(p):
    return frozenset(pow(x, 2, p) for x in range(1, p))


def interval_set(p):
    m = (p - 1) // 2
    return frozenset(range(1, m + 1))


def circulant_eigenvalues(n, S):
    omega = np.exp(2j * np.pi / n)
    return [sum(omega ** (k * s) for s in S) for k in range(n)]


def trace_power(evals, k):
    """tr(A^k) from eigenvalues."""
    return sum(e**k for e in evals).real


def cycle_count_from_trace(traces, n, max_k):
    """Extract cycle counts C_k from tr(A^k) using Newton's identity.

    tr(A^k) = sum_{d | k, d odd} d * C_d  (for tournaments, only odd d contribute)

    Actually for directed graphs: tr(A^k) counts walks of length k returning to start.
    For tournaments specifically: tr(A^k) = sum_{ℓ | k} ℓ * C_ℓ where C_ℓ = # directed ℓ-cycles.
    """
    C = {}
    for k in range(3, max_k + 1, 2):
        C[k] = traces[k] / k
        for d in range(3, k, 2):
            if k % d == 0:
                C[k] -= d * C[d] / k
    return C


def main():
    print("CROSSOVER MECHANISM: WHY INTERVAL BEATS PALEY AT LARGE p")
    print("=" * 75)

    # Part 1: Trace alternation with explicit Δ values
    print("\nPART 1: TRACE DIFFERENCES Δ_k = C_k(Paley) - C_k(Interval)")
    print("-" * 75)

    for p in [7, 11, 19, 23]:
        if p > 13 and p not in [19, 23]:
            continue
        m = (p - 1) // 2
        S_paley = paley_set(p)
        S_interval = interval_set(p)

        evals_P = circulant_eigenvalues(p, S_paley)
        evals_I = circulant_eigenvalues(p, S_interval)

        print(f"\n  p = {p}, m = {m}")
        print(f"  {'k':>4} {'k mod 4':>8} {'C_k(P)':>16} {'C_k(I)':>16} {'Δ_k':>16} {'|Δ_k|*2^k':>20}")

        cumulative_advantage_P = 0.0
        cumulative_advantage_I = 0.0

        for k in range(3, min(p + 1, 20), 2):
            tr_P = trace_power(evals_P, k)
            tr_I = trace_power(evals_I, k)
            C_P = tr_P / k
            C_I = tr_I / k
            delta = C_P - C_I
            weighted = abs(delta) * (2 ** k)

            if delta > 0.5:
                cumulative_advantage_P += delta * (2 ** k)
                marker = "P"
            elif delta < -0.5:
                cumulative_advantage_I += abs(delta) * (2 ** k)
                marker = "I"
            else:
                marker = "="

            print(f"  {k:>4} {k%4:>8} {C_P:>16.1f} {C_I:>16.1f} {delta:>+16.1f} {weighted:>20.0f} {marker}")

        print(f"  Cumulative 2^k-weighted: Paley advantage = {cumulative_advantage_P:.0f}")
        print(f"  Cumulative 2^k-weighted: Interval advantage = {cumulative_advantage_I:.0f}")
        print(f"  Net: {'PALEY' if cumulative_advantage_P > cumulative_advantage_I else 'INTERVAL'} wins by {abs(cumulative_advantage_P - cumulative_advantage_I):.0f}")

    # Part 2: The Dirichlet kernel explanation
    print(f"\n{'='*75}")
    print("PART 2: DIRICHLET KERNEL vs GAUSS SUM EIGENVALUES")
    print("-" * 75)

    for p in [7, 11, 19]:
        m = (p - 1) // 2
        S_paley = paley_set(p)
        S_interval = interval_set(p)

        evals_P = circulant_eigenvalues(p, S_paley)
        evals_I = circulant_eigenvalues(p, S_interval)

        # Paley: |λ_k| = √((p+1)/4) for all k≥1
        paley_mag = np.sqrt((p + 1) / 4)
        # Interval: |λ_k| = |sin(mπk/p) / sin(πk/p)| (Dirichlet kernel)

        print(f"\n  p = {p}")
        print(f"  Paley: all |λ_k| = {paley_mag:.4f}")
        print(f"  Interval eigenvalue magnitudes:")

        mags_I = []
        for k in range(1, m + 1):
            mag = abs(evals_I[k])
            mags_I.append(mag)

        # Sort by magnitude
        mags_sorted = sorted(mags_I, reverse=True)
        print(f"    Sorted: {[f'{m:.3f}' for m in mags_sorted]}")
        print(f"    Ratio max/min: {mags_sorted[0]/mags_sorted[-1]:.2f}")
        print(f"    Dominant eigenvalue |λ_1| = {abs(evals_I[1]):.4f} ≈ p/π = {p/np.pi:.4f}")

        # The key: how does eigenvalue concentration affect cycle counts?
        # C_k ∝ Σ|λ_j|^k cos(k·phase_j)
        # For Paley: all mags equal, phases from Gauss sum
        # For Interval: one dominant mag, rest small
        # At large k: dominant eigenvalue of interval^k >> Paley flat^k when max_I > flat

        ratio = abs(evals_I[1]) / paley_mag
        print(f"    |λ_1(I)| / |λ(P)| = {ratio:.4f}")
        print(f"    At k=p: dominant^p / flat^p ≈ {ratio**p:.2e}")

    # Part 3: Dihedral group structure
    print(f"\n{'='*75}")
    print("PART 3: DIHEDRAL GROUP INTERPRETATION")
    print("-" * 75)

    for p in [7, 11, 19]:
        m = (p - 1) // 2
        print(f"\n  p = {p}, D_{2*p} = {{rotations}} ∪ {{reflections}}")

        # Automorphism group of interval tournament
        # Rotations: always automorphisms (circulant)
        # Reflection σ: i → -i mod p sends S={1,...,m} to {-m,...,-1}={p-m,...,p-1}
        # This is the COMPLEMENT set, not S itself
        # So σ is NOT an automorphism of interval, but maps T to T̄ (same H)

        S_int = set(range(1, m + 1))
        S_comp = set(range(m + 1, p))
        print(f"    Interval S = {sorted(S_int)}")
        print(f"    σ(S) = -S mod p = {sorted(S_comp)}")
        print(f"    σ is automorphism? {S_int == S_comp}")

        # For Paley: σ(QR) = -QR
        S_paley = set(paley_set(p))
        S_paley_neg = set((p - s) % p for s in S_paley) - {0}
        print(f"    Paley S = {sorted(S_paley)}")
        print(f"    σ(S) = -S mod p = {sorted(S_paley_neg)}")
        print(f"    σ is automorphism of Paley? {S_paley == S_paley_neg}")

        # Multiplication automorphisms
        # For Paley: k*QR = QR for all k ∈ QR → |Aut| ≥ p * m = p(p-1)/2
        # For Interval: k*{1,...,m} = {k, 2k,...,mk} mod p.
        # This equals {1,...,m} only if k=1.
        aut_int = 0
        for k in range(1, p):
            kS = frozenset((k * s) % p for s in S_int)
            if kS == frozenset(S_int):
                aut_int += 1

        aut_paley = 0
        for k in range(1, p):
            kS = frozenset((k * s) % p for s in S_paley)
            if kS == frozenset(S_paley):
                aut_paley += 1

        print(f"    Multiplicative automorphisms: Interval={aut_int}, Paley={aut_paley}")
        print(f"    Full circulant |Aut|: Interval={p*aut_int}, Paley={p*aut_paley}")

    # Part 4: The key question — is interval the global max at p=19?
    # We can't do full exhaustive search, but we can check ALL 512 circulants
    # Problem: each takes ~30s at p=19. Let's check a strategic subset.

    # Part 4: At feasible sizes (p=7, 11), check interval vs all circulants
    print(f"\n{'='*75}")
    print("PART 4: INTERVAL vs ALL CIRCULANTS (exhaustive at p=7, 11)")
    print("-" * 75)

    for p in [7, 11]:
        m = (p - 1) // 2
        pairs = [(d, p - d) for d in range(1, m + 1)]

        all_H = []
        paley_H = None
        interval_H = None
        S_paley = paley_set(p)
        S_int = interval_set(p)

        for bits in range(1 << m):
            S = set()
            for i in range(m):
                if bits & (1 << i):
                    S.add(pairs[i][0])
                else:
                    S.add(pairs[i][1])
            S = frozenset(S)
            A = circulant_adj(p, S)
            H = hamiltonian_paths_dp(A, p)
            all_H.append((H, S))
            if S == S_paley or S == frozenset(set(range(1, p)) - S_paley):
                paley_H = H
            if S == S_int or S == frozenset(set(range(1, p)) - S_int):
                interval_H = H

        all_H.sort(reverse=True)
        print(f"\n  p = {p}: {2**m} circulants")
        print(f"    Top 5 H values:")
        seen_H = set()
        count = 0
        for H, S in all_H:
            if H not in seen_H:
                seen_H.add(H)
                is_paley = (S == S_paley or S == frozenset(set(range(1, p)) - S_paley))
                is_int = (S == S_int or S == frozenset(set(range(1, p)) - S_int))
                label = " ★PALEY" if is_paley else (" ★INT" if is_int else "")
                print(f"      H = {H}, S = {sorted(S)}{label}")
                count += 1
                if count >= 5:
                    break

        print(f"    Paley H = {paley_H}")
        print(f"    Interval H = {interval_H}")
        print(f"    Paley rank: {sorted(set(h for h, _ in all_H), reverse=True).index(paley_H) + 1}")
        print(f"    Interval rank: {sorted(set(h for h, _ in all_H), reverse=True).index(interval_H) + 1}")

    # Part 5: Analytical crossover formula
    print(f"\n{'='*75}")
    print("PART 5: CROSSOVER ANALYSIS — WHEN DOES INTERVAL BEAT PALEY?")
    print("-" * 75)

    # At p=7: Paley wins by 14 (189 vs 175)
    # At p=11: Paley wins by 2068 (95095 vs 93027)
    # At p=19: Interval wins by 11517077848

    # The gap as fraction of Paley H:
    data = [
        (7, 189, 175),
        (11, 95095, 93027),
        (19, 1172695746915, 1184212824763),
    ]

    print(f"\n  {'p':>4} {'H(Paley)':>20} {'H(Interval)':>20} {'Δ/H(P)':>12} {'Winner':>8}")
    for p, Hp, Hi in data:
        delta = (Hp - Hi) / Hp
        winner = "PALEY" if Hp > Hi else "INTERVAL"
        print(f"  {p:>4} {Hp:>20} {Hi:>20} {delta:>+12.6f} {winner:>8}")

    # The crossover: Paley wins by decreasing margin at p=7 (7.4%), p=11 (2.2%),
    # then loses at p=19 (-1.0%).
    # Extrapolation: crossover prime is between 11 and 19.

    print(f"\n  Paley advantage decreasing: 7.41% → 2.17% → -0.98%")
    print(f"  The crossover prime p₀ is between 11 and 19.")
    print(f"  The only p ≡ 3 mod 4 in (11, 19) is p = 15 (not prime).")
    print(f"  So p₀ = 19 is the FIRST prime where interval beats Paley.")

    # Part 6: Why the interval tournament creates flow
    print(f"\n{'='*75}")
    print("PART 6: THE FLOW STRUCTURE OF INTERVAL TOURNAMENTS")
    print("-" * 75)

    # The interval tournament C_p has the structure:
    # i → j iff 0 < (j-i) mod p ≤ m
    # This means: each vertex beats the next m vertices in cyclic order.
    # Hamiltonian paths can "ride the flow" — following the cyclic order.

    # Count how many HP follow mostly-forward steps
    for p in [7, 11]:
        m = (p - 1) // 2
        S_int = interval_set(p)
        A = circulant_adj(p, S_int)
        H = hamiltonian_paths_dp(A, p)

        # Count HPs by their "displacement pattern"
        # A step from v to w has displacement (w - v) mod p
        # Forward steps: displacement in {1,...,m}
        # All steps must be forward (by construction of interval tournament)
        # So the question is: how many permutations of Z_p have ALL
        # consecutive differences in {1,...,m}?

        # Actually in the interval tournament, i→j iff (j-i) mod p ∈ {1,...,m}
        # So EVERY edge is a "forward" step of size 1 to m.
        # The HP count is: number of permutations (v_0, v_1, ..., v_{p-1})
        # such that (v_{i+1} - v_i) mod p ∈ {1,...,m} for all i.

        print(f"\n  p = {p}: Interval tournament")
        print(f"    Every arc has displacement in {{1,...,{m}}}")
        print(f"    HP count = # permutations with all steps in {{1,...,{m}}} mod p")
        print(f"    H = {H}")

        # This is a known combinatorial object!
        # It's the number of "m-step" Hamiltonian paths on Z_p.
        # For m ≥ ⌈p/2⌉ (which it is, since m = (p-1)/2),
        # this counts permutations where no step exceeds m.

        # Compare with Paley:
        S_paley = paley_set(p)
        A_P = circulant_adj(p, S_paley)
        H_P = hamiltonian_paths_dp(A_P, p)
        print(f"    H(Paley) = {H_P}")

        # What's the step size distribution in Paley?
        print(f"    Paley connection set: {sorted(S_paley)}")
        print(f"    Step sizes: {sorted(S_paley)}")

        # The interval has consecutive steps {1,2,...,m}
        # This means paths can make small adjustments — step 1 is always available
        # Paley has GAPS in its step sizes
        paley_sorted = sorted(S_paley)
        gaps = [paley_sorted[i+1] - paley_sorted[i] for i in range(len(paley_sorted)-1)]
        print(f"    Paley gaps between consecutive steps: {gaps}")
        print(f"    Interval has NO gaps (steps 1,2,...,{m})")

    # Part 7: Connection to permanents and doubly stochastic matrices
    print(f"\n{'='*75}")
    print("PART 7: PERMANENT INTERPRETATION")
    print("-" * 75)

    # H(T) for a tournament T on Z_p can be related to a "permanent-like" quantity.
    # The adjacency matrix A of the tournament is {0,1}-valued with A + A^T = J - I.
    # H(T) = number of Hamiltonian paths = permanent of a certain matrix.
    #
    # Actually: H(T) = sum over all permutations σ of prod_{i=0}^{n-2} A[σ(i)][σ(i+1)]
    # This is NOT exactly a permanent, but a "path permanent" — the permanent
    # of a matrix where we sum over permutations that form a single path.

    print("  H(T) = Σ_σ Π_{i=0}^{n-2} A[σ(i), σ(i+1)]")
    print("  This is a 'path permanent' — sum over Hamiltonian paths, not matchings.")
    print()
    print("  Van der Waerden analogue:")
    print("  For doubly stochastic matrices: perm(J/n) ≤ perm(M) is FALSE")
    print("  (VdW says perm(J/n) ≤ perm(M) — wait, VdW says MINIMUM is at J/n)")
    print("  So the flat matrix MINIMIZES the permanent.")
    print()
    print("  For H on tournaments:")
    print("  At small p: flat spectrum (Paley) MAXIMIZES H")
    print("  At large p: concentrated spectrum (interval) MAXIMIZES H")
    print("  The crossover suggests H is NOT Schur-concave in general.")
    print("  This is consistent with THM-134 (Schur-concavity fails at p≡1(4)).")


if __name__ == '__main__':
    main()
    print("\nDONE.")
