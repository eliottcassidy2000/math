"""
dihedral_paley.py — Dihedral group structure of tournaments and H-maximization

Key insight: vertices on unit circle → D_{2n} symmetry group.
- Every circulant tournament: rotations (Z_n) are automorphisms,
  reflections (k↦-k) are anti-automorphisms
- D_{2n} = Z_n ⋊ Z_2 where Z_2 = complement duality
- Paley tournament: ADDITIONALLY has Aut ⊃ Z_p ⋊ QR_p (affine group)

This script:
1. Enumerate all circulant tournaments on Z_n for odd n
2. Compute H(T) for each via DP
3. Compute eigenvalues and spectral flatness
4. Identify the maximizer and its connection to D_{2n} structure
5. Decompose H(T) using D_{2n} representation theory

Author: opus-2026-03-12-S60
"""
import sys
import time
import numpy as np
from itertools import combinations
from collections import defaultdict
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)


def circulant_tournaments(n):
    """Enumerate all circulant tournaments on Z_n (n odd).

    A circulant tournament has connection set S ⊂ Z_n \ {0} with |S| = (n-1)/2,
    such that for each d ≠ 0, exactly one of {d, n-d} is in S.

    The constraint S ∩ (-S) = ∅ means we choose one from each pair {d, n-d}.
    For odd n, the pairs are {1, n-1}, {2, n-2}, ..., {(n-1)/2, (n+1)/2}.
    """
    k = (n - 1) // 2
    pairs = [(d, n - d) for d in range(1, k + 1)]
    # Each pair contributes one element to S: either d or n-d
    for bits in range(1 << k):
        S = set()
        for i in range(k):
            if bits & (1 << i):
                S.add(pairs[i][0])
            else:
                S.add(pairs[i][1])
        yield frozenset(S)


def adjacency_matrix(n, S):
    """Build adjacency matrix for circulant tournament on Z_n with connection set S."""
    A = np.zeros((n, n), dtype=int)
    for i in range(n):
        for s in S:
            j = (i + s) % n
            A[i, j] = 1
    return A


def hamiltonian_paths_dp(A, n):
    """Count Hamiltonian paths by bitmask DP."""
    # dp[mask][v] = number of Hamiltonian paths ending at v with visited set = mask
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
                if A[v, w]:
                    dp[mask | (1 << w)][w] += dp[mask][v]

    return sum(dp[full][v] for v in range(n))


def eigenvalues_circulant(n, S):
    """Compute eigenvalues of circulant adjacency matrix.

    For circulant with connection set S:
      λ_k = Σ_{s ∈ S} ω^{ks}, k = 0, ..., n-1
    where ω = e^{2πi/n}.
    """
    omega = np.exp(2j * np.pi / n)
    evals = np.zeros(n, dtype=complex)
    for k in range(n):
        evals[k] = sum(omega ** (k * s) for s in S)
    return evals


def spectral_flatness(evals):
    """Measure how "flat" the eigenvalue magnitudes are (excluding k=0).
    Flatness = 0 means all |λ_k| equal for k ≥ 1.
    """
    mags = np.abs(evals[1:])
    if len(mags) == 0:
        return 0.0
    return np.std(mags) / (np.mean(mags) + 1e-30)


def is_paley_set(n, S):
    """Check if S equals the set of quadratic residues mod n (n prime)."""
    if not is_prime(n):
        return False
    qr = set()
    for x in range(1, n):
        qr.add((x * x) % n)
    return S == frozenset(qr) or S == frozenset(set(range(1, n)) - qr)


def is_prime(n):
    if n < 2:
        return False
    for p in [2, 3, 5, 7, 11, 13]:
        if n == p:
            return True
        if n % p == 0:
            return False
    d = 17
    while d * d <= n:
        if n % d == 0:
            return False
        d += 2
    return True


def gauss_sum(p):
    """Compute the quadratic Gauss sum g = Σ_{a=0}^{p-1} (a/p) ω^a
    where (a/p) is the Legendre symbol and ω = e^{2πi/p}.

    For p ≡ 3 mod 4: g = i√p
    For p ≡ 1 mod 4: g = √p
    """
    omega = np.exp(2j * np.pi / p)
    g = 0
    for a in range(p):
        leg = legendre(a, p)
        g += leg * omega ** a
    return g


def legendre(a, p):
    """Legendre symbol (a/p)."""
    a = a % p
    if a == 0:
        return 0
    if pow(a, (p - 1) // 2, p) == 1:
        return 1
    return -1


def connection_set_type(n, S):
    """Classify the connection set S."""
    S_set = set(S)
    # Check if cyclic interval: S = {a, a+1, ..., a+(n-1)/2-1} mod n
    k = (n - 1) // 2
    for start in range(1, n):
        interval = frozenset((start + i) % n for i in range(k))
        # Exclude 0 from interval
        if 0 in interval:
            continue
        if frozenset(S) == interval:
            return f"interval[{start}..{(start+k-1)%n}]"

    # Check if QR
    if is_prime(n):
        qr = set()
        for x in range(1, n):
            qr.add((x * x) % n)
        if S_set == qr:
            return "QR (Paley)"
        if S_set == set(range(1, n)) - qr:
            return "NQR (Paley complement)"

    return "other"


def aut_group_order(n, S):
    """Compute |Aut(T)| for circulant tournament with connection set S on Z_n.
    Automorphisms must be permutations π with i→j iff π(i)→π(j).
    For circulant: translations x↦x+a are always automorphisms.
    Additional: x↦bx for b such that bS = S.
    """
    S_set = set(S)
    # Find multipliers b such that {b*s mod n for s in S} = S
    multipliers = 0
    for b in range(1, n):
        if all((b * s) % n in S_set for s in S_set):
            multipliers += 1
    return n * multipliers


def dihedral_character_decomposition(n, S):
    """Decompose the H(T) computation using D_{2n} representations.

    D_{2n} for odd n has irreducible representations:
    - ρ_0: trivial (dim 1)
    - ρ_1: sign representation (rotation → 1, reflection → -1) (dim 1)
    - ρ_k: 2-dim rep indexed by k=1,...,(n-1)/2

    The eigenvalues λ_k of the circulant can be organized by D_{2n} irreps:
    - ρ_0 eigenvalue: λ_0 = |S| = (n-1)/2
    - ρ_k eigenvalue pair: (λ_k, λ_{n-k}) = (λ_k, λ_k*) since λ_{n-k} = conj(λ_k)

    Key: λ_k and λ_{n-k} are paired by the reflection symmetry!
    """
    evals = eigenvalues_circulant(n, S)
    m = (n - 1) // 2

    decomp = {}
    decomp['trivial'] = {'eval': evals[0].real, 'dim': 1}

    for k in range(1, m + 1):
        # D_{2n} pairs eigenspace k with eigenspace n-k
        # The 2-dim irrep has trace = 2*Re(λ_k) for rotations
        # and trace = 0 for reflections (since λ_{n-k} = conj(λ_k))
        pair_mag = abs(evals[k])
        pair_phase = np.angle(evals[k])
        decomp[f'ρ_{k}'] = {
            'eval_k': evals[k],
            'eval_nk': evals[n - k],
            'magnitude': pair_mag,
            'phase': pair_phase,
            'real_part': evals[k].real,
        }

    return decomp


def main():
    print("DIHEDRAL GROUP STRUCTURE OF TOURNAMENTS AND H-MAXIMIZATION")
    print("=" * 75)

    for n in [3, 5, 7, 9, 11, 13]:
        print(f"\n{'='*75}")
        print(f"n = {n}  |  D_{{{2*n}}}  |  {2**((n-1)//2)} circulant tournaments")
        print(f"{'='*75}")

        is_paley_prime = is_prime(n) and n % 4 == 3
        print(f"  Paley prime (p ≡ 3 mod 4): {'YES' if is_paley_prime else 'NO'}")
        if is_prime(n):
            print(f"  n mod 4 = {n % 4}, n mod 8 = {n % 8}")

        tournaments = list(circulant_tournaments(n))
        print(f"  Number of circulants: {len(tournaments)}")

        if n > 13:
            print(f"  (Skipping DP — too large)")
            continue

        # Compute H for each
        results = []
        t0 = time.time()
        for S in tournaments:
            A = adjacency_matrix(n, S)
            H = hamiltonian_paths_dp(A, n)
            evals = eigenvalues_circulant(n, S)
            sf = spectral_flatness(evals)
            ct = connection_set_type(n, S)
            ag = aut_group_order(n, S)
            results.append({
                'S': S,
                'H': H,
                'evals': evals,
                'spectral_flatness': sf,
                'type': ct,
                'aut_order': ag,
            })
        elapsed = time.time() - t0
        print(f"  Computed in {elapsed:.2f}s")

        # Sort by H descending
        results.sort(key=lambda r: -r['H'])

        # Print ranking
        print(f"\n  --- H RANKING ---")
        max_H = results[0]['H']
        for rank, r in enumerate(results):
            marker = " ★" if r['H'] == max_H else ""
            evals_str = ", ".join(f"{abs(r['evals'][k]):.3f}" for k in range(1, min(4, n)))
            print(f"  #{rank+1}: H={r['H']:>12}  |S|={len(r['S'])}  "
                  f"|Aut|={r['aut_order']:>4}  "
                  f"SF={r['spectral_flatness']:.4f}  "
                  f"type={r['type']:<25}  "
                  f"|λ|=[{evals_str}]"
                  f"{marker}")
            if rank >= 9 and r['H'] < max_H:
                remaining = len(results) - rank - 1
                if remaining > 0:
                    print(f"  ... ({remaining} more with H < {r['H']})")
                break

        # Analyze the maximizer
        best = results[0]
        print(f"\n  --- MAXIMIZER ANALYSIS ---")
        print(f"  Connection set S = {sorted(best['S'])}")
        print(f"  H = {best['H']}")
        print(f"  |Aut(T)| = {best['aut_order']}")
        print(f"  Spectral flatness = {best['spectral_flatness']:.6f}")
        print(f"  Type: {best['type']}")

        # Eigenvalue decomposition via D_{2n}
        decomp = dihedral_character_decomposition(n, best['S'])
        print(f"\n  --- D_{{{2*n}}} EIGENVALUE DECOMPOSITION ---")
        print(f"  Trivial: λ_0 = {decomp['trivial']['eval']:.3f} (= (n-1)/2 = {(n-1)/2})")
        m = (n - 1) // 2
        for k in range(1, m + 1):
            d = decomp[f'ρ_{k}']
            print(f"  ρ_{k}: λ_{k} = {d['eval_k']:.3f}, |λ_{k}| = {d['magnitude']:.4f}, "
                  f"phase = {d['phase']*180/np.pi:.1f}°")

        if is_prime(n):
            # Gauss sum
            g = gauss_sum(n)
            print(f"\n  --- GAUSS SUM ---")
            print(f"  g(η,ψ) = {g:.4f}")
            print(f"  |g| = {abs(g):.4f} (expected √{n} = {np.sqrt(n):.4f})")
            print(f"  g²/n = {(g**2/n):.4f} (expected {'i' if n%4==3 else '1'} × sign)")

        # Compare: spectral flatness vs H correlation
        Hs = [r['H'] for r in results]
        SFs = [r['spectral_flatness'] for r in results]
        if len(set(Hs)) > 1:
            corr = np.corrcoef(Hs, SFs)[0, 1]
            print(f"\n  Correlation(H, spectral_flatness) = {corr:.4f}")
            print(f"  {'NEGATIVE = flat spectrum → more H' if corr < 0 else 'POSITIVE = peaked spectrum → more H'}")

        # Aut group order vs H
        Auts = [r['aut_order'] for r in results]
        if len(set(Auts)) > 1 and len(set(Hs)) > 1:
            corr_aut = np.corrcoef(Hs, Auts)[0, 1]
            print(f"  Correlation(H, |Aut|) = {corr_aut:.4f}")
            print(f"  {'POSITIVE = more symmetry → more H' if corr_aut > 0 else 'NEGATIVE'}")

        # The dihedral "pairing" structure
        print(f"\n  --- DIHEDRAL PAIRING: λ_k and λ_{{n-k}} ---")
        for k in range(1, m + 1):
            lk = best['evals'][k]
            lnk = best['evals'][n - k]
            print(f"  (λ_{k}, λ_{n-k}) = ({lk:.3f}, {lnk:.3f}), "
                  f"sum = {(lk+lnk):.4f} (= 2·Re(λ_{k}))")


if __name__ == '__main__':
    main()
    print("\nDONE.")
