"""
cross_field_connections.py — Creative connections to other mathematical fields

The tournament H-maximization problem has deep connections to:

1. STATISTICAL MECHANICS (Ising model)
   - The orientation cube {±1}^m IS the Ising configuration space
   - H(σ) is the partition function / energy
   - J[i,j] = degree-2 Walsh coefficients = pairwise spin couplings
   - Paley σ_P is the GROUND STATE of the 2-body Hamiltonian
   - The p≥19 crossover = phase transition where higher-body interactions dominate

2. ERROR-CORRECTING CODES (QR codes)
   - The Paley matrix is the core structure of Quadratic Residue codes
   - QR codes have optimal minimum distance (random-error correction)
   - Interval = burst-error codes (local structure)
   - H-maximization ↔ some coding-theoretic optimization?

3. ADDITIVE COMBINATORICS (sum-product phenomenon)
   - QR has LOW additive energy (Bourgain-Katz-Tao)
   - {1,...,m} has HIGH additive energy
   - The crossover: H REWARDS additive structure at large p

4. EXPANDER GRAPHS (Ramanujan property)
   - Paley = Ramanujan tournament (optimal expansion)
   - Interval = poor expansion but strong "flow"
   - H-maximization trades expansion for flow at large p

5. RANDOM MATRIX THEORY
   - Circulant = DFT-diagonalized
   - Paley spectrum ~ CUE (circular unitary ensemble)
   - Interval spectrum = rank-1 perturbation of flat

Author: opus-2026-03-12-S60b
"""
import sys
import time
import math
import numpy as np
from collections import defaultdict
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)


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


def additive_energy(S, p):
    """Additive energy E(S) = |{(a,b,c,d) ∈ S^4 : a+b = c+d mod p}|."""
    count = defaultdict(int)
    S = list(S)
    for a in S:
        for b in S:
            count[(a + b) % p] += 1
    return sum(v * v for v in count.values())


def main():
    print("CROSS-FIELD CONNECTIONS FROM TOURNAMENT H-MAXIMIZATION")
    print("=" * 75)

    # ================================================================
    # CONNECTION 1: ISING MODEL
    # ================================================================
    print(f"\n{'='*75}")
    print("1. ISING MODEL ON THE ORIENTATION CUBE")
    print("=" * 75)

    print("""
  The orientation cube {±1}^m parametrizes circulant tournaments on Z_p.

  The Walsh-Fourier expansion:
    H(σ) = H₀ + Σ_{i<j} J[i,j] σ_i σ_j + Σ_{i<j<k<l} K[i,j,k,l] σ_i σ_j σ_k σ_l + ...

  This IS the Ising model Hamiltonian on m spins with:
    - H₀ = external field (constant)
    - J[i,j] = pairwise spin coupling
    - K[i,j,k,l] = 4-body interaction

  Paley σ_P is the GROUND STATE of the 2-body term (maximizes σ^T J σ).
  At p≥19, higher-body interactions change the ground state to Interval.

  PHASE TRANSITION INTERPRETATION:
    Small p (= high temperature): 2-body dominance → Paley ground state
    Large p (= low temperature): many-body dominance → Interval ground state

  The "temperature" parameter is effectively log(p):
    - Number of spins m = (p-1)/2
    - Degree-2k Walsh term scales as ~p^{2k} (eigenvalue products)
    - Relative weight of degree-2k vs degree-2: ~p^{2(k-1)}
    - At large p, higher-degree terms inevitably dominate
    """)

    # Compute the "Ising parameters" at p=7
    p = 7
    m = (p - 1) // 2
    pairs = [(d, p - d) for d in range(1, m + 1)]

    # Compute H for all 2^m orientations
    H_values = {}
    for bits in range(1 << m):
        sigma = []
        S = set()
        for i in range(m):
            if bits & (1 << i):
                sigma.append(1)
                S.add(pairs[i][0])
            else:
                sigma.append(-1)
                S.add(pairs[i][1])
        A = circulant_adj(p, S)
        H = hamiltonian_paths_dp(A, p)
        H_values[tuple(sigma)] = H

    # Extract Walsh coefficients
    N = 1 << m
    print(f"  p={p}: Ising model on m={m} spins, {N} configurations")

    # Degree 0: mean
    H0 = sum(H_values.values()) / N
    print(f"  H₀ = {H0:.4f}")

    # Degree 2: J[i,j]
    J = np.zeros((m, m))
    for i in range(m):
        for j in range(i + 1, m):
            val = 0
            for sigma, H in H_values.items():
                val += sigma[i] * sigma[j] * H
            J[i][j] = J[j][i] = val / N
    print(f"  J matrix (pairwise couplings):")
    print(f"  {J}")

    # Eigenvalues of J
    eigenvalues = np.linalg.eigvalsh(J)
    print(f"  J eigenvalues: {eigenvalues}")

    # Paley sigma
    qr = paley_set(p)
    sigma_P = tuple(1 if pairs[i][0] in qr else -1 for i in range(m))
    print(f"  Paley σ = {sigma_P}")
    print(f"  σ^T J σ = {sum(J[i][j] * sigma_P[i] * sigma_P[j] for i in range(m) for j in range(m)):.4f}")

    # "Temperature" = ratio of degree-2 to higher-degree energy
    E2_paley = sum(J[i][j] * sigma_P[i] * sigma_P[j] for i in range(m) for j in range(m))
    E_total = H_values[sigma_P] - H0
    E_higher = E_total - E2_paley
    print(f"  E₂(Paley) = {E2_paley:.4f}")
    print(f"  E_total(Paley) = {E_total:.4f}")
    print(f"  E_higher(Paley) = {E_higher:.4f}")
    print(f"  Ratio E₂/E_total = {E2_paley/E_total:.4f}")

    # ================================================================
    # CONNECTION 2: ADDITIVE ENERGY
    # ================================================================
    print(f"\n{'='*75}")
    print("2. ADDITIVE COMBINATORICS: SUM-PRODUCT PHENOMENON")
    print("=" * 75)

    print("""
  The sum-product theorem (Bourgain-Katz-Tao, 2004):
    For S ⊂ Z_p with p^δ < |S| < p^{1-δ}:
    max(|S+S|, |S·S|) ≥ |S|^{1+ε}

  QR is "multiplicatively structured" (S·S = S) → |S+S| must be large
  → LOW additive energy

  Interval {1,...,m} is "additively structured" → HIGH additive energy

  H-maximization at large p REWARDS additive structure:
    - Hamiltonian paths follow "flows" — consecutive vertex visits
    - Additive structure = adjacent vertices have similar neighborhoods
    - This creates more path options (higher H)
    """)

    for p in [7, 11, 19, 23]:
        m = (p - 1) // 2
        S_paley = paley_set(p)
        S_interval = interval_set(p)

        E_paley = additive_energy(S_paley, p)
        E_interval = additive_energy(S_interval, p)

        # Trivial bound: E(S) ≥ |S|^3/p (random) and ≤ |S|^3 (arithmetic progression)
        trivial_lb = m**3 / p
        trivial_ub = m**3

        print(f"  p={p}: E(Paley) = {E_paley:>12}, E(Interval) = {E_interval:>12}, "
              f"ratio I/P = {E_interval/E_paley:.4f}")

    # ================================================================
    # CONNECTION 3: EXPANDER MIXING LEMMA
    # ================================================================
    print(f"\n{'='*75}")
    print("3. EXPANDER GRAPHS: SPECTRAL GAP vs HAMILTONIAN PATHS")
    print("=" * 75)

    print("""
  The Expander Mixing Lemma:
    For d-regular graph G with second eigenvalue λ₂:
    |e(S,T) - d|S||T|/n| ≤ λ₂ √(|S||T|)

  For tournaments (directed):
    - Paley: all |λ_k| = √((p+1)/4) → OPTIMAL expansion (Ramanujan)
    - Interval: dominant |λ_1| ≈ p/π → POOR expansion

  But H-maximization is NOT about expansion!
  It's about counting PATHS — a flow problem, not a mixing problem.

  The tradeoff:
    Expansion → edges go "everywhere" → many cycles but they overlap
    Flow → edges go "nearby" → fewer cycles but they're disjoint

  This is EXACTLY our OCF independence crossover:
    Expansion ↔ high α_1 (many cycles)
    Flow ↔ high α_2+ (many disjoint packings)
    """)

    for p in [7, 11, 19]:
        m = (p - 1) // 2
        # Spectral gap: λ_1 (excluding λ_0 = m)
        S_paley = paley_set(p)
        S_interval = interval_set(p)

        omega = np.exp(2j * np.pi / p)
        evals_P = [sum(omega ** (k * s) for s in S_paley) for k in range(1, p)]
        evals_I = [sum(omega ** (k * s) for s in S_interval) for k in range(1, p)]

        gap_P = max(abs(e) for e in evals_P)
        gap_I = max(abs(e) for e in evals_I)

        print(f"  p={p}: Paley λ_max = {gap_P:.4f} (≈√p/2), "
              f"Interval λ_max = {gap_I:.4f} (≈p/π)")

    # ================================================================
    # CONNECTION 4: CODING THEORY
    # ================================================================
    print(f"\n{'='*75}")
    print("4. CODING THEORY: QR CODES vs CYCLIC CODES")
    print("=" * 75)

    print("""
  Quadratic Residue (QR) codes:
    Generator polynomial g(x) = Π_{i ∈ QR} (x - α^i)
    where α is a primitive p-th root of unity in GF(2^k).

  The Paley tournament adjacency matrix IS the QR code check structure.
  QR codes achieve the highest known minimum distance for their length.

  Analogy:
    QR code: optimal error CORRECTION (max Hamming distance)
    Paley tournament: optimal path COUNT at small p

  Interval → cyclic code with generator g(x) = 1 + x + ... + x^{m-1}
    Cyclic code: optimal burst-error correction
    Interval tournament: optimal path COUNT at large p

  The crossover in H-maximization mirrors the crossover in coding theory:
    - Random errors → QR codes best (like small p = "random" structure)
    - Burst errors → cyclic codes best (like large p = "structured" flows)

  DEEP CONNECTION: Both use the SAME algebraic machinery:
    - Discrete Fourier Transform over Z_p
    - Eigenvalues = code distance spectrum
    - Weight enumerator polynomial ↔ independence polynomial I(Ω,x)
    """)

    # ================================================================
    # CONNECTION 5: TROPICAL GEOMETRY / ASSIGNMENT PROBLEM
    # ================================================================
    print(f"\n{'='*75}")
    print("5. TROPICAL GEOMETRY: PATH PERMANENT")
    print("=" * 75)

    print("""
  H(T) = Σ_{permutations π} Π_{i=0}^{n-2} A[π(i), π(i+1)]

  This is a "path permanent" — like the matrix permanent but summing
  over permutations viewed as paths, not matchings.

  The TROPICAL (min-plus) version:
    H_trop(T) = min_{π path} Σ_{i=0}^{n-2} cost(π(i), π(i+1))
    = the shortest Hamiltonian path (TSP)

  Connection:
    Ordinary permanent (Σ Π) → counts 1-factors → Brégman bound
    Path permanent (Σ Π over paths) → counts HP → Alon's bound
    Tropical permanent (min Σ) → shortest HP → TSP

  The Brégman bound that Alon uses connects permanent optimization
  to our H-maximization. Can we get a TIGHTER bound by using the
  path structure directly, without going through 1-factors?

  Conjecture: For regular tournaments on n vertices:
    H(T) ≤ (k!)^{(2n-1)/k} for k = (n-1)/2
    (tighter than Brégman for tournaments specifically)
    """)

    # ================================================================
    # CONNECTION 6: REPRESENTATION THEORY / HARMONIC ANALYSIS
    # ================================================================
    print(f"\n{'='*75}")
    print("6. REPRESENTATION THEORY: H AS A SPHERICAL FUNCTION")
    print("=" * 75)

    print("""
  H(T) for circulant T on Z_p is determined by the connection set S.
  S is a subset of (Z/pZ)* / {±1} (chord types).

  The space of such functions is the group algebra C[(Z/pZ)* / {±1}].

  For p ≡ 3 mod 4:
    (Z/pZ)* / {±1} ≅ QR_p (the quadratic residue group)
    This is a cyclic group of ODD order m = (p-1)/2.

  H is a function on the POWER SET of QR_p → integers.
  The Walsh expansion diagonalizes this as:
    H = Σ_S Ĥ(S) χ_S

  The Paley orientation σ_P = Legendre symbol = the UNIQUE
  nontrivial character of (Z/pZ)* of order 2.

  THM-137 says: σ_P is an eigenvector of J.
  This is SCHUR'S LEMMA applied to the QR action!

  DEEPER: The full Walsh expansion lives in the group algebra
  of the BOOLEAN group (Z/2Z)^m. The QR multiplication acts
  on this via permutations. The irreducible decomposition of
  H under this action determines which orientations maximize H.

  For the AFFINE GROUP Aff(Z_p) = {x → ax + b : a ∈ (Z/pZ)*, b ∈ Z_p}:
    - Paley tournament is the unique tournament invariant under Aff(Z_p)|_QR
    - Its automorphism group has order p·m = p(p-1)/2
    - This is the LARGEST possible for a circulant tournament
    """)

    # ================================================================
    # CONNECTION 7: QUANTUM INFORMATION
    # ================================================================
    print(f"\n{'='*75}")
    print("7. QUANTUM INFORMATION: MUTUALLY UNBIASED BASES")
    print("=" * 75)

    print("""
  Mutually Unbiased Bases (MUBs) in C^d:
    Two orthonormal bases B₁, B₂ are MU if |⟨b₁|b₂⟩|² = 1/d
    for all b₁ ∈ B₁, b₂ ∈ B₂.

  For d = p prime: the DFT matrix and its p-1 shifts give p+1 MUBs.
  The Paley tournament's FLAT spectrum means:
    |λ_k|² = (p+1)/4 for all k ≥ 1

  This is the tournament analogue of the MUB condition:
    All "measurements" (eigenvalue magnitudes) give equal information.

  The interval tournament CONCENTRATES information:
    |λ_1|² >> |λ_k|² for k > 1

  In quantum info terms:
    Paley = maximally uncertain (like MUB measurement)
    Interval = maximally informative (like computational basis)

  H-maximization at large p favoring CONCENTRATED information
  suggests: for counting Hamiltonian paths, it's better to have
  ONE strong "channel" than many weak ones.

  This connects to the "superadditivity" of quantum channels:
    Sometimes a single strong channel beats many independent weak ones.
    """)

    # ================================================================
    # CONNECTION 8: BOOLEAN FUNCTION COMPLEXITY
    # ================================================================
    print(f"\n{'='*75}")
    print("8. BOOLEAN FUNCTION COMPLEXITY: FOURIER DEGREE OF H")
    print("=" * 75)

    # Compute the Walsh/Fourier spectrum of H at p=7, 11
    for p in [7, 11]:
        m = (p - 1) // 2
        pairs = [(d, p - d) for d in range(1, m + 1)]

        # Compute H for all orientations
        all_H = {}
        for bits in range(1 << m):
            S = set()
            sigma = []
            for i in range(m):
                if bits & (1 << i):
                    S.add(pairs[i][0])
                    sigma.append(1)
                else:
                    S.add(pairs[i][1])
                    sigma.append(-1)
            A = circulant_adj(p, S)
            H = hamiltonian_paths_dp(A, p)
            all_H[tuple(sigma)] = H

        # Walsh transform
        N = 1 << m
        mean_H = sum(all_H.values()) / N
        total_energy = sum((H - mean_H)**2 for H in all_H.values()) / N

        # Energy by degree
        energy_by_deg = defaultdict(float)
        for size in range(m + 1):
            from itertools import combinations
            for subset in combinations(range(m), size):
                coeff = 0
                for sigma, H in all_H.items():
                    prod = 1
                    for idx in subset:
                        prod *= sigma[idx]
                    coeff += prod * H
                coeff /= N
                energy_by_deg[size] += coeff**2

        print(f"\n  p={p} (m={m}): Walsh spectrum of H")
        print(f"  Mean H = {mean_H:.2f}")
        for deg in sorted(energy_by_deg.keys()):
            if energy_by_deg[deg] > 0.001:
                frac = energy_by_deg[deg] / total_energy * 100 if total_energy > 0 else 0
                print(f"    Degree {deg}: energy = {energy_by_deg[deg]:.2f} ({frac:.1f}% of variance)")

    print("""
  KEY OBSERVATION:
  H(σ) on the Boolean cube is a LOW-DEGREE function!
  At p=7: ONLY degree 2 (quadratic)
  At p=11: degree 2 + degree 4 (quartic)
  No odd-degree terms (since H(σ) = H(-σ)).

  In Boolean function complexity:
  - Low degree → learnable by low-degree polynomial regression
  - Even-only → invariant under global flip (complement symmetry)

  The Fourier concentration at low degrees means:
  H is "smooth" on the Boolean cube — nearby orientations have similar H.
  This is JUNTA-like behavior: H depends on a few "important" chord types.

  Connection to computational learning theory:
  If H is a degree-2k polynomial on {±1}^m, it can be learned from
  O(m^{2k}) samples. For small k, this is efficient!
    """)


if __name__ == '__main__':
    main()
    print("\nDONE.")
