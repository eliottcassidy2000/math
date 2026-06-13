#!/usr/bin/env python3
"""
tasep_tournament_bridge.py — opus-2026-03-13-S67h

CREATIVE CONNECTION: Tournament Hamiltonian paths ↔ TASEP ↔ Last Passage Percolation

The Totally Asymmetric Simple Exclusion Process (TASEP) on a ring of p sites
has deep connections to:
  - KPZ universality class (H(m) scales with m^{4/3} exponent)
  - Last passage percolation (LPP) on directed graphs
  - Random matrix eigenvalue statistics (Tracy-Widom)
  - RSK correspondence / Schur functions

KEY INSIGHT: A tournament on Z_p with connection set S can be viewed as
defining a DIRECTED PERCOLATION model. A Hamiltonian path is a maximal
self-avoiding directed path — essentially a last-passage percolation path.

This script explores:
1. Encoding Hamiltonian paths as LPP paths on a torus
2. The connection to Robinson-Schensted-Knuth (RSK) correspondence
3. Tracy-Widom statistics for tournament path count fluctuations
4. The exact solvability via Bethe ansatz (circulant = translation-invariant)

Also explores:
5. Wiener-Hopf factorization of the transfer matrix
6. Szegő's theorem for Toeplitz determinants → log-gas interpretation
7. Fisher-Hartwig singularities at the phase transition
"""

import numpy as np
from math import pi, sin, cos, log, sqrt, comb, factorial
from itertools import combinations, permutations

phi = (1 + sqrt(5)) / 2

def is_prime(n):
    if n < 2: return False
    if n < 4: return n > 1
    if n % 2 == 0 or n % 3 == 0: return False
    d = 5
    while d*d <= n:
        if n % d == 0 or (n+2) % d == 0: return False
        d += 6
    return True

def count_H_fast(S, p):
    """Bitmask DP for Hamiltonian path count."""
    n = p
    adj = [[False]*n for _ in range(n)]
    for i in range(n):
        for s in S:
            j = (i + s) % n
            adj[i][j] = True
    dp = [{} for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    full = (1 << n) - 1
    count = 0
    for mask in range(1, 1 << n):
        for v in dp[mask]:
            if not dp[mask][v]: continue
            if mask == full:
                count += dp[mask][v]
                continue
            for u in range(n):
                if mask & (1 << u): continue
                if adj[v][u]:
                    if u not in dp[mask | (1 << u)]:
                        dp[mask | (1 << u)][u] = 0
                    dp[mask | (1 << u)][u] += dp[mask][v]
    return count

print("=" * 70)
print("CONNECTION: TOURNAMENT HAMILTONIAN PATHS AS LPP / RSK")
print("=" * 70)
print()

# 1. ENCODING AS LPP
print("1. ENCODING AS LAST PASSAGE PERCOLATION")
print("-" * 40)
print()
print("Given tournament T on Z_p (connection set S, |S|=m=(p-1)/2):")
print("  Define weight matrix W[i][j] = 1 if i→j, 0 otherwise")
print("  A Hamiltonian path maximizes Σ W[π(k)][π(k+1)] = p-1")
print("  (every step is a tournament edge)")
print()
print("This is a DIRECTED LPP on the complete graph K_p with")
print("weights drawn from the tournament. The 'passage time' is p-1")
print("for every Hamiltonian path (since tournaments are complete).")
print()
print("But the NUMBER of such paths (= H(T)) is the LPP MULTIPLICITY.")
print("In LPP language: all geodesics have length p-1,")
print("and H(T) counts the number of geodesics.")
print()

# 2. RSK CORRESPONDENCE
print("2. RSK CORRESPONDENCE")
print("-" * 40)
print()
print("RSK maps a permutation σ to a pair (P,Q) of standard Young tableaux.")
print("The shape λ of P,Q satisfies:")
print("  λ_1 = length of longest increasing subsequence of σ")
print("  (Vershik-Kerov / Baik-Deift-Johansson theorem: λ_1 ~ 2√n + n^{1/6}·χ_{TW})")
print()
print("For a tournament T on [p], each Hamiltonian path σ defines a permutation.")
print("The RSK shape λ(σ) captures the 'pattern structure' of the path.")
print()

# Let's actually compute RSK shapes for small tournaments
def rsk_insert(P, val):
    """RSK insertion of val into Young tableau P. Returns (P, bumped_row)."""
    P = [list(row) for row in P]
    row_idx = 0
    while True:
        if row_idx >= len(P):
            P.append([val])
            return P, row_idx
        row = P[row_idx]
        # Find first element > val
        pos = None
        for i, x in enumerate(row):
            if x > val:
                pos = i
                break
        if pos is None:
            row.append(val)
            return P, row_idx
        else:
            bumped = row[pos]
            row[pos] = val
            val = bumped
            row_idx += 1

def rsk_shape(perm):
    """Compute RSK shape (partition) for a permutation."""
    P = []
    for val in perm:
        P, _ = rsk_insert(P, val)
    return tuple(len(row) for row in P)

# For p=5 interval tournament, enumerate all Hamiltonian paths and compute RSK
for p in [5, 7]:
    m = (p-1)//2
    S = set(range(1, m+1))

    adj = [[False]*p for _ in range(p)]
    for i in range(p):
        for s in S:
            adj[i][(i+s)%p] = True

    # Enumerate Hamiltonian paths and compute RSK shapes
    shape_counts = {}
    total_paths = 0

    for start in range(p):
        # DFS
        stack = [(start, [start])]
        while stack:
            curr, path = stack.pop()
            if len(path) == p:
                total_paths += 1
                # Convert path to permutation (relative order)
                perm = path  # Already a permutation of {0,...,p-1}
                shape = rsk_shape(perm)
                shape_counts[shape] = shape_counts.get(shape, 0) + 1
                continue
            for j in range(p):
                if j not in path and adj[curr][j]:
                    stack.append((j, path + [j]))

    print(f"\np={p} Interval tournament: H={total_paths}, RSK shape distribution:")
    for shape, count in sorted(shape_counts.items(), key=lambda x: -x[1]):
        pct = count / total_paths * 100
        print(f"  {shape}: {count} paths ({pct:.1f}%)")

    # Compute average LIS (longest increasing subsequence)
    avg_lis = sum(shape[0] * count for shape, count in shape_counts.items()) / total_paths
    print(f"  Average LIS = {avg_lis:.2f} (random permutation: ~{2*sqrt(p):.2f})")

print()
print()
print("=" * 70)
print("3. TOEPLITZ DETERMINANT / SZEGŐ'S THEOREM CONNECTION")
print("=" * 70)
print()
print("The transfer matrix M(T) for a circulant tournament is a CIRCULANT")
print("matrix, hence diagonalized by the DFT. Its determinant is a product")
print("of eigenvalues — which are Fourier coefficients of the 'symbol function'.")
print()
print("Szegő's strong limit theorem: for Toeplitz matrices T_n with symbol f(θ),")
print("  det(T_n) ~ G(f)^n · E(f)")
print("where G(f) = exp((1/2π)∫ log f dθ) and E(f) depends on the Fourier")
print("coefficients of log f.")
print()
print("Our situation: the 'symbol' is f(θ) = 1 + Q(θ) where")
print("Q(θ) = Fejér kernel = sin²(mθ/2)/sin²(θ/2)")
print()
print("G(f) controls the exponential growth rate of det(T_n), analogous to")
print("how φ controls the growth of F_p.")
print()

# Compute G(f) = exp((1/2π)∫_0^{2π} log(1+F_m(θ)) dθ)
# We already computed (1/π)∫_0^π log(1+F_m) dθ ≈ 1.166 for large m.
# So G = exp(1.166/2) ≈ 1.79

for m in [10, 50, 100, 500, 1000, 5000]:
    theta = np.linspace(1e-10, pi, 100000)
    F = (np.sin(m*theta/2) / np.sin(theta/2))**2
    integral = np.trapezoid(np.log(1 + F), theta) / pi
    G = np.exp(integral / 2)
    print(f"  m={m:4d}: (1/π)∫ log(1+F_m) = {integral:.6f}, G = {G:.6f}")

print()
print("G converges to exp(1.1662/2) ≈ 1.7927 as m → ∞")
print(f"But we need G^m to match F_p ~ φ^{p}/√5")
print(f"G^m = 1.79^m while F_p = φ^(2m+1)/√5 = {phi:.4f}^(2m+1)/√5")
print()
print("The issue: Szegő gives det(T_n) ~ G^n, but our F_p = sum(det terms),")
print("not a single determinant. The permanent (sum over positive terms)")
print("grows faster than the determinant.")

print()
print()
print("=" * 70)
print("4. FISHER-HARTWIG SINGULARITIES")
print("=" * 70)
print()
print("The symbol f(θ) = 1 + sin²(mθ/2)/sin²(θ/2) has a SINGULARITY")
print("at θ = 0: f(0) = 1 + m². For Toeplitz asymptotics, this is a")
print("Fisher-Hartwig singularity of 'type 0' (jump discontinuity in log f).")
print()
print("The Fisher-Hartwig conjecture (proved by Deift-Its-Krasovsky):")
print("For symbols with r singularities at θ_j:")
print("  det(T_n) ~ C · n^{Σ(α_j²-β_j²)} · G^n · ∏ E_j")
print()
print("Our symbol: at θ=0, log f ≈ log(m²) + O(θ²/m²) → log f has a")
print("PEAK, not a jump. This gives α=0, β=0 (smooth case),")
print("consistent with standard Szegő behavior.")
print()
print("But the DISCRETE Riemann sum vs continuous integral gap is the key:")
print("the Fejér kernel peak captures O(1) sample points worth O(log m²)")
print("each, giving an EXTRA factor of m² in the discrete product vs")
print("the continuous determinant.")

print()
print()
print("=" * 70)
print("5. WIENER-HOPF FACTORIZATION OF THE TRANSFER MATRIX")
print("=" * 70)
print()
print("The Wiener-Hopf factorization: for a symbol f(θ) = f+(θ)·f-(θ)")
print("where f+ has Fourier coefficients supported on n≥0 and f- on n≤0.")
print()
print("For our symbol f(θ) = 1 + |D_m(θ)|²:")
print("  D_m(θ) = Σ_{k=1}^m e^{ikθ} = Dirichlet kernel")
print("  |D_m(θ)|² = D_m(θ)·D_m(-θ) = Σ_{j,k=1}^m e^{i(k-j)θ}")
print()
print("So f(θ) = 1 + Σ_{|n|≤m-1} (m-|n|) e^{inθ}")
print("        = 1 + m + 2 Σ_{n=1}^{m-1} (m-n) cos(nθ)")
print("        = m+1 + 2 Σ_{n=1}^{m-1} (m-n) cos(nθ)")
print()
print("This is a trigonometric polynomial of degree m-1.")
print("Its Wiener-Hopf factorization: f = g·g* where")
print("  g(θ) = Σ_{n=0}^{m-1} a_n e^{inθ}")
print("  g*(θ) = Σ_{n=0}^{m-1} ā_n e^{-inθ}")
print()
print("IMPORTANT: f(θ) = |1 + D_m(θ)|² = |Σ_{k=0}^m e^{ikθ}|²")
print("Wait, let me recheck: 1 + |D_m|² where D_m = Σ_{k=1}^m e^{ikθ}")
print("= 1 + |D_m|² ≠ |1+D_m|²  [unless D_m is real, which it's not]")
print()
print("But the F-product ∏(1+Q_k) = F_p IS computed from the eigenvalues")
print("Q_k = |D_m(2πk/p)|². The product of eigenvalues of the matrix")
print("I + DD* where D is the DFT of the indicator of S.")

# Verify this matrix interpretation
for p in [5, 7, 11]:
    m = (p-1)//2
    S = list(range(1, m+1))

    # DFT matrix approach
    omega = np.exp(2j*pi/p)
    D = np.zeros(p, dtype=complex)
    for s in S:
        D += omega**(np.arange(p)*s)

    # Eigenvalues of DD*
    Q_vals = np.abs(D)**2

    # Product 1+Q_k for k=1..p-1 (skip k=0)
    prod_val = np.prod(1 + Q_vals[1:])

    # Compare with Fibonacci
    fib_p = [0, 1]
    for _ in range(p-1):
        fib_p.append(fib_p[-1] + fib_p[-2])
    Fp = fib_p[p]

    # Also: Q_0 = m² (the DC component)
    print(f"  p={p}: ∏(1+Q_k, k=1..{p-1}) = {prod_val:.4f}, F_p = {Fp}")
    print(f"    Q_0 = {Q_vals[0]:.1f} = m² = {m**2}")
    print(f"    Product including k=0: {np.prod(1+Q_vals):.4f} = (1+m²)·F_p = {(1+m**2)*Fp}")

print()
print("Note: ∏_{k=0}^{p-1} (1+Q_k) = (1+m²)·F_p")
print("This is the determinant det(I + DD*) where D=diag(D_m(2πk/p)).")
print("In terms of the Fourier matrix: det(I + M) where M is the")
print("(p×p) circulant with first row [m, m-1, ..., 1, 0, ..., 0, 1, ..., m-1].")

print()
print()
print("=" * 70)
print("6. BETHE ANSATZ / COORDINATE BETHE ANSATZ")
print("=" * 70)
print()
print("The circulant structure of our tournament means the transfer matrix")
print("is TRANSLATION-INVARIANT on Z_p. This is the setting where the")
print("Bethe ansatz is exact!")
print()
print("In the XXZ spin chain / TASEP language:")
print("  - Sites: Z_p (the tournament vertices)")
print("  - Occupation: n_i = 1 if vertex i has been visited")
print("  - Hopping rate: c_{i→j} = 1 if (j-i) mod p ∈ S, 0 otherwise")
print("  - Boundary: periodic (ring)")
print()
print("The Bethe equations for p particles on p sites (fully packed):")
print("  e^{ik_j p} ∏_{l≠j} S(k_j, k_l) = 1  for j=1,...,p")
print()
print("where S(k,k') is the two-body S-matrix for particle scattering.")
print("For TASEP on Z_p: S(k,k') = (1-e^{ik'})/(1-e^{ik}).")
print()
print("The TOURNAMENT S-matrix is different: S encodes which direction")
print("particles can pass. The asymmetry is determined by S.")

# Let's compute something concrete: the transfer matrix eigenvalues
# for small p and compare with Bethe ansatz predictions
print()
print("Transfer matrix spectrum for Interval tournaments:")
for p in [5, 7]:
    m = (p-1)//2
    S_set = set(range(1, m+1))

    # Build transfer matrix: M[a][b] = #{Hamiltonian paths starting at a, ending at b}
    adj = [[False]*p for _ in range(p)]
    for i in range(p):
        for s in S_set:
            adj[i][(i+s)%p] = True

    # Count paths of length p-1 from a to b (Hamiltonian paths)
    M = np.zeros((p, p))
    for start in range(p):
        stack = [(start, frozenset([start]), start)]
        while stack:
            curr, visited, st = stack.pop()
            if len(visited) == p:
                M[st][curr] += 1
                continue
            for j in range(p):
                if j not in visited and adj[curr][j]:
                    stack.append((j, visited | {j}, st))

    evals = np.linalg.eigvals(M)
    evals_sorted = sorted(evals, key=lambda x: -abs(x))

    print(f"\n  p={p}: M eigenvalues (sorted by |λ|):")
    for i, ev in enumerate(evals_sorted):
        if abs(ev.imag) < 1e-8:
            print(f"    λ_{i} = {ev.real:.4f}")
        else:
            print(f"    λ_{i} = {ev.real:.4f} ± {ev.imag:.4f}i")

    print(f"    H(T) = {int(sum(M[0]))*p} (= p * row sum = {int(sum(M[0]))}*{p})")
    print(f"    M is {'circulant' if all(abs(M[0][(j-0)%p] - M[1][(j-1)%p]) < 1e-10 for j in range(p)) else 'NOT circulant'}")

print()
print()
print("=" * 70)
print("7. RANDOM MATRIX STATISTICS OF Q_k")
print("=" * 70)
print()
print("For large p, the bulk eigenvalues {Q_2,...,Q_m} should follow")
print("some universal distribution. Let's check if it's semicircular,")
print("Marchenko-Pastur, or something else.")
print()

for p in [997, 4999]:
    if not is_prime(p): continue
    m = (p-1)//2
    Q = [sin(m*pi*k/p)**2 / sin(pi*k/p)**2 for k in range(2, m+1)]

    # Normalize: Q_k for bulk modes, removing the giant Q_1
    Q_arr = np.array(Q)
    mean = np.mean(Q_arr)
    std = np.std(Q_arr)

    # Compute level spacings (for nearest-neighbor statistics)
    Q_sorted = np.sort(Q_arr)
    spacings = np.diff(Q_sorted)
    spacings_norm = spacings / np.mean(spacings)

    # Wigner surmise: P(s) = (π/2)s·exp(-πs²/4) for GUE
    # Poisson: P(s) = exp(-s)
    # Compute mean spacing ratio
    mean_s = np.mean(spacings_norm)
    var_s = np.var(spacings_norm)

    # r-ratio (adjacent gap ratio) used in RMT
    r_ratios = []
    for i in range(len(spacings)-1):
        r = min(spacings[i], spacings[i+1]) / max(spacings[i], spacings[i+1])
        r_ratios.append(r)
    mean_r = np.mean(r_ratios)

    print(f"  p={p:4d} (m={m}):")
    print(f"    Bulk stats: mean={mean:.4f}, std={std:.4f}, skew={float(np.mean(((Q_arr-mean)/std)**3)):.4f}")
    print(f"    Level spacing ratio <r> = {mean_r:.4f}")
    print(f"    (Poisson: <r>≈0.386, GOE: <r>≈0.536, GUE: <r>≈0.603)")
    print(f"    → Closest to: {'GUE' if abs(mean_r-0.603) < abs(mean_r-0.536) and abs(mean_r-0.603) < abs(mean_r-0.386) else 'GOE' if abs(mean_r-0.536) < abs(mean_r-0.386) else 'Poisson'}")

    # Tracy-Widom check: does the maximum bulk eigenvalue scale as m^{1/3}?
    Q_max_bulk = max(Q)
    # For Fejér kernel, Q_2 ≈ (p/(3π))² is O(m²/9) for the 2nd peak
    # Actually Q_2 ≈ F_m(2·2π/p) = sin²(2mπ/p)/sin²(2π/p) ≈ m²·sin²(π)/sin²(2π/p)
    # Wait, sin(2mπ/p) = sin(2·(p-1)/2·π/p) = sin(π(1-1/p)) = sin(π/p) → 0
    # So Q_2 → sin²(π/p)/sin²(2π/p) → 1/4 for large p!
    # That means Q_2 is O(1), not O(m²).
    print(f"    Q_2 (2nd mode) = {sin(m*pi*2/p)**2/sin(2*pi/p)**2:.4f} (O(1), not O(m²))")
    print(f"    Max bulk Q = {Q_max_bulk:.4f}")

print()
print()
print("=" * 70)
print("8. PARTITION FUNCTION AS FREDHOLM DETERMINANT")
print("=" * 70)
print()
print("In LPP and TASEP, the partition function can be written as a")
print("Fredholm determinant: Z = det(I - K) where K is an integral operator.")
print()
print("For our tournament: H(T) = p · ∏(1+Q_k) · A where A = amplification.")
print("The F-product ∏(1+Q_k) IS a finite Fredholm determinant:")
print("  F = det(I + DD*) / (1+m²)")
print("where D is the diagonal matrix of Dirichlet kernel values.")
print()
print("This connects to the Borodin-Okounkov identity:")
print("  det(I - K) = exp(Σ tr(K^n)/n) [when K is trace class]")
print()
print("In our case:")
print("  log F = Σ_{k=1}^{m} log(1+Q_k)")
print("  = Σ_{k=1}^m Σ_{n=1}^∞ (-1)^{n+1} Q_k^n / n")
print("  = Σ_{n=1}^∞ (-1)^{n+1} S_n / n  where S_n = Σ Q_k^n")
print()
print("This is the PLETHYSTIC LOGARITHM! And S_n are the power-sum")
print("symmetric functions of the eigenvalues.")
print()

# Compute this explicitly
for p in [7, 11, 13]:
    m = (p-1)//2
    Q = [sin(m*pi*k/p)**2 / sin(pi*k/p)**2 for k in range(1, m+1)]

    Fp_exact = 1
    for q in Q:
        Fp_exact *= (1 + q)

    # Plethystic expansion
    S = [sum(q**n for q in Q) for n in range(1, 10)]
    logF_series = sum((-1)**(n) * S[n-1] / n for n in range(1, 10))
    # Should equal log(Fp_exact) with alternating signs
    logF_exact = sum(log(1+q) for q in Q)

    print(f"\n  p={p}: F_p = {Fp_exact:.2f}")
    print(f"    log F (exact) = {logF_exact:.6f}")
    print(f"    Plethystic (9 terms) = {logF_series:.6f}")
    print(f"    Power sums S_n: {[f'{s:.2f}' for s in S[:5]]}")
    print(f"    S_1 = ΣQ = {S[0]:.4f} = m(m+1)/2 = {m*(m+1)/2}")

print()
print()
print("=" * 70)
print("SYNTHESIS: FIVE-WAY DICTIONARY")
print("=" * 70)
print()
print("Tournament    | TASEP          | LPP              | RMT            | Algebraic")
print("------------- | -------------- | ---------------- | -------------- | ---------")
print("Vertex i      | Site i         | Node (i,t)       | Eigenvalue i   | Root α_i")
print("Edge i→j      | Hop i→j        | Edge weight       | Level repul.   | Root diff")
print("Connection S  | Hop range      | Edge distribution | Potential      | Weyl group")
print("HP σ          | Trajectory     | Geodesic         | Eigenstate     | Weyl element")
print("H(T)          | Partition fn Z | # geodesics      | det(I+K)       | |W|")
print("Q_k           | Bethe root     | Shape param      | λ_k            | Character val")
print("F_p = ∏(1+Q)  | Free energy    | Weight gen fn    | Fredholm det   | Schur fn")
print("A (amplif.)   | Interaction    | Multiplicity     | Tracy-Widom    | Correction")
print("Interval S    | Nearest-neigh  | Directed lattice | GUE            | Longest element")
print("Paley S       | Long-range     | Random graph     | CUE            | Coxeter element")
print("m = (p-1)/2   | Density        | Dimension        | Matrix size    | Rank")
print("φ (growth)    | Current        | Shape limit      | Spectral edge  | Coxeter number")
print("8/π² (distrib)| Occupation     | Aspect ratio     | Bulk density   | Root density")
