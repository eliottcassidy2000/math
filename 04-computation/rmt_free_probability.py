#!/usr/bin/env python3
"""
Random matrix theory and free probability connections to H maximization.

THE SETUP (from circulant_walsh_analytic.py):
  For circulant tournament on Z_p with orientation σ ∈ {±1}^m:
    λ_t = a_t + i·b_t(σ)
  where:
    a_t = Σ_{j=1}^m cos(2πjt/p) = -1/2 for t ≥ 1  (CONSTANT!)
    b_t(σ) = Σ_{j=1}^m σ_j sin(2πjt/p)  (VARIABLE)

  The matrix B_{t,j} = sin(2πjt/p) has ALL SINGULAR VALUES EQUAL to √(p/4).
  So b = Bσ lives on a sphere of radius √(m·p/4) = m√(p/(4m)) = m/√2·√(p/(2m)).

  H is a function of |λ_t|² = 1/4 + b_t², which depends only on b_t².

KEY OBSERVATION FROM RMT:
  The CIRCULAR UNITARY ENSEMBLE (CUE) has eigenvalue distribution
  with |λ_t| = 1 (all on the unit circle). The Paley tournament's
  eigenvalues are on a circle of radius √((p+1)/4) → equidistribution
  on a circle, like CUE.

  The Interval tournament's eigenvalues are NOT equidistributed —
  they cluster near one point, like a DETERMINISTIC matrix.

  In free probability: the Paley tournament corresponds to a FREE
  random variable (semicircular distribution), while the Interval
  corresponds to a DETERMINISTIC one.

  FREE CONVEXITY: For free random variables, many functionals are
  CONVEX (entropy is concave, free energy is convex). This could
  give us the convexity we need!

THE PROOF IDEA:
  If H can be expressed as a functional of the eigenvalue distribution
  that is CONVEX under free convolution, then the deterministic
  (Interval) distribution maximizes it among all distributions with
  the same moments.

opus-2026-03-12-S66
"""

import numpy as np
from sympy.ntheory import legendre_symbol as legendre

def eigenvalues_circulant(sigma, p):
    m = (p - 1) // 2
    omega = np.exp(2j * np.pi / p)
    eigs = np.zeros(p, dtype=complex)
    for t in range(p):
        for j in range(1, m + 1):
            if sigma[j-1] == 1:
                eigs[t] += omega**(j * t)
            else:
                eigs[t] += omega**((p - j) * t)
    return eigs

# ========================================================================
print("=" * 72)
print("PART I: EIGENVALUE DISTRIBUTIONS AS PROBABILITY MEASURES")
print("=" * 72)

for p in [7, 13, 19, 29, 37]:
    m = (p - 1) // 2

    # Interval eigenvalues
    omega = np.exp(2j * np.pi / p)
    eigs_int = np.array([sum(omega**(j*t) for j in range(1, m+1)) for t in range(p)])
    lam_int = np.abs(eigs_int[1:])  # exclude t=0

    # Paley eigenvalues (if p ≡ 3 mod 4)
    qr_set = set()
    for a in range(1, p):
        qr_set.add((a*a) % p)
    S_pal = set(k for k in range(1, p) if k in qr_set)

    eigs_pal = np.array([sum(omega**(s*t) for s in S_pal) for t in range(p)])
    lam_pal = np.abs(eigs_pal[1:])

    print(f"\n  p={p}, m={m}:")

    # Moments of the eigenvalue distribution
    # μ_k = (1/(p-1)) Σ |λ_t|^k
    for k in [2, 4, 6, 8]:
        mu_int = np.mean(lam_int**k)
        mu_pal = np.mean(lam_pal**k)
        ratio = mu_int / mu_pal if mu_pal > 0 else float('inf')
        print(f"    μ_{k}: Int = {mu_int:.2f}, Pal = {mu_pal:.2f}, ratio = {ratio:.4f}")

    # Rényi entropy of |λ|² distribution
    p_int = lam_int**2 / np.sum(lam_int**2)
    p_pal = lam_pal**2 / np.sum(lam_pal**2)

    H2_int = -np.log2(np.sum(p_int**2))
    H2_pal = -np.log2(np.sum(p_pal**2))

    print(f"    H_2: Int = {H2_int:.4f} bits, Pal = {H2_pal:.4f} bits")
    print(f"    H_2 ratio: {H2_int/H2_pal:.4f}")

    # The KEY: the ratio of moments μ_k(Int)/μ_k(Pal) GROWS with k
    # This means the Interval concentrates weight on higher moments
    # → dominates in the "tail" → dominates in Z(Ω, 2) for large fugacity

# ========================================================================
print("\n" + "=" * 72)
print("PART II: MOMENT RATIO SCALING")
print("=" * 72)
print("""
For Interval: μ_k = (1/(p-1)) Σ |λ_t|^k ≈ |μ₁|^k / (p-1)
  (dominated by the top eigenvalue |μ₁| ≈ 2m/π)

For Paley: μ_k = (1/(p-1)) Σ |λ_t|^k = ((p+1)/4)^{k/2}
  (all eigenvalues have the same magnitude √((p+1)/4))

Ratio: μ_k(Int) / μ_k(Pal) ≈ (|μ₁|/√((p+1)/4))^k / (p-1)^{k/??}

Wait, more carefully:
  μ_k(Int) ≈ (2/p) |μ₁|^k  (two conjugate eigenvalues contribute)
  μ_k(Pal) = ((p+1)/4)^{k/2}

  Ratio ≈ (2/(p-1)) · (|μ₁|² / ((p+1)/4))^{k/2}

For large p: |μ₁| ≈ 2m/π ≈ p/π, so |μ₁|² ≈ p²/π²
  (p+1)/4 ≈ p/4

  Ratio ≈ (2/p) · (4p/π²)^{k/2} = (2/p) · (4/π²)^{k/2} · p^{k/2}

For k ≥ 3: this ratio → ∞ as p → ∞.
For k = 2: ratio → (2/p) · (4/π²) · p = 8/π² ≈ 0.81 (constant, < 1!)

So: μ_2(Int) < μ_2(Pal) for large p (Paley has more total "cycle energy")
But: μ_k(Int) > μ_k(Pal) for k ≥ 4 and p large enough.

This is the MOMENT CROSSOVER: higher moments favor Interval.
And H = Z(Ω, 2) depends on ALL moments (through the cycle counts).
""")

for p in [7, 13, 19, 29, 37, 53, 97]:
    m = (p - 1) // 2
    omega = np.exp(2j * np.pi / p)

    eigs_int = np.array([sum(omega**(j*t) for j in range(1, m+1)) for t in range(1, p)])
    eigs_pal_abs2 = np.full(p-1, (p+1)/4.0)  # exact for p ≡ 3 mod 4

    lam_int = np.abs(eigs_int)

    print(f"\n  p={p}:")
    for k in [2, 4, 6, 8, 10]:
        mu_int = np.mean(lam_int**k)
        mu_pal = ((p+1)/4.0)**(k/2)  # exact for flat spectrum
        ratio = mu_int / mu_pal
        print(f"    μ_{k:>2}: ratio(Int/Pal) = {ratio:.6f}")

# ========================================================================
print("\n" + "=" * 72)
print("PART III: THE CYCLE COUNT DICHOTOMY")
print("=" * 72)
print("""
From the trace formula: c_k = (1/p) Σ_t λ_t^k

For the Interval tournament (large p):
  c_k ≈ (2/p) |μ₁|^k cos(kθ₁) + (2/p) |μ₃|^k cos(kθ₃) + ...
  where θ_t is the argument of λ_t.

  Dominated by the TOP TWO eigenvalues (conjugate pair):
  c_k ≈ (2/p) |μ₁|^k cos(kθ₁)

  This OSCILLATES between positive and negative as k varies!
  For k where cos(kθ₁) > 0: Interval has MORE k-cycles
  For k where cos(kθ₁) < 0: Interval has FEWER k-cycles

For the Paley tournament (p ≡ 3 mod 4):
  c_k = (1/p) Σ_t (√((p+1)/4) · e^{iφ_t})^k
  where φ_t are the Gauss sum phases.

  All terms have the same magnitude → more "averaging" → c_k is smoother.

THE KEY: At the Hamiltonian length k=p:
  c_p(Int) > c_p(Pal) even though c_k(Int) < c_k(Pal) for most shorter k.
  (This is the "Hamiltonian cycle reversal" from synthesis_missing_links.out)

  Why? Because cos(p·θ₁) happens to be close to 1 for Interval.
  θ₁ ≈ π/(2m) → p·θ₁ ≈ p·π/(p-1) ≈ π → cos(pθ₁) ≈ -1 → c_p < 0??

  Wait, this doesn't match. Let me compute directly.
""")

for p in [7, 13, 19]:
    m = (p - 1) // 2
    omega = np.exp(2j * np.pi / p)

    eigs_int = np.array([sum(omega**(j*t) for j in range(1, m+1)) for t in range(p)])

    print(f"\n  p={p}: Top eigenvalue")
    lam1 = eigs_int[1]
    theta1 = np.angle(lam1)
    print(f"    λ₁ = {lam1:.4f}, |λ₁| = {abs(lam1):.4f}, arg = {theta1:.4f} rad = {np.degrees(theta1):.1f}°")

    # Cycle counts via trace formula
    for k in range(3, min(p+1, 20), 2):
        ck = np.real(sum(eigs_int**k)) / p
        ck_dom = 2 * abs(lam1)**k * np.cos(k * theta1) / p
        print(f"    c_{k:>2}: trace = {ck:.2f}, dominant = {ck_dom:.2f}, "
              f"cos({k}θ) = {np.cos(k*theta1):.4f}")

# ========================================================================
print("\n" + "=" * 72)
print("PART IV: THE LOG-PARTITION FUNCTION AND CONVEXITY")
print("=" * 72)
print("""
THE CENTRAL QUESTION: Why does max IPR → max H for large p?

APPROACH: Consider log H = log Z(Ω, 2) as a function of the eigenvalue
distribution. The partition function Z is a sum of exponentials (in the
independence polynomial expansion), so log Z is a "logsumexp" function.

CLAIM: log Z(Ω, λ) is CONVEX in λ for large enough λ.
  (This is because the dominant term 2^k α_k for the largest k
   grows exponentially with λ.)

But we need convexity in the eigenvalue DISTRIBUTION, not in λ.

ALTERNATIVE: Think of the eigenvalue distribution as a PROBABILITY measure
μ on the complex plane. For each μ:
  1. Compute the cycle counts c_k(μ) via the trace formula
  2. Build the overlap graph Ω(μ)
  3. Compute H(μ) = Z(Ω(μ), 2)

If H is convex in μ (in the sense of optimal transport), then
the most concentrated μ (Dirac delta → Interval) maximizes H.

This connects to the MONGE-KANTOROVICH theory:
  The optimal transport cost from Interval to Paley is related to
  the Wasserstein distance between their eigenvalue distributions.

For large p: the Wasserstein distance W₂(μ_Int, μ_Pal) → ∞
  (Interval concentrates while Paley spreads over a circle)
  So if H has any convexity, the advantage of Interval grows with p.
""")

# Compute Wasserstein-like distances
for p in [7, 13, 19, 29, 53]:
    m = (p - 1) // 2
    omega = np.exp(2j * np.pi / p)

    eigs_int = np.array([sum(omega**(j*t) for j in range(1, m+1)) for t in range(1, p)])
    lam_int = np.abs(eigs_int)

    # Sorted eigenvalue magnitudes
    lam_int_sorted = np.sort(lam_int)[::-1]
    lam_pal_sorted = np.full(p-1, np.sqrt((p+1)/4.0))

    # L2 distance between sorted eigenvalue distributions
    W2 = np.sqrt(np.mean((lam_int_sorted - lam_pal_sorted)**2))

    # L∞ distance
    Linf = np.max(np.abs(lam_int_sorted - lam_pal_sorted))

    # Concentration ratio: max|λ| / mean|λ|
    conc_int = np.max(lam_int) / np.mean(lam_int)
    conc_pal = 1.0  # all equal

    print(f"  p={p}: W₂ = {W2:.4f}, L∞ = {Linf:.4f}, "
          f"conc(Int) = {conc_int:.4f}, conc(Pal) = 1.0000")

# ========================================================================
print("\n" + "=" * 72)
print("PART V: THE b-VECTOR OPTIMIZATION LANDSCAPE")
print("=" * 72)
print("""
Since B is orthogonal (scaled), the map σ → b = Bσ is a rotation.
The hypercube {±1}^m maps to a DISCRETE SET on a sphere.

All b-vectors have the same norm (||b||² = m · p/4).
H depends on b through |λ_t|² = 1/4 + b_t².

Define: F(b) = H(eigenvalues with |λ_t|² = 1/4 + b_t²)

The optimization is: max_{b ∈ B({±1}^m)} F(b)

Since B is orthogonal, the vertices of the hypercube map to the
vertices of a ROTATED hypercube on the sphere. The question is:
which vertex maximizes F?

THE INTERVAL VERTEX: b_Int = B · (1,...,1)
  This has b_1 = Σ sin(2πj/p) ≈ cot(π/(2m+2)) (very large)
  and b_t ≈ 0 for t ≥ 2.

THE PALEY VERTEX: b_Pal = B · σ_P
  This has all |b_t| ≈ √(p/4) (equidistributed).

F is maximized at the Interval vertex for p ≥ 13.
This is equivalent to saying: the function F prefers "spike" distributions
(one large b_t, rest small) over "flat" distributions.

Is F convex on the sphere? If so, the vertices of the hypercube
(which are the "most extreme" points) would all achieve larger F
than interior points. And among the vertices, the one with the
most extreme b-distribution would win.

Actually: F is NOT convex on the sphere (we showed violations at p=13).
But it IS maximized at the most extreme vertex.
""")

# Compute F at various intermediate points (interpolations)
for p in [13]:
    m = (p - 1) // 2

    # B matrix
    B = np.zeros((m, m))
    for t in range(m):
        for j in range(m):
            B[t][j] = np.sin(2 * np.pi * (j+1) * (t+1) / p)

    sigma_int = np.ones(m)
    b_int = B @ sigma_int

    # Paley sigma
    sigma_pal = np.array([float(legendre(k, p)) for k in range(1, m+1)])
    b_pal = B @ sigma_pal

    print(f"\n  p={p}: Interpolation between Interval and Paley b-vectors")
    print(f"  b_Int: {[f'{x:.3f}' for x in b_int]}")
    print(f"  b_Pal: {[f'{x:.3f}' for x in b_pal]}")

    # Note: interpolation leaves the discrete hypercube, so we can't compute H directly.
    # But we CAN compute H at nearby hypercube vertices.
    # Let's measure the Hamming distance from Interval and see how H varies.

    print(f"\n  H vs Hamming distance from Interval:")
    from collections import defaultdict
    import time

    def orientation_to_S_local(sigma, p):
        m = (p - 1) // 2
        S = set()
        for k in range(1, m + 1):
            if sigma[k-1] == 1:
                S.add(k)
            else:
                S.add(p - k)
        return S

    def make_tournament_local(p, S):
        A = [[0]*p for _ in range(p)]
        for i in range(p):
            for d in S:
                A[i][(i+d)%p] = 1
        return A

    def count_H_local(A):
        n = len(A)
        dp = [[0]*n for _ in range(1 << n)]
        for i in range(n):
            dp[1 << i][i] = 1
        for mask in range(1, 1 << n):
            for v in range(n):
                if not (mask & (1 << v)):
                    continue
                if dp[mask][v] == 0:
                    continue
                for u in range(n):
                    if mask & (1 << u):
                        continue
                    if A[v][u]:
                        dp[mask | (1 << u)][u] += dp[mask][v]
        full = (1 << n) - 1
        return sum(dp[full][v] for v in range(n))

    hamming_H = defaultdict(list)
    N = 2**m

    for bits in range(N):
        sigma = tuple(1 if (bits >> k) & 1 else -1 for k in range(m))
        # Hamming distance from all-1
        d = sum(1 for k in range(m) if sigma[k] == -1)
        S = orientation_to_S_local(sigma, p)
        A = make_tournament_local(p, S)
        H = count_H_local(A)
        hamming_H[d].append(H)

    print(f"  {'Hamming d':>10} {'#orient':>8} {'min H':>12} {'max H':>12} {'mean H':>12}")
    for d in sorted(hamming_H.keys()):
        vals = hamming_H[d]
        print(f"  {d:>10} {len(vals):>8} {min(vals):>12} {max(vals):>12} {np.mean(vals):>12.1f}")

# ========================================================================
print("\n" + "=" * 72)
print("PART VI: THE PROOF PATH — PUTTING IT ALL TOGETHER")
print("=" * 72)
print("""
FROM ALL THE EVIDENCE, the proof structure is:

THEOREM (Interval H-Maximization, p ≥ 13):
  Among all circulant tournaments on Z_p (p ≥ 13 prime), the Interval
  tournament S = {1, 2, ..., (p-1)/2} maximizes the Hamiltonian path count H.

PROOF ARCHITECTURE:

  1. ALGEBRAIC REFORMULATION:
     H(σ) = Σ_{S⊂[m]} ĥ[S] · Π_{k∈S} σ_k    (Walsh expansion)
     H(σ) depends on σ only through b = Bσ where B is the DFT sine matrix
     Since B has constant singular values, H is a function on a sphere

  2. THE WALSH DIFFERENCE FORMULA (proved for p ≡ 3 mod 4):
     H(1,...,1) - H(σ_P) = 2 · Σ_{NQR products} ĥ[S]

  3. EXHAUSTIVE VERIFICATION (computed):
     p = 13: H(Interval) = 3,711,175 = max over all 64 orientations
     p = 17: H(Interval) = 13,689,269,499 = max over all 256 orientations

  4. ASYMPTOTIC SCALING (from THM-142 + entropy analysis):
     The disjointness excess Δα₂ grows as O(p⁴) [THM-142]
     The cycle count advantage Δα₁ grows slower (relative to α₁)
     The H correlation with IPR INCREASES with p (0.51 → 0.71 → ...)
     The Rényi entropy ratio decreases: H₂(Int)/H₂(Pal) → 0

  5. FOR A COMPLETE PROOF: need one of:
     (a) Exhaustive verification at p = 19, 23 (computationally hard)
     (b) Analytic bound: H(Int) - H(σ) ≥ 0 for all σ using Walsh structure
     (c) Induction: if true at p, then true at p' > p
         (using the increasing correlation as evidence)

  The STRONGEST current approach is (b): the Walsh expansion gives
  H(Int) - H(σ) = 2·Σ_{ψ(S,σ)=-1} ĥ[S]
  where ψ(S,σ) = Π_{k∈S} σ_k.

  For ANY σ, the set {S : ψ(S,σ) = -1} contains exactly 2^{m-1} subsets.
  If Σ_{NQR(σ)} ĥ[S] ≥ 0 for ALL σ (not just Paley), then Interval wins.

  This is a HYPERPLANE INEQUALITY on the Walsh spectrum:
  for each σ, the half-space defined by σ must have non-negative
  projection of ĥ.

  At p=13: this holds (verified by exhaustive search).
  At p=17: this holds (verified by exhaustive search).

  The question is: does it hold for all p ≥ 13?
""")

print("\nDONE.")
