#!/usr/bin/env python3
"""
kam_stability_connection.py — opus-2026-03-12-S67d

CREATIVE CROSS-FIELD CONNECTION: KAM Theory & Golden Ratio Stability

The golden ratio φ is the "most irrational" number — its continued fraction
[1;1,1,1,...] has the slowest convergence of any CF. In KAM theory,
dynamical systems with frequency ratio φ are the MOST STABLE against
perturbation.

KEY OBSERVATIONS:
1. Interval tournament eigenvalues: Q_k^{1/m} → φ² as m → ∞
2. The ratio B_m/B_{m-1} → φ² = [3;3,3,...] (periodic CF)
3. Interval has the MOST CONCENTRATED Q-spectrum (max Q_1, min others)
4. This concentration is analogous to a system "locked" to one frequency

ANALOGY TABLE:
  KAM theory              | Tournament theory
  ========================|========================
  Frequency ratio ω       | Q_1/sum(Q) (spectral concentration)
  Golden ratio φ           | φ² = B_m/B_{m-1} limit
  Perturbation ε          | Flipping σ_i → -σ_i
  Invariant torus          | Maximum-H configuration
  Resonance destruction    | H drops when moving away from Interval
  Small divisors 1/(qω-p) | 1/(Q_k - 1) corrections in the OCF

ADDITIONAL CONNECTIONS:
A. Penrose tilings & quasicrystals
B. Zeckendorf representation (every integer as sum of non-consecutive Fibonacci)
C. Beatty sequences and the Wythoff game
D. The Fibonacci sequence in biology (phyllotaxis)
"""

import numpy as np
from math import comb, gcd
from fractions import Fraction

def morgan_voyce_B(m, x):
    return sum(comb(m+j, 2*j) * x**j for j in range(m+1))

def fib(n):
    a, b = 0, 1
    for _ in range(n):
        a, b = b, a+b
    return a

def interval_Q(p):
    m = (p-1)//2
    return np.array([
        (np.sin(m*np.pi*k/p)/np.sin(np.pi*k/p))**2
        for k in range(1, m+1)
    ])

phi = (1 + np.sqrt(5)) / 2
psi = (1 - np.sqrt(5)) / 2

print("=" * 70)
print("KAM THEORY & GOLDEN RATIO STABILITY — opus-2026-03-12-S67d")
print("=" * 70)
print()

# ============================================================
# PART 1: SPECTRAL CONCENTRATION AND THE GOLDEN RATIO
# ============================================================

print("PART 1: SPECTRAL CONCENTRATION → φ²")
print("=" * 70)
print()

print("The Interval tournament Q-spectrum concentrates:")
print("  Q_1 = m²(1 - O(1/p²)), Q_k << Q_1 for k ≥ 2")
print()
print("Spectral concentration ratio R = Q_1 / sum(Q_k):")
print("  sum(Q_k) = m(m+1)/2 (constant for all circulant tournaments)")
print()

for p in [5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47]:
    m = (p-1)//2
    Q = interval_Q(p)
    Q1 = Q[0]
    sum_Q = np.sum(Q)
    R = Q1 / sum_Q

    # Also compute the Herfindahl index (IPR)
    IPR = np.sum(Q**2) / sum_Q**2

    # Effective number of frequencies
    N_eff = 1 / IPR

    # Entropy of Q/sum(Q) distribution
    probs = Q / sum_Q
    entropy = -np.sum(probs * np.log(probs))

    print(f"  p={p:2d}, m={m:2d}: Q_1/Σ = {R:.6f}, IPR = {IPR:.6f}, "
          f"N_eff = {N_eff:.3f}, S = {entropy:.4f}")

print()
print(f"As m→∞: Q_1/Σ → 2/m · (1/sin²(π/p))/(m+1) → ...")
print("The concentration R → 2/3 for large m (since Q_1 ≈ m² and Σ ≈ m²/2 + m/3)")
print()

# ============================================================
# PART 2: CONTINUED FRACTION STABILITY
# ============================================================

print("=" * 70)
print("PART 2: CONTINUED FRACTION AND STABILITY")
print("=" * 70)
print()

# In KAM theory, the stability of an orbit with frequency ω depends on
# how well ω resists rational approximation. The measure is:
# |ω - p/q| ≥ C/q^τ for all p/q
# φ has τ = 2 (best possible) and C = 1/√5 (also optimal among τ=2)

# For our tournament, the analogous quantity is:
# How well does the Q-spectrum resist being "rationalized" (made uniform)?

# The Paley tournament has UNIFORM Q_k = m for all k
# The Interval tournament has the MOST NON-UNIFORM Q_k

# CF of Q_1/Q_2 for Interval:
print("Continued fraction of Q_1/Q_2 (Interval):")
for p in [7, 11, 13, 17, 23, 29, 37, 41, 47]:
    m = (p-1)//2
    Q = interval_Q(p)
    ratio = Q[0] / Q[1]

    # Compute CF of ratio
    cf = []
    x = ratio
    for _ in range(8):
        a = int(x)
        cf.append(a)
        if abs(x - a) < 1e-10:
            break
        x = 1 / (x - a)

    print(f"  p={p:2d}: Q_1/Q_2 = {ratio:10.4f}, CF = {cf}")

print()

# ============================================================
# PART 3: ZECKENDORF REPRESENTATION
# ============================================================

print("=" * 70)
print("PART 3: ZECKENDORF REPRESENTATION OF H AND RELATED QUANTITIES")
print("=" * 70)
print()
print("Every positive integer has a unique Zeckendorf representation")
print("as a sum of non-consecutive Fibonacci numbers.")
print()

def zeckendorf(n):
    """Compute Zeckendorf representation of n."""
    if n <= 0:
        return []
    fibs = [1, 2]
    while fibs[-1] <= n:
        fibs.append(fibs[-1] + fibs[-2])

    rep = []
    remaining = n
    for f in reversed(fibs):
        if f <= remaining:
            rep.append(f)
            remaining -= f
    return rep

# Zeckendorf of H/p values
H_int = {3: 3, 5: 25, 7: 175, 11: 93027, 13: 3711175, 17: 13689269499}

print("Zeckendorf representations of H/p:")
for p in sorted(H_int.keys()):
    h = H_int[p] // p
    zeck = zeckendorf(h)
    print(f"  p={p:2d}: H/p = {h}, Zeckendorf = {zeck} ({len(zeck)} terms)")

print()

# Zeckendorf of F_p related quantities
print("Zeckendorf of prod(1+Q_k) = F_p:")
for p in [5, 7, 11, 13, 17, 23, 29]:
    Fp = fib(p)
    zeck = zeckendorf(Fp)
    print(f"  p={p:2d}: F_p = {Fp}, Zeckendorf = {zeck}")

print()

# ============================================================
# PART 4: BEATTY SEQUENCES AND TOURNAMENT STRUCTURE
# ============================================================

print("=" * 70)
print("PART 4: BEATTY SEQUENCES — φ AND √2 PARTITIONS")
print("=" * 70)
print()
print("Beatty's theorem: For irrational α > 1, define β s.t. 1/α + 1/β = 1")
print("Then B_α = {floor(nα)} and B_β = {floor(nβ)} partition Z+")
print()
print("For α = φ²: β = φ²/(φ²-1) = φ²·φ/(φ+1)·... = φ²·(φ-1) = φ")
print("Wait, 1/φ² + 1/β = 1 ⟹ β = φ²/(φ²-1) = (3+√5)/2 / ((1+√5)/2)")
print(f"= {phi**2 / (phi**2 - 1):.6f}")
print(f"1/α + 1/β = {1/phi**2 + (phi**2 - 1)/phi**2:.6f}")
print()

# Actually: α = φ, β = φ² (since 1/φ + 1/φ² = (φ+1)/φ² = φ²/φ² = 1)
print(f"For α = φ = {phi:.6f}, β = φ² = {phi**2:.6f}:")
print(f"  1/φ + 1/φ² = {1/phi + 1/phi**2:.6f} = 1 ✓")
print()

# The Beatty sequences for φ and φ²:
B_phi = [int(n * phi) for n in range(1, 20)]
B_phi2 = [int(n * phi**2) for n in range(1, 20)]
print(f"B_φ  = {B_phi}")
print(f"B_φ² = {B_phi2}")
print()

# These partition the natural numbers!
# B_phi = {1,3,4,6,8,9,11,12,14,...} (lower Wythoff sequence)
# B_phi2 = {2,5,7,10,13,15,18,20,...} (upper Wythoff sequence)

# For primes p, does S_int = {1,...,m} relate to a Beatty sequence mod p?
print("Connection to tournament: for prime p, the Interval set is {1,...,m}")
print("Does this relate to Beatty sequences?")
print()

for p in [7, 11, 13, 17, 23]:
    m = (p-1)//2
    # Beatty sequence mod p
    B_mod = sorted(set(int(n * phi) % p for n in range(1, p)))
    B2_mod = sorted(set(int(n * phi**2) % p for n in range(1, p)))
    S_int = list(range(1, m+1))
    overlap = len(set(S_int) & set(B_mod))

    print(f"  p={p}: S_int = {S_int}")
    print(f"    B_φ mod p = {B_mod[:m+2]}... (first {m+2} elements)")
    print(f"    Overlap |S_int ∩ B_φ mod p| = {overlap}/{m}")
    print()

# ============================================================
# PART 5: THE PENROSE TILING CONNECTION
# ============================================================

print("=" * 70)
print("PART 5: PENROSE TILINGS AND TOURNAMENT Q-SPECTRUM")
print("=" * 70)
print()

# Penrose tilings have diffraction pattern with:
# - Peaks at positions related to the Fibonacci lattice
# - Intensities proportional to |sin(πmτ)/sin(πτ)|² where τ = 1/φ
# - This is EXACTLY the Dirichlet kernel at the golden ratio!

# Our Q_k = |D_m(2πk/p)|² = |sin(mπk/p)/sin(πk/p)|²
# For the Fibonacci lattice, the diffraction peaks are at
# k = a + bφ for integers a, b
# And the peak intensity is |D_m(2π(a+bφ))|²

print("The Penrose tiling diffraction pattern has peaks at:")
print("  positions q = h + k·τ where τ = 1/φ")
print("  intensities ∝ |sin(mπτ)/sin(πτ)|²")
print()
print("This is EXACTLY our Q_k formula with k → τ = 1/φ!")
print()

# So: Q_k for the Interval tournament samples the Penrose diffraction
# pattern at the discrete frequencies k/p.

# The peak frequency in the Penrose pattern is at q = 0 (DC)
# The next peak is at q ≈ τ = 1/φ ≈ 0.618...
# For our tournament, k/p ≈ 1/p (small) gives Q_1 ≈ m² (large)
# and k/p ≈ 0.5 gives small Q_k

# Compare: k*/p closest to 1/φ
print("Which k/p is closest to 1/φ = 0.618...?")
for p in [7, 11, 13, 17, 23, 29, 37, 41, 47]:
    m = (p-1)//2
    Q = interval_Q(p)

    # k/p closest to 1/phi
    k_star = round(p / phi)
    if k_star > m:
        k_star = p - k_star  # use the "folded" version

    Q_star = Q[k_star - 1] if k_star <= m else None

    print(f"  p={p:2d}: k*/p = {k_star}/{p} = {k_star/p:.4f}, Q_{k_star} = {Q_star:.4f}" if Q_star else
          f"  p={p:2d}: k*/p = {k_star}/{p} (folded)")

print()

# ============================================================
# PART 6: THE PHYLLOTAXIS / SUNFLOWER CONNECTION
# ============================================================

print("=" * 70)
print("PART 6: PHYLLOTAXIS — BIOLOGICAL FIBONACCI")
print("=" * 70)
print()

# In phyllotaxis, successive leaves/seeds are placed at angles
# separated by the golden angle 2π(2-φ) ≈ 137.5°.
# This creates the optimal packing (no pattern of gaps forms).
#
# For a tournament on p vertices arranged on a circle,
# the Interval tournament S = {1,...,m} connects each vertex to
# its m nearest neighbors (one direction).
# The Paley tournament connects based on quadratic residues —
# a seemingly "random" pattern.
#
# The Fibonacci permutation σ: k ↦ k·g mod p where g ≈ p/φ²
# maps the consecutive Interval set to a "scattered" set.

print("Phyllotaxis analogy:")
print("  Sunflower seeds → tournament vertices (on circle)")
print("  Golden angle → Fibonacci permutation")
print("  Optimal packing → maximum H")
print()
print("  Interval: connect to NEAREST neighbors → CONCENTRATED")
print("  Paley: connect via QR → SCATTERED (Ramanujan expander)")
print()
print("  Phyllotaxis says: golden-angle placement is optimal")
print("  Tournaments say: Interval (consecutive) is optimal for H")
print()
print("  KEY DIFFERENCE: phyllotaxis optimizes PACKING,")
print("  tournaments optimize HAMILTONIAN PATHS")
print("  But both relate to the golden ratio!")
print()

# ============================================================
# PART 7: MARKOV SPECTRUM AND WORST-CASE APPROXIMATION
# ============================================================

print("=" * 70)
print("PART 7: MARKOV SPECTRUM — WORST APPROXIMABLE NUMBERS")
print("=" * 70)
print()

# The Markov spectrum describes how well numbers can be approximated:
# For irrational α, define μ(α) = liminf_{q→∞} q|qα - p|
# The Markov spectrum is the set of all μ(α).
# The worst approximable number is φ with μ(φ) = 1/√5
# Next: √2 with μ = 1/(2√2), then (9+√221)/10 with μ = 1/√221...

# Our Interval tournament has B_m/B_{m-1} → φ²
# And at x=2: B_m(2)/B_{m-1}(2) → (4+√20)/2 = 2+√5 (Pell-related)
# And at x=4: B_m(4)/B_{m-1}(4) → 3+2√2 (silver ratio squared)

print("Metallic means and Morgan-Voyce:")
for x in [1, 2, 3, 4, 5]:
    # B_m(x)/B_{m-1}(x) → fixed point of t = (2+x) - 1/t
    # t² - (2+x)t + 1 = 0
    # t = ((2+x) + sqrt((2+x)²-4))/2 = ((2+x) + sqrt(x²+4x))/2
    disc = x**2 + 4*x
    ratio_limit = ((2+x) + np.sqrt(disc)) / 2

    # What is this algebraically?
    # x=1: (3+√5)/2 = φ²
    # x=2: (4+√12)/2 = 2+√3
    # x=3: (5+√21)/2
    # x=4: (6+√32)/2 = 3+2√2 = (1+√2)² = silver ratio²

    print(f"  x={x}: limit = {ratio_limit:.6f}, disc = {disc}, √disc = {np.sqrt(disc):.6f}")

    # Check if this is a metallic ratio squared
    if x == 1:
        print(f"         = φ² = ((1+√5)/2)² (golden)")
    elif x == 4:
        print(f"         = (1+√2)² = δ² (silver)")
    elif x == 9:
        print(f"         = (3+√13)/2)² (bronze)")

print()

# The sequence: at x = n²-4 (perfect square discriminant):
# x=0: disc=0, degenerate
# x=5: disc=45, not perfect square
# x=12: disc=192, not perfect square
# x=21: disc=525, not perfect square

# For x=n²: disc = n⁴ + 4n² = n²(n²+4), limit = (2+n² + n√(n²+4))/2
# These are NOT metallic ratios in general.

# ============================================================
# PART 8: THE FREE ENERGY CONNECTION
# ============================================================

print("=" * 70)
print("PART 8: FREE ENERGY — THERMODYNAMIC ANALOGY")
print("=" * 70)
print()

# In statistical mechanics, the partition function Z determines
# the free energy F = -kT log Z.
#
# For our system:
# Z_OCF(β) = I(Ω(T), e^β) = sum over independent sets of Ω with activity e^β
# Z_Fib(β) = prod(1 + e^β Q_k)
#
# At β = log(2): Z_OCF = H(T), the Hamiltonian path count
# At β = 0: Z_OCF = number of independent sets of Ω
# As β → ∞: Z_OCF ~ e^{β·α(Ω)} where α = independence number
#
# The "temperature" log(2) is SPECIAL — it's where the physics
# gives us the tournament invariant H.

print("Thermodynamic analogy:")
print("  Lattice gas on Ω(T) at activity z = 2:")
print("  Z(z) = I(Ω, z) = Σ z^|S| over independent sets S")
print("  At z=2: Z = H (Hamiltonian path count)")
print("  At z=1: Z = number of independent sets")
print("  At z→∞: Z ~ z^α(Ω) (dominated by max independent set)")
print()

# The free energy per site:
# f(z) = (1/|Ω|) log Z(z)
# For the Interval tournament, |Ω| = number of odd cycles

# Let's compute I(Ω, z) for various z at p=5
# p=5 is small enough to do exhaustively
p = 5
m = (p-1)//2
S_int = set(range(1, m+1))

# Build tournament
A = [[0]*p for _ in range(p)]
for i in range(p):
    for j in range(p):
        if i != j and (j-i) % p in S_int:
            A[i][j] = 1

# Find 3-cycles (only odd cycles at p=5)
three_cycles = set()
for i in range(p):
    for j in range(p):
        if i == j or not A[i][j]: continue
        for k in range(p):
            if k in (i,j): continue
            if A[j][k] and A[k][i]:
                three_cycles.add(tuple(sorted([i,j,k])))

five_cycles = set()
# At p=5, a 5-cycle uses all vertices
# Check if 0→1→2→3→4→0 exists
for perm in [[0,1,2,3,4]]:  # just check identity
    pass
# Actually enumerate all 5-cycles
from itertools import permutations
for perm in permutations(range(5)):
    if perm[0] != 0: continue  # fix first vertex to avoid counting rotations
    is_cycle = True
    for i in range(5):
        if not A[perm[i]][perm[(i+1)%5]]:
            is_cycle = False
            break
    if is_cycle:
        five_cycles.add(tuple(sorted(perm)))

all_cycles = list(three_cycles) + list(five_cycles)
n_cyc = len(all_cycles)

print(f"p=5 Interval: {len(three_cycles)} 3-cycles, {len(five_cycles)} 5-cycles")
print(f"  Total odd cycles: {n_cyc}")
print()

# Build Ω adjacency
omega_adj = [[0]*n_cyc for _ in range(n_cyc)]
for a in range(n_cyc):
    for b in range(a+1, n_cyc):
        if set(all_cycles[a]) & set(all_cycles[b]):
            omega_adj[a][b] = omega_adj[b][a] = 1

# Compute I(Ω, z) for various z
print(f"I(Ω, z) for p=5 Interval:")
for z in [0.5, 1.0, 1.5, 2.0, 3.0, 5.0, 10.0]:
    I_z = 0
    for mask in range(1 << n_cyc):
        vertices = [i for i in range(n_cyc) if mask & (1 << i)]
        indep = True
        for i in range(len(vertices)):
            for j in range(i+1, len(vertices)):
                if omega_adj[vertices[i]][vertices[j]]:
                    indep = False
                    break
            if not indep: break
        if indep:
            I_z += z ** len(vertices)

    Fp = fib(p)
    Q = interval_Q(p)
    prod_z = np.prod(z + Q) if z > 0 else 0

    print(f"  z={z:5.1f}: I(Ω,z) = {I_z:10.2f}, prod(z+Q) = {prod_z:10.2f}, "
          f"ratio = {I_z/prod_z:.6f}" if prod_z > 0 else
          f"  z={z:5.1f}: I(Ω,z) = {I_z:10.2f}")

print()

# ============================================================
# PART 9: SYNTHESIS — THE KAM-FIBONACCI-TOURNAMENT TRIANGLE
# ============================================================

print("=" * 70)
print("SYNTHESIS: THE KAM-FIBONACCI-TOURNAMENT TRIANGLE")
print("=" * 70)
print()
print("THREE PILLARS connect through the golden ratio:")
print()
print("1. FIBONACCI / NUMBER THEORY:")
print("   F_p = prod(1+Q_k) for Interval tournament")
print("   B_m(x) = Morgan-Voyce polynomial = Chebyshev at golden argument")
print("   Cassini: B_{m-1}·B_{m+1} - B_m² = x")
print()
print("2. KAM / DYNAMICAL SYSTEMS:")
print("   B_m/B_{m-1} → φ² (golden ratio squared, 'most irrational' fixed point)")
print("   Interval Q-spectrum = most concentrated ↔ most 'locked' orbit")
print("   Phase transition at p~13: like onset of chaos")
print("   Metallic ratios at special x: B_m(4)/B_{m-1}(4) → (1+√2)² (silver)")
print()
print("3. STATISTICAL MECHANICS:")
print("   H = I(Ω, 2) = lattice gas partition function at z=2")
print("   F_p = I(∅, Q) = unconstrained partition function")
print("   H ≥ F_p for all circulant tournaments (conjectured)")
print("   Free energy analogy: z=2 is the 'physical temperature'")
print()
print("THE UNIFYING THEME:")
print("The Interval tournament maximizes H because its Q-spectrum is")
print("maximally concentrated, which is the spectral fingerprint of the")
print("golden ratio's irrationality. Just as φ is hardest to approximate")
print("by rationals (most stable in KAM theory), the Interval spectrum")
print("is hardest to 'flatten' (most resistant to perturbation).")
print()
print("CONCRETE PROOF DIRECTION:")
print("The Morgan-Voyce polynomial B_m(x) encodes the Interval spectrum.")
print("By Cassini: B_{m-1}B_{m+1} = B_m² + x, providing a RECURSIVE")
print("stability condition. If we can show this cascade makes the Interval")
print("the unique maximizer of I(Ω, 2) among all circulant tournaments,")
print("we'd have the proof.")
print()
print("OPEN QUESTIONS:")
print("Q1: Is H ≥ F_p a theorem? (Would give: max H ≥ prod(1+Q_k) for Interval)")
print("Q2: Does the Cassini identity B_{m-1}B_{m+1} - B_m² = x give a")
print("    FUNCTIONAL EQUATION for H(p) in terms of H(p-2)?")
print("Q3: Can the KAM small-divisor technique bound the Walsh coefficients?")
print("Q4: Is there a 'Kolmogorov theorem' for tournament optimization?")
print()
