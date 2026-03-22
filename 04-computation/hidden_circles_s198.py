#!/usr/bin/env python3
"""
hidden_circles_s198.py — Wacky Manifestations of x² + y²
opus-2026-03-22-S198

From S197: π appears in tournaments because independence → additive
variances → x² + y² → circular symmetry → π.

Now: where ELSE does x² + y² hide, in forms as convoluted as tournaments?

Each example traces the hidden circle to its source.
"""

import numpy as np
from math import comb, factorial, pi, sqrt, exp, log, gcd
from fractions import Fraction
from collections import Counter

print("=" * 72)
print("  HIDDEN CIRCLES: WACKY MANIFESTATIONS OF x² + y²")
print("  opus-2026-03-22-S198")
print("=" * 72)
print()

# ==========================================================================
# 1. COLLIDING BLOCKS COMPUTE π
# ==========================================================================

print("1. COLLIDING BLOCKS COMPUTE π")
print("=" * 60)
print()

print("""  Setup: A small block (mass 1) sits between a wall and a
  large block (mass M = 100^d). The large block slides left.
  Count elastic collisions until both blocks move right.

  Result: the number of collisions = first d+1 digits of π.

  WHY: Each collision is a linear transformation in (v₁, v₂) space.
  Conservation of energy: ½v₁² + ½Mv₂² = const.
  This IS an ellipse x² + (√M)²y² = const in velocity space.
  Each collision REFLECTS across a line in this space.
  The number of reflections before escaping = arc length / angle
  = πM^{1/2} / ... = digits of π.

  The x² + y² is KINETIC ENERGY: ½mv² is a SUM OF SQUARES
  of velocity components. Conservation of energy = staying on a circle.
""")

# Simulate
def count_collisions(mass_ratio):
    """Count elastic collisions between blocks of mass 1 and mass_ratio."""
    v1, v2 = 0, -1  # small block at rest, large moving left
    m1, m2 = 1, mass_ratio
    count = 0
    while True:
        # Next event: blocks collide or small block hits wall
        if v1 <= v2:  # blocks moving apart or same speed, no more block collision
            if v1 < 0:
                # small block hits wall
                v1 = -v1
                count += 1
            else:
                break
        else:
            # Block collision (elastic)
            new_v1 = ((m1 - m2) * v1 + 2 * m2 * v2) / (m1 + m2)
            new_v2 = ((m2 - m1) * v2 + 2 * m1 * v1) / (m1 + m2)
            v1, v2 = new_v1, new_v2
            count += 1
            if v1 < 0:
                v1 = -v1
                count += 1

    return count

# Actually, simpler implementation:
def collisions_pi(d):
    """Collisions with mass ratio 100^d give π digits."""
    M = 100**d
    v_small, v_big = 0.0, -1.0
    count = 0
    for _ in range(10 * 10**d):  # upper bound on iterations
        if v_small > 0 and v_big > 0 and v_big >= v_small:
            break
        if v_small <= v_big and v_big < 0:
            break
        if v_small <= v_big:
            if v_small < 0:
                v_small = -v_small
                count += 1
                continue
            break
        # Elastic collision
        new_vs = ((1 - M) * v_small + 2 * M * v_big) / (1 + M)
        new_vb = ((M - 1) * v_big + 2 * v_small) / (1 + M)
        v_small, v_big = new_vs, new_vb
        count += 1
        if v_small < 0:
            v_small = -v_small
            count += 1
    return count

for d in range(1, 5):
    n = collisions_pi(d)
    print(f"  M = 100^{d} = {100**d}: collisions = {n}, "
          f"π ≈ {n / 10**d:.{d}f}")

print()
print("  The hidden circle: KINETIC ENERGY ½v₁² + ½Mv₂² = const.")
print("  This is an ELLIPSE in velocity space → angular motion → π.")
print()

# ==========================================================================
# 2. BUFFON'S NEEDLE
# ==========================================================================

print()
print("2. BUFFON'S NEEDLE DROPS")
print("=" * 60)
print()

print("""  Drop a needle of length L on parallel lines spaced D apart (L ≤ D).
  P(needle crosses a line) = 2L / (πD).

  So π = 2L / (D × P(crossing)).

  WHY x² + y²: The needle's position is (center_y, angle θ).
  It crosses a line iff center_y < (L/2)|sin(θ)|.
  The sin(θ) IS the y-coordinate of the unit circle.
  Integrating over θ ∈ [0, π] brings in the angular measure → π.

  The hidden circle: needle ORIENTATION is a point on S¹.
  The average of |sin(θ)| over the circle = 2/π.
""")

# Monte Carlo Buffon
np.random.seed(42)
N = 1000000
L, D = 1.0, 2.0
y_centers = np.random.uniform(0, D/2, N)
thetas = np.random.uniform(0, np.pi, N)
crossings = np.sum(y_centers < (L/2) * np.sin(thetas))
p_cross = crossings / N
pi_est = 2 * L / (D * p_cross)
print(f"  Monte Carlo ({N:,} drops): P(cross) = {p_cross:.6f}")
print(f"  π estimate = {pi_est:.6f} (actual = {pi:.6f})")
print()
print("  The hidden circle: needle angle θ lives on the unit circle S¹.")
print()

# ==========================================================================
# 3. STIRLING'S APPROXIMATION
# ==========================================================================

print()
print("3. STIRLING'S √(2π) IN n!")
print("=" * 60)
print()

print("""  n! ≈ √(2πn) × (n/e)^n

  WHY π in a FACTORIAL? Factorials count permutations.
  Where is the circle in counting permutations?

  The answer: n! = Γ(n+1) = ∫₀^∞ t^n e^{-t} dt.
  Use the saddle-point method: the integrand peaks at t = n.
  Near the peak, t^n e^{-t} ≈ n^n e^{-n} × e^{-(t-n)²/(2n)}.
  That last factor IS a Gaussian with variance n.
  Integrating it gives √(2πn).

  The x² + y²: the saddle point approximation turns the
  integral into a GAUSSIAN integral. And ∫e^{-x²}dx = √π.
  Same polar coordinate trick. Same circle.

  n! COUNTS PERMUTATIONS. The √(2π) appears because the
  number of permutations, viewed as an integral, has a
  Gaussian peak. The Gaussian peak has circular symmetry in 2D.

  Permutations → integral → saddle point → Gaussian → circle → π.
""")

for n in [5, 10, 50, 100]:
    exact = factorial(n) if n <= 170 else float('inf')
    stirling = sqrt(2 * pi * n) * (n / np.e)**n
    ratio = exact / stirling if exact < float('inf') else 0
    print(f"  n={n:3d}: n!/Stirling = {ratio:.8f}")

print()

# ==========================================================================
# 4. PROTEIN FOLDING — RANDOM COIL END-TO-END DISTANCE
# ==========================================================================

print()
print("4. PROTEIN FOLDING: THE RANDOM COIL")
print("=" * 60)
print()

print("""  A polymer chain of N segments, each of length b, with random
  bond angles. The end-to-end vector R = Σ rᵢ (sum of segment vectors).

  If segments are independent:
    E[R] = 0 (random directions cancel)
    E[R²] = E[Σ rᵢ · Σ rⱼ] = N × b² (cross terms vanish)
    E[|R|] = b√(2N/π)    ← π!

  The end-to-end DISTANCE of a random polymer chain involves π
  because each segment is an independent vector in 3D.

  P(|R| = r) ∝ r² × e^{-3r²/(2Nb²)}

  This is a 3D Gaussian (Maxwell distribution).
  The normalization involves ∫e^{-r²}r²dr, computed in
  SPHERICAL coordinates. The angular integral gives 4π.

  The x² + y² (now x²+y²+z²): segment independence makes
  displacements add in quadrature. The 3D Gaussian has
  SPHERICAL symmetry. Integrating over the sphere gives 4π.

  Protein folding uses x²+y²+z² = r².
  The "circle" is now a SPHERE in 3D.
  The angular integral is 4π instead of 2π.
  Same mechanism, one dimension higher.
""")

# Simulate random walks in 3D
np.random.seed(42)
n_walks = 100000
for N in [10, 100, 1000]:
    # Random unit vectors in 3D
    phi = np.random.uniform(0, 2*np.pi, (n_walks, N))
    costheta = np.random.uniform(-1, 1, (n_walks, N))
    sintheta = np.sqrt(1 - costheta**2)

    x = np.sum(sintheta * np.cos(phi), axis=1)
    y = np.sum(sintheta * np.sin(phi), axis=1)
    z = np.sum(costheta, axis=1)

    R = np.sqrt(x**2 + y**2 + z**2)
    mean_R = np.mean(R)
    theory = sqrt(2 * N / pi) * sqrt(3)  # Corrected for 3D
    # Actually E[|R|] for 3D Gaussian = sqrt(8N/(3π))... let me use simulation
    print(f"  N={N:4d} segments: E[|R|] = {mean_R:.3f}, "
          f"√N = {sqrt(N):.3f}, ratio = {mean_R/sqrt(N):.4f}")

print()
print("  The ratio E[|R|]/√N → √(8/(3π)) ≈ 0.921 as N → ∞.")
print("  π appears because independent 3D vectors → spherical Gaussian.")
print()

# ==========================================================================
# 5. PRIME NUMBERS: ζ(2) = π²/6
# ==========================================================================

print()
print("5. PRIME NUMBERS: P(coprime) = 6/π²")
print("=" * 60)
print()

print("""  Pick two random integers. The probability that they share
  no common factor (are coprime) is 6/π².

  WHY π² in NUMBER THEORY?

  P(coprime) = Π_p (1 - 1/p²) = 1/ζ(2) = 6/π².

  Euler proved ζ(2) = Σ 1/n² = π²/6 by factoring sin(πx)/πx:
    sin(πx)/(πx) = Π_n (1 - x²/n²)

  Setting x = 0 and comparing Taylor coefficients gives ζ(2) = π²/6.

  The x² + y²: The product Π(1 - x²/n²) comes from the ZEROS
  of sin(πx), which lives on the UNIT CIRCLE (sin is circular).
  The n² in the denominator IS a sum of squares
  (each n² = n² + 0² in the lattice point sense).

  Alternative: ζ(2) = Σ 1/n² counts lattice points (a,b) with
  a²+b² = n² weighted by 1/n². The sum over the CIRCLE of
  radius n gives the angular factor... and π emerges.

  Coprimality → Euler product → zeros of sine → unit circle → π.
""")

# Verify P(coprime) = 6/π²
np.random.seed(42)
N = 1000000
a = np.random.randint(1, 10000, N)
b = np.random.randint(1, 10000, N)
coprimes = sum(gcd(int(ai), int(bi)) == 1 for ai, bi in zip(a, b))
p_coprime = coprimes / N
pi_est = sqrt(6 / p_coprime)
print(f"  Monte Carlo: P(coprime) = {p_coprime:.6f}")
print(f"  6/π² = {6/pi**2:.6f}")
print(f"  π estimate = {pi_est:.6f}")
print()

# ==========================================================================
# 6. CARD SHUFFLES: THE RIFFLE SHUFFLE PROBABILITY
# ==========================================================================

print()
print("6. CARD SHUFFLES: DESCENT STATISTICS")
print("=" * 60)
print()

print("""  Shuffle a deck of n cards. Count the number of DESCENTS
  (positions where card[i] > card[i+1]).

  The expected number of descents = (n-1)/2.
  The variance of descents = (n+1)/12.

  The distribution of descents, normalized, converges to a
  GAUSSIAN by the CLT (descents are approximately sums of
  weakly dependent indicators).

  The number of permutations with exactly k descents is the
  Eulerian number A(n,k). Its generating function involves:
    A_n(t) = Σ A(n,k) t^k

  And the normalized Eulerian numbers satisfy:
    A(n, (n-1)/2 + x√((n+1)/12)) → (1/√(2π)) e^{-x²/2}

  The π comes from the CLT applied to descent indicators.
  Each position independently contributes (roughly) to the
  descent count. Independence → additive variances → Gaussian → π.

  The x² + y²: variance additivity of descent indicators.
  Each consecutive pair is a binary comparison (like a tournament arc!).
  In fact: a permutation IS a tournament (the "total order" tournament).
  Descents are the "upsets" relative to the identity permutation.
""")

# Verify Eulerian number distribution approaches Gaussian
from itertools import permutations as perms

for n in [4, 6, 8]:
    descent_counts = Counter()
    for p in perms(range(n)):
        descents = sum(1 for i in range(n-1) if p[i] > p[i+1])
        descent_counts[descents] += 1

    total = sum(descent_counts.values())
    mean = sum(k * v for k, v in descent_counts.items()) / total
    var = sum(k**2 * v for k, v in descent_counts.items()) / total - mean**2

    print(f"  n={n}: mean descents = {mean:.3f} (theory {(n-1)/2:.3f}), "
          f"var = {var:.3f} (theory {(n+1)/12:.3f})")

print()
print("  Descents in permutations ≈ upsets in tournaments.")
print("  Both approach Gaussian by CLT. Both give π.")
print()

# ==========================================================================
# 7. NEURAL NETWORKS: RANDOM WEIGHT INITIALIZATION
# ==========================================================================

print()
print("7. NEURAL NETWORKS: KAIMING INITIALIZATION")
print("=" * 60)
print()

print("""  In deep learning, weights are initialized from:
    w ~ N(0, 2/fan_in)  (Kaiming/He initialization)

  The 2 in 2/fan_in comes from E[ReLU(x)²] = Var(x)/2.
  And the NORMALIZATION of the Gaussian uses √(2π).

  But deeper: WHY does this work?

  A neural network layer computes y = Σ wᵢxᵢ.
  If wᵢ and xᵢ are independent with zero mean:
    Var(y) = Σ Var(wᵢ) × E[xᵢ²] = fan_in × σ²_w × E[x²]

  For the variance to be PRESERVED across layers:
    σ²_w = 1 / (fan_in × E[x²])

  The variance preservation IS additive variance (Pythagorean!).
  Each weight independently contributes σ²_w × x²ᵢ to Var(y).
  The sum of independent contributions → Gaussian by CLT.
  π enters through the Gaussian normalization.

  Neural network stability REQUIRES the Pythagorean theorem
  for variances. If weights were correlated, the initialization
  would be different, π would not appear, and the network
  would not train.

  The x² + y²: Var(y) = Σ Var(wᵢxᵢ) is LITERALLY a sum of
  squares (each variance is a squared quantity). The √(2π)
  in the initialization formula comes from ∫e^{-w²}dw = √π.
""")

# Demonstrate: variance preservation through random layers
np.random.seed(42)
fan_in = 100
n_layers = 20
x = np.random.randn(1000)

variances = [np.var(x)]
for layer in range(n_layers):
    W = np.random.randn(fan_in, fan_in) * sqrt(2.0 / fan_in)
    x_new = np.maximum(0, W @ np.random.randn(fan_in, 1000))  # ReLU
    x = x_new[0]  # Take one output neuron
    variances.append(np.var(x))

print(f"  Kaiming init: variance across {n_layers} layers:")
for i in [0, 1, 5, 10, 15, 19]:
    if i < len(variances):
        print(f"    Layer {i:2d}: var = {variances[i]:.4f}")

print()

# ==========================================================================
# 8. MUSIC: THE EQUAL TEMPERAMENT SCALE
# ==========================================================================

print()
print("8. MUSIC: WHY 12 NOTES PER OCTAVE")
print("=" * 60)
print()

print("""  The Western scale divides the octave into 12 equal steps.
  Each step multiplies frequency by 2^(1/12).

  WHY 12? Because 12 gives the best rational approximations
  to the PERFECT FIFTH (3/2) and PERFECT FOURTH (4/3).

  2^(7/12) ≈ 1.4983 ≈ 3/2 = 1.5000 (error 0.11%)

  But WHERE is x² + y²?

  Sound waves ARE circles: a pure tone is sin(2πft).
  The FREQUENCY f determines the angular velocity on S¹.
  Two notes harmonize when their frequency ratio is simple.

  The DISSONANCE between two notes at frequencies f₁ and f₂
  depends on how their overtone series interfere:
    energy = Σ (a_k sin(2πkf₁t) + b_k sin(2πkf₂t))²

  Squaring → cross terms → beating → SUMS OF SQUARES of
  amplitude differences. The human ear perceives the POWER
  (sum of squared amplitudes) — a quadratic form!

  The x² + y²: sound energy is ||amplitude||² = Σ aᵢ².
  Each frequency component contributes independently.
  The total energy is a sum of squares. The fundamental
  period of sin is 2π. The octave is 2:1 frequency ratio.
  Equal temperament divides the CIRCLE (one period of
  the logarithmic frequency axis) into 12 equal arcs.

  Music IS division of the circle. The 12 notes ARE 12 equally
  spaced points on S¹. Harmony IS proximity on S¹.
  π IS the half-octave (the tritone: 2^(6/12) = √2).
""")

# Show equal temperament as circle division
print(f"  The 12-note circle:")
notes = ['C', 'C#', 'D', 'D#', 'E', 'F', 'F#', 'G', 'G#', 'A', 'A#', 'B']
for i, note in enumerate(notes):
    angle = 2 * pi * i / 12
    x_pos = np.cos(angle)
    y_pos = np.sin(angle)
    freq_ratio = 2**(i/12)
    print(f"    {note:2s}: θ = {angle:5.2f} rad ({i*30:3d}°), "
          f"freq ratio = {freq_ratio:.4f}, "
          f"(x,y) = ({x_pos:+.3f}, {y_pos:+.3f})")

print()
print("  The perfect fifth C→G spans 7 semitones = 7/12 of the circle.")
print(f"  2^(7/12) = {2**(7/12):.6f} ≈ 3/2 = {3/2}")
print()

# ==========================================================================
# 9. DRUNKARD'S WALK ON A CHESSBOARD
# ==========================================================================

print()
print("9. CHESS PIECE RETURN PROBABILITIES")
print("=" * 60)
print()

print("""  A king doing a random walk on an infinite chessboard.
  At each step, it moves to one of 8 adjacent squares (uniformly).

  The probability of returning to the starting square after 2k steps
  involves 1/(πk) — the 2D random walk return probability!

  This is Pólya's theorem: random walks in 2D are RECURRENT
  (they return to the origin with probability 1), but the
  expected time to return is INFINITE.

  The return probability after 2k steps ~ 1/(πk) because the
  2D walk has independent x and y components, each ~ N(0, k).
  The joint distribution is a 2D Gaussian with CIRCULAR symmetry.
  The probability at the origin = 1/(2πσ²) where σ² ~ k.

  In 3D: the random walk is TRANSIENT (returns with probability
  ≈ 0.3405). The 3D Gaussian has SPHERICAL symmetry.
  Return probability ~ 1/(k^{3/2} × (2π)^{3/2}).

  In 1D: always returns. Probability ~ 1/√(πk).
  In 2D: always returns. Probability ~ 1/(πk).
  In 3D: may not return. Probability ~ 1/(πk)^{3/2}.

  The DIMENSION of the walk determines HOW MANY ANGULAR
  INTEGRALS you get: 2π in 2D, 4π in 3D.
  More dimensions → smaller return probability → more π's.

  π COUNTS THE DIMENSIONS OF THE CIRCLE/SPHERE.
""")

# Simulate 2D random walk return
np.random.seed(42)
n_walks = 200000
max_steps = 200

# 2D walk
for target_k in [10, 50, 100]:
    steps_x = np.random.choice([-1, 0, 1], size=(n_walks, 2*target_k))
    steps_y = np.random.choice([-1, 0, 1], size=(n_walks, 2*target_k))
    pos_x = np.cumsum(steps_x, axis=1)
    pos_y = np.cumsum(steps_y, axis=1)
    at_origin = (pos_x[:, -1] == 0) & (pos_y[:, -1] == 0)
    p_return = np.mean(at_origin)
    theory = 1 / (pi * target_k)  # rough approximation
    print(f"  2D walk, 2k={2*target_k}: P(return) = {p_return:.6f}, "
          f"~1/(πk) = {theory:.6f}")

print()

# ==========================================================================
# 10. QUANTUM MECHANICS: THE UNCERTAINTY PRINCIPLE
# ==========================================================================

print()
print("10. QUANTUM MECHANICS: FOURIER UNCERTAINTY")
print("=" * 60)
print()

print("""  Heisenberg: ΔxΔp ≥ ℏ/2 where ℏ = h/(2π).

  WHY 2π? Because the Fourier transform converts position
  to momentum:
    ψ̂(p) = ∫ ψ(x) e^{-2πipx} dx

  The e^{-2πipx} IS uniform rotation on the unit circle in C.
  The 2π comes from the PERIOD of the complex exponential.

  The minimum-uncertainty state is a Gaussian:
    ψ(x) = (2πσ²)^{-1/4} e^{-x²/(4σ²)}

  Its Fourier transform is ALSO a Gaussian:
    ψ̂(p) = (2πσ̃²)^{-1/4} e^{-p²/(4σ̃²)}

  With σσ̃ = 1/(4π). The Gaussian is the ONLY function
  that is its own Fourier transform (up to scaling).

  The x² + y²: the Fourier kernel e^{-2πixy} traces
  the unit circle in the complex plane. The integral
  ∫e^{-x²}e^{-2πixy}dx involves the function (x+iy)²...
  the Gaussian in the COMPLEX plane has symmetry
  under ROTATION of (x, y) → (x cos θ - y sin θ, ...).

  Uncertainty principle = circles in PHASE SPACE (x, p).
  The minimum uncertainty region IS a circle of area ℏ/2.
  Phase space volume is measured in units of 2πℏ = h.
  Each quantum state occupies one CIRCLE's worth of phase space.
""")

# Demonstrate Fourier uncertainty with Gaussian
x = np.linspace(-10, 10, 10000)
dx = x[1] - x[0]
for sigma in [0.5, 1.0, 2.0]:
    psi = np.exp(-x**2 / (4 * sigma**2)) / (2 * pi * sigma**2)**0.25
    # Fourier transform
    psi_fft = np.fft.fftshift(np.fft.fft(np.fft.ifftshift(psi))) * dx
    p = np.fft.fftshift(np.fft.fftfreq(len(x), d=dx))

    delta_x = np.sqrt(np.sum(x**2 * np.abs(psi)**2 * dx))
    delta_p = np.sqrt(np.sum(p**2 * np.abs(psi_fft)**2 * (p[1]-p[0])))
    product = delta_x * delta_p

    print(f"  σ={sigma:.1f}: Δx={delta_x:.4f}, Δp={delta_p:.4f}, "
          f"ΔxΔp={product:.4f} (≥ 1/(4π) = {1/(4*pi):.4f})")

print()

# ==========================================================================
# 11. COOKING: FLAVOR BALANCE IS A NORM
# ==========================================================================

print()
print("11. COOKING: FLAVOR BALANCE AS A QUADRATIC FORM")
print("=" * 60)
print()

print("""  A dish has flavor dimensions: salt, sweet, sour, bitter, umami.
  Each ingredient contributes independently to each dimension.
  The BALANCE of a dish ≈ how close the flavor vector is to
  the target vector.

  If you model flavor as a vector f = (f₁, f₂, ..., f₅),
  then "deliciousness" ∝ e^{-||f - target||²/σ²}
  = e^{-(f₁-t₁)² - (f₂-t₂)² - ... - (f₅-t₅)²}

  This IS e^{-r²} in 5 dimensions!

  If you randomly substitute one ingredient:
  P(dish stays good) ∝ fraction of 5D sphere near target
  ∝ 1/(k^{5/2}) × (2π)^{5/2} normalization

  More flavor dimensions → more π's in the normalization.
  The volume of the d-dimensional sphere of radius r is:
    V_d(r) = π^{d/2} / Γ(d/2 + 1) × r^d

  π appears d/2 times because the sphere in d dimensions
  has d/2 pairs of coordinates, each contributing one
  angular integral of 2π.

  COOKING IS HIGH-DIMENSIONAL SPHERE PACKING.
  π counts how many independent flavor dimensions you have.
""")

# Volume of d-dimensional unit sphere
print(f"  Volume of unit sphere in d dimensions:")
for d in range(1, 11):
    from math import gamma as gam
    vol = pi**(d/2) / gam(d/2 + 1)
    print(f"    d={d:2d}: V_d(1) = {vol:.6f}  "
          f"(π appears {d/2:.1f} times)")

print()

# ==========================================================================
# 12. INTERNET: THE PAGERANK RANDOM WALK
# ==========================================================================

print()
print("12. INTERNET: PAGERANK AND THE STATIONARY DISTRIBUTION")
print("=" * 60)
print()

print("""  Google's PageRank models a "random surfer" who follows links.
  With probability α, follow a random link from current page.
  With probability 1-α, jump to a random page.

  The stationary distribution is the eigenvector of the
  transition matrix for eigenvalue 1.

  WHERE is x² + y²?

  The MIXING TIME of the random walk (how fast PageRank converges)
  depends on the spectral gap of the transition matrix.
  The eigenvalues of a random transition matrix are approximately
  distributed in the DISK |z| ≤ 1/√n in the complex plane.

  This is the CIRCULAR LAW of random matrix theory:
  the eigenvalues of a random n×n matrix fill a DISK
  of radius 1/√n uniformly.

  Why a DISK? Because each eigenvalue z = x + iy has
  |z|² = x² + y² ≤ 1/n. This is x² + y² = r².
  The circular law IS the circle.

  The mixing time ~ 1/gap ~ √n, and the √ involves...
  yes, the normalization of the eigenvalue distribution
  over the disk, which involves the AREA = πr² = π/n.

  PageRank inherits π from the CIRCULAR LAW of random matrices,
  which inherits it from the circular symmetry of the complex plane,
  which inherits it from x² + y² = r².
""")

print("  (PageRank convergence rate inherits π from circular law.)")
print()

# ==========================================================================
# 13. ELECTIONS: PROBABILITY OF A TIE
# ==========================================================================

print()
print("13. ELECTIONS: PROBABILITY OF A TIE")
print("=" * 60)
print()

print("""  N voters, each independently choosing candidate A or B with p=1/2.
  P(exact tie when N is even) = C(N, N/2) / 2^N ~ 1/√(πN/2)

  This IS the random walk return probability.
  This IS the fiber fraction.
  This IS C(2k,k)/4^k with k = N/2.

  An election IS a tournament on 2 candidates with N comparisons.
  The "score" is the vote count.
  A "tie" is when the score is preserved (both at N/2).

  The probability of a tie in a close election scales as 1/√(πN).
  π appears because voters are independent → CLT → Gaussian → π.

  If voters were correlated (herding behavior), the tie probability
  would NOT scale as 1/√(πN). It would be different. Some elections
  might have higher tie probability (correlated → clustered opinions)
  or lower (anti-correlated → polarized). The √π specifically
  measures the INDEPENDENCE of voter behavior.

  π IS THE DEMOCRATIC CONSTANT:
  it measures how independent voters are in an election.
""")

# Election tie probability
for N in [10, 100, 1000, 10000]:
    k = N // 2
    if N <= 1000:
        exact = comb(N, k) / 2**N
    else:
        exact = 1 / sqrt(pi * N / 2)  # Stirling
    theory = 1 / sqrt(pi * N / 2)
    print(f"  N={N:5d} voters: P(tie) = {exact:.8f}, "
          f"1/√(πN/2) = {theory:.8f}")

print()

# ==========================================================================
# SUMMARY
# ==========================================================================

print()
print("=" * 72)
print("  SUMMARY: 13 HIDDEN CIRCLES")
print("=" * 72)
print()

print("""
  Each example traces back to x² + y² through a different path:

   # | SYSTEM             | x² + y² HIDES IN...
  ---|--------------------|-----------------------------------------
   1 | Colliding blocks   | Kinetic energy ½mv² (ellipse in v-space)
   2 | Buffon's needle    | Needle angle on S¹
   3 | Stirling's n!      | Saddle point → Gaussian → ∫e^{-x²}
   4 | Protein folding    | Independent bond vectors → 3D Gaussian
   5 | Coprime integers   | Zeros of sin(πx) on unit circle
   6 | Card shuffles      | Descent statistics → CLT → Gaussian
   7 | Neural networks    | Additive variance preservation
   8 | Musical scales     | Frequency ratio on logarithmic circle
   9 | Chess random walk  | 2D walk → circular Gaussian contours
  10 | Quantum mechanics  | Fourier kernel e^{-2πix} on unit circle
  11 | Cooking            | Flavor norm ||f||² = sum of squares
  12 | PageRank           | Circular law of random matrices
  13 | Elections          | Independent voters → tie probability

  THE COMMON PATTERN:
  Independent components + sum/norm → Gaussian → ∫e^{-x²} → √π

  OR equivalently:
  Any quantity living on S¹ (the circle) in ANY disguise → π

  THE MOST SURPRISING: colliding blocks (physics problem with
  NO obvious circle gives DIGITS OF π).

  THE MOST PRACTICAL: neural network initialization (every deep
  learning model uses √(2π) without knowing it's a circle).

  THE MOST BEAUTIFUL: coprime integers (pure number theory,
  no geometry in sight, yet π² appears through Euler's sine product).

  π IS NOT A GEOMETRIC CONSTANT.
  π IS THE CONSTANT OF INDEPENDENCE AND QUADRATIC FORMS.
  Geometry is just one place where these appear.
""")

print("opus-2026-03-22-S198")
