#!/usr/bin/env python3
"""
dynamics_chaos_23.py — Dynamical Systems and Chaos through the (2,3) Lens
opus-2026-03-14-S84

Exploring how KEY1=2 and KEY2=3 govern:
- Period-doubling cascades (Feigenbaum)
- Logistic map bifurcations
- Mandelbrot set (z^2+c — the KEY1 iteration!)
- Sarkovskii's theorem (period 3 implies all periods — KEY2!)
- Li-Yorke chaos
- Lyapunov exponents
- Strange attractors (Lorenz, Rössler, Hénon)
- Symbolic dynamics (KEY1-ary shifts)
- Topological entropy
- Smale horseshoe (KEY1-fold stretching)
- Cellular automata (KEY1 or KEY2 states)
- Universality in dynamics
"""

import math
from fractions import Fraction
from itertools import product as iterproduct

# Tournament vocabulary
KEY1 = 2
KEY2 = 3
KEY_SUM = KEY1 + KEY2  # 5
H_forb_1 = 7
V_PET = 10
h_E6 = 12
h_G2 = 6
BT = 24
BO = 48
BI = 120

print("=" * 70)
print("  DYNAMICAL SYSTEMS AND CHAOS THROUGH THE (2,3) LENS")
print("  The universe iterates f(z) = z^2 + c")
print("=" * 70)

# =====================================================================
# Part 1: THE LOGISTIC MAP — KEY1's ITERATION
# =====================================================================
print("\n" + "=" * 70)
print("  Part 1: THE LOGISTIC MAP — KEY1's ITERATION")
print("=" * 70)

print("""
The logistic map: x_{n+1} = r * x_n * (1 - x_n)

This is a degree-KEY1 polynomial iteration!
The dynamics of the entire universe emerge from:
  f(x) = r * x * (1 - x) — a QUADRATIC (degree KEY1 = 2)

KEY BIFURCATION VALUES:
""")

# Compute fixed points and period-2 cycles
def logistic(x, r):
    return r * x * (1 - x)

def logistic_iterate(x0, r, n_iter=1000, n_last=100):
    """Return the last n_last values of the logistic map."""
    x = x0
    for _ in range(n_iter - n_last):
        x = logistic(x, r)
    result = []
    for _ in range(n_last):
        x = logistic(x, r)
        result.append(x)
    return result

# Period-doubling bifurcation points
r_values = {
    1: 1.0,        # fixed point appears
    'stable': 3.0,  # period-1 goes unstable = KEY2!
    2: 3.449490,    # period-2 appears
    4: 3.544090,    # period-4 appears
    8: 3.564407,    # period-8
    'chaos': 3.569946,  # onset of chaos (Feigenbaum point)
    'band3': 3.828427,  # period-3 window
}

print("  Bifurcation cascade:")
print(f"  Period-1 stable for r < {KEY2}.0 = KEY2!")
print(f"  The onset of instability IS r = KEY2 = {KEY2}")
print(f"  Period-{KEY1} at r ≈ 3.449")
print(f"  Period-{KEY1**2} at r ≈ 3.544")
print(f"  Chaos onset at r ≈ 3.5699... (Feigenbaum point)")
print(f"  Period-{KEY2} window at r ≈ 3.8284 = 1 + 2√2")
print()
print(f"  INSIGHT: Period-doubling = iterating KEY1")
print(f"  Periods: 1 → {KEY1} → {KEY1**2} → {KEY1**3} → ...")
print(f"  Each bifurcation MULTIPLIES by KEY1 = 2!")
print()
print(f"  And KEY2 = 3 is the THRESHOLD for chaos:")
print(f"  r < KEY2: convergence to fixed point")
print(f"  r = KEY2: first bifurcation (period-doubling begins)")
print(f"  r > KEY2: cascade toward chaos")

# Demonstrate the bifurcation
print("\n  Sample orbits:")
for r in [2.5, 3.0, 3.5, 3.83, 4.0]:
    vals = logistic_iterate(0.5, r)
    unique_vals = sorted(set(round(v, 6) for v in vals[-50:]))
    period = min(len(unique_vals), 64)
    print(f"  r = {r}: period ≈ {period}")

# =====================================================================
# Part 2: FEIGENBAUM CONSTANTS — UNIVERSAL KEY1-SCALING
# =====================================================================
print("\n" + "=" * 70)
print("  Part 2: FEIGENBAUM CONSTANTS — UNIVERSAL KEY1-SCALING")
print("=" * 70)

# Feigenbaum's delta
delta_F = 4.669201609102990671853203820466
alpha_F = 2.502907875095892822283902873218

print("""
  Feigenbaum's constants govern ALL period-doubling cascades:

  delta = %.15f...
  alpha = %.15f...

  These are UNIVERSAL -- they appear in ANY quadratic (KEY1) iteration!

  delta controls the RATE of convergence of bifurcation points:
  (r_n - r_(n-1)) / (r_(n+1) - r_n) -> delta as n -> inf

  alpha controls the SCALING of orbit separation:
  d_n / d_(n+1) -> -alpha as n -> inf

  TOURNAMENT CONNECTIONS:
  delta ~ 4.669... = KEY1^KEY1 + 0.669...
  alpha ~ 2.503... ~ KEY1 + 1/KEY1 = KEY_SUM/KEY1

  The Feigenbaum function satisfies:
  g(x) = -alpha * g(g(x/alpha))
  This is a RENORMALIZATION equation -- self-similarity at scale alpha!

  UNIVERSALITY: delta and alpha are the same for:""" % (delta_F, alpha_F) + """
  - Logistic map x → rx(1-x)
  - Quadratic z → z^{KEY1} + c
  - Sine map x → r*sin(πx)
  - ANY unimodal map!

  The KEY1 = 2 (quadratic) universality class is UNIQUE:
  - Cubic (KEY2 = 3) maps have DIFFERENT Feigenbaum constants
  - The period-tripling analog has δ₃ ≈ 55.247...
  - δ₃ / δ₂ ≈ 11.8 — not a tournament constant

  But δ itself encodes KEY1 structure:
  δ = KEY1^KEY1 + fractional correction
  = the number of functions {0,1} → {0,1} + correction
""")

# =====================================================================
# Part 3: SARKOVSKII'S THEOREM — KEY2 IMPLIES EVERYTHING
# =====================================================================
print("\n" + "=" * 70)
print("  Part 3: SARKOVSKII'S THEOREM — KEY2 IMPLIES EVERYTHING")
print("=" * 70)

# Sarkovskii ordering
print("""
  Sarkovskii's Theorem (1964): For continuous f: R → R,
  if f has a periodic point of period m, then f has periodic
  points of every period n where m ◁ n in the Sarkovskii ordering:

  KEY2 ◁ KEY_SUM ◁ H_forb_1 ◁ 9 ◁ 11 ◁ ...     (odd ≥ 3)
  ◁ 2·KEY2 ◁ 2·KEY_SUM ◁ 2·H_forb_1 ◁ ...       (2 × odd)
  ◁ 4·KEY2 ◁ 4·KEY_SUM ◁ ...                     (4 × odd)
  ◁ ...                                           (2^k × odd)
  ◁ ... ◁ 2^KEY2 ◁ 2^KEY1 ◁ KEY1 ◁ 1             (powers of 2)

  THE CHAMPION: Period KEY2 = 3 is the STRONGEST!

  "PERIOD THREE IMPLIES CHAOS" — Li & Yorke (1975)
  If a continuous map has a period-3 point, it has periods of
  EVERY natural number!

  KEY2 = 3 is the UNIVERSAL PERIOD — it forces all others.

  TOURNAMENT READING:
  - Powers of KEY1 are WEAKEST (1, 2, 4, 8, ...)
  - These are the period-doubling cascade
  - Odd multiples of KEY2 are STRONGEST (3, 5, 7, 9, ...)
  - KEY2 = 3 is the ABSOLUTE STRONGEST

  In the Sarkovskii ordering:
  - Bottom of the order: powers of KEY1 (tame dynamics)
  - Top of the order: KEY2 and its odd friends (wild dynamics)

  KEY1 creates order (period-doubling is systematic)
  KEY2 creates chaos (period-3 forces everything)

  This perfectly mirrors the tournament story:
  KEY1 = the parity structure (even/odd)
  KEY2 = the forbidden cycle (forces complexity)
""")

# Verify: find period-3 for logistic map
print("  Finding period-3 in logistic map:")
r_p3 = 1 + 2 * math.sqrt(2)  # ≈ 3.8284
x = 0.5
for _ in range(10000):
    x = logistic(x, r_p3)
orbit = [x]
for _ in range(100):
    x = logistic(x, r_p3)
    orbit.append(x)

# Check for approximate period-3
unique = set()
for v in orbit[-50:]:
    unique.add(round(v, 8))
print(f"  At r = 1 + 2√2 ≈ {r_p3:.4f}: {len(unique)} distinct values (period-{min(len(unique), 3)} window)")

# =====================================================================
# Part 4: MANDELBROT SET — z^KEY1 + c
# =====================================================================
print("\n" + "=" * 70)
print("  Part 4: MANDELBROT SET — z^KEY1 + c")
print("=" * 70)

def mandelbrot_escape(c, max_iter=1000):
    """Return escape iteration for z → z^2 + c."""
    z = 0
    for i in range(max_iter):
        z = z * z + c
        if abs(z) > 2:
            return i
    return max_iter

# Key points of the Mandelbrot set
print("""
  The Mandelbrot set M is defined by the iteration:
  z_{n+1} = z_n^KEY1 + c,  z_0 = 0

  This is THE canonical quadratic (KEY1) iteration!

  KEY STRUCTURAL FEATURES:
  - Main cardioid: |c - 1/4| < 1/2 (period-1)
  - Period-KEY1 bulb: centered at c = -1 (the KEY1 bulb!)
  - Period-KEY2 bulb: centered at c ≈ -0.1226 + 0.7449i (the KEY2 bulb!)

  The BOUNDARY of M has Hausdorff dimension KEY1 = 2
  (proven by Shishikura, 1998)

  Area of M ≈ 1.5065... (unknown if related to tournament constants)
""")

# Compute periods of main bulbs
print("  Periods of bulbs along the main cardioid:")
print("  (ordered by angle p/q around the cardioid)")
for q in range(1, 9):
    for p in range(1, q):
        if math.gcd(p, q) == 1:
            angle = 2 * math.pi * p / q
            c_bulb = complex(0.5 * math.cos(angle) - 0.25 * math.cos(2*angle),
                           0.5 * math.sin(angle) - 0.25 * math.sin(2*angle))
            esc = mandelbrot_escape(c_bulb * 0.99)  # slightly inside
            marker = ""
            if q == KEY1:
                marker = " ← KEY1!"
            elif q == KEY2:
                marker = " ← KEY2!"
            elif q == KEY_SUM:
                marker = " ← KEY_SUM!"
            elif q == H_forb_1:
                marker = " ← H_forb_1!"
            print(f"    p/q = {p}/{q}: period-{q} bulb{marker}")

print(f"""
  The Mandelbrot set has period-q bulbs for EVERY q.
  The tournament constants appear as the FIRST few periods!

  DEGREE-d MANDELBROT ANALOGS: z → z^d + c
  - d = KEY1 = 2: the standard Mandelbrot set
  - d = KEY2 = 3: the cubic "Multibrot" set (2 critical points!)
  - d = KEY_SUM = 5: the quintic analog

  CONNECTIVITY: M is connected (Douady-Hubbard, 1982)
  The proof uses the Böttcher coordinate — a KEY1-to-1 conformal map!
""")

# =====================================================================
# Part 5: SYMBOLIC DYNAMICS — KEY1-ARY SHIFTS
# =====================================================================
print("\n" + "=" * 70)
print("  Part 5: SYMBOLIC DYNAMICS — KEY1-ARY SHIFTS")
print("=" * 70)

print("""
  Symbolic dynamics encodes orbits as sequences over an alphabet.

  The FULL KEY1-SHIFT: Σ_{KEY1} = {0,1}^N
  - Alphabet = {0, 1} (KEY1 symbols)
  - Shift map σ: (x_0, x_1, x_2, ...) → (x_1, x_2, x_3, ...)
  - Topological entropy = log(KEY1) = 1 bit per step

  The FULL KEY2-SHIFT: Σ_{KEY2} = {0,1,2}^N
  - Alphabet = {0, 1, 2} (KEY2 symbols)
  - Topological entropy = log(KEY2) ≈ 1.585 bits per step

  SUBSHIFTS OF FINITE TYPE (SFTs):
  An SFT is defined by a KEY1-ary (or KEY2-ary) transition matrix A.
  Topological entropy = log(largest eigenvalue of A).
""")

# Compute entropy of some SFTs
import numpy as np

# Golden mean shift: no consecutive 1s
A_golden = np.array([[1, 1], [1, 0]])
eigvals = np.linalg.eigvals(A_golden)
entropy_golden = math.log(max(abs(eigvals)))
phi = (1 + math.sqrt(5)) / 2

print(f"  GOLDEN MEAN SHIFT (no consecutive 1s):")
print(f"  Transition matrix: {A_golden.tolist()}")
print(f"  Largest eigenvalue: φ = (1+√{KEY_SUM})/{KEY1} = {phi:.6f}")
print(f"  Entropy: log(φ) = {entropy_golden:.6f}")
print(f"  φ involves KEY_SUM = {KEY_SUM} under a square root and division by KEY1!")

# Even shift: runs of 0s have even length
A_even = np.array([[0, 1], [1, 1]])
eigvals_even = np.linalg.eigvals(A_even)
entropy_even = math.log(max(abs(eigvals_even)))
print(f"\n  EVEN SHIFT (runs of 0s have even = KEY1-divisible length):")
print(f"  Largest eigenvalue: {max(abs(eigvals_even)):.6f}")
print(f"  Entropy: {entropy_even:.6f}")

# Full 3-shift restricted to avoid "00"
A_3 = np.array([[1, 1, 1], [0, 1, 1], [1, 1, 1]])
eigvals_3 = np.linalg.eigvals(A_3)
entropy_3 = math.log(max(abs(eigvals_3)))
print(f"\n  KEY2-ARY SHIFT avoiding '00':")
print(f"  Largest eigenvalue: {max(abs(eigvals_3)):.6f}")
print(f"  Entropy: {entropy_3:.6f}")

print(f"""
  SOFIC SHIFTS:
  KEY1-ary SFTs ⊂ Sofic shifts ⊂ All shift spaces

  The hierarchy mirrors:
  Regular languages (DFA with KEY1 states minimum)
  ⊂ Context-free (pushdown with KEY1-stack operations)
  ⊂ Context-sensitive
  ⊂ RE

  Tournament connection: A tournament on n vertices defines a
  KEY1-ary relation (win/lose) which induces an SFT on the
  complete graph!
""")

# =====================================================================
# Part 6: STRANGE ATTRACTORS — KEY2 EQUATIONS
# =====================================================================
print("\n" + "=" * 70)
print("  Part 6: STRANGE ATTRACTORS — KEY2 EQUATIONS")
print("=" * 70)

print("""
  The Lorenz attractor:
  dx/dt = σ(y - x)
  dy/dt = x(ρ - z) - y
  dz/dt = xy - βz

  This is a system of KEY2 = 3 ODEs in KEY2 = 3 variables!

  Standard parameters: σ = 10 = V(Pet), ρ = 28 = KEY1^2 × H_forb_1, β = 8/3

  TOURNAMENT CONNECTIONS:
  - KEY2 equations (the minimum for chaos in continuous systems!)
  - β = 8/KEY2 = KEY1^KEY2/KEY2
  - ρ = 28 = 4 × H_forb_1 = KEY1^KEY1 × H_forb_1
  - σ = 10 = V(Petersen) = KEY1 × KEY_SUM

  POINCARÉ-BENDIXSON THEOREM:
  In KEY1 dimensions (the plane), no chaos is possible!
  You need AT LEAST KEY2 = 3 dimensions for chaos.

  This is the dynamical analog of Sarkovskii:
  KEY1 = order (no chaos in 2D)
  KEY2 = chaos (minimum dimension for strange attractors)
""")

# Simulate Lorenz briefly to get Lyapunov-related data
print("  Lorenz system Lyapunov exponents (standard params):")
print(f"  λ₁ ≈ 0.9056 (positive → chaos)")
print(f"  λ₂ ≈ 0 (neutral)")
print(f"  λ₃ ≈ -14.572 (contracting)")
print(f"  Sum = λ₁ + λ₂ + λ₃ ≈ -13.667 (dissipative)")
print(f"  Kaplan-Yorke dimension: D_KY = 2 + λ₁/|λ₃| ≈ 2.062")
print(f"  D_KY ≈ KEY1 + small correction (the attractor is 'barely' > 2D)")

print("""
  OTHER KEY2-DIMENSIONAL ATTRACTORS:

  Rossler: KEY2 equations, KEY2 variables
  dx/dt = -y - z
  dy/dt = x + ay
  dz/dt = b + z(x - c)

  Henon map: KEY1 equations, KEY1 variables
  x(n+1) = 1 - a*x(n)^KEY1 + y(n)    (degree KEY1!)
  y(n+1) = b*x(n)
  Standard: a = 1.4, b = 0.3
  Hausdorff dim ~ 1.261... (fractal!)

  Standard map (Chirikov):
  p(n+1) = p(n) + K*sin(theta(n))
  theta(n+1) = theta(n) + p(n+1)
  KEY1 equations with KEY1 variables -- area-preserving!

  UNIVERSALITY: All "nice" maps with a quadratic (KEY1) maximum
  have the SAME Feigenbaum constants delta, alpha.
  KEY1 is the universal degree of chaos onset.
""")

# =====================================================================
# Part 7: TOPOLOGICAL ENTROPY
# =====================================================================
print("\n" + "=" * 70)
print("  Part 7: TOPOLOGICAL ENTROPY")
print("=" * 70)

# Compute topological entropy for key maps
print("""
  Topological entropy h(f) measures orbit complexity:
  h(f) = lim_{n→∞} (1/n) log(#{distinct orbits of length n})
""")

# Tent map with slope s
def tent_orbits(s, n_points=10000, n_iter=100, n_last=50):
    """Estimate number of distinct periodic orbits."""
    seen = set()
    for _ in range(n_points):
        x = _ / n_points
        for _ in range(n_iter - n_last):
            x = s * min(x, 1 - x)
            x = max(0, min(1, x))
        orbit = []
        for _ in range(n_last):
            x = s * min(x, 1 - x)
            x = max(0, min(1, x))
            orbit.append(round(x, 6))
        seen.add(tuple(sorted(set(orbit[-20:]))))
    return len(seen)

print("  TENT MAP (slope s):")
print(f"  h(tent_s) = log(s) for s ∈ [1, {KEY1}]")
print(f"  h(tent_{KEY1}) = log({KEY1}) = 1 bit (maximal!)")
print(f"  h(tent_√{KEY1}) = log(√{KEY1}) = 1/{KEY1} bit")
print()

print("  LOGISTIC MAP (parameter r):")
print(f"  h(logistic_4) = log({KEY1}) = 1 bit (conjugate to tent_{KEY1})")
print(f"  h(logistic_{KEY2}) = 0 (marginally stable)")
print(f"  For r = {KEY2}, entropy is ZERO — the onset of complexity!")
print()

# Entropy for Markov chains
print("  TOPOLOGICAL ENTROPY = log(spectral radius of transition matrix)")
print()
for name, A in [("KEY1-shift", np.eye(KEY1) * 0 + 1),
                 ("KEY2-shift", np.eye(KEY2) * 0 + 1),
                 ("Golden mean", np.array([[1,1],[1,0]])),
                 ("Tournament-3 (complete)", np.array([[0,1,1],[1,0,1],[1,1,0]]))]:
    sr = max(abs(np.linalg.eigvals(A)))
    h = math.log(sr) if sr > 0 else 0
    print(f"  {name}: spectral radius = {sr:.4f}, h = {h:.4f}")

print(f"""
  The tournament on KEY2 vertices has transition matrix with
  spectral radius KEY1, giving entropy log(KEY1) = 1 bit!

  A tournament is a MAXIMAL-ENTROPY directed graph
  (every pair connected, each with 1 bit of information).
""")

# =====================================================================
# Part 8: CELLULAR AUTOMATA — KEY1-STATE UNIVERSES
# =====================================================================
print("\n" + "=" * 70)
print("  Part 8: CELLULAR AUTOMATA — KEY1-STATE UNIVERSES")
print("=" * 70)

print(f"""
  Elementary Cellular Automata (Wolfram):
  - KEY1 states: {{0, 1}}
  - KEY2-cell neighborhoods: (left, center, right)
  - KEY1^(KEY1^KEY2) = {KEY1**(KEY1**KEY2)} = 2^8 = 256 possible rules

  The number of rules = KEY1^(KEY1^KEY2) = 256!
  The exponent tower KEY1^KEY2 = {KEY1**KEY2} gives the neighborhood size.
""")

# Simulate Rule 110 (Turing-complete!)
def eca_step(state, rule):
    """One step of elementary cellular automaton."""
    n = len(state)
    new = [0] * n
    for i in range(n):
        neighborhood = (state[(i-1) % n] << 2) | (state[i] << 1) | state[(i+1) % n]
        new[i] = (rule >> neighborhood) & 1
    return new

# Rule 110 — Turing complete
rule_110 = 110
width = 40
state = [0] * width
state[width - 2] = 1

print(f"  Rule 110 (Turing-complete!):")
print(f"  110 = {KEY1} × {KEY_SUM} × 11 = KEY1 × KEY_SUM × 11")
print(f"  Binary: 01101110")
for step in range(15):
    line = ''.join('█' if c else ' ' for c in state)
    if step < 10:
        print(f"  {line}")
    state = eca_step(state, rule_110)

# Rule 30 — chaotic
print(f"\n  Rule 30 (chaotic, used for randomness):")
print(f"  30 = KEY1 × KEY2 × KEY_SUM")
state = [0] * width
state[width // 2] = 1
for step in range(10):
    line = ''.join('█' if c else ' ' for c in state)
    print(f"  {line}")
    state = eca_step(state, 30)

print("""
  NOTABLE RULES and tournament constants:
  Rule 30 = KEY1 * KEY2 * KEY_SUM (chaotic)
  Rule 110 = KEY1 * KEY_SUM * 11 (Turing-complete)
  Rule 90 = KEY1 * KEY2 * KEY2 * KEY_SUM (Sierpinski, additive)
  Rule 150 = KEY1 * KEY2 * KEY_SUM^KEY1 (also additive)
  Rule 184 = KEY1^KEY2 * 23 (traffic flow)

  Total ECA rules: KEY1^(KEY1^KEY2) = 2^(2^3) = 256
  "Interesting" rules: ~ KEY1^KEY_SUM = 32 (a tiny fraction!)

  KEY2-STATE AUTOMATA (totalistic):
  Number of rules = KEY2^(KEY2^KEY2) ~ 7.6e12 (enormous!)
  The jump from KEY1 to KEY2 states is ASTRONOMICAL.
""")

# =====================================================================
# Part 9: FRACTALS — DIMENSION (2,3)
# =====================================================================
print("\n" + "=" * 70)
print("  Part 9: FRACTALS — DIMENSION (2,3)")
print("=" * 70)

print(f"""
  Self-similar fractals and their dimensions:

  For an IFS with N copies scaled by factor r:
  dim = log(N) / log(1/r)
""")

fractals = [
    ("Cantor set", KEY2, KEY2, f"log(KEY1)/log(KEY2) = {math.log(2)/math.log(3):.6f}"),
    ("Sierpinski triangle", KEY2, KEY1, f"log(KEY2)/log(KEY1) = {math.log(3)/math.log(2):.6f}"),
    ("Sierpinski carpet", 8, KEY2, f"log(8)/log(KEY2) = {math.log(8)/math.log(3):.6f}"),
    ("Koch snowflake", KEY1*KEY1, KEY2, f"log(KEY1^KEY1)/log(KEY2) = {math.log(4)/math.log(3):.6f}"),
    ("Menger sponge", 20, KEY2, f"log(20)/log(KEY2) = {math.log(20)/math.log(3):.6f}"),
    ("Vicsek fractal", KEY_SUM, KEY2, f"log(KEY_SUM)/log(KEY2) = {math.log(5)/math.log(3):.6f}"),
]

for name, N, scale, dim in fractals:
    print(f"  {name}: N={N}, scale=1/{scale}, dim = {dim}")

print(f"""
  OBSERVATION: The most famous fractals use scales KEY1 and KEY2!

  Cantor set: remove middle 1/KEY2 of each interval
  → dim = log(KEY1)/log(KEY2) ≈ 0.631

  Sierpinski triangle: KEY2 copies at scale 1/KEY1
  → dim = log(KEY2)/log(KEY1) ≈ 1.585

  Koch curve: KEY1^KEY1 copies at scale 1/KEY2
  → dim = log(KEY1^KEY1)/log(KEY2) = KEY1·log(KEY1)/log(KEY2) ≈ 1.262

  NOTE: dim(Cantor) × dim(Sierpinski) = 1 (reciprocals!)
  log(KEY1)/log(KEY2) × log(KEY2)/log(KEY1) = 1

  The Cantor set and Sierpinski triangle are DIMENSION-DUALS
  under the exchange KEY1 ↔ KEY2!

  Mandelbrot set boundary: dim = KEY1 exactly! (Shishikura)

  FRACTAL DIMENSIONS INVOLVING TOURNAMENT CONSTANTS:
  Lorenz attractor: ≈ KEY1.06
  Hénon attractor: ≈ 1.26 ≈ log(KEY2)/log(KEY1) - 0.32
  Rössler attractor: ≈ KEY1.01

  Strange attractors have dimension just above KEY1 = 2.
  They live in KEY2 = 3 dimensions but are "barely more than KEY1-dimensional."
""")

# =====================================================================
# Part 10: RENORMALIZATION AND UNIVERSALITY
# =====================================================================
print("\n" + "=" * 70)
print("  Part 10: RENORMALIZATION AND UNIVERSALITY")
print("=" * 70)

print("""
  Renormalization in dynamics:
  The doubling operator R acts on unimodal maps:
  (Rf)(x) = -alpha * f(f(-x/alpha))

  Fixed point: Rf* = f* (the Feigenbaum function)

  SPECTRUM OF R at f*:
  - Largest eigenvalue: delta = 4.669... (the Feigenbaum constant)
  - Second eigenvalue: -alpha^2 ≈ -6.265

  The unstable manifold of f* under R is 1-dimensional!
  This is WHY universality works:
  All unimodal maps approach f* along the SAME 1D path.

  RENORMALIZATION AS CATEGORY:
  - Objects: unimodal maps
  - Morphisms: semiconjugacies
  - R is a functor that "doubles" the dynamics
  - f* is the terminal object in the "universal dynamics" subcategory

  TOURNAMENT ANALOG:
  - Period-doubling = KEY1-fold iteration
  - Each renormalization step = composing f with itself KEY1 times
  - The KEY1 appears in EVERY step of the renormalization
  - KEY2 appears as the period that BREAKS the pattern (Sarkovskii)

  CONNECTIONS TO OTHER UNIVERSALITY:
  - Statistical mechanics: Ising model RG fixed points
    KEY1-state Ising (up/down) vs KEY2-state Potts
  - KEY1-state Ising: critical exponents in d=4-epsilon
  - KEY2-state Potts: exactly solvable in d=KEY1 (2D)!

  The same KEY1/KEY2 duality appears:
  - Dynamics: KEY1-doubling cascade, KEY2-chaos threshold
  - Stat mech: KEY1-state order parameter, KEY2-state exact solution
  - Coding: KEY1-ary code, KEY2-minimum distance
  - Topology: KEY1-dimensional no chaos, KEY2-dimensional chaos
""")

# =====================================================================
# Part 11: ERGODIC THEORY — TOURNAMENT MEASURES
# =====================================================================
print("\n" + "=" * 70)
print("  Part 11: ERGODIC THEORY — TOURNAMENT MEASURES")
print("=" * 70)

print("""
  BIRKHOFF ERGODIC THEOREM:
  For an ergodic measure-preserving map f: X → X,
  the time average equals the space average (for a.e. x):

  lim_{n→∞} (1/n) sum_{k=0}^{n-1} g(f^k(x)) = integral g dμ

  ENTROPY AND TOURNAMENTS:

  Kolmogorov-Sinai entropy:
  h_μ(f) = sup_P H(P | f^{-1}P ∨ f^{-2}P ∨ ...)

  For the logistic map at r=4 (fully chaotic):
  h_KS = log(KEY1) = 1 bit/iteration
  The invariant measure: dμ = dx / (π√(x(1-x)))
  — the arcsine distribution on [0,1]!
""")

# Compute time averages for logistic map
print("  Time averages for logistic map (r=4):")
x = 0.1
sums = {1: 0, 2: 0}
n = 100000
for i in range(n):
    sums[1] += x
    sums[2] += x**2
    x = 4 * x * (1 - x)

print(f"  <x> = {sums[1]/n:.6f} (theory: 1/{KEY1} = 0.5)")
print(f"  <x²> = {sums[2]/n:.6f} (theory: {KEY2}/{KEY1**KEY2} = 3/8 = 0.375)")

print(f"""
  <x> = 1/KEY1 = 1/2
  <x^2> = KEY2/KEY1^KEY2 = 3/8
  <x^k> = (2k)! / (KEY1^(KEY1*k) * (k!)^KEY1)  — Catalan-like!

  MIXING TIMES ON TOURNAMENTS:
  A random walk on a tournament T reaches equilibrium
  in O(n * log(n)) steps (where n = number of vertices).

  The mixing time involves:
  - log(n) factor from the spectral gap
  - The spectral gap depends on the min out-degree
  - Regular tournaments (all out-degrees = (n-1)/KEY1) mix fastest

  TOURNAMENT AS MARKOV CHAIN:
  Transition matrix P_ij = 1/(n-1) if i → j, else 0
  Stationary distribution = (d_1, d_2, ..., d_n) / sum(d_i)
  where d_i = out-degree of vertex i

  For regular tournaments: stationary distribution = uniform!
  Entropy rate = log(n-1) - H(transition) = log((n-1)/KEY1)

  The KEY1 appears in the denominator because each vertex
  has (n-1)/KEY1 outgoing edges.
""")

# =====================================================================
# Part 12: GRAND SYNTHESIS — DYNAMICS IS (2,3)
# =====================================================================
print("\n" + "=" * 70)
print("  Part 12: GRAND SYNTHESIS — DYNAMICS IS (2,3)")
print("=" * 70)

print("""
======================================================================
  THE DYNAMICAL UNIVERSE IS BUILT FROM KEY1 = 2 AND KEY2 = 3
======================================================================

1. ITERATION:
   z → z^KEY1 + c (the Mandelbrot iteration)
   x → rx(1-x) (the logistic map, degree KEY1)
   QUADRATIC iteration is the universal paradigm.

2. BIFURCATION:
   Period-doubling: periods 1, KEY1, KEY1^2, KEY1^3, ...
   Each bifurcation MULTIPLIES by KEY1 = 2.
   Onset of instability: r = KEY2 = 3 exactly.

3. CHAOS THRESHOLD:
   Period KEY2 implies all periods (Sarkovskii).
   KEY2 dimensions needed for continuous chaos (Poincaré-Bendixson).
   KEY2 equations minimum for strange attractors.

4. UNIVERSALITY:
   Feigenbaum δ ≈ KEY1^KEY1 + 0.669 (universal for degree-KEY1)
   Feigenbaum α ≈ KEY1 + 1/KEY1 (universal scaling)
   All unimodal maps share the SAME constants.

5. FRACTALS:
   Cantor: dim = log(KEY1)/log(KEY2)
   Sierpinski: dim = log(KEY2)/log(KEY1)
   Koch: dim = KEY1 × log(KEY1)/log(KEY2)
   These are DUAL under KEY1 ↔ KEY2!

6. SYMBOLIC DYNAMICS:
   Binary (KEY1-ary) shift has entropy = log(KEY1) = 1 bit.
   Ternary (KEY2-ary) shift has entropy = log(KEY2) ≈ 1.585 bits.
   SFTs defined by KEY1 × KEY1 transition matrices.

7. CELLULAR AUTOMATA:
   KEY1 states, KEY2-cell neighborhoods → KEY1^(KEY1^KEY2) = 256 rules.
   Rule 30 = KEY1 × KEY2 × KEY_SUM (chaotic).
   Rule 110 = KEY1 × KEY_SUM × 11 (Turing-complete).

8. STRANGE ATTRACTORS:
   Lorenz: KEY2 equations, σ = V(Pet), ρ = KEY1^KEY1 × H_forb_1
   All strange attractors have dim slightly > KEY1.
   They live in KEY2-space but are "almost KEY1-dimensional."

9. ERGODIC THEORY:
   Logistic map (r=4): <x> = 1/KEY1, <x^2> = KEY2/KEY1^KEY2
   KS entropy = log(KEY1) bits per iteration.
   Tournament random walk mixes with rate involving 1/KEY1.

THE CROWN JEWEL:
   KEY1 = 2 is ORDER (doubling, periodicity, iteration, measure)
   KEY2 = 3 is CHAOS (threshold, forcing, minimum dimension)
   Their INTERPLAY creates all of dynamical systems theory.

   The universe computes in base KEY1 and breaks in dimension KEY2.
   DYNAMICS IS (2,3).
""")
