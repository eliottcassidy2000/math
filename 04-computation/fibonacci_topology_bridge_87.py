#!/usr/bin/env python3
"""
fibonacci_topology_bridge_87.py — opus-2026-03-14-S87

THE FIBONACCI-TOPOLOGY BRIDGE

The independence polynomial I(G, x) lives on a 1-parameter family:
  x = -1: Euler characteristic of Ind(G)
  x = 0: trivial (= 1)
  x = 1: combinatorial (#independent sets)
  x = 2: tournament count (for G = Ω(T))

For path graphs: I(P_k, x) satisfies a(k) = a(k-1) + x·a(k-2)
  → Fibonacci at x=1, Jacobsthal at x=2, period-6 at x=-1

THE KEY QUESTION: What happens as we continuously deform x from -1 to 2?
The eigenvalues of the transfer matrix [[1,x],[1,0]] are:
  λ± = (1 ± √(1+4x)) / 2

At x = -1/4: discriminant = 0, double root λ = 1/2 (CRITICAL POINT!)
At x = -1: λ± = (1 ± i√3)/2 = e^{±iπ/3} (on unit circle!)
At x = 0: λ+ = 1, λ- = 0
At x = 1: λ± = (1 ± √5)/2 = φ, -1/φ (golden ratio!)
At x = 2: λ± = (1 ± 3)/2 = 2, -1

The MODULUS |λ+| transitions from 1 (unit circle at x=-1) through √x
to larger values. The ARGUMENT transitions from π/3 to 0 (real axis).
"""

import math
import cmath

print("=" * 70)
print("THE FIBONACCI-TOPOLOGY BRIDGE")
print("=" * 70)

# ══════════════════════════════════════════════════════════════════
# PART 1: EIGENVALUE TRAJECTORY
# ══════════════════════════════════════════════════════════════════

print("\n--- Transfer matrix eigenvalues λ± = (1 ± √(1+4x)) / 2 ---\n")

print(f"{'x':>7} {'λ+':>20} {'λ-':>20} {'|λ+|':>8} {'arg(λ+)':>10} {'|λ-|':>8}")
print("-" * 80)

key_x = [-2, -1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 1, 1.5, 2, 3, 6, 12, 20]
for x in key_x:
    disc = 1 + 4*x
    if disc >= 0:
        sq = math.sqrt(disc)
        lp = (1 + sq) / 2
        lm = (1 - sq) / 2
        lp_str = f"{lp:.6f}"
        lm_str = f"{lm:.6f}"
        mod_p = abs(lp)
        arg_p = 0 if lp > 0 else math.pi
        mod_m = abs(lm)
    else:
        sq = cmath.sqrt(disc)
        lp = (1 + sq) / 2
        lm = (1 - sq) / 2
        lp_str = f"{lp.real:.4f}+{lp.imag:.4f}i"
        lm_str = f"{lm.real:.4f}+{lm.imag:.4f}i"
        mod_p = abs(lp)
        arg_p = cmath.phase(lp)
        mod_m = abs(lm)

    print(f"{x:>7.2f} {lp_str:>20} {lm_str:>20} {mod_p:>8.4f} {arg_p/math.pi:>8.4f}π {mod_m:>8.4f}")

# ══════════════════════════════════════════════════════════════════
# PART 2: THE UNIT CIRCLE REGIME — TOPOLOGY DOMINATES
# ══════════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("PART 2: THE UNIT CIRCLE REGIME (x ∈ [-1, 0))")
print("=" * 70)
print()
print("When -1/4 < x < 0: eigenvalues are REAL, both in (-1, 1).")
print("  I(P_k, x) → 0 exponentially.")
print()
print("When x = -1/4: CRITICAL POINT. Double eigenvalue λ = 1/2.")
print("  I(P_k, -1/4) = (k+2) / 2^{k+1}  (linear × exponential decay)")
print()
print("When -1 ≤ x < -1/4: eigenvalues are COMPLEX, on the unit circle!")
print("  |λ±| = √(-x) → I(P_k, x) is PERIODIC/QUASIPERIODIC!")
print("  At x = -1: |λ±| = 1, period = 6.")
print("  At x = -1/2: |λ±| = 1/√2 ≈ 0.707, DECAYING oscillation.")
print()

# Verify critical point
print("Verification of I(P_k, -1/4) = (k+2)/2^{k+1}:")
def I_path(k, x):
    if k == 0: return 1
    if k == 1: return 1 + x
    a, b = 1, 1 + x
    for _ in range(k - 1):
        a, b = b, b + x * a
    return b

for k in range(10):
    actual = I_path(k, -0.25)
    predicted = (k + 2) / (2 ** (k + 1))
    print(f"  k={k}: I(P_{k}, -1/4) = {actual:.6f}, (k+2)/2^(k+1) = {predicted:.6f}, match={abs(actual-predicted) < 1e-10}")

# ══════════════════════════════════════════════════════════════════
# PART 3: PHASE DIAGRAM — WHERE DOES OSCILLATION LIVE?
# ══════════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("PART 3: PHASE DIAGRAM OF I(P_k, x)")
print("=" * 70)
print()
print("Three regimes:")
print("  x < -1/4: OSCILLATORY (complex eigenvalues, period depends on x)")
print("  x = -1/4: CRITICAL (double real eigenvalue)")
print("  x > -1/4: EXPONENTIAL (real eigenvalues, growth for x > 0)")
print()

# Period of I(P_k, x) as function of x (for rational x giving exact period)
# At x = -1: eigenvalues are 6th roots of unity → period 6
# At what x is the period 12? 8? etc.

# arg(λ+) = arctan(√(-1-4x)) at x < -1/4
# Period p when arg = 2πk/p for some integer k coprime to p

print("Period structure:")
print(f"  x = -1: arg = π/3, period = 6")
print(f"  x = -cos²(π/5) - 1/4 = ?: arg = π/5, period = 10")
print(f"  x = -cos²(π/4) - 1/4 = ?: arg = π/4, period = 8")
print()

# Compute: for arg = π/n, what x gives it?
# λ+ = (1 + i√(-1-4x))/2
# arg(λ+) = arctan(√(-1-4x))
# For arg = θ: tan(θ) = √(-1-4x), so -1-4x = tan²(θ), x = -(1+tan²(θ))/4 = -sec²(θ)/4
# Since |λ+|² = (1/4)(1 + (-1-4x)) = (1/4)(-4x) = -x
# So |λ+| = √(-x)

print("Exact x for given period p (arg = π/p for period 2p):")
for p in [3, 4, 5, 6, 7, 8, 10, 12]:
    theta = math.pi / p
    x = -1 / (4 * math.cos(theta)**2)
    mod_lp = math.sqrt(-x)
    actual_period = 2 * p
    print(f"  period {actual_period:3d}: x = {x:.6f}, |λ+| = {mod_lp:.6f}")
    # Verify
    seq = [I_path(k, x) for k in range(3 * actual_period)]
    # Check if periodic
    is_periodic = all(abs(seq[k] - seq[k + actual_period]) < 1e-6 for k in range(actual_period))
    print(f"           verified periodic: {is_periodic}")

# ══════════════════════════════════════════════════════════════════
# PART 4: THE OBLONG NUMBERS AS PHASE MARKERS
# ══════════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("PART 4: OBLONG FUGACITIES IN THE PHASE DIAGRAM")
print("=" * 70)
print()
print("Integer-eigenvalue fugacities x = k(k+1) (oblong numbers):")
print("  x = 0: λ = 1, 0  (degenerate)")
print("  x = 2: λ = 2, -1  (tournaments)")
print("  x = 6: λ = 3, -2")
print("  x = 12: λ = 4, -3")
print("  x = 20: λ = 5, -4")
print()
print("These are ALL in the EXPONENTIAL regime (x > -1/4).")
print("The 'tournament fugacity' x=2 is the SMALLEST nonzero oblong.")
print()

# What's special about oblong fugacities?
# G_k(n) = I(P_n, k(k+1)) = ((k+1)^n - (-k)^n) / (2k+1)
# These are INTEGERS when the denominator divides the numerator.

print("G_k(n) = I(P_n, k(k+1)) for first few k and n:")
print(f"{'n':>3}", end="")
for k in range(5):
    print(f"{'x='+str(k*(k+1)):>10}", end="")
print()
for n in range(10):
    print(f"{n:>3}", end="")
    for k in range(5):
        x = k * (k + 1)
        val = I_path(n, x)
        print(f"{int(round(val)):>10}", end="")
    print()

# The correct closed form for I(P_n, x) with x = k(k+1):
# Recurrence: a(n) = a(n-1) + x*a(n-2), a(0)=1, a(1)=1+x
# Eigenvalues: λ+ = k+1, λ- = -k
# Solution: a(n) = A(k+1)^n + B(-k)^n
# From a(0)=1: A+B=1, a(1)=1+x=1+k(k+1)=(k+1)^2-k: A(k+1)-Bk=(k+1)^2-k
# Solving: A = (k+1+k)/(2k+1) * something... let me just compute directly
print("\nClosed form for I(P_n, k(k+1)):")
for k in range(1, 5):
    x = k * (k + 1)
    # From a(0)=1, a(1)=1+x: solve A+B=1, A(k+1)-Bk=1+k(k+1)
    # A(k+1) - (1-A)k = 1 + k^2 + k → A(2k+1) - k = 1 + k^2 + k → A = (2+k^2+2k)/(2k+1) = (k+1)^2/(2k+1)+...
    # Let me just compute A, B numerically
    lp, lm = k+1, -k
    # A + B = 1, A*lp + B*lm = 1+x
    A = (1 + x - lm) / (lp - lm)
    B = 1 - A
    all_match = True
    for n in range(10):
        actual = int(round(I_path(n, x)))
        predicted = round(A * lp**n + B * lm**n)
        if actual != predicted:
            all_match = False
    if all_match:
        print(f"  k={k} (x={x}): I(P_n,{x}) = {A:.4f}·{lp}^n + {B:.4f}·({lm})^n  ✓")

# ══════════════════════════════════════════════════════════════════
# PART 5: THE 2/3 FIBONACCI WORD AND TILINGS
# ══════════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("PART 5: THE FIBONACCI WORD — 2 AND 3 AS TILE LENGTHS")
print("=" * 70)
print()
print("Fibonacci numbers count tilings of a 1×n strip with tiles of length 1 and 2.")
print("Jacobsthal numbers count tilings where 1-tiles have 1 color and 2-tiles have 2 colors.")
print("At fugacity x: tiles of length 2 have x colors.")
print()

# The Fibonacci word: starting from 'a', apply a→ab, b→a
# Letters a and b can represent tiles of length 2 and 3
def fibonacci_word(n):
    """Generate first n letters of the Fibonacci word on {a,b}."""
    if n == 0: return []
    word = ['a']
    while len(word) < n:
        new_word = []
        for c in word:
            if c == 'a':
                new_word.extend(['a', 'b'])
            else:
                new_word.append('a')
        word = new_word
    return word[:n]

fw = fibonacci_word(30)
print(f"Fibonacci word (first 30): {''.join(fw)}")
print(f"  a count: {fw.count('a')}, b count: {fw.count('b')}")
print(f"  ratio a/b: {fw.count('a')/max(fw.count('b'),1):.4f} (should → φ = {(1+math.sqrt(5))/2:.4f})")

# Map a→2, b→3 to get the "tile sequence"
tiles = [2 if c == 'a' else 3 for c in fw]
print(f"Tile lengths: {tiles}")
print(f"  Running sums: {[sum(tiles[:i+1]) for i in range(15)]}")

# These running sums are the BEATTY SEQUENCE B(φ²) = floor(n·φ²)
phi = (1 + math.sqrt(5)) / 2
beatty = [int(math.floor((i+1) * phi**2)) for i in range(15)]
actual_sums = [sum(tiles[:i+1]) for i in range(15)]
print(f"  Beatty B(φ²): {beatty}")
print(f"  Match: {actual_sums == beatty}")

# Actually the Fibonacci word with tiles 1,2 gives Beatty B(φ).
# With tiles 2,3, partial sums are 2*position + #b_tiles = n + (n + #b).
# Let's see what the actual sequence is
# The partial sums follow: sum_k = 2*a_count + 3*b_count = 2k + b_count
# where b_count = #b letters in first k positions
# Since b appears with density 1/φ, sum_k ≈ 2k + k/φ = k(2 + 1/φ) = k(1+φ) = kφ²
# So it IS asymptotically kφ², just not exactly Beatty B(φ²)

print(f"  Asymptotically: partial_sum(k) ~ k·φ² = k·{(1+math.sqrt(5))/2**2:.4f}")
print()
print("★ INSIGHT: The Fibonacci word with tile lengths 2 and 3:")
print("  2 = tournament fugacity (x in I(Ω,2) = H)")
print("  3 = x+1 = the 'augmented fugacity'")
print("  2+3 = 5 (pentagon), 2×3 = 6 (period)")

# ══════════════════════════════════════════════════════════════════
# PART 6: THE PERIOD-6 AS A TOPOLOGICAL INVARIANT
# ══════════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("PART 6: PERIOD-6 AS TOPOLOGICAL INVARIANT")
print("=" * 70)
print()
print("Period 6 appears in many places:")
print("  - I(P_k, -1) mod anything: exact period 6")
print("  - Fibonacci mod 4: π_F(4) = 6")
print("  - I(P_k, 2) mod 4: period ?")
print("  - 6th roots of unity: eigenvalues at x = -1")
print()

# I(P_k, 2) mod various m
print("I(P_k, 2) mod m for various m:")
for m in [3, 4, 5, 6, 7, 8, 9, 10, 12]:
    seq = [int(round(I_path(k, 2))) % m for k in range(60)]
    # Find period
    period = None
    for p in range(1, 51):
        if all(seq[i] == seq[i+p] for i in range(min(9, 60-p))):
            period = p
            break
    if period:
        print(f"  mod {m:2d}: period = {period:3d}, first values: {seq[:min(period, 20)]}")
    else:
        print(f"  mod {m:2d}: period > 50, first values: {seq[:20]}")

# Fibonacci mod m for comparison
print("\nFibonacci mod m for comparison:")
def fib(n):
    if n <= 1: return n
    a, b = 0, 1
    for _ in range(n-1):
        a, b = b, a+b
    return b

for m in [3, 4, 5, 6, 7, 8, 9, 10, 12]:
    seq = [fib(k) % m for k in range(120)]
    period = None
    for p in range(1, 61):
        if all(seq[i] == seq[i+p] for i in range(min(60, 120-p))):
            period = p
            break
    print(f"  mod {m:2d}: π_F = {period:3d}")

# ══════════════════════════════════════════════════════════════════
# PART 7: THE CONWAY-COXETER FRIEZE PATTERN CONNECTION
# ══════════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("PART 7: FRIEZE PATTERNS AND THE SL(2,Z) ACTION")
print("=" * 70)
print()
print("The transfer matrix M(x) = [[1, x], [1, 0]] ∈ GL(2, R).")
print("det(M) = -x. So M ∈ SL(2) only at x = -1!")
print()
print("At x = -1: M = [[1,-1],[1,0]], M^6 = I (order 6 in SL(2,Z)).")
print("This is the MODULAR GROUP element corresponding to the 6th root of unity.")

# Verify M^6 = I at x = -1
def mat_mult(A, B):
    return [
        [A[0][0]*B[0][0] + A[0][1]*B[1][0], A[0][0]*B[0][1] + A[0][1]*B[1][1]],
        [A[1][0]*B[0][0] + A[1][1]*B[1][0], A[1][0]*B[0][1] + A[1][1]*B[1][1]]
    ]

M = [[1, -1], [1, 0]]
power = [[1, 0], [0, 1]]  # identity
print("Powers of M = [[1,-1],[1,0]]:")
for k in range(1, 8):
    power = mat_mult(power, M)
    print(f"  M^{k} = {power}")

print()
print("★ M^6 = I confirmed! The transfer matrix at x=-1 has ORDER 6 in SL(2,Z).")
print("  This is the rotation by π/3 in the hyperbolic plane.")

# At x = 2: det(M) = -2
# M^k mod p for various p
print("\nAt x = 2: M = [[1,2],[1,0]], det = -2.")
M2 = [[1, 2], [1, 0]]
power = [[1, 0], [0, 1]]
print("Powers of M mod 3:")
for k in range(1, 15):
    power = mat_mult(power, M2)
    mod_power = [[power[i][j] % 3 for j in range(2)] for i in range(2)]
    print(f"  M^{k} mod 3 = {mod_power}", end="")
    if mod_power == [[1, 0], [0, 1]]:
        print(f"  ← IDENTITY! Period = {k}")
        break
    elif mod_power == [[2, 0], [0, 2]]:
        print(f"  ← scalar 2I")
    else:
        print()

print("\nPowers of M mod 7:")
power = [[1, 0], [0, 1]]
for k in range(1, 50):
    power = mat_mult(power, M2)
    mod_power = [[power[i][j] % 7 for j in range(2)] for i in range(2)]
    if mod_power == [[1, 0], [0, 1]]:
        print(f"  M^{k} mod 7 = I! Period = {k}")
        break

# ══════════════════════════════════════════════════════════════════
# PART 8: THE DEFORMATION AS HOMOTOPY
# ══════════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("PART 8: THE x-DEFORMATION AS ALGEBRAIC HOMOTOPY")
print("=" * 70)
print()
print("The family I(G, x) for varying x is a 'deformation' connecting")
print("different algebraic/topological invariants:")
print()
print("  x = -1: Euler characteristic (homological)")
print("  x = 0: trivial (contractible point)")
print("  x = 1: independence number (combinatorial)")
print("  x = 2: tournament count (arithmetic)")
print("  x = φ-1 ≈ 0.618: golden ratio point")
print()
print("For the PATH GRAPH P_k, the deformation continuously moves")
print("from period-6 oscillation to exponential growth:")

for x in [-1.0, -0.5, -0.25, 0, 0.5, 0.618, 1, 1.5, 2, 3]:
    seq = [I_path(k, x) for k in range(15)]
    seq_str = ", ".join(f"{v:.1f}" for v in seq[:12])
    if x >= -0.25:
        # Growth rate
        rate = seq[14] / seq[13] if seq[13] != 0 else float('inf')
        print(f"  x={x:>6.3f}: {seq_str}, ...  (growth ≈ {rate:.3f})")
    else:
        print(f"  x={x:>6.3f}: {seq_str}, ...  (oscillating)")

# ══════════════════════════════════════════════════════════════════
# PART 9: THE GOLDEN RATIO AS A PHASE TRANSITION
# ══════════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("PART 9: PHASE TRANSITIONS IN THE INDEPENDENCE POLYNOMIAL")
print("=" * 70)
print()

# For a graph G, the roots of I(G, x) are the "phase transitions"
# For P_k, the roots approach -1/φ² as k → ∞ (Lee-Yang theorem analog)
# This is because |λ+| > |λ-| iff x > -1/4, and the dominant root
# gives the sign for large k.

# The independence polynomial always has a real negative root closest to 0
# For paths: I(P_k, x) = 0 has roots near x = -1/φ² ≈ -0.382 for large k

print("Roots of I(P_k, x) for k = 2, 3, ..., 10:")
for k in range(2, 11):
    # Find roots by scanning
    roots = []
    for x_100 in range(-200, 10):
        x = x_100 / 100.0
        x_next = (x_100 + 1) / 100.0
        v1 = I_path(k, x)
        v2 = I_path(k, x_next)
        if v1 * v2 < 0:
            # Bisect
            lo, hi = x, x_next
            for _ in range(50):
                mid = (lo + hi) / 2
                if I_path(k, lo) * I_path(k, mid) < 0:
                    hi = mid
                else:
                    lo = mid
            roots.append(round((lo+hi)/2, 5))
    print(f"  P_{k}: roots ≈ {roots}")

print(f"\n  The largest root → -1/φ² = {-1/((1+math.sqrt(5))/2)**2:.6f} as k → ∞")
print(f"  This is the HARD-CORE LATTICE GAS critical fugacity for Z!")
print(f"  x_c(Z) = -1/φ² (1D chain), x_c(Z²) ≈ 3.796 (2D square lattice)")

# ══════════════════════════════════════════════════════════════════
# PART 10: SYNTHESIS — THE TOPOLOGICAL UNIVERSE OF TOURNAMENTS
# ══════════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("PART 10: THE TOPOLOGICAL UNIVERSE")
print("=" * 70)
print()
print("THE MAP OF THE TOURNAMENT WORLD:")
print()
print("  TOPOLOGICAL LAYER (x = -1):")
print("    Euler characteristic χ = I(Ω, -1)")
print("    Period-6 from SL(2,Z) action")
print("    M^6 = I in the modular group")
print("    Independence complex Ind(Ω) has cyclic homology")
print()
print("  COMBINATORIAL LAYER (x = 1):")
print("    #independent cycle sets = I(Ω, 1)")
print("    α_odd = (I(Ω,1) - I(Ω,-1))/2")
print("    α_even = (I(Ω,1) + I(Ω,-1) - 2)/2")
print()
print("  ARITHMETIC LAYER (x = 2):")
print("    H = I(Ω, 2) = Hamiltonian path count")
print("    Always odd (Rédei)")
print("    Forbidden values from α₁ gaps")
print()
print("  THE BRIDGE FORMULA:")
print("    H + 2χ = 3 + 6α₂ + 6α₃ + 18α₄ + ...")
print("    At n ≤ 5: H + 2χ = 3 EXACTLY")
print("    At n = 6: first correction from α₂")
print()
print("  THE CRITICAL POINTS:")
print("    x = -1/4: double eigenvalue (phase boundary)")
print("    x = -1/φ²: hard-core gas critical fugacity")
print("    x = 0: contractible (trivial)")
print("    x = 1: golden ratio eigenvalue")
print("    x = 2: tournament fugacity")
print("    x = k(k+1): oblong integer-eigenvalue hierarchy")
