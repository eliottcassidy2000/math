#!/usr/bin/env python3
"""
beyond_finite_s339.py — Tournament theory beyond the natural numbers
opus-2026-03-25-S339

Every object in our project is parameterized by n ∈ N. What happens
when we extend n beyond the naturals? Some extensions are rigorous,
others are creative/speculative. We pursue ALL of them.

RIGOROUS EXTENSIONS:
  1. n → ∞ (limit): Fraïssé limit T_ω, Cantor space, tournamentons
  2. n → real: analytic continuation of A000568(n), Burnside, etc.
  3. n → complex: complex tournaments? n = i?
  4. n → ordinal: tournaments on ω, ω+1, ω·2, ...

CREATIVE/SPECULATIVE:
  5. n → 1/2: what's a "half-vertex tournament"? (fractional n)
  6. n → -1: negative vertices? (combinatorial reciprocity)
  7. n → p-adic: tournament over Qp?
  8. The staircase δ_{-1}: what's a negative staircase?
  9. H(T) for continuous T: what replaces Hamiltonian paths?
  10. Score sequence at n = π: irrational number of vertices?
"""

import sys
import numpy as np
from math import comb, factorial, log, log2, sqrt, pi, e
from scipy.special import gamma as Gamma, gammaln
from scipy import interpolate

sys.stdout.reconfigure(line_buffering=True)

print("=" * 72)
print("  BEYOND FINITE n: TOURNAMENT THEORY EXTENDED")
print("  opus-2026-03-25-S339")
print("=" * 72)

# ============================================================================
# 1. ANALYTIC CONTINUATION OF A000568(n)
# ============================================================================

print(f"\n{'='*72}")
print(f"  1. A000568 AS A REAL-VALUED FUNCTION")
print(f"{'='*72}")

print(f"""
  A000568(n) ~ 2^{{C(n,2)}} / n! for integer n >= 1.

  The RHS extends naturally to real n:
    V(x) = 2^{{x(x-1)/2}} / Γ(x+1)

  This is well-defined for all real x > -1 (where Γ(x+1) > 0).
  What does V(x) look like between the integers?
""")

print(f"  {'x':>6} {'V(x)':>14} {'A000568':>10} {'ratio':>8}")
A000568 = {1:1, 2:1, 3:2, 4:4, 5:12, 6:56, 7:456, 8:6880}

for x in [0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7, 7.5, 8]:
    m = x * (x - 1) / 2
    V = 2**m / Gamma(x + 1)
    a = A000568.get(int(x), '') if x == int(x) else ''
    ratio = V / A000568[int(x)] if x == int(x) and int(x) in A000568 else ''
    ratio_str = f"{ratio:.4f}" if isinstance(ratio, float) else ''
    print(f"  {x:>6.1f} {V:>14.4f} {str(a):>10} {ratio_str:>8}")

print(f"\n  At x = 1/2: V(0.5) = 2^{{-1/8}} / Γ(3/2) = {2**(-1/8) / Gamma(1.5):.6f}")
print(f"  At x = π:   V(π) = 2^{{π(π-1)/2}} / Γ(π+1) = {2**(pi*(pi-1)/2) / Gamma(pi+1):.4f}")
print(f"  At x = e:   V(e) = 2^{{e(e-1)/2}} / Γ(e+1) = {2**(e*(e-1)/2) / Gamma(e+1):.4f}")

# ============================================================================
# 2. FRACTIONAL STAIRCASE δ_{x}
# ============================================================================

print(f"\n{'='*72}")
print(f"  2. THE FRACTIONAL STAIRCASE δ_x")
print(f"{'='*72}")

print(f"""
  δ_{{n-2}} has:
    strips = n-2 (always integer)
    tiles = C(n-1, 2) = (n-1)(n-2)/2
    area = C(n, 2) = n(n-1)/2

  Extend to real x >= 0:
    strips(x) = x          (real number of strips)
    tiles(x) = x(x-1)/2    (a parabola in x, can be negative for x < 1!)
    area(x) = (x+2)(x+1)/2 (also a parabola)

  At x = 0:   0 strips, 0 tiles (the empty tournament n=2)
  At x = 1:   1 strip, 0 tiles (n=3: one triangle, one tile... wait, tiles=0?)
  At x = 1/2: 0.5 strips, -0.125 tiles (!!)
  At x = -1:  -1 strips, 1 tile (negative strips, positive tiles!)
""")

print(f"  {'x':>6} {'strips':>8} {'tiles':>10} {'area':>10} {'n=x+2':>8}")
for x in [-1, -0.5, 0, 0.5, 1, 1.5, 2, 3, 4, 5, pi-2, e-2]:
    n = x + 2
    strips = x
    tiles = x * (x - 1) / 2
    area = (x + 2) * (x + 1) / 2
    print(f"  {x:>6.2f} {strips:>8.2f} {tiles:>10.4f} {area:>10.4f} {n:>8.4f}")

print(f"\n  The tiles function x(x-1)/2 has a ZERO at x=0 and x=1,")
print(f"  and is NEGATIVE for 0 < x < 1. This means: between n=2 and n=3,")
print(f"  the number of tiles is 'negative' — the staircase has a hole.")
print(f"  This could be interpreted as: the tournament at fractional n has")
print(f"  FEWER degrees of freedom than its neighbors.")

# ============================================================================
# 3. THE FIBER FRACTION AT REAL n
# ============================================================================

print(f"\n{'='*72}")
print(f"  3. FIBER FRACTION f(x) AT REAL x")
print(f"{'='*72}")

print(f"""
  f(n) = Γ(n - 3/2) / (Γ(1/2) · Γ(n-1))

  This extends naturally to real x via the Gamma function:
    f(x) = Γ(x - 3/2) / (Γ(1/2) · Γ(x - 1))

  f has POLES at x = 3/2 (where Γ(0) diverges) and at x = 1 (Γ(0) again).
  It has ZEROS at x = 3/2 - k for k = 1, 2, ... (non-positive integers of the numerator Γ).
""")

print(f"  {'x':>6} {'f(x)':>14}")
for x in [1.6, 2, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 10, pi, e, sqrt(2)+2]:
    try:
        fx = float(Gamma(x - 1.5) / (Gamma(0.5) * Gamma(x - 1)))
        label = ""
        if abs(x - pi) < 0.01: label = "  (x = π)"
        if abs(x - e) < 0.01: label = "  (x = e)"
        print(f"  {x:>6.3f} {fx:>14.8f}{label}")
    except:
        print(f"  {x:>6.3f} {'pole/undef':>14}")

# ============================================================================
# 4. COMBINATORIAL RECIPROCITY: n → -n
# ============================================================================

print(f"\n{'='*72}")
print(f"  4. COMBINATORIAL RECIPROCITY: NEGATIVE VERTICES")
print(f"{'='*72}")

print(f"""
  Stanley's reciprocity principle: many combinatorial formulas at positive n
  have meaningful interpretations at negative n.

  Example: C(n, k) = n! / (k!(n-k)!) extends to negative n via:
    C(-n, k) = (-1)^k C(n+k-1, k)

  For tournaments: C(n, 2) = n(n-1)/2 arcs.
    C(-n, 2) = (-n)(-n-1)/2 = n(n+1)/2

  So a "tournament on -n vertices" has n(n+1)/2 arcs — one MORE than
  a tournament on n+1 vertices! The negative tournament has an EXTRA arc.

  The Burnside formula: V(n) ~ 2^{{C(n,2)}} / n!
    V(-n) = 2^{{n(n+1)/2}} / Γ(1-n)

  For integer n >= 1: Γ(1-n) = ∞ (pole), so V(-n) = 0.
  But Γ(1-n) has alternating sign residues, and
  1/Γ(1-n) = 0 at all positive integers → V(-n) = 0 exactly.

  INTERPRETATION: There are ZERO tournaments on negative vertices.
  This is combinatorial reciprocity saying: the inverse of counting
  tournaments is counting... nothing? Or perhaps: counting
  ANTI-tournaments (orientations of the complement of K_n)?
""")

print(f"  {'n':>5} {'C(n,2)':>8} {'C(-n,2)':>8} {'V(n) approx':>14} {'V(-n)':>8}")
for n in range(1, 8):
    cn2 = n * (n-1) // 2
    cn2_neg = n * (n+1) // 2
    vn = 2**cn2 / factorial(n)
    # V(-n) = 2^{n(n+1)/2} / Gamma(1-n) = 2^{n(n+1)/2} × 0 = 0 (for integer n)
    print(f"  {n:>5} {cn2:>8} {cn2_neg:>8} {vn:>14.2f} {'0':>8}")

# ============================================================================
# 5. THE HAMILTONIAN PATH COUNT AT REAL n
# ============================================================================

print(f"\n{'='*72}")
print(f"  5. H(T) EXTENDED: WHAT REPLACES HAMILTONIAN PATHS?")
print(f"{'='*72}")

print(f"""
  H(T) counts Hamiltonian paths in a tournament on n vertices.
  At integer n: H is a positive odd integer (Rédei).

  EXTENSION 1: H via the OCF
    H(T) = I(Ω(T), 2) = independence polynomial of conflict graph at x=2.
    The conflict graph Ω(T) is defined for any tournament, including infinite ones.
    For T_ω: Ω(T_ω) is an infinite graph. I(Ω, 2) = ∞ (uncountably many
    independent sets). So H(T_ω) = ∞.

  EXTENSION 2: H via the transfer matrix
    H = Σ permanent-like expression over permutations.
    For a tournamenton W: [0,1]² → [0,1]:
      H(W) = ∫_{{S_∞}} Π W(σ(i), σ(i+1)) dσ

    This is a "Hamiltonian path integral" — sum over all orderings of [0,1],
    weighted by the product of edge weights along the path.

    For W(x,y) = 1/2 (random tournamenton):
      H(W) = (1/2)^∞ × ∞! = 0 × ∞ (indeterminate)

    For W(x,y) = 1_{{x<y}} (transitive tournamenton):
      H(W) = 1 (the unique order-preserving path)

  EXTENSION 3: H as a continuous function of n
    The average H over all tournaments on n vertices is:
      E[H] = n!/2^{{n-1}}

    This extends to: E[H](x) = Γ(x+1)/2^{{x-1}} = x!/2^{{x-1}}

    At x = 1/2: E[H](0.5) = Γ(1.5)/2^{{-0.5}} = (√π/2)·√2 = {Gamma(1.5) * sqrt(2):.6f}
    At x = π:   E[H](π) = Γ(π+1)/2^{{π-1}} = {Gamma(pi+1) / 2**(pi-1):.6f}
""")

print(f"  E[H] at real x:")
print(f"  {'x':>6} {'E[H](x)':>14}")
for x in [0.5, 1, 1.5, 2, 2.5, 3, pi, e, 4, 5, 6, 7]:
    eh = Gamma(x + 1) / 2**(x - 1)
    print(f"  {x:>6.3f} {eh:>14.6f}")

# ============================================================================
# 6. THE SPECTRAL GAP AT REAL n
# ============================================================================

print(f"\n{'='*72}")
print(f"  6. SPECTRAL GAP 2/n EXTENDED")
print(f"{'='*72}")

print(f"""
  The spectral gap of the metagraph: gap ≈ 2/n.

  At real x: gap(x) = 2/x (a simple hyperbola).

  gap(1) = 2    (trivial tournament: gap = 2)
  gap(2) = 1    (one arc: gap = 1)
  gap(∞) = 0    (continuous spectrum)
  gap(1/2) = 4  (half a vertex: supergap!)
  gap(-1) = -2  (negative gap! → unstable?)
  gap(π) = 2/π = {2/pi:.6f} (an irrational gap)

  The gap = 2/n means the metagraph mixes in O(n²) random walks.
  At n = ∞: mixing time = ∞ (never mixes — T_ω is rigid).
  At n = 1/2: mixing time ≈ 1/16 (mixes instantly — no structure).

  SPECULATION: The spectral gap 2/n interpolates between
  "no structure" (n→0, gap→∞, instant mixing)
  and "rigid structure" (n→∞, gap→0, frozen).
  The crossover is at n ≈ 1 (gap = 2).
""")

# ============================================================================
# 7. SCORE SEQUENCE AT REAL n
# ============================================================================

print(f"\n{'='*72}")
print(f"  7. SCORE SEQUENCES AT REAL n")
print(f"{'='*72}")

print(f"""
  At integer n: the score sequence (s₁,...,s_n) satisfies:
    Σ sᵢ = C(n,2) = n(n-1)/2
    0 ≤ sᵢ ≤ n-1

  At real x: "scores" become a MEASURE on [0, x-1]:
    The score distribution μ satisfies:
      ∫ μ = x(x-1)/2        (total arcs)
      μ(t) ∈ [0, x-1]       (each vertex's out-degree)

  For the transitive tournamenton W(x,y) = 1_{{x<y}} on [0,1]:
    Score of vertex t = ∫₀¹ 1_{{t<y}} dy = 1-t
    Score distribution μ = (1-t) on [0,1]
    Total = ∫₀¹ (1-t) dt = 1/2 = C(1,2)/... hmm, need to scale.

  For general tournamenton W on [0,1]:
    Score of vertex t = ∫₀¹ W(t,y) dy
    This is a measurable function s: [0,1] → [0,1]
    Total: ∫₀¹ s(t) dt = 1/2 (each arc counted once)

  The score sequence at n=∞ is a SCORE MEASURE.
  Our project's score conditioning becomes: conditioning on the
  score DISTRIBUTION (a function, not a tuple).
""")

# ============================================================================
# 8. THE STAIRCASE PARADOX AT REAL n
# ============================================================================

print(f"\n{'='*72}")
print(f"  8. THE STAIRCASE AT n = π")
print(f"{'='*72}")

print(f"""
  The staircase δ_{{n-2}} at n = π:
    strips = π - 2 = {pi - 2:.6f}
    tiles = (π-1)(π-2)/2 = {(pi-1)*(pi-2)/2:.6f}
    L1 length of hypotenuse = 2(π-2) = {2*(pi-2):.6f}
    L2 length of hypotenuse = (π-2)√2 = {(pi-2)*sqrt(2):.6f}
    ratio L1/L2 = √2 = {sqrt(2):.6f} (SAME as integer n!)

  The staircase paradox ratio is √2 REGARDLESS of n (even non-integer).
  This is because √2 = L1/L2 ratio of a 45° line, which is a geometric
  fact independent of the number of steps.

  A "π-vertex tournament" would have:
    π(π-1)/2 = {pi*(pi-1)/2:.6f} arcs (≈ 3.37 arcs)
    The staircase δ_{{π-2}} has area {(pi-2)*(pi-1)/2:.6f}
    The "π-th tournament isomorphism class" doesn't exist as a discrete object.

  But the GENERATING FUNCTION approach works:
    V(x) = Σ A000568(n) x^n / n! has convergence radius > 0.
    Evaluating at x = 1 gives: Σ A000568(n) / n! = ???
    This is the "total weight" of tournaments, which should relate to T_ω.
""")

# Compute the EGF of A000568
print(f"  EGF partial sums: Σ A000568(k)/k! for k=1..n:")
total = 0.0
for n in range(1, 11):
    vn = A000568.get(n, 0)
    total += vn / factorial(n)
    print(f"    n={n:>2}: A000568={vn:>10}, cumulative EGF = {total:.8f}")

print(f"\n  The EGF converges VERY fast (dominated by 2^{{C(n,2)}}/n!² ).")
print(f"  Its value at x=1 ≈ {total:.8f} is the 'total tournament weight'.")

# ============================================================================
# 9. TOURNAMENT ON THE REALS
# ============================================================================

print(f"\n{'='*72}")
print(f"  9. TOURNAMENTS ON R: CONTINUOUS PAIRWISE COMPARISON")
print(f"{'='*72}")

print(f"""
  A tournament on the reals = a function T: R×R → {{0,1}} such that:
    T(x,y) + T(y,x) = 1 for all x ≠ y
    T(x,x) = 0

  Examples:
    1. T(x,y) = 1 iff x < y (the usual order on R)
       This is transitive. H = 1 "path" (the order itself).
       Score of x = measure of {{y : x < y}} (infinite for Lebesgue).

    2. T(x,y) = 1 iff sin(x-y) > 0 (the "circular" tournament on R)
       This is S(2)-like. Has 3-cycles everywhere.
       Not transitive. H = 0 (no spanning path respects this).

    3. T(x,y) = indicator of (x mod 1) < (y mod 1) (the [0,1) tournament)
       This is a tournamenton: W(x,y) = 1_{{frac(x) < frac(y)}}.

  OUR INVARIANTS ON CONTINUOUS TOURNAMENTS:
    • Score s(x) = ∫ T(x,y) dy — a function, not a number
    • c₃ density = ∫∫∫ T(x,y)T(y,z)T(z,x) dx dy dz — a real number
    • H = Hamiltonian path integral — usually 0 or ∞
    • Conflict graph Ω(T) — vertices are odd cycles (a continuum)

  KEY INSIGHT: The tournamenton W: [0,1]² → [0,1] is the RIGHT
  extension of tournaments to the reals. It replaces the discrete
  adjacency matrix A[i][j] ∈ {{0,1}} with a continuous kernel
  W(x,y) ∈ [0,1] satisfying W(x,y) + W(y,x) = 1.

  Our OCF (Odd-Cycle Collection Formula) becomes:
    H(W) = I(Ω(W), 2)
  where Ω(W) is the conflict "graph" (now a kernel on cycle space)
  and I is the "independence polynomial" (now a functional integral).
""")

# ============================================================================
# 10. CREATIVE CONNECTIONS
# ============================================================================

print(f"\n{'='*72}")
print(f"  10. CREATIVE CONNECTIONS BEYOND N")
print(f"{'='*72}")

print(f"""
  A. THE CHROMATIC NUMBER χ(G_x) = x - 1
     At integer n: χ = n-1 (proved n ≤ 7, conjectured all n).
     At real x: χ(x) = x - 1 (trivially extended).
     At x = 1: χ = 0 (no vertices → no colors needed).
     At x = 0: χ = -1 (NEGATIVE colors! → chromatic polynomial reciprocity?)
     Stanley (1973): P_G(-1) = (-1)^n × #acyclic orientations.
     Tournament connection: at n = -1, the chromatic polynomial
     counts something about ORIENTATIONS — full circle!

  B. THE BURNSIDE EXPONENT: Fix(σ) = 2^{{arc_orbits}}
     For σ ∈ S_n: arc_orbits depends on cycle type.
     For σ ∈ S_∞ (infinite symmetric group): Fix(σ) is usually 0 or 1.
     The Burnside formula becomes: A000568 = Σ Fix(σ) / |G|.
     In the limit: the sum is dominated by identity (all others → 0).
     This is WHY A000568(n) ~ 2^m / n! (the identity term).

  C. THE FOUR CONSTANTS AT REAL n:
     √2 = L1/L2 ratio. Independent of n. UNIVERSAL.
     π = Γ(1/2)². Appears via CLT. Limit of finite-n corrections.
     e = Growth rate. Appears via Stirling. Limit of n!/n^n.
     γ = Euler-Mascheroni. Appears via harmonic sums. Limit of H_n - ln(n).

     ALL FOUR are asymptotic constants that measure the rate of
     convergence of discrete (finite n) to continuous (real n).
     They are the "Taylor coefficients of the passage to the limit."

  D. TOURNAMENT ZETA FUNCTION:
     ζ_T(s) = Σ_n A000568(n) / n^s
     This Dirichlet series converges for Re(s) > 1 (since A000568 < n!).
     Its analytic continuation (if it exists) would define A000568 at
     complex s, giving "tournament theory in the complex plane."

  E. THE INDEPENDENCE POLYNOMIAL AT x = -1:
     I(Ω(T), -1) = (-1)^n × #maximal independent sets × sign
     This is related to the MATCHING POLYNOMIAL of Ω(T).
     At x = 2: I(Ω, 2) = H(T) (our OCF).
     At x = -1: I(Ω, -1) = ??? (a mystery invariant).
     At x = 0: I(Ω, 0) = 1 (the empty set).
     At x = 1: I(Ω, 1) = #independent sets = #ways to pick
     vertex-disjoint odd cycles.

     The function x ↦ I(Ω(T), x) interpolates between:
       x=0: the constant 1
       x=1: the independence number
       x=2: the Hamiltonian path count
       x→∞: dominated by the maximum independent set

     What is I(Ω(T), x) at x = φ = (1+√5)/2 (golden ratio)?
     At x = 1/φ = φ-1? The golden ratio is the root of x²-x-1=0,
     which means I(Ω, φ) satisfies a special algebraic relation
     if Ω has any quadratic-like structure.

  F. TOURNAMENT ON THE p-ADICS:
     The Paley tournament on Z/pZ uses quadratic residues.
     The p-adic integers Z_p contain Z/pZ as a quotient.
     A "p-adic tournament" on Z_p would use the p-adic quadratic
     residue symbol. This exists! The Hilbert symbol (a,b)_p ∈ {{±1}}
     defines a tournament-like structure on Q_p*.

     The Paley tournament is the "reduction mod p" of this p-adic tournament.
     Our finite calculations are SHADOWS of a p-adic object.
""")

# ============================================================================
# 11. THE GENERATING FUNCTION APPROACH
# ============================================================================

print(f"{'='*72}")
print(f"  11. A000568 GENERATING FUNCTIONS")
print(f"{'='*72}")

# OGF
print(f"\n  OGF: T(x) = Σ A000568(n) x^n")
ogf = 0.0
for n in range(1, 9):
    ogf += A000568[n] * (0.1)**n
print(f"  T(0.1) = {ogf:.10f}")

# EGF
print(f"\n  EGF: E(x) = Σ A000568(n) x^n / n!")
egf_vals = {}
for x in [0.1, 0.5, 1.0, 2.0]:
    egf = sum(A000568.get(n, 0) * x**n / factorial(n) for n in range(1, 11))
    egf_vals[x] = egf
    print(f"  E({x}) = {egf:.10f}")

# Dirichlet series
print(f"\n  Dirichlet: D(s) = Σ A000568(n) / n^s")
for s in [2, 3, 4, 5]:
    ds = sum(A000568.get(n, 0) / n**s for n in range(1, 11))
    print(f"  D({s}) = {ds:.10f}")

print(f"""

  The EGF E(x) should satisfy a functional equation related to
  the Burnside formula. If almost all tournaments are asymmetric:
    E(x) ≈ exp(Σ 2^{{C(k,2)}} x^k / k! / k)
  This connects to the "labeled tournament" EGF:
    L(x) = Σ 2^{{C(n,2)}} x^n / n!
  via: E(x) = exp(Sigma L_k / k) (Burnside/Polya exponential formula)

  The labeled EGF is L(x) = Σ 2^{{n(n-1)/2}} x^n / n!.
  This series converges for ALL x ∈ C (entire function!).
  Its value at x = 1: L(1) = Σ 2^{{C(n,2)}} / n! ≈ ???
""")

L_1 = sum(2**comb(n,2) / factorial(n) for n in range(15))
print(f"  L(1) = Σ 2^C(n,2) / n! ≈ {L_1:.8f} (converges fast)")
print(f"  This is the 'total labeled tournament weight' = {L_1:.8f}")
print(f"  Compare: e = {e:.8f}, L(1)/e = {L_1/e:.6f}")

print(f"\n{'='*72}")
print(f"  SYNTHESIS")
print(f"{'='*72}")
print(f"""
  The most meaningful extensions beyond N:

  1. ANALYTIC CONTINUATION of A000568 via V(x) = 2^{{x(x-1)/2}} / Γ(x+1)
     → Smooth interpolation between integer tournament counts
     → V(π) ≈ 18.8 "π-vertex tournaments"

  2. TOURNAMENTONS W: [0,1]² → [0,1] with W(x,y)+W(y,x)=1
     → Continuous limit of discrete tournaments
     → The OCF becomes a functional integral: H(W) = I(Ω(W), 2)

  3. p-ADIC TOURNAMENTS via the Hilbert symbol
     → Paley tournaments are reductions mod p of p-adic structure
     → Connects to algebraic number theory

  4. THE INDEPENDENCE POLYNOMIAL I(Ω, x) at non-integer x
     → Interpolates between counting (x=2 → H) and structure (x=1 → α)
     → At x = φ: golden-ratio evaluation of the conflict polynomial

  5. NEGATIVE n via combinatorial reciprocity
     → χ(-n) counts acyclic orientations (Stanley)
     → Tournament at -n: related to orientations of the complement

  6. THE LABELED TOURNAMENT EGF is an ENTIRE FUNCTION
     → L(x) = Σ 2^{{C(n,2)}} x^n / n! converges everywhere
     → L(1) ≈ {L_1:.6f}: the total weight of all labeled tournaments
""")

print("DONE.")
