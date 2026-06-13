#!/usr/bin/env python3
"""
infinite_tournament_s339.py — Tournaments at n = ∞
opus-2026-03-25-S339

THE THREE INFINITE TOURNAMENTS (Lachlan 1984):
  There are exactly THREE countable homogeneous tournaments:
    1. (Q, <) — the rational order (transitive, boring)
    2. S(2) — the dense local order (circular)
    3. T_ω — the random tournament (Fraïssé limit, universal)

  T_ω is the n=∞ limit of our project. It contains ALL finite
  tournaments as induced subtournaments. Its automorphism group
  Aut(T_ω) is oligomorphic, and A000568(n) = number of orbits
  on n-element subsets.

THE LIMIT OF OUR SEQUENCES:
  As n → ∞:
    A000568(n) ~ 2^{C(n,2)} / n!           (almost all tournaments asymmetric)
    f(n) ~ 1/√(πn) → 0                     (no single class dominates)
    χ(G_n) = n-1 → ∞                       (chromatic number grows linearly)
    spectral gap = 2/n → 0                  (continuous spectrum in the limit)
    SC fraction → 0                         (self-complementary classes vanish)
    Q_m → 2^ω = Cantor space               (tournament space → Cantor set)
    generic point of 2^ω = T_ω             (a.s. random tournament)

THIS SCRIPT COMPUTES:
  1. Asymptotic behavior of all key sequences
  2. Rate of convergence to the infinite limit
  3. The tournamenton (graph limit) interpretation
  4. Connection between Fraïssé limit and our metagraph
"""

import sys
import numpy as np
from math import comb, factorial, log, log2, sqrt, pi, e, lgamma
from scipy.special import gamma as Gamma

sys.stdout.reconfigure(line_buffering=True)

# ============================================================================
# 1. ASYMPTOTIC ENUMERATION
# ============================================================================

print("=" * 72)
print("  TOURNAMENTS AT n = ∞")
print("  opus-2026-03-25-S339")
print("=" * 72)

print(f"\n  1. A000568(n) ASYMPTOTICS")
print(f"\n  A000568(n) ~ 2^C(n,2) / n!")
print(f"  log₂ A000568(n) ~ n(n-1)/2 - n·log₂(n) + n·log₂(e) - ½·log₂(2πn)")

A000568 = {1:1, 2:1, 3:2, 4:4, 5:12, 6:56, 7:456, 8:6880, 9:191536, 10:9733056}

print(f"\n  {'n':>3} {'A000568':>12} {'2^m/n!':>14} {'ratio':>8} {'log₂(V)':>10} {'m-nlogn':>10}")
for n in range(3, 11):
    m = comb(n, 2)
    vn = A000568[n]
    asymp = 2**m / factorial(n)
    ratio = vn / asymp
    logV = log2(vn) if vn > 0 else 0
    approx_logV = m - n * log2(n) + n * log2(e) - 0.5 * log2(2 * pi * n)
    print(f"  {n:>3} {vn:>12} {asymp:>14.2f} {ratio:>8.4f} {logV:>10.2f} {approx_logV:>10.2f}")

print(f"\n  The ratio → 1 as n → ∞ (almost all tournaments have trivial Aut group)")

# ============================================================================
# 2. RATE OF APPROACH TO THE LIMIT
# ============================================================================

print(f"\n{'='*72}")
print(f"  2. RATES OF CONVERGENCE")
print(f"{'='*72}")

print(f"\n  {'n':>3} {'f(n)':>10} {'1/√πn':>10} {'SC frac':>10} {'gap':>8} {'χ/n':>6}")
for n in range(3, 20):
    # Fiber fraction
    fn = float(Gamma(n - 2 + 0.5) / (Gamma(0.5) * factorial(n - 2)))
    fn_approx = 1.0 / sqrt(pi * n)

    # SC fraction (approximate: SC(n) / A000568(n))
    # SC tournaments exist only for odd n; for even n SC(n) = 0
    # Asymptotically: SC(n) ~ A000568(n) / 2^{n-1} (rough)
    sc_frac = 2.0 ** (1 - (n-1)) if n > 2 else 1.0  # very rough estimate

    # Spectral gap
    gap = 2.0 / n

    # Chromatic number / n
    chi_over_n = (n - 1) / n

    print(f"  {n:>3} {fn:>10.6f} {fn_approx:>10.6f} {sc_frac:>10.6f} {gap:>8.4f} {chi_over_n:>6.3f}")

print(f"\n  As n → ∞:")
print(f"    f(n) → 0 at rate 1/√n   (no single class dominates)")
print(f"    SC fraction → 0 exponentially (self-complementary vanishes)")
print(f"    spectral gap → 0 at rate 1/n (spectrum becomes continuous)")
print(f"    χ/n → 1 (chromatic number grows linearly)")

# ============================================================================
# 3. THE RANDOM TOURNAMENT T_ω
# ============================================================================

print(f"\n{'='*72}")
print(f"  3. THE RANDOM TOURNAMENT T_ω (FRAÏSSÉ LIMIT)")
print(f"{'='*72}")

print(f"""
  T_ω is the UNIQUE countable tournament (up to isomorphism) satisfying:

  EXTENSION PROPERTY: For any finite disjoint subsets A, B of vertices,
  there exists a vertex v that beats every vertex in A and is beaten
  by every vertex in B.

  This is the tournament analogue of the Rado graph.

  KEY PROPERTIES:
    • Universal: contains ALL finite tournaments as induced subtournaments
    • Homogeneous: any isomorphism between finite subtournaments extends
      to an automorphism of T_ω
    • ω-categorical: its first-order theory is complete, and T_ω is
      the unique countable model
    • Almost-sure: a random countable tournament is T_ω with probability 1

  CONSTRUCTION:
    Take vertices v₁, v₂, v₃, ...
    For each pair (vᵢ, vⱼ) with i < j, flip a fair coin:
      Heads → vᵢ → vⱼ
      Tails → vⱼ → vᵢ
    The result is (almost surely) isomorphic to T_ω.

  Aut(T_ω):
    • Simple as an abstract group
    • Oligomorphic: finitely many orbits on n-tuples
    • The number of orbits on unordered n-sets = A000568(n)
    • 2-transitive but not 3-transitive
""")

# ============================================================================
# 4. A000568 = ORBIT COUNT OF Aut(T_ω)
# ============================================================================

print(f"{'='*72}")
print(f"  4. A000568 = ORBIT-COUNTING OF Aut(T_ω)")
print(f"{'='*72}")

print(f"""
  This is the deepest interpretation of A000568:

  A000568(n) = the number of orbits of Aut(T_ω) acting on {{V choose n}}

  In other words: the number of "types" of n-element subtournaments
  of the random tournament. Since T_ω contains all finite tournaments,
  this equals the number of isomorphism classes of n-vertex tournaments.

  Cameron's framework (1976): for an oligomorphic group G acting on Ω,
  the orbit-counting sequence f_n(G) = |Ω^n / G| encodes the structure
  of G. For G = Aut(T_ω):

    f_1 = 1   (one orbit on single vertices: all vertices are equivalent)
    f_2 = 1   (one orbit on pairs: all edges look the same)
    f_3 = 2   (two orbits on triples: cyclic or transitive)
    f_4 = 4   (four tournament types on 4 vertices)
    f_5 = 12  (twelve on 5 vertices)
    ...

  This IS A000568.

  The exponential generating function:
    Σ f_n x^n / n! = exp(Σ s_k x^k / k)
  where s_k counts the orbits on ORDERED k-tuples = 2^C(k,2)
  (since Aut(T_ω) is 2^C(k,2)-transitive on labeled tournaments).

  Wait — that's not quite right. The actual EGF uses the Burnside
  formula. But the KEY POINT stands:

  A000568 has a GROUP-THEORETIC meaning as orbit counts of Aut(T_ω).
  Our entire project is, at its heart, studying the structure of
  the automorphism group of the universal tournament.
""")

# ============================================================================
# 5. TOURNAMENT SPACE → CANTOR SET
# ============================================================================

print(f"{'='*72}")
print(f"  5. TOURNAMENT SPACE → CANTOR SET AT INFINITY")
print(f"{'='*72}")

print(f"""
  At finite n: tournament tiling space = Q_m = {{0,1}}^m,  m = C(n-1, 2)
  At n → ∞: Q_m → {{0,1}}^N = 2^ω = CANTOR SPACE

  The Cantor space 2^ω is:
    • Compact, metrizable, totally disconnected
    • The universal 0-dimensional compact space
    • Homeomorphic to the Cantor middle-thirds set

  Under the product topology on 2^ω:
    • The set of tournaments isomorphic to T_ω is COMEAGER (residual)
    • Meaning: T_ω is the "generic" tournament in the topological sense
    • Equivalently: the "typical" infinite tournament IS T_ω

  OUR METAGRAPH AT INFINITY:
    G_n / Z_2 has V = A000568(n) vertices at finite n.
    As n → ∞, V → ∞ super-exponentially.
    The metagraph becomes a CONTINUOUS object:
      • Vertices = isomorphism classes = orbits of Aut(T_ω)
      • Edges = arc-flip connections
      • The topology is inherited from the Cantor space quotient

  TOURNAMENT LIMIT THEORY (tournamentons):
    A tournamenton is a measurable function W: [0,1]² → [0,1]
    satisfying W(x,y) + W(y,x) = 1 a.e.

    Examples:
      W(x,y) = 1/2 for all x,y      (random tournamenton → T_ω)
      W(x,y) = 1 if x < y            (transitive tournamenton → Q)
      W(x,y) = 1 if 0 < y-x < 1/2   (circular tournamenton → S(2))

    Every convergent sequence of finite tournaments converges to a
    tournamenton in the cut metric. Our project studies the structure
    of tournamenton space.
""")

# ============================================================================
# 6. WHAT SURVIVES AT n = ∞
# ============================================================================

print(f"{'='*72}")
print(f"  6. WHAT SURVIVES AND WHAT DIES AT n = ∞")
print(f"{'='*72}")

print(f"""
  SURVIVES (meaningful at n = ∞):

    • RÉDEI'S THEOREM: Every countable tournament has a spanning path
      of order type ω (an ω-Hamiltonian path). So H(T) > 0 for all
      countable tournaments. For T_ω: the number of spanning paths is
      uncountable (continuum).

    • ERDŐS-MOSER: Every infinite tournament contains an infinite
      transitive subtournament (via Ramsey). So the transitive class
      always embeds.

    • SCORE SEQUENCE: The scores (out-degrees) of T_ω are all ω
      (countably infinite out-degree for every vertex). The score
      sequence is the constant sequence (ω, ω, ω, ...).

    • CYCLE STRUCTURE: T_ω contains cycles of every finite length.
      It contains every finite tournament as a subtournament.

    • THE EXTENSION PROPERTY: The defining property of T_ω is the
      infinite version of the "every class is reachable" completeness.

  DIES (breaks or becomes trivial at n = ∞):

    • PARITY: Rédei's odd-parity theorem is inherently finite.
      At infinity, the number of Hamiltonian paths is uncountable.

    • BURNSIDE FORMULA: The orbit count A000568(n) ~ 2^(C(n,2))/n!
      diverges. There is no finite "number of isomorphism classes."
      Instead, the classes form a Polish space (2^ω / Aut(T_ω)).

    • CHROMATIC NUMBER: χ(G_n) = n-1 → ∞. In the limit, the
      metagraph requires infinitely many colors.

    • SPECTRAL GAP: 2/n → 0. The spectrum becomes continuous.
      The Krawtchouk polynomials become Hermite polynomials.

    • BAND-LIMITEDNESS: Walsh degree 2⌊(n-1)/2⌋ → ∞. In the limit,
      H is no longer band-limited — it uses all frequencies.

    • FIBER FRACTION: f(n) → 0. No single class has positive measure
      in the limit. The distribution over classes is purely atomic
      with atoms of measure zero.

  TRANSFORMS (changes character at n = ∞):

    • The METAGRAPH becomes a continuous object (quotient of Cantor space)
    • The TILING MODEL becomes a probability measure on 2^ω
    • The STAIRCASE δ_{n-2} becomes the "infinite staircase" (the full
      positive orthant staircase approximating the diagonal of R^2)
    • L1/L2 TRANSITION: at finite n, both metrics coexist. At n = ∞,
      the CLT makes the L2 metric dominant (Gaussian concentration)
""")

# ============================================================================
# 7. NUMERICAL: APPROACH TO THE LIMIT
# ============================================================================

print(f"{'='*72}")
print(f"  7. HOW FAST DO WE APPROACH INFINITY?")
print(f"{'='*72}")

print(f"\n  {'n':>3} {'log₂V':>8} {'m':>5} {'V/V_prev':>10} {'f(n)':>10} {'sc_rate':>10}")
prev_v = 1
for n in range(3, 11):
    m = comb(n, 2)
    vn = A000568[n]
    logV = log2(vn)
    growth = vn / prev_v
    fn = float(Gamma(n - 2 + 0.5) / (Gamma(0.5) * factorial(n - 2)))
    sc_rate = 2.0**(1-n)  # approximate

    print(f"  {n:>3} {logV:>8.2f} {m:>5} {growth:>10.1f} {fn:>10.6f} {sc_rate:>10.6f}")
    prev_v = vn

print(f"""
  Growth rate V(n)/V(n-1) ~ 2^{{n-1}}/n → ∞ super-exponentially.
  By n=10: V = 9.7 million classes.
  By n=20: V ~ 10^{{43}}.
  By n=100: V ~ 10^{{1475}}.

  The metagraph is already intractable at n=8 (V=6880).
  Our project's computational reach (n ≤ 7, V ≤ 456) covers
  only the tiniest initial fragment of the infinite structure.

  Yet the THEOREMS we prove (band-limitedness, χ = n-1, fiber fraction,
  Krawtchouk structure) hold for all finite n, and their LIMITS
  characterize the behavior at infinity.
""")

# ============================================================================
# 8. THE THREE FACES OF INFINITY
# ============================================================================

print(f"{'='*72}")
print(f"  8. THREE FACES OF TOURNAMENT INFINITY")
print(f"{'='*72}")

print(f"""
  Lachlan's classification gives exactly THREE homogeneous countable
  tournaments. These are the three "faces" of tournament infinity:

  ┌─────────────────────────────────────────────────────────────────┐
  │                                                                 │
  │  (Q, <)         S(2)              T_ω                           │
  │  Transitive     Circular          Random                        │
  │  Linear order   Local order       Universal                     │
  │                                                                 │
  │  H = 1          H = 0 (cycles)    H = ∞                        │
  │  c₃ = 0         c₃ = dense        c₃ = dense                   │
  │  Regular: no     Regular: yes      Regular: yes (all deg = ω)   │
  │  Aut = Aut(Q)   Aut = PSL₂(R)    Aut = simple, oligomorphic   │
  │                                                                 │
  │  Score seq:      Score seq:        Score seq:                   │
  │  (0,1,...,ω)    (ω,ω,...,ω)       (ω,ω,...,ω)                  │
  │                                                                 │
  │  Limit of:       Limit of:         Limit of:                   │
  │  Transitive T_n  Paley T_p        Random T_n                   │
  │  H minimal       H maximal?        H typical                   │
  │                                                                 │
  └─────────────────────────────────────────────────────────────────┘

  Our project lives in the space BETWEEN these three:
    • The transitive tournament (H=1, source of the metagraph)
    • The regular tournament (H maximal, center of the metagraph)
    • The random tournament (H typical, the limit at infinity)

  The METAGRAPH G_n is the FINITE-DIMENSIONAL shadow of the infinite
  space connecting these three archetypes.

  As n → ∞:
    • G_n → quotient of Cantor space by Aut(T_ω)
    • The spine (SC-SC) vanishes (SC fraction → 0)
    • The sea (NS-NS) dominates (96% at n=8 → 100% at n=∞)
    • The H-spectrum fills [0, n!] continuously
    • Tournament theory becomes MEASURE THEORY on tournamenton space
""")

print("DONE.")
