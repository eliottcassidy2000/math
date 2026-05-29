#!/usr/bin/env python3
"""
sixteen_dimensions_s93d.py — The 16-dimensional tournament decomposition of transformers
opus-2026-03-20-S93d

THE KEY INSIGHT:

gl(4,R) = 4×4 matrices = 16 dimensions.
This decomposes as:
  ANTISYMMETRIC (skew) 4×4 matrices: dim = C(4,2) = 6  → ACTIVE (tournament)
  SYMMETRIC 4×4 matrices:            dim = C(5,2) = 10 → DARK (cooperation)

The 6 active modes ARE the tournament structure (who beats whom).
The 10 dark modes ARE the cooperation structure (mutual strength).

If the model's output is driven by COMPETITION (tournament/antisymmetric)
but its self-knowledge lives in COOPERATION (symmetric/dark),
this explains Napolitano's finding AND connects to our framework.
"""

from math import comb, sqrt, pi, tanh, atanh
import numpy as np

# =============================================================================
# PART 1: THE DECOMPOSITION gl(4,R) = TOURNAMENT ⊕ COOPERATION
# =============================================================================

print("=" * 70)
print("THE 16-DIMENSIONAL TOURNAMENT DECOMPOSITION")
print("=" * 70)
print()

print("""
  A 4×4 real matrix M has 16 entries. It decomposes UNIQUELY as:

    M = A + S   where  A = (M - M^T)/2  (antisymmetric)
                        S = (M + M^T)/2  (symmetric)

  Dimensions:
    A (antisymmetric 4×4): C(4,2) = 6 independent entries.
    S (symmetric 4×4):     C(4+1,2) = 10 independent entries (including diagonal).

    6 + 10 = 16 = dim(gl(4,R)). ✓

  TOURNAMENT INTERPRETATION:
    A[i][j] = -A[j][i]: the COMPETITION between objects i and j.
      If A[i][j] > 0: object i "beats" j. If < 0: j beats i.
      This IS a weighted tournament on 4 vertices.

    S[i][j] = S[j][i]: the COOPERATION between objects i and j.
      S[i][j] > 0: i and j reinforce each other.
      S[i][i] = self-confidence of object i.
      This IS a weighted undirected graph with self-loops.

  NAPOLITANO'S FINDING (if correct):
    The 6 antisymmetric (tournament) modes → output head → ACTIVE.
    The 10 symmetric (cooperation) modes → self-knowledge → DARK.

    The model COMPETES to generate tokens (tournament: who wins?).
    But it COOPERATES to know if it's right (symmetric: mutual support).
""")

# Demonstrate the decomposition
print("  Example: random 4×4 matrix decomposed into tournament + cooperation:")
print()

np.random.seed(42)
M = np.random.randn(4, 4)
A = (M - M.T) / 2  # antisymmetric (tournament)
S = (M + M.T) / 2  # symmetric (cooperation)

print("  M (full):")
for row in M: print(f"    [{', '.join(f'{x:+.2f}' for x in row)}]")
print()
print("  A (antisymmetric = TOURNAMENT):")
for row in A: print(f"    [{', '.join(f'{x:+.2f}' for x in row)}]")
print()
print("  S (symmetric = COOPERATION):")
for row in S: print(f"    [{', '.join(f'{x:+.2f}' for x in row)}]")
print()

# Verify
recon = A + S
error = np.max(np.abs(recon - M))
print(f"  Reconstruction error: {error:.2e} (should be ~0)")
print(f"  A is antisymmetric: {np.allclose(A, -A.T)}")
print(f"  S is symmetric: {np.allclose(S, S.T)}")

# =============================================================================
# PART 2: WHY 4? THE TOURNAMENT SIZE OF ATTENTION
# =============================================================================

print()
print("=" * 70)
print("WHY 4 OBJECTS? THE NATURAL TOURNAMENT SIZE OF ATTENTION")
print("=" * 70)
print()

print("""
  In a transformer with k normalization layers per block:
    Each normalization creates a 2-dimensional internal space.
    With k=2 norms (standard): 2 × 2 = 4 internal objects.
    gl(2k, R) = gl(4, R) for k=2.

  These 4 objects are:
    1. Query representation (what am I looking for?)
    2. Key representation (what do I offer?)
    3. Value representation (what do I provide?)
    4. Output representation (what do I produce?)

  The TOURNAMENT among them:
    Q beats K: query dominates key selection (attention is sharp)
    K beats V: key selection dominates value retrieval
    V beats O: value retrieval dominates output
    O beats Q: output feeds back to modify the query

  This 4-tournament has C(4,2) = 6 arcs = 6 active modes.
  It controls HOW attention flows through the transformer.

  The COOPERATION among them:
    Q·K: query-key similarity (attention weights)
    Q·V: query-value alignment (direct retrieval)
    Q·O: query-output consistency (did we answer the question?)
    K·V: key-value coherence (is the stored information consistent?)
    K·O: key-output relevance (did we use the right keys?)
    V·O: value-output fidelity (was the value accurately transmitted?)
    Q·Q, K·K, V·V, O·O: self-consistency of each representation

  10 symmetric interactions = 10 dark modes.
  These track WHETHER the computation is internally consistent.
""")

# =============================================================================
# PART 3: THE FORMAL GROUP ON THE 4-TOURNAMENT
# =============================================================================

print("=" * 70)
print("THE FORMAL GROUP ON THE 4-TOURNAMENT")
print("=" * 70)
print()

print("""
  The 4-tournament has T(4) = 4 isomorphism classes:
    1. Transitive (score 3,2,1,0): one clear winner
    2. Near-transitive A (score 3,1,1,1): strong winner, tied rest
    3. Near-transitive B (score 2,2,2,0): no clear winner, one clear loser
    4. Diamond (score 2,2,1,1): balanced competition

  The formal group F(x,y) = (x+y)/(1+xy) acts on the COUPLING SPACE
  of the tournament. Each arc has coupling x_ij ∈ (-1, 1).

  The RAPIDITY of each arc: ρ_ij = arctanh(x_ij).
  The total rapidity of the tournament: Σ ρ_ij.

  For the 4-tournament with 6 arcs:
    Total rapidity space = R^6 (one rapidity per arc).
    S_4 acts on this space (permuting vertex labels).
    The quotient R^6 / S_4 has dimension... complex to compute.

  The KEY: the formal group linearizes the tournament dynamics.
  In the linearized (rapidity) space, the 6 tournament modes
  become ADDITIVE. Evidence aggregation = addition of rapidities.

  This is why FormalRank works: arctanh converts the 6 active modes
  from multiplicative (tournament) to additive (rapidity) space.
""")

# =============================================================================
# PART 4: WHAT THE 10 DARK MODES COMPUTE
# =============================================================================

print("=" * 70)
print("WHAT THE 10 DARK MODES COMPUTE")
print("=" * 70)
print()

print("""
  The 10 symmetric modes S[i][j] = S[j][i] form a METRIC on the 4 objects.

  A metric tells you DISTANCES between objects:
    S[i][j] large: i and j are "close" (cooperating, reinforcing)
    S[i][j] small: i and j are "far" (independent, unrelated)
    S[i][j] negative: i and j are "opposing" (contradicting)

  The diagonal S[i][i] = NORM of each object:
    S[i][i] large: object i is "strong" (high self-confidence)
    S[i][i] small: object i is "weak" (low self-confidence)

  THE METRIC DETECTS INCONSISTENCY:
    If Q and K are "close" (S[Q][K] large) but V and O are "far" (S[V][O] small):
    the model is paying attention to the right things but not producing
    the right output. This is a HALLUCINATION signature.

  TOURNAMENT INTERPRETATION:
    The 6 active modes say: "here is the ranking" (who beats whom).
    The 10 dark modes say: "here is the confidence" (how reliable is the ranking).

    The model can produce a ranking (via tournament/active modes) that it
    KNOWS is wrong (via metric/dark modes). The dark modes are the model's
    UNCERTAINTY ESTIMATE, orthogonal to the ranking itself.

  THIS EXPLAINS NAPOLITANO'S FINDING:
    Dark modes alone achieve 93.86% on ARC-Challenge because:
    - ARC asks "is this answer correct?"
    - The answer's correctness lives in the METRIC (dark modes)
    - The actual answer lives in the TOURNAMENT (active modes)
    - You don't need to know the answer to know if you're right!
""")

# =============================================================================
# PART 5: THE 6+10 DECOMPOSITION IN OUR THREE-CHANNEL FRAMEWORK
# =============================================================================

print("=" * 70)
print("THE 6+10 DECOMPOSITION vs OUR THREE CHANNELS")
print("=" * 70)
print()

print("""
  Our BoostRanker decomposes each comparison into:
    INERT (2):    who wins? → binary ranking
    RAMIFIED (3): by how much? → confidence/margin
    SPLIT (7):    surprise? → anomaly detection

  The 6+10 decomposition:
    ACTIVE (6):   the tournament → who wins (INERT + RAMIFIED)
    DARK (10):    the metric → surprise and self-knowledge (SPLIT + more)

  MAPPING:
    Active = INERT + RAMIFIED (the competitive channels)
    Dark = SPLIT + METRIC (the cooperative channels)

    The SPLIT channel detects intransitivity (cycles in the tournament).
    The METRIC detects inconsistency (distance violations in the cooperation).
    Both are "dark" in the sense that they're invisible to the output head.

  DIMENSIONAL ACCOUNTING:
    INERT:    1 dimension (binary who-wins)
    RAMIFIED: 5 dimensions (margin for each of C(4,2)=6 pairs minus 1 for binary)
    SPLIT:    4 dimensions (surprise for each of the 4 objects)
    METRIC:   6 dimensions (pairwise distances for C(4,2)=6 pairs)

    1 + 5 + 4 + 6 = 16 ✓

    So: INERT(1) + RAMIFIED(5) = 6 active.
        SPLIT(4) + METRIC(6) = 10 dark.

  THIS IS A REFINEMENT OF NAPOLITANO:
    His 6+10 is our (INERT+RAMIFIED) + (SPLIT+METRIC).
    Our decomposition is FINER: 1+5+4+6 vs his 6+10.
    Each of our four channels has a specific mathematical meaning.
""")

# =============================================================================
# PART 6: PRACTICAL IMPLICATIONS
# =============================================================================

print("=" * 70)
print("PRACTICAL IMPLICATIONS")
print("=" * 70)
print()

print("""
  IF the 16-dimensional decomposition is correct:

  1. HALLUCINATION DETECTION (improved):
     Instead of just tracking the logit gap (our TournamentConfidence),
     decompose the hidden state into tournament (6) + metric (10).
     Check if the METRIC is inconsistent with the TOURNAMENT.
     If yes: the model knows it's wrong but is generating anyway.

  2. SELF-CONSISTENCY TRAINING:
     Add a loss term that aligns dark modes with active modes.
     The dark metric should AGREE with the tournament ranking:
     if the tournament says A > B, the metric should say d(A,answer) < d(B,answer).
     Misalignment = the model is learning to hallucinate.

  3. DARK MODE EXTRACTION (zero-parameter):
     For any 4-object interaction in the transformer:
       M = hidden state reshaped as 4×4 (or projected to 4 dims)
       A = (M - M^T)/2  → tournament modes (6 dims, active)
       S = (M + M^T)/2  → metric modes (10 dims, dark)
     Check S for inconsistency. No probes needed!

  4. ARCHITECTURE REDESIGN:
     Current: one attention head computes QK^T (tournament) and uses V (cooperation).
     Better: SEPARATE tournament and metric computation.
       Tournament head: computes antisymmetric interactions (who beats whom).
       Metric head: computes symmetric interactions (how consistent am I?).
     The metric head outputs a CONFIDENCE SCORE, not a token prediction.
     This naturally produces calibrated uncertainty.

  5. ADAPTIVE NORM COUNT:
     Napolitano says: k normalizations → gl(2k, R).
     k=2 (standard): gl(4,R), 16 dims, 6+10.
     k=3 (his proposal): gl(6,R), 36 dims, 15+21.
     k=4 (further): gl(8,R), 64 dims, 28+36.

     More norms → more dark modes → better self-knowledge.
     The COST: more parameters (+72.8M for k=3 per his estimate).
     The BENEFIT: predicted 98.7% on ARC-Challenge.

     OUR PERSPECTIVE: for gl(2k,R), the dark sector has
       D_k = C(2k+1, 2) = k(2k+1)
     symmetric modes. Increasing the norm count k -> k+1 changes
       gl(2k,R) -> gl(2k+2,R)
     and adds
       D_{k+1} - D_k = C(2k+3, 2) - C(2k+1, 2) = 4k + 3
     new dark modes. For k=2 -> 3 this is 11, matching 10 -> 21.

     The MARGINAL VALUE of each dark mode should decrease
     (like the layer-d correction in the k-periodicity tower).

     PREDICTION: accuracy scales as 1 - C × q^N where N = dark modes
     and q = 1 - 1/N ≈ 0.881 (from the scaling law).
     This gives diminishing returns: going from 10→21 dark modes
     is much more valuable than from 21→36.
""")

# Compute predicted accuracy for different numbers of dark modes
print("  Predicted accuracy by dark mode count:")
print(f"  {'k (norms)':>9} | {'algebra':>8} | {'dark modes':>10} | {'predicted acc':>13}")
print("  " + "-" * 50)

for k in [1, 2, 3, 4, 5]:
    n = 2*k
    total = n*n
    antisym = comb(n, 2)
    sym = total - antisym
    # Scaling: errors = 209 * 0.881^N
    dark = sym
    errors = 209 * 0.881**dark
    acc = max(0, (1172 - errors) / 1172 * 100)  # ARC has 1172 questions
    print(f"  {k:9d} | gl({n},R) | {dark:10d} | {acc:12.1f}%")

# =============================================================================
# PART 7: THE NUMBER 16 IN OUR FRAMEWORK
# =============================================================================

print()
print("=" * 70)
print("THE NUMBER 16 IN TOURNAMENT THEORY")
print("=" * 70)
print()

print("""
  16 appears in our framework at several points:

  1. 16 = 2^4 = the number of labeled tournaments on n=4 vertices.
     Wait: 2^{C(4,2)} = 2^6 = 64, not 16. So this doesn't work directly.

  2. 16 = C(4,2) + C(5,2) = 6 + 10 = arcs at n=4 + arcs at n=5.
     The 16 dimensions span TWO tournament levels simultaneously.

  3. 16 = the dimension of the Clifford algebra Cl(4) representation.
     Cl(4) = M(2, H) = 2×2 quaternionic matrices.
     This connects to Bott periodicity at period 8:
     Cl(8) = M(16, R) (16×16 real matrices).

  4. 15 = 16 - 1 = C(6, 2) = arcs at the PIVOT n=6.
     The pivot tournament has 15 arcs and eigenvalue denominator 15 = 3×5.
     16 = C(6,2) + 1 = arcs + identity.

  5. The SEDENION boundary: 16-dimensional hypercomplex numbers (sedenions).
     At 16-ons: the division algebra property FAILS.
     No norms, no nice multiplication. The formal group "dies" here.
     2→4→8→16: real→complex→quaternion→octonion→sedenion(broken).

  THE DEEPEST CONNECTION:
    The 16-dimensional structure of gl(4,R) = tournament + cooperation
    might be the ALGEBRAIC reason why transformers with 2 norms per block
    sit at the SEDENION BOUNDARY — the edge of where nice algebraic
    structure (division algebras) can exist.

    Going to gl(6,R) (3 norms, 36 dims) crosses into the SEDENION regime
    where the formal group structure is richer but more fragile.
    This might explain why adding the third norm helps accuracy
    (more dark modes = more self-knowledge) but with diminishing returns
    (the sedenion boundary imposes structural limits).
""")

print("opus-2026-03-20-S93d: THE 16-DIMENSIONAL TOURNAMENT DECOMPOSITION")
