#!/usr/bin/env python3
"""
erdos_gallai_s188.py — Erdős-Gallai and Landau in Tournament Theory
opus-2026-03-22-S188

Erdős-Gallai characterizes graphic sequences (degree sequences of graphs).
Landau characterizes score sequences (out-degree sequences of tournaments).

Both are FEASIBILITY conditions on sequences. Together they control:
  - Which nodes can exist in G_n (Landau → score classes → iso classes)
  - Which edges can exist in G_n (Erdős-Gallai on G_n's own degree sequence)
  - The structure of the score simplex (Landau polytope)
"""

from math import comb, log2
from itertools import combinations
from collections import Counter

print("=" * 72)
print("  ERDŐS-GALLAI AND LANDAU IN TOURNAMENT THEORY")
print("  opus-2026-03-22-S188")
print("=" * 72)
print()

# ==========================================================================
# PART 1: LANDAU'S THEOREM — WHICH SCORE SEQUENCES ARE VALID?
# ==========================================================================

print("PART 1: LANDAU'S THEOREM")
print("-" * 60)
print()

def is_landau(s):
    """Check if sorted sequence s is a valid tournament score sequence.
    Landau: Σ_{i=1}^k s_i ≥ C(k,2) for all k, with equality at k=n."""
    n = len(s)
    s_sorted = sorted(s)
    for k in range(1, n+1):
        if sum(s_sorted[:k]) < comb(k, 2):
            return False
    return sum(s_sorted) == comb(n, 2)

def enumerate_score_sequences(n):
    """Generate all valid tournament score sequences on n vertices."""
    # Scores are non-negative integers summing to C(n,2) = n(n-1)/2
    # Each score is between 0 and n-1
    # Must satisfy Landau's conditions
    target = comb(n, 2)
    valid = []

    def backtrack(remaining, min_val, current, depth):
        if depth == n:
            if remaining == 0:
                if is_landau(current):
                    valid.append(tuple(current))
            return
        for v in range(min_val, n):
            if remaining - v < 0:
                break
            if remaining - v > (n - depth - 1) * (n - 1):
                continue
            current.append(v)
            backtrack(remaining - v, v, current, depth + 1)
            current.pop()

    backtrack(target, 0, [], 0)
    return valid

print(f"  Landau's theorem: sequence (s₁ ≤ ... ≤ sₙ) is a tournament score iff")
print(f"    Σ{{i=1}}^k s_i ≥ C(k,2) for k=1..n, with Σ s_i = C(n,2)")
print()

for n in range(2, 9):
    seqs = enumerate_score_sequences(n)
    total_seqs = len(seqs)
    # How many integer partitions of C(n,2) into n parts ≤ n-1?
    # vs how many satisfy Landau?
    print(f"  n={n}: {total_seqs} valid score sequences, "
          f"C(n,2)={comb(n,2)}")
    if n <= 6:
        for s in seqs:
            print(f"    {s}")
    print()

# This IS the A000571 sequence
print(f"  Score sequence counts: {[len(enumerate_score_sequences(n)) for n in range(2, 9)]}")
print(f"  = A000571: 1, 2, 4, 9, 22, 59, 167")
print()

# ==========================================================================
# PART 2: THE LANDAU POLYTOPE
# ==========================================================================

print()
print("PART 2: THE LANDAU POLYTOPE")
print("-" * 60)
print()

print("""  The set of valid score sequences forms a POLYTOPE in R^n:
    P_n = {s ∈ Z^n : s sorted, Σs = C(n,2), Landau conditions}

  The VERTICES of this polytope are the extreme score sequences.
  The INTERIOR points are the "generic" score sequences.

  At n=5:
    Vertices of P_5: (0,1,2,3,4) and (2,2,2,2,2)
    These are the TRANSITIVE and REGULAR score sequences.
    Everything else is "between" these extremes.

  The Landau polytope IS the score simplex from S174.
  Its integer points = the valid score sequences.
  The number of integer points = A000571(n).
""")

# Compute the "distance" of each score from the regular score
n = 5
seqs = enumerate_score_sequences(n)
regular = tuple([n//2] * n) if n % 2 == 1 else None

print(f"  SCORE SEQUENCES AT n={n} WITH DISTANCE FROM REGULAR:")
print(f"  {'score':>16s}  {'S₂':>4s}  {'dist_to_reg':>11s}  {'H range':>10s}")
print(f"  {'─'*16}  {'─'*4}  {'─'*11}  {'─'*10}")

for s in seqs:
    s2 = sum(v*v for v in s)
    if n % 2 == 1:
        reg = (n-1)//2
        dist = sum(abs(v - reg) for v in s)
    else:
        dist = 0  # harder for even n
    # H range from our formula H = 1 + n(n-1)(2n-1)/6 - S₂ (+ correction)
    base = 1 + n*(n-1)*(2*n-1)//6
    h_approx = base - s2

    print(f"  {str(s):>16s}  {s2:4d}  {dist:11d}  ~{h_approx}")

print()

# ==========================================================================
# PART 3: ERDŐS-GALLAI ON G_n's DEGREE SEQUENCE
# ==========================================================================

print()
print("PART 3: ERDŐS-GALLAI ON G_n's OWN DEGREE SEQUENCE")
print("-" * 60)
print()

def is_erdos_gallai(d):
    """Check if non-increasing sequence d is graphic."""
    n = len(d)
    d_sorted = sorted(d, reverse=True)
    if sum(d_sorted) % 2 != 0:
        return False
    for k in range(1, n+1):
        lhs = sum(d_sorted[:k])
        rhs = k * (k-1) + sum(min(d_sorted[i], k) for i in range(k, n))
        if lhs > rhs:
            return False
    return True

# Degree sequences of G_n (from S177)
gn_degrees = {
    3: [1, 1],
    4: [2, 2, 3, 3],
    5: [2, 3, 3, 3, 4, 6, 6, 6, 6, 7, 7, 7],
}

for n_g in [3, 4, 5]:
    d = sorted(gn_degrees[n_g], reverse=True)
    eg = is_erdos_gallai(d)
    print(f"  G_{n_g} degree sequence: {d}")
    print(f"    Erdős-Gallai valid: {eg}")
    print(f"    Sum = {sum(d)} (must be even: {sum(d) % 2 == 0})")
    print(f"    Max degree = {d[0]}, min = {d[-1]}")
    print()

# ==========================================================================
# PART 4: THE LANDAU-ERDŐS-GALLAI CONNECTION
# ==========================================================================

print()
print("PART 4: THE DEEP CONNECTION")
print("-" * 60)
print()

print("""  Landau conditions on TOURNAMENT SCORES control which
  VERTICES exist in G_n (each vertex = one score class or iso class).

  Erdős-Gallai conditions on G_n's DEGREE SEQUENCE control
  which EDGES can exist.

  TOGETHER they constrain the entire meta-graph G_n.

  THE PARALLEL:

  TOURNAMENT LEVEL:
    Landau: which score sequences are valid?
    Answer: those satisfying Σ_{i=1}^k s_i ≥ C(k,2).
    This is a CONVEX polytope condition.

  META-GRAPH LEVEL:
    Erdős-Gallai: which degree sequences for G_n are valid?
    Answer: those satisfying the EG conditions.
    G_n's degree sequence IS graphic (it comes from a real graph).

  THE HIERARCHY:
    Level 0: Landau conditions → valid tournaments → nodes of G_n
    Level 1: Flip structure → edges of G_n
    Level 2: Erdős-Gallai → G_n's degree sequence is consistent

  Each level constrains the next. Landau determines the nodes;
  the flip structure determines the edges; EG verifies consistency.
""")

# ==========================================================================
# PART 5: LANDAU GAPS AND FORBIDDEN SCORES
# ==========================================================================

print()
print("PART 5: LANDAU GAPS — FORBIDDEN SCORE SEQUENCES")
print("-" * 60)
print()

# Which score-like sequences FAIL Landau?
print(f"  WHICH PARTITIONS OF C(n,2) INTO n PARTS FAIL LANDAU?")
print()

for n in range(3, 7):
    target = comb(n, 2)
    all_partitions = []
    def gen_parts(remaining, min_val, current, depth, max_val):
        if depth == n:
            if remaining == 0:
                all_partitions.append(tuple(current))
            return
        for v in range(min_val, min(max_val + 1, remaining + 1)):
            if remaining - v < 0: break
            current.append(v)
            gen_parts(remaining - v, v, current, depth + 1, max_val)
            current.pop()

    gen_parts(target, 0, [], 0, n-1)
    valid = [p for p in all_partitions if is_landau(p)]
    invalid = [p for p in all_partitions if not is_landau(p)]

    print(f"  n={n}: {len(all_partitions)} partitions, "
          f"{len(valid)} valid (Landau), {len(invalid)} invalid")
    if invalid and n <= 5:
        for p in invalid[:5]:
            # Which Landau condition fails?
            p_sorted = sorted(p)
            for k in range(1, n+1):
                if sum(p_sorted[:k]) < comb(k, 2):
                    print(f"    INVALID: {p} — fails at k={k}: "
                          f"Σ={{sum(p_sorted[:k])}} < C({k},2)={comb(k,2)}")
                    break
    print()

# ==========================================================================
# PART 6: THE SCORE SEQUENCE AS A CODE
# ==========================================================================

print()
print("PART 6: SCORE SEQUENCES AS CODES")
print("-" * 60)
print()

print("""  The Landau conditions define a LINEAR CODE over Z:
    {s ∈ Z^n : s sorted, Σs = C(n,2), Σ_{i=1}^k s_i ≥ C(k,2) ∀k}

  This is a POLYTOPE CODE: the codewords are the integer lattice
  points inside the Landau polytope P_n.

  The CODE RATE = log₂(#codewords) / n:
    n=3: log₂(2)/3 = 0.33 bits/symbol
    n=4: log₂(4)/4 = 0.50 bits/symbol
    n=5: log₂(9)/5 = 0.63 bits/symbol
    n=6: log₂(22)/6 = 0.74 bits/symbol
    n=7: log₂(59)/7 = 0.84 bits/symbol

  The rate INCREASES with n → the Landau polytope becomes
  RELATIVELY LARGER as n grows (fewer constraints binding).

  The MINIMUM DISTANCE of this code:
  Two score sequences differ by at least one score value.
  Min distance = minimum L₁ distance between valid sequences.
""")

for n in range(3, 8):
    seqs = enumerate_score_sequences(n)
    k = len(seqs)
    rate = (log2(k) if k > 0 else 0) / n if n > 0 else 0
    print(f"  n={n}: {k} codewords, rate = {rate:.3f} bits/symbol")

print()

# Minimum L1 distance between score sequences
for n in range(3, 7):
    seqs = enumerate_score_sequences(n)
    min_dist = float('inf')
    for i in range(len(seqs)):
        for j in range(i+1, len(seqs)):
            dist = sum(abs(a-b) for a, b in zip(seqs[i], seqs[j]))
            min_dist = min(min_dist, dist)
    print(f"  n={n}: min L₁ distance between score sequences = {min_dist}")

print()

# ==========================================================================
# PART 7: ERDŐS-GALLAI PREDICTS G_n STRUCTURE
# ==========================================================================

print()
print("PART 7: WHAT ERDŐS-GALLAI PREDICTS FOR G_n")
print("-" * 60)
print()

print("""  For a graph on N vertices with degree sequence d₁ ≥ ... ≥ d_N:
    EG condition: Σ_{i=1}^k d_i ≤ k(k-1) + Σ_{i=k+1}^N min(d_i, k)

  For G_n, the degree sequence encodes how CONNECTED each iso class is.
  HIGH degree → many structural transitions available.
  LOW degree → few transitions → structurally "rigid."

  THE DEGREE SEQUENCE OF G_n IS A TOURNAMENT INVARIANT:
    It's determined by the tournament structure (iso classes + flip graph).
    It satisfies EG (it's a real graph).
    Its shape tells us about the TOPOLOGY of tournament space.

  PREDICTIONS FROM THE DEGREE SEQUENCE:
    Degree sequence of G_5: [7,7,7,6,6,6,6,4,3,3,3,2]
    Sum = 60 = 2 × 30 (30 edges).
    Mean degree = 5.0.
    The sequence is UNIMODAL: most classes have degree 6-7,
    a few have degree 2-3 (the high-symmetry extremes).
""")

# Check: does the degree sequence of G_n satisfy any nice formula?
d5 = sorted(gn_degrees[5], reverse=True)
print(f"  G_5 degree sequence: {d5}")
print(f"  Sorted scores of G_5's degree-tournament would need Landau check.")
print()

# Does the degree sequence relate to tournament H values?
# From S176: corr(degree, |Aut|) ≈ -0.94
# From S177: max_degree = 7, achieved by |Aut|=1 classes
print(f"  DEGREE vs |Aut| at n=5:")
# Class data from S170
class_info_5 = [
    (1, 1, 120, 6), (3, 3, 40, 3), (3, 3, 40, 4), (3, 3, 40, 3),
    (5, 1, 120, 7), (5, 1, 120, 7), (9, 1, 120, 6), (9, 1, 120, 7),
    (11, 1, 120, 6), (13, 1, 120, 6), (15, 3, 40, 3), (15, 5, 24, 2),
]
# (H, |Aut|, size, degree)

print(f"  {'H':>3s}  {'|Aut|':>5s}  {'size':>5s}  {'degree':>6s}  {'size×degree':>11s}")
print(f"  {'─'*3}  {'─'*5}  {'─'*5}  {'─'*6}  {'─'*11}")
for h, aut, size, deg in class_info_5:
    print(f"  {h:3d}  {aut:5d}  {size:5d}  {deg:6d}  {size*deg:11d}")

print()
print(f"  NOTE: size × degree is ROUGHLY proportional to class_size × C(n,2)")
print(f"  because total outgoing transitions = size × C(n,2) (each tournament")
print(f"  has C(n,2) possible flips), and degree = # distinct external classes reached.")
print()

# ==========================================================================
# SUMMARY
# ==========================================================================

print()
print("=" * 72)
print("  SUMMARY: ERDŐS-GALLAI AND LANDAU UNIFIED")
print("=" * 72)
print()

from math import log2

print(f"""
  TWO THEOREMS, ONE FRAMEWORK:

  LANDAU (tournaments):
    Characterizes valid score sequences.
    Controls which NODES exist in G_n.
    A000571(n) = number of valid score sequences.
    The Landau polytope is a convex body in R^n.
    Rate: {log2(len(enumerate_score_sequences(5)))/5:.2f} bits/symbol at n=5.

  ERDŐS-GALLAI (graphs):
    Characterizes valid degree sequences.
    Applied to G_n: verifies its degree sequence is consistent.
    G_n's degree sequence: {sorted(gn_degrees[5], reverse=True)}
    The degree sequence encodes the topology of tournament space.

  THE HIERARCHY:
    Landau → valid scores → nodes of G_n → degree sequence → Erdős-Gallai
    Each level constrains the next.

  THE CODE-THEORETIC VIEW:
    Score sequences = codewords of the Landau code.
    Tournaments = messages in the tournament code.
    Syndrome = score sequence (Landau codeword).
    The OCR (97%) = syndrome decoding success rate.

  FOR G_n:
    The degree of a class ≈ C(n,2) × (1 - self-loop fraction) / (# neighbors)
    High |Aut| → high self-loop → low degree (formula from S182: 2/n at source)
    Low |Aut| → low self-loop → high degree

  THE LANDAU POLYTOPE IS THE "CONFIGURATION SPACE" OF G_n's NODES.
  The Erdős-Gallai conditions are the "CONSISTENCY CHECK" on G_n's edges.
  Together they define the FEASIBLE REGION for the meta-graph.
""")

print("opus-2026-03-22-S188")
