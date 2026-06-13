#!/usr/bin/env python3
"""
imaginary_units_s138.py — What Imaginary Units Actually Are
opus-2026-03-21-S138

The user said: "please understand what imaginary units are."

This is a call to stop making analogies and START UNDERSTANDING.

An imaginary unit e is an element of an algebra satisfying e² = -1.
What does this MEAN? What does it DO? And what is its precise
relationship to tournaments?

Let me think carefully.
"""

from math import comb, sqrt, pi, cos, sin
from itertools import permutations

def tournament_adj(n, bits):
    adj = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx):
                adj[i][j] = 1
            else:
                adj[j][i] = 1
            idx += 1
    return adj

def H_dp(adj, n):
    dp = [0]*((1<<n)*n)
    for v in range(n): dp[(1<<v)*n+v] = 1
    for S in range(1, 1<<n):
        for v in range(n):
            if not (S&(1<<v)): continue
            val = dp[S*n+v]
            if val == 0: continue
            for u in range(n):
                if S&(1<<u): continue
                if adj[v][u]: dp[(S|(1<<u))*n+u] += val
    full = (1<<n)-1
    return sum(dp[full*n+v] for v in range(n))

print("=" * 72)
print("  WHAT IMAGINARY UNITS ACTUALLY ARE")
print("  opus-2026-03-21-S138")
print("=" * 72)
print()

# ==========================================================================
# PART 1: WHAT e² = -1 MEANS
# ==========================================================================

print("PART 1: WHAT e² = -1 MEANS")
print("-" * 60)
print()

print("""  An imaginary unit e satisfies e² = -1.

  This is NOT mysterious. It means:

  APPLYING e TWICE NEGATES.

  If you think of e as an OPERATION (a transformation),
  then doing it twice gives you the OPPOSITE of where you started.

  In C: multiplying by i rotates the plane by 90°.
        Doing it twice rotates by 180° = negation.
        That's ALL that i² = -1 says: half of negation.

  In H: multiplying by i, j, or k EACH rotates in a different
        plane of 4D space by 90°. Each is "half a negation"
        in its own plane.

  The REAL CONTENT of e² = -1 is:

    e is a SQUARE ROOT OF NEGATION.

  Negation is the simplest nontrivial operation (order 2).
  An imaginary unit is its SQUARE ROOT (order 4).
  That's why i⁴ = (i²)² = (-1)² = 1: four applications
  of the half-negation return to the identity.
""")

# ==========================================================================
# PART 2: HOW IMAGINARY UNITS MULTIPLY — THE 3-CYCLE
# ==========================================================================

print()
print("PART 2: HOW IMAGINARY UNITS MULTIPLY — THE 3-CYCLE")
print("-" * 60)
print()

print("""  In the QUATERNIONS, the three imaginary units i, j, k satisfy:

    ij = k,   jk = i,   ki = j     (CYCLIC: i → j → k → i)
    ji = -k,  kj = -i,  ik = -j    (ANTI-CYCLIC: reverse = negate)

  The cyclic rule ij = k, jk = i, ki = j IS a directed 3-cycle.
  The direction determines the sign of the product.

  Going WITH the cycle: product is POSITIVE.
  Going AGAINST the cycle: product is NEGATIVE.

  This is EXACTLY how a tournament 3-cycle works:
    a → b → c → a   (directed cycle)
    Going with: a beats b, b beats c, c beats a
    Going against: each relation is reversed

  BUT THERE IS A CRUCIAL DIFFERENCE:

  In quaternions: ij = +k and ji = -k. BOTH products EXIST.
    The sign tells you the DIRECTION along the cycle.

  In tournaments: a→b is an ARC. b→a is the OPPOSITE arc.
    They can't BOTH exist. One or the other, never both.

  So the tournament 3-cycle is the ABSOLUTE VALUE of the
  quaternion product rule. Tournaments see the CYCLE STRUCTURE
  (which triples form a cycle) but not the SIGN (which direction).

  The sign IS the imaginary unit. The cycle IS the tournament.
  The imaginary unit ADDS the sign to the cycle.
""")

# Demonstrate with explicit quaternion multiplication
print("  QUATERNION MULTIPLICATION TABLE:")
print()
# Represent quaternions as (a, b, c, d) for a + bi + cj + dk
# Multiplication: (a,b,c,d)(e,f,g,h)
# = (ae-bf-cg-dh, af+be+ch-dg, ag-bh+ce+df, ah+bg-cf+de)

units = {'1': (1,0,0,0), 'i': (0,1,0,0), 'j': (0,0,1,0), 'k': (0,0,0,1)}

def qmul(p, q):
    a,b,c,d = p
    e,f,g,h = q
    return (a*e-b*f-c*g-d*h,
            a*f+b*e+c*h-d*g,
            a*g-b*h+c*e+d*f,
            a*h+b*g-c*f+d*e)

def qname(q):
    a,b,c,d = q
    parts = []
    if a == 1: return "+1"
    if a == -1: return "-1"
    if b == 1: parts.append("+i")
    elif b == -1: parts.append("-i")
    if c == 1: parts.append("+j")
    elif c == -1: parts.append("-j")
    if d == 1: parts.append("+k")
    elif d == -1: parts.append("-k")
    if not parts: return "0"
    result = "".join(parts)
    if result.startswith("+"): result = " " + result
    return result

print(f"     {'1':>4s}  {'i':>4s}  {'j':>4s}  {'k':>4s}")
for name1, u1 in [('1', units['1']), ('i', units['i']),
                   ('j', units['j']), ('k', units['k'])]:
    row = f"  {name1}  "
    for name2, u2 in [('1', units['1']), ('i', units['i']),
                       ('j', units['j']), ('k', units['k'])]:
        prod = qmul(u1, u2)
        row += f"{qname(prod):>4s}  "
    print(row)

print()
print(f"  THE DIRECTED 3-CYCLE: i → j → k → i")
print(f"    ij = +k  (WITH the cycle)")
print(f"    ji = -k  (AGAINST the cycle)")
print(f"    The SIGN encodes the DIRECTION.")
print()

# ==========================================================================
# PART 3: THE OCTONION IMAGINARY UNITS — 7 = FANO
# ==========================================================================

print()
print("PART 3: THE SEVEN OCTONION IMAGINARY UNITS")
print("-" * 60)
print()

# Octonion multiplication table (one standard convention)
# 7 imaginary units: e1, e2, ..., e7
# 7 Fano lines: each gives a quaternionic subalgebra
# Standard Fano lines: {1,2,4}, {2,3,5}, {3,4,6}, {4,5,7}, {5,6,1}, {6,7,2}, {7,1,3}
fano_lines = [
    (1, 2, 4),
    (2, 3, 5),
    (3, 4, 6),
    (4, 5, 7),
    (5, 6, 1),
    (6, 7, 2),
    (7, 1, 3),
]

print(f"  The octonions O have 7 imaginary units e₁, ..., e₇.")
print(f"  Each satisfies eᵢ² = -1.")
print()
print(f"  Their multiplication is governed by the FANO PLANE PG(2,F₂):")
print(f"  7 points, 7 lines, 3 points per line, 3 lines per point.")
print()
print(f"  Each Fano line {{a, b, c}} gives a DIRECTED 3-CYCLE:")
print(f"    eₐ · eᵦ = eᶜ  (with the cycle)")
print(f"    eᵦ · eₐ = -eᶜ (against the cycle)")
print()

print(f"  FANO LINES (= directed 3-cycles of imaginary units):")
for i, (a, b, c) in enumerate(fano_lines):
    print(f"    Line {i+1}: {{e{a}, e{b}, e{c}}} → e{a}·e{b} = e{c}")

print()
print(f"  WHAT MAKES OCTONIONS NON-ASSOCIATIVE:")
print(f"    In quaternions: i(jk) = i·i = -1 = (ij)k = k·k = -1 ✓")
print(f"    In octonions: e₁(e₂·e₃) ≠ (e₁·e₂)·e₃ for some triples")
print(f"    because the Fano plane doesn't have a CONSISTENT global order.")
print()
print(f"  A quaternionic subalgebra (one Fano line) IS associative.")
print(f"  But TWO DIFFERENT Fano lines INTERFERE: the order of operations")
print(f"  matters because there's no consistent way to parenthesize.")
print()
print(f"  This is EXACTLY the tournament situation:")
print(f"  One 3-cycle is fine (associative within the triple).")
print(f"  Two 3-cycles that SHARE A VERTEX create interference")
print(f"  (the shared vertex sees conflicting orderings).")
print(f"  This is the CONFLICT GRAPH Ω(T).")
print()

# ==========================================================================
# PART 4: WHAT THE IMAGINARY UNIT DOES TO A TOURNAMENT
# ==========================================================================

print()
print("PART 4: WHAT AN IMAGINARY UNIT DOES TO A TOURNAMENT")
print("-" * 60)
print()

print("""  An imaginary unit e satisfies e² = -1.
  It is a HALF-NEGATION: a 90° rotation.

  In tournament theory, NEGATION is COMPLEMENTATION:
    T^c = the tournament with all arcs reversed.
    (T^c)^c = T  (negation squared = identity)
    H(T^c) = H(T)  (negation preserves H)

  So the tournament complement T → T^c is multiplication by -1.

  QUESTION: Is there a tournament operation that is "half" of complement?
  An operation σ such that σ(σ(T)) = T^c?

  This would be the TOURNAMENT IMAGINARY UNIT.

  For it to exist: σ² = complement, σ⁴ = identity.
  The operation σ has ORDER 4, just like i (since i⁴ = 1).

  CANDIDATE: The Vitali atom!
    The Vitali atom reverses 4 arcs.
    Applying it again reverses 4 more arcs.
    After two applications: 8 arcs reversed... but C(n,2) total arcs.
    This is NOT complement unless 8 = C(n,2), i.e. n=4 only.

  So the Vitali atom is NOT literally σ² = complement.
  But it IS a square root of SOMETHING:
""")

# Let's check: at n=5 (10 arcs), complement reverses all 10 arcs.
# A "half-complement" would reverse 5 arcs.
# At n=3 (3 arcs), complement reverses all 3.
# A "half-complement" would reverse 1.5 arcs — impossible!

# BUT: at n=4 (6 arcs), complement reverses all 6.
# A "half-complement" reverses 3.

n = 4
print(f"  At n={n}: C({n},2) = {comb(n,2)} arcs.")
print(f"  Complement reverses all {comb(n,2)} arcs.")
print(f"  A 'half-complement' reverses {comb(n,2)//2} arcs.")
print()

# Can we find σ with σ² = complement at n=4?
# That means: σ reverses 3 specific arcs, and σ(σ(T)) reverses ALL 6.
# The 3 arcs reversed by σ and the 3 reversed by σ again must
# be disjoint and cover all 6. So σ reverses a set S of arcs,
# and σ applied to the result reverses the COMPLEMENT of S.
# But σ is a FIXED operation — it always reverses the same arcs.
# So σ(T) reverses S, and σ(σ(T)) reverses S again = back to T.
# That gives σ² = IDENTITY, not complement!

# The issue: a fixed set of arc reversals squares to identity, not complement.
# To get σ² = complement, σ must be CONTEXT-DEPENDENT.

print(f"  PROBLEM: A fixed set of arc reversals squares to identity,")
print(f"  not complement. A 'tournament imaginary unit' cannot be")
print(f"  a simple arc-reversal operation.")
print()
print(f"  BUT: rotation by 90° in the COMPLEX PLANE is not a flip.")
print(f"  It's a ROTATION. The tournament analogue of 90° rotation")
print(f"  is a cyclic relabeling that acts on the VERTEX SET:")
print()

# At n=3: the cyclic permutation (0 1 2) acts on tournaments.
# It relabels vertices: 0→1, 1→2, 2→0.
# This is order 3, not order 4.

# At n=4: the cyclic permutation (0 1 2 3) has order 4!
# σ: relabel 0→1, 1→2, 2→3, 3→0
# σ²: relabel 0→2, 1→3, 2→0, 3→1 (swap pairs)
# σ⁴: identity

# Does σ² = complement? No, σ² is a relabeling, not a reversal.
# But σ² acts on TOURNAMENTS by permuting vertices.

print(f"  At n=4: the cyclic permutation σ = (0 1 2 3) has ORDER 4.")
print(f"  σ⁴ = identity, σ² = (0 2)(1 3) (swap pairs).")
print(f"  This is ORDER 4 — the same order as the imaginary unit i.")
print()

# Let's see what σ does to a specific n=4 tournament
adj_orig = [[0,1,1,1],[0,0,1,0],[0,0,0,1],[0,1,0,0]]  # A specific tournament
h_orig = H_dp(adj_orig, 4)

# Apply σ: relabel 0→1, 1→2, 2→3, 3→0
# New adj[σ(i)][σ(j)] = old adj[i][j]
perm = [1, 2, 3, 0]
adj_sigma = [[0]*4 for _ in range(4)]
for i in range(4):
    for j in range(4):
        adj_sigma[perm[i]][perm[j]] = adj_orig[i][j]
h_sigma = H_dp(adj_sigma, 4)

# Apply σ²
perm2 = [2, 3, 0, 1]
adj_sigma2 = [[0]*4 for _ in range(4)]
for i in range(4):
    for j in range(4):
        adj_sigma2[perm2[i]][perm2[j]] = adj_orig[i][j]
h_sigma2 = H_dp(adj_sigma2, 4)

print(f"  Example at n=4:")
print(f"    T:      H = {h_orig}")
print(f"    σ(T):   H = {h_sigma}")
print(f"    σ²(T):  H = {h_sigma2}")
print(f"    T^c:    H = {H_dp([[1-adj_orig[i][j] if i!=j else 0 for j in range(4)] for i in range(4)], 4)}")
print()

# ==========================================================================
# PART 5: THE REAL MEANING — ROTATION VS REFLECTION
# ==========================================================================

print()
print("PART 5: THE REAL MEANING — ROTATION VS REFLECTION")
print("-" * 60)
print()

print("""  The DEEP answer to "what is an imaginary unit":

  An imaginary unit is a ROTATION, not a REFLECTION.

  REFLECTION (order 2): negation, complement, swap.
    T → T^c  (reverse all arcs)
    (-1)     (negate)
    Applying twice: identity.

  ROTATION (order 4): quarter-turn, imaginary unit.
    i: rotate 90° in a plane
    Applying twice: 180° rotation = negation = -1
    Applying four times: 360° = identity

  The DIFFERENCE:
    Reflections COMMUTE (all reflections form Z/2 or products thereof)
    Rotations DON'T commute in general (SO(n) is non-abelian for n≥3)

  In tournaments:
    COMPLEMENT T → T^c is a REFLECTION (order 2, commutes with everything)
    A 3-CYCLE is a ROTATION of order 3 (NOT order 4!)

  This is the KEY DISCREPANCY: tournament 3-cycles have order 3,
  but imaginary units have order 4.

  HOW THEY RELATE:
  The QUATERNION imaginary unit i has i⁴ = 1 and i² = -1.
  The CYCLIC PERMUTATION (0 1 2) has order 3.

  They're related through the COVERING MAP:
    The quaternion group Q₈ = {±1, ±i, ±j, ±k} has order 8.
    The rotation group of the tetrahedron has order 12.
    The DOUBLE COVER of the rotation group IS the quaternion group:
      SO(3) → Q₈/±1

  The "±" (the sign) is what the imaginary unit adds.
  The UNSIGNED rotation (the 3-cycle) is the tournament.
  The SIGNED rotation (the imaginary unit) is the algebra.

  TOURNAMENTS SEE THE UNSIGNED PART.
  IMAGINARY UNITS CARRY THE SIGN.

  This is why:
    7 = dim(O) - 1: the 7 unsigned rotations (Fano lines)
    8 = dim(O):     the 7 signed rotations + 1 identity
    The forbidden H = 7 is the UNSIGNED count.
    The tournament space at n=3 has 8 = SIGNED count.
""")

# ==========================================================================
# PART 6: THE SIGN AND THE ARC DIRECTION
# ==========================================================================

print()
print("PART 6: THE SIGN IS THE ARC DIRECTION")
print("-" * 60)
print()

print("""  In quaternions: ij = +k, ji = -k.
  The SIGN distinguishes the two orderings of a pair.

  In tournaments: arc (i→j) vs arc (j→i).
  The DIRECTION distinguishes the two orderings of a pair.

  These are THE SAME THING seen differently:
    Sign ±1 ↔ arc direction (forward/backward)
    Product eᵢeⱼ ↔ arc (i,j) in the tournament
    POSITIVE product ↔ arc goes WITH the Fano line's orientation
    NEGATIVE product ↔ arc goes AGAINST the Fano line's orientation

  So a tournament IS an assignment of SIGNS to the imaginary unit products.
  Each arc (i,j) chooses: does eᵢeⱼ = +eₖ or -eₖ?

  A tournament on n=3 with 3 arcs IS an assignment of 3 signs
  to the 3 products ij, jk, ki of a quaternionic triple.
  There are 2³ = 8 sign assignments = 8 tournaments.

  A tournament on n=7 with 21 arcs IS an assignment of signs
  to the products of the 7 octonion imaginary units.
  But there are only 7 independent products (one per Fano line),
  plus 14 = 21 - 7 products involving pairs NOT on the same line.

  THE CRITICAL POINT: only 7 of the 21 products are INDEPENDENT
  (determined by the Fano structure). The other 14 are DERIVED
  from the 7 via the octonion multiplication rule.

  A tournament on 7 vertices makes 21 INDEPENDENT binary choices.
  The octonion multiplication table makes only 7 independent choices.
  The tournament has 21 - 7 = 14 EXTRA degrees of freedom
  beyond the octonionic structure.

  This is WHY the Fano plane embeds in Paley T₇ but doesn't
  determine it: the Fano gives 7 of 21 arc directions;
  the remaining 14 are additional tournament structure.
""")

# Verify: 21 arcs, 7 Fano lines, 14 extra
n = 7
arcs = comb(n, 2)
fano_lines_count = 7
extra = arcs - fano_lines_count
print(f"  At n=7: {arcs} arcs, {fano_lines_count} Fano lines, {extra} extra degrees of freedom")
print(f"  Ratio: {fano_lines_count}/{arcs} = {fano_lines_count/arcs:.3f} = 1/3")
print(f"  The Fano structure determines 1/3 of the tournament!")
print(f"  The remaining 2/3 is non-octonionic structure.")
print()

# ==========================================================================
# PART 7: THE COMPLETE PICTURE — WHAT e² = -1 MEANS FOR TOURNAMENTS
# ==========================================================================

print()
print("=" * 72)
print("  PART 7: WHAT IMAGINARY UNITS ARE, FOR TOURNAMENTS")
print("=" * 72)
print()

print("""
  An imaginary unit is a SIGNED ROTATION.

  It has two components:
    1. A ROTATION (order 3 for quaternions, order 7 for octonions)
       This is the 3-CYCLE in tournament theory.
    2. A SIGN (±1)
       This is the ARC DIRECTION in tournament theory.

  The imaginary unit COMBINES the cycle structure with the sign:
    eᵢeⱼ = ±eₖ means: arcs i→j and j→k compose to give ±(arc i→k)
    The ± encodes whether the composition follows or opposes the cycle.

  A TOURNAMENT is an assignment of signs to ALL pairs.
  An OCTONION is an assignment of signs to the 7 Fano-line pairs.

  The tournament has MORE information than the octonion:
    21 signs (all pairs) vs 7 signs (Fano lines only).
    The extra 14 signs are the non-octonionic content.

  The OCF says: H = I(Ω, 2) counts independent sets of cycles.
  Each cycle IS a signed rotation (an imaginary unit product).
  The independence polynomial counts how many COMPATIBLE signed
  rotations can coexist.

  AT THE DEEPEST LEVEL:

  The imaginary unit is NOT a "number."
  It is a HALF-NEGATION: an operation that, applied twice, reverses.

  In tournament theory:
    The complement T^c is FULL negation (reverse everything).
    The 3-cycle is a THIRD of negation (rotate 120°, not 90°).
    There is NO natural "half-negation" in tournaments
    because 3-cycles have ODD order (3), not EVEN order (4).

  This is WHY tournament theory is TYPE {3,∞} and not {4,∞}:
    The face is a TRIANGLE (order 3), not a SQUARE (order 4).
    The imaginary unit lives in the SQUARE world (order 4).
    Tournaments live in the TRIANGLE world (order 3).

  The Cayley-Dickson tower connects them through the COVERING:
    Triangle rotations (3-cycles) COVER quaternion rotations
    via the binary polyhedral groups (double covers).
    The ± sign that the imaginary unit carries is the KERNEL
    of this covering map.

  The tournament sees the rotation but not the sign.
  The algebra sees both. The sign is the imaginary unit.
  The gap between them is the OCR residual.
""")

print("opus-2026-03-21-S138")
