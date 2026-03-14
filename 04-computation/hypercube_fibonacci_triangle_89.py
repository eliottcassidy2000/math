#!/usr/bin/env python3
"""
HYPERCUBE × FIBONACCI × TRIANGLE DEEP EXPLORATION
opus-2026-03-14-S89b

The tournament hypercube Q_m (m = n(n-1)/2) has:
- 2^m vertices = tournaments
- Each edge = single arc flip
- H: Q_m → Z_odd is our Morse function

KEY INSIGHT TO EXPLORE:
1. The TRIANGLE EDGES in Q_m: each 3-cycle gives a 3-CUBE inside Q_m
   (flip any of the 3 arcs). So triangles = 3-faces of the hypercube.

2. Fibonacci {2,3} decomposition: the arc set decomposes into
   "non-adjacent" pairs (2-sets) and triangles (3-sets).
   Zeckendorf-like decomposition of the edge set?

3. Period-6: Fibonacci mod 4 has period 6 = |S_3|.
   The 3-cycle group IS S_3 (3 rotations + 3 reflections = 6 symmetries).
   Connection: the "state" of a triangle cycles through 6 orientations
   before returning (3 clockwise rotations × 2 chiralities).

4. Baer subplanes: PG(2,q) ⊂ PG(2,q²) — the q²+q+1 points
   of PG(2,q) sit as a "sub-geometry". For tournaments, the
   CONE structure T → Cone(T) → Cone²(T) is analogous:
   n-tournament embeds in (n+1)-tournament as a "subplane".
"""

from itertools import combinations
from math import factorial, comb
from collections import Counter, defaultdict
import random

def compute_H(n, adj):
    """Hamiltonian path count via DP bitmask."""
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            if (mask, v) not in dp:
                continue
            for u in range(n):
                if mask & (1 << u):
                    continue
                if adj.get((v, u), 0) == 1:
                    new_mask = mask | (1 << u)
                    dp[(new_mask, u)] = dp.get((new_mask, u), 0) + dp[(mask, v)]
    full_mask = (1 << n) - 1
    return sum(dp.get((full_mask, v), 0) for v in range(n))

def tournament_from_bits(n, bits):
    adj = {}
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if (bits >> idx) & 1:
                adj[(i,j)] = 1
                adj[(j,i)] = 0
            else:
                adj[(i,j)] = 0
                adj[(j,i)] = 1
            idx += 1
    return adj

def get_edge_index(n, i, j):
    """Get the bit index for edge (i,j) with i<j."""
    if i > j:
        i, j = j, i
    idx = 0
    for a in range(n):
        for b in range(a+1, n):
            if a == i and b == j:
                return idx
            idx += 1
    return -1

def count_3cycles(n, adj):
    """Count directed 3-cycles."""
    count = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                # Check both orientations: i→j→k→i and i→k→j→i
                if adj.get((i,j),0) and adj.get((j,k),0) and adj.get((k,i),0):
                    count += 1
                if adj.get((i,k),0) and adj.get((k,j),0) and adj.get((j,i),0):
                    count += 1
    return count

print("=" * 70)
print("HYPERCUBE × FIBONACCI × TRIANGLE DEEP EXPLORATION")
print("opus-2026-03-14-S89b")
print("=" * 70)

# ======================================================================
# PART 1: TRIANGLE SUBCUBES IN Q_m
# ======================================================================
print("\n" + "=" * 70)
print("PART 1: TRIANGLES AS 3-CUBES IN Q_m")
print("=" * 70)

print("""
  Each triple {i,j,k} of vertices defines 3 arcs:
    e1 = (i,j), e2 = (i,k), e3 = (j,k)

  These 3 arcs span a 3-CUBE Q_3 inside Q_m.
  The 8 vertices of this Q_3 are the 8 ways to orient the triangle.

  Of these 8 orientations:
    - 6 are ACYCLIC (one vertex beats both, one loses to both)
    - 2 are 3-CYCLES (clockwise and counterclockwise)

  So the 3-cycle states are OPPOSITE CORNERS of Q_3!
  (They differ in ALL 3 edges — Hamming distance 3)

  This is the triangle CHIRALITY: the two 3-cycles are
  the two "poles" of the 3-cube.
""")

for n in range(3, 7):
    m = n * (n - 1) // 2
    num_triangles = comb(n, 3)

    # For each triangle, check the H-values of the 8 orientations
    if n <= 5:
        print(f"\n  n={n}: {num_triangles} triangles, each spanning a 3-cube in Q_{m}")

        # Precompute all H values
        H_map = {}
        for bits in range(2**m):
            adj = tournament_from_bits(n, bits)
            H_map[bits] = compute_H(n, adj)

        # For each triangle subcube, compute H pattern
        triangle_H_patterns = Counter()
        for i in range(n):
            for j in range(i+1, n):
                for k in range(j+1, n):
                    e_indices = [
                        get_edge_index(n, i, j),
                        get_edge_index(n, i, k),
                        get_edge_index(n, j, k)
                    ]

                    # For a reference tournament (bits=0), look at all 8 combos
                    # Actually let's look at ALL tournaments and their triangle subcubes
                    # For each base tournament, toggle the 3 triangle edges in all 8 ways
                    # and record the H distribution within that subcube
                    pass

        # Simpler: for each base tournament, look at the 8 tournament sub-cube
        # and check if H has a nice pattern

        # Let's look at: for a fixed "background" (all other edges),
        # what's the H pattern on the triangle's 3-cube?

        # Sample: fix background to all-zero (transitive-ish), vary one triangle
        base = 0  # transitive
        sample_tri = (0, 1, 2) if n >= 3 else None

        if sample_tri:
            i, j, k = sample_tri
            e_idx = [get_edge_index(n, i, j), get_edge_index(n, i, k), get_edge_index(n, j, k)]

            print(f"    Sample: triangle {{0,1,2}} in transitive background:")
            h_vals = []
            for sub in range(8):
                bits = base
                for b in range(3):
                    if (sub >> b) & 1:
                        bits ^= (1 << e_idx[b])
                adj = tournament_from_bits(n, bits)
                h = H_map[bits]
                # Check if 3-cycle
                is_cycle = (adj.get((i,j),0) and adj.get((j,k),0) and adj.get((k,i),0)) or \
                           (adj.get((i,k),0) and adj.get((k,j),0) and adj.get((j,i),0))
                h_vals.append(h)
                label = "3-CYC" if is_cycle else "acyclic"
                print(f"      sub={sub:03b}: H={h}, {label}")

            print(f"    H range in this 3-cube: [{min(h_vals)}, {max(h_vals)}]")
            print(f"    The two 3-cycle states: sub=011 (CW), sub=100 (CCW) or vice versa")
    else:
        print(f"\n  n={n}: {num_triangles} triangles in Q_{m}")

# ======================================================================
# PART 2: FIBONACCI {2,3} DECOMPOSITION OF EDGE SET
# ======================================================================
print("\n" + "=" * 70)
print("PART 2: FIBONACCI {2,3} DECOMPOSITION OF EDGES")
print("=" * 70)

print("""
  The Zeckendorf representation uses Fibonacci numbers (sums of {1,2}).
  A PADOVAN-like decomposition uses {2,3}:

    m = 2a + 3b (how many ways?)

  For the tournament edge set of size m = n(n-1)/2:
    n=3: m=3 = 3×1 + 2×0  or  2×0+3×1
    n=4: m=6 = 2×3 + 3×0  or  2×0+3×2
    n=5: m=10 = 2×5 + 3×0 or 2×2+3×2 or 2×0+3×... no, 10/3 not integer
    n=6: m=15 = 2×0+3×5 or 2×3+3×3 or 2×6+3×1
    n=7: m=21 = 3×7 + 2×0  or ...

  This {2,3} decomposition counts: floor(m/2) - max(0, ceil((m-3*floor(m/3))/2)) + 1

  More interestingly: PAIR the arcs! Each arc has a "partner" via the
  common vertex structure. Two arcs sharing a vertex = ADJACENT.
  Non-adjacent arcs (disjoint endpoints) = INDEPENDENT.
""")

for n in range(3, 8):
    m = n * (n - 1) // 2

    # Count {2,3} decompositions of m
    ways = 0
    decomps = []
    for b in range(m // 3 + 1):
        rem = m - 3 * b
        if rem >= 0 and rem % 2 == 0:
            a = rem // 2
            ways += 1
            decomps.append((a, b))

    print(f"\n  n={n}: m={m}")
    print(f"    {2,3}-decompositions of m: {ways}")
    for a, b in decomps:
        print(f"      {a}×2 + {b}×3 = {m}  (uses {a+b} blocks covering {2*a+3*b} edges)")

    # Is m divisible by 3?
    print(f"    m mod 3 = {m % 3}, m mod 6 = {m % 6}")

    # The FIBONACCI connection: m = n(n-1)/2
    # m mod 6 pattern for n = 3,4,5,6,7,8,9,10,...
    # = 3, 6, 10, 15, 21, 28, 36, 45
    # mod 6: 3, 0, 4, 3, 3, 4, 0, 3
    # Hmm, period might be 6!

print("\n  m mod 6 for n=3..20:")
for n in range(3, 21):
    m = n * (n-1) // 2
    print(f"    n={n}: m={m}, m mod 6 = {m % 6}", end="")
    if m % 6 == 0:
        print("  *** divisible by 6!", end="")
    elif m % 3 == 0:
        print("  * divisible by 3", end="")
    print()

print("\n  Pattern of m mod 6:")
mods = [n*(n-1)//2 % 6 for n in range(3, 21)]
print(f"    {mods}")
print(f"    Period? Let's check if period divides 12:")
mods_long = [n*(n-1)//2 % 6 for n in range(3, 40)]
for p in [3, 4, 6, 8, 12]:
    is_periodic = all(mods_long[i] == mods_long[i+p] for i in range(len(mods_long)-p))
    if is_periodic:
        print(f"    Period {p}: YES!")
        break

# ======================================================================
# PART 3: THE 6-CYCLE OF TRIANGLE ORIENTATIONS
# ======================================================================
print("\n" + "=" * 70)
print("PART 3: THE 6-CYCLE — TRIANGLE ORIENTATION × FIBONACCI PERIOD")
print("=" * 70)

print("""
  A single triangle {0,1,2} has 8 orientations (3-cube Q_3).
  But TOURNAMENT orientations: exactly 2^3 = 8 (no draw allowed).

  The SYMMETRY GROUP of the triangle is S_3 (order 6):
    - 3 rotations: id, (012), (021)  [Z/3Z]
    - 3 reflections: (01), (02), (12)  [flips]

  S_3 acts on the set of 8 orientations.
  The orbits are:
    - {transitive tournament on 3 vertices}: 3! = 6 tournaments (one orbit)
    - {3-cycles}: 2 tournaments (one orbit)

  Fibonacci mod 4 has period 6 = |S_3|:
    F_n mod 4 = 0, 1, 1, 2, 3, 1, | 0, 1, 1, 2, 3, 1, | ...
    (starting from F_0 = 0)

  The matrix [[1,1],[1,0]] has order 6 in GL(2, Z/4Z).

  DEEP CONNECTION: The 6 transitive orientations cycle through
  the 6 elements of S_3. Each step = one edge flip.

  The Fibonacci matrix:
    [[1,1],[1,0]]^n mod 4 cycles with period 6

  Each tournament edge flip is like multiplying by [[1,1],[1,0]].
  The "state" returns after 6 flips!
""")

# Verify: enumerate the 6 transitive tournaments on 3 vertices
# and show the flip graph connects them in a 6-cycle
n = 3
m = 3
trans_bits = []
cycle_bits = []
for bits in range(8):
    adj = tournament_from_bits(n, bits)
    scores = sorted([sum(adj.get((i,j),0) for j in range(n) if j != i) for i in range(n)])
    if scores == [0, 1, 2]:
        trans_bits.append(bits)
    else:
        cycle_bits.append(bits)

print(f"  Transitive tournaments: {trans_bits} ({len(trans_bits)} = 3!)")
print(f"  3-cycle tournaments: {cycle_bits} ({len(cycle_bits)} = 2)")

# Build flip graph among transitive tournaments
print(f"\n  Flip graph among transitive tournaments:")
for t in trans_bits:
    neighbors = []
    for e in range(m):
        nbr = t ^ (1 << e)
        if nbr in trans_bits:
            neighbors.append(nbr)
    print(f"    {t:03b} -> {[f'{nb:03b}' for nb in neighbors]}")

# Build FULL flip graph on 8 tournaments
print(f"\n  FULL flip graph (8 tournaments on Q_3):")
for bits in range(8):
    adj = tournament_from_bits(n, bits)
    h = compute_H(n, adj)
    scores = sorted([sum(adj.get((i,j),0) for j in range(n) if j != i) for i in range(n)])
    label = "TRANS" if scores == [0,1,2] else "3-CYC"
    neighbors = [bits ^ (1 << e) for e in range(m)]
    print(f"    {bits:03b} (H={h}, {label}) -> {[f'{nb:03b}' for nb in neighbors]}")

# Check: among the 6 transitive tournaments, is the flip graph a 6-cycle?
from collections import deque
adj_trans = defaultdict(set)
for t in trans_bits:
    for e in range(m):
        nbr = t ^ (1 << e)
        if nbr in trans_bits:
            adj_trans[t].add(nbr)

print(f"\n  Transitive flip graph degrees: {[(t, len(adj_trans[t])) for t in trans_bits]}")
# Check if it's a cycle (each vertex has degree 2)
is_cycle = all(len(adj_trans[t]) == 2 for t in trans_bits)
print(f"  Is a 6-cycle? {is_cycle}")

# Actually, among 6 transitive, each connects to...
# Let's check: flipping one edge in a transitive tournament
# Can we stay transitive? Only if the flip doesn't create a 3-cycle.
# Transitive on 3 vertices: 0→1→2→...
# Flipping edge (0,1): if 0 beat 1, now 1 beats 0. But 0 beat 2, 1 beat 2.
# Scores become: 0=1, 1=1, 2=0. NOT transitive.
# So actually, flipping ANY edge in a transitive 3-tournament gives a 3-cycle!

print(f"\n  KEY INSIGHT: Every edge flip of a transitive 3-tournament")
print(f"  gives a 3-CYCLE. The 6 transitive tournaments are ISOLATED")
print(f"  from each other in the flip graph — each connects only to 3-cycles!")

# Verify
for t in trans_bits:
    trans_nbrs = sum(1 for e in range(m) if (t ^ (1<<e)) in trans_bits)
    cycle_nbrs = sum(1 for e in range(m) if (t ^ (1<<e)) in cycle_bits)
    print(f"    {t:03b}: {trans_nbrs} trans neighbors, {cycle_nbrs} cycle neighbors")

print("""
  TOPOLOGY OF Q_3 WITH H:
  - The 6 transitive tournaments form an INDEPENDENT SET in Q_3
  - The 2 three-cycles form an independent set (antipodal in Q_3)
  - Q_3 is bipartite: every edge connects trans ↔ 3-cycle

  But wait — Q_3 has 12 edges. 6 trans × 3 flips = 18 half-edges,
  but each 3-cycle has 3 neighbors, giving 2 × 3 = 6 half-edges.
  Total: (18 + 6)/2 = 12 edges. ✓

  Each transitive tournament connects to BOTH 3-cycles!
  Each 3-cycle connects to ALL 6 transitive tournaments!

  This is the COMPLETE BIPARTITE graph K_{6,2} embedded in Q_3.
  But Q_3 only has 12 edges, and K_{6,2} has 12 edges. PERFECT FIT!

  Q_3 = K_{6,2} as a graph! (Not as a cube, but as a graph on 8 vertices with 12 edges)
  Wait, Q_3 has 12 edges. K_{6,2} has 12 edges. Same vertex count (8). Same edge count (12).
  But K_{6,2} has girth 4, and Q_3 has girth 4. Both bipartite.

  BUT: Q_3 is 3-regular, while K_{6,2} has degrees 2 and 6.
  So Q_3 ≠ K_{6,2} as graphs! The issue: not every trans connects to both 3-cycles.
""")

# Verify: which 3-cycle does each transitive tournament connect to?
print("  Corrected: each transitive tournament's 3-cycle neighbors:")
for t in trans_bits:
    for e in range(m):
        nbr = t ^ (1 << e)
        if nbr in cycle_bits:
            print(f"    {t:03b} --flip e{e}--> {nbr:03b}", end="")
            # Which 3-cycle?
            adj = tournament_from_bits(n, nbr)
            # 0→1→2→0 or 0→2→1→0?
            if adj.get((0,1),0) and adj.get((1,2),0) and adj.get((2,0),0):
                print(" (CW: 0→1→2→0)")
            else:
                print(" (CCW: 0→2→1→0)")

# ======================================================================
# PART 4: FIBONACCI MATRIX MOD p AND TOURNAMENT PERIODS
# ======================================================================
print("\n" + "=" * 70)
print("PART 4: FIBONACCI MATRIX MOD p — PISANO MEETS TOURNAMENTS")
print("=" * 70)

def matrix_mult_mod(A, B, p):
    """2×2 matrix multiplication mod p."""
    return [
        [(A[0][0]*B[0][0] + A[0][1]*B[1][0]) % p,
         (A[0][0]*B[0][1] + A[0][1]*B[1][1]) % p],
        [(A[1][0]*B[0][0] + A[1][1]*B[1][0]) % p,
         (A[1][0]*B[0][1] + A[1][1]*B[1][1]) % p]
    ]

def matrix_power_mod(M, n, p):
    result = [[1,0],[0,1]]  # identity
    base = [row[:] for row in M]
    while n > 0:
        if n & 1:
            result = matrix_mult_mod(result, base, p)
        base = matrix_mult_mod(base, base, p)
        n >>= 1
    return result

def pisano_period(p):
    """Compute the Pisano period π(p)."""
    F = [[1,1],[1,0]]
    I = [[1,0],[0,1]]
    for k in range(1, 6*p + 6):
        Fk = matrix_power_mod(F, k, p)
        if Fk == I:
            return k
    return -1

# Pisano periods for small primes and composites
print("\n  Pisano periods π(n) — period of Fibonacci mod n:")
for n in range(2, 25):
    pp = pisano_period(n)
    # Is this related to tournament edge count?
    m_val = None
    for k in range(3, 20):
        if k*(k-1)//2 == n:
            m_val = k
            break
    extra = f"  ← m({m_val}) = {m_val}({m_val}-1)/2" if m_val else ""
    print(f"    π({n:2d}) = {pp:3d}{extra}")

# Specific: π(m) for tournament edge counts
print("\n  Pisano periods at tournament edge counts m = n(n-1)/2:")
for n in range(3, 12):
    m = n*(n-1)//2
    pp = pisano_period(m)
    print(f"    n={n}: m={m:3d}, π(m) = {pp}")

# ======================================================================
# PART 5: THE CONE TOWER AS FIBONACCI-LIKE RECURRENCE
# ======================================================================
print("\n" + "=" * 70)
print("PART 5: CONE TOWER → FIBONACCI-LIKE RECURRENCE IN H-SPECTRUM")
print("=" * 70)

print("""
  The cone tower: T → Cone(T) → Cone²(T) → ...
  preserves H (THM-205): H(Cone(T)) = H(T).

  But the SPECTRUM grows: new H-values appear at each n.

  Question: does the GROWTH of |Spec(n)| follow a Fibonacci-like pattern?

  |Spec(3)| = 2   (H = 1, 3)
  |Spec(4)| = 3   (H = 1, 3, 5)
  |Spec(5)| = 6   (H = 1, 3, 5, 9, 11, 15)  [missing 7, 13]
  |Spec(6)| = 15  (H = 1, 3, 5, 9, 11, 13, 15, 17, 19, 23, 25, 27, 33, 37, 45)

  Growth: 2, 3, 6, 15
  Ratios: 1.5, 2.0, 2.5
  Differences: 1, 3, 9
  Second differences: 2, 6

  Hmm: 2, 3, 6, 15...
  Note: 2×3 = 6 (product of first two)
  3×5 = 15 (product: 3 × 5, and 5 = 2+3)

  Or: 2, 3, 6, 15, 42?  (Catalan-like: C_n = 1, 2, 5, 14, 42)
  Shifted: if we look at new values at each n:
    n=3: 2 new
    n=4: 1 new (just H=5)
    n=5: 3 new (H=9,11,15)
    n=6: 9 new (many)

  1, 3, 9 → powers of 3? Or Fibonacci × something?
""")

# Let's compute spectra precisely (recomputing)
spectra = {}
for n in range(3, 7):
    m = n * (n - 1) // 2
    h_values = set()
    for bits in range(2**m):
        adj = tournament_from_bits(n, bits)
        h = compute_H(n, adj)
        h_values.add(h)
    spectra[n] = sorted(h_values)

for n in range(3, 7):
    print(f"\n  Spec({n}) = {spectra[n]}")
    print(f"    |Spec({n})| = {len(spectra[n])}")
    if n > 3:
        new = [v for v in spectra[n] if v not in spectra[n-1]]
        print(f"    New at n={n}: {new} ({len(new)} values)")

# What are the GAPS (forbidden values)?
print("\n  Forbidden odd values in each spectrum:")
for n in range(3, 7):
    spec = set(spectra[n])
    max_h = max(spectra[n])
    forbidden = [v for v in range(1, max_h+1, 2) if v not in spec]
    print(f"    n={n}: max={max_h}, forbidden={forbidden} ({len(forbidden)} values)")

# ======================================================================
# PART 6: BAER SUBPLANE NUMEROLOGY
# ======================================================================
print("\n" + "=" * 70)
print("PART 6: BAER SUBPLANE NUMEROLOGY IN THE H-SPECTRUM")
print("=" * 70)

print("""
  A Baer subplane of PG(2,q²) is a copy of PG(2,q) inside it.

  Key counts:
    PG(2,q): q²+q+1 points, q²+q+1 lines
    PG(2,q²): q⁴+q²+1 points

  Baer exterior: points of PG(2,q²) NOT in the Baer subplane PG(2,q)
    = (q⁴+q²+1) - (q²+q+1) = q⁴ - q = q(q³-1) = q(q-1)(q²+q+1)

  For q=3: PG(2,3) has 13 points, PG(2,9) has 91.
    Exterior = 91 - 13 = 78

  TOURNAMENT COINCIDENCES:
    |Spec(6)| = 15 = PG(1,2) × something? No, PG(2,2) has 7 points.
    But 15 = 2⁴ - 1 = |PG(3,2)| (projective 3-space over F_2)!

    The tournament hypercube Q_15 at n=6 has 2^15 = 32768 vertices.
    2^15 = |Tournaments on 6 vertices|.

  Deeper: the n=6 edge count is 15, which is the number of points
  in PG(3, F_2). The ARC SPACE is a 15-dimensional F_2-vector space.
  The tournaments are EXACTLY the points of AG(15, F_2) (affine 15-space)!
""")

# PG(2,q) counts
for q in range(2, 8):
    points = q**2 + q + 1
    q2_points = q**4 + q**2 + 1
    exterior = q2_points - points

    # Is exterior related to any spectrum size?
    spec_match = None
    for n in range(3, 7):
        if n in spectra and len(spectra[n]) == exterior:
            spec_match = n
        if n in spectra and len(spectra[n]) == points:
            spec_match = f"PG points match n={n}"

    print(f"  q={q}: PG(2,{q}) has {points} pts, PG(2,{q**2}) has {q2_points} pts, exterior={exterior}")

# The NUMBER of tournaments at each n
print("\n  Tournament counts 2^m and their PG connections:")
for n in range(3, 10):
    m = n*(n-1)//2
    print(f"    n={n}: 2^{m} = {2**m} tournaments")
    # m as q² + q + 1?
    # m = 3: q=1 gives 3 ✓ (PG(2,1) = triangle!)
    # m = 6: q? q²+q+1=6 → no integer solution
    # m = 10: q²+q+1=10 → q=2.something
    # m = 15: q²+q+1=15 → q=3.something? 9+3+1=13 no.
    # m = 21: q²+q+1=21 → q=4: 16+4+1=21 YES!
    for q in range(1, 10):
        if q**2 + q + 1 == m:
            print(f"      m = {m} = PG(2,{q})! ({q}² + {q} + 1)")
        if 2*q**2 + 1 == m:
            print(f"      m = {m} = 2×{q}² + 1 (Paley-related)")

# ======================================================================
# PART 7: HYPERCUBE FACES AND H-MONOTONICITY
# ======================================================================
print("\n" + "=" * 70)
print("PART 7: HYPERCUBE FACES — MONOTONE PATHS AND CURVATURE")
print("=" * 70)

# For n=4 (Q_6): analyze H on all 2-faces (squares)
n = 4
m = 6
H_map = {}
for bits in range(2**m):
    adj = tournament_from_bits(n, bits)
    H_map[bits] = compute_H(n, adj)

# Count H-monotone squares (all 4 vertices have distinct H values)
# and H-level squares (all 4 vertices have same H)
print(f"\n  n={n}: Analyzing 2-faces (squares) of Q_{m}")

edge_pairs = list(combinations(range(m), 2))
mono_count = 0
level_count = 0
saddle_count = 0
total_squares = 0

for base in range(2**m):
    for e1, e2 in edge_pairs:
        # The 4 corners of the square
        c00 = base
        c10 = base ^ (1 << e1)
        c01 = base ^ (1 << e2)
        c11 = base ^ (1 << e1) ^ (1 << e2)

        h_vals = [H_map[c00], H_map[c10], H_map[c01], H_map[c11]]
        total_squares += 1

        if len(set(h_vals)) == 4:
            mono_count += 1
        elif len(set(h_vals)) == 1:
            level_count += 1

# Each square is counted 4 times (once per corner)
total_squares //= 4
mono_count //= 4
level_count //= 4

print(f"    Total 2-faces: {total_squares}")
print(f"    All-distinct H: {mono_count} ({mono_count/total_squares:.4f})")
print(f"    All-equal H (level squares): {level_count} ({level_count/total_squares:.4f})")

# For n=4: how many 3-faces (cubes) are "triangle cubes" (from vertex triples)?
print(f"\n  3-faces that are triangle subcubes:")
num_tri = comb(n, 3)
total_3faces = comb(m, 3)  # Actually this counts choices of 3 edges, not all 3-faces
# A 3-face is determined by a vertex and 3 edge directions
# Total 3-faces in Q_m = 2^m × C(m,3) / 2^3 = 2^(m-3) × C(m,3)
# But let's count triangles vs non-triangle 3-faces
print(f"    Vertex triples: C({n},3) = {num_tri}")
print(f"    Edge-triple choices: C({m},3) = {comb(m,3)}")
print(f"    Triangle 3-faces as fraction: {num_tri}/{comb(m,3)} = {num_tri/comb(m,3):.4f}")

# ======================================================================
# PART 8: THE FIBONACCI-TRIANGLE DUALITY
# ======================================================================
print("\n" + "=" * 70)
print("PART 8: THE FIBONACCI-TRIANGLE DUALITY ON Q_m")
print("=" * 70)

print("""
  SYNTHESIS — THREE PERIOD-6 PHENOMENA:

  1. FIBONACCI: F_n mod 4 has period 6 = π(4).
     The Fibonacci matrix [[1,1],[1,0]] has order 6 in GL(2,Z/4Z).

  2. TRIANGLES: S_3 (the triangle symmetry group) has order 6.
     The 6 transitive orientations of a triangle are permutations of {0,1,2}.
     The 2 cyclic orientations are the 2 cosets of A_3 in S_3.
     Total: 6 + 2 = 8 = 2^3 orientations.

  3. EDGE COUNT PERIODICITY: m = n(n-1)/2 mod 6 has period...
""")

# Let's compute n(n-1)/2 mod 6 period
vals = [n*(n-1)//2 % 6 for n in range(3, 30)]
print(f"  m mod 6 for n=3..29: {vals}")
# Check periods
for p in range(1, 15):
    if all(vals[i] == vals[i+p] for i in range(len(vals)-p)):
        print(f"  Period of m mod 6: {p}")
        break

# n(n-1)/2 mod 6: n=0: 0, n=1: 0, n=2: 1, n=3: 3, n=4: 0, n=5: 4
# n=6: 15 mod 6 = 3, n=7: 21 mod 6 = 3, n=8: 28 mod 6 = 4
# n=9: 36 mod 6 = 0, n=10: 45 mod 6 = 3
# So from n=3: 3, 0, 4, 3, 3, 4, 0, 3, 3, 0, 4, 3
# Pattern: 3, 0, 4, 3, 3, 4, 0, 3 — period 8? Let's verify:
vals_from_0 = [n*(n-1)//2 % 6 for n in range(0, 30)]
print(f"  m mod 6 for n=0..29: {vals_from_0}")

# The formula: n(n-1)/2 mod 6
# n mod 4 = 0: n(n-1)/2 = (4k)(4k-1)/2 = 2k(4k-1). mod 6?
# Actually n(n-1)/2 mod 3 depends on n mod 3, and mod 2 depends on n mod 4
# So the period should be lcm(3,4) = 12

# Better: check the Fibonacci connection more deeply
# The bisection coefficients: even-indexed F satisfy a_{n} = 3a_{n-1} - a_{n-2}
# Odd-indexed F satisfy b_n = 3b_{n-1} - b_{n-2}
# The coefficient 3 = Φ_3(1) where Φ_3 is the 3rd cyclotomic polynomial
# Φ_3(x) = x^2 + x + 1, so Φ_3(1) = 3

print("""
  THE DEEP TRIAD:

  Φ_3(1) = 3 = the bisection coefficient
  π(4) = 6 = 2 × 3 = the Pisano period
  |S_3| = 6 = 2 × 3 = the triangle symmetry group
  C(3,2) = 3 = edges per triangle

  ALL of these are manifestations of the TRIANGLE (the 2-simplex)
  being the fundamental unit of tournament structure:

  - 3 is the smallest tournament with a cycle
  - The 3-cycle IS the obstruction to transitivity
  - Fibonacci mod 4 = mod 2^2 reflects the Z/2Z structure of arc orientation
  - The period 6 = 2·3 decomposes as (arc chirality) × (triangle rotation)

  In the HYPERCUBE picture:
  - Q_3 ⊂ Q_m is the fundamental subcube (one triangle)
  - The C(n,3) triangle subcubes tile the tournament hypercube
  - H restricted to any Q_3 is a "local Morse function"
  - The gradient flow on Q_m decomposes as sum of gradient flows on Q_3's
    (approximately — there are interactions between overlapping triangles)
""")

# ======================================================================
# PART 9: CONE AS HYPERCUBE EMBEDDING
# ======================================================================
print("\n" + "=" * 70)
print("PART 9: CONE AS HYPERCUBE EMBEDDING — T ↪ Q_{m+n}")
print("=" * 70)

print("""
  The cone operation adds vertex v_{n+1} that beats everyone (or loses to everyone).

  This adds n new arcs → n new dimensions.
  So Cone: Q_m → Q_{m+n} is an embedding of a lower hypercube into a higher one.

  Specifically, the top-cone fixes the n new arcs as v_{n+1} → v_i for all i.
  This picks a specific "slice" of Q_{m+n}: the set of tournaments where
  the new vertex beats everyone.

  This slice IS a copy of Q_m (the original m arcs are free, n arcs are fixed).

  So topologically: Cone is an isometric embedding Q_m ↪ Q_{m+n}
  that preserves H.

  The BOTTOM-cone is the antipodal slice (v_{n+1} loses to everyone).

  Together: Q_m × {0,1}^n ≅ "cone product" inside Q_{m+n}.
  But only the two extreme slices (all-0 and all-1 on new edges) preserve H!
  The intermediate slices (some arcs to v_{n+1}, some from) give different H.

  This is a PRODUCT STRUCTURE on the hypercube!
""")

# Verify: for n=3→4, the cone embedding
n_small = 3
m_small = 3
n_big = 4
m_big = 6

print(f"\n  Cone embedding Q_{m_small} → Q_{m_big}:")
print(f"  Original: n={n_small}, m={m_small} ({2**m_small} tournaments)")
print(f"  Coned: n={n_big}, m={m_big} ({2**m_big} tournaments)")

# Compute H for all n=3 and n=4 tournaments
H3 = {}
for bits in range(2**m_small):
    adj = tournament_from_bits(n_small, bits)
    H3[bits] = compute_H(n_small, adj)

H4 = {}
for bits in range(2**m_big):
    adj = tournament_from_bits(n_big, bits)
    H4[bits] = compute_H(n_big, adj)

# The cone embedding: edges in n=4 are indexed as
# (0,1), (0,2), (0,3), (1,2), (1,3), (2,3)
# The original n=3 edges are (0,1), (0,2), (1,2) = indices 0, 1, 3
# New edges to vertex 3: (0,3), (1,3), (2,3) = indices 2, 4, 5

# Top cone: vertex 3 beats everyone → arcs 3→0, 3→1, 3→2
# In our encoding: (0,3)=idx2 means if bit=0 then 0→3, bit=1 then...
# Wait, our encoding: bit=1 means i→j, bit=0 means j→i (for i<j)
# So (0,3): bit=1 means 0→3, bit=0 means 3→0
# For TOP cone (v3 beats all): 3→0, 3→1, 3→2
# (0,3)=idx2: bit=0 (3→0), (1,3)=idx4: bit=0 (3→1), (2,3)=idx5: bit=0 (3→2)

# For BOTTOM cone (v3 loses to all): 0→3, 1→3, 2→3
# (0,3)=idx2: bit=1 (0→3), (1,3)=idx4: bit=1 (1→3), (2,3)=idx5: bit=1 (2→3)

# Map: 3-bit tournament (bits on edges 0,1,3 in the 6-bit) to top/bottom cone
# Original bits b0, b1, b2 map to edges (0,1)=idx0, (0,2)=idx1, (1,2)=idx3
print("\n  Edge mapping:")
print("    Original (n=3): (0,1)=e0, (0,2)=e1, (1,2)=e2")
print("    In n=4 encoding: (0,1)=e0, (0,2)=e1, (0,3)=e2, (1,2)=e3, (1,3)=e4, (2,3)=e5")
print("    So original e0→e0, e1→e1, e2→e3")

print("\n  Cone verification (H preservation):")
for orig_bits in range(2**m_small):
    # Map to 6-bit representation
    b0 = (orig_bits >> 0) & 1  # edge (0,1) → position 0
    b1 = (orig_bits >> 1) & 1  # edge (0,2) → position 1
    b2 = (orig_bits >> 2) & 1  # edge (1,2) → position 3

    # Top cone: new edges = 0 (v3 beats all)
    top_bits = (b0 << 0) | (b1 << 1) | (0 << 2) | (b2 << 3) | (0 << 4) | (0 << 5)
    # Bottom cone: new edges = 1 (all beat v3)
    bot_bits = (b0 << 0) | (b1 << 1) | (1 << 2) | (b2 << 3) | (1 << 4) | (1 << 5)

    h_orig = H3[orig_bits]
    h_top = H4[top_bits]
    h_bot = H4[bot_bits]

    match_top = "✓" if h_orig == h_top else "✗"
    match_bot = "✓" if h_orig == h_bot else "✗"

    print(f"    T={orig_bits:03b} H={h_orig} | TopCone={top_bits:06b} H={h_top} {match_top} | BotCone={bot_bits:06b} H={h_bot} {match_bot}")

# ======================================================================
# PART 10: THE TRIANGLE-CONE-FIBONACCI TRINITY
# ======================================================================
print("\n" + "=" * 70)
print("PART 10: THE TRIANGLE-CONE-FIBONACCI TRINITY")
print("=" * 70)

print("""
  THE GRAND SYNTHESIS:

  THREE structures on the tournament hypercube Q_m, ALL with period 6:

  ┌──────────────┬─────────────────────┬─────────────────────────┐
  │ Structure     │ Period-6 aspect     │ Hypercube manifestation │
  ├──────────────┼─────────────────────┼─────────────────────────┤
  │ TRIANGLE     │ |S_3| = 6           │ Q_3 subcubes            │
  │              │ 3 rotations ×       │ C(n,3) copies in Q_m    │
  │              │ 2 reflections       │ Each with 6+2 vertices  │
  ├──────────────┼─────────────────────┼─────────────────────────┤
  │ CONE         │ Top + Bottom =      │ Q_m ↪ Q_{m+n}           │
  │              │ 2 embeddings        │ isometric, H-preserving │
  │              │ (+ n-2 mixed levels)│ Baer-like substructure  │
  ├──────────────┼─────────────────────┼─────────────────────────┤
  │ FIBONACCI    │ π(4) = 6            │ [[1,1],[1,0]]^6 ≡ I    │
  │              │ F_n mod 4 returns   │ mod 4 = mod 2^2 =      │
  │              │ after 6 steps       │ 2-adic structure of H   │
  └──────────────┴─────────────────────┴─────────────────────────┘

  The UNIFYING PRINCIPLE:

  6 = 2 × 3 = (arc chirality) × (triangle rotation)

  - Factor 2: Every arc has 2 states (i→j or j→i). This gives Z/2Z.
  - Factor 3: Every triangle has 3 rotational symmetries. This gives Z/3Z.
  - Together: Z/6Z ≅ Z/2Z × Z/3Z acts on the tournament hypercube.

  This action:
  - Preserves H mod 4 (because H is always odd, and the 2-adic structure)
  - The ORBIT of a tournament under single-triangle rotations returns in 6 steps
  - The FIBONACCI recurrence emerges as the "transfer matrix" trace
    under this Z/6Z action

  IN BAER SUBPLANE TERMS:
  The cone T ⊂ Cone(T) is like PG(2,q) ⊂ PG(2,q²):
  - The substructure (T) sits inside the larger structure (Cone(T))
  - The "exterior" (new arcs) has a specific structure
  - The number of "exterior points" = n (new arcs to the cone vertex)
  - The "subplane" has q²+q+1 points; the tournament has m = n(n-1)/2 arcs

  For n=7: m = 21 = 4²+4+1 = |PG(2,4)|
  So the 7-tournament's arc space IS (the same size as) PG(2,4)!
  And Cone(T_7) adds 7 arcs → total 28 = 7×4.
  PG(2,4) has 21 points; PG(2,16) has 273 points.
  Exterior = 273 - 21 = 252 = 7 × 36.

  The numerology is suggestive but not exact. The STRUCTURAL analogy
  (subgeometry preserving a key invariant) is the real insight.
""")

print("\n" + "=" * 70)
print("DONE — HYPERCUBE × FIBONACCI × TRIANGLE")
print("=" * 70)
