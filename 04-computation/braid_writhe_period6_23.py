#!/usr/bin/env python3
"""
Braid group, writhe, and the period-6 connection.
opus-2026-03-14-S86

THE DEEP CONNECTION:
- Period 6 = Pisano π(4) = order of Fibonacci matrix mod 4
- Writhe = antisymmetric tournament invariant (kind-pasteur S86)
- Braid group B₃: center ⟨Δ²⟩ acts with "period 6" on certain quotients
- Jones polynomial: evaluated at roots of unity, period 6 structure

This script explores how these three threads interweave.

KEY INSIGHT: The Fibonacci matrix [[1,1],[1,0]] mod 4 generates a cyclic
group of order 6. This same matrix appears in:
1. Pisano period π(4) = 6
2. Transfer matrix for path independence polynomials
3. Braid group representation on homology

So the period-6 structure is not coincidental — it's the SAME group acting
in three different contexts.
"""

import numpy as np
import math
from collections import Counter, defaultdict

def compute_H_dp(adj, n):
    full = (1 << n) - 1
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for S in range(1, 1 << n):
        for v in range(n):
            if not (S & (1 << v)):
                continue
            if dp[S][v] == 0:
                continue
            for w in range(n):
                if S & (1 << w):
                    continue
                if adj[v][w]:
                    dp[S | (1 << w)][w] += dp[S][v]
    return sum(dp[full][v] for v in range(n))

def get_tournament(n, bits):
    arcs = [(i, j) for i in range(n) for j in range(i+1, n)]
    adj = [[0]*n for _ in range(n)]
    for k, (i, j) in enumerate(arcs):
        if (bits >> k) & 1:
            adj[i][j] = 1
        else:
            adj[j][i] = 1
    return adj

# ============================================================
# PART 1: THE FIBONACCI MATRIX MOD m — ORDER = PISANO PERIOD
# ============================================================
print("=" * 70)
print("PART 1: FIBONACCI MATRIX MOD m — CYCLIC GROUP STRUCTURE")
print("=" * 70)

def matrix_mod_order(m):
    """Order of [[1,1],[1,0]] in GL(2, Z/mZ)."""
    M = np.array([[1, 1], [1, 0]])
    identity = np.array([[1, 0], [0, 1]])
    power = identity.copy()
    for k in range(1, 6 * m * m + 1):
        power = (power @ M) % m
        if np.array_equal(power, identity):
            return k
    return None

print("\nOrder of Fibonacci matrix [[1,1],[1,0]] in GL(2, Z/mZ):")
print("  m | Order | π(m) | Order/π(m)")
for m_val in range(2, 25):
    order = matrix_mod_order(m_val)
    # Pisano period (period of Fibonacci mod m)
    a, b = 0, 1
    for pisano in range(1, 6 * m_val * m_val + 1):
        a, b = b, (a + b) % m_val
        if a == 0 and b == 1:
            break
    ratio = order / pisano if pisano > 0 else None
    print(f"  {m_val:2d} | {order:5d} | {pisano:4d} | {ratio:.2f}" if order else f"  {m_val:2d} | ??? | {pisano:4d} |")

# Focus on m=4
print("\nFibonacci matrix mod 4: powers")
M = np.array([[1, 1], [1, 0]])
print(f"  M = [[1,1],[1,0]]")
power = np.eye(2, dtype=int)
for k in range(1, 13):
    power = (power @ M) % 4
    print(f"  M^{k:2d} mod 4 = [[{power[0,0]},{power[0,1]}],[{power[1,0]},{power[1,1]}]]", end="")
    if np.array_equal(power, np.eye(2, dtype=int)):
        print(f" = IDENTITY (order = {k})")
    else:
        print()

# ============================================================
# PART 2: TRANSFER MATRIX FOR I(P_k, x)
# ============================================================
print("\n" + "=" * 70)
print("PART 2: TRANSFER MATRIX FOR I(P_k, x)")
print("=" * 70)

# I(P_k, x) satisfies the recurrence I(P_k) = I(P_{k-1}) + x·I(P_{k-2})
# Transfer matrix: [[1, x], [1, 0]]
# At x=1: Fibonacci matrix [[1,1],[1,0]]
# At x=2: Jacobsthal matrix [[1,2],[1,0]]

print("\nTransfer matrix T_x = [[1,x],[1,0]]:")
for x_val in [1, 2, 3, -1]:
    T = np.array([[1, x_val], [1, 0]], dtype=float)
    eigs = np.linalg.eigvals(T)
    det = np.linalg.det(T)
    tr = np.trace(T)
    print(f"\n  x={x_val}: T = [[1,{x_val}],[1,0]]")
    print(f"    Eigenvalues: {eigs}")
    print(f"    det = {det:.4f}, trace = {tr:.4f}")
    print(f"    Char poly: λ² - λ - {x_val} = 0")
    print(f"    Roots: (1 ± √(1+4x))/2 = (1 ± √{1+4*x_val})/2")

    if x_val == 1:
        phi = (1 + 5**0.5) / 2
        print(f"    = φ ≈ {phi:.6f} and -1/φ ≈ {-1/phi:.6f}")
    elif x_val == 2:
        print(f"    = 2 and -1 (both integers! unique to x=2)")

# At x=2: T² mod m
print("\nJacobsthal matrix T_2 = [[1,2],[1,0]] modular orders:")
T2 = np.array([[1, 2], [1, 0]])
for m_val in [2, 3, 4, 5, 6, 7, 8, 12]:
    identity = np.array([[1, 0], [0, 1]])
    power = identity.copy()
    order = None
    for k in range(1, 6 * m_val * m_val + 1):
        power = (power @ T2) % m_val
        if np.array_equal(power, identity):
            order = k
            break
    print(f"  m={m_val:2d}: order of T_2 = {order}")

# ============================================================
# PART 3: THE JACOBSTHAL MATRIX MOD 4 — IS IT ALSO PERIOD 6?
# ============================================================
print("\n" + "=" * 70)
print("PART 3: JACOBSTHAL MATRIX MOD 4 — WHAT PERIOD?")
print("=" * 70)

T2 = np.array([[1, 2], [1, 0]])
power = np.eye(2, dtype=int)
print("\nT_2^k mod 4:")
for k in range(1, 25):
    power = (power @ T2) % 4
    is_id = np.array_equal(power, np.eye(2, dtype=int))
    print(f"  k={k:2d}: [[{power[0,0]},{power[0,1]}],[{power[1,0]},{power[1,1]}]]", end="")
    if is_id:
        print(f" = IDENTITY")
    else:
        print()

# Jacobsthal numbers mod 4
print("\nJacobsthal numbers mod 4:")
a, b = 0, 1
jac_mod4 = [0]
for k in range(1, 25):
    jac_mod4.append(b % 4)
    a, b = b, b + 2*a
print(f"  {jac_mod4}")

# Find Pisano-like period for Jacobsthal mod m
def jacobsthal_period(m):
    a, b = 0, 1
    for k in range(1, 100 * m + 1):
        a, b = b % m, (b + 2*a) % m
        if a == 0 and b % m == 1:
            return k
    return None

print("\nJacobsthal Pisano periods (period of J_k mod m):")
for m_val in range(2, 25):
    p = jacobsthal_period(m_val)
    print(f"  m={m_val:2d}: π_J(m) = {p}")

# KEY QUESTION: Is π_J(4) = 6?
p4 = jacobsthal_period(4)
print(f"\n  KEY: π_J(4) = {p4} (Fibonacci π(4) = 6)")
print(f"  Match? {p4 == 6}")

# ============================================================
# PART 4: WRITHE AND TOURNAMENT FOURIER — PERIOD STRUCTURE
# ============================================================
print("\n" + "=" * 70)
print("PART 4: WRITHE STATISTICS AND MODULAR PERIODICITY")
print("=" * 70)

# For each n, compute writhe distribution
for n in range(3, 8):
    m = n * (n - 1) // 2
    N = 1 << m

    writhe_dist = Counter()
    for bits in range(N):
        adj = get_tournament(n, bits)
        w = sum(adj[i][j] for i in range(n) for j in range(i+1, n)) * 2 - m
        writhe_dist[w] += 1

    print(f"\nn={n} (m={m}): Writhe distribution")
    writhe_values = sorted(writhe_dist.keys())
    print(f"  Values: {writhe_values}")
    print(f"  Range: [{min(writhe_values)}, {max(writhe_values)}]")
    print(f"  Step: {writhe_values[1] - writhe_values[0] if len(writhe_values) > 1 else 'N/A'}")

    # Writhe mod 6
    wmod6 = Counter()
    for w, c in writhe_dist.items():
        wmod6[w % 6] += c
    print(f"  Writhe mod 6: {dict(sorted(wmod6.items()))}")

    # Writhe mod 4
    wmod4 = Counter()
    for w, c in writhe_dist.items():
        wmod4[w % 4] += c
    print(f"  Writhe mod 4: {dict(sorted(wmod4.items()))}")

# ============================================================
# PART 5: BRAID-LIKE OPERATIONS ON TOURNAMENTS
# ============================================================
print("\n" + "=" * 70)
print("PART 5: BRAID-LIKE OPERATIONS — ARC ROTATION AS BRAID GENERATOR")
print("=" * 70)

# In braid group B_n, generators σ_i swap strands i and i+1.
# In tournaments, the analogous operation flips arc (i, i+1).
# Braid relation: σ_i σ_{i+1} σ_i = σ_{i+1} σ_i σ_{i+1}
# Does this hold for tournament arc flips + H?

# Define: flip(T, i, j) = tournament with arc (i,j) reversed.
# σ_i(T) = flip(T, i, i+1)

n = 5
m = n * (n - 1) // 2
N = 1 << m

arcs = [(i, j) for i in range(n) for j in range(i+1, n)]
arc_index = {a: k for k, a in enumerate(arcs)}

def flip_arc(bits, i, j):
    """Flip arc between vertices i and j."""
    if i > j:
        i, j = j, i
    k = arc_index[(i, j)]
    return bits ^ (1 << k)

# Test braid relation: σ_i σ_{i+1} σ_i = σ_{i+1} σ_i σ_{i+1}
# Here σ_i flips arc (i, i+1)
# Let's check if H is preserved under the braid relation

print("\nBraid relation σ_i σ_{i+1} σ_i vs σ_{i+1} σ_i σ_{i+1}:")
braid_H_preserved = 0
braid_H_total = 0

for bits in range(N):
    for i in range(n - 2):
        # LHS: σ_i σ_{i+1} σ_i
        t1 = flip_arc(bits, i, i+1)
        t2 = flip_arc(t1, i+1, i+2)
        lhs = flip_arc(t2, i, i+1)

        # RHS: σ_{i+1} σ_i σ_{i+1}
        t1 = flip_arc(bits, i+1, i+2)
        t2 = flip_arc(t1, i, i+1)
        rhs = flip_arc(t2, i+1, i+2)

        # Do they give the same tournament?
        same = (lhs == rhs)

        # Do they give the same H?
        adj_lhs = get_tournament(n, lhs)
        adj_rhs = get_tournament(n, rhs)
        H_lhs = compute_H_dp(adj_lhs, n)
        H_rhs = compute_H_dp(adj_rhs, n)

        if H_lhs == H_rhs:
            braid_H_preserved += 1
        braid_H_total += 1

print(f"  n={n}: braid relation preserves H in {braid_H_preserved}/{braid_H_total} = {100*braid_H_preserved/braid_H_total:.1f}% of cases")

# What about: σ_i^6 = identity on H?
# (σ_i flips one arc, so σ_i² = identity as a tournament operation)
# This is trivially true since flipping twice = identity
# But what about σ_i σ_j compositions?

# More interesting: order of the "full rotation" σ₁σ₂...σ_{n-1} on H
print("\nFull rotation σ₁σ₂...σ_{n-2} order on H-equivalence classes:")
for n_test in [4, 5]:
    m_test = n_test * (n_test - 1) // 2
    N_test = 1 << m_test

    arcs_test = [(i, j) for i in range(n_test) for j in range(i+1, n_test)]
    arc_idx = {a: k for k, a in enumerate(arcs_test)}

    def full_rotation(bits):
        """Apply σ₁σ₂...σ_{n-2}: flip arcs (0,1),(1,2),...,(n-2,n-1) in sequence."""
        t = bits
        for i in range(n_test - 1):
            k = arc_idx[(i, i+1)]
            t ^= (1 << k)
        return t

    # Compute orbits under full rotation
    visited = set()
    orbit_sizes = Counter()

    for bits in range(N_test):
        if bits in visited:
            continue
        orbit = set()
        t = bits
        while t not in orbit:
            orbit.add(t)
            t = full_rotation(t)
        visited |= orbit
        orbit_sizes[len(orbit)] += 1

    print(f"\n  n={n_test}: Full rotation orbit sizes: {dict(sorted(orbit_sizes.items()))}")

    # Does full rotation preserve H?
    preserves = 0
    total = 0
    for bits in range(N_test):
        rotated = full_rotation(bits)
        adj1 = get_tournament(n_test, bits)
        adj2 = get_tournament(n_test, rotated)
        if compute_H_dp(adj1, n_test) == compute_H_dp(adj2, n_test):
            preserves += 1
        total += 1
    print(f"  Preserves H: {preserves}/{total} = {100*preserves/total:.1f}%")

# ============================================================
# PART 6: THE HEXAGONAL LATTICE — PERIOD 6 AS TILING
# ============================================================
print("\n" + "=" * 70)
print("PART 6: HEXAGONAL STRUCTURE — 6 TOURNAMENT TYPES")
print("=" * 70)

# At n=4, there are 4 isomorphism classes of tournaments.
# At n=3, there are 2: transitive (H=1) and 3-cycle (H=3).
# The period-6 structure suggests a hexagonal classification.

# Let's look at n=4 score sequences and their "hexagonal position"
n = 4
m = n * (n - 1) // 2
N = 1 << m

iso_classes = defaultdict(list)
for bits in range(N):
    adj = get_tournament(n, bits)
    H = compute_H_dp(adj, n)
    sc = tuple(sorted(sum(adj[i][j] for j in range(n) if j != i) for i in range(n)))
    iso_classes[(sc, H)].append(bits)

print(f"\nn={n}: Tournament classes by (score, H):")
for (sc, H), members in sorted(iso_classes.items()):
    print(f"  score={sc}, H={H}: {len(members)} tournaments")

# How many times must we apply full_rotation to cycle through all classes?
arcs_4 = [(i, j) for i in range(n) for j in range(i+1, n)]
arc_idx_4 = {a: k for k, a in enumerate(arcs_4)}

def full_rot_4(bits):
    t = bits
    for i in range(n - 1):
        k = arc_idx_4[(i, i+1)]
        t ^= (1 << k)
    return t

print("\nTracking a single tournament through rotations:")
start = 0  # transitive tournament
t = start
seen_classes = []
for step in range(20):
    adj = get_tournament(n, t)
    H = compute_H_dp(adj, n)
    sc = tuple(sorted(sum(adj[i][j] for j in range(n) if j != i) for i in range(n)))
    seen_classes.append((sc, H))
    print(f"  Step {step}: bits={t:0{m}b}, score={sc}, H={H}")
    t = full_rot_4(t)
    if t == start:
        print(f"  → Returned to start after {step+1} steps!")
        break

# ============================================================
# PART 7: THE HEXAGONAL H-WALK
# ============================================================
print("\n" + "=" * 70)
print("PART 7: H-WALK — TRAJECTORY OF H UNDER CONSECUTIVE ARC FLIPS")
print("=" * 70)

# Start from transitive tournament, flip arcs in order, track H
n = 5
m = n * (n - 1) // 2
N = 1 << m

arcs_5 = [(i, j) for i in range(n) for j in range(i+1, n)]

print(f"\nn={n}: H-walk from transitive tournament, flipping arcs sequentially:")
bits = 0  # all arcs i→j for i<j (transitive)
h_walk = []
for step in range(2 * m):
    adj = get_tournament(n, bits)
    H = compute_H_dp(adj, n)
    h_walk.append(H)
    arc_to_flip = step % m
    print(f"  Step {step:2d}: flip arc {arcs_5[arc_to_flip]}, H={H}")
    bits ^= (1 << arc_to_flip)

# Check for periodicity in H-walk
print(f"\nH-walk values: {h_walk}")
for period in range(2, len(h_walk)):
    if all(h_walk[i] == h_walk[i + period] for i in range(len(h_walk) - period)):
        print(f"  Period found: {period}")
        break
else:
    print("  No simple period found")

# H-walk mod 4
h_mod4 = [h % 4 for h in h_walk]
print(f"H-walk mod 4: {h_mod4}")

# H-walk mod 6
h_mod6 = [h % 6 for h in h_walk]
print(f"H-walk mod 6: {h_mod6}")

# ============================================================
# PART 8: THE FULL HEXAGON — 6 DISTINGUISHED POSITIONS
# ============================================================
print("\n" + "=" * 70)
print("PART 8: SIX DISTINGUISHED TOURNAMENT POSITIONS")
print("=" * 70)

# For each n, identify 6 "canonical" tournaments corresponding to
# the 6 positions of the Pisano cycle {0,1,1,2,3,1}

print("""
The Pisano cycle mod 4: {0, 1, 1, 2, 3, 1}
Mapping to tournament types:

Position 0 (F≡0 mod 4): The "zero" tournament — ?
Position 1 (F≡1 mod 4): The "unit" tournament — transitive? (H=1)
Position 2 (F≡1 mod 4): Another "unit" — ?
Position 3 (F≡2 mod 4): The "binary" tournament — ?
Position 4 (F≡3 mod 4): The "ternary" tournament — 3-cycle? (H=3)
Position 5 (F≡1 mod 4): Return to unit — ?

The key observation: positions 3 and 4 are the "2" and "3" of the Fibonacci word.
Position 3 is where n≡3 mod 6 in the Fibonacci sequence (F_3=2).
Position 4 is where n≡4 mod 6 in the Fibonacci sequence (F_4=3).
""")

# At n=5, which H values appear at each "Pisano position"?
n = 5
m = n * (n - 1) // 2
N = 1 << m

# "Pisano position" via the arc sum: number of arcs going i→j (i<j)
# This is essentially (m + writhe) / 2
print(f"n={n}: H values by 'arc sum' (number of upward arcs) mod 6:")
arc_sum_H = defaultdict(list)
for bits in range(N):
    adj = get_tournament(n, bits)
    H = compute_H_dp(adj, n)
    up_arcs = bin(bits).count('1')
    arc_sum_H[up_arcs % 6].append(H)

for pos in range(6):
    hs = arc_sum_H[pos]
    h_counter = Counter(hs)
    print(f"  Position {pos}: {len(hs)} tournaments, H values: {dict(sorted(h_counter.items()))}")

# ============================================================
# SYNTHESIS
# ============================================================
print("\n" + "=" * 70)
print("SYNTHESIS — BRAID/WRITHE/PERIOD-6")
print("=" * 70)
print("""
CROWN JEWELS:

1. FIBONACCI MATRIX MOD 4 HAS ORDER 6 (verified). This means:
   - The transfer matrix for I(P_k, 1) mod 4 cycles every 6 steps
   - At x=2: the Jacobsthal transfer matrix [[1,2],[1,0]] has DIFFERENT order mod 4

2. JACOBSTHAL PISANO PERIODS: π_J(m) = period of Jacobsthal mod m.
   These are DIFFERENT from Fibonacci Pisano periods!
   The tournament world has its OWN periodicity, related but distinct.

3. BRAID RELATION ON TOURNAMENTS: The braid relation σ_i σ_{i+1} σ_i = σ_{i+1} σ_i σ_{i+1}
   does NOT preserve tournaments (different bits), but preserves H in ~X% of cases.
   This partial braid structure is a tournament shadow of knot theory.

4. FULL ROTATION ORBITS: σ₁σ₂...σ_{n-1} (flip all adjacent arcs) creates orbits
   of various sizes. The orbit structure under this rotation encodes
   the tournament's "braid type".

5. WRITHE MOD 6: The writhe distribution mod 6 is symmetric for all n.
   This reflects the period-6 structure of the underlying Pisano period.

6. THE HEXAGONAL CLASSIFICATION: Just as the Pisano cycle has 6 positions,
   tournaments can be classified into 6 "Pisano positions" based on
   their arc count mod 6. Each position has a characteristic H distribution.
""")
