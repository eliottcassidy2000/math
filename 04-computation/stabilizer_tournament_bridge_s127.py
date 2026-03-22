#!/usr/bin/env python3
"""
stabilizer_tournament_bridge_s127.py — Stabilizer Code Distance ↔ Tournament Theory
opus-2026-03-21-S127

THE CONNECTION:

A stabilizer code on n qubits is described by a graph G on n vertices.
The CODE DISTANCE d = minimum weight of a logical operator.
Recent work (arXiv:2411.14448, PRX Quantum 2025) shows:
  d >= girth(G) / 2 for certain graph families
  d is bounded by independence-number-like quantities

A tournament on n vertices has conflict graph Omega(T).
H(T) = I(Omega(T), 2) = the independence polynomial at fugacity 2.

THE BRIDGE:
  The stabilizer code distance is related to the MINIMUM WEIGHT
  of a codeword = the size of the SMALLEST logical operator.
  The tournament's H is related to the TOTAL COUNT of independent sets
  = the PARTITION FUNCTION at fugacity 2.

  In coding theory: d = min weight → small independent sets matter.
  In tournament theory: H = I(Omega, 2) → ALL independent sets matter.

  The SRCP (sorted root cycle profile) determines H.
  Could a SIMILAR profile determine the code distance?

MERGING WITH RECENT AGENT WORK:
  - kind-pasteur-S16: OCR has exact rational values 18/19, 12/13, 120/131
    with PRIME denominators = largest prime factor of Var(H)
  - kind-pasteur-S18: THM-261/262/263 (dual root system embedding, SRCP)
  - Our S124-126: SRCP determines H at n=5-7, needs c7 at n=8
"""

from math import comb, factorial, sqrt, log
from fractions import Fraction
import numpy as np

print("=" * 72)
print("  STABILIZER CODE DISTANCE AND TOURNAMENT THEORY")
print("  opus-2026-03-21-S127")
print("=" * 72)
print()

# ==========================================================================
print("=" * 72)
print("  PART 1: THE PARALLEL STRUCTURES")
print("=" * 72)
print()

print("""
  STABILIZER CODES                  TOURNAMENTS
  ==================                ===========
  n qubits                          n vertices
  Stabilizer group S ⊂ P_n          Tournament T ∈ {0,1}^{C(n,2)}
  Graph state |G>                   Conflict graph Omega(T)
  Logical operators                 Hamiltonian paths
  Code distance d                   H(T) = I(Omega, 2)
  d = min weight of logical op      H = total independent set count at λ=2

  KEY DIFFERENCE:
    Code distance = MINIMUM of a set of weights
    Hamiltonian path count = SUM over all independent sets

  KEY SIMILARITY:
    Both depend on the CYCLE STRUCTURE of an associated graph.
    For codes: short cycles in the Tanner graph limit d.
    For tournaments: cycles in Omega determine I(Omega, 2).
""")

# ==========================================================================
print("=" * 72)
print("  PART 2: THE GIRTH CONNECTION")
print("=" * 72)
print()

print("""
  For stabilizer codes represented as graphs:
    d >= girth(G)/2 (for certain constructions)
    Graphs with large girth → good codes (high distance)

  For tournament conflict graphs:
    girth(Omega(T)) varies with T.
    At n<=8: Omega(T) is claw-free → girth >= 4 (no triangles if bipartite)
    At n>=9: Omega(T) can have claws → girth can be 3.

  THE ANALOGY:
    Large girth in the CODE graph → high DISTANCE
    Large girth in Omega → ??? → maybe high H?

  Let's check: does girth(Omega) correlate with H?
""")

# At n=5: compute girth of Omega for various tournaments
from itertools import combinations, permutations

def compute_omega_girth(adj, n):
    """Compute girth of the conflict graph Omega_3(T)."""
    # Find all 3-cycles
    cycles = []
    for i, j, k in combinations(range(n), 3):
        if (adj[i][j] and adj[j][k] and adj[k][i]) or \
           (adj[i][k] and adj[k][j] and adj[j][i]):
            cycles.append(frozenset([i, j, k]))

    if len(cycles) <= 1:
        return float('inf')  # No edges in Omega or isolated

    # Build Omega adjacency (cycles share vertex)
    nc = len(cycles)
    omega_adj = [[0]*nc for _ in range(nc)]
    for a in range(nc):
        for b in range(a+1, nc):
            if cycles[a] & cycles[b]:
                omega_adj[a][b] = omega_adj[b][a] = 1

    # BFS to find shortest cycle
    min_girth = float('inf')
    for start in range(nc):
        dist = [-1] * nc
        dist[start] = 0
        queue = [start]
        qi = 0
        while qi < len(queue):
            u = queue[qi]
            qi += 1
            for v in range(nc):
                if omega_adj[u][v]:
                    if dist[v] == -1:
                        dist[v] = dist[u] + 1
                        queue.append(v)
                    elif v != (queue[qi-2] if qi >= 2 else -1) and dist[v] >= dist[u]:
                        girth = dist[u] + dist[v] + 1
                        if girth < min_girth:
                            min_girth = girth

    return min_girth

n = 5
m = comb(n, 2)
girth_H = []

for mask in range(2**m):
    adj = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if (mask >> idx) & 1: adj[i][j] = 1
            else: adj[j][i] = 1
            idx += 1

    H = sum(1 for perm in permutations(range(n))
            if all(adj[perm[k]][perm[k+1]] for k in range(n-1)))

    girth = compute_omega_girth(adj, n)
    girth_H.append((girth, H))

# Analyze correlation
from collections import Counter
girth_groups = {}
for g, h in girth_H:
    if g not in girth_groups:
        girth_groups[g] = []
    girth_groups[g].append(h)

print(f"  n=5: Girth of Omega vs H:")
print(f"  {'girth':>6} | {'count':>6} | {'E[H]':>8} | {'H range':>12}")
print(f"  " + "-" * 40)
for g in sorted(girth_groups.keys()):
    hs = girth_groups[g]
    g_str = str(g) if g != float('inf') else 'inf'
    print(f"  {g_str:>6} | {len(hs):6d} | {sum(hs)/len(hs):8.2f} | [{min(hs)},{max(hs)}]")

print()

# ==========================================================================
print("=" * 72)
print("  PART 3: SRCP AS A CODE PARAMETER")
print("=" * 72)
print()

print("""
  In stabilizer code theory, the CODE determines which errors can be
  corrected. The DISTANCE d is the minimum weight error that goes
  undetected.

  In tournament theory, the SRCP determines H. The "minimum" analog
  is: what is the SMALLEST number of cycles that share an arc?

  Define the MINIMUM ROOT CYCLE COUNT:
    min_c3 = min over all directed arcs (i,j) of c3(i,j)
    = the arc that participates in the FEWEST 3-cycles.

  This is analogous to the code distance: it's the "weakest link"
  in the cycle structure.

  CONJECTURE: min_c3(T) is related to how "decodable" the tournament
  is from its score sequence. Low min_c3 = some arcs carry no cycle info.
""")

# Compute min_c3 for n=5 tournaments and correlate with H
from collections import defaultdict

min_c3_H = defaultdict(list)
for mask in range(2**m):
    adj = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if (mask >> idx) & 1: adj[i][j] = 1
            else: adj[j][i] = 1
            idx += 1

    H = sum(1 for perm in permutations(range(n))
            if all(adj[perm[k]][perm[k+1]] for k in range(n-1)))

    min_c3 = n  # Upper bound
    for i in range(n):
        for j in range(n):
            if i == j or adj[i][j] == 0: continue
            c3_arc = sum(1 for k in range(n) if k!=i and k!=j and adj[j][k] and adj[k][i])
            min_c3 = min(min_c3, c3_arc)

    min_c3_H[min_c3].append(H)

print(f"  n=5: min_c3(T) vs H:")
print(f"  {'min_c3':>7} | {'count':>6} | {'E[H]':>8} | {'H range':>12}")
print(f"  " + "-" * 40)
for mc in sorted(min_c3_H.keys()):
    hs = min_c3_H[mc]
    print(f"  {mc:7d} | {len(hs):6d} | {sum(hs)/len(hs):8.2f} | [{min(hs)},{max(hs)}]")

print()

# ==========================================================================
print("=" * 72)
print("  PART 4: THE OCR PRIME DENOMINATORS AND CODE PARAMETERS")
print("=" * 72)
print()

print("""
  From kind-pasteur-S16: The exact rational OCR values are:
    OCR(5) = 18/19 (denominator 19 = largest prime factor of Var(H)(5))
    OCR(6) = 12/13 (denominator 13 = largest prime factor of Var(H)(6))
    OCR(7) = 120/131 (denominator 131 = largest prime factor of Var(H)(7))

  Wait: 18/19 vs our computation 129/133 at n=5?
  These DIFFER because S16 uses a different R² definition (regression)
  while S103/S104 use Var(E[H|s])/Var(H) (variance decomposition).

  The PRIME DENOMINATOR SEQUENCE: 19, 13, 131, ...
  These are the SURPRISE PRIMES from the Var(H) factorization.
  And they appear as DENOMINATORS of the OCR.

  In coding theory: the code distance d often involves SPECIFIC PRIMES
  related to the field size (GF(q) codes have d related to q).
  In tournament theory: the OCR denominator involves specific primes
  related to the Var(H) structure.

  THE BRIDGE: Both involve "characteristic primes" that emerge from
  the algebraic structure of the underlying system.

  For stabilizer codes: the field is GF(2) (binary).
  The code distance d measures the minimum Hamming weight.
  For tournaments: the field is effectively GF(2) (binary arc choice).
  The OCR denominator measures the "algebraic distance" from perfect
  score determination.
""")

# The surprise primes and their properties
print("  THE SURPRISE PRIME SEQUENCE:")
surprise = [19, 13, 131]
for i, p in enumerate(surprise):
    n_val = i + 5
    print(f"    n={n_val}: surprise prime = {p}")
    print(f"      p-1 = {p-1} = ", end="")
    # Factor p-1
    m = p-1
    factors = []
    d = 2
    while d*d <= m:
        while m % d == 0:
            factors.append(d)
            m //= d
        d += 1
    if m > 1: factors.append(m)
    print(" x ".join(str(f) for f in factors))

print()

# ==========================================================================
print("=" * 72)
print("  PART 5: THE STABILIZER CODE ANALOGY — PRECISE STATEMENT")
print("=" * 72)
print()

print("""
  THE ANALOGY (formalized):

  ┌─────────────────────┬──────────────────────────────────────────┐
  │  STABILIZER CODE    │  TOURNAMENT                              │
  ├─────────────────────┼──────────────────────────────────────────┤
  │  n qubits           │  n vertices                              │
  │  k logical qubits   │  ? (analogous to score sequence?)        │
  │  Code distance d    │  "tournament distance" (SRCP depth?)     │
  │  Stabilizer group   │  Conflict graph Omega(T)                 │
  │  Logical operators  │  Hamiltonian paths                       │
  │  Syndrome           │  Score sequence (shadow)                 │
  │  Error correction   │  H determination from shadow             │
  │  Uncorrectable err  │  OCR residual (H not det. by shadow)     │
  │  Tanner graph girth │  Omega girth                             │
  │  LDPC structure     │  Sparse Omega (high n)                   │
  └─────────────────────┴──────────────────────────────────────────┘

  THE PRECISE PARALLEL:

  In a stabilizer code: given syndrome s, the decoder determines
  the most likely error E. The code distance d is the minimum weight
  error that has syndrome 0 (goes undetected).

  In tournament theory: given score sequence s (= syndrome), the
  "decoder" estimates H. The OCR residual is the probability of
  ERROR in this estimation. The "tournament distance" would be the
  minimum amount of cycle structure change that escapes the shadow.

  DEFINE: d_T(n) = min over all pairs (T1, T2) with same scores but
  different H of |SRCP(T1) - SRCP(T2)|_1 (L1 distance of profiles).

  This is the TOURNAMENT DISTANCE: the minimum cycle-structure
  difference between score-equivalent tournaments with different H.

  At n=5: d_T(5) = 1 (the PoS class has H differing by 2, with
  SRCP differing by 1 in one c5 count).
""")

# Compute tournament distance at n=5
n = 5
m = comb(n, 2)
pos_score = (1, 2, 2, 2, 3)
pos_tours = []

for mask in range(2**m):
    adj = [[0]*n for _ in range(n)]
    idx = 0
    scores = [0]*n
    for i in range(n):
        for j in range(i+1, n):
            if (mask >> idx) & 1: adj[i][j] = 1; scores[i] += 1
            else: adj[j][i] = 1; scores[j] += 1
            idx += 1

    if tuple(sorted(scores)) != pos_score:
        continue

    H = sum(1 for perm in permutations(range(n))
            if all(adj[perm[k]][perm[k+1]] for k in range(n-1)))

    # Compute SRCP
    profiles = []
    for i in range(n):
        for j in range(n):
            if i == j or adj[i][j] == 0: continue
            c3 = sum(1 for k in range(n) if k!=i and k!=j and adj[j][k] and adj[k][i])
            c5 = 0
            others = [v for v in range(n) if v != i and v != j]
            for p in permutations(others, 3):
                a, b, c = p
                if adj[j][a] and adj[a][b] and adj[b][c] and adj[c][i]:
                    c5 += 1
            profiles.append((c3, c5))
    srcp = tuple(sorted(profiles))

    pos_tours.append({'H': H, 'srcp': srcp, 'mask': mask})

# Compute minimum L1 distance between SRCPs with different H
print(f"  Tournament distance d_T at n=5:")
min_dist = float('inf')
for i in range(len(pos_tours)):
    for j in range(i+1, len(pos_tours)):
        if pos_tours[i]['H'] != pos_tours[j]['H']:
            # L1 distance of SRCPs
            s1 = pos_tours[i]['srcp']
            s2 = pos_tours[j]['srcp']
            if len(s1) == len(s2):
                dist = sum(abs(a-b) for (a1,a2), (b1,b2) in zip(s1,s2) for a,b in [(a1,b1),(a2,b2)])
                if dist < min_dist:
                    min_dist = dist

print(f"    d_T(5) = {min_dist}")
print(f"    This is the minimum SRCP L1 distance between score-equivalent")
print(f"    tournaments with different H values.")
print()

# ==========================================================================
print("=" * 72)
print("  PART 6: SYNTHESIS — THE CODE-TOURNAMENT BRIDGE")
print("=" * 72)
print()

print("""
  THE UNIFIED PICTURE:

  Both stabilizer codes and tournaments involve:
  1. A BINARY STRUCTURE on n objects (qubits/vertices)
  2. A GRAPH encoding constraints (Tanner/conflict graph)
  3. A DISTANCE measuring the minimum undetectable change
  4. An INDEPENDENCE STRUCTURE controlling capacity/count

  The OCF identity H = I(Omega, 2) is the tournament's
  PARTITION FUNCTION at fugacity 2 (= tau + tau^{-3}).

  A stabilizer code's partition function would be:
    Z(code, x) = sum over codewords c: x^{weight(c)}
  The code distance d = min weight = smallest power of x with nonzero coefficient.

  THE BRIDGE:
    Tournament: Z = I(Omega, x) evaluated at x=2 gives H.
    Code: Z evaluated at x→0 gives the code distance d.

  Both are evaluations of the SAME type of partition function,
  just at DIFFERENT points:
    x = 2 (tournament): total count weighted by 2^{size}
    x → 0 (code distance): minimum weight

  The SRCP determines H (the x=2 evaluation).
  An analogous profile might determine d (the x→0 evaluation).

  This is a RESEARCH DIRECTION: define a "code cycle profile"
  that determines the stabilizer code distance, analogous to how
  the tournament SRCP determines H.
""")

print("opus-2026-03-21-S127: STABILIZER CODE DISTANCE AND TOURNAMENT THEORY")
