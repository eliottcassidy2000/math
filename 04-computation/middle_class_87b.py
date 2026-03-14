#!/usr/bin/env python3
"""
middle_class_87b.py — opus-2026-03-14-S87b

The "middle class" of regular tournaments on n=7.

Three classes:
  QR₇:    α₂=7,  H=189, |Aut|=21 (circulant, doubly-regular)
  Middle:  α₂=10, H=171, |Aut|=3  (non-circulant!) ← THIS ONE
  AP₇:    α₂=14, H=175, |Aut|=7  (circulant, max disjoint)

The middle class is the ONLY non-circulant regular tournament on n=7.
It has |Aut|=3, orbit=1680=7!/3.

Questions:
1. What does the middle class look like?
2. Is its automorphism Z₃?
3. Why is α₂=10 (not 7 or 14)?
4. H=171 = 9 × 19 — interesting factorization
5. Is it related to the Steiner system or some other design?
"""

from itertools import combinations, permutations
from collections import Counter, defaultdict

def find_3cycles(adj, n):
    cycles = []
    for a, b, c in combinations(range(n), 3):
        if adj[a][b] and adj[b][c] and adj[c][a]:
            cycles.append(((a,b,c), frozenset({a,b,c})))
        elif adj[a][c] and adj[c][b] and adj[b][a]:
            cycles.append(((a,c,b), frozenset({a,b,c})))
    return cycles

def compute_alpha2(cycles):
    nc = len(cycles)
    if nc < 2:
        return 0
    count = 0
    for i in range(nc):
        for j in range(i+1, nc):
            if not (cycles[i][1] & cycles[j][1]):
                count += 1
    return count

def compute_H_dp(adj, n):
    full = (1 << n) - 1
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for S in range(1, 1 << n):
        for v in range(n):
            if not (S & (1 << v)): continue
            if dp[S][v] == 0: continue
            for w in range(n):
                if S & (1 << w): continue
                if adj[v][w]:
                    dp[S | (1 << w)][w] += dp[S][v]
    return sum(dp[full][v] for v in range(n))

n = 7
edges = [(i,j) for i in range(n) for j in range(i+1, n)]
m = len(edges)

# ══════════════════════════════════════════════════════════════════
# PART 1: Find a middle-class representative
# ══════════════════════════════════════════════════════════════════

print("=" * 70)
print("PART 1: FINDING A MIDDLE-CLASS REPRESENTATIVE")
print("=" * 70)

middle_adj = None
for bits in range(1 << m):
    adj = [[0]*n for _ in range(n)]
    for k, (i,j) in enumerate(edges):
        if bits & (1 << k):
            adj[i][j] = 1
        else:
            adj[j][i] = 1

    if any(sum(adj[i]) != 3 for i in range(n)):
        continue

    cycles = find_3cycles(adj, n)
    a2 = compute_alpha2(cycles)
    if a2 == 10:
        middle_adj = adj
        break

print("Middle-class tournament adjacency:")
for i in range(n):
    print(f"  {''.join(str(middle_adj[i][j]) for j in range(n))}")

# ── Automorphism group ─────────────────────────────────────────

print("\nAutomorphism group:")
auts = []
for perm in permutations(range(n)):
    if all(middle_adj[a][b] == middle_adj[perm[a]][perm[b]]
           for a in range(n) for b in range(n)):
        auts.append(perm)

print(f"  |Aut| = {len(auts)}")
print(f"  Automorphisms:")
for perm in auts:
    cycle_str = perm_to_cycles(perm) if False else str(perm)
    print(f"    {perm}")

# Express as cycle notation
def perm_cycles(perm):
    n = len(perm)
    visited = [False] * n
    cycles = []
    for i in range(n):
        if visited[i]:
            continue
        cycle = []
        j = i
        while not visited[j]:
            visited[j] = True
            cycle.append(j)
            j = perm[j]
        if len(cycle) > 1:
            cycles.append(tuple(cycle))
    return cycles if cycles else [()]

for perm in auts:
    cycles = perm_cycles(perm)
    if not cycles or cycles == [()]:
        cycle_str = "id"
    else:
        cycle_str = " ".join(f"({' '.join(map(str, c))})" for c in cycles)
    print(f"    {cycle_str}")

# Is it Z_3?
if len(auts) == 3:
    # Check if the non-identity automorphisms have order 3
    for perm in auts:
        if perm == tuple(range(n)):
            continue
        # Apply perm 3 times
        p = list(perm)
        p2 = [p[p[i]] for i in range(n)]
        p3 = [p[p2[i]] for i in range(n)]
        is_id = all(p3[i] == i for i in range(n))
        print(f"    {perm} has order {'3' if is_id else '?'}")
    print(f"  Aut ≅ Z₃")

# ── Score sequence and local structure ────────────────────────

print("\nScore sequence (sorted):", sorted([sum(middle_adj[i]) for i in range(n)], reverse=True))

# For each vertex, its "outneighborhood"
print("\nOut-neighborhoods:")
for v in range(n):
    out = [j for j in range(n) if middle_adj[v][j]]
    inn = [j for j in range(n) if middle_adj[j][v]]
    print(f"  v={v}: out={out}, in={inn}")

# ── Cycle structure ────────────────────────────────────────────

cycles = find_3cycles(middle_adj, n)
print(f"\n3-cycles ({len(cycles)}):")
cycle_sets = [c[1] for c in cycles]
for c in cycles:
    print(f"  {sorted(c[1])}")

# Disjoint pairs
pairs_7 = []
for A in combinations(range(7), 3):
    remaining = [x for x in range(7) if x not in A]
    for B in combinations(remaining, 3):
        pair = tuple(sorted([frozenset(A), frozenset(B)], key=lambda s: min(s)))
        if pair not in [(p[0], p[1]) for p in pairs_7]:
            pairs_7.append(pair)

both_cyclic = []
for A, B in pairs_7:
    if A in cycle_sets and B in cycle_sets:
        both_cyclic.append((A, B))

print(f"\nBoth-cyclic pairs (α₂={len(both_cyclic)}):")
uncov_count = Counter()
for A, B in both_cyclic:
    uncov = [v for v in range(7) if v not in A and v not in B]
    uncov_count[uncov[0]] += 1
    print(f"  {sorted(A)} | {sorted(B)} → uncov {uncov[0]}")

print(f"\nUncovered vertex distribution: {dict(sorted(uncov_count.items()))}")
# Is it uniform?
if len(set(uncov_count.values())) == 1:
    print("  UNIFORM (all vertices equally likely uncovered)")
else:
    print(f"  NOT uniform — some vertices more covered than others")

# ══════════════════════════════════════════════════════════════════
# PART 2: The orbit structure under Z₃
# ══════════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("PART 2: Z₃ ORBIT STRUCTURE")
print("=" * 70)

# The Z₃ automorphism fixes some vertices and permutes others
# Find the generator
gen = None
for perm in auts:
    if perm != tuple(range(n)):
        gen = perm
        break

if gen:
    cycles = perm_cycles(gen)
    print(f"Generator: {gen}")
    print(f"Cycle decomposition: {cycles}")

    # Fixed points
    fixed = [i for i in range(n) if gen[i] == i]
    print(f"Fixed points: {fixed}")

    # 3-cycles
    three_cycles = [c for c in cycles if len(c) == 3]
    print(f"3-cycles: {three_cycles}")

    # The vertex set splits into: fixed points + orbits of size 3
    # With |Aut|=3: must have exactly one fixed point and two 3-orbits
    # (1 + 2×3 = 7)
    print(f"\nVertex partition under Z₃:")
    print(f"  Fixed: {fixed}")
    for c in three_cycles:
        print(f"  Orbit: {c}")

    # The fixed vertex is special — it's the "axis" of the Z₃ symmetry
    if fixed:
        v0 = fixed[0]
        out = [j for j in range(n) if middle_adj[v0][j]]
        inn = [j for j in range(n) if middle_adj[j][v0]]
        print(f"\n  Fixed vertex {v0}:")
        print(f"    Out-neighbors: {out}")
        print(f"    In-neighbors:  {inn}")

        # How does the fixed vertex relate to the two orbits?
        for c in three_cycles:
            out_in_orbit = [v for v in c if v in out]
            in_in_orbit = [v for v in c if v in inn]
            print(f"    Orbit {c}: {len(out_in_orbit)} out, {len(in_in_orbit)} in")

# ══════════════════════════════════════════════════════════════════
# PART 3: Compare all three classes
# ══════════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("PART 3: THE THREE CLASSES COMPARED")
print("=" * 70)

# QR₇
qr7 = {1, 2, 4}
adj_qr = [[0]*n for _ in range(n)]
for i in range(n):
    for d in qr7:
        adj_qr[i][(i+d)%n] = 1

# AP₇
ap7 = {1, 2, 3}
adj_ap = [[0]*n for _ in range(n)]
for i in range(n):
    for d in ap7:
        adj_ap[i][(i+d)%n] = 1

classes = [
    ("QR₇", adj_qr),
    ("Middle", middle_adj),
    ("AP₇", adj_ap),
]

for name, adj in classes:
    cycles = find_3cycles(adj, n)
    cycle_sets = [c[1] for c in cycles]
    a2 = compute_alpha2(cycles)
    H = compute_H_dp(adj, n)

    # Compute doubly-regular parameter
    dr_counts = Counter()
    for i in range(n):
        for j in range(n):
            if i == j: continue
            if adj[i][j]:
                common = sum(1 for k in range(n) if k != i and k != j
                            and adj[i][k] and adj[k][j])
                dr_counts[common] += 1

    # Check vertex-3-cycle incidence
    vertex_cycle_count = [sum(1 for cs in cycle_sets if v in cs) for v in range(n)]

    print(f"\n  {name}:")
    print(f"    α₂={a2}, H={H}")
    print(f"    Doubly-regular parameter: {dict(dr_counts)}")
    print(f"    Vertex-cycle counts: {vertex_cycle_count}")
    print(f"    H = {H} = ", end="")
    # Factor H
    h = H
    factors = []
    for p in [2,3,5,7,11,13,17,19,23]:
        while h % p == 0:
            factors.append(p)
            h //= p
    if h > 1: factors.append(h)
    print(" × ".join(map(str, factors)))

# ══════════════════════════════════════════════════════════════════
# PART 4: What makes the middle class special?
# ══════════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("PART 4: WHAT MAKES THE MIDDLE CLASS SPECIAL?")
print("=" * 70)

print("""
THE THREE REGULAR TOURNAMENT CLASSES ON n=7:

  QR₇ (Multiplicative):
    - Connection set {1,2,4} = quadratic residues
    - α₂ = 7 (minimum), H = 189 = 27 × 7 = 3³ × 7
    - |Aut| = 21 = 3 × 7 (LARGEST symmetry)
    - Doubly-regular with λ=1
    - Vertex-cycle incidence: uniform (each vertex in 6 cycles)
    - Cycles spread maximally → fewest disjoint pairs

  Middle (Non-circulant):
    - NOT invariant under any rotation
    - α₂ = 10 (middle), H = 171 = 9 × 19
    - |Aut| = 3 (Z₃ — SMALLEST symmetry)
    - NOT doubly-regular
    - 1 fixed point + 2 orbits of 3 under Z₃
    - BREAKS the circulant pattern

  AP₇ (Additive):
    - Connection set {1,2,3} = consecutive residues
    - α₂ = 14 (maximum), H = 175 = 25 × 7 = 5² × 7
    - |Aut| = 7 (Z₇ — pure rotation)
    - Vertex-cycle incidence: uniform (each vertex in 6 cycles)
    - Cycles clustered → most disjoint pairs

THE PATTERN:
  |Aut|: 21 > 7 > 3
  α₂:     7 < 10 < 14
  H:     189 > 175 > 171

  More symmetry → fewer disjoint cycles → higher H
  (Symmetry "forces" cycles to overlap, reducing independence)

  H factorizations:
    189 = 3³ × 7  (pure Mersenne structure)
    175 = 5² × 7  (Fibonacci × Mersenne)
    171 = 3² × 19 (p=3 × p=19, Jacobsthal primes)

  All three H values share factor 7... wait, 171/7 = 24.43, not integer.
  171 = 9 × 19. 189 = 27 × 7. 175 = 25 × 7.
  Only QR and AP share the factor 7!
  Middle has H = 171 = 9 × 19, with NO factor of 7.

  The FORBIDDEN value 7 divides H for QR and AP but NOT for Middle.
  The middle class AVOIDS the forbidden number!
""")
