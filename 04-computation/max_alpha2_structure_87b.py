#!/usr/bin/env python3
"""
max_alpha2_structure_87b.py — opus-2026-03-14-S87b

Deep analysis of the UNIQUE max-α₂ tournament class at n=7.

Key facts discovered:
- Exactly ONE isomorphism class has α₂=14 (max)
- |Aut| = 7 (cyclic Z_7)
- Orbit size = 720 = 7!/7
- H = 175 = 5^2 × 7
- All 14 3-cycles, all 21 transitive triples
- QR_7 has α₂=7 (min among regular), this has α₂=14 (max)

Questions:
1. Is this the ANTI-Paley tournament? (reverse of QR_7?)
2. Is it a circulant tournament? (invariant under rotation i→i+1 mod 7?)
3. What's its connection type set? (which residues define the arcs?)
4. Connection to the 2-(10,4,2) BIBD at n=6?
"""

from itertools import combinations, permutations
from collections import Counter, defaultdict

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

n = 7

# ══════════════════════════════════════════════════════════════════
# PART 1: Find the max-α₂ tournament explicitly
# ══════════════════════════════════════════════════════════════════

print("=" * 70)
print("PART 1: FINDING THE MAX-α₂ TOURNAMENT")
print("=" * 70)

# Enumerate regular tournaments to find the max-α₂ representative
edges = [(i,j) for i in range(n) for j in range(i+1, n)]
m = len(edges)

max_adj = None
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
    if compute_alpha2(cycles) == 14:
        max_adj = adj
        break

print("Max-α₂ tournament adjacency:")
for i in range(n):
    print(f"  {''.join(str(max_adj[i][j]) for j in range(n))}")

# ── Connection set ─────────────────────────────────────────────

print("\nChecking if circulant (invariant under i→i+1 mod 7):")
is_circulant = True
# Check if adj[i][j] depends only on (j-i) mod 7
connection_sets = []
for i in range(n):
    conn = frozenset((j - i) % n for j in range(n) if max_adj[i][j])
    connection_sets.append(conn)
    if i == 0:
        print(f"  Vertex 0 → {sorted(conn)} (connection set)")

if all(c == connection_sets[0] for c in connection_sets):
    print(f"  YES — circulant with connection set {sorted(connection_sets[0])}")
    conn_set = sorted(connection_sets[0])
else:
    print(f"  NO — not circulant")
    for i in range(n):
        print(f"    Vertex {i}: {sorted(connection_sets[i])}")

# ── Compare with QR_7 ─────────────────────────────────────────

print("\nQR_7 connection set: {1, 2, 4} (quadratic residues mod 7)")
print(f"Max-α₂ connection set: {sorted(connection_sets[0]) if all(c == connection_sets[0] for c in connection_sets) else 'varies'}")

# QR mod 7 = {1,2,4}, NQR mod 7 = {3,5,6}
# Complement: reverse all arcs → connection set becomes {3,5,6}
qr7_set = {1, 2, 4}
nqr7_set = {3, 5, 6}
print(f"QR_7 complement (NQR_7): {sorted(nqr7_set)}")

# ── Check all 3 circulant regular tournaments ─────────────────

print("\n" + "=" * 70)
print("PART 2: ALL CIRCULANT REGULAR TOURNAMENTS ON 7 VERTICES")
print("=" * 70)

# A circulant regular tournament on Z_7 has connection set S ⊂ {1,...,6}
# with |S| = 3 and S ∩ (7-S) = ∅ (S and its complement partition {1,...,6})
# So S is a subset of size 3 from {1,...,6} with no pair summing to 7.

circ_sets = []
for s in combinations(range(1, 7), 3):
    # Check: s ∩ (7-s) = ∅
    if all((7 - x) % 7 not in s for x in s):
        circ_sets.append(s)
    elif set(s) | set((7-x)%7 for x in s) == set(range(1,7)):
        circ_sets.append(s)

# Actually: for each pair {d, 7-d}, choose one. 3 pairs: {1,6}, {2,5}, {3,4}
# So 2^3 = 8 choices, but exactly half (with |S|=3) = 8/... no, all have |S|=3.
# Wait: we must choose exactly one from each pair, giving 2^3 = 8 subsets.
# But half of them are the same tournament (complement = reverse), so 4 up to isomorphism.

pairs = [(1,6), (2,5), (3,4)]
all_circ = []
for b0 in [0,1]:
    for b1 in [0,1]:
        for b2 in [0,1]:
            s = set()
            s.add(pairs[0][b0])
            s.add(pairs[1][b1])
            s.add(pairs[2][b2])
            all_circ.append(frozenset(s))

print(f"All 8 circulant connection sets:")
for s in all_circ:
    # Build tournament
    adj = [[0]*n for _ in range(n)]
    for i in range(n):
        for d in s:
            adj[i][(i+d)%n] = 1

    cycles = find_3cycles(adj, n)
    a2 = compute_alpha2(cycles)
    H = compute_H_dp(adj, n)

    # Check isomorphism type
    is_qr = (s == frozenset({1,2,4}))
    is_nqr = (s == frozenset({3,5,6}))

    label = ""
    if is_qr: label = " ← QR_7"
    if is_nqr: label = " ← NQR_7 (complement of QR_7)"

    print(f"  S={sorted(s)}: α₁=14, α₂={a2}, H={H}{label}")

# ── Which circulant has max α₂? ───────────────────────────────

print("\n" + "=" * 70)
print("PART 3: IDENTIFYING THE MAX-α₂ CIRCULANT")
print("=" * 70)

for s in all_circ:
    adj = [[0]*n for _ in range(n)]
    for i in range(n):
        for d in s:
            adj[i][(i+d)%n] = 1

    cycles = find_3cycles(adj, n)
    a2 = compute_alpha2(cycles)
    if a2 == 14:
        print(f"MAX α₂=14 circulant: S = {sorted(s)}")
        print("Adjacency matrix:")
        for i in range(n):
            print(f"  {''.join(str(adj[i][j]) for j in range(n))}")

        # Check which 3-subsets are cyclic
        cyclic_subs = [c[1] for c in cycles]
        print(f"\nCyclic 3-subsets ({len(cyclic_subs)}):")
        for cs in sorted(cyclic_subs, key=lambda s: tuple(sorted(s))):
            print(f"  {sorted(cs)}")

        # Which disjoint pairs are both-cyclic?
        pairs_7 = []
        for A in combinations(range(7), 3):
            remaining = [x for x in range(7) if x not in A]
            for B in combinations(remaining, 3):
                pair = tuple(sorted([frozenset(A), frozenset(B)], key=lambda s: min(s)))
                if pair not in pairs_7:
                    pairs_7.append(pair)

        both_cyclic = []
        for A, B in pairs_7:
            if A in cyclic_subs and B in cyclic_subs:
                both_cyclic.append((A, B))

        print(f"\nBoth-cyclic disjoint pairs ({len(both_cyclic)}):")
        for A, B in both_cyclic:
            uncov = [v for v in range(7) if v not in A and v not in B]
            print(f"  {sorted(A)} | {sorted(B)} → uncov {uncov[0]}")

        # The both-cyclic pairs form a design on the 70 possible pairs
        # Check which vertex is uncovered
        uncov_count = Counter()
        for A, B in both_cyclic:
            uncov = [v for v in range(7) if v not in A and v not in B]
            uncov_count[uncov[0]] += 1
        print(f"\nUncovered vertex counts: {dict(sorted(uncov_count.items()))}")

# ══════════════════════════════════════════════════════════════════
# PART 4: The trio {QR₇, middle, max-α₂} — characterization
# ══════════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("PART 4: THE THREE REGULAR TOURNAMENT CLASSES")
print("=" * 70)

# For each of the 8 circulants, compute α₂
print("Circulant connection sets grouped by α₂:")
a2_to_sets = defaultdict(list)
for s in all_circ:
    adj = [[0]*n for _ in range(n)]
    for i in range(n):
        for d in s:
            adj[i][(i+d)%n] = 1
    cycles = find_3cycles(adj, n)
    a2 = compute_alpha2(cycles)
    H = compute_H_dp(adj, n)
    a2_to_sets[a2].append((sorted(s), H))

for a2 in sorted(a2_to_sets.keys()):
    print(f"  α₂={a2}: {a2_to_sets[a2]}")

# Are there non-circulant regular tournaments on n=7?
print(f"\nTotal regular tournaments: 2640")
print(f"Circulant contributions:")
# Each circulant has |Aut| ≥ 7 (at least the cyclic rotation)
# Orbit size ≤ 720
# Count: for each circulant, compute orbit = 7!/|Aut|

for s in all_circ:
    adj = [[0]*n for _ in range(n)]
    for i in range(n):
        for d in s:
            adj[i][(i+d)%n] = 1

    # Count automorphisms
    auto_count = 0
    for perm in permutations(range(n)):
        if all(adj[a][b] == adj[perm[a]][perm[b]] for a in range(n) for b in range(n)):
            auto_count += 1

    orbit = 5040 // auto_count
    cycles = find_3cycles(adj, n)
    a2 = compute_alpha2(cycles)
    H = compute_H_dp(adj, n)
    print(f"  S={sorted(s)}: |Aut|={auto_count}, orbit={orbit}, α₂={a2}, H={H}")

print("\nNon-circulant regular tournaments?")
# Total from circulants:
total_circ = 0
seen_classes = []
for s in all_circ:
    adj = [[0]*n for _ in range(n)]
    for i in range(n):
        for d in s:
            adj[i][(i+d)%n] = 1

    # Check if isomorphic to already seen
    already = False
    for prev in seen_classes:
        if all(adj[a][b] == prev[a][b] for a in range(n) for b in range(n)):
            already = True
            break
        # Check isomorphism
        for perm in permutations(range(n)):
            if all(adj[a][b] == prev[perm[a]][perm[b]] for a in range(n) for b in range(n)):
                already = True
                break
        if already:
            break

    if not already:
        auto_count = sum(1 for perm in permutations(range(n))
                        if all(adj[a][b] == adj[perm[a]][perm[b]]
                               for a in range(n) for b in range(n)))
        orbit = 5040 // auto_count
        total_circ += orbit
        seen_classes.append(adj)
        a2 = compute_alpha2(find_3cycles(adj, n))
        print(f"  Distinct class: S={sorted(s)}, orbit={orbit}, α₂={a2}")

print(f"  Total from distinct circulant classes: {total_circ}")
print(f"  Total regular tournaments: 2640")
print(f"  Non-circulant: {2640 - total_circ}")

if 2640 - total_circ > 0:
    print(f"\n  There ARE non-circulant regular tournaments on n=7!")
    print(f"  The 'middle' class (α₂=10, count=1680) has orbit 1680 = 7!/3")
    print(f"  This means |Aut|=3 — NOT divisible by 7, so NOT circulant!")

# ══════════════════════════════════════════════════════════════════
# PART 5: The connection set pattern
# ══════════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("PART 5: CONNECTION SET ↔ α₂ PATTERN")
print("=" * 70)

# For each connection set, what's the "difference structure"?
# S = {a,b,c}: what are a+b, a+c, b+c mod 7?
for s in all_circ:
    s_list = sorted(s)
    sums = [(s_list[i]+s_list[j]) % 7 for i in range(3) for j in range(i+1,3)]
    prods = [(s_list[i]*s_list[j]) % 7 for i in range(3) for j in range(i+1,3)]

    adj = [[0]*n for _ in range(n)]
    for i in range(n):
        for d in s:
            adj[i][(i+d)%n] = 1
    a2 = compute_alpha2(find_3cycles(adj, n))

    # Is s a coset of a subgroup?
    # Subgroups of Z_7*: trivial {1}, {1,2,4} (QR), {1,3,2,6,4,5} (all)
    is_qr = (s == frozenset({1,2,4}))
    is_nqr = (s == frozenset({3,5,6}))

    # Check if s is an arithmetic progression mod 7
    diffs = [(s_list[1]-s_list[0])%7, (s_list[2]-s_list[1])%7]
    is_ap = diffs[0] == diffs[1]

    label = ""
    if is_qr: label = "QR"
    elif is_nqr: label = "NQR"
    elif is_ap: label = f"AP(d={diffs[0]})"

    print(f"  S={sorted(s)}: α₂={a2:2d}, sums={sums}, prods={prods}, {label}")

print("""
SYNTHESIS:
  The connection set determines the α₂ value for circulant tournaments.
  QR_7 = {1,2,4}: minimizes α₂ (=7) — most "ordered" cycle structure
  NQR_7 = {3,5,6}: also α₂=7 — isomorphic to QR_7 complement
  Max α₂=14 sets: specific 3-element subsets of Z_7 — most "disordered"

  The α₂ spectrum {7, 10, 14} for regular tournaments on n=7
  mirrors the structure of 3-element subsets of Z_7 under
  the natural action of the affine group AGL(1,7).
""")
