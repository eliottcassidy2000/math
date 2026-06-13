#!/usr/bin/env python3
"""
paley_max_alpha2_87b.py — opus-2026-03-14-S87b

The Paley tournament QR_7 and max-α₂ tournaments at n=7.

Key findings to investigate:
1. ALL max-α₂=14 tournaments at n=7 are regular (scores all 3)
2. 14 cyclic + 21 transitive = 35 = C(7,3), and 21 = H_forb_2!
3. Is QR_7 among the max-α₂ tournaments?
4. How many non-isomorphic regular tournaments on 7 vertices?
5. The number 14 = 2×7... is this 2n for n=7?

Also: connection to the 21 = H_forb_2 mystery.
The max-α₂ tournament has exactly 21 transitive triples.
Tournament theory says: a tournament is LOCALLY ALMOST REGULAR
iff it has the minimum number of transitive triples.
"""

from itertools import combinations, permutations
from collections import Counter, defaultdict
import sys

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

def is_isomorphic(adj1, adj2, n):
    """Check if two tournaments are isomorphic (brute force for small n)."""
    for perm in permutations(range(n)):
        match = True
        for i in range(n):
            for j in range(n):
                if adj1[i][j] != adj2[perm[i]][perm[j]]:
                    match = False
                    break
            if not match:
                break
        if match:
            return True
    return False

# ══════════════════════════════════════════════════════════════════
# PART 1: Construct QR_7 (Paley tournament on 7 vertices)
# ══════════════════════════════════════════════════════════════════

print("=" * 70)
print("PART 1: THE PALEY TOURNAMENT QR_7")
print("=" * 70)

# QR_7: vertices {0,...,6}, arc i→j iff j-i is a QR mod 7
# QR mod 7: {1, 2, 4} (since 1²=1, 2²=4, 3²=2, etc.)
qr7 = {1, 2, 4}
print(f"Quadratic residues mod 7: {sorted(qr7)}")

n = 7
adj_qr7 = [[0]*n for _ in range(n)]
for i in range(n):
    for j in range(n):
        if i != j and (j - i) % 7 in qr7:
            adj_qr7[i][j] = 1

print("QR_7 adjacency matrix:")
for i in range(n):
    print(f"  {''.join(str(adj_qr7[i][j]) for j in range(n))}")

scores = [sum(adj_qr7[i]) for i in range(n)]
print(f"Scores: {scores} (all 3 → regular)")

cycles_qr7 = find_3cycles(adj_qr7, n)
a2_qr7 = compute_alpha2(cycles_qr7)
H_qr7 = compute_H_dp(adj_qr7, n)

print(f"\nα₁ (3-cycles) = {len(cycles_qr7)}")
print(f"α₂ (disjoint pairs) = {a2_qr7}")
print(f"H (Hamiltonian paths) = {H_qr7}")

# Count transitive vs cyclic triples
total_triples = len(list(combinations(range(n), 3)))
print(f"\nTransitive triples: {total_triples - len(cycles_qr7)}")
print(f"Cyclic triples: {len(cycles_qr7)}")
print(f"Total: {total_triples} = C(7,3)")

# ══════════════════════════════════════════════════════════════════
# PART 2: Enumerate ALL regular tournaments on 7 vertices
# ══════════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("PART 2: ALL REGULAR TOURNAMENTS ON n=7")
print("=" * 70)

# Regular tournament: every vertex has score (n-1)/2 = 3
# Enumerate by trying all tournaments and filtering

edges = [(i,j) for i in range(n) for j in range(i+1, n)]
m = len(edges)  # 21

print(f"Total tournaments: 2^{m} = {1 << m}")
print("Enumerating regular tournaments (score 3 for all vertices)...")

regular_tournaments = []
regular_count = 0

for bits in range(1 << m):
    adj = [[0]*n for _ in range(n)]
    for k, (i,j) in enumerate(edges):
        if bits & (1 << k):
            adj[i][j] = 1
        else:
            adj[j][i] = 1

    # Quick check: all scores = 3
    ok = True
    for i in range(n):
        if sum(adj[i]) != 3:
            ok = False
            break
    if not ok:
        continue

    regular_count += 1
    regular_tournaments.append(adj)

    if regular_count % 10000 == 0:
        print(f"  ... found {regular_count} so far")

print(f"\nTotal regular tournaments on n=7: {regular_count}")

# ── Classify by (α₁, α₂, H) ──────────────────────────────────

print("\nClassifying regular tournaments by (α₁, α₂, H)...")

class_counter = Counter()
iso_classes = []

for idx, adj in enumerate(regular_tournaments):
    cycles = find_3cycles(adj, n)
    a1 = len(cycles)
    a2 = compute_alpha2(cycles)
    H = compute_H_dp(adj, n)
    class_counter[(a1, a2, H)] += 1

    if idx < 100 or a2 == 14:
        # Check isomorphism with known classes
        found_iso = False
        for cls_adj, cls_key in iso_classes:
            if cls_key == (a1, a2, H):
                if is_isomorphic(adj, cls_adj, n):
                    found_iso = True
                    break
        if not found_iso:
            iso_classes.append((adj, (a1, a2, H)))

    if (idx + 1) % 10000 == 0:
        print(f"  ... classified {idx+1}/{regular_count}")

print(f"\n{'α₁':>4} {'α₂':>4} {'H':>6} {'count':>8}")
for key in sorted(class_counter.keys()):
    a1, a2, H = key
    print(f"{a1:>4} {a2:>4} {H:>6} {class_counter[key]:>8}")

print(f"\nRegular tournaments with α₂=14 (max): {sum(v for (a1,a2,H), v in class_counter.items() if a2 == 14)}")
print(f"Regular tournaments with α₂=13: {sum(v for (a1,a2,H), v in class_counter.items() if a2 == 13)}")

# ── Check QR_7 classification ─────────────────────────────────

print(f"\nQR_7 has (α₁,α₂,H) = ({len(cycles_qr7)}, {a2_qr7}, {H_qr7})")

# Is QR_7 the unique tournament with max α₂?
max_a2_vals = [k for k, v in class_counter.items() if k[1] == 14]
if max_a2_vals:
    print(f"Max α₂=14 classes: {max_a2_vals}")

# ══════════════════════════════════════════════════════════════════
# PART 3: The 21 = H_forb_2 connection
# ══════════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("PART 3: THE 21 TRANSITIVE TRIPLES = H_forb_2")
print("=" * 70)

# In a regular tournament on n=7:
# Score sequence [3,3,3,3,3,3,3]
# Number of transitive triples = C(7,3) - (number of 3-cycles)
# For a regular tournament: the number of 3-cycles = ?

# Formula: number of 3-cycles = C(n,3) - sum_i C(s_i, 2)
# where s_i are scores.
# For regular: s_i = 3 for all i, so sum C(3,2) = 7 × 3 = 21
# 3-cycles = C(7,3) - 21 = 35 - 21 = 14

print("For ANY regular tournament on n=7:")
print(f"  3-cycles = C(7,3) - n × C((n-1)/2, 2)")
print(f"           = 35 - 7 × C(3,2)")
print(f"           = 35 - 7 × 3")
print(f"           = 35 - 21 = 14")
print(f"  Transitive triples = 21 = H_forb_2")
print()
print("This is NOT a coincidence!")
print(f"  21 = n × C((n-1)/2, 2) for n=7")
print(f"     = 7 × C(3,2) = 7 × 3 = 21")
print(f"     = C(7,2) = 21")
print()
print("In general for regular tournaments on n vertices (n odd):")
print(f"  transitive triples = n × C((n-1)/2, 2)")
print(f"  For n=3: 3 × C(1,2) = 0 (every triple is a 3-cycle!)")
print(f"  For n=5: 5 × C(2,2) = 5")
print(f"  For n=7: 7 × C(3,2) = 21")
print(f"  For n=9: 9 × C(4,2) = 54")
print()

# But wait: is α₁=14 the SAME for all regular tournaments on n=7?
# Yes! Every regular tournament on 7 has exactly 14 3-cycles.
# So α₁ is a CONSTANT for regular tournaments.
# But α₂ varies!

all_a2 = set()
for key in class_counter:
    all_a2.add(key[1])
print(f"α₂ values among regular n=7 tournaments: {sorted(all_a2)}")
print(f"So α₁=14 is constant but α₂ varies — the TOPOLOGY varies!")

# ══════════════════════════════════════════════════════════════════
# PART 4: The doubly-regular QR_7
# ══════════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("PART 4: QR_7 SPECIAL PROPERTIES")
print("=" * 70)

# QR_7 is doubly-regular: for any two vertices i,j with i→j,
# the number of vertices k with i→k and k→j is constant.

# Check:
in_common = Counter()
for i in range(n):
    for j in range(n):
        if i == j:
            continue
        if adj_qr7[i][j]:  # i→j
            common = sum(1 for k in range(n) if k != i and k != j
                        and adj_qr7[i][k] and adj_qr7[k][j])
            in_common[common] += 1

print(f"Common out-neighbor counts for arcs i→j in QR_7: {dict(in_common)}")
if len(in_common) == 1:
    val = list(in_common.keys())[0]
    print(f"  QR_7 is doubly-regular with parameter λ = {val}")
else:
    print(f"  QR_7 is NOT doubly-regular")

# Automorphism group of QR_7
auto_count = 0
for perm in permutations(range(n)):
    match = True
    for i in range(n):
        for j in range(n):
            if adj_qr7[i][j] != adj_qr7[perm[i]][perm[j]]:
                match = False
                break
        if not match:
            break
    if match:
        auto_count += 1

print(f"|Aut(QR_7)| = {auto_count}")
print(f"  Compare: 7 × 3 = 21 (cyclic × QR translations)")
print(f"  Compare: |GL(3,F_2)| = 168 = 8 × 21")

# ══════════════════════════════════════════════════════════════════
# PART 5: Isomorphism classes
# ══════════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("PART 5: ISOMORPHISM CLASSIFICATION OF max-α₂ TOURNAMENTS")
print("=" * 70)

# Collect all max-α₂ tournaments
max_a2_tourns = [adj for adj in regular_tournaments
                 if compute_alpha2(find_3cycles(adj, n)) == 14]
print(f"Max-α₂=14 tournaments: {len(max_a2_tourns)}")

# Check if QR_7 is among them
qr7_is_max = False
for adj in max_a2_tourns:
    if is_isomorphic(adj, adj_qr7, n):
        qr7_is_max = True
        break
print(f"QR_7 has max α₂: {qr7_is_max}")

# Find isomorphism classes among max-α₂ tournaments
iso_classes_max = []
for adj in max_a2_tourns:
    found = False
    for cls_adj in iso_classes_max:
        if is_isomorphic(adj, cls_adj, n):
            found = True
            break
    if not found:
        iso_classes_max.append(adj)
        if len(iso_classes_max) % 5 == 0:
            print(f"  ... found {len(iso_classes_max)} classes so far")

print(f"\nIsomorphism classes among max-α₂=14 tournaments: {len(iso_classes_max)}")

for i, adj in enumerate(iso_classes_max):
    H = compute_H_dp(adj, n)
    # Count automorphisms
    auto = 0
    for perm in permutations(range(n)):
        if all(adj[a][b] == adj[perm[a]][perm[b]] for a in range(n) for b in range(n)):
            auto += 1
    orbit_size = 5040 // auto  # |S_7| / |Aut|
    print(f"  Class {i}: H={H}, |Aut|={auto}, orbit size={orbit_size}")

total_in_orbits = sum(5040 // sum(1 for perm in permutations(range(n))
                      if all(adj[a][b] == adj[perm[a]][perm[b]]
                             for a in range(n) for b in range(n)))
                      for adj in iso_classes_max)
print(f"  Total from orbit sizes: {total_in_orbits} (should be {len(max_a2_tourns)})")
