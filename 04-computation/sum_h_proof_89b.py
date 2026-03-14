#!/usr/bin/env python3
"""
PROOF OF Sum H = n! · 2^{m-n+1} AND DEEPER PALEY CONNECTIONS
opus-2026-03-14-S89b

The identity Sum_{T} H(T) = n! · 2^{m-n+1} where m = C(n,2)
has a beautiful proof:

PROOF:
  Sum_T H(T) = Sum_T Sum_π Π_{i=1}^{n-1} A_{π(i),π(i+1)}

  Swap summation order:
  = Sum_π Sum_T Π_{i=1}^{n-1} A_{π(i),π(i+1)}

  For a fixed permutation π, the product Π A_{π(i),π(i+1)} is 1
  iff the tournament has arcs π(1)→π(2)→...→π(n).

  This constrains n-1 arcs out of m = C(n,2) arcs. The remaining
  m - (n-1) arcs are free. So:

  Sum_T Π A = 2^{m-n+1}

  And there are n! permutations π:
  Sum_T H(T) = n! · 2^{m-n+1}. □

This is elegant and completely independent of the IP formula!

Now: what does this mean for E[i_k]?
  E[H] = n!/2^{n-1} = n!·2^{m-n+1}/2^m

Using H = Σ_k 2^k · i_k:
  E[H] = Σ_k 2^k · E[i_k]
  n!/2^{n-1} = 1 + 2·E[total cycles] + 4·E[disjoint pairs] + ...

We know:
  E[t₃] = C(n,3)/4
  E[t₅] = C(n,5) · (3/4)  ← where does 3/4 come from?

Let me verify and find exact formulas for E[t_k] and E[d_ij].
"""

from math import factorial, comb
from itertools import combinations, permutations
from collections import Counter

def compute_H(n, adj):
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)): continue
            if (mask, v) not in dp: continue
            val = dp[(mask, v)]
            for u in range(n):
                if mask & (1 << u): continue
                if adj[v][u]:
                    key = (mask | (1 << u), u)
                    dp[key] = dp.get(key, 0) + val
    full = (1 << n) - 1
    return sum(dp.get((full, v), 0) for v in range(n))

def tournament_adj(n, bits):
    adj = [[False]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if (bits >> idx) & 1:
                adj[i][j] = True
            else:
                adj[j][i] = True
            idx += 1
    return adj

print("="*70)
print("Sum H = n! · 2^{m-n+1}: PROOF AND CONSEQUENCES")
print("="*70)

# Verify the identity
print("\nVerification:")
for n in range(3, 8):
    m = n*(n-1)//2
    expected = factorial(n) * (1 << (m - n + 1))
    if n <= 6:
        actual = sum(compute_H(n, tournament_adj(n, b)) for b in range(1 << m))
        print(f"  n={n}: m={m}, Sum H = {actual}, n!·2^(m-n+1) = {expected}, match = {actual==expected}")
    else:
        print(f"  n={n}: m={m}, n!·2^(m-n+1) = {expected} (not verified exhaustively)")

# ===== E[t_k] formulas =====
print("\n" + "="*70)
print("EXACT E[t_k] FORMULAS")
print("="*70)

# E[t_3] = C(n,3) · p(3-cycle)
# p(3-cycle on 3 vertices) = 2/2^3 = 1/4
# Because: 2 of the 8 tournaments on 3 vertices are 3-cycles

# E[t_5] = C(n,5) · p(5-cycle on 5 vertices)
# p(5-cycle on 5 vertices) = (number of directed 5-cycles) / 2^10

# Count directed 5-cycles on 5 labeled vertices
n5 = 5
m5 = n5*(n5-1)//2
total_5cycles = 0
for bits in range(1 << m5):
    adj = tournament_adj(n5, bits)
    count = 0
    for perm in permutations(range(n5)):
        ok = True
        for idx in range(5):
            if not adj[perm[idx]][perm[(idx+1)%5]]:
                ok = False
                break
        if ok:
            count += 1
    total_5cycles += count // 5

avg_5_per_tournament = total_5cycles / (1 << m5)
prob_5cycle = avg_5_per_tournament  # for the single 5-subset {0,1,2,3,4}

print(f"\nFor a random tournament on 5 vertices:")
print(f"  Total directed 5-cycles across all tournaments: {total_5cycles}")
print(f"  Average t₅ per tournament: {avg_5_per_tournament}")
print(f"  Expected: ?")

# The number of labeled directed 5-cycles on 5 vertices
# = (5-1)! = 24 (fix a start, permute the rest)
# But each tournament supports only some of these
# The probability that a SPECIFIC directed 5-cycle exists in a random tournament:
# = Pr(all 5 arcs go in the cycle direction) = (1/2)^5 = 1/32
# But a vertex set of 5 has C(5,2)=10 arcs, and a 5-cycle constrains only 5 of them
# So Pr(specific directed 5-cycle) = (1/2)^5 = 1/32
# Number of labeled directed 5-cycles = (5-1)!/2? No...

# Actually: number of distinct directed 5-cycles on {0,1,2,3,4}:
# = 4!/2 = 12? Let me count.
# Fix direction: cyclic permutations of (a,b,c,d,e) give the same cycle.
# So number = 5!/5 = 24. But these include both clockwise and counterclockwise.
# In a tournament, both directions can't coexist (arc i->j is fixed).
# So the number of POTENTIAL directed 5-cycles = 24.
# Each has probability (1/2)^5 = 1/32 of existing.

# Expected t₅ on 5 vertices = 24 · (1/2)^5 = 24/32 = 3/4 = 0.75 ✓

print(f"\n  Derivation:")
print(f"    Directed 5-cycles on 5 labeled vertices: (5-1)! = {factorial(4)} = 24")
print(f"    Each constrains 5 arcs, probability = (1/2)^5 = 1/32")
print(f"    E[t₅ on 5 verts] = 24/32 = 3/4 = {24/32}")
print(f"    Matches: {abs(avg_5_per_tournament - 24/32) < 0.001}")

print(f"\n  For general n:")
print(f"    E[t₃] = C(n,3) · 2/2³ = C(n,3)/4")
print(f"    E[t₅] = C(n,5) · 24/2^10 = C(n,5) · 3/128")
print(f"    E[t₇] = C(n,7) · 720/2^21 = C(n,7) · 720/2097152")

# Verify E[t₅] for n=6
n6 = 6
m6 = n6*(n6-1)//2
total_t5_n6 = 0
for bits in range(1 << m6):
    adj = tournament_adj(n6, bits)
    t5 = 0
    for combo in combinations(range(n6), 5):
        for perm in permutations(combo):
            ok = True
            for idx in range(5):
                if not adj[perm[idx]][perm[(idx+1)%5]]:
                    ok = False
                    break
            if ok:
                t5 += 1
    total_t5_n6 += t5 // 5

avg_t5_n6 = total_t5_n6 / (1 << m6)
expected_t5_n6 = comb(6, 5) * 24 / 2**10
print(f"\n  n=6: E[t₅] = {avg_t5_n6}, expected = C(6,5)·24/2^10 = {expected_t5_n6}")
print(f"  Match: {abs(avg_t5_n6 - expected_t5_n6) < 0.001}")

# ===== THE MASTER EQUATION =====
print("\n" + "="*70)
print("THE MASTER EQUATION")
print("="*70)

print("""
From E[H] = n!/2^{n-1} and H = Σ_k 2^k · i_k:

  n!/2^{n-1} = 1 + 2·E[Σ t_odd] + 4·E[Σ d_ij] + ...

For small n:
  n=3: 3/2 = 1 + 2·(1/4) = 1 + 1/2 = 3/2 ✓
  n=4: 3 = 1 + 2·1 = 3 ✓  (E[t3] = C(4,3)/4 = 1)
  n=5: 15/2 = 1 + 2·(5/2 + 3/4) = 1 + 2·13/4 = 1 + 13/2 = 15/2 ✓
""")

# Verify the master equation for n=3,4,5
for n in [3, 4, 5]:
    m = n*(n-1)//2
    mean_H = factorial(n) / 2**(n-1)

    # E[t3]
    et3 = comb(n, 3) / 4

    # E[t5] (only for n >= 5)
    et5 = comb(n, 5) * 24 / 2**10 if n >= 5 else 0

    # E[d33] for n=5: impossible (need 6 vertices)
    ed33 = 0  # for n <= 5

    rhs = 1 + 2*(et3 + et5) + 4*ed33
    print(f"  n={n}: E[H] = {mean_H}, 1 + 2·(E[t3]+E[t5]) + 4·E[d33] = 1 + 2·({et3}+{et5}) + 0 = {rhs}")
    print(f"    Match: {abs(mean_H - rhs) < 0.001}")

# For n=6, we need E[d33]
print(f"\n  For n=6:")
n = 6
m = n*(n-1)//2
mean_H = factorial(n) / 2**(n-1)
et3 = comb(n, 3) / 4
et5 = comb(n, 5) * 24 / 2**10

# Compute E[d33] directly
total_d33 = 0
for bits in range(1 << m):
    adj = tournament_adj(n, bits)
    triples = []
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if (adj[i][j] and adj[j][k] and adj[k][i]) or \
                   (adj[i][k] and adj[k][j] and adj[j][i]):
                    triples.append(frozenset([i,j,k]))
    d33 = 0
    for i in range(len(triples)):
        for j in range(i+1, len(triples)):
            if len(triples[i] & triples[j]) == 0:
                d33 += 1
    total_d33 += d33

ed33 = total_d33 / (1 << m)
rhs6 = 1 + 2*(et3 + et5) + 4*ed33
print(f"  E[t3] = {et3}, E[t5] = {et5}, E[d33] = {ed33}")
print(f"  1 + 2·({et3}+{et5}) + 4·{ed33} = {rhs6}")
print(f"  E[H] = {mean_H}")
print(f"  Match: {abs(mean_H - rhs6) < 0.001}")

# Can we derive E[d33] from combinatorics?
# E[d33] = E[number of disjoint 3-cycle pairs]
# = C(n,3) · C(n-3,3) / 2 · Pr(both triples are 3-cycles)
# = [C(n,3)·C(n-3,3)/2] · (1/4)^2 = [C(n,3)·C(n-3,3)/2] · 1/16

# For n=6: C(6,3)·C(3,3)/2 · 1/16 = 20·1/2 · 1/16 = 10/16 = 5/8
expected_ed33 = comb(6,3) * comb(3,3) / 2 * (1/4)**2
print(f"\n  E[d33] derivation: C(6,3)·C(3,3)/2 · (1/4)² = {expected_ed33}")
print(f"  Computed: {ed33}")
print(f"  Match: {abs(ed33 - expected_ed33) < 0.001}")

# ===== PALEY TOURNAMENT DEEP DIVE =====
print("\n" + "="*70)
print("PALEY TOURNAMENT DEEP DIVE")
print("="*70)

# P_7: d33 = 7. The 7 disjoint pairs correspond to the 7 partitions
# of {0,1,2,3,4,5,6} into two disjoint triples + one leftover vertex?
# No, 3+3 = 6, leaving 1 vertex. So each disjoint pair uses 6 of 7 vertices.

# Which vertex is left out in each of the 7 disjoint pairs?
qr7 = {1, 2, 4}
n = 7
adj_paley = [[False]*n for _ in range(n)]
for i in range(n):
    for j in range(n):
        if i != j and (j - i) % n in qr7:
            adj_paley[i][j] = True

triples_paley = []
for i in range(n):
    for j in range(i+1, n):
        for k in range(j+1, n):
            if (adj_paley[i][j] and adj_paley[j][k] and adj_paley[k][i]) or \
               (adj_paley[i][k] and adj_paley[k][j] and adj_paley[j][i]):
                triples_paley.append(frozenset([i,j,k]))

print(f"\nPaley P₇: {len(triples_paley)} 3-cycles")

# Find all 7 disjoint pairs
disjoint_pairs = []
for i in range(len(triples_paley)):
    for j in range(i+1, len(triples_paley)):
        if len(triples_paley[i] & triples_paley[j]) == 0:
            disjoint_pairs.append((triples_paley[i], triples_paley[j]))

print(f"Disjoint 3-cycle pairs: {len(disjoint_pairs)}")
for pair in disjoint_pairs:
    leftover = set(range(7)) - pair[0] - pair[1]
    print(f"  {sorted(pair[0])} + {sorted(pair[1])}, leftover = {sorted(leftover)}")

# Count how many times each vertex is the "leftover"
leftover_count = Counter()
for pair in disjoint_pairs:
    leftover = set(range(7)) - pair[0] - pair[1]
    for v in leftover:
        leftover_count[v] += 1

print(f"\nLeftover vertex distribution: {dict(sorted(leftover_count.items()))}")
print(f"Each vertex is leftover the same number of times: {len(set(leftover_count.values())) == 1}")

# ===== CONNECTION TO STEINER SYSTEMS =====
print("\n" + "="*70)
print("STEINER SYSTEM CONNECTION")
print("="*70)

print("""
The 14 three-cycles of P₇ form a structure where:
- 7 vertices
- 14 triples (twice the Fano plane)
- Each pair of vertices in exactly 2 triples
- This is a 2-(7, 3, 2) design!

The disjoint pairs (d₃₃ = 7) correspond to:
- Each disjoint pair covers 6 vertices, leaving 1
- Each vertex is left out exactly once
- So the 7 pairs partition into 7 "complementary triads"

This is the COMPLEMENT of the Fano plane:
- Fano: 7 blocks, each pair in 1 block
- Our structure: 14 blocks, each pair in 2 blocks = 2× Fano

The 7 disjoint pairs, one per leftover vertex, form a
RESOLUTION of the 2-(7,3,2) design!
""")

# Check: is each pair in exactly 2 triples?
from collections import defaultdict
pair_in_triples = defaultdict(int)
for t in triples_paley:
    for pair in combinations(t, 2):
        pair_in_triples[frozenset(pair)] += 1

print(f"Pair multiplicities: {sorted(set(pair_in_triples.values()))}")
print(f"All pairs in exactly 2 triples: {all(v == 2 for v in pair_in_triples.values())}")

print("\n" + "="*70)
print("DONE — Sum H proof and Paley connections")
print("="*70)
