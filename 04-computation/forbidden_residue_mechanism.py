#!/usr/bin/env python3
"""
forbidden_residue_mechanism.py — opus-2026-03-14-S74

WHY does α₁+2α₂ ≡ 3 (mod 7) fail at n≤6?
WHY does H ≡ 2 (mod 5) fail at n≤5?

The answers must lie in:
1. What (α₁, α₂) pairs are achievable
2. What structural constraints the tournament imposes
3. How the cyclotomic period connects to the lifting point

Also: the surprising α₂=0 at n=5 needs explanation.
"""

from itertools import combinations, permutations
from math import comb
from collections import Counter, defaultdict
import random

def count_hamiltonian_paths(adj, n):
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            if dp[mask][v] == 0:
                continue
            for u in range(n):
                if mask & (1 << u):
                    continue
                if adj[v][u]:
                    dp[mask | (1 << u)][u] += dp[mask][v]
    full = (1 << n) - 1
    return sum(dp[full][v] for v in range(n))

def count_directed_3cycles(adj, n):
    count = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if adj[i][j] and adj[j][k] and adj[k][i]:
                    count += 1
                if adj[i][k] and adj[k][j] and adj[j][i]:
                    count += 1
    return count

def count_directed_5cycles(adj, n):
    count = 0
    for verts in combinations(range(n), 5):
        v = list(verts)
        for perm in permutations(v):
            is_cycle = True
            for i in range(5):
                if not adj[perm[i]][perm[(i+1) % 5]]:
                    is_cycle = False
                    break
            if is_cycle:
                count += 1
    return count // 5

def get_score_sequence(adj, n):
    scores = sorted([sum(adj[i]) for i in range(n)])
    return tuple(scores)

# ====================================================================
# PART 1: WHY α₂ = 0 AT n=5
# ====================================================================
print("=" * 70)
print("PART 1: WHY α₂ = 0 AT n=5")
print("=" * 70)

print("""
  α₂ = number of vertex-disjoint PAIRS of directed odd cycles.

  At n=5, the smallest pair needs:
  - Two 3-cycles: needs 6 vertices (3+3). But n=5 < 6. IMPOSSIBLE!
  - A 3-cycle + 5-cycle: needs 3+5=8 vertices. IMPOSSIBLE!

  Therefore α₂ = 0 at n=5. Always.

  This means I(x) = 1 + α₁·x at n=5 (linear!).
  And H = 1 + 2α₁, so α₁ = (H-1)/2.

  THE REASON: n=5 is too small for vertex-disjoint pairs.
  First α₂ > 0 at n=6 (two disjoint 3-cycles).
""")

# Verify
n = 5
edges = [(i, j) for i in range(n) for j in range(i+1, n)]
all_H = []
for bits in range(2**len(edges)):
    adj = [[0]*n for _ in range(n)]
    for idx, (i, j) in enumerate(edges):
        if bits & (1 << idx):
            adj[i][j] = 1
        else:
            adj[j][i] = 1
    H = count_hamiltonian_paths(adj, n)
    all_H.append(H)

H_set = sorted(set(all_H))
print(f"  H values at n=5: {H_set}")
print(f"  α₁ values: {sorted(set((H-1)//2 for H in H_set))}")
print(f"  H mod 5: {sorted(set(H % 5 for H in H_set))}")
print(f"  Missing mod 5: {[r for r in range(5) if r not in set(H % 5 for H in H_set)]}")

# ====================================================================
# PART 2: ACHIEVABLE α₁ VALUES AT n=5 AND THE GAP
# ====================================================================
print("\n" + "=" * 70)
print("PART 2: THE α₁ GAP — WHY 3 IS MISSING")
print("=" * 70)

# α₁ = dc3 + dc5 at n=5
# We know dc3 ∈ {0,1,2,3,4,5}, dc5 ∈ {0,1,2,3}
# But (dc3, dc5) pairs are constrained

n = 5
edges = [(i, j) for i in range(n) for j in range(i+1, n)]
pair_counts = Counter()
score_analysis = defaultdict(list)

for bits in range(2**len(edges)):
    adj = [[0]*n for _ in range(n)]
    for idx, (i, j) in enumerate(edges):
        if bits & (1 << idx):
            adj[i][j] = 1
        else:
            adj[j][i] = 1
    dc3 = count_directed_3cycles(adj, n)
    dc5 = count_directed_5cycles(adj, n)
    scores = get_score_sequence(adj, n)
    pair_counts[(dc3, dc5)] += 1
    score_analysis[scores].append((dc3, dc5))

print(f"  Achievable (dc3, dc5) pairs at n=5:")
print(f"  {'dc3':>4} {'dc5':>4} {'α₁':>4} {'H':>4} {'H%5':>4} {'H%7':>4} {'count':>6}")
for (dc3, dc5) in sorted(pair_counts.keys()):
    a1 = dc3 + dc5
    H = 1 + 2*a1
    print(f"  {dc3:>4} {dc5:>4} {a1:>4} {H:>4} {H%5:>4} {H%7:>4} {pair_counts[(dc3,dc5)]:>6}")

print(f"\n  α₁ values that DON'T appear: ", end="")
achievable_a1 = set(dc3+dc5 for dc3, dc5 in pair_counts.keys())
max_a1 = max(achievable_a1)
missing = [a for a in range(max_a1+1) if a not in achievable_a1]
print(missing)

print(f"\n  WHY α₁=3 is impossible:")
print(f"  Would need (dc3+dc5)=3, which requires one of:")
print(f"    (3,0): dc3=3 forces dc5≥1 (proven in S73)")
print(f"    (2,1): dc3=2 with dc5=1? Let's check...")

# Check which score sequences allow each dc3 value
print(f"\n  Score sequence → (dc3, dc5) map:")
for scores in sorted(score_analysis.keys()):
    pairs = set(score_analysis[scores])
    # dc3 is determined by score sequence
    dc3_vals = set(p[0] for p in pairs)
    dc5_vals = set(p[1] for p in pairs)
    a1_vals = set(p[0]+p[1] for p in pairs)
    print(f"    {scores}: dc3={dc3_vals}, dc5={dc5_vals}, α₁={a1_vals}")

# ====================================================================
# PART 3: THE SCORE-SEQUENCE CONSTRAINT
# ====================================================================
print("\n" + "=" * 70)
print("PART 3: dc3 IS DETERMINED BY SCORE SEQUENCE")
print("=" * 70)

print("""
  At n=5, the number of directed 3-cycles is:
    dc3 = C(5,3) - Σᵢ C(sᵢ, 2)
  where s₁,...,s₅ is the score sequence (sorted).

  This is because each triple of vertices either:
  - Has a 3-cycle (if no vertex dominates the other two)
  - Has a transitive triple (if some vertex beats both others)

  The transitive triples are counted by Σ C(sᵢ, 2).
  So dc3 = C(5,3) - Σ C(sᵢ, 2) = 10 - Σ C(sᵢ, 2).
""")

# Verify
for scores in sorted(set(get_score_sequence(adj, n)
                         for bits in range(1)  # just use last adj
                         for adj in [[[0]*n for _ in range(n)]])):
    pass  # just to get the structure

print(f"  {'score seq':>20} {'Σ C(sᵢ,2)':>10} {'dc3=10-Σ':>10}")
for scores in sorted(score_analysis.keys()):
    sum_c = sum(s*(s-1)//2 for s in scores)
    dc3_formula = 10 - sum_c
    dc3_actual = list(set(p[0] for p in score_analysis[scores]))[0]
    match = "✓" if dc3_formula == dc3_actual else "✗"
    print(f"  {str(scores):>20} {sum_c:>10} {dc3_formula:>10} {match}")

# ====================================================================
# PART 4: dc5 IS NOT DETERMINED BY SCORE — WHAT DETERMINES IT?
# ====================================================================
print("\n" + "=" * 70)
print("PART 4: WHAT DETERMINES dc5?")
print("=" * 70)

# For each score sequence that allows dc5 variation, analyze
for scores in sorted(score_analysis.keys()):
    dc5_vals = set(p[1] for p in score_analysis[scores])
    if len(dc5_vals) > 1:
        print(f"\n  Score {scores}: dc5 ∈ {sorted(dc5_vals)}")
        # Count
        dc5_dist = Counter(p[1] for p in score_analysis[scores])
        for d5, cnt in sorted(dc5_dist.items()):
            print(f"    dc5={d5}: {cnt} tournaments")

# ====================================================================
# PART 5: EXHAUSTIVE n=6 — ACHIEVABLE (α₁, α₂)
# ====================================================================
print("\n" + "=" * 70)
print("PART 5: ACHIEVABLE (α₁, α₂) AT n=6")
print("=" * 70)

import time
n = 6
num_edges = comb(n, 2)
edges = [(i, j) for i in range(n) for j in range(i+1, n)]

print(f"  Computing (dc3, dc5, H) for all {2**num_edges} tournaments at n=6...")
t0 = time.time()

a1_a2_pairs = Counter()
H_mod7_analysis = defaultdict(list)
a1_2a2_vals = set()

for bits in range(2**num_edges):
    adj = [[0]*n for _ in range(n)]
    for idx, (i, j) in enumerate(edges):
        if bits & (1 << idx):
            adj[i][j] = 1
        else:
            adj[j][i] = 1

    H = count_hamiltonian_paths(adj, n)
    dc3 = count_directed_3cycles(adj, n)
    dc5 = count_directed_5cycles(adj, n)
    alpha1 = dc3 + dc5
    alpha2 = (H - 1 - 2*alpha1) // 4

    a1_a2_pairs[(alpha1, alpha2)] += 1
    a1_2a2_vals.add(alpha1 + 2*alpha2)
    H_mod7_analysis[H % 7].append((alpha1, alpha2, H))

elapsed = time.time() - t0
print(f"  Done in {elapsed:.1f}s")

print(f"\n  Achievable (α₁, α₂) pairs at n=6 (showing first 20):")
print(f"  {'α₁':>5} {'α₂':>5} {'H':>5} {'α₁+2α₂':>8} {'%7':>3} {'count':>6}")
for (a1, a2) in sorted(a1_a2_pairs.keys())[:30]:
    H = 1 + 2*a1 + 4*a2
    v = a1 + 2*a2
    print(f"  {a1:>5} {a2:>5} {H:>5} {v:>8} {v%7:>3} {a1_a2_pairs[(a1,a2)]:>6}")

print(f"\n  α₁+2α₂ values mod 7: {sorted(set(v % 7 for v in a1_2a2_vals))}")
print(f"  Missing mod 7: {[r for r in range(7) if r not in set(v % 7 for v in a1_2a2_vals)]}")

print(f"\n  Range of α₁+2α₂: [{min(a1_2a2_vals)}, {max(a1_2a2_vals)}]")
print(f"  All achievable values: {sorted(a1_2a2_vals)}")

# ====================================================================
# PART 6: WHAT α₁+2α₂ VALUES ARE ACHIEVABLE?
# ====================================================================
print("\n" + "=" * 70)
print("PART 6: STRUCTURE OF ACHIEVABLE α₁+2α₂")
print("=" * 70)

# α₁+2α₂ = dc3+dc5 + 2·(H-1-2(dc3+dc5))/4 = dc3+dc5 + (H-1-2dc3-2dc5)/2
# = dc3+dc5 + H/2 - 1/2 - dc3 - dc5 = H/2 - 1/2 = (H-1)/2
# WAIT: α₁+2α₂ = (H-1)/2 ??

# Let's check: H = 1+2α₁+4α₂, so α₁+2α₂ = (H-1)/2.
# Since H is always odd, (H-1)/2 is always an integer.
# So α₁+2α₂ = (H-1)/2 ALWAYS!

print("""
  CRITICAL INSIGHT:

  α₁+2α₂ = (H-1)/2  ALWAYS!

  Proof: H = 1 + 2α₁ + 4α₂ = 1 + 2(α₁ + 2α₂)
  So α₁ + 2α₂ = (H-1)/2.

  Therefore:
    α₁+2α₂ ≡ 3 (mod 7)  ⟺  (H-1)/2 ≡ 3 (mod 7)  ⟺  H ≡ 7 ≡ 0 (mod 7)

  So H ≡ 0 (mod 7) IS EQUIVALENT TO α₁+2α₂ ≡ 3 (mod 7).
  This is just a restatement, not a deeper constraint!

  Similarly:
    H ≡ 2 (mod 5) ⟺ (H-1)/2 ≡ (2-1)/2... wait, need to be careful.
    H ≡ 2 (mod 5): H = 5k+2, (H-1)/2 = (5k+1)/2.
    This needs 5k+1 even, so k odd. H = 10m+7. (H-1)/2 = 5m+3.
    So (H-1)/2 ≡ 3 (mod 5).

    Or: H ≡ 2 (mod 5) means (H-1)/2 ≡ (H-1)/2 mod 5.
    H=2: (H-1)/2 = 1/2 ≡ 3 (mod 5) [since 2·3=6≡1]
    H=7: (H-1)/2 = 3 ≡ 3 (mod 5) ✓
    H=12: (H-1)/2 = 11/2... H=12 is even, impossible (Rédei).

  Actually H is always odd, so H ≡ 2 (mod 5) means H ∈ {7, 17, 27, ...}.
  (H-1)/2 ∈ {3, 8, 13, ...} ≡ 3 (mod 5).

  So: H ≡ 0 (mod 7) ⟺ (H-1)/2 ≡ 3 (mod 7)
      H ≡ 2 (mod 5) ⟺ (H-1)/2 ≡ 3 (mod 5)

  BOTH forbidden residues correspond to (H-1)/2 ≡ 3 mod p!
  (for p=7 and p=5 respectively)
""")

# Verify
print("  Verification:")
for H in [1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21]:
    half = (H-1)//2
    print(f"    H={H:>3}: (H-1)/2={half:>3}, mod 5={half%5}, mod 7={half%7}")

print(f"\n  The forbidden residue (H-1)/2 ≡ 3 (mod p) for p=5:")
print(f"  Achievable (H-1)/2 mod 5 at n=5: ", end="")
r5 = sorted(set((H-1)//2 % 5 for H in all_H))
print(r5, "missing:", [r for r in range(5) if r not in r5])

print(f"  Achievable (H-1)/2 mod 7 at n=5: ", end="")
r7 = sorted(set((H-1)//2 % 7 for H in all_H))
print(r7, "missing:", [r for r in range(7) if r not in r7])

# ====================================================================
# PART 7: THE DEEP QUESTION — WHY IS (H-1)/2 ≡ 3 AVOIDED?
# ====================================================================
print("\n" + "=" * 70)
print("PART 7: WHY (H-1)/2 ≡ 3 IS SPECIAL")
print("=" * 70)

print("""
  At n=5: (H-1)/2 ∈ {0, 1, 2, 4, 5, 6, 7} — skips 3!

  H ∈ {1, 3, 5, 9, 11, 13, 15}
  (H-1)/2 ∈ {0, 1, 2, 4, 5, 6, 7}

  The gap at 3 means H=7 is impossible at n=5.
  7 = Φ₃(2)! The cyclotomic prime itself is the forbidden H value!

  WHY is H=7 impossible at n=5?
  H=7 would need α₁=3 (since α₂=0, H=1+2α₁).
  α₁=3 is impossible because dc3+dc5=3 is impossible.

  At n=6: what H values near 7 are achievable?
""")

# Check H=7 at n=6
H_at_6 = set()
for bits in range(2**comb(6,2)):
    adj = [[0]*6 for _ in range(6)]
    edges6 = [(i,j) for i in range(6) for j in range(i+1,6)]
    for idx, (i,j) in enumerate(edges6):
        if bits & (1 << idx):
            adj[i][j] = 1
        else:
            adj[j][i] = 1
    H = count_hamiltonian_paths(adj, 6)
    H_at_6.add(H)

H_at_6 = sorted(H_at_6)
print(f"  H values at n=6: {H_at_6[:20]}... (showing first 20)")
print(f"  Is H=7 achievable? {7 in H_at_6}")
print(f"  H mod 7 at n=6: {sorted(set(H % 7 for H in H_at_6))}")
print(f"  Missing mod 7: {[r for r in range(7) if r not in set(H%7 for H in H_at_6)]}")
print(f"  Is H=14 achievable? {14 in H_at_6}")
print(f"  Is H=21 achievable? {21 in H_at_6}")
print(f"  Is H=35 achievable? {35 in H_at_6}")

# Which H ≡ 0 mod 7 are achievable at n=6?
h_mod7_0 = sorted([H for H in H_at_6 if H % 7 == 0])
print(f"  H ≡ 0 (mod 7) at n=6: {h_mod7_0}")

# ====================================================================
# PART 8: THE MECHANISM — SMALL n RESTRICTIONS ON (H-1)/2
# ====================================================================
print("\n" + "=" * 70)
print("PART 8: GROWTH OF (H-1)/2 RANGE")
print("=" * 70)

for nn in [3, 4, 5, 6]:
    edges_n = [(i,j) for i in range(nn) for j in range(i+1,nn)]
    H_vals = set()
    for bits in range(2**len(edges_n)):
        adj = [[0]*nn for _ in range(nn)]
        for idx, (i,j) in enumerate(edges_n):
            if bits & (1 << idx):
                adj[i][j] = 1
            else:
                adj[j][i] = 1
        H_vals.add(count_hamiltonian_paths(adj, nn))

    H_sorted = sorted(H_vals)
    half_vals = sorted(set((H-1)//2 for H in H_sorted))
    max_half = max(half_vals)

    missing_all = [v for v in range(max_half+1) if v not in half_vals]

    print(f"  n={nn}: H ∈ {H_sorted}")
    print(f"       (H-1)/2 ∈ {half_vals}")
    print(f"       missing: {missing_all if missing_all else 'none'}")
    print(f"       max (H-1)/2 = {max_half}")
    print(f"       missing mod 5: {[r for r in range(5) if r not in set(v%5 for v in half_vals)]}")
    print(f"       missing mod 7: {[r for r in range(7) if r not in set(v%7 for v in half_vals)]}")
    print()

# ====================================================================
# PART 9: WHEN DOES THE GAP CLOSE?
# ====================================================================
print("=" * 70)
print("PART 9: THE PARITY OF (H-1)/2 AND TOURNAMENT STRUCTURE")
print("=" * 70)

print("""
  At n=5:  (H-1)/2 ∈ {0,1,2,4,5,6,7} — gap at 3
  At n=6:  need to check if gap at (H-1)/2 ≡ 3 (mod 7) persists

  The gap MUST close eventually because:
  - H grows exponentially with n (max H ~ n!/e)
  - The achievable H values become denser
  - The arithmetic progressions fill up

  The question is: WHEN does each gap close?
  Answer from previous results:
    H≡2 mod 5: lifts at n=6 (first H with (H-1)/2 ≡ 3 mod 5)
    H≡0 mod 7: lifts at n=7 (first H with (H-1)/2 ≡ 3 mod 7)
""")

# What is the smallest H ≡ 0 (mod 7) at n=7?
random.seed(42)
n = 7
h_mod7_0_at_7 = set()
for trial in range(50000):
    adj = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5:
                adj[i][j] = 1
            else:
                adj[j][i] = 1
    H = count_hamiltonian_paths(adj, n)
    if H % 7 == 0:
        h_mod7_0_at_7.add(H)

if h_mod7_0_at_7:
    print(f"  Smallest H ≡ 0 (mod 7) found at n=7: {min(h_mod7_0_at_7)}")
    print(f"  All H ≡ 0 (mod 7) found: {sorted(h_mod7_0_at_7)[:10]}...")
else:
    print(f"  No H ≡ 0 (mod 7) found in 50000 samples at n=7!")

# ====================================================================
# PART 10: SYNTHESIS — THE FORBIDDEN 3
# ====================================================================
print("\n" + "=" * 70)
print("PART 10: SYNTHESIS — THE FORBIDDEN 3")
print("=" * 70)

print("""
  THE PATTERN: (H-1)/2 ≡ 3 (mod p) is the forbidden residue for
  BOTH p=5 and p=7 (and implicitly for p=3 it's (H-1)/2 ≡ 1 mod 2,
  which would mean H ≡ 3 mod 4, which IS achievable).

  This is because:
    H ≡ 0 (mod p) ⟺ (H-1)/2 ≡ (p-1)/2 (mod p)

  For p=3: (p-1)/2 = 1. H≡0 mod 3 ⟺ (H-1)/2 ≡ 1 mod 3.
  For p=5: (p-1)/2 = 2. H≡0 mod 5 ⟺ (H-1)/2 ≡ 2 mod 5.
  For p=7: (p-1)/2 = 3. H≡0 mod 7 ⟺ (H-1)/2 ≡ 3 mod 7.

  Wait — H≡2 mod 5, not H≡0 mod 5. Let me recheck.

  H≡2 mod 5: (H-1)/2 = (5k+1)/2. For H odd: 5k+2 odd → 5k odd → k odd.
  H=7: (H-1)/2=3≡3 mod 5. H=17: (H-1)/2=8≡3 mod 5.
  So H≡2 mod 5 ⟹ (H-1)/2 ≡ 3 mod 5? Let me check more carefully.

  Actually:
  (H-1)/2 mod 5:
  H=1: 0
  H=3: 1
  H=5: 2
  H=7: 3  ← H≡2 mod 5, (H-1)/2≡3 mod 5
  H=9: 4
  H=11: 0
  H=13: 1
  H=15: 2
  H=17: 3  ← H≡2 mod 5, (H-1)/2≡3 mod 5 ✓

  Yes! H≡2 mod 5 ⟺ (H-1)/2 ≡ 3 mod 5.
  And H≡0 mod 7 ⟺ (H-1)/2 ≡ 3 mod 7.

  BOTH forbidden conditions give (H-1)/2 ≡ 3!

  The number 3 itself is the universal forbidden residue for (H-1)/2.

  WHY 3? Because at n=5, α₁=3 is the impossible value (dc3+dc5 can't sum to 3).
  And α₁ = (H-1)/2 when α₂=0.

  The impossibility of α₁=3 at n=5 creates the same gap mod 5 and mod 7!
  It lifts at different n values because:
  - Mod 5: the gap fills when α₂ > 0 compensates (n=6)
  - Mod 7: the gap persists because α₁+2α₂ ≡ 3 (mod 7) is ALSO impossible at n=6
""")

# Verify the universal forbidden 3
print("  Universal pattern verification:")
for p in [3, 5, 7, 11, 13]:
    # What H mod p gives (H-1)/2 ≡ 3 mod p?
    # (H-1)/2 ≡ 3 mod p → H ≡ 2·3+1 ≡ 7 mod 2p... no.
    # More carefully: H-1 ≡ 6 mod 2p if p is odd
    # H ≡ 7 mod 2p → H mod p = 7 mod p
    forbidden_H_mod_p = (7 % p) if p > 7 else (7 % p)
    # Actually (H-1)/2 ≡ 3 mod p means H = 2(pk+3)+1 = 2pk+7
    # So H ≡ 7 mod 2p. And H mod p ≡ 7 mod p.
    h_mod_p = 7 % p
    # But H is odd, so we also need 7 mod 2 = 1 ✓ (7 is odd)
    print(f"  p={p}: (H-1)/2 ≡ 3 mod {p} ⟺ H ≡ {h_mod_p} mod {p}", end="")
    # Check if this is indeed a forbidden residue at n=5
    h5_mods = set(H % p for H in all_H)
    forbidden = h_mod_p not in h5_mods
    print(f"  {'← FORBIDDEN at n=5' if forbidden else '← achievable at n=5'}")

print("\nDone.")
