#!/usr/bin/env python3
"""
betti_criterion.py - Find a simple criterion for β₁ > 0 and β₃ > 0.

KEY DATA from exhaustive n=5:
  t₃=0: β₁=0 always (120/120)
  t₃=1: β₁=0 always (120/120)
  t₃=2: β₁=0 always (240/240)
  t₃=3: β₁=1 for 120/240 (50%)
  t₃=4: β₁=1 for 160/280 (57%)
  t₃=5: β₁=1 for 24/24 (100%)

So β₁ is NOT just a function of t₃. What else matters?

HYPOTHESIS 1: β₁ = 1 iff Ω(T) contains a specific substructure.
HYPOTHESIS 2: β₁ relates to the "homological complexity" of the 3-cycle graph.
HYPOTHESIS 3: β₁ = 1 iff the tournament "wraps around" in a topological sense.

Let me check: at n=4, β₁=1 iff t₃=2. At n=3, β₁=1 iff t₃=1.
The threshold is t₃ = C(n,3) - sum C(s_i,2) for... the regular tournament?

Actually at n=4: t₃ can be 0 (transitive, 24 tournaments) or 1 (almost-transitive, 16)
or 2 (the two "cycle-rich" ones, 24). Wait, but earlier data shows:
  n=4: t₃=0: {β₁=0: 24}, t₃=1: {β₁=0: 16}, t₃=2: {β₁=1: 24}
So β₁=1 iff t₃=2 iff the tournament has TWO directed 3-cycles.

At n=3: β₁=1 iff t₃=1 (the 3-cycle). Max t₃ at n=3 is 1.
At n=4: β₁=1 iff t₃=2. Max t₃ at n=4 is 4... wait no.
  n=4 has C(4,3)=4 triples. t₃ = 4 - sum C(s_i,2).
  Score (0,1,2,3): sum C(s_i,2) = 0+0+1+3 = 4, t₃=0.
  Score (1,1,1,3): sum = 0+0+0+3 = 3, t₃=1.
  Score (1,1,2,2): sum = 0+0+1+1 = 2, t₃=2.
  Max t₃ at n=4 = 2.

So at n=4, β₁=1 iff t₃ = max(t₃) = 2. At n=3, β₁=1 iff t₃=max=1.
But at n=5, β₁=1 doesn't require max t₃=5.

Let me investigate what property WITHIN each t₃ class determines β₁.

Author: opus-2026-03-07-S46e
"""
import sys
sys.path.insert(0, '/Users/e/Documents/GitHub/math/04-computation')
from path_homology_v2 import path_betti_numbers
from itertools import combinations, permutations
from collections import Counter, defaultdict

def count_t3(A, n):
    return sum(1 for i,j,k in combinations(range(n), 3)
               if (A[i][j] and A[j][k] and A[k][i]) or
                  (A[i][k] and A[k][j] and A[j][i]))

def score_seq(A, n):
    return tuple(sorted(sum(A[i][j] for j in range(n)) for i in range(n)))

def is_strongly_connected(A, n):
    """Check if tournament is strongly connected (has a directed Hamiltonian cycle)."""
    # BFS reachability from vertex 0
    visited = {0}
    queue = [0]
    while queue:
        v = queue.pop(0)
        for u in range(n):
            if A[v][u] and u not in visited:
                visited.add(u)
                queue.append(u)
    if len(visited) < n:
        return False
    # Reverse reachability
    visited = {0}
    queue = [0]
    while queue:
        v = queue.pop(0)
        for u in range(n):
            if A[u][v] and u not in visited:
                visited.add(u)
                queue.append(u)
    return len(visited) == n

# === Exhaustive n=5 ===
print("=" * 60)
print("n=5: WHAT DETERMINES β₁ WITHIN EACH t₃ CLASS?")
print("=" * 60)

n = 5
m = n*(n-1)//2
data = []

for bits in range(1 << m):
    A = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if (bits >> idx) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    beta = path_betti_numbers(A, n)
    b1 = int(beta[1]) if len(beta) > 1 else 0
    t3 = count_t3(A, n)
    sc = is_strongly_connected(A, n)
    score = score_seq(A, n)
    data.append({'b1': b1, 't3': t3, 'sc': sc, 'score': score})

# β₁ vs strong connectivity
print("\nβ₁ vs strong connectivity:")
tab = defaultdict(Counter)
for d in data:
    tab[d['sc']][d['b1']] += 1
for sc in [True, False]:
    print(f"  SC={sc}: β₁ distribution = {dict(tab[sc])}")

# β₁ vs (t₃, SC)
print("\nβ₁ vs (t₃, SC):")
tab2 = defaultdict(Counter)
for d in data:
    tab2[(d['t3'], d['sc'])][d['b1']] += 1
for key in sorted(tab2.keys()):
    total = sum(tab2[key].values())
    b1_1 = tab2[key].get(1, 0)
    print(f"  t₃={key[0]}, SC={key[1]}: β₁=1 in {b1_1}/{total} ({100*b1_1/total:.0f}%)")

# β₁ vs score sequence
print("\nβ₁ vs score sequence:")
tab3 = defaultdict(Counter)
for d in data:
    tab3[d['score']][d['b1']] += 1
for score in sorted(tab3.keys()):
    total = sum(tab3[score].values())
    b1_1 = tab3[score].get(1, 0)
    if b1_1 > 0 or total > 50:
        print(f"  score={score}: β₁=1 in {b1_1}/{total} ({100*b1_1/total:.0f}%)")

# KEY QUESTION: is β₁ = 1 iff strongly connected?
b1_sc = sum(1 for d in data if d['b1'] == 1 and d['sc'] == True)
b1_nsc = sum(1 for d in data if d['b1'] == 1 and d['sc'] == False)
sc_total = sum(1 for d in data if d['sc'] == True)
nsc_total = sum(1 for d in data if d['sc'] == False)

print(f"\nIs β₁ = 1 iff SC?")
print(f"  SC + β₁=1: {b1_sc}/{sc_total}")
print(f"  NSC + β₁=1: {b1_nsc}/{nsc_total}")
if b1_nsc == 0:
    print("  β₁ > 0 implies SC (necessary)")
if b1_sc == sc_total:
    print("  SC implies β₁ > 0 (sufficient)")

# === Exhaustive n=4 ===
print("\n" + "=" * 60)
print("n=4: VERIFICATION")
print("=" * 60)

n = 4
m = n*(n-1)//2
data4 = []

for bits in range(1 << m):
    A = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if (bits >> idx) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    beta = path_betti_numbers(A, n)
    b1 = int(beta[1]) if len(beta) > 1 else 0
    t3 = count_t3(A, n)
    sc = is_strongly_connected(A, n)
    data4.append({'b1': b1, 't3': t3, 'sc': sc})

tab4 = defaultdict(Counter)
for d in data4:
    tab4[(d['t3'], d['sc'])][d['b1']] += 1
print("n=4: β₁ vs (t₃, SC)")
for key in sorted(tab4.keys()):
    total = sum(tab4[key].values())
    b1_1 = tab4[key].get(1, 0)
    print(f"  t₃={key[0]}, SC={key[1]}: β₁=1 in {b1_1}/{total}")

# At n=4, SC iff t₃ >= 1 (no, SC iff there's a Hamiltonian cycle)
# Actually a tournament is SC iff it's not a "star" (has more than one SCC)
# A 4-tournament is SC iff not transitive => has at least one 3-cycle
# Wait, is that right? Score (1,1,1,3) has t₃=1 and vertex 3 beats everyone.
# Can we reach 3 from other vertices? 0→? 1→? 2→? They all point somewhere.
# Actually a tournament is SC iff every vertex can reach every other.

sc4 = sum(1 for d in data4 if d['sc'] == True)
print(f"\nn=4: {sc4}/64 are SC")

# === Exhaustive n=3 ===
print("\n" + "=" * 60)
print("n=3: VERIFICATION")
print("=" * 60)

n = 3
data3 = []
for bits in range(1 << 3):
    A = [[0]*3 for _ in range(3)]
    idx = 0
    for i in range(3):
        for j in range(i+1, 3):
            if (bits >> idx) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    beta = path_betti_numbers(A, 3)
    b1 = int(beta[1]) if len(beta) > 1 else 0
    t3 = count_t3(A, 3)
    sc = is_strongly_connected(A, 3)
    data3.append({'b1': b1, 't3': t3, 'sc': sc})
    print(f"  t₃={t3}, SC={sc}, β₁={b1}")

print(f"\nn=3: β₁=1 iff t₃=1 iff SC (the 3-cycle)")

print("\n" + "=" * 60)
print("SUMMARY")
print("=" * 60)
