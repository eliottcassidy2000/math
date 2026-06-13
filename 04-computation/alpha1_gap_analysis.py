#!/usr/bin/env python3
"""
Full analysis of the alpha_1=3 gap.

Key findings so far:
- n<=4: max alpha_1 = 2, so alpha_1=3 impossible trivially
- n=5: alpha_1=3 impossible (c3=3 forces c5>=1, and c3<=2 with c5<=1 can't reach 3)
- n=6: alpha_1=3 impossible (exhaustive)
- n=7: alpha_1=3 IS POSSIBLE (3360 tournaments)

Questions for this script:
1. At n=5: why exactly does c3=3 force c5>=1?
   Answer: the 3 triples always share a common vertex, span exactly 5 vertices,
   and the induced sub-tournament on those 5 vertices has score (1,1,2,3,3) which
   always has a Hamiltonian cycle = 5-cycle.

2. At n=6: why does alpha_1=3 fail?
   Two sub-questions:
   a) c3=3, c5=0: impossible because triples still span 5 vertices (exhaustive),
      and the 5-vertex sub-tournament still has c5>=1.
   b) c3<=2 with remaining cycles summing to 3-c3: need to check.

3. At n=7: why does it suddenly work?
   Because the triples can span 7 vertices (a disjoint pair), avoiding the
   5-vertex confinement that forces c5>=1.

4. For n>=7: is alpha_1=3 always achievable?

kind-pasteur-2026-03-06
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '03-artifacts', 'code'))
from tournament_lib import tournament_from_bits
from itertools import combinations

def c3_from_scores(scores):
    n = len(scores)
    return (n*(n-1)*(n-2)//6) - sum(s*(s-1)//2 for s in scores)

def count_directed_cycles_on_subset(T, verts):
    k = len(verts)
    if k < 3:
        return 0
    dp = {}
    dp[(1, 0)] = 1
    for mask in range(1, 1 << k):
        if not (mask & 1):
            continue
        for vi in range(k):
            if not (mask & (1 << vi)):
                continue
            c = dp.get((mask, vi), 0)
            if c == 0:
                continue
            for ui in range(k):
                if mask & (1 << ui):
                    continue
                if T[verts[vi]][verts[ui]]:
                    key = (mask | (1 << ui), ui)
                    dp[key] = dp.get(key, 0) + c
    full = (1 << k) - 1
    total = 0
    for vi in range(1, k):
        c = dp.get((full, vi), 0)
        if c > 0 and T[verts[vi]][verts[0]]:
            total += c
    return total

def count_ck(T, n, k):
    total = 0
    for combo in combinations(range(n), k):
        total += count_directed_cycles_on_subset(T, list(combo))
    return total

def get_cyclic_triples(T, n):
    result = []
    for combo in combinations(range(n), 3):
        a, b, c = combo
        if T[a][b] and T[b][c] and T[c][a]:
            result.append(frozenset(combo))
        elif T[a][c] and T[c][b] and T[b][a]:
            result.append(frozenset(combo))
    return result

# ============================================================
# n=5: Full (c3, c5) table
# ============================================================
print("=== n=5: All (c3, c5) pairs ===")
n = 5
m = n*(n-1)//2
pairs = {}
for bits in range(1 << m):
    T = tournament_from_bits(n, bits)
    scores = [sum(T[i]) for i in range(n)]
    c3 = c3_from_scores(scores)
    c5 = count_ck(T, n, 5)
    pair = (c3, c5)
    pairs[pair] = pairs.get(pair, 0) + 1

print(f"  (c3, c5) -> count [alpha_1 = c3+c5]:")
for (c3, c5), cnt in sorted(pairs.items()):
    a1 = c3 + c5
    marker = " <-- alpha_1=3" if a1 == 3 else ""
    print(f"  ({c3}, {c5}): {cnt} tournaments, alpha_1={a1}{marker}")

# ============================================================
# n=6: Full (c3, c5) table
# ============================================================
print(f"\n=== n=6: All (c3, c5) pairs ===")
n = 6
m = n*(n-1)//2
pairs = {}
for bits in range(1 << m):
    T = tournament_from_bits(n, bits)
    scores = [sum(T[i]) for i in range(n)]
    c3 = c3_from_scores(scores)
    if c3 > 5:  # alpha_1 >= c3 > 3, skip
        continue
    c5 = count_ck(T, n, 5)
    pair = (c3, c5)
    pairs[pair] = pairs.get(pair, 0) + 1

print(f"  (c3, c5) -> count [only c3 <= 5 shown]:")
for (c3, c5), cnt in sorted(pairs.items()):
    a1 = c3 + c5
    marker = " <-- alpha_1=3!" if a1 == 3 else ""
    print(f"  ({c3}, {c5}): {cnt} tournaments, alpha_1>={a1}{marker}")

# ============================================================
# n=5: Prove c3=3 => c5>=1 structurally
# ============================================================
print(f"\n=== PROOF: c3=3 forces c5>=1 at n=5 ===")
n = 5
m = n*(n-1)//2
print("Examining all n=5 tournaments with c3=3:")

for bits in range(1 << m):
    T = tournament_from_bits(n, bits)
    scores_list = [sum(T[i]) for i in range(n)]
    c3 = c3_from_scores(scores_list)
    if c3 != 3:
        continue
    scores = tuple(sorted(scores_list))
    triples = get_cyclic_triples(T, n)
    common = triples[0] & triples[1] & triples[2]
    c5 = count_ck(T, n, 5)

    # Find the two vertices shared by all triples
    if len(common) == 2:
        cv = sorted(common)
        print(f"  bits={bits}: scores={scores}, common={cv}, c5={c5}")
        print(f"    T[{cv[0]}][{cv[1]}]={T[cv[0]][cv[1]]}, T[{cv[1]}][{cv[0]}]={T[cv[1]][cv[0]]}")
        break  # Just one example

# ============================================================
# KEY INSIGHT: at n<=6 with c3=3, triples span exactly 5 vertices.
# The 5-vertex sub-tournament has score (1,1,2,3,3).
# Does score (1,1,2,3,3) ALWAYS have c5>=1?
# ============================================================
print(f"\n=== Score (1,1,2,3,3) at n=5: always c5>=1? ===")
n = 5
m = n*(n-1)//2
for bits in range(1 << m):
    T = tournament_from_bits(n, bits)
    scores = tuple(sorted(sum(T[i]) for i in range(n)))
    if scores != (1, 1, 2, 3, 3):
        continue
    c5 = count_ck(T, n, 5)
    c3 = c3_from_scores(list(scores))
    if c5 == 0:
        print(f"  COUNTEREXAMPLE: bits={bits} has score (1,1,2,3,3) but c5=0!")
        break
else:
    print(f"  CONFIRMED: every n=5 tournament with score (1,1,2,3,3) has c5>=1")
    print(f"  (These are exactly the c3=3 tournaments)")

# But wait: at n=6, the c3=3 tournament has 5 active vertices.
# Do those 5 vertices also have score (1,1,2,3,3)?
print(f"\n=== n=6, c3=3: score of the 5-vertex span ===")
n = 6
m = n*(n-1)//2
checked = 0
for bits in range(1 << m):
    T = tournament_from_bits(n, bits)
    scores = [sum(T[i]) for i in range(n)]
    c3 = c3_from_scores(scores)
    if c3 != 3:
        continue
    triples = get_cyclic_triples(T, n)
    all_verts = set()
    for t in triples:
        all_verts.update(t)
    span_verts = sorted(all_verts)
    # Induced score on span
    induced_scores = []
    for v in span_verts:
        s = sum(1 for u in span_verts if u != v and T[v][u])
        induced_scores.append(s)
    induced = tuple(sorted(induced_scores))
    # c5 on those 5 vertices
    c5_span = count_directed_cycles_on_subset(T, span_verts)
    if checked < 3:
        print(f"  bits={bits}, span={span_verts}, induced_scores={induced}, c5_on_span={c5_span}")
    checked += 1
    if induced != (1, 1, 2, 3, 3):
        print(f"  DIFFERENT SCORE: {induced}")
    if c5_span == 0:
        print(f"  c5=0 on span!")

print(f"  Checked {checked} tournaments. All have induced score (1,1,2,3,3) on span? (no error = yes)")

print("\nDone.")
