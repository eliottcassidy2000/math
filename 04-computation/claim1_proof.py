#!/usr/bin/env python3
"""
Prove Claim 1: c3 <= 2 implies c5 = 0.

Approach: A 5-vertex subset has a directed Hamiltonian cycle (contributing to c5)
iff the induced sub-tournament is strongly connected.
A strongly connected tournament on 5 vertices has c3 >= 3.
So if the global tournament has c3 <= 2, no 5-subset is strongly connected,
hence c5 = 0.

Verify: every strongly connected tournament on 5 vertices has c3 >= 3.

kind-pasteur-2026-03-06
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '03-artifacts', 'code'))
from tournament_lib import tournament_from_bits
from itertools import combinations

def c3_from_scores(scores):
    n = len(scores)
    return (n*(n-1)*(n-2)//6) - sum(s*(s-1)//2 for s in scores)

def is_strongly_connected(T, n):
    """Check if tournament T on n vertices is strongly connected."""
    # BFS from 0
    visited = {0}
    queue = [0]
    while queue:
        v = queue.pop(0)
        for u in range(n):
            if u not in visited and T[v][u]:
                visited.add(u)
                queue.append(u)
    if len(visited) != n:
        return False
    # BFS on reverse graph from 0
    visited = {0}
    queue = [0]
    while queue:
        v = queue.pop(0)
        for u in range(n):
            if u not in visited and T[u][v]:
                visited.add(u)
                queue.append(u)
    return len(visited) == n

# Check: at n=5, strong connectivity => c3 >= ?
print("=== n=5: Minimum c3 among strongly connected tournaments ===")
n = 5
m = n*(n-1)//2
min_c3_sc = 100
sc_by_c3 = {}
for bits in range(1 << m):
    T = tournament_from_bits(n, bits)
    if not is_strongly_connected(T, n):
        continue
    scores = [sum(T[i]) for i in range(n)]
    c3 = c3_from_scores(scores)
    if c3 < min_c3_sc:
        min_c3_sc = c3
    sc_by_c3[c3] = sc_by_c3.get(c3, 0) + 1

print(f"  Minimum c3 among SC tournaments: {min_c3_sc}")
print(f"  c3 distribution among SC tournaments: {dict(sorted(sc_by_c3.items()))}")

# Also check at n=4
print("\n=== n=4: SC tournaments ===")
n = 4
m = n*(n-1)//2
min_c3_sc = 100
for bits in range(1 << m):
    T = tournament_from_bits(n, bits)
    if not is_strongly_connected(T, n):
        continue
    scores = [sum(T[i]) for i in range(n)]
    c3 = c3_from_scores(scores)
    if c3 < min_c3_sc:
        min_c3_sc = c3
print(f"  Minimum c3 among SC tournaments on 4 vertices: {min_c3_sc}")

# Theoretical argument: Landau's strong connectivity condition
# A tournament on n vertices is SC iff for all k=1,...,n-1:
#   sum of k smallest scores > C(k,2)
# i.e., sum_{i=1}^k s_{(i)} > k(k-1)/2

# For n=5 SC tournament: c3 = 10 - sum C(s_i, 2)
# Minimize sum C(s_i,2) subject to SC constraint and sum s_i = 10.
# SC constraint: s_{(1)} >= 1, s_{(1)}+s_{(2)} >= 2,
#   s_{(1)}+s_{(2)}+s_{(3)} >= 4, s_{(1)}+s_{(2)}+s_{(3)}+s_{(4)} >= 7.

# With integer scores 0 <= s_i <= 4, sum = 10:
# To minimize sum C(s_i,2) = maximize c3, we want scores as equal as possible.
# Regular: (2,2,2,2,2), sum C(s_i,2) = 5*1 = 5, c3 = 10-5 = 5.

# To maximize sum C(s_i,2) = minimize c3 among SC:
# We want scores as spread as possible while maintaining SC.
# The SC constraint forces partial sums > C(k,2).
# k=1: s_{(1)} >= 1
# k=2: s_{(1)} + s_{(2)} >= 2 -> both >= 1
# k=3: sum of 3 smallest >= 4. With two >= 1, third >= 2.
# k=4: sum of 4 smallest >= 7. With (1,1,2,...), 4th >= 3. So (1,1,2,3,...), 5th = 10-7 = 3.
# Score: (1,1,2,3,3). sum C(s_i,2) = 0+0+1+3+3 = 7. c3 = 10-7 = 3.

# Can we do (1,1,2,2,4)? sum = 10. Partial sums: 1,2,4,6. k=4: 6 >= 7? NO. Not SC.
# Can we do (1,1,1,3,4)? sum = 10. k=3: 1+1+1=3 >= 4? NO. Not SC.
# Can we do (1,1,2,3,3)? k=1:1>0, k=2:2>1, k=3:4>3, k=4:7>6. YES, SC.
# Can we do (1,2,2,2,3)? k=1:1>0, k=2:3>1, k=3:5>3, k=4:7>6. YES, SC.
#   sum C(s_i,2) = 0+1+1+1+3 = 6. c3 = 10-6 = 4.

# So the minimum c3 for SC tournament on 5 vertices is 3, achieved by (1,1,2,3,3).
print("\n=== Theoretical: minimum c3 for SC tournament on 5 vertices ===")
print("Score (1,1,2,3,3): c3 = 10 - (0+0+1+3+3) = 10 - 7 = 3")
print("This is the minimum because any more spread score violates Landau's SC condition.")
print("Therefore: SC on 5 vertices => c3 >= 3.")
print("Contrapositive: c3 <= 2 => no 5-vertex subset is SC => c5 = 0.")

# Now verify this extends to n=6: if total c3 <= 2, no 5-subset is SC
print("\n=== n=6: If c3(T) <= 2, is every 5-subset non-SC? ===")
n = 6
m = n*(n-1)//2
violation = False
for bits in range(1 << m):
    T = tournament_from_bits(n, bits)
    scores = [sum(T[i]) for i in range(n)]
    c3 = c3_from_scores(scores)
    if c3 > 2:
        continue
    # Check each 5-subset
    for combo in combinations(range(n), 5):
        verts = list(combo)
        # Build induced tournament
        T5 = [[0]*5 for _ in range(5)]
        for i in range(5):
            for j in range(5):
                if i != j:
                    T5[i][j] = T[verts[i]][verts[j]]
        if is_strongly_connected(T5, 5):
            sc_scores = [sum(T5[i]) for i in range(5)]
            sc_c3 = c3_from_scores(sc_scores)
            print(f"  VIOLATION: n=6 bits={bits}, c3={c3}, 5-subset {verts} is SC with c3={sc_c3}")
            violation = True
            break
    if violation:
        break

if not violation:
    print("  CONFIRMED: c3 <= 2 at n=6 => every 5-subset is non-SC.")

# The argument is purely about 5-vertex subsets:
# If the GLOBAL tournament has c3 <= 2, then EVERY 5-vertex sub-tournament
# has c3(sub) <= c3(global) <= 2 < 3. Hence no 5-subset is SC.
# Hence c5 = 0.
#
# Wait: c3 of a sub-tournament can be LESS than or equal to c3 of the global
# tournament, since the sub-tournament only counts triples within those 5 vertices.
# So c3(sub) <= c3(T). If c3(T) <= 2, then c3(sub) <= 2 < 3.
# Since SC on 5 vertices requires c3(sub) >= 3, we get: no 5-subset is SC.
# Hence c5 = 0.

print("\n=== Complete proof of Claim 1 ===")
print("Claim: For any tournament T on n >= 5 vertices, c3(T) <= 2 => c5(T) = 0.")
print()
print("Proof:")
print("  Let S be any 5-vertex subset. The induced sub-tournament T[S] has")
print("  c3(T[S]) <= c3(T) <= 2.")
print("  The minimum c3 for a strongly connected tournament on 5 vertices is 3")
print("  (achieved by score (1,1,2,3,3); any lower c3 violates Landau's SC condition).")
print("  So T[S] is not strongly connected.")
print("  A non-SC tournament has no Hamiltonian cycle.")
print("  So T[S] contributes 0 to c5.")
print("  Since S was arbitrary, c5(T) = 0. QED.")
print()
print("This argument works for ALL n >= 5, not just n <= 6!")

# Verify the claim at n=7 too
print("\n=== n=7 verification: c3 <= 2 => c5 = 0 ===")
n = 7
m = 21
total = 1 << m
checked = 0
violation = False
for bits in range(total):
    T = tournament_from_bits(n, bits)
    scores = [sum(T[i]) for i in range(n)]
    c3 = c3_from_scores(scores)
    if c3 > 2:
        continue
    checked += 1
    # Only need to check c5
    for combo in combinations(range(n), 5):
        verts = list(combo)
        T5 = [[0]*5 for _ in range(5)]
        for i in range(5):
            for j in range(5):
                if i != j:
                    T5[i][j] = T[verts[i]][verts[j]]
        sc_scores = [sum(T5[i]) for i in range(5)]
        sc_c3 = c3_from_scores(sc_scores)
        if sc_c3 >= 3:
            print(f"  Sub-c3 violation at bits={bits}")
            violation = True
            break
    if violation:
        break

print(f"  Checked {checked} tournaments with c3 <= 2 at n=7.")
if not violation:
    print("  CONFIRMED: no 5-subset has c3 >= 3.")

print("\nDone.")
