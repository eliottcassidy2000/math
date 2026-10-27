"""Brute-force check of MISTAKE-683 (2): S_4 => S_4 (S_4 the strong 4-tournament) is rigid and self-converse."""
from itertools import permutations
# strong 4-tournament: 0->1->2->3->0, 0->2, 1->3
S4 = {(0, 1), (1, 2), (2, 3), (3, 0), (0, 2), (1, 3)}
arcs = set()
for (a, b) in S4:
    arcs.add((a, b)); arcs.add((a + 4, b + 4))
for a in range(4):
    for b in range(4, 8):
        arcs.add((a, b))   # first copy beats second
n = 8
assert all(((i, j) in arcs) != ((j, i) in arcs) for i in range(n) for j in range(n) if i != j)
scores = [sum((i, j) in arcs for j in range(n)) for i in range(n)]
aut = anti = 0
for p in permutations(range(n)):
    if all(((p[i], p[j]) in arcs) == ((i, j) in arcs) for i in range(n) for j in range(n) if i != j):
        aut += 1
    if all(((p[j], p[i]) in arcs) == ((i, j) in arcs) for i in range(n) for j in range(n) if i != j):
        anti += 1
# is S_4 strong and unique? count strong labelled 4-tournaments up to iso by score check
print("scores", scores, "|Aut| =", aut, "#anti-automorphisms =", anti)
