"""Wave 44 side facts.  (1) H(T'-A) for the 7-21 notes/22 window: B with its out-set on the 7 window sites (head sites
are sources of T'-A and do not matter).  (2) P_7 = QR_7: H = 189, Aut = {x -> a x + b, a in {1,2,4}} of order 21 acts
freely on Hamiltonian paths, so H = 21 * (number of orbits)."""
import itertools
es = [68, 98, 8, 8]
def ham(n, adj):
    return sum(1 for p in itertools.permutations(range(n)) if all(adj[p[i]][p[i + 1]] for i in range(n - 1)))
for c in range(4):
    n = 8; B = 7; adj = [[0] * n for _ in range(n)]
    for i in range(7):
        for j in range(i + 1, 7): adj[i][j] = 1
        if es[c] >> i & 1: adj[B][i] = 1
        else: adj[i][B] = 1
    print(f"window setting {c}: H(T'-A) = H_B = {ham(n, adj)}")
Q = {1, 2, 4}
adj = [[1 if (j - i) % 7 in Q else 0 for j in range(7)] for i in range(7)]
paths = [p for p in itertools.permutations(range(7)) if all(adj[p[i]][p[i + 1]] for i in range(6))]
aut = [(a, b) for a in (1, 2, 4) for b in range(7)]
orbits = set()
for p in paths:
    orbits.add(min(tuple((a * v + b) % 7 for v in p) for a, b in aut))
print(f"P_7: H = {len(paths)}, |Aut| = {len(aut)}, orbits of Hamiltonian paths = {len(orbits)}, paths from vertex 0 = {sum(1 for p in paths if p[0] == 0)}")
