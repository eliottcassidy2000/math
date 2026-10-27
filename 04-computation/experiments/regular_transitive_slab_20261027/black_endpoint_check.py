"""Check of the black endpoint lemma (THM-4819 (D), Lemma B) in its strong form: for a regular tournament T on 2m+1 >= 7 vertices,
a vertex v and ANY Hamilton cycle C of T - v, at least one of the m paths (succ_C(u), ..., u, v), u in N^-(v), is black
(its positional reversal is not an anti-automorphism).  Instances: all circulants and random regular tournaments
(3-cycle reversal walks) of orders 7..15; every vertex; several random Hamilton cycles of T - v per vertex.
Also counts, per instance, how many of the m rotation paths are blue (the lemma says fewer than m)."""
import random, sys
from itertools import product
def circulant(n, S): return [[1 if (j - i) % n in S else 0 for j in range(n)] for i in range(n)]
def random_regular(n, rng):
    m = (n - 1) // 2; A = circulant(n, set(range(1, m + 1)))
    for _ in range(40 * n * n):
        a, b, c = rng.sample(range(n), 3)
        if A[a][b] and A[b][c] and A[c][a]:
            A[a][b] = A[b][c] = A[c][a] = 0; A[b][a] = A[c][b] = A[a][c] = 1
    return A
def ham_cycle(A, verts, rng, tries=2000):
    """random Hamilton cycle of the sub-tournament on verts (randomised DFS)"""
    verts = list(verts); k = len(verts)
    for _ in range(tries):
        s = rng.choice(verts); path = [s]; used = {s}
        def dfs():
            if len(path) == k: return A[path[-1]][s] == 1
            nb = [w for w in verts if w not in used and A[path[-1]][w]]; rng.shuffle(nb)
            for w in nb:
                path.append(w); used.add(w)
                if dfs(): return True
                path.pop(); used.discard(w)
            return False
        if dfs(): return path
    return None
def is_blue(A, P):
    n = len(P); r = {P[i]: P[n - 1 - i] for i in range(n)}
    return all(A[r[y]][r[x]] == A[x][y] for x in range(n) for y in range(n) if x != y)
rng = random.Random(int(sys.argv[1]) if len(sys.argv) > 1 else 1)
worst = {}; inst = 0
for n in (7, 9, 11, 13, 15):
    m = (n - 1) // 2; Ts = []
    for mask in range(1 << m):
        S = {d if (mask >> (d - 1)) & 1 else n - d for d in range(1, m + 1)}; Ts.append(circulant(n, S))
    Ts += [random_regular(n, rng) for _ in range(30)]
    for A in Ts:
        for v in range(n):
            rest = [w for w in range(n) if w != v]
            for _ in range(3):
                C = ham_cycle(A, rest, rng)
                assert C is not None
                k = len(C); blue = 0
                for j in range(k):
                    u = C[j]
                    if not A[u][v]: continue
                    P = [C[(j + 1 + i) % k] for i in range(k)] + [v]
                    assert all(A[P[i]][P[i+1]] for i in range(n - 1))
                    blue += is_blue(A, P)
                inst += 1; worst[n] = max(worst.get(n, 0), blue)
                if blue == m: print("COUNTEREXAMPLE", n, v, C); sys.exit(1)
print(f"instances (tournament, vertex, Hamilton cycle of T - v): {inst}; no cycle has all m rotation paths blue; "
      f"max number of blue rotation paths per order: {worst}")
