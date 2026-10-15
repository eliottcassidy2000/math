# Independent check of the (C) counterexample found by tbal5 at n = 9: enumerate every Hamiltonian path and every
# Hamiltonian cycle by plain DFS (no subset DP), and print the cycle-arc probabilities.
from fractions import Fraction as Fr
s = "110111111101001111111111111111111101"; n = 9
A = [[False] * n for _ in range(n)]; k = 0
for i in range(n):
    for j in range(i + 1, n):
        if s[k] == '1': A[i][j] = True
        else: A[j][i] = True
        k += 1
paths = []
def dfs(p, used):
    if len(p) == n: paths.append(tuple(p)); return
    for u in range(n):
        if not used >> u & 1 and A[p[-1]][u]: dfs(p + [u], used | 1 << u)
for v in range(n): dfs([v], 1 << v)
H = len(paths); pos = [{v: i for i, v in enumerate(P)} for P in paths]
p = lambda x, y: Fr(sum(1 for q in pos if q[x] < q[y]), H)
print("H =", H, "; scores", [sum(A[i]) for i in range(n)])
print("backward arcs (j -> i with j > i):", [(j, i) for i in range(n) for j in range(i + 1, n) if A[j][i]])
cycles = []
def dfc(c, used):
    if len(c) == n:
        if A[c[-1]][0]: cycles.append(c)
        return
    for u in range(1, n):
        if not used >> u & 1 and A[c[-1]][u]: dfc(c + [u], used | 1 << u)
dfc([0], 1)
bal = lambda q: Fr(1, 3) <= q <= Fr(2, 3)
bad = 0
for C in cycles:
    vals = [p(C[i], C[(i + 1) % n]) for i in range(n)]
    if not any(bal(q) for q in vals):
        bad += 1
        print("Ham cycle with NO balanced arc:", C)
        print("   p(c_i before c_(i+1)):", [str(q) for q in vals])
        print("   decimal:", [round(float(q), 4) for q in vals])
print("Ham cycles:", len(cycles), "; without a balanced arc:", bad)
d = max(min(p(x, y), 1 - p(x, y)) for x in range(n) for y in range(x + 1, n)); print("delta =", d, "=", round(float(d), 4))
print("vertex-wise best partner balance:", [str(max(min(p(v, w), 1 - p(v, w)) for w in range(n) if w != v)) for v in range(n)])
