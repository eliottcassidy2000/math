"""Independent confirmation that (62, 0) on L_2 (Q = u^2 + 2v^2, u = v mod 2) has g = 964 = 2 floor(3844/8) + 4,
i.e. excess 4 at the prime p = 31 = 7 (mod 8), refuting the 'iff p = 5 (mod 8)' clause of HYP-9428 (3).
Lower bound: the window lemma (Q <= 4C + 7 = 3859 for C = 963) with an exhaustive window computation W = 3860.
Upper bound: an explicit tree extracted from the solver and re-verified from scratch (values, leaves, beta_+)."""
import sys
from knuth_debt import knuth_debt, INF
Qf = lambda u, v: u*u + 2*v*v
for W in (3860, 6000):
    R = int(W ** 0.5) + 2
    pts = [(u, v) for u in range(-R, R + 1) for v in range(-R, R + 1) if (u - v) % 2 == 0 and Qf(u, v) <= W]
    g, Q, par = knuth_debt(pts, [[1, 0], [0, 2]], 1, record_parent=True)
    idx = {p: k for k, p in enumerate(pts)}
    k = idx[(62, 0)]
    print(f"W = {W}: g(62,0) = {int(g[k])}, 2 floor(Q/8) = {2*(3844//8)}, certified (4g+7 <= W): {4*int(g[k])+7 <= W}")
# extract the tree from the last run and verify it independently
sys.setrecursionlimit(100000)
memo = {}
def build(k):
    if k in memo: return memo[k]
    if k not in par:
        assert Qf(*pts[k]) == 3, pts[k]; memo[k] = ('leaf', pts[k]); return memo[k]
    i, j = par[k]; memo[k] = ('m', build(i), build(j)); return memo[k]
tree = build(idx[(62, 0)])
cost_memo = {}
def check(t):
    """returns (value, cost, leaves) recomputed from the leaves; memoised on object identity of shared subtrees"""
    key = id(t)
    if key in cost_memo: return cost_memo[key]
    if t[0] == 'leaf':
        r = (t[1], 0, 1)
    else:
        (y, cy, ly), (z, cz, lz) = check(t[1]), check(t[2])
        x = (y[0] + z[0], y[1] + z[1])
        beta = Qf(*x) - Qf(*y) - Qf(*z)
        assert beta == 2*(y[0]*z[0] + 2*y[1]*z[1])
        r = (x, cy + cz + max(beta, 0), ly + lz)
    cost_memo[key] = r
    return r
val, cost, leaves = check(tree)
print(f"extracted tree: root {val}, total cost {cost}, {leaves} unit leaves (all of norm 3)")
