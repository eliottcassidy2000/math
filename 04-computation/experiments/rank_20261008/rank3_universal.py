#!/usr/bin/env python3
"""Universal rank-3 families for translation-only maps with three independent non-unit positions on Z_d (d prime).
Lag-b covariance (scale-free part):  K = 2I - A, A = adjacency of the lag-b cycle restricted to the 3 positions.
  AP type (positions p, p+c, p+2c):  K in { 2I - A(P_3) [lags +-c], 2I - E_13 [lags +-2c], 2I [other lags] }
  non-AP type (pairwise differences distinct up to sign): K in { 2I - E_12, 2I - E_23, 2I - E_13, 2I }
(for d = 5 every 3-set is an AP; for d = 7 the non-AP type has no 2I lag; '2I' occurs iff d - 1 > #edge lags).
Balanced:  (1/2) tr(K Q^-1) Q - K > 0  for every K in the family.  Verify explicit integer Q exactly, and check
that Q = I fails (equality) so that a non-standard form is genuinely needed.  Also: rank >= 4 standard-form check."""
from fractions import Fraction as Fr
from balanced_lamperti import matinv, posdef_exact
import itertools, random, math
def K_from_edges(n, edges):
    K = [[Fr(2) if i == j else Fr(0) for j in range(n)] for i in range(n)]
    for i, j in edges: K[i][j] -= 1; K[j][i] -= 1
    return K
def balanced(Q, Ks):
    if not posdef_exact(Q): return False
    Qi = matinv(Q)
    for K in Ks:
        n = len(K); t = sum(K[a][c] * Qi[c][a] for a in range(n) for c in range(n))
        S = [[t / 2 * Q[a][c] - K[a][c] for c in range(n)] for a in range(n)]
        if not posdef_exact(S): return False
    return True
I3 = [[Fr(int(i == j)) for j in range(3)] for i in range(3)]
AP = [K_from_edges(3, [(0, 1), (1, 2)]), K_from_edges(3, [(0, 2)]), K_from_edges(3, [])]     # ends 0, 2; middle 1
NAP = [K_from_edges(3, [(0, 1)]), K_from_edges(3, [(1, 2)]), K_from_edges(3, [(0, 2)]), K_from_edges(3, [])]
NAP7 = NAP[:3]                                                                                 # d = 7: no edgeless lag
for name, fam in (('AP', AP), ('AP without 2I (d = 5)', AP[:2]), ('non-AP (d >= 11)', NAP), ('non-AP (d = 7)', NAP7)):
    print(f"{name}: Q = I balanced? {balanced(I3, fam)}")
    found = []
    rnd = random.Random(4)
    for _ in range(200000):
        a, b, c = (rnd.randint(1, 9) for _ in range(3)); x, y, z = (rnd.randint(-5, 5) for _ in range(3))
        Q = [[Fr(a), Fr(x), Fr(y)], [Fr(x), Fr(b), Fr(z)], [Fr(y), Fr(z), Fr(c)]]
        if balanced(Q, fam): found.append((a + b + c, [[int(v) for v in row] for row in Q]))
    found.sort()
    print(f"   integer balanced forms found: {len(found)}; smallest-trace examples: {[f[1] for f in found[:3]]}")
# rank >= 4: standard form works for every union-of-paths family (prime d): check all path partitions of rho <= 8 vertices
ok = True
for rho in range(4, 9):
    # every induced subgraph of a cycle on rho vertices after deleting >= 1 vertex is a disjoint union of paths;
    # worst case lambda_max(2I - A) is the longest path P_rho: 2 + 2 cos(pi/(rho+1)) < 4 <= rho = tr/2 of K
    lam = 2 + 2 * math.cos(math.pi / (rho + 1)); ok &= lam < rho
    print(f"rank {rho}: max eigenvalue of 2I - A(paths) <= {lam:.4f} < tr/2 = {rho}  -> standard form balanced: {lam < rho}")
print("rank >= 4 standard-form claim:", ok)
