"""Probe for a linear identity among single-bin coefficients c_1 (coefficient of H(M) in the P-form, i.e. of 1^j
in the exponential polynomial of j = |TT_j| at u) of H(B), U_B(e), U_B(f), W_ef over random rooted (B, u, e, f)."""
import random, sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from tern_check import blowup
from tail_c2_multi import responses
from fractions import Fraction
def c1(vals):
    # vals[j-1] for j = 1..L ; fit sum_{l=1}^{L} c_l l^j exactly (Vandermonde), return c_1
    L = len(vals)
    import sympy as sp
    M = sp.Matrix([[sp.Integer(l)**j for l in range(1, L+1)] for j in range(1, L+1)])
    sol = M.LUsolve(sp.Matrix(vals))
    return sol[0]
rng = random.Random(7); rows = []
for t in range(40):
    n = rng.randint(5, 7)
    adj = [0]*n
    for i in range(n):
        for j in range(i+1, n):
            if rng.random() < 0.5: adj[i] |= 1 << j
            else: adj[j] |= 1 << i
    u = rng.randrange(n); e, f = rng.sample([x for x in range(n) if x != u], 2)
    H = []; Ue = []; Uf = []; W = []
    for j in range(1, 8):
        A2, n2, block, nid = blowup(adj, n, u, j)
        h, U, Wd = responses(A2, n2)
        H.append(h); Ue.append(U[nid[e]]); Uf.append(U[nid[f]]); W.append(Wd[(min(nid[e], nid[f]), max(nid[e], nid[f]))])
    rows.append((c1(H), c1(Ue), c1(Uf), c1(W)))
import sympy as sp
A = sp.Matrix([[r[1] + r[2], r[0], 1] for r in rows]); b = sp.Matrix([r[3] for r in rows])
sol = (A.T*A).LUsolve(A.T*b); res = A*sol - b
print('fit c1(W) = a(c1(Ue)+c1(Uf)) + b c1(H) + d :', [sp.nsimplify(x) for x in sol], ' max residual', max(abs(x) for x in res))
print('sample rows (c1 H, c1 Ue, c1 Uf, c1 W):', rows[:8])
