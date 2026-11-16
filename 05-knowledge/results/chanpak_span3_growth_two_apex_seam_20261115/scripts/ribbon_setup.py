#!/usr/bin/env python3
"""ribbon_setup.py: the literal span-3 transfer (119 profile states, THM-4880) via ribbon.transfer_system(3);
exports the eight letter matrices (mats3.txt), computes the reachable span and its annihilator, and solves the
annihilator LP for the bounded-growth certificate  h A_l - psi_l <= K h  (psi_l in the annihilator of the reachable span)."""
import sys, os, numpy as np
sys.path.insert(0, '/tmp/math-wt-ribbon-20261115/04-computation/experiments/ribbon_spectrum_20261104')
from ribbon import transfer_system
order, states, mats, final = transfer_system(3)
N = len(order)
def lt(b): return frozenset(g for g in (1, 2, 3) if (b >> (g - 1)) & 1)
A = []
for b in range(8):
    M = np.zeros((N, N), dtype=np.int64)
    for (j, i), m in mats[lt(b)].items(): M[j, i] += m
    A.append(M)
e0 = np.zeros(N, dtype=np.int64); e0[0] = 1
v0 = A[0] @ e0
h = np.zeros(N, dtype=np.int64)
for j in final: h[j] = 1
with open('mats3.txt', 'w') as f:
    f.write(f"{N}\n")
    for b in range(8):
        for i in range(N): f.write(" ".join(map(str, A[b][i])) + "\n")
    f.write(" ".join(map(str, v0)) + "\n" + " ".join(map(str, h)) + "\n")
print("states", N, "complete states (h=1):", int(h.sum()), "v0 support", int((v0 != 0).sum()))
# sanity: H of a few words by matrices vs brute force
from ribbon import tournament_from_word, ham_count_bruteforce
import random
rng = random.Random(1)
for trial in range(5):
    n = rng.randint(4, 9); word = [lt(rng.randrange(8)) for _ in range(n)]
    adj = tournament_from_word(3, word); Hb = ham_count_bruteforce(adj)
    x = e0.copy()
    for a in word: x = A[sum(1 << (g - 1) for g in a)] @ x
    assert int(h @ x) == Hb, (word, int(h @ x), Hb)
print("matrix counts agree with brute force on 5 random words")
# reachable span (BFS, exact rank over Q via fractions-free modular test)
P1 = 2147483629
def rank_mod(rows, p=P1):
    M = [[int(x) % p for x in r] for r in rows]; rk = 0; cols = len(M[0]) if M else 0
    for c in range(cols):
        piv = next((i for i in range(rk, len(M)) if M[i][c]), None)
        if piv is None: continue
        M[rk], M[piv] = M[piv], M[rk]; inv = pow(M[rk][c], p - 2, p); M[rk] = [x * inv % p for x in M[rk]]
        for i in range(len(M)):
            if i != rk and M[i][c]:
                f = M[i][c]; M[i] = [(x - f * y) % p for x, y in zip(M[i], M[rk])]
        rk += 1
    return rk
B = [list(v0)]; frontier = [v0]
while frontier:
    nf = []
    for v in frontier:
        for b in range(8):
            w = A[b] @ v
            if not w.any(): continue
            if rank_mod(B + [list(w)]) > len(B): B.append(list(w)); nf.append(w)
    frontier = nf
r = len(B); print("reachable span dimension", r)
Rb = np.array(B, dtype=np.int64).T  # N x r
# annihilator: left null space of Rb (functionals f with f Rb = 0), exact via sympy
import sympy
Rs = sympy.Matrix(Rb.tolist())
ann = Rs.T.nullspace()   # vectors f (N) with Rb^T f = 0  <=> f^T Rb = 0
print("annihilator dimension", len(ann))
Ann = np.array([[float(x) for x in v] for v in ann])   # d x N
np.save('ann.npy', Ann); np.save('Rb.npy', Rb); np.save('A.npy', np.array(A)); np.save('h.npy', h); np.save('v0.npy', v0)
# LP: minimise K s.t. for each letter l and state s: (h A_l)_s - (c_l . Ann)_s <= K h_s
from scipy.optimize import linprog
d = Ann.shape[0]; nl = 8
nvar = 1 + nl * d  # K, then c_l
A_ub = []; b_ub = []
for b in range(8):
    hA = h @ A[b]
    for s in range(N):
        row = np.zeros(nvar); row[0] = -h[s]
        row[1 + b * d: 1 + (b + 1) * d] = -Ann[:, s]
        A_ub.append(row); b_ub.append(-hA[s])
c = np.zeros(nvar); c[0] = 1
res = linprog(c, A_ub=np.array(A_ub), b_ub=np.array(b_ub), bounds=[(None, None)] * nvar, method='highs')
print("annihilator LP:", res.status, res.message, "K* =", res.fun if res.status == 0 else None)
if res.status == 0:
    K = res.x[0]
    for b in range(8):
        hA = h @ A[b]; psi = res.x[1 + b * d: 1 + (b + 1) * d] @ Ann
        slack = K * h - hA + psi
        print(f"  letter {b}: min slack {slack.min():.3e}; states where (hA)_s > 0 and h_s = 0: {int(((hA > 0) & (h == 0)).sum())}")
