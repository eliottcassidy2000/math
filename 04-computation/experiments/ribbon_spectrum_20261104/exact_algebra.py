#!/usr/bin/env python3
"""exact_algebra.py -- exact (integer) verification that the letters of span 2 generate the full matrix
algebra on the 13-dimensional span of nonempty-word states (irreducibility), using python-flint."""
import itertools, numpy as np
import flint
from intrep import build
from ribbon import letters
w = 2
basis, A, y0, hcoef, order, states, mats, final = build(w, verbose=False)
r = len(basis); L = letters(w)
Am = {a: flint.fmpz_mat(A[a]) for a in L}
# integer basis of the span of nonempty-word states: states of words of length 1..5, echelon over Z via
# flint hnf of the stacked vectors (rows)
y = flint.fmpz_mat([[v] for v in y0])
vecs = []
frontier = [y]
for n in range(1, 6):
    nf = []
    for v in frontier:
        for a in L:
            nv = Am[a] * v
            nf.append(nv)
            vecs.append([int(nv[i, 0]) for i in range(r)])
    frontier = nf
S = flint.fmpz_mat(vecs)
rk = S.rank()
print("rank of the span of nonempty-word states:", rk)
# pick an integer basis: rows of the HNF
Hn = S.hnf()
rows = [[int(Hn[i, j]) for j in range(r)] for i in range(Hn.nrows()) if any(int(Hn[i, j]) for j in range(r))]
assert len(rows) == rk
Sb = flint.fmpz_mat(rows).transpose()   # r x rk, columns = basis vectors
# restricted matrices over Q: B_a = solve S B = A S  (S has full column rank); use fmpq
Sq = flint.fmpq_mat(Sb)
def restrict(M):
    AS = flint.fmpq_mat(M) * Sq
    # least squares exact: (S^T S)^{-1} S^T (A S)
    StS = Sq.transpose() * Sq
    B = StS.inv() * (Sq.transpose() * AS)
    assert Sq * B == AS, "not invariant"
    return B
Bq = {a: restrict(Am[a]) for a in L}
k = rk
def algebra_rank(gens):
    vecs = [flint.fmpq_mat(k, k, [1 if i == j else 0 for i in range(k) for j in range(k)])]
    flat = [[int(x.p) if x.q == 1 else None for x in [v[i, j] for i in range(k) for j in range(k)]] for v in vecs]
    frontier = list(vecs)
    cur = flint.fmpq_mat([[vecs[0][i, j] for i in range(k) for j in range(k)]])
    while frontier:
        nf = []
        for X in frontier:
            for G in gens:
                Y = G * X
                cand = flint.fmpq_mat([[cur[i, j] for j in range(k * k)] for i in range(cur.nrows())] + [[Y[i, j] for i in range(k) for j in range(k)]])
                if cand.rank() > cur.rank():
                    cur = cand; nf.append(Y)
        frontier = nf
    return cur.rank()
print("exact algebra dimension, all four letters:", algebra_rank([Bq[a] for a in L]), "of", k * k)
strong = [frozenset([1]), frozenset([2]), frozenset([1, 2])]
print("exact algebra dimension, three nonempty letters:", algebra_rank([Bq[a] for a in strong]), "of", k * k)
