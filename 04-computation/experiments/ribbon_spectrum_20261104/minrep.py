#!/usr/bin/env python3
"""minrep.py -- minimal linear representation of the span-w ribbon series (exact, over Q):
forward reachable span, backward co-reachable span, Hankel rank; characteristic polynomials of the
letters on the reachable space (sympy on the reduced matrices only)."""
import sys
from fractions import Fraction
import sympy as sp
from ribbon import transfer_system, letters

def rref_insert(basis, v):
    """basis: list of (pivot, row) with rows reduced; insert v (list of Fractions) if independent."""
    v = v[:]
    for piv, row in basis:
        if v[piv] != 0:
            f = v[piv] / row[piv]
            v = [a - f * b for a, b in zip(v, row)]
    for i, a in enumerate(v):
        if a != 0:
            f = a
            v = [x / f for x in v]
            # reduce existing rows
            newbasis = []
            for piv, row in basis:
                if row[i] != 0:
                    g = row[i]
                    row = [x - g * y for x, y in zip(row, v)]
                newbasis.append((piv, row))
            newbasis.append((i, v))
            return newbasis, True
    return basis, False

def matvec(M, v, d):
    out = [Fraction(0)] * d
    for (j, i), m in M.items():
        if v[i] != 0:
            out[j] += m * v[i]
    return out

def vecmat(f, M, d):
    out = [Fraction(0)] * d
    for (j, i), m in M.items():
        if f[j] != 0:
            out[i] += m * f[j]
    return out

for w in [int(x) for x in (sys.argv[1:] or ['1', '2', '3'])]:
    order, states, mats, final = transfer_system(w)
    d = len(order)
    L = letters(w)
    init = [Fraction(0)] * d; init[0] = Fraction(1)
    fin = [Fraction(0)] * d
    for j in final: fin[j] = Fraction(1)
    # forward
    basis, _ = rref_insert([], init)
    fwd = [init]; frontier = [init]
    while frontier:
        nf = []
        for v in frontier:
            for a in L:
                nv = matvec(mats[a], v, d)
                basis, new = rref_insert(basis, nv)
                if new:
                    fwd.append(nv); nf.append(nv)
        frontier = nf
    fr = len(basis)
    cob, _ = rref_insert([], fin)
    bwd = [fin]; frontier = [fin]
    while frontier:
        nf = []
        for f in frontier:
            for a in L:
                nfv = vecmat(f, mats[a], d)
                cob, new = rref_insert(cob, nfv)
                if new:
                    bwd.append(nfv); nf.append(nfv)
        frontier = nf
    br = len(cob)
    # Hankel: pairings of bwd functionals with fwd vectors
    Hk = sp.Matrix([[sum(f[i] * v[i] for i in range(d)) for v in fwd] for f in bwd])
    hr = Hk.rank()
    print(f"w={w}: states {d}, forward span {fr}, backward span {br}, Hankel rank (minimal dimension) {hr}", flush=True)
    if w <= 2:
        B = sp.Matrix([[v[i] for v in fwd] for i in range(d)])  # d x fr
        rows = B.T.rref()[1]
        Bsq = B.extract(list(rows), list(range(fr)))
        Binv = Bsq.inv()
        for a in L:
            MB = sp.Matrix([[sum(mats[a].get((j, i), 0) * B[i, c] for i in range(d)) for c in range(fr)] for j in range(d)])
            A = Binv * MB.extract(list(rows), list(range(fr)))
            assert MB == B * A
            print(f"  letter {sorted(a)}: charpoly on reachable space = {sp.factor(A.charpoly().as_expr())}", flush=True)
