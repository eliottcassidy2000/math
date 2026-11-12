#!/usr/bin/env python3
"""Lane "digits": minimal linear representation of a transfer system and its exact affine classes.

A transfer system: letters l with integer matrices T_l (D x D), initial column v0, final row phi;
H(word l_1 ... l_n) = phi T_(l_n) ... T_(l_1) v0.
Minimal representation (Schuetzenberger): d = rank of the Hankel pairing.  We compute
  R  = span of the reachable columns  T_w v0           (basis by Krylov),
  Ob = span of the observable rows    phi T_w  restricted to R,
and coordinates  x(s) = (rho_1 s, ..., rho_d s)  with rho_1 = phi and rho_k = phi T_(u_k) for explicit words u_k.
In these coordinates the letters act by rational d x d matrices C_l, x0 = x(v0), readout e_1.
"""
from fractions import Fraction

from affine import RowSpace, matmul, matvec, rowmat


def transpose(M):
    return [list(r) for r in zip(*M)]


def krylov_basis(seed, mats, D, left=False):
    """independent vectors spanning the smallest space containing seed and closed under the mats;
    left = False: columns v -> M v; left = True: rows r -> r M.  Returns (vectors, words)."""
    V = RowSpace(D)
    vecs = []
    words = []
    queue = []
    if V.add(seed):
        vecs.append([Fraction(x) for x in seed])
        words.append(())
        queue.append(0)
    qi = 0
    while qi < len(queue):
        idx = queue[qi]
        qi += 1
        v = vecs[idx]
        for li, M in enumerate(mats):
            w = rowmat(v, M) if left else matvec(M, v)
            if V.add(w):
                vecs.append(w)
                words.append(words[idx] + (li,))
                queue.append(len(vecs) - 1)
    return vecs, words


def solve_coords(basis, target):
    """coefficients c with sum c_k basis[k] = target (basis independent), exact; None if not in the span."""
    d = len(basis)
    n = len(target)
    A = [[Fraction(basis[k][i]) for k in range(d)] + [Fraction(target[i])] for i in range(n)]
    row = 0
    piv = []
    for c in range(d):
        p = next((r for r in range(row, n) if A[r][c] != 0), None)
        if p is None:
            return None
        A[row], A[p] = A[p], A[row]
        pv = A[row][c]
        A[row] = [x / pv for x in A[row]]
        for r in range(n):
            if r != row and A[r][c] != 0:
                f = A[r][c]
                A[r] = [x - f * y for x, y in zip(A[r], A[row])]
        piv.append(c)
        row += 1
    for r in range(row, n):
        if A[r][d] != 0:
            return None
    return [A[k][d] for k in range(d)]


def minimal_representation(mats, v0, phi):
    D = len(v0)
    # reachable space: columns
    Rvecs, Rwords = krylov_basis(v0, mats, D, left=False)
    dR = len(Rvecs)
    # restrict to R: coordinates of T_l on the basis Rvecs
    A = []
    for M in mats:
        cols = [solve_coords(Rvecs, matvec(M, v)) for v in Rvecs]
        A.append([[cols[j][i] for j in range(dR)] for i in range(dR)])
    a0 = solve_coords(Rvecs, v0)
    f = [sum(Fraction(phi[i]) * v[i] for i in range(D)) for v in Rvecs]
    # observable rows inside R
    Ovecs, Owords = krylov_basis(f, A, dR, left=True)
    d = len(Ovecs)
    # coordinates x = Ovecs . a  (a in R-coordinates)
    C = []
    for M in A:
        rows = []
        for rho in Ovecs:
            rows.append(solve_coords(Ovecs, rowmat(rho, M)))
        C.append(rows)
    x0 = [sum(rho[i] * a0[i] for i in range(dR)) for rho in Ovecs]
    return {"dR": dR, "d": d, "C": C, "x0": x0, "obs_words": Owords}


def to_int(M):
    """common-denominator integer version (matrix, denominator)."""
    den = 1
    for row in M:
        for x in row:
            x = Fraction(x)
            g = den
            b = x.denominator
            while b:
                g, b = b, g % b
            den = den * x.denominator // g
    return [[int(Fraction(x) * den) for x in row] for row in M], den
