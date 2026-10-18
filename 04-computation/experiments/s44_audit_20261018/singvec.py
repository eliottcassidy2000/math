#!/usr/bin/env python3
"""Singular vectors of the Virasoro vacuum module V^c at weight (p-1)(q-1) (S44 audit; independent of vir_ss.py).
Exact over Q (flint.fmpq_mat) or mod a prime P (flint.nmod_mat).  Ascending PBW basis (see vcore.py)."""
import flint
from fractions import Fraction
from vcore import parts, pos_mono_cached, central_charge


def _equations(D, c):
    basis = list(parts(D))
    eqs = {}
    for n in (1, 2):
        for j, mono in enumerate(basis):
            for mo, co in pos_mono_cached(n, mono, c):
                eqs.setdefault((n, mo), {})
                eqs[(n, mo)][j] = eqs[(n, mo)].get(j, 0) + co
    return basis, list(eqs.values())


def singular_vector_exact(p, q):
    c = central_charge(p, q)
    D = (p - 1) * (q - 1)
    N = D // 2
    basis, rows = _equations(D, c)
    M = flint.fmpq_mat(len(rows), len(basis), [flint.fmpq(Fraction(r.get(j, 0)).numerator, Fraction(r.get(j, 0)).denominator)
                                               for r in rows for j in range(len(basis))])
    # nullspace via rref of M
    R, rank = M.rref()
    piv = []
    for i in range(rank):
        for j in range(len(basis)):
            if R[i, j] != 0:
                piv.append(j)
                break
    free = [j for j in range(len(basis)) if j not in piv]
    assert len(free) == 1, ('nullity', len(free))
    fj = free[0]
    vec = [flint.fmpq(0)] * len(basis)
    vec[fj] = flint.fmpq(1)
    for i, pj in enumerate(piv):
        vec[pj] = -R[i, fj]
    top = basis.index((2,) * N)
    lead = vec[top]
    assert lead != 0, 'coefficient of L_-2^N vanishes'
    out = {}
    for j, x in enumerate(vec):
        if x != 0:
            y = x / lead
            out[basis[j]] = Fraction(int(y.p), int(y.q))
    return out


def singular_vector_modp(p, q, P):
    c = central_charge(p, q)
    D = (p - 1) * (q - 1)
    N = D // 2
    basis, rows = _equations(D, c)

    def red(x):
        x = Fraction(x)
        return x.numerator % P * pow(x.denominator % P, -1, P) % P
    M = flint.nmod_mat(len(rows), len(basis), [red(r.get(j, 0)) for r in rows for j in range(len(basis))], P)
    X, nul = M.nullspace()
    assert nul == 1, ('nullity mod P', nul)
    col = [int(X[i, 0]) for i in range(len(basis))]
    lead = col[basis.index((2,) * N)]
    assert lead, 'L_-2^N coefficient vanishes mod P'
    inv = pow(lead, -1, P)
    return {basis[i]: col[i] * inv % P for i in range(len(basis)) if col[i]}


if __name__ == '__main__':
    import sys
    v = singular_vector_exact(3, 4)
    print('Ising singular vector (ascending order):', {k: str(x) for k, x in v.items()})
