#!/usr/bin/env python3
"""S44 audit, HYP-9316: build the Cramer (determinantal) syzygy for a monomial ideal Q = (h_1..h_K) containing the
jets F_0..F_K, check it is a nonzero Koszul 1-cycle, and test whether it is a nonzero class in H_1 at its bidegree
(i.e. not a Koszul boundary); also report dim H_1 there.  Jets F_n = [t^n] f(x(t), y(t)) (see curves_audit.py).
Usage: python cramer_audit.py <curve> <gens>   gens like x0^2,y0^1  (x_j^e or y_j^e)
"""
import sys, itertools, json, time
import flint
from curves_audit import Curve, CURVES, P


def poly_mul(p1, p2):
    out = {}
    for m1, c1 in p1.items():
        for m2, c2 in p2.items():
            k = tuple(sorted(m1 + m2))
            out[k] = (out.get(k, 0) + c1 * c2) % P
    return {k: v for k, v in out.items() if v}


def poly_add(p1, p2, s=1):
    out = dict(p1)
    for k, v in p2.items():
        out[k] = (out.get(k, 0) + s * v) % P
    return {k: v for k, v in out.items() if v}


def det(Mx):
    """Determinant of a square matrix of polynomials (Laplace along first row)."""
    n = len(Mx)
    if n == 0:
        return {(): 1}
    if n == 1:
        return Mx[0][0]
    res = {}
    for j in range(n):
        if not Mx[0][j]:
            continue
        minor = [row[:j] + row[j + 1:] for row in Mx[1:]]
        term = poly_mul(Mx[0][j], det(minor))
        res = poly_add(res, term, 1 if j % 2 == 0 else -1)
    return res


def parse_gens(s):
    gens = []
    for g in s.split(','):
        var, e = g.split('^')
        v, j = var[0], int(var[1:])
        code = 2 * j + (0 if v == 'x' else 1)
        gens.append((code, int(e)))
    return gens


def main():
    name, gs = sys.argv[1], sys.argv[2]
    a, b, terms = CURVES[name]
    gens = parse_gens(gs)
    K = len(gens)
    d = a * terms[0][1] + b * terms[0][2]
    qdeg = lambda code: a if code % 2 == 0 else b
    idx_ = lambda code: code // 2
    hq = sum(e * qdeg(c) for c, e in gens)
    hi = sum(e * idx_(c) for c, e in gens)
    Dq = (K + 1) * d - hq
    I = sum(range(K + 1)) - hi
    W = Dq + I
    C = Curve(a, b, terms, max(W, d + K) + 1)
    F = C.F
    # coefficient matrix: assign each monomial of F_l to the first generator dividing it
    A = [[{} for _ in range(K)] for _ in range(K + 1)]
    for l in range(K + 1):
        for mono, co in F[l].items():
            placed = False
            for j, (code, e) in enumerate(gens):
                if mono.count(code) >= e:
                    rest = list(mono)
                    for _ in range(e):
                        rest.remove(code)
                    key = tuple(rest)
                    A[l][j][key] = (A[l][j].get(key, 0) + co) % P
                    placed = True
                    break
            assert placed, ('F_%d not in Q' % l, mono)
    # Cramer cycle z = sum_l (-1)^l det(A without row l) e_l
    z = {}
    for l in range(K + 1):
        sub = [A[i] for i in range(K + 1) if i != l]
        dl = det(sub)
        if dl:
            z[l] = {k: (v if l % 2 == 0 else (-v) % P) for k, v in dl.items()}
    assert z, 'Cramer cycle is zero (rank < K)'
    # check cycle
    tot = {}
    for l, g in z.items():
        tot = poly_add(tot, poly_mul(g, F[l]))
    assert not tot, 'not a cycle'
    # bidegree checks
    for l, g in z.items():
        for mono in g:
            assert sum(qdeg(c) for c in mono) + d == Dq and sum(idx_(c) for c in mono) + l == I
    # Koszul complex at (Dq, I)
    B1 = [(n, m) for n in range(0, I + 1) for m in C.monos(Dq - d, I - n)]
    i1 = {x: i for i, x in enumerate(B1)}
    M0 = C.monos(Dq, I)
    i0 = {m: i for i, m in enumerate(M0)}
    rows1 = []
    for (n, m) in B1:
        row = [0] * len(M0)
        for mono, co in F[n].items():
            row[i0[tuple(sorted(mono + m))]] = (row[i0[tuple(sorted(mono + m))]] + co) % P
        rows1.append(row)
    rk1 = flint.nmod_mat(len(rows1), len(M0), [x for r in rows1 for x in r], P).rank() if rows1 and M0 else 0
    dimZ = len(B1) - rk1
    brows = []
    for n1 in range(0, I + 1):
        for n2 in range(n1 + 1, I + 1 - n1):
            for m in C.monos(Dq - 2 * d, I - n1 - n2):
                row = [0] * len(B1)
                for mono, co in F[n1].items():
                    j = i1[(n2, tuple(sorted(mono + m)))]
                    row[j] = (row[j] + co) % P
                for mono, co in F[n2].items():
                    j = i1[(n1, tuple(sorted(mono + m)))]
                    row[j] = (row[j] - co) % P
                brows.append(row)
    rkB = flint.nmod_mat(len(brows), len(B1), [x for r in brows for x in r], P).rank() if brows else 0
    zrow = [0] * len(B1)
    for l, g in z.items():
        for mono, co in g.items():
            zrow[i1[(l, mono)]] = (zrow[i1[(l, mono)]] + co) % P
    rkBz = flint.nmod_mat(len(brows) + 1, len(B1), [x for r in brows + [zrow] for x in r], P).rank()
    print('%s Q=(%s): Cramer cycle bidegree (q-deg %d, index %d), weight %d; nonzero cycle; dim K_1=%d dim Z_1=%d dim B_1=%d '
          'dim H_1=%d; Cramer cycle is a boundary: %s'
          % (name, gs, Dq, I, W, len(B1), dimZ, rkB, dimZ - rkB, rkBz == rkB))
    json.dump(dict(curve=name, gens=gs, qdeg=Dq, index=I, weight=W, dimZ=dimZ, dimB=rkB, dimH1=dimZ - rkB,
                   cramer_is_boundary=(rkBz == rkB)), open('cramer_audit_%s.json' % name, 'w'), indent=1)


if __name__ == '__main__':
    main()
