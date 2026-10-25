"""Literal physical lines versus periodic parameters; no producer imports."""
from itertools import product
from math import comb
from pathlib import Path
import json


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def reduce(rows, width, ell):
    rows = [[v % ell for v in r] for r in rows]
    pivots = []
    for c in range(width):
        i = next((i for i in range(len(pivots), len(rows)) if rows[i][c]), None)
        if i is None:
            continue
        k = len(pivots)
        rows[k], rows[i] = rows[i], rows[k]
        inv = pow(rows[k][c], -1, ell)
        rows[k] = [v*inv % ell for v in rows[k]]
        for j in range(len(rows)):
            if j != k and rows[j][c]:
                v = rows[j][c]
                rows[j] = [(a-v*b) % ell for a,b in zip(rows[j], rows[k])]
        pivots.append(c)
    return rows[:len(pivots)], pivots


def kernel(rows, width, ell):
    red, pivots = reduce(rows, width, ell)
    basis = []
    for c in range(width):
        if c not in pivots:
            v = [0]*width
            v[c] = 1
            for row, pivot in zip(red, pivots):
                v[pivot] = -row[c] % ell
            basis.append(v)
    return basis


def physical(p, L):
    addresses = list(product(range(2), range(2), range(L)))
    points = [((-a-3*b) % p, (j-9*b) % p) for a,b,j in addresses]
    need(len(set(points)) == len(points), ('point collision', p, L))
    rows = []
    for dx,dy in ((1,0), (0,1), (3,9)):
        labels = [(x*dy-y*dx) % p for x,y in points]
        for label in sorted(set(labels)):
            rows.append([int(v == label) for v in labels])
    return addresses, rows


def periodic(p, L, ell):
    I = set(range(L))
    E = [j for j in range(p) if j in I and (j+9) % p in I]
    rows = [[1]*len(E)+[0], [0]*len(E)+[p]]
    for j in range(L, p):
        rows.append([sum(int(t == (j-s) % p) for s in (3,6,9)) for t in E]+[-1])
    parameters = kernel(rows, len(E)+1, ell)
    addresses, physical_rows = physical(p, L)

    def array(h, r, lam):
        P = [(-sum(r[(j-s) % p] for s in (3,6,9))+lam) % ell for j in range(p)]
        S = [-r[(j-9) % p] % ell for j in range(p)]
        return [(h[j], P[j]-h[j], r[j]-h[j], S[j]-P[j]+h[j])[2*a+b] % ell
                for a,b,j in addresses]

    modes = []
    for j in range(L-1):
        h = [0]*p
        h[j], h[j+1] = 1, -1
        modes.append(array(h, [0]*p, 0))
    for v in parameters:
        r = [0]*p
        for j, t in enumerate(E):
            r[t] = v[j]
        modes.append(array([0]*p, r, v[-1]))
    actual = 4*L-len(reduce(physical_rows, 4*L, ell)[1])
    need(actual == len(modes), ('dimension', p,L,ell,actual,len(modes)))
    need(len(reduce(modes, 4*L, ell)[1]) == actual, ('independence',p,L,ell))
    need(all(sum(a*b for a,b in zip(row,v)) % ell == 0
             for row in physical_rows for v in modes), ('line sum',p,L,ell))
    if L <= p-9:
        need(actual == L-1+max(0,L-10), ('inherited no-wrap',p,L,ell))
    if p >= 19 and p-L <= 8:
        d = p-L
        c = max(0,9-2*d,6-d)
        need(actual == 2*L-11+c+int(ell == p and d <= 6), ('wrap formula',p,L,ell))
        need(any(v[-1] for v in parameters) == (ell == p and d <= 6),
             ('invariant mode',p,L,ell))
    return dict(p=p,L=L,coefficient_prime=ell,J=actual,
                lambda_mode=any(v[-1] for v in parameters))


def difference_power(p, m):
    A = [[0]*p for _ in range(p)]
    for j in range(p):
        for s in range(m+1):
            A[j][(j-3*s) % p] += (-1)**s*comb(m,s)
    return A


def determinant(A):
    A = [r[:] for r in A]
    previous = 1
    sign = 1
    for j in range(len(A)-1):
        if A[j][j] == 0:
            i = next(i for i in range(j+1,len(A)) if A[i][j])
            A[i],A[j] = A[j],A[i]
            sign *= -1
        pivot = A[j][j]
        for i in range(j+1,len(A)):
            for k in range(j+1,len(A)):
                v = A[i][k]*pivot-A[i][j]*A[j][k]
                need(v % previous == 0, 'Bareiss divisibility')
                A[i][k] = v//previous
            A[i][j] = 0
        previous = pivot
    return sign*A[-1][-1]


records = []
for p in (5,7,11,13,17,19,23,29):
    lengths = range(2,p+1) if p != 29 else range(p-8,p+1)
    for L in lengths:
        for ell in sorted({2,3,5,p}):
            records.append(periodic(p,L,ell))
depth_checks = 0
cofactor_checks = 0
hasse_checks = 0
for p in (5,7,11,19):
    jets = [sum(comb(j,i) for j in range(i,p)) % p for i in range(p)]
    need(jets == [0]*(p-1)+[1], ('full-fibre Hasse mode',p))
    hasse_checks += 1
    for m in range(1,p+3):
        D = difference_power(p,m)
        for ell in (2,3,p):
            expected = min(m,p) if ell == p else 1
            need(p-len(reduce(D,p,ell)[1]) == expected, ('deep rank',p,m,ell))
            depth_checks += 1
        if p <= 11 and m <= p:
            need(abs(determinant([r[1:] for r in D[1:]])) == p**(m-1),
                 ('integral cofactor',p,m))
            cofactor_checks += 1
receipt = dict(status='PASS', physical_cases=len(records), records=records,
               depth_cases=depth_checks, integer_cofactors=cofactor_checks,
               full_fibre_hasse_checks=hasse_checks,
               wrap_c=[max(0,9-2*d,6-d) for d in range(9)],
               scope='literal physical incidence; coefficient field distinguished from physical F_p')
print(json.dumps(receipt,indent=2,sort_keys=True))
