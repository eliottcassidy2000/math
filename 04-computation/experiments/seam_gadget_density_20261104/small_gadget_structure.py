"""Exploratory structure of small gadget matrices (the session note, section 4):
(1) rational eigen-directions shared by two or more gadget matrices of order <= 5 (right and left);
(2) words of length <= 3 over the gadget matrices of order <= 4 whose characteristic polynomial has |det| = 1 and
    |trace|, |c2| <= 3 (the only candidates for unipotent or finite-order elements).
Reads mats1.txt ... mats5.txt (gadget_matrices.c output)."""
import itertools
from fractions import Fraction as F
from math import gcd
from collections import defaultdict
def load(ks):
    out = []
    for k in ks:
        for line in open(f'mats{k}.txt'):
            a = list(map(int, line.split())); out.append((k, [a[0:3], a[3:6], a[6:9]]))
    return out
def det3(m):
    return (m[0][0]*(m[1][1]*m[2][2]-m[1][2]*m[2][1]) - m[0][1]*(m[1][0]*m[2][2]-m[1][2]*m[2][0]) + m[0][2]*(m[1][0]*m[2][1]-m[1][1]*m[2][0]))
def kernel(A):
    m = [[F(x) for x in r] for r in A]; piv = []; r = 0
    for c in range(3):
        p = next((i for i in range(r, 3) if m[i][c] != 0), None)
        if p is None: continue
        m[r], m[p] = m[p], m[r]; pv = m[r][c]; m[r] = [x / pv for x in m[r]]
        for i in range(3):
            if i != r and m[i][c] != 0:
                f = m[i][c]; m[i] = [a - f * b for a, b in zip(m[i], m[r])]
        piv.append(c); r += 1
    vecs = []
    for fc in [c for c in range(3) if c not in piv]:
        v = [F(0)] * 3; v[fc] = F(1)
        for i, pc in enumerate(piv): v[pc] = -m[i][fc]
        vecs.append(v)
    return vecs
def norm(v):
    den = 1
    for x in v: den = den * x.denominator // gcd(den, x.denominator)
    w = [int(x * den) for x in v]; g = 0
    for x in w: g = gcd(g, abs(x))
    w = [x // g for x in w]
    s = next(x for x in w if x != 0)
    return tuple(-x for x in w) if s < 0 else tuple(w)
right = defaultdict(set); left = defaultdict(set)
G5 = load(range(1, 6))
for idx, (k, M) in enumerate(G5):
    D = det3(M); cands = {0} | {d for j in range(1, abs(D) + 1) if D % j == 0 for d in (j, -j)}
    for lam in cands:
        A = [[M[i][j] - (lam if i == j else 0) for j in range(3)] for i in range(3)]
        if det3(A) != 0: continue
        for v in kernel(A): right[norm(v)].add((idx, lam))
        for v in kernel([[A[j][i] for j in range(3)] for i in range(3)]): left[norm(v)].add((idx, lam))
def shared(dct):
    return {v: sorted(L) for v, L in dct.items() if len({i for i, l in L if l != 0}) >= 2}
print('order <= 5:', len(G5), 'matrices')
print('right eigen-directions shared by >= 2 matrices (nonzero eigenvalue):', {v: [l for i, l in L if l] for v, L in shared(right).items()})
print('left eigen-functionals shared by >= 2 matrices (nonzero eigenvalue):', {v: [l for i, l in L if l] for v, L in shared(left).items()})
G4 = load(range(1, 5))
def mm(A, B): return [[sum(A[i][k] * B[k][j] for k in range(3)) for j in range(3)] for i in range(3)]
found = defaultdict(int)
for L in (1, 2, 3):
    for w in itertools.product(range(len(G4)), repeat=L):
        M = G4[w[0]][1]
        for i in w[1:]: M = mm(G4[i][1], M)
        tr = M[0][0] + M[1][1] + M[2][2]
        c2 = (M[0][0]*M[1][1]-M[0][1]*M[1][0]) + (M[0][0]*M[2][2]-M[0][2]*M[2][0]) + (M[1][1]*M[2][2]-M[1][2]*M[2][1])
        D = det3(M)
        if abs(D) == 1 and abs(tr) <= 3 and abs(c2) <= 3: found[(L, tr, c2, D)] += 1
print('words over order <= 4 (', len(G4), 'matrices ), lengths 1-3, with |det| = 1, |tr|, |c2| <= 3: (length, tr, c2, det) -> count', dict(found))
print('char poly t^3 - tr t^2 + c2 t - det; (1,-1,1): t^3-t^2-t-1 (singleton, tribonacci); (3,-1,1): t^3-3t^2-t-1 (TT_2 source->sink); neither is unipotent')
