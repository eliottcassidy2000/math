"""Exploratory structure of small gadget matrices (note section 4 (b)-(c)); corrected after audit S11 (M2).

(1) Rational eigen-directions shared by two or more DISTINCT gadget matrices of order <= 5 with nonzero eigenvalue
    (right eigenvectors and left eigen-functionals). Eigenvalues are the integer roots of the characteristic
    polynomial t^3 - tr t^2 + c2 t - det (for det = 0 the roots of the quadratic cofactor t^2 - tr t + c2 are also
    tried; the first version tried only divisors of det and so missed every nonzero eigenvalue of a singular matrix).
    Shared directions are found by intersecting eigenspaces exactly, pair by pair.
(2) Words of length <= 3 over the distinct gadget matrices of order <= 4 with |det| = 1 and |tr|, |c2| <= 3, and
    words whose characteristic polynomial has (t - 1)^2 as a factor (a double eigenvalue 1).
Reads mats1.txt ... mats5.txt (gadget_matrices.c output)."""
import itertools
from fractions import Fraction as F
from math import gcd
from collections import defaultdict
def load(ks):
    out = []
    for k in ks:
        for line in open(f'mats{k}.txt'):
            a = tuple(map(int, line.split())); out.append(a)
    return sorted(set(out))                      # distinct matrices only
def M3(a): return [list(a[0:3]), list(a[3:6]), list(a[6:9])]
def det3(m):
    return (m[0][0]*(m[1][1]*m[2][2]-m[1][2]*m[2][1]) - m[0][1]*(m[1][0]*m[2][2]-m[1][2]*m[2][0]) + m[0][2]*(m[1][0]*m[2][1]-m[1][1]*m[2][0]))
def cp(m):
    tr = m[0][0]+m[1][1]+m[2][2]
    c2 = (m[0][0]*m[1][1]-m[0][1]*m[1][0]) + (m[0][0]*m[2][2]-m[0][2]*m[2][0]) + (m[1][1]*m[2][2]-m[1][2]*m[2][1])
    return tr, c2, det3(m)
def int_roots(tr, c2, d):
    f = lambda t: t**3 - tr*t**2 + c2*t - d
    cands = set()
    if d != 0:
        for j in range(1, abs(d) + 1):
            if d % j == 0: cands |= {j, -j}
    else:
        cands.add(0)
        if c2 != 0:
            for j in range(1, abs(c2) + 1):
                if c2 % j == 0: cands |= {j, -j}
        else:
            cands.add(tr)
    return sorted(t for t in cands if f(t) == 0)
def nullspace(A):
    m = [[F(x) for x in r] for r in A]; piv = []; r = 0; n = len(A[0])
    for c in range(n):
        p = next((i for i in range(r, len(m)) if m[i][c] != 0), None)
        if p is None: continue
        m[r], m[p] = m[p], m[r]; pv = m[r][c]; m[r] = [x / pv for x in m[r]]
        for i in range(len(m)):
            if i != r and m[i][c] != 0:
                f = m[i][c]; m[i] = [a - f * b for a, b in zip(m[i], m[r])]
        piv.append(c); r += 1
    vecs = []
    for fc in [c for c in range(n) if c not in piv]:
        v = [F(0)] * n; v[fc] = F(1)
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
def eigenspaces(m, left=False):
    A0 = [[m[j][i] for j in range(3)] for i in range(3)] if left else m
    tr, c2, d = cp(m); out = []
    for lam in int_roots(tr, c2, d):
        if lam == 0: continue
        B = [[A0[i][j] - (lam if i == j else 0) for j in range(3)] for i in range(3)]
        out.append((lam, nullspace(B)))
    return out
def intersect(B1, B2):
    # basis of span(B1) ∩ span(B2): solve sum x_i B1_i = sum y_j B2_j
    cols = [list(v) for v in B1] + [[-x for x in v] for v in B2]
    A = [[cols[c][r] for c in range(len(cols))] for r in range(3)]
    sol = nullspace(A); out = []
    for s in sol:
        v = [sum(s[i] * B1[i][r] for i in range(len(B1))) for r in range(3)]
        if any(x != 0 for x in v): out.append(v)
    return out
G5 = load(range(1, 6))
print('order <= 5: distinct matrices', len(G5))
for side in ('right', 'left'):
    E = [(idx, lam, B) for idx, a in enumerate(G5) for lam, B in eigenspaces(M3(a), left=(side == 'left'))]
    shared = defaultdict(set)
    for (i, l1, B1), (j, l2, B2) in itertools.combinations(E, 2):
        if i == j: continue
        I = intersect(B1, B2)
        if len(I) == 1:
            v = norm(I[0]); shared[v] |= {(i, l1), (j, l2)}
        elif len(I) >= 2:                       # whole plane shared: record its coordinate-basis directions
            for w in I: v = norm(w); shared[v] |= {(i, l1), (j, l2)}
    rows = []
    for v, S in shared.items():
        mats = {i for i, l in S}; lams = sorted({l for i, l in S})
        rows.append((v, len(mats), lams, any(abs(l) > 1 for l in lams)))
    rows.sort(key=lambda r: -r[1])
    print(f'{side}: shared rational eigen-directions (nonzero eigenvalue, >= 2 distinct matrices): {len(rows)}; expanding: {sum(r[3] for r in rows)}')
    for v, nm, lams, exp in rows: print('   ', v, 'matrices', nm, 'eigenvalues', lams)
G4 = load(range(1, 5))
def mm(A, B): return [[sum(A[i][k] * B[k][j] for k in range(3)) for j in range(3)] for i in range(3)]
found = defaultdict(int); double1 = 0
for L in (1, 2, 3):
    for w in itertools.product(range(len(G4)), repeat=L):
        M = M3(G4[w[0]])
        for i in w[1:]: M = mm(M3(G4[i]), M)
        tr, c2, D = cp(M)
        if abs(D) == 1 and abs(tr) <= 3 and abs(c2) <= 3: found[(L, tr, c2, D)] += 1
        # (t-1)^2 | charpoly  <=>  f(1) = 0 and f'(1) = 0
        if 1 - tr + c2 - D == 0 and 3 - 2 * tr + c2 == 0: double1 += 1
print('words over the', len(G4), 'distinct matrices of order <= 4, lengths 1-3, |det| = 1, |tr|, |c2| <= 3: (length, tr, c2, det) -> count', dict(found))
print('words of length <= 3 with a double eigenvalue 1 ((t-1)^2 | charpoly):', double1)
