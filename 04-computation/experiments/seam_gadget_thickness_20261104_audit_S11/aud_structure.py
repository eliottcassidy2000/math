"""aud_structure.py -- audit S11 of the note's section 4 (b) and (c), with a COMPLETE rational eigenvalue search.

(c) For every distinct gadget matrix of order <= 5 (audit files aud_mats1..5.txt, 264 distinct; the note says 281,
    which is the sum of the per-order counts with cross-order repeats), find ALL integer eigenvalues (rational
    eigenvalues of an integer matrix are integers; for det = 0 the quadratic cofactor is searched too), then for every
    pair of distinct matrices and nonzero eigenvalues compute the exact intersection of the eigenspaces (right and left).
    Also reports what the audited search would have seen (it tried only divisors of det, i.e. only lambda = 0 when
    det = 0, and compared only kernel BASIS vectors).
(b) Words of length <= 3 over the distinct matrices of order <= 4: characteristic polynomials with |det| = 1 and
    |tr|, |c2| <= 3 (the audited filter), plus all words with eigenvalue 1 (p(1) = 0) and with (t-1)^2 | p.
Usage: python3 aud_structure.py DIR"""
import sys, itertools
from fractions import Fraction as F
from math import gcd
from collections import defaultdict
D = sys.argv[1]

def load(ks):
    seen = {}
    for k in ks:
        for line in open(f'{D}/aud_mats{k}.txt'):
            a = tuple(map(int, line.split()))
            seen.setdefault(a, k)
    return [((a[0:3], a[3:6], a[6:9]), k) for a, k in seen.items()]

def det3(m):
    return (m[0][0]*(m[1][1]*m[2][2]-m[1][2]*m[2][1]) - m[0][1]*(m[1][0]*m[2][2]-m[1][2]*m[2][0])
            + m[0][2]*(m[1][0]*m[2][1]-m[1][1]*m[2][0]))
def charpoly(m):
    tr = m[0][0] + m[1][1] + m[2][2]
    c2 = (m[0][0]*m[1][1]-m[0][1]*m[1][0]) + (m[0][0]*m[2][2]-m[0][2]*m[2][0]) + (m[1][1]*m[2][2]-m[1][2]*m[2][1])
    return tr, c2, det3(m)
def divisors(n):
    n = abs(n); return [d for d in range(1, n + 1) if n % d == 0]
def int_roots(tr, c2, dt):
    p = lambda t: t**3 - tr*t*t + c2*t - dt
    roots = set()
    if dt != 0:
        for d in divisors(dt):
            for s in (d, -d):
                if p(s) == 0: roots.add(s)
    else:
        roots.add(0)
        # t^2 - tr t + c2
        q = lambda t: t*t - tr*t + c2
        cands = divisors(c2) if c2 != 0 else [abs(tr)] if tr != 0 else []
        for d in cands:
            for s in (d, -d):
                if q(s) == 0: roots.add(s)
    return roots

def rref_kernel(A):
    """exact kernel basis of a rational matrix (list of rows)"""
    m = [[F(x) for x in r] for r in A]; rows, cols = len(m), len(m[0]); piv = []; r = 0
    for c in range(cols):
        p = next((i for i in range(r, rows) if m[i][c] != 0), None)
        if p is None: continue
        m[r], m[p] = m[p], m[r]; pv = m[r][c]; m[r] = [x / pv for x in m[r]]
        for i in range(rows):
            if i != r and m[i][c] != 0:
                f = m[i][c]; m[i] = [a - f * b for a, b in zip(m[i], m[r])]
        piv.append(c); r += 1
        if r == rows: break
    basis = []
    for fc in [c for c in range(cols) if c not in piv]:
        v = [F(0)] * cols; v[fc] = F(1)
        for i, pc in enumerate(piv): v[pc] = -m[i][fc]
        basis.append(v)
    return basis
def prim(v):
    den = 1
    for x in v: den = den * x.denominator // gcd(den, x.denominator)
    w = [int(x * den) for x in v]; g = 0
    for x in w: g = gcd(g, abs(x))
    w = [x // g for x in w]; s = next(x for x in w if x != 0)
    return tuple(-x for x in w) if s < 0 else tuple(w)
def sub_lam(m, lam, left=False):
    A = [[m[i][j] - (lam if i == j else 0) for j in range(3)] for i in range(3)]
    return [[A[j][i] for j in range(3)] for i in range(3)] if left else A

G5 = load(range(1, 6))
print(f"distinct gadget matrices of order <= 5: {len(G5)}")
eig = []   # (index, lambda, right basis, left basis)
n_sing_nonzero = 0
for idx, (m, k) in enumerate(G5):
    tr, c2, dt = charpoly(m)
    for lam in int_roots(tr, c2, dt):
        rb = rref_kernel(sub_lam(m, lam)); lb = rref_kernel(sub_lam(m, lam, True))
        eig.append((idx, lam, rb, lb))
        if dt == 0 and lam != 0: n_sing_nonzero += 1
print(f"integer eigenvalues found: {len(eig)} (matrix, lambda) pairs; of which nonzero eigenvalues of SINGULAR matrices "
      f"(not examined by the audited search): {n_sing_nonzero}")
for side, pos in (("right", 2), ("left", 3)):
    shared = defaultdict(set)   # direction -> set of (idx, lam)
    planes = []
    nz = [e for e in eig if e[1] != 0]
    for (i, li, *bi), (j, lj, *bj) in itertools.combinations(nz, 2):
        if i == j: continue
        Bi = bi[pos - 2]; Bj = bj[pos - 2]
        # intersection of span(Bi) and span(Bj): solve sum a_s Bi_s = sum b_t Bj_t
        if not Bi or not Bj: continue
        cols = [list(v) for v in Bi] + [[-x for x in v] for v in Bj]
        A = [[cols[c][r] for c in range(len(cols))] for r in range(3)]
        ker = rref_kernel(A)
        if not ker: continue
        vecs = [[sum(c[s] * Bi[s][r] for s in range(len(Bi))) for r in range(3)] for c in ker]
        vecs = [v for v in vecs if any(x != 0 for x in v)]
        dim = len(rref_kernel([list(map(list, zip(*vecs)))[r] for r in range(3)])) if False else None
        if len(vecs) == 1:
            shared[prim(vecs[0])].update({(i, li), (j, lj)})
        else:
            planes.append((i, li, j, lj, len(vecs)))
    print(f"{side} eigen-directions shared by two distinct matrices with nonzero eigenvalues: {len(shared)}")
    for v, S in sorted(shared.items(), key=lambda kv: -len(kv[1])):
        lams = sorted({l for _, l in S}); mats = sorted({i for i, _ in S})
        expanding = [l for l in lams if abs(l) > 1]
        print(f"   {v}: {len(mats)} matrices, eigenvalues {lams}{'   <-- EXPANDING' if expanding else ''}")
    if planes: print(f"   pairs with a common 2-dim eigenspace: {len(planes)} (first {planes[:5]})")

# what the audited procedure could see: only lambda | det for det != 0 (lambda = 0 for det = 0), kernel basis vectors
aud = defaultdict(set)
for idx, lam, rb, lb in eig:
    m = G5[idx][0]
    if det3(m) == 0 and lam != 0: continue
    for v in rb: aud[('R', prim(v))].add((idx, lam))
    for v in lb: aud[('L', prim(v))].add((idx, lam))
print("audited-style search (det-divisor eigenvalues only, basis vectors only), shared with nonzero eigenvalue:",
      {k: sorted({l for _, l in S if l}) for k, S in aud.items() if len({i for i, l in S if l != 0}) >= 2})

# (b) short words
G4 = [m for m, k in load(range(1, 5))]
print(f"distinct gadget matrices of order <= 4: {len(G4)}")
def mm(A, B): return tuple(tuple(sum(A[i][t] * B[t][j] for t in range(3)) for j in range(3)) for i in range(3))
filt = defaultdict(int); ev1 = defaultdict(int); jord = defaultdict(int); examples = {}
for L in (1, 2, 3):
    for w in itertools.product(range(len(G4)), repeat=L):
        M = G4[w[0]]
        for i in w[1:]: M = mm(G4[i], M)
        tr, c2, dt = charpoly(M)
        if abs(dt) == 1 and abs(tr) <= 3 and abs(c2) <= 3: filt[(L, tr, c2, dt)] += 1
        p1 = 1 - tr + c2 - dt; dp1 = 3 - 2 * tr + c2
        if p1 == 0:
            ev1[(L, tr, c2, dt)] += 1
            if dp1 == 0: jord[(L, tr, c2, dt)] += 1; examples.setdefault((tr, c2, dt), (w, M))
print("words with |det| = 1, |tr|, |c2| <= 3 (audited filter): (length, tr, c2, det) -> count", dict(filt))
print("words with eigenvalue 1 (any det):", sum(ev1.values()), "; polynomials:", sorted({k[1:] for k in ev1}))
print("words with (t-1)^2 dividing the characteristic polynomial (any det):", sum(jord.values()), dict(jord))
for key, (w, M) in list(examples.items())[:5]:
    print("   example word", w, "matrix", M, "charpoly (tr,c2,det)", key)
