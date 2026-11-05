"""aud_affine_lines.py -- audit S11, follow-up to the note's section 4 (c).
An affine line l = v0 + R e (v0 not in R e) is invariant under a linear map M iff M e = lam e and the induced map on
the quotient V / R e fixes the class of v0; M then acts on l by s -> lam s + mu (an affine digit).
For every rational eigen-direction e shared by >= 2 distinct gadget matrices of order <= 5 (complete search, as in
aud_structure.py), list the matrices with M e = lam e, lam != 0, compute their quotient maps on V / R e, and report
pairs of DISTINCT matrices with a common nonzero fixed class (a common invariant affine line not through 0), with the
eigenvalues and digits.  Usage: python3 aud_affine_lines.py DIR"""
import sys, itertools
from fractions import Fraction as F
D = sys.argv[1]
mats = {}
for k in range(1, 6):
    for line in open(f'{D}/aud_mats{k}.txt'):
        a = tuple(map(int, line.split())); mats.setdefault(a, k)
M3 = [((a[0:3], a[3:6], a[6:9]), k) for a, k in mats.items()]
def mul(m, v): return [sum(F(m[i][j]) * v[j] for j in range(3)) for i in range(3)]
dirs = [(1, 1, 1), (1, 0, 2), (1, 0, 4), (3, 1, 3), (5, 1, 5), (1, 0, 8), (2, 1, -1), (9, 1, 9), (5, 2, 6), (3, 0, -2)]
def complement_basis(e):
    # two standard vectors completing e to a basis
    for i, j in itertools.combinations(range(3), 2):
        b1 = [0, 0, 0]; b1[i] = 1; b2 = [0, 0, 0]; b2[j] = 1
        det = (e[0] * (b1[1] * b2[2] - b1[2] * b2[1]) - e[1] * (b1[0] * b2[2] - b1[2] * b2[0]) + e[2] * (b1[0] * b2[1] - b1[1] * b2[0]))
        if det != 0: return b1, b2
def coords(v, e, b1, b2):
    # solve v = x e + y b1 + z b2 exactly (Cramer)
    A = [[F(e[r]), F(b1[r]), F(b2[r])] for r in range(3)]
    def det(m): return m[0][0]*(m[1][1]*m[2][2]-m[1][2]*m[2][1]) - m[0][1]*(m[1][0]*m[2][2]-m[1][2]*m[2][0]) + m[0][2]*(m[1][0]*m[2][1]-m[1][1]*m[2][0])
    d = det(A); out = []
    for c in range(3):
        B = [row[:] for row in A]
        for r in range(3): B[r][c] = F(v[r])
        out.append(det(B) / d)
    return out
total_pairs = 0
for e in dirs:
    b1, b2 = complement_basis(e)
    fam = []
    for m, k in M3:
        w = mul(m, e)
        lam = w[0] / e[0] if e[0] else (w[1] / e[1] if e[1] else w[2] / e[2])
        if lam != 0 and all(w[i] == lam * e[i] for i in range(3)):
            # quotient map on V / R e in basis (b1, b2): columns = coords of M b1, M b2 dropping the e-part
            c1 = coords(mul(m, b1), e, b1, b2); c2 = coords(mul(m, b2), e, b1, b2)
            Q = [[c1[1], c2[1]], [c1[2], c2[2]]]
            fam.append((m, k, lam, Q))
    # fixed vectors of each Q (eigenvalue 1); common nonzero fixed vectors of two distinct matrices
    pairs = []
    for (m1, k1, l1, Q1), (m2, k2, l2, Q2) in itertools.combinations(fam, 2):
        A = [[Q1[0][0] - 1, Q1[0][1]], [Q1[1][0], Q1[1][1] - 1], [Q2[0][0] - 1, Q2[0][1]], [Q2[1][0], Q2[1][1] - 1]]
        # nonzero x with A x = 0 iff rank(A) < 2
        r = 0
        rows = [row[:] for row in A]
        for c in range(2):
            p = next((i for i in range(r, 4) if rows[i][c] != 0), None)
            if p is None: continue
            rows[r], rows[p] = rows[p], rows[r]
            for i in range(4):
                if i != r and rows[i][c] != 0:
                    f = rows[i][c] / rows[r][c]; rows[i] = [a - f * b for a, b in zip(rows[i], rows[r])]
            r += 1
        if r < 2: pairs.append((k1, l1, k2, l2))
    total_pairs += len(pairs)
    print(f"direction {e}: {len(fam)} matrices with M e = lam e (lam != 0), eigenvalues {sorted({f[2] for f in fam})}, "
          f"rank-1 among them {sum(1 for f in fam if all(x == 0 for row in f[3] for x in row))}; "
          f"pairs with a common invariant affine line not through 0: {len(pairs)} {pairs[:6]}")
print("total pairs of distinct order<=5 gadgets sharing an invariant affine line not through the origin:", total_pairs)
