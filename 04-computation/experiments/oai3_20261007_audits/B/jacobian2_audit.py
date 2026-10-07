import sympy as sp, numpy as np, itertools
from collections import Counter
x, y, z = sp.symbols('x y z')
u = 1 + x*y
F1 = sp.expand(u**3*z + y**2*u*(4 + 3*x*y)); F2 = sp.expand(y + 3*x*u**2*z + 3*x*y**2*(4 + 3*x*y)); F3 = sp.expand(2*x - 3*x**2*y - x**3*z)
J = sp.Matrix([F1, F2, F3]).jacobian([x, y, z])
print("det JF =", sp.factor(J.det()))
R = sp.Poly(sp.expand(F2*y**2 + F2**3 + F1**2*F3), x, y, z)
print("F2 y^2 + F2^3 + F1^2 F3: all coefficients even:", all(c % 2 == 0 for c in R.coeffs()), " (nonzero over Z:", not R.is_zero, ")")
# mod-2 field theory: x^2, y^2, z^2 in F_2(F) explicitly, and rank of JF mod 2 = 2
P2 = lambda e: sp.Poly(sp.expand(e), x, y, z, modulus=2)
# identities (cleared denominators) mod 2:
#   y^2 F2 = F2^3 + F1^2 F3 ;   x^2 (F2 + y^2 F3) = F3 ;   x^4 z^2 F3^0 ... z^2 x^6 = F3^2 + x^4 y^2
print("y^2*F2 == F2^3 + F1^2*F3 (mod 2):", P2(y**2*F2 - F2**3 - F1**2*F3).is_zero)
print("x^2*(F2 + y^2*F3) == F3 (mod 2)   [so x^2 in F_2(F,y^2)]:", P2(x**2*(F2 + y**2*F3) - F3).is_zero)
print("x^6*z^2 == F3^2 + x^4*y^2 (mod 2) [so z^2 in F_2(F,x^2,y^2)]:", P2(x**6*z**2 - F3**2 - x**4*y**2).is_zero)
minors = [J.extract(r, c).det() for r in itertools.combinations(range(3), 2) for c in itertools.combinations(range(3), 2)]
print("some 2x2 minor of JF nonzero mod 2 (rank 2 generically):", any(not P2(mm).is_zero for mm in minors))
# ---------- numeric: image of (Z/2^k)^3 ----------
def Fmod(X, Y, Z, m):
    U = (1 + X*Y) % m
    f1 = (U*U % m * U % m * Z + Y*Y % m * U % m * ((4 + 3*X*Y) % m)) % m
    f2 = (Y + 3*X % m * (U*U % m) % m * Z + 3*X % m * (Y*Y % m) % m * ((4 + 3*X*Y) % m)) % m
    f3 = (2*X - 3*X*X % m * Y - X*X % m * X % m * Z) % m
    return f1, f2, f3
for k in range(1, 8):
    m = 2**k; r = np.arange(m, dtype=np.int64)
    X, Y, Z = (a.ravel() for a in np.meshgrid(r, r, r, indexing='ij'))
    f1, f2, f3 = Fmod(X, Y, Z, m)
    code = (f1*m + f2)*m + f3
    uq, cnt = np.unique(code, return_counts=True)
    print(f"k={k}: #image/8^k = {len(uq)}/{m**3} = {sp.Rational(len(uq), m**3)}; fibre-size histogram {dict(sorted(Counter(cnt.tolist()).items()))}")
# ---------- rigorous N(w) via level-2 Hensel: F(c+4Z^3) = F(c) + 4 J(c) Z^3 ----------
Jf = sp.lambdify((x, y, z), J.tolist(), 'math')
Ff = sp.lambdify((x, y, z), [F1, F2, F3], 'math')
classes = list(itertools.product(range(4), repeat=3))
def colspace_mod2(M):  # set of vectors in F_2^3 spanned by columns
    cols = [tuple(int(M[i][j]) % 2 for i in range(3)) for j in range(3)]
    S = {(0,0,0)}
    for c in cols:
        S |= {tuple((a+b) % 2 for a, b in zip(s, c)) for s in S}
    return S
info = []
for c in classes:
    Fc = [int(t) % 8 for t in Ff(*c)]
    Jc = [[int(Jf(*c)[i][j]) for j in range(3)] for i in range(3)]
    S = colspace_mod2(Jc)
    assert len(S) == 4
    info.append((c, Fc, S))
def preimage_classes(w):  # classes C mod 4 whose image F(C) contains the 2-adic points of w mod 8
    out = []
    for c, Fc, S in info:
        d = [(w[i] - Fc[i]) % 8 for i in range(3)]
        if all(t % 4 == 0 for t in d) and tuple((t//4) % 2 for t in d) in S:
            out.append(c)
    return out
Nw = {}
for w in itertools.product(range(8), repeat=3):
    Nw[w] = preimage_classes(w)
hist = Counter(len(v) for v in Nw.values())
print("N(w) for w mod 8 (Hensel level 2):", dict(sorted(hist.items())))
img = sum(1 for v in Nw.values() if v)
print("mu(F(Z_2^3)) =", sp.Rational(img, 512), "; integral of N =", sp.Rational(sum(len(v) for v in Nw.values()), 512))
# partner law: P Haar-random; F(P) has density 2*N(w) w.r.t. Haar on target (|det|_2 = 1/2)
law = Counter()
for w, v in Nw.items():
    if v: law[len(v) - 1] += sp.Rational(2*len(v), 512)
print("P(#other preimages = j):", dict(sorted(law.items())))
# in-disc partners by residue class mod 2 of P: P mod 8 determines F(P) mod 8 and P's class mod 4
disc = {}
for P in itertools.product(range(8), repeat=3):
    w = tuple(int(t) % 8 for t in Ff(*P))
    pcs = Nw[w]; C0 = tuple(t % 4 for t in P)
    assert C0 in pcs
    same = sum(1 for c in pcs if c != C0 and all((a - b) % 2 == 0 for a, b in zip(c, C0)))
    other = sum(1 for c in pcs if any((a - b) % 2 for a, b in zip(c, C0)))
    disc.setdefault(tuple(t % 2 for t in P), Counter())[(same, other)] += 1
for k in sorted(disc):
    print("class mod 2", k, " (#in-disc partners, #other partners) over the 64 subclasses mod 8:", dict(disc[k]))
