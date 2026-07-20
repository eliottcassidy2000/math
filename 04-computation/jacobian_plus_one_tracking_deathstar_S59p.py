#!/usr/bin/env python3
"""
death-star-2026-07-19-S59p (HYP-8120, THM-1325) -- where the +1 goes under
Druzkowski reduction.  Three parts, exact.

D1: G = L^{-1} o F is Yagzhev-normalized (G = X + H, H no constant/linear):
    the unit constants A(0), C(0), E0(0) ARE the linear part L (antidiagonal).
    So THM-1320's factorization c0(0) = -E0(0)A(0)C(0) is det JF(0) -- true,
    engine-verified, but CLASSICALLY TRIVIAL (Keller => linear part invertible).
    The +1 does not hide: it is the Yagzhev X.
D2: torus-blindness -- equivariant cubic-linear maps are triangular over the
    weight-0 core (THM-1320 P1); back-substitution makes the triangular layer
    injective, so ALL collisions sit in the core; in ambient dim <= 3 with a
    nontrivial torus the core is <= 2-dim = shears/triangular (P2) = injective:
    NO equivariant cubic-linear map in dim <= 3 is non-injective, and in ANY
    dimension the torus direction is collision-free.  Machine spot-checks.
D3: ONE VERIFIED STABILIZATION STEP: dim 4, u's cube argument becomes the
    AFFINE LINEAR form 1 + w; the 4th component P is solved exactly (linear
    system via row-4 multilinearity of det) so that det == const; collision
    pair transported.  Degree drops 7 -> 4.  The +1 = the affine part of ell.
"""
from fractions import Fraction as Fr
from itertools import product

def pmul(a, b):
    r = {}
    for ka, ca in a.items():
        for kb, cb in b.items():
            k = tuple(p+q for p, q in zip(ka, kb))
            v = r.get(k, 0) + ca*cb
            if v: r[k] = v
            elif k in r: del r[k]
    return r
def padd(*ps):
    r = {}
    for p in ps:
        for k, c in p.items():
            v = r.get(k, 0) + c
            if v: r[k] = v
            elif k in r: del r[k]
    return r
def psc(p, s): return {k: c*s for k, c in p.items() if c*s != 0}
def pdiff(p, i):
    r = {}
    for k, c in p.items():
        if k[i] > 0:
            k2 = list(k); k2[i] -= 1
            r[tuple(k2)] = c*k[i]
    return r
def peval(p, pt):
    s = Fr(0)
    for k, c in p.items():
        m = Fr(c)
        for xi, e in zip(pt, k):
            m *= Fr(xi)**e
        s += m
    return s

# ---------------- D1 ----------------
print("=== D1: the +1 IS the Yagzhev X ===")
NV = 3
X3 = {(1,0,0):1}; Y3 = {(0,1,0):1}; Z3 = {(0,0,1):1}; ONE3 = {(0,0,0):1}
U = padd(ONE3, pmul(X3,Y3)); U2 = pmul(U,U); U3_ = pmul(U2,U)
W3_ = padd(psc(ONE3,4), psc(pmul(X3,Y3),3))
F1 = padd(pmul(U3_, Z3), pmul(pmul(pmul(Y3,Y3),U), W3_))
F2 = padd(Y3, psc(pmul(pmul(X3,U2),Z3),3), psc(pmul(pmul(X3,pmul(Y3,Y3)),W3_),3))
F3 = padd(psc(X3,2), psc(pmul(pmul(X3,X3),Y3),-3), psc(pmul(pmul(pmul(X3,X3),X3),Z3),-1))
J0 = [[peval(pdiff(Fi, j), (0,0,0)) for j in range(3)] for Fi in (F1,F2,F3)]
print("  JF(0) =", J0, " (antidiag(A(0)=1 | C(0)=1 | E0(0)=2))")
# G = L^{-1} o F with L(x,y,z) = (z, y, 2x): L^{-1}(a,b,c) = (c/2, b, a)
G1 = psc(F3, Fr(1,2)); G2 = F2; G3 = F1
JG = [[pdiff(Gi, j) for j in range(3)] for Gi in (G1,G2,G3)]
detG = padd(
    pmul(JG[0][0], padd(pmul(JG[1][1], JG[2][2]), psc(pmul(JG[1][2], JG[2][1]), -1))),
    psc(pmul(JG[0][1], padd(pmul(JG[1][0], JG[2][2]), psc(pmul(JG[1][2], JG[2][0]), -1))), -1),
    pmul(JG[0][2], padd(pmul(JG[1][0], JG[2][1]), psc(pmul(JG[1][1], JG[2][0]), -1))))
lin_ok = all(peval(pdiff(Gi, j), (0,0,0)) == (1 if i == j else 0)
             for i, Gi in enumerate((G1,G2,G3)) for j in range(3))
const_ok = all(peval(Gi, (0,0,0)) == 0 for Gi in (G1,G2,G3))
print(f"  G = L^-1 F: G(0) = 0: {const_ok}; JG(0) = I: {lin_ok}; det JG = {detG} (expect 1)")
pts = [(0,0,Fr(-1,4)), (1,Fr(-3,2),Fr(13,2)), (-1,Fr(3,2),Fr(13,2))]
imgs = [tuple(peval(Gi, pt) for Gi in (G1,G2,G3)) for pt in pts]
print(f"  triple collision transported to G: {imgs[0] == imgs[1] == imgs[2]} (image {imgs[0]})")
print("  => the unit constants ARE L; THM-1320 P3 factorization = det JF(0): true, trivial classically.")

# ---------------- D2 ----------------
print("\n=== D2: torus-blindness spot-checks ===")
# triangular-over-shear-core sample in dim 4, weights (3,1,0,0):
# F = (x1 + c*x2^3, x2, x3 + l^3, x4 + mu^3 l^3), l = a x3 + b x4? weights: w1=3, w2=1: H1 = c x2^3 ok (3 = 3*1)
# core (x3,x4) weight 0: planar shear (P2 family: l2 = mu l1, coefficient a = -mu^3 b).
a, b, mu, c = Fr(2), Fr(1), Fr(-2), Fr(5)
av = -mu**3 * b
import random
random.seed(3)
def Fmap(p):
    x1, x2, x3, x4 = p
    l = av*x3 + b*x4
    return (x1 + c*x2**3, x2, x3 + l**3, x4 + mu**3 * l**3)
inj_fail = 0
seen = {}
for _ in range(4000):
    p = tuple(Fr(random.randint(-6,6), random.randint(1,3)) for _ in range(4))
    im = Fmap(p)
    if im in seen and seen[im] != p: inj_fail += 1
    seen[im] = p
print(f"  sample triangular-over-shear-core map: injectivity violations in 4000 random pts: {inj_fail}")
print("  (proof in THM-1325: back-substitution up the triangle; collisions must sit in the core;")
print("   in dim <= 3 the core is <= 2-dim = shears = injective => equivariant cubic-linear maps")
print("   in dim <= 3 are ALL injective; in any dim, the torus direction is collision-free.)")

# ---------------- D3 ----------------
print("\n=== D3: one verified stabilization step (dim 4, degree 7 -> 4) ===")
NV = 4
def V4(i):
    k = [0,0,0,0]; k[i] = 1
    return {tuple(k): 1}
x, y, z, w = (V4(i) for i in range(4))
ONE4 = {(0,0,0,0): 1}
uw = padd(ONE4, w)                       # 1 + w  (the cube's new AFFINE LINEAR argument)
ww = padd(psc(ONE4,4), psc(w,3))         # 4 + 3w
Ft1 = padd(pmul(pmul(pmul(uw,uw),uw), z), pmul(pmul(pmul(y,y),uw), ww))
Ft2 = padd(y, psc(pmul(pmul(x,pmul(uw,uw)),z),3), psc(pmul(pmul(x,pmul(y,y)),ww),3))
Ft3 = padd(psc(x,2), psc(pmul(pmul(x,x),y),-3), psc(pmul(pmul(pmul(x,x),x),z),-1))
# unknown 4th component P: all monomials of degree <= 4 in (x,y,z,w)
monos = [k for k in product(range(5), repeat=4) if sum(k) <= 4]
NU = len(monos)
# det is LINEAR in row 4 = grad P: det = sum_j dP/dx_j * C4j, cofactors from rows 1-3
Jt = [[pdiff(Fi, j) for j in range(4)] for Fi in (Ft1, Ft2, Ft3)]
def minor3(cols):
    (c0, c1, c2) = cols
    return padd(
        pmul(Jt[0][c0], padd(pmul(Jt[1][c1], Jt[2][c2]), psc(pmul(Jt[1][c2], Jt[2][c1]), -1))),
        psc(pmul(Jt[0][c1], padd(pmul(Jt[1][c0], Jt[2][c2]), psc(pmul(Jt[1][c2], Jt[2][c0]), -1))), -1),
        pmul(Jt[0][c2], padd(pmul(Jt[1][c0], Jt[2][c1]), psc(pmul(Jt[1][c1], Jt[2][c0]), -1))))
cof = []
for j in range(4):
    cols = [cc for cc in range(4) if cc != j]
    sign = (-1)**(4+1+j+1)   # cofactor C_{4j}: (-1)^{4+j} with 1-indexing
    cof.append(psc(minor3(cols), (-1)**(4+ (j+1))))
# det(P) = sum_j (dP/dx_j) * cof[j]  -- build rows: for each unknown monomial m,
# contribution sum_j d(m)/dx_j * cof[j]; RHS: det must equal const (free)
rows = {}
for ui, m in enumerate(monos):
    pm_ = {m: 1}
    contrib = padd(*[pmul(pdiff(pm_, j), cof[j]) for j in range(4)])
    for k, cc in contrib.items():
        rows.setdefault(k, [0]*(NU))
        rows[k][ui] = rows[k].get if False else rows[k][ui] + cc
M = []; keys = []
for k, rr in rows.items():
    if k == (0,0,0,0): continue
    M.append([Fr(v) for v in rr]); keys.append(k)
# solve homogeneous M u = 0, then pick u with: nonzero det-const (contribution at (0,0,0,0)),
# gradient nondegenerate, collision-4th-coords equal.
aug = [row[:] for row in M]
r = 0; piv = []
for cc in range(NU):
    pv = next((i for i in range(r, len(aug)) if aug[i][cc] != 0), None)
    if pv is None: continue
    aug[r], aug[pv] = aug[pv], aug[r]
    inv = aug[r][cc]
    aug[r] = [xx/inv for xx in aug[r]]
    for i in range(len(aug)):
        if i != r and aug[i][cc] != 0:
            f = aug[i][cc]
            aug[i] = [aa - f*bb for aa, bb in zip(aug[i], aug[r])]
    r += 1; piv.append(cc)
free = [cc for cc in range(NU) if cc not in piv]
print(f"  P-solution space: {NU} unknowns, rank {r}, dim {len(free)}")
basis = []
for fc in free:
    v = [Fr(0)]*NU; v[fc] = Fr(1)
    for i, cc in enumerate(piv): v[cc] = -aug[i][fc]
    basis.append(v)
# const-det functional: contribution of each monomial to det at key (0,0,0,0)
const_row = rows.get((0,0,0,0), [0]*NU)
# find basis members with nonzero const: 
good = [(vb, sum(Fr(const_row[i])*vb[i] for i in range(NU))) for vb in basis]
good = [(vb, cv) for vb, cv in good if cv != 0]
print(f"  directions with nonzero det-const: {len(good)}")
# among solutions prefer SIMPLE: check the natural candidate P = w (linear):
w_idx = monos.index((0,0,0,1))
in_space = None
for vb, cv in good:
    supp = [monos[i] for i in range(NU) if vb[i] != 0]
    if supp == [(0,0,0,1)]:
        in_space = (vb, cv); break
if in_space:
    print("  P = w works: det =", in_space[1])
    Psol = {(0,0,0,1): 1}
else:
    vb, cv = good[0]
    Psol = {monos[i]: vb[i] for i in range(NU) if vb[i] != 0}
    print(f"  simplest found direction: {len(Psol)} terms, det-const {cv}; support sample:",
          dict(list(Psol.items())[:6]))
# verify chosen P: det == const exactly
Ft4 = Psol
Jt4 = [pdiff(Ft4, j) for j in range(4)]
det4 = padd(*[pmul(Jt4[j], cof[j]) for j in range(4)])
const = det4.get((0,0,0,0), 0)
print(f"  VERIFY: det J Ftilde = {dict(det4) if len(det4) < 4 else 'nonconstant!'} (constant {const})")
# collisions: transported pair on the slice w = xy at the orbit points (w* = -3/2)
Pt = (1, Fr(-3,2), Fr(13,2), Fr(-3,2))
Qt = (-1, Fr(3,2), Fr(13,2), Fr(-3,2))
imP = tuple(peval(Fi, Pt) for Fi in (Ft1, Ft2, Ft3, Ft4))
imQ = tuple(peval(Fi, Qt) for Fi in (Ft1, Ft2, Ft3, Ft4))
print(f"  Ftilde(P~) = {imP}")
print(f"  Ftilde(Q~) = {imQ}")
print(f"  COLLISION IN DIM 4: {imP == imQ}")
deg = max(sum(k) for Fi in (Ft1,Ft2,Ft3,Ft4) for k in Fi)
print(f"  degree of Ftilde: {deg} (witness was 7)")
print("  => the u-cube's argument is now the AFFINE LINEAR form 1 + w: after one")
print("     stabilization the +1 sits as the affine part of ell -- exactly what")
print("     Druzkowski's homogeneous-linear ell forbids; translating w' = 1 + w")
print("     moves it to a base-point constant (F(0) != 0), never to zero.")
