#!/usr/bin/env python3
"""
jacobian_weyl_endo_and_structure_boxeph_S141.py  (HYP-8080)

DEEP SESSION on the verified JC counterexample F (HYP-8065/8070/8075):
 (A) the EXPLICIT non-surjective Weyl endomorphism phi of A_3, with every Weyl
     relation verified symbolically: phi(x_i) = F_i, phi(d_j) = sum_k G_kj d_k,
     G = (JF)^{-1} = -1/2 adj(JF).  Well-definedness needs [D_i, D_j] = 0 for
     the vector fields D_j = sum_k G_kj d_k -- 9 polynomial identities, checked
     exactly.  Non-surjectivity: a preimage of x_i would give a polynomial
     left-inverse P with P(F) = x, contradicting the verified 3-point collision.
 (B) STRUCTURE: (i) z-LINEARITY  F = a(x,y) + z*b(x,y), b = (u^3, 3xu^2, -x^3),
     u = 1+xy, and 4+3xy = 1+3u (Pascal-row fingerprints);
     (ii) C*-EQUIVARIANCE  F(t^{-1}x, ty, t^2 z) = (t^2 F1, t F2, t^{-1} F3)
     (weights (-1,1,2) -> (2,1,-1); sigma/tau = the t = -1 element)  =>
     the 3-point collision is a slice of a COLLISION CURVE: for every t != 0,
     {(0,0,-t^2/4), (1/t, -3t/2, 13t^2/2), (-1/t, 3t/2, 13t^2/2)} maps to
     (-t^2/4, 0, 0) -- verified exactly at t = 2, 3, 1/2;
     (iii) the DESCENT: on invariants s = xy, r = x^2 z (source) and
     vw, uw^2 (target), F induces a PLANE map (P(s,r), Q(s,r)) -- computed.
 (C) SPORADIC vs FAMILY, first exact datum: the EQUIVARIANT KELLER TANGENT at F
     -- all weight-respecting, degree-bounded perturbations G with
     d/de det J(F + eG) = 0, i.e. tr(adj(JF) . JG) = 0 -- exact linear algebra
     over Q.  Compare with the visible trivial directions (torus, target scaling,
     source/target equivariant triangular moves).

Pure Python, exact.  boxeph-2026-07-19-S141.
"""

from fractions import Fraction as Fr
from itertools import product

# ---------- poly arithmetic over Q[x,y,z] ----------
def pmul(a, b):
    out = {}
    for (i, j, k), c in a.items():
        for (l, m, n), d in b.items():
            key = (i + l, j + m, k + n)
            out[key] = out.get(key, 0) + c * d
    return {k: v for k, v in out.items() if v}

def padd(*ps):
    out = {}
    for p in ps:
        for k, v in p.items():
            out[k] = out.get(k, 0) + v
    return {k: v for k, v in out.items() if v}

def psc(p, c):
    return {k: v * c for k, v in p.items() if v * c}

def pdiff(p, var):
    out = {}
    for e, c in p.items():
        if e[var]:
            k = list(e); k[var] -= 1
            out[tuple(k)] = out.get(tuple(k), 0) + c * e[var]
    return out

def pev(p, x, y, z):
    return sum(c * x**i * y**j * z**k for (i, j, k), c in p.items())

X, Y, Z, ONE = {(1,0,0):1}, {(0,1,0):1}, {(0,0,1):1}, {(0,0,0):1}
U = padd(ONE, pmul(X, Y))
A43 = padd(psc(ONE, 4), psc(pmul(X, Y), 3))
f1 = padd(pmul(pmul(pmul(U,U),U), Z), pmul(pmul(pmul(Y,Y),U), A43))
f2 = padd(Y, psc(pmul(pmul(X,pmul(U,U)),Z),3), psc(pmul(pmul(X,pmul(Y,Y)),A43),3))
f3 = padd(psc(X,2), psc(pmul(pmul(X,X),Y),-3), psc(pmul(pmul(pmul(X,X),X),Z),-1))
F = [f1, f2, f3]
J = [[pdiff(f, v) for v in range(3)] for f in F]

def det3(M):
    return padd(
        pmul(M[0][0], padd(pmul(M[1][1],M[2][2]), psc(pmul(M[1][2],M[2][1]),-1))),
        psc(pmul(M[0][1], padd(pmul(M[1][0],M[2][2]), psc(pmul(M[1][2],M[2][0]),-1))),-1),
        pmul(M[0][2], padd(pmul(M[1][0],M[2][1]), psc(pmul(M[1][1],M[2][0]),-1))))

assert det3(J) == {(0,0,0): -2}
print("[re-verified] det JF = -2")

# ---------------- (A) the Weyl endomorphism --------------------------------------
# adjugate: adj[i][j] = cofactor C_ji
def cof(M, i, j):
    r = [a for a in range(3) if a != i]
    c = [a for a in range(3) if a != j]
    d = padd(pmul(M[r[0]][c[0]], M[r[1]][c[1]]), psc(pmul(M[r[0]][c[1]], M[r[1]][c[0]]), -1))
    return d if (i + j) % 2 == 0 else psc(d, -1)

ADJ = [[cof(J, j, i) for j in range(3)] for i in range(3)]     # adj(J)_{ij} = C_ji
G = [[psc(ADJ[i][j], Fr(-1, 2)) for j in range(3)] for i in range(3)]   # J^{-1}
# check J.G = I
IDp = [[{(0,0,0):1} if i == j else {} for j in range(3)] for i in range(3)]
JG = [[padd(*[pmul(J[i][k], G[k][j]) for k in range(3)]) for j in range(3)] for i in range(3)]
assert all(JG[i][j] == ({(0,0,0):1} if i == j else {}) for i in range(3) for j in range(3))
print("[A] J . J^{-1} = I  (J^{-1} = -1/2 adj J, polynomial entries)")

# [D_i, D_j] = 0: for all l:  sum_k ( G_ki dk(G_lj) - G_kj dk(G_li) ) = 0
ok = True
for i in range(3):
    for j in range(i + 1, 3):
        for l in range(3):
            s = {}
            for k in range(3):
                s = padd(s, pmul(G[k][i], pdiff(G[l][j], k)),
                            psc(pmul(G[k][j], pdiff(G[l][i], k)), -1))
            if s != {}:
                ok = False
assert ok
print("[A] [D_i, D_j] = 0 for all i<j (9 identities, exact)  =>  phi IS a Weyl endomorphism:")
print("    phi(x_i) = F_i;  phi(d_j) = sum_k G_kj d_k;  [phi d_j, phi x_i] = delta_ij verified via J.G=I.")
print("    Injective (A_3 simple); NOT surjective: phi onto would give P with P(F)=x,")
print("    contradicting the verified 3-point collision.  ==> explicit proper self-embedding of A_3.")
sizes = [[len(G[i][j]) for j in range(3)] for i in range(3)]
degs = [[max((sum(e) for e in G[i][j]), default=0) for j in range(3)] for i in range(3)]
print("    G entry sizes (#monomials):", sizes, " degrees:", degs)
print("    G row 1 sample G[0][0] leading terms:", dict(list(sorted(G[0][0].items()))[:4]), "...")

# ---------------- (B)(i) z-linearity + Pascal ------------------------------------
bz = [pdiff(f, 2) for f in F]
az = [padd(f, psc(pmul(bz_i, Z), -1)) for f, bz_i in zip(F, bz)]
assert all(all(e[2] == 0 for e in b) for b in bz) and all(all(e[2] == 0 for e in a) for a in az)
U3 = pmul(pmul(U, U), U)
assert bz[0] == U3 and bz[1] == psc(pmul(X, pmul(U, U)), 3) and bz[2] == psc(pmul(pmul(X,X),X), -1)
assert A43 == padd(ONE, psc(U, 3))
print("\n[B-i] F = a(x,y) + z b(x,y); b = (u^3, 3x u^2, -x^3), u = 1+xy; and 4+3xy = 1+3u.")

# ---------------- (B)(ii) C*-equivariance + collision curve ----------------------
def weight_ok(p, w):
    return all(-e[0] + e[1] + 2 * e[2] == w for e in p)
assert weight_ok(f1, 2) and weight_ok(f2, 1) and weight_ok(f3, -1)
print("[B-ii] C*-equivariance verified: F(t^-1 x, t y, t^2 z) = (t^2 F1, t F2, t^-1 F3);")
print("       sigma/tau (S140) = the t = -1 element of this torus.")
for t in (Fr(2), Fr(3), Fr(1, 2)):
    pts = [(Fr(0), Fr(0), -t*t/4), (1/t, -3*t/2, 13*t*t/2), (-1/t, 3*t/2, 13*t*t/2)]
    imgs = [tuple(pev(f, *p) for f in F) for p in pts]
    assert imgs[0] == imgs[1] == imgs[2] == (-t*t/4, Fr(0), Fr(0)), (t, imgs)
print("       COLLISION CURVE verified at t = 2, 3, 1/2:  fiber over (-t^2/4, 0, 0) contains")
print("       (0,0,-t^2/4), (1/t, -3t/2, 13t^2/2), (-1/t, 3t/2, 13t^2/2)  -- a 1-parameter")
print("       family of triple collisions; the 3 original points are the t=1 slice.")

# ---------------- (B)(iii) the descended plane map -------------------------------
# invariants: source s = xy, r = x^2 z ; target v*w = F2*F3, u*w^2 = F1*F3^2
def to_sr(p):
    """express a weight-0 polynomial in x,y,z as a poly in s=xy, r=x^2 z (unique)."""
    out = {}
    for (i, j, k), c in p.items():
        # x^i y^j z^k = s^j r^k requires i = j + 2k
        assert i == j + 2 * k, ("not in C[s,r]!", (i, j, k))
        out[(j, k)] = out.get((j, k), 0) + c
    return out
P = to_sr(pmul(f2, f3))
Q = to_sr(pmul(f1, pmul(f3, f3)))
print("[B-iii] descended plane map on (s,r) = (xy, x^2 z):")
print("        P(s,r) = F2 F3   =", dict(sorted(P.items())))
print("        Q(s,r) = F1 F3^2 =", dict(sorted(Q.items())))
def sr_diff(p, var):
    out = {}
    for e, c in p.items():
        if e[var]:
            k = list(e); k[var] -= 1
            out[tuple(k)] = out.get(tuple(k), 0) + c * e[var]
    return out
def sr_mul(a, b):
    out = {}
    for (i, j), c in a.items():
        for (l, m), d in b.items():
            out[(i + l, j + m)] = out.get((i + l, j + m), 0) + c * d
    return {k: v for k, v in out.items() if v}
def sr_add(a, b, sgn=1):
    out = dict(a)
    for k, v in b.items():
        out[k] = out.get(k, 0) + sgn * v
    return {k: v for k, v in out.items() if v}
jac2 = sr_add(sr_mul(sr_diff(P, 0), sr_diff(Q, 1)),
              sr_mul(sr_diff(P, 1), sr_diff(Q, 0)), -1)
print("        plane Jacobian d(P,Q)/d(s,r) =", dict(sorted(jac2.items())))
print("        (NOT constant -- the descent breaks Keller-ness at the torus-fixed loci,")
print("         which is exactly why this construction lives in dim 3 and evades JC_2.)")

# ---------------- (C) equivariant Keller tangent ----------------------------------
# basis: per component i with weight w_i in (2,1,-1), monomials x^a y^b z^c with
# -a+b+2c = w_i and total degree <= deg(f_i) (7, 6, 4)
wts = [2, 1, -1]
degcap = [7, 6, 4]
basis = []
for i in range(3):
    for a in range(0, 8):
        for b in range(0, 8):
            for c in range(0, 4):
                if -a + b + 2 * c == wts[i] and a + b + c <= degcap[i]:
                    basis.append((i, (a, b, c)))
print("\n[C] equivariant tangent: basis size =", len(basis))
rows = {}          # monomial -> {col: coeff}   from tr(adj(J) . J_G) = 0
for col, (i, e) in enumerate(basis):
    mono = {e: Fr(1)}
    for v in range(3):
        d = pdiff(mono, v)
        if not d: continue
        # contribution to tr(adj . JG): sum_{v} adj[v][i]-row? tr(ADJ . JG) = sum_{a,b} ADJ[a][b] JG[b][a]
        # JG[b][a] = d(G_b)/dx_a contribution: G component b = i only
        term = pmul(ADJ[v][i], d)
        for m, cf in term.items():
            d = rows.setdefault(m, {})
            d[col] = d.get(col, 0) + cf
mat = [[Fr(r.get(c, 0)) for c in range(len(basis))] for r in rows.values()]
# rank via exact Gauss
rank = 0
nc = len(basis)
pivcols = []
for c in range(nc):
    piv = None
    for r in range(rank, len(mat)):
        if mat[r][c] != 0: piv = r; break
    if piv is None: continue
    mat[rank], mat[piv] = mat[piv], mat[rank]
    pr = mat[rank]
    inv = 1 / pr[c]
    for r in range(len(mat)):
        if r != rank and mat[r][c] != 0:
            f = mat[r][c] * inv
            for cc in range(c, nc):
                mat[r][cc] -= f * pr[cc]
    pivcols.append(c); rank += 1
dimT = len(basis) - rank
print("    rank of the linearized Keller condition =", rank, " =>  equivariant tangent dim =", dimT)
print("    visible trivial directions (equivariant): source torus (1) + target weighted scalings (3)")
print("    + F itself... compare: if dim > trivial count, GENUINE deformation directions exist.")
free = [c for c in range(nc) if c not in pivcols]
print("    free columns (candidate deformation directions), first 8:")
for c in free[:8]:
    i, e = basis[c]
    print("      component f%d  +  eps * x^%d y^%d z^%d" % (i + 1, *e))
print("DONE-CORE.")
