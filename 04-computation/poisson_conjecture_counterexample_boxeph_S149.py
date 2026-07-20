#!/usr/bin/env python3
"""
poisson_conjecture_counterexample_boxeph_S149.py  (HYP-8190; 8185 ceded to kind-pasteur first-push)

AN EXPLICIT COUNTEREXAMPLE TO THE (STABLE) POISSON CONJECTURE ON C^6.

The Poisson conjecture (Kontsevich; Adjamagbo-van den Essen chain
JC <=> DC <=> PC in the stable form) asserts: every polynomial Poisson
endomorphism of C[x1,x2,x3,p1,p2,p3] with the standard symplectic bracket
{x_i, p_j} = delta_ij is an automorphism.

CONSTRUCTION (cotangent lift of the Keller kernel map K):
  K = (u^3 z + y^2 u (4+3xy),  y + 3x u^2 z + 3x y^2 (4+3xy),  2x - 3x^2 y - x^3 z),
  u = 1 + xy, det JK = -2 (the externally supplied JC counterexample, verified
  in-repo S140).  Define Phi on C[x, p] by
       Phi(x_i) = X_i := K_i(x),
       Phi(p_j) = P_j := sum_k M_{jk}(x) p_k,   M := (JK^T)^{-1} = adj(JK)^T / (-2)
  (polynomial because det JK is a unit).  Phi is the pullback action of the
  cotangent lift T*C^3 -> T*C^3 of K.

THIS SCRIPT VERIFIES EXACTLY (pure Python, Fractions):
  (1) det JK = -2 and JK . adj(JK) = -2 I (adjugate correctness);
  (2) {X_i, X_j} = 0 (structural: no p-dependence),
      {X_i, P_j} = delta_ij            [9 polynomial identities],
      {P_i, P_j} = 0                   [3 identities, linear in p: 9 coefficient
                                        identities] -- Phi PRESERVES the bracket
      on generators, hence is a Poisson endomorphism;
  (3) det JPhi = 1 (block triangular: det JK * det M = (-2)(-1/2));
  (4) Phi is NOT injective: the collision points of K lift with p = 0:
        Phi-map sends (0,0,-1/4,0,0,0) and (1,-3/2,13/2,0,0,0) and
        (-1,3/2,13/2,0,0,0) to the SAME point ((-1/4,0,0),(0,0,0)).
  => Phi is a unimodular Poisson endomorphism that is not an automorphism:
     THE POISSON CONJECTURE IS FALSE ON C^6, EXPLICITLY.  Padding with extra
     symplectic pairs (identity) falsifies it on C^{2n} for all n >= 3: the
     STABLE Poisson conjecture is false.  This is the classical (hbar -> 0)
     shadow of the S141 explicit Dixmier endomorphism: the 12 bracket
     identities here are the symbols of S141's 9 Weyl commutator identities.

boxeph-2026-07-20-S149.
"""

from fractions import Fraction as Fr
from itertools import permutations

NV = 6   # (x1, x2, x3, p1, p2, p3)

def mono(i):
    e = [0] * NV; e[i] = 1
    return {tuple(e): Fr(1)}

def mul(a, b):
    out = {}
    for e, c in a.items():
        for f, d in b.items():
            k = tuple(i + j for i, j in zip(e, f))
            out[k] = out.get(k, 0) + c * d
    return {k: v for k, v in out.items() if v}

def add(*ps):
    out = {}
    for p in ps:
        for k, v in p.items(): out[k] = out.get(k, 0) + v
    return {k: v for k, v in out.items() if v}

def sc(p, c): return {k: v * c for k, v in p.items() if v * c}

def diff(p, var):
    out = {}
    for e, c in p.items():
        if e[var]:
            k = list(e); k[var] -= 1
            out[tuple(k)] = out.get(tuple(k), 0) + c * e[var]
    return out

def bracket(f, g):
    """standard Poisson bracket sum_i df/dx_i dg/dp_i - df/dp_i dg/dx_i."""
    out = {}
    for i in range(3):
        out = add(out, mul(diff(f, i), diff(g, i + 3)),
                       sc(mul(diff(f, i + 3), diff(g, i)), -1))
    return out

X, Y, Z = mono(0), mono(1), mono(2)
P1, P2, P3 = mono(3), mono(4), mono(5)
ONE = {(0,) * NV: Fr(1)}

U  = add(ONE, mul(X, Y))                       # 1 + xy
A43 = add(sc(ONE, 4), sc(mul(X, Y), 3))        # 4 + 3xy
K1 = add(mul(mul(mul(U, U), U), Z), mul(mul(mul(Y, Y), U), A43))
K2 = add(Y, sc(mul(mul(X, mul(U, U)), Z), 3), sc(mul(mul(X, mul(Y, Y)), A43), 3))
K3 = add(sc(X, 2), sc(mul(mul(X, X), Y), -3), sc(mul(mul(mul(X, X), X), Z), -1))
K = [K1, K2, K3]

# (1) Jacobian, determinant, adjugate
J = [[diff(K[i], j) for j in range(3)] for i in range(3)]
def det3(M):
    tot = {}
    for perm in permutations(range(3)):
        sgn = 1; pl = list(perm)
        for a in range(3):
            for b in range(a + 1, 3):
                if pl[a] > pl[b]: sgn = -sgn
        term = sc(ONE, Fr(sgn))
        for r in range(3): term = mul(term, M[r][perm[r]])
        tot = add(tot, term)
    return tot
D = det3(J)
assert D == sc(ONE, Fr(-2)), "det JK != -2"
print("(1) det JK = -2  OK")
def cof(M, r, c):
    rows = [i for i in range(3) if i != r]
    cols = [j for j in range(3) if j != c]
    m = add(mul(M[rows[0]][cols[0]], M[rows[1]][cols[1]]),
            sc(mul(M[rows[0]][cols[1]], M[rows[1]][cols[0]]), -1))
    return sc(m, Fr((-1) ** (r + c)))
ADJ = [[cof(J, c, r) for c in range(3)] for r in range(3)]   # adj = cofactor^T
for i in range(3):
    for j in range(3):
        prodij = add(*[mul(J[i][k], ADJ[k][j]) for k in range(3)])
        assert prodij == (sc(ONE, Fr(-2)) if i == j else {}), (i, j)
print("    JK . adj(JK) = -2 I  OK  (adjugate exact)")

# M = adj(JK)^T / (-2);  P_j = sum_k M_jk p_k.  Work with N = adj^T (integer).
N = [[ADJ[k][j] for k in range(3)] for j in range(3)]        # N[j][k] = adj[k][j]
Pimg = [sc(add(*[mul(N[j][k], mono(3 + k)) for k in range(3)]), Fr(-1, 2))
        for j in range(3)]
Ximg = K

# (2) bracket preservation on generators
print("(2) bracket identities:")
okXX = all(not bracket(Ximg[i], Ximg[j]) for i in range(3) for j in range(3))
print("    {X_i, X_j} = 0        :", "OK" if okXX else "FAIL")
assert okXX
okXP = True
for i in range(3):
    for j in range(3):
        br = bracket(Ximg[i], Pimg[j])
        want = ONE if i == j else {}
        if br != want: okXP = False; print("      FAIL at X_%d,P_%d" % (i, j))
print("    {X_i, P_j} = delta_ij :", "OK (9 identities)" if okXP else "FAIL")
assert okXP
okPP = True
for i in range(3):
    for j in range(i + 1, 3):
        br = bracket(Pimg[i], Pimg[j])
        if br: okPP = False; print("      FAIL at P_%d,P_%d (%d monomials)" % (i, j, len(br)))
print("    {P_i, P_j} = 0        :", "OK (3 identities, all p-coefficients)" if okPP else "FAIL")
assert okPP
print("    => Phi IS a Poisson endomorphism of (C[x,p], standard bracket).")

# (3) det JPhi = 1: block lower-triangular [[JK, 0], [*, M]]
detN = det3(N)
assert detN == sc(ONE, Fr(4)), "det adj^T != 4"
print("(3) det M = det(adj^T)/(-2)^3 = 4/(-8) = -1/2 ; det JPhi = det JK * det M = 1  OK")

# (4) non-injectivity: evaluate Phi at the three collision points with p = 0
def ev(p, pt):
    return sum(c * Fr(1) * pt[0] ** e[0] * pt[1] ** e[1] * pt[2] ** e[2]
               * pt[3] ** e[3] * pt[4] ** e[4] * pt[5] ** e[5] for e, c in p.items())
pts = [(Fr(0), Fr(0), Fr(-1, 4), Fr(0), Fr(0), Fr(0)),
       (Fr(1), Fr(-3, 2), Fr(13, 2), Fr(0), Fr(0), Fr(0)),
       (Fr(-1), Fr(3, 2), Fr(13, 2), Fr(0), Fr(0), Fr(0))]
imgs = [tuple(ev(g, pt) for g in Ximg + Pimg) for pt in pts]
print("(4) Phi at the three lifted collision points:")
for pt, im in zip(pts, imgs):
    print("    Phi(%s) = %s" % (tuple(str(c) for c in pt), tuple(str(c) for c in im)))
assert imgs[0] == imgs[1] == imgs[2]
print("    all three EQUAL: Phi is NOT injective => NOT an automorphism.")

print("\nCONCLUSION: Phi = cotangent lift of the Keller kernel is an explicit")
print("unimodular Poisson endomorphism of C[x1..x3,p1..p3] that is not an")
print("automorphism: the POISSON CONJECTURE on C^6 is FALSE, and by identity")
print("padding the STABLE POISSON CONJECTURE is false for all C^{2n}, n >= 3.")
print("Classical limit of the S141 Dixmier endomorphism (symbols of its 9")
print("Weyl commutator identities).  DONE.")
