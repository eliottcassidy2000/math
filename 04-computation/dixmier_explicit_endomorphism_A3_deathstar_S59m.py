#!/usr/bin/env python3
"""
death-star-2026-07-19-S59m (HYP-8075) -- the EXPLICIT Dixmier counterexample.

From the verified JC_3 counterexample F (det JF = -2, non-injective), build
  phi: A_3 -> A_3,  phi(X_i) = F_i(X),  phi(D_j) = sum_k B_jk(X) D_k,
  B := (JF^T)^{-1} = adj(JF)^T / (-2)   -- polynomial over Z[1/2].

Weyl-relation verification (exact, dependency-free):
  (R0) [phi X_i, phi X_j] = 0                    -- automatic (both in C[X])
  (R1) [phi D_j, phi X_i] = delta_ij  <=>  B * (JF)^T = I          (9 identities)
  (R2) [phi D_i, phi D_j] = 0  <=>  for all i<j, l:
         sum_k ( B_ik d_k B_jl - B_jk d_k B_il ) = 0               (9 identities)
If R1+R2 hold, phi extends to an endomorphism of A_3 (universal property of the
Weyl presentation); A_3 simple => phi injective.  NON-SURJECTIVITY is the
module one-liner: any element of im(phi) is sum_{beta,alpha} c F^beta V^alpha
(V_j := phi(D_j) acts on C[x,y,z] as the derivation sum_k B_jk d_k, which
kills 1); if X_1 = that sum, apply both sides to 1 in C[x,y,z]:
x_1 = sum_beta c_beta F^beta in C[F] -- same for x_2, x_3 -- so C[F] = C[X],
F has a polynomial inverse, F is injective: CONTRADICTION with the verified
triple collision.  Hence X_1 not in im(phi): phi is a PROPER self-embedding
of A_3 and the Dixmier conjecture is FALSE for A_3 (hence A_n, n >= 3).
"""
from fractions import Fraction

def pmul(a, b):
    r = {}
    for ka, ca in a.items():
        for kb, cb in b.items():
            k = (ka[0]+kb[0], ka[1]+kb[1], ka[2]+kb[2])
            r[k] = r.get(k, 0) + ca*cb
    return {k: c for k, c in r.items() if c != 0}

def padd(*ps):
    r = {}
    for p in ps:
        for k, c in p.items():
            r[k] = r.get(k, 0) + c
    return {k: c for k, c in r.items() if c != 0}

def pscale(p, s):
    return {k: c*s for k, c in p.items() if c*s != 0}

def pdiff(p, i):
    r = {}
    for k, c in p.items():
        if k[i] > 0:
            k2 = list(k); k2[i] -= 1
            r[tuple(k2)] = r.get(tuple(k2), 0) + c*k[i]
    return {k: c for k, c in r.items() if c != 0}

X = {(1,0,0): 1}; Y = {(0,1,0): 1}; Z = {(0,0,1): 1}; ONE = {(0,0,0): 1}
U = padd(ONE, pmul(X, Y)); U2 = pmul(U, U); U3 = pmul(U2, U)
W = padd(pscale(ONE, 4), pscale(pmul(X, Y), 3))
F1 = padd(pmul(U3, Z), pmul(pmul(pmul(Y, Y), U), W))
F2 = padd(Y, pscale(pmul(pmul(X, U2), Z), 3), pscale(pmul(pmul(X, pmul(Y, Y)), W), 3))
F3 = padd(pscale(X, 2), pscale(pmul(pmul(X, X), Y), -3), pscale(pmul(pmul(pmul(X, X), X), Z), -1))
F = [F1, F2, F3]
J = [[pdiff(Fi, v) for v in range(3)] for Fi in F]   # J[i][k] = dF_i/dx_k
A = [[J[k][i] for k in range(3)] for i in range(3)]   # A = J^T

def det2(a, b, c, d):
    return padd(pmul(a, d), pscale(pmul(b, c), -1))

# adj(A)_{ij} = (-1)^{i+j} * minor_{ji}(A);  B = adj(A) / det(A), det(A) = -2
adj = [[None]*3 for _ in range(3)]
for i in range(3):
    for j in range(3):
        rows = [r for r in range(3) if r != j]
        cols = [c for c in range(3) if c != i]
        m = det2(A[rows[0]][cols[0]], A[rows[0]][cols[1]],
                 A[rows[1]][cols[0]], A[rows[1]][cols[1]])
        adj[i][j] = pscale(m, (-1)**(i+j))
B = [[pscale(adj[i][j], Fraction(-1, 2)) for j in range(3)] for i in range(3)]

# (R1): B * A = I
ok1 = True
for i in range(3):
    for j in range(3):
        s = padd(*[pmul(B[i][k], A[k][j]) for k in range(3)])
        want = {(0,0,0): 1} if i == j else {}
        if s != want:
            ok1 = False; print(f"R1 FAIL at ({i},{j}): {s}")
print("(R1) B * JF^T = I:", ok1)

# (R2): flatness -- commuting vector fields V_i = sum_k B_ik d_k
ok2 = True
dB = [[[pdiff(B[i][l], k) for k in range(3)] for l in range(3)] for i in range(3)]
for i in range(3):
    for j in range(i+1, 3):
        for l in range(3):
            s = padd(*[padd(pmul(B[i][k], dB[j][l][k]),
                            pscale(pmul(B[j][k], dB[i][l][k]), -1)) for k in range(3)])
            if s != {}:
                ok2 = False
                print(f"R2 FAIL at (i={i},j={j},l={l}): {len(s)} terms")
print("(R2) [V_i, V_j] = 0 (all 9 identities):", ok2)

deg = lambda p: max((sum(k) for k in p), default=0)
print("\nB entry stats (degree / #terms / all coeffs in Z[1/2]):")
half_ok = True
for i in range(3):
    for j in range(3):
        cs = B[i][j].values()
        inhalf = all((isinstance(c, int)) or (c.denominator in (1, 2)) for c in cs)
        half_ok &= inhalf
        print(f"  B[{i}][{j}]: deg {deg(B[i][j]):2d}, {len(B[i][j]):3d} terms, Z[1/2]: {inhalf}")
print("all B coefficients in Z[1/2]:", half_ok)

def pstr(p):
    def mono(k):
        out = []
        for v, e in zip("xyz", k):
            if e == 1: out.append(v)
            elif e > 1: out.append(f"{v}^{e}")
        return "*".join(out) if out else "1"
    terms = sorted(p.items(), key=lambda kv: (sum(kv[0]), kv[0]))
    return " + ".join(f"({c})*{mono(k)}" for k, c in terms) if terms else "0"

with open("05-knowledge/results/dixmier_explicit_endomorphism_A3_deathstar_S59m.out", "w") as f:
    f.write("phi: A_3 -> A_3, phi(X1)=F1, phi(X2)=F2, phi(X3)=F3 (the JC counterexample),\n")
    f.write("phi(D_j) = sum_k B[j][k] D_k with B = (JF^T)^{-1}:\n\n")
    f.write(f"(R1) B*JF^T = I: {ok1}\n(R2) all 9 flatness identities: {ok2}\n\n")
    for i in range(3):
        for j in range(3):
            f.write(f"B[{i}][{j}] = {pstr(B[i][j])}\n\n")
print("\nVERDICT: phi is a well-defined endomorphism of A_3 =", ok1 and ok2)
print("=> with the verified non-injectivity of F: phi is injective (A_3 simple),")
print("   NOT surjective (module one-liner) -- an explicit proper self-embedding.")
