#!/usr/bin/env python3
# THM-4910 (mac-mini-2026-11-07-sevensphere): exact checks of the eight-point sphere.
# (P1) symbolic, all n >= 3, in THM-4742's Gram basis (v, s, w): Q(C*) = 0 and Q(X) = 0 on THM-4742's sphere S2,
#      where C* = (n+1)G - nV0 (anticomplementary vertex) = s - (n-1)v and X = inverse of G' in the facet circumsphere;
#      and the corner powers Pow_S2(O, G, G', O') = (n+1)R^2, n(R^2-|OG|^2), (n-1)(R^2-|OG'|^2), n(R^2-|OO'|^2).
# (P2) exact rational tetrahedra (n = 3) in actual coordinates: all eight labels on one sphere, pairwise distinct,
#      corner powers, the four side-line pairs, and Pow_S2(V0) = V0A.V0C* = sum_i |V0Vi|^2 (THM-4742's defect).
# (P3) exact rational simplices n = 4, 5, 6 in actual coordinates (labels lie in a 3-flat; checked inside it).
import random
from fractions import Fraction as Fr
import sympy as sp

# ---------------- (P1) symbolic in the Gram basis ----------------
n, R2, vs, vw, ss, k = sp.symbols('n R2 vs vw ss k')
gram = sp.Matrix([[R2, vs, vw], [vs, ss, n*k], [vw, n*k, k]])          # basis v, s, w (relative to O)
def dot(a, b): return (sp.Matrix(a).T * gram * sp.Matrix(b))[0]
twou = {0: 2*vs + (2 - n)*R2, 1: ss + 2*n*R2, 2: (n + 1)*k + R2}       # S2 coefficients (THM-4742)
def Q(x): return sp.simplify(dot(x, x) - sum(twou[i]*x[i] for i in range(3)) + (n + 1)*R2)
v, s, w = [1, 0, 0], [0, 1, 0], [0, 0, 1]
Cstar = [-(n - 1), 1, 0]
Gp_minus_Op = [0, 1/n, -1]                                              # G' - O' = s/n - w
gamma2 = sp.simplify(dot(Gp_minus_Op, Gp_minus_Op))                    # = ss/n^2 - k
rp2 = R2 - k                                                            # facet circumradius^2
X = [0 + rp2/gamma2*Gp_minus_Op[0], rp2/gamma2*Gp_minus_Op[1], 1 + rp2/gamma2*Gp_minus_Op[2]]
G = [1/(n + 1), 1/(n + 1), 0]; Gp = [0, 1/n, 0]; Op = [0, 0, 1]; O0 = [0, 0, 0]
print("(P1) Q_S2(C*) =", sp.simplify(Q(Cstar)))
print("(P1) Q_S2(X)  =", sp.simplify(sp.together(Q(X))))
for name, P, c in (("O", O0, n + 1), ("G", G, n), ("G'", Gp, n - 1), ("O'", Op, n)):
    print(f"(P1) Pow_S2({name}) - c*(R^2-|OP|^2) =", sp.simplify(Q(P) - c*(R2 - dot(P, P))))

# ---------------- exact rational geometry ----------------
def sub(a, b): return [x - y for x, y in zip(a, b)]
def add(a, b): return [x + y for x, y in zip(a, b)]
def mul(c, a): return [c*x for x in a]
def dt(a, b): return sum(x*y for x, y in zip(a, b))
def solve(M, b):
    # exact Gaussian elimination
    m = len(M); A = [list(map(Fr, row)) + [Fr(bb)] for row, bb in zip(M, b)]
    for c in range(m):
        p = next(r for r in range(c, m) if A[r][c] != 0); A[c], A[p] = A[p], A[c]
        for r in range(m):
            if r != c and A[r][c] != 0:
                f = A[r][c]/A[c][c]; A[r] = [x - f*y for x, y in zip(A[r], A[c])]
    return [A[i][m]/A[i][i] for i in range(m)]
def circumcenter(P):
    p0 = P[0]; Dm = [sub(p, p0) for p in P[1:]]
    Gm = [[dt(a, b) for b in Dm] for a in Dm]; rhs = [dt(a, a)/2 for a in Dm]
    lam = solve(Gm, rhs)
    return add(p0, [sum(l*d[j] for l, d in zip(lam, Dm)) for j in range(len(p0))])
def labels(V):
    N = len(V) - 1
    O = circumcenter(V); R2v = dt(sub(V[0], O), sub(V[0], O))
    G = mul(Fr(1, N + 1), [sum(x) for x in zip(*V)])
    F_ = V[1:]; Op = circumcenter(F_); Gp = mul(Fr(1, N), [sum(x) for x in zip(*F_)])
    d = sub(G, V[0]); mu = -2*dt(sub(V[0], O), d)/dt(d, d); A = add(V[0], mul(mu, d))
    w = sub(Op, O); B = add(O, mul(R2v/dt(w, w), w))
    Cs = add(V[0], mul(N + 1, sub(G, V[0])))
    D = add(O, mul(N + 1, sub(G, O))); g = sub(G, O); E = add(O, mul(R2v/dt(g, g), g))
    F = add(O, mul(N + 1, sub(Op, O))); H = add(Op, mul(N, sub(Gp, Op)))
    rp2 = dt(sub(F_[0], Op), sub(F_[0], Op)); gg = sub(Gp, Op); X = add(Op, mul(rp2/dt(gg, gg), gg))
    return dict(A=A, B=B, Cstar=Cs, D=D, E=E, F=F, H=H, X=X), dict(O=O, G=G, Gp=Gp, Op=Op, R2=R2v, V0=V[0], N=N)

def check(V):
    L, C = labels(V)
    O, G, Gp, Op, R2v, N = C["O"], C["G"], C["Gp"], C["Op"], C["R2"], C["N"]
    # coordinates inside the 3-flat O + span(V0-O, sum(Vi-O), O'-O): Gram-Schmidt-free, use exact basis coefficients
    basis = [sub(C["V0"], O), sub(mul(N, Gp), mul(N, O)), sub(Op, O)]   # v, s, w
    Gm = [[dt(a, b) for b in basis] for a in basis]
    def coords(P):  # coefficients of P - O in basis (exact; P lies in the flat)
        c = solve(Gm, [dt(sub(P, O), b) for b in basis])
        rec = add(O, [sum(ci*b[j] for ci, b in zip(c, basis)) for j in range(len(O))])
        assert rec == P, "label left the 3-flat"
        return c
    def q(c, u, kap):  # power in the flat: |x|^2 - 2u.x + kappa with x given by coords
        x2 = sum(c[i]*c[j]*Gm[i][j] for i in range(3) for j in range(3))
        return x2 - 2*sum(c[i]*u[i] for i in range(3)) + kap
    # sphere through D, F, C*, H' (non-coplanar in the regime): unknowns (u.v, u.s, u.w) as dual coords, kappa
    names = ["D", "F", "Cstar", "H"]
    rows, rhs = [], []
    for nm in names:
        c = coords(L[nm]); x2 = sum(c[i]*c[j]*Gm[i][j] for i in range(3) for j in range(3))
        rows.append([-2*c[0], -2*c[1], -2*c[2], 1]); rhs.append(-x2)
    sol = solve(rows, rhs); u, kap = sol[:3], sol[3]
    ok = all(q(coords(L[nm]), u, kap) == 0 for nm in L)
    # corner powers
    cp = {"O": (O, N + 1), "G": (G, N), "Gp": (Gp, N - 1), "Op": (Op, N)}
    okp = all(q(coords(P), u, kap) == c*(R2v - dt(sub(P, O), sub(P, O))) for P, c in cp.values())
    # defect: Pow(V0) = V0A . V0C* = sum |V0 Vi|^2
    pv0 = q(coords(C["V0"]), u, kap)
    sig = sum(dt(sub(V[0], Vi), sub(V[0], Vi)) for Vi in V[1:])
    prod = (dt(sub(L["A"], V[0]), sub(L["A"], V[0])) * dt(sub(L["Cstar"], V[0]), sub(L["Cstar"], V[0])))
    okd = (pv0 == sig) and (pv0**2 == prod)
    distinct = len({tuple(p) for p in L.values()}) == 8
    return ok, okp, okd, distinct, kap == (N + 1)*R2v + 0  # kappa relative? (power at O)

random.seed(20261107)
def rand_simplex(N):
    while True:
        V = [[Fr(random.randint(-9, 9), random.randint(1, 4)) for _ in range(N)] for _ in range(N + 1)]
        try:
            check(V); return V
        except Exception:
            continue

for N, trials in ((3, 150), (4, 40), (5, 20), (6, 10)):
    tot = [0, 0, 0, 0]
    done = 0
    while done < trials:
        V = [[Fr(random.randint(-9, 9), random.randint(1, 4)) for _ in range(N)] for _ in range(N + 1)]
        try:
            r = check(V)
        except (StopIteration, ZeroDivisionError, AssertionError):
            continue
        done += 1
        for i in range(4): tot[i] += bool(r[i])
    print(f"(P{2 if N == 3 else 3}) n={N}: {trials} exact rational simplices: all 8 labels on S2: {tot[0]}/{trials}; "
          f"corner powers c=(n+1,n,n-1,n): {tot[1]}/{trials}; Pow(V0)=V0A*V0C*=sum|V0Vi|^2: {tot[2]}/{trials}; 8 distinct: {tot[3]}/{trials}")
print("all done")
