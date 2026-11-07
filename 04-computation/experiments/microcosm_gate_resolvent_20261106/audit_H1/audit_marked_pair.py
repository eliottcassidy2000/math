# audit_H1: independent exact checks of THM-4900 (1)-(5): gate/padding transfer on (N,P) (N via det(zI-M+J)-det(zI-M)), eigen-coordinates, nested-gate closed form, phi forms/bounds/doubling invariance, D/E formulas, series law, golden map, Fibonacci orbit, D-Gamma interaction.
import random
import sys
from fractions import Fraction
import sympy as sp
from audit_tournament_lib import (rand_tour, classes, transitive, dominance, scores, gate, pad, double, join,
                                  phi, phi_gram, F_at)

z, s, w = sp.symbols('z s w')
rng = random.Random(20261106)
fails = []


def check(cond, msg):
    if not cond:
        fails.append(msg)
        print("FAIL:", msg)


def marked_pair(A):
    """(N, P) with P = det(zI - M) and N = u^T adj(zI - M) u = det(zI - M + J) - det(zI - M) (matrix determinant lemma)."""
    n = len(A)
    M = sp.Matrix(dominance(A))
    P = sp.expand(M.charpoly(z).as_expr())
    J = sp.ones(n, n)
    PJ = sp.expand((M - J).charpoly(z).as_expr())   # det(zI - (M - J)) = det(zI - M + J)
    return sp.expand(PJ - P), P


def all_labelled(n):
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    for bits in range(1 << len(pairs)):
        A = [[0] * n for _ in range(n)]
        for k, (i, j) in enumerate(pairs):
            if (bits >> k) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
        yield A


# ---------------- (1) gate and padding transfers ----------------
cores = []
for n in range(1, 5):
    cores += list(all_labelled(n))          # 1 + 2 + 8 + 64 = 75 labelled cores
for n in (5, 6):
    cores += classes(n)                     # 12 + 56 classes
cnt = 0
for A in cores:
    N0, P0 = marked_pair(A)
    G = gate(A, rng)
    E = pad(A, rng)
    NG, PG = marked_pair(G)
    NE, PE = marked_pair(E)
    check(sp.expand(NG - ((z**2 + 9) * N0 + 2 * z * P0)) == 0, f"N_Gamma n={len(A)}")
    check(sp.expand(PG - (2 * z * N0 + (z**2 + 1) * P0)) == 0, f"P_Gamma n={len(A)}")
    check(sp.expand(NE - ((z**2 + 1) * N0 + 2 * z * P0)) == 0, f"N_E n={len(A)}")
    check(sp.expand(PE - PG) == 0, f"P_E = P_Gamma n={len(A)}")
    cnt += 1
print(f"(1) gate transfer [[z^2+9,2z],[2z,z^2+1]] and padding [[z^2+1,2z],[2z,z^2+1]] on (N,P), P_Gamma = P_E: "
      f"{cnt} cores (all labelled n<=4, all classes n=5,6), randomly relabelled shells")

# eigen-structure, symbolic with s^2 = z^2 + 4
Gm = sp.Matrix([[z**2 + 9, 2 * z], [2 * z, z**2 + 1]])
def red(e):
    e = sp.expand(e)
    return sp.expand(sp.rem(sp.Poly(e, s), sp.Poly(s**2 - z**2 - 4, s)).as_expr())
for sign in (1, -1):
    v = sp.Matrix([z, -2 + sign * s])
    lam = (s + sign) ** 2
    diff = Gm * v - lam * v
    check(all(red(d) == 0 for d in diff), f"eigenvector sign {sign}")
check(red(Gm.det() - ((s + 1) ** 2 * (s - 1) ** 2)) == 0, "det G = (s^2-1)^2")
print("(1) eigenvectors (z, -2 +- s) with eigenvalues (s +- 1)^2 (s^2 = z^2+4): checked symbolically")

# ---------------- (2) nested gates closed form ----------------
def closed_form_P(N0, P0, k):
    Psip = z * N0 + (-2 + s) * P0
    Psim = z * N0 + (-2 - s) * P0
    num = (s + 1) ** (2 * k) * Psip - (s - 1) ** (2 * k) * Psim
    # divide by 2s: num is odd in s after reduction
    r = red(num)
    poly = sp.Poly(r, s)
    check(poly.degree() <= 1 and poly.coeff_monomial(1) == 0, "closed form numerator not s*f(z)")
    return sp.expand(poly.coeff_monomial(s) / 2)

ccount = 0
for n in range(1, 6):
    for trial in range(3):
        A = rand_tour(n, rng)
        N0, P0 = marked_pair(A)
        X = A
        for k in range(1, 5):
            X = gate(X, rng)
            if len(X) > 13:
                break
            NX, PX = marked_pair(X)
            check(sp.expand(PX - closed_form_P(N0, P0, k)) == 0, f"closed form n={n} k={k}")
            ccount += 1
print(f"(2) P_(Gamma^k T) = [(s+1)^2k Psi+ - (s-1)^2k Psi-]/(2s): {ccount} (core, k) pairs, exact")

# one-vertex formula as polynomial identity, k = 1..8
for k in range(1, 9):
    X = [[0]]
    for _ in range(k):
        X = gate(X)
    NX, PX = marked_pair(X)
    br = red((s + 1) ** (2 * k - 1) + (s - 1) ** (2 * k - 1))
    pol = sp.Poly(br, s)
    form = sp.expand(z * (z**2 + 3) * pol.coeff_monomial(s) / 2)
    check(pol.coeff_monomial(1) == 0 and sp.expand(PX - form) == 0, f"one-vertex formula k={k}")
print("(2) P_(Gamma^k K1) = z(z^2+3)[(s+1)^(2k-1) + (s-1)^(2k-1)]/(2s): k = 1..8 exact")

# ---------------- (3) phi ----------------
nreg = 0
tot = 0
for n in range(1, 9):
    cl = classes(n)
    for A in cl:
        f = phi(A)
        g = phi_gram(A)
        check(f == g, f"phi forms differ n={n}")
        check(0 < f <= n, f"phi bound n={n}")
        regular = all(2 * x == n - 1 for x in scores(A))
        check((f == n) == regular, f"equality iff regular n={n}")
        nreg += regular
        tot += 1
print(f"(3) phi = F_T(1) = u^T(I+M^TM)^(-1)u, 0 < phi <= n, phi = n iff regular: all {tot} classes of orders 1..8 "
      f"({nreg} regular)")
# phi equals THM-1950's total inverse-response s(T) = 1^T (I+K)^(-1) 1 with K = A - A^T = M
from audit_tournament_lib import solve_frac
same = 0
for n in range(1, 8):
    for A in classes(n):
        M = dominance(A)
        IK = [[(1 if i == j else 0) + M[i][j] for j in range(n)] for i in range(n)]
        sval = sum(solve_frac(IK, [1] * n))
        check(sval == phi(A), "phi != THM-1950 s")
        same += 1
print(f"(3) phi(T) = THM-1950's s(T) = 1^T(I+K)^(-1)1: {same} classes of orders 1..7")

dcount = 0
for n in range(1, 9):
    for trial in range(4):
        A = rand_tour(n, rng)
        D = double(A)
        check(phi(D) == phi(A), f"phi(DT) n={n}")
        # P_D(1) = 2^n P(1)^2 exactly
        PD1 = sp.Matrix(dominance(D)).charpoly(z).as_expr().subs(z, 1)
        P1 = sp.Matrix(dominance(A)).charpoly(z).as_expr().subs(z, 1)
        check(PD1 == 2**n * P1**2, f"P_D(1) n={n}")
        dcount += 1
print(f"(3) phi(DT) = phi(T) and P_D(1) = 2^n P(1)^2: {dcount} random cores of orders 1..8, exact")

# Session35 D formulas as identities (independent derivation check), small n
for n in range(1, 6):
    A = rand_tour(n, rng)
    N0, P0 = marked_pair(A)
    D = double(A)
    ND, PD = marked_pair(D)
    # P_D(z) = 2^n P(w)^2 with w^2 = (z^2+1)/2 ; N_D = 2^n (z/w) N(w) P(w)
    Pw2 = sp.expand(P0.subs(z, w) ** 2)
    Pw2 = sp.Poly(Pw2, w)
    # even polynomial in w: substitute w^2
    sub = sum(c * ((z**2 + 1) / 2) ** (m[0] // 2) for m, c in zip(Pw2.monoms(), Pw2.coeffs()) if m[0] % 2 == 0)
    check(all(m[0] % 2 == 0 for m in Pw2.monoms()), "P(w)^2 not even")
    check(sp.expand(PD - 2**n * sub) == 0, f"P_D formula n={n}")
    NP = sp.Poly(sp.expand(N0.subs(z, w) * P0.subs(z, w)), w)   # odd in w
    check(all(m[0] % 2 == 1 for m in NP.monoms()), "N(w)P(w) not odd")
    sub2 = sum(c * ((z**2 + 1) / 2) ** ((m[0] - 1) // 2) for m, c in zip(NP.monoms(), NP.coeffs()))
    check(sp.expand(ND - 2**n * z * sub2) == 0, f"N_D formula n={n}")
print("(3)/(5) Session35 D formulas P_D = 2^n P(w)^2, N_D = 2^n (z/w) N(w) P(w): exact, n = 1..5")

# ---------------- (4) series law, golden map, Fibonacci ----------------
for trial in range(25):
    A = rand_tour(rng.randint(1, 5), rng)
    B = rand_tour(rng.randint(1, 5), rng)
    NA, PA = marked_pair(A)
    NB, PB = marked_pair(B)
    NS, PS = marked_pair(join(A, B))
    check(sp.expand(PS - (PA * PB + NA * NB)) == 0, "P join")
    check(sp.expand(NS - (NA * PB + PA * NB)) == 0, "N join")
print("(4) series law: P_(A=>B) = P_A P_B + N_A N_B, N_(A=>B) = N_A P_B + P_A N_B (=> relativistic addition at every z): 25 random pairs")
C3 = gate([[0]])
check(phi(join(C3, C3)) == Fraction(3, 5), "phi(C3=>C3)")
for n in range(1, 8):
    for trial in range(4):
        A = rand_tour(n, rng)
        f = phi(A)
        check(phi(gate(A, rng)) == (5 * f + 1) / (f + 1), "golden map")
        check(phi(pad(A, rng)) == 1, "phi(E T) = 1")
        check(phi(join([[0]], A)) == 1 and phi(join(A, [[0]])) == 1, "K1 absorbing")
print("(4) phi(Gamma T) = (5 phi+1)/(phi+1), phi(E T) = 1, K1 absorbing, phi(C3=>C3) = 3/5: checked")
fib = [0, 1]
while len(fib) < 40:
    fib.append(fib[-1] + fib[-2])
def F(i):
    return 1 if i == -1 else fib[i]
X = [[0]]
orbit = []
for k in range(0, 11):
    if k > 0:
        X = gate(X)
    f = phi(X)
    orbit.append(f)
    check(f == Fraction(F(2 * k + 2), F(2 * k - 1)), f"Fibonacci k={k}")
print("(4) phi(Gamma^k K1), k=0..10:", [str(x) for x in orbit])

# ---------------- (5) D and Gamma ----------------
bad_gate = 0
g = lambda x: (5 * x + 1) / (x + 1)
for n in range(1, 6):
    for A in classes(n):
        DG = double(gate(A))
        GD = gate(double(A))
        GGD = gate(gate(double(A)))
        f = phi(A)
        check(phi(DG) == phi(GD) == g(f), "commutation at z=1")
        check(phi(DG) != phi(GGD), "BS relation unexpectedly holds")
        if any(x == 1 for x in scores(DG)):
            bad_gate += 1
        if n <= 4:
            PDG = sp.Matrix(dominance(DG)).charpoly(z).as_expr()
            PDE = sp.Matrix(dominance(double(pad(A)))).charpoly(z).as_expr()
            check(sp.expand(PDG - PDE) == 0, "P_(D Gamma T) = P_(D E T)")
print(f"(5) phi(D Gamma T) = phi(Gamma D T) = g(phi T) != phi(Gamma^2 D T), all classes n<=5; "
      f"D(Gamma T) with an out-degree-1 vertex: {bad_gate}; P_(D Gamma T) = P_(D E T) for n <= 4")

print("FAILURES:", len(fails))
sys.exit(1 if fails else 0)
