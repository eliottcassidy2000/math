"""A failed proof strategy for the A_24 deep-hole theta identity (session opus-2026-10-21-S52; HYP-9363).

Claim:  theta_(Lambda + c)(tau) - theta_(N + rho_N/25)(tau) = 25 Delta(tau),
where Lambda is the Leech lattice, c its A_24 deep hole, N = N(A_24) and rho_N its Weyl vector (h = 25).

Realisation in II_{25,1} (tournament coordinates of THM-4752):
- rho = (0,1,...,24 | 70) is the Leech cusp and nu = (1/2^25 | 5/2) the A_24 cusp, with rho . nu = -25.
- M = {rho, nu}^perp is positive definite of rank 24 and determinant 625; M^*/M = D is (Z/25)^2.
- Lambda = rho^perp/rho and N = nu^perp/nu are realised inside M (x) Q as the lifts orthogonal to the other cusp:
    Lambda_M = { x + (x . nu/25) rho : x in rho^perp },   N_M = { y + (y . rho/25) nu : y in nu^perp }.
  Their images in D are two isotropic subgroups H_L, H_N of order 25.
- The deep hole is c_M = nu/25 + rho - rho' (rho' = rho + e_0 - e_1).
- rho_N/25 is lifted as w_M = rho_N/25 + ((rho_N/25) . rho/25) nu, with rho_N = rho + (rho . nu') nu - 25 nu' and
  nu' = nu - e_0 - e_1.
- Write A = c_M + H_L and B = w_M + H_N, as sets of classes in D.

theta_A - theta_B = sum_gamma f_gamma theta_(M + gamma), with f = 1_A - 1_B.

Theta_M = sum_gamma theta_(M+gamma) e_gamma is a weight-12 modular form for the Weil representation rho_M (Borcherds
1998, Thm 4.1; sig M = 24 so the eighth-root phase is 1). Hence <Theta_M, f> is a modular form of weight 12 for
SL_2(Z) as soon as f is fixed by the transposed generators:
- T:  e(q(gamma)) f_gamma = f_gamma;
- S:  (1/25) sum_delta e(-(gamma, delta)) f_delta = f_gamma.
These are exact finite checks over 625 classes, done in Z[zeta_25]. If they held, the difference would be a level-1
cusp form, i.e. a multiple of Delta, and its q^1 coefficient (the 25 hole vertices) would make it 25 Delta.

RESULT: they do NOT hold. A and B are lines of D = (Z/25)^2 (u = -1 and v = -1 in coordinates (x.rho, x.nu) mod 25)
meeting in one class, so f is supported on non-isotropic classes. The identity verified numerically
(hyp9363_holes.py) therefore needs relations among different coset theta series:
    theta_(M + (-1, v)) = theta_(M + (v, -1)) for v != 0,
plus one integral relation at v = 0.
- For units v the first relation would follow if the torus (u, v) -> (a u, v/a) of O(D) were induced by isometries
  of M.
- For v in 5Z/25, and for the integral relation, O(D) gives nothing: its other elements contain the swap
  (u, v) -> (v, u), which would force theta_Lambda = theta_N.
This script certifies the structural facts and records that the naive sufficient condition fails.
Needs cypari.  Prints ALL CHECKS PASSED (the expected failures are asserted as such)."""
import itertools
from fractions import Fraction as Fr
from cypari import pari

pari.allocatemem(10 ** 9)
FAIL = []


def check(c, msg):
    print(('  ok   ' if c else '  FAIL ') + msg, flush=True)
    if not c:
        FAIL.append(msg)


n = 25


def lor(x, y):
    return sum(Fr(a) * Fr(b) for a, b in zip(x[:-1], y[:-1])) - Fr(x[-1]) * Fr(y[-1])


H = Fr(1, 2)
rho = [Fr(i) for i in range(n)] + [Fr(70)]
nu = [H] * n + [5 * H]
rhop = list(rho); rhop[0] += 1; rhop[1] -= 1
nup = list(nu); nup[0] -= 1; nup[1] -= 1
check(lor(rho, rho) == 0 and lor(nu, nu) == 0 and lor(rhop, rhop) == 0 and lor(nup, nup) == 0 and lor(rho, rhop) == -1
      and lor(nu, nup) == -1 and lor(rho, nu) == -25,
      "rho, rho' = rho + e_0 - e_1, nu, nu' = nu - e_0 - e_1 null; rho.rho' = nu.nu' = -1; rho.nu = -25")

# II_{25,1} basis (doubled coordinates, HNF)
gens = []
for i in range(n - 1):
    v = [0] * (n + 1); v[i] = 2; v[i + 1] = -2; gens.append(v)
v = [0] * (n + 1); v[n - 2] = 2; v[n - 1] = 2; gens.append(v)
gens.append([1] * (n + 1))
v = [0] * (n + 1); v[n - 1] = 2; v[n] = 2; gens.append(v)
v = [0] * (n + 1); v[0] = 4; gens.append(v)
Am = pari.matrix(n + 1, len(gens), [gens[j][i] for i in range(n + 1) for j in range(len(gens))])
Hb = pari.mathnf(Am)
basis = [[Fr(int(Hb[i, j]), 2) for i in range(n + 1)] for j in range(n + 1)]


def kernel_basis(forms):
    F = pari.matrix(len(forms), n + 1, [int(lor(b, f)) for f in forms for b in basis])
    K = pari.matkerint(F)
    d = int(pari.matsize(K)[1])
    return [[sum(int(K[i, j]) * basis[i][t] for i in range(n + 1)) for t in range(n + 1)] for j in range(d)]


Mb = kernel_basis([rho, nu])                       # basis of M
GM = [[lor(a, b) for b in Mb] for a in Mb]
GMp = pari.matrix(24, 24, [int(x) for row in GM for x in row])
check(len(Mb) == 24 and int(pari.matdet(GMp)) == 625 and all(GM[i][i] % 2 == 0 for i in range(24)),
      'M = {rho, nu}^perp: even, positive definite rank 24, det 625')
GMinv = pari.matsolve(GMp, pari.matid(24))


def coords(x):
    """Coordinates of x in M (x) Q (x orthogonal to rho and nu) w.r.t. the basis Mb."""
    g = pari.matrix(24, 1, [lor(b, x) for b in Mb])
    a = GMinv * g
    return [Fr(str(a[i, 0])) for i in range(24)]


def vec(a):
    return [sum(a[j] * Mb[j][t] for j in range(24)) for t in range(n + 1)]


# D = M^*/M: elements are coordinate vectors a with G a integral, modulo Z^24.  Generate via the dual basis.
dual = [[Fr(str(GMinv[i, j])) for i in range(24)] for j in range(24)]   # columns of G^-1 = dual basis coords


def red(a):
    return tuple(t - (t.numerator // t.denominator) for t in a)


classes = {red([Fr(0)] * 24)}
frontier = list(classes)
while frontier:
    nxt = []
    for g in frontier:
        for d in dual:
            h = red([x + y for x, y in zip(g, d)])
            if h not in classes:
                classes.add(h); nxt.append(h)
    frontier = nxt
D = sorted(classes)
check(len(D) == 625, f'|D| = |M^*/M| = {len(D)}')


def qf(a):
    v = sum(a[i] * sum(GM[i][j] * a[j] for j in range(24)) for i in range(24)) / 2
    return v - (v.numerator // v.denominator)


def bf(a, b):
    v = sum(a[i] * sum(GM[i][j] * b[j] for j in range(24)) for i in range(24))
    return v - (v.numerator // v.denominator)


# the two isotropic subgroups
def lift_L(x):            # x in rho^perp: x + (x.nu/25) rho  (orthogonal to nu)
    return [t + lor(x, nu) / 25 * r for t, r in zip(x, rho)]


def lift_N(y):            # y in nu^perp: y + (y.rho/25) nu  (orthogonal to rho)
    return [t + lor(y, rho) / 25 * s for t, s in zip(y, nu)]


rho_perp = kernel_basis([rho])
nu_perp = kernel_basis([nu])
HL = {red([Fr(0)] * 24)}
for x in rho_perp:
    a = red(coords(lift_L(x)))
    new = set()
    for g in HL:
        h = g
        for _ in range(25):
            h = red([p + q for p, q in zip(h, a)])
            new.add(h)
    HL |= new
HN = {red([Fr(0)] * 24)}
for y in nu_perp:
    a = red(coords(lift_N(y)))
    new = set()
    for g in HN:
        h = g
        for _ in range(25):
            h = red([p + q for p, q in zip(h, a)])
            new.add(h)
    HN |= new
check(len(HL) == 25 and len(HN) == 25 and all(qf(g) == 0 for g in HL | HN) and len(HL & HN) == 1,
      'Lambda_M/M and N_M/M are isotropic subgroups H_L, H_N of order 25 with H_L n H_N = 0 (D = H_L + H_N)')

cM = [a / 25 + b - c for a, b, c in zip(nu, rho, rhop)]
rhoN = [r + lor(rho, nup) * s - 25 * t for r, s, t in zip(rho, nu, nup)]
y0 = [t / 25 for t in rhoN]
wM = [t + lor(y0, rho) / 25 * s for t, s in zip(y0, nu)]
check(lor(cM, rho) == 0 and lor(cM, nu) == 0 and lor(cM, cM) == 2 and lor(wM, rho) == 0 and lor(wM, nu) == 0
      and lor(wM, wM) == Fr(52, 25) and lor(rhoN, rhoN) == 1300,
      'c_M (|c|^2 = 2) and w_M = lift of rho_N/25 (|rho_N|^2 = 1300 = 2h(h+1), |w_M|^2 = 52/25) lie in M (x) Q')
ac, aw = red(coords(cM)), red(coords(wM))
check(all((sum(GM[i][j] * ac[j] for j in range(24))).denominator == 1 for i in range(24))
      and all((sum(GM[i][j] * aw[j] for j in range(24))).denominator == 1 for i in range(24)),
      'c_M and w_M lie in M^*')
A = {red([p + q for p, q in zip(ac, g)]) for g in HL}
B = {red([p + q for p, q in zip(aw, g)]) for g in HN}
negA = {red([-t for t in g]) for g in A}
negB = {red([-t for t in g]) for g in B}
print(f'     |A| = {len(A)}, |B| = {len(B)}, |A n B| = {len(A & B)}, A = -A: {A == negA}, B = -B: {B == negB}')
check(len(A) == 25 and len(B) == 25 and len(A & B) == 1, 'A = c_M + H_L and B = w_M + H_N are lines of D meeting in exactly one class')
# f symmetrised: theta_(-gamma) = theta_gamma
f = {g: Fr(0) for g in D}
for g in A:
    f[g] += H; f[red([-t for t in g])] += H
for g in B:
    f[g] -= H; f[red([-t for t in g])] -= H
supp = [g for g in D if f[g] != 0]
check(not all(qf(g) == 0 for g in supp), f'naive T-condition FAILS as expected: f = 1_A - 1_B (symmetrised) is supported on {len(supp)} classes, not all isotropic')
# S: (1/25) sum_delta e(-(g, delta)) f_delta = f_g, exactly in Z[zeta_25]: reduce exponent vectors mod Phi_25
Phi25 = [1 if i % 5 == 0 else 0 for i in range(21)]   # 1 + x^5 + ... + x^20


def reduce_cyc(vecz):
    v = list(vecz)
    for d in range(len(v) - 1, 19, -1):
        c_ = v[d]
        if c_:
            for i in range(21):
                v[d - 20 + i] -= c_ * Phi25[i]
    return v[:20]


okS = True
for g in D:
    acc = [Fr(0)] * 25
    for dlt in supp:
        e = (-bf(g, dlt) * 25) % 25
        acc[int(e)] += f[dlt]
    lhs = reduce_cyc([x / 25 for x in acc])
    rhs = [f[g]] + [Fr(0)] * 19
    okS &= lhs == rhs
check(not okS, 'naive S-condition FAILS as expected (exact in Z[zeta_25]): the identity is not a consequence of Weil-invariance of f alone')
print()
print('ALL CHECKS PASSED' if not FAIL else f'FAILED: {FAIL}')
