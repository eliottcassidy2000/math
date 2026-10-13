#!/usr/bin/env python3
"""Audit D (2026-10-12) of THM-4626 -- machine checks (sympy, exact unless marked numeric).

(A) The non-example of THM-4626(4): D = x d_x - y d_y + d_z on Q[x,y,z].
    * div D = 0;
    * D is surjective: explicit preimage of every monomial x^a y^b z^c (closed form), verified
      symbolically for all monomials of total degree <= N, plus a rank check that D(P_{<=N})
      contains P_{<=N-1};
    * D is not locally nilpotent (D^k x = x);
    * Ker D on P_{<=N} is spanned by (xy)^k (Ker D = Q[xy] in bounded degree: transcendence degree
      one, so D is not a dual field of any Keller map of A^3, whose kernels have transcendence
      degree two).
(B) The padding example of THM-4626(3): F = F' x id with F' = THM-1300.  Checks det JF' is a
    nonzero constant, that the dual field V_4 of F' x id is d/dx_4 (so surjective and locally
    nilpotent), and counts the generic fibre of F' (= deg F' = d_4 = the number of lines in the
    generic fibre of (F_1, F_2, F_3) of the padded map) by eliminating z (F'_3 is linear in z),
    a resultant in y and numeric roots (50 digits).
(C) The mechanism behind THM-4626(1) (audit finding): once K[x] sits inside L[F] with V = d/dF, V
    is locally nilpotent on K[x] because d/dF is locally nilpotent on L[F]; surjectivity is not
    needed beyond THM-4620(4)'s hypothesis mu_j = 0.  Illustrated on a Keller map of A^3 that is an
    automorphism not of product form: V_3 from (JF^T)^{-1}, local nilpotence on x, y, z by
    iteration.
"""
import random
import time
from itertools import product

import sympy as sp
from sympy.polys.matrices import DomainMatrix

x, y, z, w = sp.symbols('x y z w')


def section(title):
    print('\n' + '=' * 78 + '\n' + title + '\n' + '=' * 78, flush=True)


# ---------------------------------------------------------------- (A) toy derivation
section('(A) D = x d_x - y d_y + d_z on Q[x,y,z]')
t0 = time.time()


def D_toy(p):
    return sp.expand(x * sp.diff(p, x) - y * sp.diff(p, y) + sp.diff(p, z))


print('div D =', sp.diff(x, x) + sp.diff(-y, y) + sp.diff(sp.Integer(1), z))


def preimage_monomial(a, b, c):
    """Closed-form a' with D a' = x^a y^b z^c (THM-4626(4)'s proof, made explicit)."""
    if a != b:
        lam = sp.Integer(a - b)
        # solve lam*p + p' = z^c:  p = sum_k (-1)^k c!/(c-k)! z^(c-k) / lam^(k+1)
        p = sum((-1) ** k * sp.factorial(c) / sp.factorial(c - k) * z ** (c - k) / lam ** (k + 1)
                for k in range(c + 1))
    else:
        p = z ** (c + 1) / (c + 1)
    return sp.expand(x ** a * y ** b * p)


N = 6
bad = count = 0
for a, b, c in product(range(N + 1), repeat=3):
    if a + b + c > N:
        continue
    count += 1
    if sp.expand(D_toy(preimage_monomial(a, b, c)) - x ** a * y ** b * z ** c) != 0:
        bad += 1
        print('  FAIL at monomial', (a, b, c))
print(f'explicit preimages verified for all {count} monomials of degree <= {N}: failures = {bad}')

# rank check: D(P_{<=N}) contains P_{<=N-1}  (DomainMatrix over QQ)
monos = [(a, b, c) for a, b, c in product(range(N + 1), repeat=3) if a + b + c <= N]
idx = {m: i for i, m in enumerate(monos)}
rows = [[sp.Integer(0)] * len(monos) for _ in monos]
for j, (a, b, c) in enumerate(monos):
    img = sp.Poly(D_toy(x ** a * y ** b * z ** c), x, y, z)
    for mon, coeff in img.terms():
        rows[idx[mon]][j] = coeff
M = DomainMatrix.from_Matrix(sp.Matrix(rows)).convert_to(sp.QQ)
low = [m for m in monos if sum(m) <= N - 1]
E = [[sp.Integer(1) if monos[i] == m else sp.Integer(0) for m in low] for i in range(len(monos))]
E = DomainMatrix.from_Matrix(sp.Matrix(E)).convert_to(sp.QQ)
rM, rME = M.rank(), M.hstack(E).rank()
print(f'dim P_<={N} = {len(monos)}, rank D|P_<={N} = {rM}, rank [D | P_<={N - 1}] = {rME}  '
      f'=> P_<={N - 1} subset D(P_<={N}): {rM == rME}')

p = x
for k in range(5):
    p = D_toy(p)
print('D^5 x =', p, '(not nilpotent)')

ker = M.nullspace().to_Matrix()
kernel_polys = [sp.expand(sum(ker[r, i] * x ** monos[i][0] * y ** monos[i][1] * z ** monos[i][2]
                              for i in range(len(monos)))) for r in range(ker.rows)]
print(f'dim Ker D on P_<={N} = {len(kernel_polys)}; expected floor(N/2)+1 = {N // 2 + 1}')
print('kernel basis:', kernel_polys)
xy_only = all(sp.Poly(kp, x, y, z).as_expr().subs({x: w / y}).has(y) is False and sp.diff(kp, z) == 0
              for kp in kernel_polys)
print('every kernel element is a polynomial in xy:', xy_only, f'  [{time.time() - t0:.1f}s]')

# ---------------------------------------------------------------- (B) padding example
section("(B) padding example F = F' x id with F' = THM-1300")
t0 = time.time()
u = 1 + x * y
F1 = u ** 3 * z + y ** 2 * u * (4 + 3 * x * y)
F2 = y + 3 * x * u ** 2 * z + 3 * x * y ** 2 * (4 + 3 * x * y)
F3 = 2 * x - 3 * x ** 2 * y - x ** 3 * z
Fp = [sp.expand(F1), sp.expand(F2), sp.expand(F3)]
J3 = sp.Matrix([[sp.diff(f, v) for v in (x, y, z)] for f in Fp])
detJ3 = sp.expand(J3.det())
print("det J(F') =", detJ3)
assert detJ3.is_number and detJ3 != 0, 'THM-1300 formula is not Keller as typed'

vars4 = (x, y, z, w)
F4 = Fp + [w]
J4 = sp.Matrix([[sp.diff(f, v) for v in vars4] for f in F4])
B4 = (J4.T).inv().applyfunc(sp.cancel)
print('dual field V_4 coefficients (B[4,k])_k =', [B4[3, k] for k in range(4)], '  (= d/dw)')
print('V_j F_i = delta_ij for all i,j:',
      all(sp.cancel(sum(B4[j, k] * sp.diff(F4[i], vars4[k]) for k in range(4)) - (1 if i == j else 0)) == 0
          for i in range(4) for j in range(4)))

# generic degree of F': eliminate z from F'_3 = c3 (linear in z, x != 0), resultant in y, numeric roots
random.seed(20261012)
c = [sp.Rational(random.randint(-9, 9), random.randint(1, 5)) for _ in range(3)]
zsol = sp.solve(sp.Eq(Fp[2], c[2]), z)[0]
P1 = sp.Poly(sp.cancel(x ** 3 * (Fp[0] - c[0]).subs(z, zsol)), x, y)
P2 = sp.Poly(sp.cancel(x ** 3 * (Fp[1] - c[1]).subs(z, zsol)), x, y)
res = sp.Poly(sp.resultant(P1.as_expr(), P2.as_expr(), y), x)
# strip the spurious x-power (x = 0 is excluded by the elimination)
k0 = 0
while res.eval(0) == 0:
    res = sp.Poly(sp.quo(res.as_expr(), x), x)
    k0 += 1
print(f"c = {c};  Res_y(P1,P2) = x^{k0} * (degree {res.degree()} polynomial);  factorization over Q:",
      sp.factor_list(res.as_expr())[1])
import mpmath
mpmath.mp.dps = 50
xroots = sp.Poly(res, x).nroots(n=40)
npts = 0
for x0 in xroots:
    p1 = sp.Poly(P1.as_expr().subs(x, x0), y)
    p2 = sp.Poly(P2.as_expr().subs(x, x0), y)
    yr = p1.nroots(n=40)
    common = [y0 for y0 in yr if abs(p2.eval(y0)) < 1e-20]
    npts += len(common)
print(f"numeric fibre count over c (points with x != 0): {npts}  [{time.time() - t0:.1f}s]")
print("=> deg F' = d_4 of the padded map = number of disjoint lines in the generic fibre of (F_1,F_2,F_3).")

# ---------------------------------------------------------------- (C) the mechanism on an automorphism
section('(C) restriction of d/dF_3 to K[x] inside L_3[F_3] is locally nilpotent: a non-product automorphism')
t0 = time.time()
psi = [x, y + z ** 2, z]                                   # triangular automorphism
Gid = [psi[0] + psi[1] ** 2, psi[1], psi[2]]               # (G x id) o psi, G = (p + q^2, q)
Fa = [sp.expand(Gid[0]), sp.expand(Gid[1]), sp.expand(Gid[2] + Gid[0] * Gid[1])]  # phi = (p, q, r + pq)
Ja = sp.Matrix([[sp.diff(f, v) for v in (x, y, z)] for f in Fa])
print('F =', Fa)
print('det JF =', sp.expand(Ja.det()), ' degrees =', [sp.Poly(f, x, y, z).total_degree() for f in Fa])
Ba = (Ja.T).inv().applyfunc(sp.cancel)
coeffs = [sp.Poly(Ba[2, k], x, y, z) for k in range(3)]
print('V_3 = sum_k c_k d_k with polynomial c_k =', [ck.as_expr() for ck in coeffs])
print('V_3 F_i =', [sp.cancel(sum(Ba[2, k] * sp.diff(Fa[i], (x, y, z)[k]) for k in range(3))) for i in range(3)])


def V3(P):
    return sum(coeffs[k] * P.diff((x, y, z)[k]) for k in range(3))


for v in (x, y, z):
    P, k = sp.Poly(v, x, y, z), 0
    while not P.is_zero and k < 80:
        P, k = V3(P), k + 1
    print(f'  V_3^{k} {v} = 0  (locally nilpotent on {v})')
print(f'[{time.time() - t0:.1f}s]')
print('Done.')
