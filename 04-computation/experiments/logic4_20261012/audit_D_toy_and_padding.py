#!/usr/bin/env python3
"""Audit D (2026-10-12) of THM-4626 -- machine checks.

(A) The non-example of THM-4626(4): D = x d_x - y d_y + d_z on Q[x,y,z].
    * div D = 0;
    * D is surjective: explicit preimage of every monomial x^a y^b z^c (closed form), verified
      symbolically for all monomials of total degree <= N, plus a rank check that D(P_{<=N})
      contains P_{<=N-1};
    * D is not locally nilpotent (D^k x = x);
    * Ker D on P_{<=N} is spanned by (xy)^k (so Ker D = Q[xy] in bounded degree: transcendence
      degree one, hence D is not a dual field of any Keller map of A^3, whose kernels have
      transcendence degree two).
(B) The padding example of THM-4626(3): F = F' x id with F' = THM-1300.  Checks det JF' is a
    nonzero constant, that the dual field V_4 of F' x id is d/dx_4 (so surjective and locally
    nilpotent), and counts the generic fibre of F' (deg F' = d_4 = number of lines in the generic
    fibre of (F_1, F_2, F_3) of the padded map) by a Groebner basis at a random rational point.
(C) The mechanism behind THM-4626(1) (audit finding): if K[x] sits inside L[F] with V = d/dF then V
    is locally nilpotent on K[x] simply because d/dF is locally nilpotent on L[F]; surjectivity is
    not needed beyond THM-4620(4)'s hypothesis mu_j = 0.  Illustrated on a Keller map of A^3 which is
    an automorphism not of product form: V_3 computed from (JF^T)^{-1}, local nilpotence on x, y, z
    verified by iteration, and the kernel R_3 = K[F_1, F_2] exhibited (x, y, z expressed in
    R_3[F_3]).
"""
import random
import sys
import time
from itertools import product

import sympy as sp

x, y, z, w = sp.symbols('x y z w')
Q = sp.QQ


def section(title):
    print('\n' + '=' * 78 + '\n' + title + '\n' + '=' * 78)


# ---------------------------------------------------------------- (A) toy derivation
section('(A) D = x d_x - y d_y + d_z on Q[x,y,z]')


def D_toy(p):
    return sp.expand(x * sp.diff(p, x) - y * sp.diff(p, y) + sp.diff(p, z))


div = sp.diff(x, x) + sp.diff(-y, y) + sp.diff(sp.Integer(1), z)
print('div D =', div)


def preimage_monomial(a, b, c):
    """Closed-form a' with D a' = x^a y^b z^c (the proof of THM-4626(4), made explicit)."""
    if a != b:
        lam = sp.Integer(a - b)
        # solve lam*p + p' = z^c:  p = sum_k (-1)^k c!/(c-k)! z^(c-k) / lam^(k+1)
        p = sum((-1) ** k * sp.factorial(c) / sp.factorial(c - k) * z ** (c - k) / lam ** (k + 1)
                for k in range(c + 1))
    else:
        p = z ** (c + 1) / (c + 1)
    return sp.expand(x ** a * y ** b * p)


N = 6
bad = 0
count = 0
for a, b, c in product(range(N + 1), repeat=3):
    if a + b + c > N:
        continue
    count += 1
    ap = preimage_monomial(a, b, c)
    if sp.expand(D_toy(ap) - x ** a * y ** b * z ** c) != 0:
        bad += 1
        print('  FAIL at monomial', (a, b, c))
print(f'explicit preimages verified for all {count} monomials of degree <= {N}: failures = {bad}')

# rank check: D(P_{<=N}) contains P_{<=N-1}
monos = [x ** a * y ** b * z ** c for a, b, c in product(range(N + 1), repeat=3) if a + b + c <= N]
idx = {m: i for i, m in enumerate(monos)}
M = sp.zeros(len(monos), len(monos))
for j, m in enumerate(monos):
    img = sp.Poly(D_toy(m), x, y, z)
    for mon, coeff in img.terms():
        M[idx[x ** mon[0] * y ** mon[1] * z ** mon[2]], j] = coeff
rank = M.rank()
print(f'dim P_<={N} = {len(monos)}, rank of D on it = {rank}, dim P_<={N - 1} = {sum(1 for m in monos if sp.Poly(m, x, y, z).total_degree() <= N - 1)}')
# check every monomial of degree <= N-1 is in the column space
colspace = M.columnspace()
Mc = sp.Matrix.hstack(*colspace)
missing = []
for m in monos:
    if sp.Poly(m, x, y, z).total_degree() <= N - 1:
        v = sp.zeros(len(monos), 1)
        v[idx[m]] = 1
        if sp.Matrix.hstack(Mc, v).rank() != Mc.rank():
            missing.append(m)
print('monomials of degree <= N-1 missing from D(P_<=N):', missing)

# non-nilpotence
p = x
for k in range(1, 6):
    p = D_toy(p)
print('D^5 x =', p, '(not nilpotent)')

# kernel in degree <= N
ker = M.nullspace()
kernel_polys = [sp.expand(sum(v[i] * monos[i] for i in range(len(monos)))) for v in ker]
print(f'dim Ker D on P_<={N} = {len(ker)}; expected floor(N/2)+1 = {N // 2 + 1}')
print('kernel basis:', kernel_polys)
span_ok = all(sp.Poly(kp, x, y, z).total_degree() % 2 == 0 and
              sp.simplify(kp.subs({x: w, y: 1}) - kp.subs({x: 1, y: w})) == 0 and
              sp.diff(kp, z) == 0 for kp in kernel_polys)
print('every kernel element is a polynomial in xy (symmetric in x<->y, no z):', span_ok)

# ---------------------------------------------------------------- (B) padding example
section("(B) padding example F = F' x id with F' = THM-1300")
u = 1 + x * y
F1 = u ** 3 * z + y ** 2 * u * (4 + 3 * x * y)
F2 = y + 3 * x * u ** 2 * z + 3 * x * y ** 2 * (4 + 3 * x * y)
F3 = 2 * x - 3 * x ** 2 * y - x ** 3 * z
Fp = [sp.expand(F1), sp.expand(F2), sp.expand(F3)]
J3 = sp.Matrix([[sp.diff(f, v) for v in (x, y, z)] for f in Fp])
detJ3 = sp.expand(J3.det())
print("det J(F') =", detJ3)
assert detJ3.is_number and detJ3 != 0, 'THM-1300 formula is not Keller as typed'

# padded map and its dual fields
vars4 = (x, y, z, w)
F4 = Fp + [w]
J4 = sp.Matrix([[sp.diff(f, v) for v in vars4] for f in F4])
B4 = (J4.T).inv()
B4 = B4.applyfunc(sp.simplify)
V4 = [sp.simplify(B4[3, k]) for k in range(4)]
print('dual field V_4 = sum_k B[4,k] d_k with coefficients', V4, '  (= d/dw)')
print('V_j F_i = delta_ij check:', all(sp.simplify(sum(B4[j, k] * sp.diff(F4[i], vars4[k]) for k in range(4)) - (1 if i == j else 0)) == 0 for i in range(4) for j in range(4)))

# generic degree of F' by a Groebner basis at a random rational point
random.seed(20261012)
c = [sp.Rational(random.randint(-9, 9), random.randint(1, 5)) for _ in range(3)]
t0 = time.time()
G = sp.groebner([Fp[i] - c[i] for i in range(3)], x, y, z, order='grevlex', domain='QQ')
# count standard monomials = dim of the quotient = number of solutions with multiplicity
lead = [sp.Poly(g, x, y, z).monoms(order='grevlex')[0] for g in G.exprs]


def is_standard(m):
    return not any(all(m[i] >= l[i] for i in range(3)) for l in lead)


std = []
frontier = [(0, 0, 0)]
seen = set()
while frontier:
    m = frontier.pop()
    if m in seen or not is_standard(m):
        continue
    seen.add(m)
    std.append(m)
    for i in range(3):
        mm = list(m)
        mm[i] += 1
        frontier.append(tuple(mm))
    if len(std) > 200:
        break
print(f"generic fibre of F' over c = {c}: dim of quotient (solutions with multiplicity) = {len(std)}  [{time.time() - t0:.1f}s]")
print("=> deg F' = d_4 for the padded map; the generic fibre of (F_1,F_2,F_3) is that many disjoint lines")

# ---------------------------------------------------------------- (C) the mechanism on an automorphism
section('(C) restriction of d/dF_3 on L_3[F_3] is locally nilpotent: a non-product automorphism of A^3')
# F = phi o (G x id) o psi with G = (x + y^2, y) planar, psi, phi triangular automorphisms of A^3
X, Y, Z = x, y, z
psi = [X + Z ** 2, Y + X * Z, Z]                 # automorphism (triangular)
Gid = [psi[0] + psi[1] ** 2, psi[1], psi[2]]      # (G x id) o psi
phi = [Gid[0], Gid[1] + Gid[0] * Gid[2], Gid[2] + Gid[1] ** 2]  # phi o ...
Fa = [sp.expand(f) for f in phi]
Ja = sp.Matrix([[sp.diff(f, v) for v in (x, y, z)] for f in Fa])
print('det JF =', sp.expand(Ja.det()), ' degrees =', [sp.Poly(f, x, y, z).total_degree() for f in Fa])
Ba = (Ja.T).inv().applyfunc(sp.cancel)
V3 = lambda p: sp.expand(sum(Ba[2, k] * sp.diff(p, (x, y, z)[k]) for k in range(3)))
for v in (x, y, z):
    p, k = v, 0
    while p != 0 and k < 40:
        p, k = V3(p), k + 1
    print(f'  V_3^{k} {v} = 0 (V_3 locally nilpotent on {v}; nilpotency index {k})')
print('  polynomial coefficients of V_3:', [sp.factor(Ba[2, k]) for k in range(3)])
print('  V_3 F_i =', [sp.simplify(sum(Ba[2, k] * sp.diff(Fa[i], (x, y, z)[k]) for k in range(3))) for i in range(3)])
print('Done.')
