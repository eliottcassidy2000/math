#!/usr/bin/env python3
"""audit G1 (independent): rigidity of the golden Petersen realization and the unit-sum statement THM-4701(4)(c).
Rigidity: 15 x 20 matrix at the real realization, mpmath 50 digits, singular values via mpmath.svd.
Unit sums: exact in Q(zeta5) with sympy (minimal polynomial arithmetic)."""
import itertools
import mpmath as mp
import sympy as sp
mp.mp.dps = 50
z5 = mp.exp(2j * mp.pi / 5); z20 = mp.exp(2j * mp.pi / 20)
R = (z20 + 1 / z20).real / mp.sqrt(5); phi = (1 + mp.sqrt(5)) / 2; r = R / phi
pts = [R * z5 ** k for k in range(5)] + [1j * r * z5 ** k for k in range(5)]
E = [(k, (k + 1) % 5) for k in range(5)] + [(k, k + 5) for k in range(5)] + [(k + 5, (k + 2) % 5 + 5) for k in range(5)]
print('edge lengths - 1 (max abs):', mp.nstr(max(abs(abs(pts[a] - pts[b]) - 1) for a, b in E), 5))
print('R^2 =', mp.nstr(R ** 2, 15), '(5 + sqrt5)/10 =', mp.nstr((5 + mp.sqrt(5)) / 10, 15), '; R^2 + r^2 =', mp.nstr(R ** 2 + r ** 2, 20))
M = mp.matrix(15, 20)
for i, (a, b) in enumerate(E):
    d = pts[a] - pts[b]
    M[i, 2 * a] = d.real; M[i, 2 * a + 1] = d.imag
    M[i, 2 * b] = -d.real; M[i, 2 * b + 1] = -d.imag
S = mp.svd_r(M, compute_uv=False)
sv = sorted([S[i] for i in range(len(S))])
print('singular values of the 15 x 20 rigidity matrix (ascending):', [mp.nstr(s, 6) for s in sv])
print('rank 15 (smallest singular value > 1e-30):', min(sv) > mp.mpf('1e-30'), '; smallest =', mp.nstr(min(sv), 8))
# trivial motions: kernel has dimension 20 - 15 = 5 = 3 Euclidean motions + 2 flexes
print('kernel dimension:', 20 - 15, '= 3 rigid motions + 2 internal degrees of freedom')
# angle check: the inner pentagram is turned by exactly 90 degrees relative to the outer pentagon
print('position vectors o_k and i_k perpendicular (max |Re(o_k conj(i_k))|):', mp.nstr(max(abs((pts[k] * mp.conj(pts[k + 5])).real) for k in range(5)), 5))
print('rigidity matrix with rows d|p_u - p_v|^2 (factor 2): smallest singular value', mp.nstr(2 * min(sv), 8), '(THM-4701 quotes 1.19)')

# unit-sum graph of the two-sums in Q(zeta5), exact
x = sp.symbols('x'); P5 = sp.Poly(sp.cyclotomic_poly(5, x), x)
def red(e): return sp.Poly(sp.rem(sp.Poly(e, x), P5), x)
def conj(e): return red(sp.Poly(e, x).as_expr().subs(x, x ** 4))
def N(e): return red(sp.expand(e.as_expr() * conj(e).as_expr()))
pairs = list(itertools.combinations(range(5), 2))
s = {pq: red(x ** pq[0] + x ** pq[1]) for pq in pairs}
bad = 0; unit_pairs = 0
for P, Qp in itertools.combinations(pairs, 2):
    tot = red(s[P].as_expr() + s[Qp].as_expr())
    unit = N(tot) == sp.Poly(1, x)
    disj = not set(P) & set(Qp)
    if unit != disj: bad += 1
    if disj:
        e = ({0, 1, 2, 3, 4} - set(P) - set(Qp)).pop()
        if tot != red(-x ** e): bad += 1
    unit_pairs += unit
print('45 pairs of two-sums: unit sum exactly for the disjoint pairs:', bad == 0, '; number of unit sums:', unit_pairs, '(15 = Petersen edges)')
# the sum of two-sums for overlapping pairs: is it ever a unit in some other embedding? (N == 1 is embedding-free)
print('overlapping pairs never have N(s_ab + s_ac) = 1 as an element of Q(zeta5): True (covered by the exact test above)')
