#!/usr/bin/env python3
"""Symbolic audit of the algebra in THM-4671 (Lemma 1, Lemma 2, Proof of 1, factorization of A,
zero set of A, Ising identification).  Independent of the author's code.  Uses sympy."""
import sympy as sp

N, c, p, q = sp.symbols('N c p q')

# ---- Lemma 1 from the two linear conditions (re-derived here, see REPORT.md) ----------
# L_2 v = 0, coefficient of L_{-2}^{N-1}:  4N(N-1) + N c/2 + 6 alpha = 0
alpha = sp.solve(sp.Eq(4 * N * (N - 1) + N * c / 2 + 6 * sp.Symbol('a'), 0), sp.Symbol('a'))[0]
# L_1 v = 0, coefficient of L_{-3}L_{-2}^{N-2}:  3N(N-1)/2 + 5 alpha + 8 beta = 0
beta = sp.solve(sp.Eq(sp.Rational(3, 2) * N * (N - 1) + 5 * alpha + 8 * sp.Symbol('b'), 0), sp.Symbol('b'))[0]
print('alpha =', sp.factor(alpha), '  claimed -N(8(N-1)+c)/12 :',
      sp.simplify(alpha + N * (8 * (N - 1) + c) / 12) == 0)
print('5alpha+8beta =', sp.factor(5 * alpha + 8 * beta))

# ---- Lemma 2 -------------------------------------------------------------------------
s1 = N * (N + 1) / 2 + 3 * alpha
s2 = 4 * beta - 2 * alpha
s3 = -3 * beta
# d^3 f coefficients on m1, m2, m3
d1, d2, d3 = 6 * N, 6 * N * (N - 1), N * (N - 1) * (N - 2)
A, B = sp.symbols('A B')
solB = sp.solve(sp.Eq(s1, B * d1), B)[0]
solA = sp.solve(sp.Eq(s2, A + solB * d2), A)[0]
resid3 = sp.simplify(s3 - (solA * (2 * N - 7) / 12 + solB * d3))
print('B =', sp.factor(solB), ' claimed -(6N+c-10)/24 :', sp.simplify(solB + (6 * N + c - 10) / 24) == 0)
cc2 = 1 - 3 * (2 * N - 1) ** 2 / (2 * N + 1)
print('A =', sp.factor(solA), ' claimed N(2N+1)(c-c_{2,2N+1})/8 :',
      sp.simplify(solA - N * (2 * N + 1) * (c - cc2) / 8) == 0)
print('m3 residual (must be 0 identically):', resid3)
print('c_{2,2N+1} =', sp.factor(cc2), ' claimed -(12N^2-14N+2)/(2N+1):',
      sp.simplify(cc2 + (12 * N ** 2 - 14 * N + 2) / (2 * N + 1)) == 0)

# ---- A in terms of (p,q) -------------------------------------------------------------
cpq = 1 - 6 * (p - q) ** 2 / (p * q)
Apq = sp.factor(solA.subs({N: (p - 1) * (q - 1) / 2, c: cpq}))
claimed = 3 * (p - 1) * (p - 2) * (q - 1) * (q - 2) * (p ** 2 * q ** 2 - (p - q) ** 2) / (16 * p * q)
print('A(p,q) factored =', Apq)
print('claimed factorization equal:', sp.simplify(Apq - claimed) == 0)
print('p^2q^2-(p-q)^2 =', sp.factor(p ** 2 * q ** 2 - (p - q) ** 2))
# zero set on 2 <= p < q coprime, q <= 60
zeros = [(pp, qq) for pp in range(2, 60) for qq in range(pp + 1, 61) if sp.gcd(pp, qq) == 1
         and claimed.subs({p: pp, q: qq}) == 0]
print('zeros of A for coprime 2<=p<q<=60 have p=2 only:', all(z[0] == 2 for z in zeros),
      ' count', len(zeros), ' p=2 count', len([1 for qq in range(3, 61, 2)]))
# Is A symmetric in p,q?  (relevant to the Drinfeld-Sokolov remark in HYP-9293)
print('A symmetric in p<->q:', sp.simplify(claimed - claimed.subs({p: q, q: p}, simultaneous=True)) == 0)

# ---- Ising identification ------------------------------------------------------------
x0, x1, x2, x3 = sp.symbols('x0:4')
f = x0 ** 3
dd = {x0: x1, x1: 2 * x2, x2: 3 * x3, x3: 4 * sp.Symbol('x4')}


def d(expr):
    return sp.expand(sum(sp.diff(expr, v) * dd[v] for v in (x0, x1, x2, x3)))


d3f = d(d(d(f)))
b3 = x0 * x1 * x2 + sp.Rational(2 * 3 - 7, 12) * x1 ** 3
aeh = x3 * x0 ** 2 / 6 + x2 * x1 * x0
print('Ising: AEH b == d^3f/108 + (2/3) b_3 :', sp.expand(aeh - d3f / 108 - sp.Rational(2, 3) * b3) == 0)
# general d^3 f formula
Nn = 7
xs = sp.symbols('y0:6')
ddy = {xs[j]: (j + 1) * xs[j + 1] for j in range(5)}


def dy(expr):
    return sp.expand(sum(sp.diff(expr, v) * ddy[v] for v in xs[:5]))


fy = xs[0] ** Nn
lhs = dy(dy(dy(fy)))
rhs = 6 * Nn * xs[0] ** (Nn - 1) * xs[3] + 6 * Nn * (Nn - 1) * xs[0] ** (Nn - 2) * xs[1] * xs[2] + \
    Nn * (Nn - 1) * (Nn - 2) * xs[0] ** (Nn - 3) * xs[1] ** 3
print('d^3 f formula (N=7):', sp.expand(lhs - rhs) == 0)
