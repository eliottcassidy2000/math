#!/usr/bin/env python3
"""Audit of THM-4870 (A): exact mean of the run length r over UH classes of order n.

Independent re-derivation from the generating function of THM-4868 (H):
    G(x,t) = t^2 x^4 (1 + (t^2 - t) x) / (1 - (1 + 2t) x + t x^2).
We (1) expand N(n,r) by an explicit recurrence obtained from the denominator,
(2) compute E[r] exactly as a rational, (3) compare with the claimed closed form
    ((5 + sqrt5) n + 12 - 8 sqrt5)/10,
(4) derive the exact error term by partial fractions and test the claim O(phi^(-4n)).
"""
from fractions import Fraction
import sympy as sp
import mpmath as mp

mp.mp.dps = 80
x, t = sp.symbols('x t')
G = t**2 * x**4 * (1 + (t**2 - t) * x) / (1 - (1 + 2*t) * x + t * x**2)

# --- (1) coefficients by series expansion in x (polynomials in t)
NMAX = 80
ser = sp.series(G, x, 0, NMAX + 1).removeO()
ser = sp.expand(ser)
rows = {}
for n in range(4, NMAX + 1):
    c = sp.expand(ser.coeff(x, n))
    poly = sp.Poly(c, t)
    rows[n] = {int(m[0]): int(coef) for m, coef in zip(poly.monoms(), poly.coeffs())}

# Compare with the published table of THM-4868 note section 4 (rows n = 5..13)
published = {
    5: [1, 1, 1], 6: [1, 2, 3, 2], 7: [1, 3, 6, 7, 4], 8: [1, 4, 10, 16, 16, 8],
    9: [1, 5, 15, 30, 41, 36, 16], 10: [1, 6, 21, 50, 85, 102, 80, 32],
    11: [1, 7, 28, 77, 155, 231, 248, 176, 64],
    12: [1, 8, 36, 112, 259, 456, 608, 592, 384, 128],
    13: [1, 9, 45, 156, 406, 819, 1289, 1560, 1392, 832, 256],
}
ok_rows = True
for n, lst in published.items():
    mine = [rows[n].get(r, 0) for r in range(2, n)]
    if mine != lst:
        ok_rows = False
        print("ROW MISMATCH n=%d mine=%s published=%s" % (n, mine, lst))
print("published N(n,r) rows n=5..13 reproduced from the GF:", ok_rows)
print("row n=4:", rows[4])

# Fibonacci check and extremes
def fib(k):
    a, b = 0, 1
    for _ in range(k):
        a, b = b, a + b
    return a
ok_fib = all(sum(rows[n].values()) == fib(2 * n - 6) for n in range(4, NMAX + 1))
print("sum_r N(n,r) = F(2n-6) for n=4..%d:" % NMAX, ok_fib)
ok_ext = all(rows[n].get(2, 0) == 1 and rows[n].get(n - 1, 0) == 2 ** (n - 5) for n in range(5, NMAX + 1))
print("N(n,2)=1 and N(n,n-1)=2^(n-5) for n=5..%d:" % NMAX, ok_ext)
ok_range = all(min(rows[n]) >= 2 and max(rows[n]) <= n - 1 for n in range(5, NMAX + 1))
print("support of r is within [2, n-1] for n=5..%d:" % NMAX, ok_range)

# --- (2) exact means
def mean_r(n):
    tot = sum(rows[n].values())
    s = sum(r * c for r, c in rows[n].items())
    return Fraction(s, tot)

s5 = mp.sqrt(5)
phi = (1 + s5) / 2
def claimed(n):
    return ((5 + s5) * n + 12 - 8 * s5) / 10

# --- (3)/(4) exact error, scaled by phi^(4n) and by n*phi^(-4n)
print("\n n   E[r] (exact rational -> float)   claimed          err=E-claimed      err*phi^(4n)       err*phi^(4n)/n")
errs = []
for n in list(range(4, 21)) + [30, 40, 50, 60, 70, 80]:
    m = mean_r(n)
    mv = mp.mpf(m.numerator) / m.denominator
    e = mv - claimed(n)
    sc = e * phi ** (4 * n)
    errs.append((n, e, sc))
    print("%3d  %-30s %-16s %-18s %-18s %-18s" % (
        n, mp.nstr(mv, 15), mp.nstr(claimed(n), 15), mp.nstr(e, 8), mp.nstr(sc, 10), mp.nstr(sc / n, 10)))

# --- exact closed form via partial fractions at t=1
dG = sp.simplify(sp.diff(G, t).subs(t, 1))
G1 = sp.simplify(G.subs(t, 1))
print("\nG(x,1)        =", sp.factor(G1))
print("d/dt G(x,1)   =", sp.factor(dG))
claimed_dG = x**4 * (x**3 - 2*x**2 - 3*x + 2) / (1 - 3*x + x**2)**2
print("matches the claimed x^4 (x^3 - 2x^2 - 3x + 2)/(1 - 3x + x^2)^2:", sp.simplify(dG - claimed_dG) == 0)

# closed form of a_n = [x^n] dG and b_n = [x^n] G1 via sympy rsolve-free approach:
# a_n = (A n + B) phi^(2n) + (A' n + B') phi^(-2n) for n >= n0; fit exactly in Q(sqrt5)
n_ = sp.symbols('n', integer=True)
S5 = sp.sqrt(5)
PHI = (1 + S5) / 2
A, B, A2, B2 = sp.symbols('A B A2 B2')
def a_exact(n):
    return sum(r * c for r, c in rows[n].items())
eqs = [sp.Eq((A * k + B) * PHI**(2*k) + (A2 * k + B2) * PHI**(-2*k), a_exact(k)) for k in range(10, 14)]
sol = sp.solve(eqs, [A, B, A2, B2], dict=True)[0]
sol = {kk: sp.nsimplify(sp.radsimp(sp.simplify(vv)), [S5]) for kk, vv in sol.items()}
print("a_n = (A n + B) phi^(2n) + (A2 n + B2) phi^(-2n) with", sol)
ok_cf = all(sp.simplify((sol[A]*k + sol[B]) * PHI**(2*k) + (sol[A2]*k + sol[B2]) * PHI**(-2*k) - a_exact(k)) == 0
            for k in range(6, 40))
print("closed form for a_n valid n=6..39:", ok_cf)
# b_n = F(2n-6) = (phi^(2n-6) - phi^(-(2n-6)))/sqrt5
C1 = PHI**(-6) / S5
C2 = -PHI**(6) / S5
# E[r] = [(A n + B) + (A2 n + B2) q] / [C1 + C2 q],  q = phi^(-4n)
q = sp.symbols('q')
Er = ((sol[A]*n_ + sol[B]) + (sol[A2]*n_ + sol[B2]) * q) / (C1 + C2 * q)
lin = sp.simplify(sp.series(Er, q, 0, 2).removeO())
lead = sp.radsimp(sp.simplify(lin.subs(q, 0)))
first = sp.radsimp(sp.simplify(sp.diff(lin, q)))
print("leading term (q^0) =", sp.simplify(sp.expand(lead)))
print("claimed            =", sp.simplify(((5 + S5) * n_ + 12 - 8 * S5) / 10))
print("difference         =", sp.simplify(sp.expand(lead - ((5 + S5) * n_ + 12 - 8 * S5) / 10)))
print("coefficient of q = phi^(-4n) in E[r] (exact):", sp.simplify(sp.expand(first)))
