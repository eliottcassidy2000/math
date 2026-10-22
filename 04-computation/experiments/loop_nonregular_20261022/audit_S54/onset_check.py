"""Auditor: HYP-9369 section 5 'onset 19' algebra (conditional on C = 0), recomputed with sympy from the formula
G = z^5/2 (1/p^2 - 1/u) + R6/u + z^10 (1 + z + z^3)/u^2, p = 1 - z - z^3, u = p(z^2), R6 = z^5 + z^6 + z^7."""
import sympy as sp
z = sp.symbols('z')
p = 1 - z - z**3
u = p.subs(z, z**2)
R6 = z**5 + z**6 + z**7
G = z**5 / 2 * (1 / p**2 - 1 / u) + R6 / u + z**10 * (1 + z + z**3) / u**2
G = sp.cancel(sp.together(G))
num, den = sp.fraction(G)
num = sp.Poly(sp.expand(num), z); den = sp.Poly(sp.expand(den), z)
c = den.coeffs()[-1] if den.eval(0) != 0 else 1
lc = den.eval(0)
num = sp.Poly(num.as_expr() / lc, z); den = sp.Poly(den.as_expr() / lc, z)
print('reduced numerator:', sp.expand(num.as_expr()))
print('reduced denominator degree', den.degree(), '; equals p^2 u^2:', sp.expand(den.as_expr() - sp.expand(p**2 * u**2)) == 0)
print('gcd(num, den) =', sp.gcd(num, den))
ser = sp.series(G, z, 0, 60).removeO()
g = [ser.coeff(z, k) for k in range(60)]
print('g_5..g_16 =', g[5:17])
D = sp.Poly(sp.expand(p**2 * u**2), z).all_coeffs()[::-1]   # D(z) = sum d_j z^j, d_0 = 1
def resid(k):
    return sum(D[j] * g[k - j] for j in range(len(D)) if k - j >= 0)
print('residual of the D-recurrence at k = 15..25:', [resid(k) for k in range(15, 26)])
