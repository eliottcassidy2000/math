"""Independent replay of the deficit-six degree-three-node series (Session30 INTEGRATION.md section 3), conditional on
the user-reported closure C = 0, i.e. R_6 = z^5 + z^6 + z^7.  Session opus-2026-10-22-S54.
G = z^5/2 (1/p^2 - 1/u) + R_6/u + z^10 (1 + z + z^3)/u^2,  p = 1 - z - z^3,  u = p(z^2).
Checks: reduced denominator, numerator, the sharp onset of the denominator recurrence, first coefficients.
Usage: python deficit_six_onset.py"""
import sympy as sp

z = sp.symbols('z')
p = 1 - z - z ** 3
u = sp.expand(p.subs(z, z ** 2))
R6 = z ** 5 + z ** 6 + z ** 7
G = sp.Rational(1, 2) * z ** 5 * (1 / p ** 2 - 1 / u) + R6 / u + z ** 10 * (1 + z + z ** 3) / u ** 2
N, D = sp.fraction(sp.cancel(sp.together(G)))
N, D = sp.Poly(N, z), sp.Poly(D, z)
lc = D.coeffs()[-1]                                  # normalise D(0) = 1
N, D = N * (1 / D.eval(0)), D * (1 / D.eval(0))
print('deg N =', N.degree(), ' deg D =', D.degree(), ' gcd =', sp.gcd(N.as_expr(), D.as_expr()))
print('D == p^2 u^2:', sp.expand(D.as_expr() - sp.expand(p ** 2 * u ** 2)) == 0)
print('N =', sp.expand(N.as_expr()))
M = 80
g = sp.Poly(sp.series(G, z, 0, M).removeO(), z)
coef = [g.coeff_monomial(z ** k) for k in range(M)]
d = [D.coeff_monomial(z ** k) for k in range(D.degree() + 1)]
bad = [k for k in range(D.degree(), M) if sum(d[j] * coef[k - j] for j in range(len(d))) != 0]
print('coefficients g_5..g_24 =', coef[5:25])
print('denominator recurrence sum_j d_j g_(k-j) = 0 fails at k =', bad, '-> holds from k =', (max(bad) + 1) if bad else D.degree())
