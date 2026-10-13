"""Audit C of THM-4625: degrees d = 7, 8 (beyond the session's d <= 6). Exact symbolic identities (cheap expansions) for symbolic p;
the squarefreeness of h = b^2 - 8aP/w^2 is certified by a NONZERO discriminant at a random integer specialization of (a, b, p_i)
(a polynomial that is nonzero at one point is nonzero), which is what 'generic squarefreeness' needs. Also the generic-injectivity
hypothesis of (5): resultant nonzero at the same specialization.
"""
import sympy as sp, random
w, s, x, a, b, r1, r2 = sp.symbols('w s x a b r1 r2')
random.seed(7)
for d in (7, 8):
    ps = sp.symbols('p1:%d' % (d + 1))
    p = sum(ps[i - 1] * w**i for i in range(1, d + 1))
    q = sp.expand(sp.integrate((s / 2) * sp.diff(p, w).subs(w, s), (s, 0, w)))
    P = sp.expand(sp.integrate(p.subs(w, s), (s, 0, w)))
    l = ps[-1]
    F1 = (q + w * (sp.Symbol('c') * b - p) / 2) / sp.Symbol('c')**2
    c = sp.Symbol('c')
    ok1 = sp.expand(sp.diff(F1, w) - (c * b - p) / (2 * c**2)) == 0 and sp.simplify(sp.Poly(sp.expand(F1), w).LC() + l / (2 * (d + 1) * c**2)) == 0
    F2 = (p + 2 * (c**2 * a - q) / w) / c
    ok2 = sp.cancel(sp.diff(F2, w) + 2 * (c**2 * a - q) / (c * w**2)) == 0
    quad = sp.expand(q * (b * x - 2)**2 + w * p * (b * x - 2) - a * p**2 * x**2)
    A, B, C = [sp.expand(quad.coeff(x, k)) for k in (2, 1, 0)]
    ok3 = sp.expand(B**2 - 4 * A * C - p**2 * (b**2 * w**2 - 8 * a * (w * p - 2 * q))) == 0 and sp.expand(w * p - 2 * q - P) == 0
    h = sp.expand(b**2 - 8 * a * P / w**2)
    spec = {a: random.randint(1, 9), b: random.randint(1, 9)}; spec.update({pi: random.randint(-9, 9) or 1 for pi in ps}); spec[ps[-1]] = random.choice([1, 2, 3, -1])
    hs = sp.Poly(h.subs(spec), w)
    D = sp.discriminant(hs, w)
    ps_ = p.subs(spec); qs_ = q.subs(spec)
    e1 = sp.cancel((ps_.subs(w, r1) - ps_.subs(w, r2)) / (r1 - r2)); e2 = sp.cancel((qs_.subs(w, r1) - qs_.subs(w, r2)) / (r1 - r2))
    res = sp.resultant(sp.Poly(e1, r2), sp.Poly(e2, r2))
    print("d=%d: (1) %s (2) %s (3) disc identity & wp-2q=P %s | deg h = %d, disc(h) at a random integer specialization = %s (nonzero: %s) => generic squarefreeness, genus %d | resultant (5) at the specialization nonzero: %s"
          % (d, ok1, ok2, ok3, hs.degree(), D, D != 0, (hs.degree() - 1) // 2, sp.expand(res) != 0), flush=True)
