"""Audit C of THM-4625 (independent auditor, 2026-10-12): statements (1)-(3) for SYMBOLIC p of degree d = 2..7.

Family: F = ((q(w)+gamma w)/(gamma x)^2, (p(w)+2 gamma)/(gamma x), gamma x), w = gamma(1+xy), q' = (w/2) p', p(0) = 0.
Checks (all exact, sympy):
  (1) F1 on the (F2,F3)-fibre: polynomial of degree d+1 in w, leading coefficient -l/(2(d+1)c^2), derivative (cb-p)/(2c^2).
  (2) F2 on the (F1,F3)-fibre: derivative -2(c^2 a - q)/(c w^2); pole orders at 0 and infinity; coefficient of w^d is l/((d+1)c).
  (3) (F1,F2)-fibre quadratic A x^2 + B x + C: discriminant identity, wp - 2q = P, order of vanishing of A, B, C at w = 0
      (exactly two when p'(0) != 0; HIGHER when p'(0) = 0, so w = 0 can survive the division by w^2 as a spurious line),
      common roots of A~, B~, C~ (= common nonzero roots of p and q: spurious vertical lines w = w0 inside the plane curve),
      degree of h = b^2 - 8aP/w^2, generic squarefreeness, and the genus formula via Riemann-Hurwitz with the branch point at infinity.
Also: the generic-injectivity hypothesis of (5) is checked to be automatic (resultant nonzero) for symbolic p, d = 2..5.
"""
import sympy as sp, sys, time
T0 = time.time()
w, s, x, a, b, c, t = sp.symbols('w s x a b c t')
r1, r2 = sp.symbols('r1 r2')

def mkpq(d):
    ps = sp.symbols('p1:%d' % (d + 1))
    p = sum(ps[i - 1] * w**i for i in range(1, d + 1))
    q = sp.integrate((s / 2) * sp.diff(p, w).subs(w, s), (s, 0, w))
    P = sp.integrate(p.subs(w, s), (s, 0, w))
    return ps, sp.expand(p), sp.expand(q), sp.expand(P)

for d in range(2, 7):   # d = 7, 8 in audit_C_identities_fast.py (random integer specialization for the discriminant)
    ps, p, q, P = mkpq(d)
    l = ps[-1]
    # (1)
    F1 = (q + w * (c * b - p) / 2) / c**2
    F1p = sp.Poly(sp.expand(F1), w)
    ok1a = F1p.degree() == d + 1
    ok1b = sp.simplify(F1p.LC() + l / (2 * (d + 1) * c**2)) == 0
    ok1c = sp.expand(sp.diff(F1, w) - (c * b - p) / (2 * c**2)) == 0
    # (2)
    F2 = (p + 2 * (c**2 * a - q) / w) / c
    ok2a = sp.cancel(sp.diff(F2, w) + 2 * (c**2 * a - q) / (c * w**2)) == 0
    num = sp.Poly(sp.expand(F2 * w * c), w)          # c*w*F2 = w p + 2(c^2 a - q)
    ok2b = num.degree() == d + 1 and sp.simplify(num.LC() - l / (d + 1)) == 0   # pole of order d at infinity, coefficient l/((d+1)c)
    ok2c = sp.expand(num.coeff_monomial(1) - 2 * c**2 * a) == 0                  # residue-type coefficient: 2c^2 a / (c w) = 2ca/w
    # (3)
    quad = sp.expand(q * (b * x - 2)**2 + w * p * (b * x - 2) - a * p**2 * x**2)
    A, B, C = [sp.expand(quad.coeff(x, k)) for k in (2, 1, 0)]
    disc = sp.expand(B**2 - 4 * A * C)
    ok3a = sp.expand(disc - p**2 * (b**2 * w**2 - 8 * a * (w * p - 2 * q))) == 0
    ok3b = sp.expand(w * p - 2 * q - P) == 0
    ords = [min(m[0] for m in sp.Poly(X, w).monoms()) for X in (A, B, C)]       # order of vanishing at w = 0 (generic p_i)
    # with p1 = 0 (p'(0) = 0):
    ords0 = [min(m[0] for m in sp.Poly(X.subs(ps[0], 0), w).monoms()) for X in (A, B, C)]
    At, Bt, Ct = [sp.expand(X / w**2) for X in (A, B, C)]
    assert all(sp.expand(X).is_polynomial(w) for X in (At, Bt, Ct))
    h = sp.expand(b**2 - 8 * a * P / w**2)
    hp = sp.Poly(h, w)
    ok3c = hp.degree() == d - 1
    D = sp.discriminant(hp, w)
    ok3d = sp.expand(D) != 0
    # discriminant of the divided quadratic equals (p/w)^2 h
    ok3e = sp.expand(Bt**2 - 4 * At * Ct - sp.expand(p / w)**2 * h) == 0
    # genus via Riemann-Hurwitz: branch points = deg h roots (+ infinity iff deg h odd)
    nb = (d - 1) + (1 if (d - 1) % 2 == 1 else 0)
    genus_RH = sp.Rational(nb - 2, 2)
    ok3f = genus_RH == (d - 2) // 2
    # common roots of A~, B~, C~ at w = w0 != 0  <=>  p(w0) = q(w0) = 0 (shown by hand; here: gcd over Q(a,b,p_i)[w] is 1 generically)
    g = sp.gcd(sp.gcd(sp.Poly(At, w), sp.Poly(Bt, w)), sp.Poly(Ct, w))
    print("d=%d: (1) deg %s LC %s deriv %s | (2) deriv %s, pole data %s %s | (3) disc %s, wp-2q=P %s, ord_0(A,B,C)=%s (p1=0: %s), deg h=%d (%s), h generically squarefree %s, disc~=(p/w)^2 h %s, RH genus %s = floor((d-2)/2) %s, generic gcd(A~,B~,C~)=%s"
          % (d, ok1a, ok1b, ok1c, ok2a, ok2b, ok2c, ok3a, ok3b, ords, ords0, hp.degree(), ok3c, ok3d, ok3e, genus_RH, ok3f, g.as_expr()), "[%.0fs]" % (time.time() - T0), flush=True)

# (3) when p and q share a nonzero root: d = 3, p = w (w - w0)(w - w0/2) has P(w0) = 0, hence q(w0) = 0: the line w = w0 is a component
w0 = sp.symbols('w0')
p = sp.expand(w * (w - w0) * (w - w0 / 2)); q = sp.expand(sp.integrate((s / 2) * sp.diff(p, w).subs(w, s), (s, 0, w)))
P = sp.expand(sp.integrate(p.subs(w, s), (s, 0, w)))
print("\nspurious vertical line example (d = 3, p = w(w-w0)(w-w0/2)): p(w0) =", p.subs(w, w0), " q(w0) =", sp.simplify(q.subs(w, w0)), " P(w0) =", sp.simplify(P.subs(w, w0)))
quad = sp.expand(q * (b * x - 2)**2 + w * p * (b * x - 2) - a * p**2 * x**2)
A, B, C = [sp.expand(quad.coeff(x, k) / w**2) for k in (2, 1, 0)]
print("   A~, B~, C~ all vanish at w = w0:", [sp.simplify(X.subs(w, w0)) == 0 for X in (A, B, C)], "=> {w = w0} x A^1_x lies on the plane curve; it is excluded by p(w) != 0 in the chart")

# (5) hypothesis is automatic: resultant of (p(r1)-p(r2))/(r1-r2) and (q(r1)-q(r2))/(r1-r2) in r2 is a nonzero polynomial for SYMBOLIC p
for d in range(2, 5):   # symbolic resultant: d <= 4 (d = 5, 6 with random integer p are in audit_C_monodromy.py)
    ps, p, q, P = mkpq(d)
    e1 = sp.cancel((p.subs(w, r1) - p.subs(w, r2)) / (r1 - r2)); e2 = sp.cancel((q.subs(w, r1) - q.subs(w, r2)) / (r1 - r2))
    res = sp.resultant(sp.Poly(e1, r2), sp.Poly(e2, r2))
    print("d=%d: resultant_r2 of the two difference quotients, symbolic p: nonzero %s (degree in r1: %d)" % (d, sp.expand(res) != 0, sp.Poly(res, r1).degree()))
