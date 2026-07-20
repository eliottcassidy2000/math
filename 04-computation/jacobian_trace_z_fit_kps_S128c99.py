#!/usr/bin/env python3
"""
kind-pasteur-2026-07-20-S128c99 -- Tr(z) by exact specialization + weight-basis fit.
Weight/parity constraints: wt(Tr z) = -2 under (a,b,g) -> (l^-2 a, l^-1 b, l g),
tau-even under (b,g) -> (-b,-g).  Candidate basis: {a, b^2, abg, b^3g, a^2g^2,
b^4g^2, ab^2g^2?...}; fallback: L*Tr(z) polynomial (pole case).
Known anchor: a-axis Tr(z) = -51a? NO WAIT: collision fiber z-sum = -1/4+13/2+13/2
= 51/4 at a = -1/4  =>  if linear in a on the axis: Tr(z) = -51a there.
"""
import sympy as sp
from sympy import Rational, symbols, Poly
import random

x, s, a, b, g = symbols('x s a b g')
L = 27*a**2*g**2 - 18*a*b*g + 16*a + b**3*g - b**2
N = sp.expand(L*x**3 + (4 - 3*b*g)*x - 2*g)
C1 = 3*a*x**2 - b*s*x - b*x + s**2 + s
C2 = sp.expand(a*x**3 + g*(1+s)**3 - x*(s**2 + 3*s + 2))
prs = sp.subresultants(Poly(C1, s), Poly(C2, s))
lin1 = [q for q in prs if sp.degree(q, s) == 1]
q1 = sp.Poly(lin1[-1], s)
s_of_x = sp.cancel(-q1.all_coeffs()[1] / q1.all_coeffs()[0])

# z = (a - B1)/u^3 with B1 = y^2 u (4+3xy) = (s^2/x^2)(1+s)(4+3s)
B1x = (s**2/x**2)*(1+s)*(4+3*s)
z_of = sp.cancel(sp.together((a - B1x.subs(s, s_of_x)) / (1 + s_of_x)**3))
zn, zd = sp.fraction(z_of)

random.seed(11)
samples = []
tries = 0
while len(samples) < 16 and tries < 200:
    tries += 1
    va = Rational(random.randint(-9, 9), random.choice([1, 2, 3, 4]))
    vb = Rational(random.randint(-9, 9), random.choice([1, 2, 3]))
    vg = Rational(random.randint(-9, 9), random.choice([1, 2, 3]))
    if va == 0 or vg == 0: continue
    sub = {a: va, b: vb, g: vg}
    Lv = L.subs(sub)
    if Lv == 0: continue
    Nv = Poly(N.subs(sub), x, domain='QQ')
    znv = Poly(sp.expand(zn.subs(sub)), x, domain='QQ')
    zdv = Poly(sp.expand(zd.subs(sub)), x, domain='QQ')
    if sp.degree(sp.gcd(zdv, Nv), x) > 0: continue   # denominator not invertible
    ginv = sp.invert(zdv, Nv)
    red = (znv * ginv) % Nv
    cs = red.all_coeffs()[::-1]
    while len(cs) < 3: cs.append(0)
    # Newton: p1 = 0, p2 = -2*(4-3bg)/L
    p2 = -2*(4 - 3*vb*vg)/Lv
    tr = 3*cs[0] + cs[2]*p2
    samples.append((va, vb, vg, sp.nsimplify(tr)))
print(f"collected {len(samples)} samples", flush=True)
for va, vb, vg, t in samples[:4]:
    print(f"  target ({va},{vb},{vg}): Tr(z) = {t}", flush=True)

# fit polynomial: weight -2, tau-even basis
basis = [a, b**2, a*b*g, b**3*g, a**2*g**2, b**4*g**2, a*b**2*g**2]
coeffs = sp.symbols('c0:%d' % len(basis))
eqs = []
for va, vb, vg, t in samples:
    eqs.append(sum(c*m.subs({a: va, b: vb, g: vg}) for c, m in zip(coeffs, basis)) - t)
sol = sp.solve(eqs, coeffs, dict=True)
if sol:
    expr = sum(sol[0].get(c, c)*m for c, m in zip(coeffs, basis))
    print(f"\nPOLYNOMIAL FIT SUCCEEDS: Tr(z) = {sp.expand(expr)}", flush=True)
    resid = [sp.simplify(expr.subs({a: va, b: vb, g: vg}) - t) for va, vb, vg, t in samples]
    print(f"  residuals all zero: {all(r == 0 for r in resid)}", flush=True)
else:
    print("\npolynomial fit FAILS -- trying L*Tr(z) polynomial (pole case), wider basis", flush=True)
    # L*Tr(z): weight -2 + wt(L); wt(L): monomials a^2g^2 wt -2, abg -2, a -2, b^3g -2, b^2 -2 => L wt -2 => product wt -4, tau-even
    basis2 = [a**2, a*b**2, b**4, a**2*b*g, a*b**3*g, b**5*g, a**3*g**2, a**2*b**2*g**2, a*b**4*g**2, b**6*g**2, a**3*b*g**3, a**2*b**3*g**3, a**4*g**4]
    co2 = sp.symbols('d0:%d' % len(basis2))
    eqs2 = []
    for va, vb, vg, t in samples:
        Lv = L.subs({a: va, b: vb, g: vg})
        eqs2.append(sum(c*m.subs({a: va, b: vb, g: vg}) for c, m in zip(co2, basis2)) - t*Lv)
    sol2 = sp.solve(eqs2, co2, dict=True)
    if sol2:
        expr2 = sum(sol2[0].get(c, c)*m for c, m in zip(co2, basis2))
        print(f"L*Tr(z) = {sp.expand(expr2)}", flush=True)
        print(f"=> Tr(z) = [{sp.factor(expr2)}] / L   -- L-POLE CONFIRMED", flush=True)
        resid = [sp.simplify(expr2.subs({a: va, b: vb, g: vg}) - t*L.subs({a: va, b: vb, g: vg})) for va, vb, vg, t in samples]
        print(f"  residuals all zero: {all(r == 0 for r in resid)}", flush=True)
    else:
        print("both fits fail -- need bigger basis / check", flush=True)
print("\nDONE.", flush=True)
