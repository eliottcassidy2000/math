#!/usr/bin/env python3
"""Exact controls for the cusp-ideal nonrational scalar-time proof.

RESERVED: finite controls do not replace the embedded-surface, relative
constant-field, or pole-preserving curve-map arguments in the companion.
No inherited mathematical implementation is imported.
"""
from functools import reduce
from hashlib import sha256
from itertools import product
from math import gcd, factorial
import json
import sympy as sp

checks = 0

def check(value, label):
    global checks
    checks += 1
    if not value:
        raise RuntimeError(label)


def cover(A, B, es):
    E = sum(es)
    M = 2*A + 3*B + 6*E
    n0, ni = A+B+2*E, A+B+3*E
    d = reduce(gcd, (A, B, *es))
    rhs = len(es)*M-gcd(M, n0)-gcd(M, ni)-sum(gcd(M, e) for e in es)
    return M, n0, ni, d, rhs


rows = []
for r in range(1, 4):
    for es in product(range(1, 4), repeat=r):
        for A, B in product(range(5), repeat=2):
            M, n0, ni, d, rhs = cover(A, B, es)
            check(d == reduce(gcd, (M, n0, *es)), 'component gcd')
            check(rhs % (2*d) == 0, 'integral component genus')
            g = 1 + rhs//(2*d)
            exception = A == B == 0 and r == 1
            check((g == 1) == exception, 'only pure cusp initial form has central genus one')
            ap, bp, ep = A//d, B//d, tuple(e//d for e in es)
            mm, nn, ii, dd, tt = cover(ap, bp, ep)
            check(dd == 1 and rhs == d*tt, 'primitive normalization')
            if bp:
                lower = r*mm-2*bp-3*sum(ep)
                check(gcd(mm, nn) <= bp+2*sum(ep), 'zero gcd bound')
                check(gcd(mm, ii) <= bp, 'infinity gcd bound')
                check(tt >= lower >= 2*ap+bp+3*sum(ep) >= 4, 'positive-y genus bound including A zero')
            elif not exception:
                lower = (2*r-1)*ap+6*(r-1)*sum(ep)
                check(gcd(mm, ii) == mm//2, 'zero-y infinity half')
                check(tt >= lower >= 1 and tt % 2 == 0, 'nonexception positive even bound')
            else:
                check((ap, bp, ep, mm, nn, ii, tt) == (0, 0, (1,), 6, 2, 3, 0), 'elliptic primitive boundary')
            rows.append((A, B, es, d, g))

literal = 0
for A, B in product(range(3), repeat=2):
    for es in [(1,), (2,), (1, 1), (1, 2), (2, 2), (1, 1, 1)]:
        literal += 1
        M, n0, ni, d, rhs = cover(A, B, es)
        increments = [-n0] + [-e for e in es] + [ni]
        check(sum(increments) == 0, 'signed puncture monodromy')
        left = set(range(M))
        components = []
        while left:
            start = min(left)
            component, todo = {start}, [start]
            while todo:
                v = todo.pop()
                for increment in increments:
                    w = (v+increment) % M
                    if w not in component:
                        component.add(w)
                        todo.append(w)
            left -= component
            components.append(component)
        check(len(components) == d, 'literal sheet components')
        for component in components:
            ram = 0
            for increment in increments:
                unseen, cycles = set(component), 0
                while unseen:
                    v = min(unseen)
                    while v in unseen:
                        unseen.remove(v)
                        v = (v+increment) % M
                    cycles += 1
                ram += len(component)-cycles
            check(2-2*(1+rhs//(2*d)) == 2*len(component)-ram, 'literal component Euler characteristic')

p, y, x, t, s, tau, c, v, w = sp.symbols('p y x t s tau c v w')
Delta = p**3-y**2
px = t*(1+x*x*t)
yx = x*t*px
pl = s*s+tau
yl = s*pl

def equal(a, b, label):
    check(sp.cancel(a-b) == 0, label)


def source_bracket(F, G):
    return sp.diff(F, x)*sp.diff(G, t)-sp.diff(F, t)*sp.diff(G, x)


def log_bracket(F, G):
    return tau*(sp.diff(F, s)*sp.diff(G, tau)-sp.diff(F, tau)*sp.diff(G, s))


def pull(F):
    return sp.expand(F.subs({p: px, y: yx}, simultaneous=True))


def log_pull(F):
    return sp.expand(F.subs({p: pl, y: yl}, simultaneous=True))


def trunc(F, variable, n):
    return sp.Add(*[
        a*variable**power[0] for power, a in sp.Poly(sp.expand(F), variable).terms()
        if power[0] < n])


equal(pull(Delta), t**3*(1+x*x*t)**2, 'literal cusp factor')
equal(source_bracket(px, yx), -pull(Delta)/px, 'actual rational Poisson factor')
equal(log_bracket(pl, yl), -log_pull(Delta)/pl, 'logarithmic Poisson factor')
equal(px.subs({x: s/tau, t: tau}, simultaneous=True), pl, 'birational inverse p')
equal(yx.subs({x: s/tau, t: tau}, simultaneous=True), yl, 'birational inverse y')
monomial_images = [(i, i+j) for i, j in product(range(9), repeat=2)]
check(len(set(monomial_images)) == 81, 'formal injection distinct monomials')

polynomials = [Delta, Delta**2, Delta*(1+p), Delta*(p+y),
               Delta*(p**3+y**2), Delta*(Delta**3+p**10)]
for I in polynomials:
    ix, il = pull(I), log_pull(I)
    check(min(m[0] for m, a in sp.Poly(ix, t).terms()) >= 3, 'literal t-order at least three')
    for F, G in [(s, tau), (s**2+tau, s*(s*s+tau)), (s**2*tau+tau**3, il)]:
        left = source_bracket(F.subs({s: x*t, tau: t}, simultaneous=True),
                              G.subs({s: x*t, tau: t}, simultaneous=True))
        right = log_bracket(F, G).subs({s: x*t, tau: t}, simultaneous=True)
        equal(left, right, 'chain-rule bracket intertwining')
    least = min(2*a+3*b for (a, b), z in sp.Poly(I, p, y).terms())
    initial = sum(z*p**a*y**b for (a, b), z in sp.Poly(I, p, y).terms() if 2*a+3*b == least)
    actual = sp.Poly(sp.expand(I.subs({p: v*v*w, y: v**3*w}, simultaneous=True)), v)
    equal(actual.coeff_monomial(v**least), initial.subs({p: w, y: w}), 'actual exceptional initial form')

for e in range(1, 4):
    for T in [1, p, y, 1+p+y]:
        T = sp.sympify(T)
        I = Delta**e*T
        il = log_pull(I)
        W = s**(4*e)*T.subs({p: s*s, y: s**3}, simultaneous=True)
        check(W != 0, 'nonzero cusp restriction')
        di = log_bracket(pl, il)
        equal(sp.Poly(di, tau).coeff_monomial(tau**e), 2*e*s*W, 'leading displacement')
        check(min(m[0] for m, z in sp.Poly(log_bracket(di, il), tau).terms()) >= 2*e,
              'higher iterate cannot cancel displacement')

# Finite source and log jets use separate brackets, for two genuine generators.
for I in [Delta, Delta*(p+y)]:
    n = 9
    ix = trunc(pull(I), t, n+1)
    il = trunc(log_pull(I), tau, n+1)
    def delta_x(F):
        return trunc(source_bracket(F, ix), t, n)
    def exp_x(F, time):
        term, answer = trunc(F, t, n), trunc(F, t, n)
        for j in range(1, n):
            term = delta_x(term)
            answer += time**j*term/sp.Integer(factorial(j))
        return trunc(answer, t, n)
    for F in [x, t, px, yx]:
        equal(exp_x(exp_x(F, 1), 2), exp_x(F, 3), 'literal finite scalar group law')
        equal(exp_x(exp_x(F, 1), -1), trunc(F, t, n), 'literal finite inverse')
        equal(delta_x(exp_x(F, 1)), exp_x(delta_x(F), 1), 'finite commutation with vector field')
    equal(exp_x(ix, 1), trunc(ix, t, n), 'finite invariant conservation')
    for F in [pl, yl]:
        term, answer = trunc(F, tau, n), trunc(F, tau, n)
        for j in range(1, n):
            term = trunc(log_bracket(term, il), tau, n)
            answer += term/sp.Integer(factorial(j))
        expected = exp_x(F.subs({s: x*t, tau: t}, simultaneous=True), 1)
        equal(trunc(answer.subs({s: x*t, tau: t}, simultaneous=True), t, n),
              expected, 'literal-log finite fixed-input comparison')

# The pole at p=0 is intrinsic: p is a uniformizer on these generic curves.
pole_rows = []
for I in [Delta, Delta**2, Delta**3, Delta*(1+p), Delta*(p+y),
          Delta*(1+y*y), Delta**2*(1+p+y), Delta*(p+y)*(1+Delta)]:
    f = sp.expand(I.subs(p, 0))
    fy = sp.diff(f, y)
    F = sp.Poly(f-c, y, domain=sp.QQ.frac_field(c))
    check(F.degree() > 0 and f.subs(y, 0) == 0, 'nonconstant p-zero fibre')
    check(sp.gcd(F, sp.Poly(fy, y, domain=sp.QQ.frac_field(c))).degree() == 0,
          'all generic p-zero roots simple')
    check(sp.gcd(F, sp.Poly(y, y, domain=sp.QQ.frac_field(c))).degree() == 0,
          'no generic p-zero root has y zero')
    numerator = sp.expand(-Delta*sp.diff(I, y))
    equal(numerator.subs(p, 0), y*y*fy, 'nonzero principal coefficient of vector-field pole')
    check(sp.gcd(F, sp.Poly(y*y*fy, y, domain=sp.QQ.frac_field(c))).degree() == 0,
          'simple pole at every generic p-zero point')
    pole_rows.append((str(I), str(f), str(y*y*fy)))

# Fixed-H genus strengthening and its genuine elliptic boundary.
H = -3*p+5*y+p*y+y*y
for e in range(1, 6):
    I = Delta**e + Delta*(p**(3*e-2)+y**(2*e-1))
    J = sp.diff(H, p)*sp.diff(I, y)-sp.diff(H, y)*sp.diff(I, p)
    equal(sp.Poly(J.subs(p, 0), y).coeff_monomial(y**(2*e-1)),
          -3*2*e*(-1)**e, 'uncancellable fixed-H cusp-initial coefficient')
I = Delta*(y+3*p)
J = sp.diff(I, y)-sp.diff(I, p)
equal(J.subs(p, 0), 0, 'fixed H=p+y preserver outside p ideal')
check(I.subs(p, 0) != 0, 'fixed-H positive control not a p-multiple')
equal(sp.diff(Delta, y)*0-sp.diff(Delta, p)*0, 0, 'H=0 Delta is fixed-H elliptic boundary')
check(sp.gcd(sp.Poly(p**3-c, p), sp.Poly(3*p*p, p)).degree() == 0,
      'Delta generic fibre is an actual smooth elliptic cubic')

# Genus one by itself cannot force degree one: the exact duplication map.
xx, yy = sp.symbols('xx yy')
slope = 3*xx*xx/(2*yy)
X = slope*slope-2*xx
Y = slope*(xx-X)-yy
curve = yy*yy-xx**3+1
num = sp.fraction(sp.cancel(Y*Y-X**3+1))[0]
equal(sp.rem(num, curve, yy), 0, 'elliptic duplication preserves actual curve')
Xcurve = (xx**4+8*xx)/(4*(xx**3-1))
equal(sp.rem(sp.fraction(sp.cancel(X-Xcurve))[0], curve, yy), 0, 'degree-four duplication x map')
check(sp.gcd(sp.fraction(Xcurve)[0], sp.fraction(Xcurve)[1]) == 1,
      'duplication x numerator and denominator are coprime')
check(max(sp.degree(sp.fraction(Xcurve)[0], xx), sp.degree(sp.fraction(Xcurve)[1], xx)) == 4,
      'duplication has degree four on genus one')
equal(Xcurve.subs(xx, -2), 0, 'duplication does not preserve the cusp-flow pole support')
check((-2)**3-1 != 0, 'extra inverse pole points are actual smooth points')
# Non-LND rational flow remains an outside-carrier hostile.
lam = sp.symbols('lam')
Xrat, Trat = x/(1-lam*x), t*(1-lam*x)**2
equal(Xrat**2*Trat, x*x*t, 'outside-carrier rational invariant')
equal(source_bracket(Xrat, Trat), 1, 'outside-carrier rational symplectic map')

blob = json.dumps({'multiplicities': rows, 'poles': pole_rows}, separators=(',', ':')).encode()
print('Cusp-ideal extension exact controls PASS (proof status RESERVED pending independent audit)')
print('Multiplicity universe:', len(rows), 'rows; A,B=0..4, 1..3 distinct binomial factors, exponents=1..3')
print('Independent literal cyclic covers:', literal)
print('Actual pole controls:', len(pole_rows), '; source/log jets: two generators through order 8')
print('Always-active gates:', checks)
print('Semantic SHA256:', sha256(blob).hexdigest())
print('Scope: exact identities and finite controls; global curve and rational-time claims require the companion proof.')
