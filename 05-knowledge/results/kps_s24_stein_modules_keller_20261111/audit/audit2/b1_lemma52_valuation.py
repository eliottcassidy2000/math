"""b1_lemma52_valuation.py - second audit, Lemma 5.2 (README 5.2).

Claim: every extension w to K = Frac(B) of nu = -deg_(u,v) satisfies w(d_u phi) >= w(phi) + 1, w(d_v phi) >= w(phi) + 1.

Tests
 (i)   double cover (u, v) = (x^2, y/x) on C^* x C: the extension is the monomial valuation w(x) = -1/2,
       w(y) = -3/2 (x = pi^-1, y = s pi^-3, pi^2 = tau; e = 2, f = 1, c = 1).
 (i')  double cover (u, v) = (u, eta^2) on C x C^*: pi = 1/eta, pi^2 = tau / s, so c = 1/s is NOT constant
       (tests the term (j/e)(c'/c) of the derivation D~_2 in the proof).
 (ii)  Keller automorphism (x + y^2, y + (x + y^2)^3): K = k, w = nu.
 (iii) the slice of README 6.2: see b1b_lemma52_slice_places.py.
"""
import itertools
import random
from fractions import Fraction
import sympy as sp

random.seed(52)

# ------------------------------------------------------------------ (i) and (i'): monomial valuations
def test_double_cover():
    x, y = sp.symbols("x y")
    wx, wy = Fraction(-1, 2), Fraction(-3, 2)

    def wval(expr):                       # monomial valuation of a Laurent polynomial
        expr = sp.expand(expr)
        if expr == 0:
            return None
        vals = []
        for term in sp.Add.make_args(expr):
            k = sp.degree(sp.together(term * x ** 40).as_numer_denom()[0], x) - 40
            l = sp.degree(term * 1, y) if term.has(y) else 0
            vals.append(k * wx + l * wy)
        return min(vals)

    du = lambda h: sp.expand(sp.diff(h, x) / (2 * x) + y * sp.diff(h, y) / (2 * x ** 2))
    dv = lambda h: sp.expand(x * sp.diff(h, y))
    u, v = x ** 2, y / x
    assert du(u) == 1 and dv(v) == 1 and du(v) == 0 and dv(u) == 0
    # w restricts to -deg on C[u, v]
    for i, j in itertools.product(range(4), range(4)):
        assert wval(u ** i * v ** j) == -(i + j)
    bad = tot = 0
    eq_u = eq_v = 0
    tests = [x ** k * y ** l for k in range(-5, 6) for l in range(0, 5)]
    for _ in range(60):
        tests.append(sum(random.randint(-4, 4) * x ** random.randint(-5, 5) * y ** random.randint(0, 4) for _ in range(4)))
    for h in tests:
        if sp.expand(h) == 0:
            continue
        w0 = wval(h)
        for d, name in ((du, "u"), (dv, "v")):
            w1 = wval(d(h))
            tot += 1
            if w1 is not None and w1 < w0 + 1:
                bad += 1
            if w1 is not None and w1 == w0 + 1:
                if name == "u":
                    eq_u += 1
                else:
                    eq_v += 1
    print("(i)  double cover (x^2, y/x), w(x) = -1/2, w(y) = -3/2: %d inequalities tested, %d violations; equality in %d (d_u), %d (d_v) cases" % (tot, bad, eq_u, eq_v))

    # (i') (u, v) = (U, eta^2) on C x C^*: w(U) = -1, w(eta) = -1/2
    U, e = sp.symbols("U eta")
    wU, we = Fraction(-1), Fraction(-1, 2)

    def wval2(expr):
        expr = sp.expand(expr)
        if expr == 0:
            return None
        vals = []
        for term in sp.Add.make_args(expr):
            k = sp.degree(term, U) if term.has(U) else 0
            l = sp.degree(sp.together(term * e ** 40).as_numer_denom()[0], e) - 40
            vals.append(k * wU + l * we)
        return min(vals)

    du2 = lambda h: sp.expand(sp.diff(h, U))
    dv2 = lambda h: sp.expand(sp.diff(h, e) / (2 * e))
    assert dv2(e ** 2) == 1 and du2(e ** 2) == 0
    bad = tot = 0
    tests = [U ** k * e ** l for k in range(0, 5) for l in range(-5, 6)]
    for _ in range(60):
        tests.append(sum(random.randint(-4, 4) * U ** random.randint(0, 4) * e ** random.randint(-5, 5) for _ in range(4)))
    for h in tests:
        if sp.expand(h) == 0:
            continue
        w0 = wval2(h)
        for d in (du2, dv2):
            w1 = wval2(d(h))
            tot += 1
            if w1 is not None and w1 < w0 + 1:
                bad += 1
    print("(i') double cover (U, eta^2), w(U) = -1, w(eta) = -1/2 (pi^2 = tau/s): %d inequalities tested, %d violations" % (tot, bad))
    # the derivation D~_2 of the proof on pi = 1/eta: D~_2 pi = (1/e)(c'/c) pi with c = 1/s, e = 2  =>  d_v eta = tau D~_2 (pi^-1)
    s_, tau_ = sp.symbols("s tau", positive=True)
    pi_ = sp.sqrt(tau_ / s_)
    D2pi = sp.Rational(1, 2) * (sp.diff(1 / s_, s_) / (1 / s_)) * pi_
    lhs = tau_ * (-pi_ ** -2 * D2pi)                # tau * D~_2(pi^-1)
    rhs = (1 / pi_) / (2 * (s_ / tau_))             # d_v eta = 1/(2 eta) = eta/(2 v), v = s/tau, eta = 1/pi
    print("     proof formula d_v = tau D~_2 on eta = 1/pi with c = 1/s: difference = %s" % sp.simplify(lhs - rhs))


# ------------------------------------------------------------------ (ii) Keller automorphism
def test_automorphism():
    x, y, U, V = sp.symbols("x y U V")
    f = x + y ** 2
    g = y + f ** 3
    # inverse: y = V - U^3, x = U - y^2
    yy = V - U ** 3
    xx = U - yy ** 2
    assert sp.expand(f.subs({x: xx, y: yy}, simultaneous=True) - U) == 0 and sp.expand(g.subs({x: xx, y: yy}, simultaneous=True) - V) == 0
    bad = tot = 0
    for _ in range(40):
        h = sum(random.randint(-3, 3) * x ** random.randint(0, 3) * y ** random.randint(0, 3) for _ in range(3))
        H = sp.expand(h.subs({x: xx, y: yy}, simultaneous=True))
        if H == 0:
            continue
        d0 = sp.Poly(H, U, V).total_degree()
        for var in (U, V):
            dH = sp.expand(sp.diff(H, var))
            tot += 1
            if dH != 0 and sp.Poly(dH, U, V).total_degree() > d0 - 1:
                bad += 1
    print("(ii) automorphism (x + y^2, y + (x + y^2)^3): w = -deg in (u, v); %d inequalities tested, %d violations" % (tot, bad))


if __name__ == "__main__":
    test_double_cover()
    test_automorphism()
