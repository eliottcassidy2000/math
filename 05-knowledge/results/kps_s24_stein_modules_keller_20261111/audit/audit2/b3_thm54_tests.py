"""b3_thm54_tests.py - second audit, Theorem 5.4 (h != 0 defined over k_0 is not in M_(a,b) for (a, b) algebraically
independent over k_0): independent tests with symbolic centres.

 (1) double cover (x^2, y/x): residue functionals of README Prop. 6.1 with symbolic (alpha, b), alpha = sqrt(a):
       lambda_+-(x^k y^l) = l! Res_(x = +-alpha) x^(k+l+1) e^(b(x^2 - a)) / (x^2 - a)^(l+1)
                          = (+-1)^(k+l) (1/2) sum_j C(l, j) ff((k+l)/2, j) b^(l-j) alpha^(k+l-2j)      (closed form)
     - closed form against sympy residues; both kill the generators of M_z (symbolic a, b);
     - no monomial has lambda_+ identically zero; the leading terms b^l alpha^(k+l) are distinct, so lambda_+ is
       injective on C[x, 1/x, y] (this proves Theorem 5.4 on the double cover independently of Lemmas 5.2-5.3);
     - zero loci of three examples, and a machine check of two of them with the engine.
 (2) standard pair: h = (e^(d_x d_y) h)(E_a, E_b) 1, so [h]_z = (e^(d_x d_y) h)(a, b) [1].
 (3) open immersion C^* x C -> C^2: L(h) = (e^(d_x d_y) h)(a, b) kills M_z (a != 0), L(1/x) = 1/a; at a = 0, S_z = 0.
 (4) Corollary 5.6, last sentence: an element h NOT defined over k_0 with [h]_z = 0 at a generic centre.
"""
import random
from math import comb
import sympy as sp
import eng2
from eng2 import RingL, add, scal, Echelon, inv

random.seed(54)
x, y, a, b, al = sp.symbols("x y a b alpha")


def ff(z, j):
    out = sp.Integer(1)
    for i in range(j):
        out *= (z - i)
    return out


def lam_closed(k, l, sign=1):
    return sp.Rational(1, 2) * sign ** (k + l) * sum(comb(l, j) * ff(sp.Rational(k + l, 2), j) * b ** (l - j) * al ** (k + l - 2 * j) for j in range(l + 1))


def lam_residue(k, l, sign=1):
    expr = x ** (k + l + 1) * sp.exp(b * (x ** 2 - al ** 2)) / (x ** 2 - al ** 2) ** (l + 1)
    # residue at x = sign * alpha: write x = sign*alpha + e
    e = sp.symbols("e")
    ex = expr.subs(x, sign * al + e)
    ser = sp.series(ex * e ** (l + 1), e, 0, l + 1).removeO()          # regular part times e^(l+1)
    return sp.simplify(sp.factorial(l) * sp.Poly(sp.expand(ser), e).coeff_monomial(e ** l))


def lam(h, sign=1):
    """lambda_sign of a Laurent polynomial in x, polynomial in y."""
    h = sp.expand(h)
    tot = 0
    for term in sp.Add.make_args(h):
        c, rest = term.as_independent(x, y)
        l = sp.degree(rest, y) if rest.has(y) else 0
        rest_x = sp.simplify(rest / y ** l)
        k = sp.degree(sp.together(rest_x * x ** 60).as_numer_denom()[0], x) - 60 if rest_x != 1 else 0
        assert sp.simplify(rest - x ** k * y ** l) == 0
        tot += c * lam_closed(k, l, sign)
    return sp.expand(tot)


def part1():
    bad = 0
    for k in range(-4, 5):
        for l in range(0, 4):
            for sgn in (1, -1):
                if sp.simplify(lam_closed(k, l, sgn) - lam_residue(k, l, sgn)) != 0:
                    bad += 1
    print("(1) double cover: closed form of lambda_+- against the residue, |k| <= 4, l <= 3: %d mismatches of 72" % bad)
    W = (x ** 2 - al ** 2) * (y / x - b)
    Wx, Wy = sp.diff(W, x), sp.diff(W, y)
    bad = tot = 0
    for k in range(-6, 7):
        for l in range(0, 6):
            m = x ** k * y ** l
            g1 = sp.expand(sp.diff(m, x) - Wx * m)
            g2 = sp.expand(sp.diff(m, y) - Wy * m)
            for g in (g1, g2):
                for sgn in (1, -1):
                    tot += 1
                    if sp.simplify(lam(g, sgn)) != 0:
                        bad += 1
    print("    lambda_+- kill the generators (d/dx - W_x) m, (d/dy - W_y) m, |k| <= 6, l <= 5, symbolic (alpha, b): %d failures of %d" % (bad, tot))
    # leading terms
    lead = {}
    zero = 0
    for k in range(-8, 9):
        for l in range(0, 7):
            f = sp.Poly(sp.expand(lam_closed(k, l) * al ** 40), al, b)
            if f.is_zero:
                zero += 1
                continue
            top = max(f.monoms(), key=lambda mm: (mm[1], mm[0]))        # largest power of b, then of alpha
            lead[(k, l)] = (top[0] - 40, top[1])
    assert all(lead[(k, l)] == (k + l, l) for (k, l) in lead)
    print("    lambda_+(x^k y^l), |k| <= 8, l <= 6: identically zero for %d monomials; leading term (largest power of b) = b^l alpha^(k+l)/2 in all %d cases, pairwise distinct: %s"
          % (zero, len(lead), len(set(lead.values())) == len(lead)))
    print("    => lambda_+ is injective on C[x, 1/x, y] (triangular in the order (l, k)); so (lambda_+(h), lambda_-(h)) is never identically zero for h != 0.")
    for name, h in (("y", y), ("x y - 1", x * y - 1), ("x^2 - y^2", x ** 2 - y ** 2), ("x^-2 - 1/a", x ** -2 - 1 / al ** 2)):
        lp, lm = sp.factor(lam(h, 1)), sp.factor(lam(h, -1))
        print("    h = %-10s lambda_+ = %s ;  lambda_- = %s" % (name, lp, lm))
    print("    zero loci (a != 0): [y]_z = 0 iff 2ab + 1 = 0;  [xy - 1]_z = 0 iff ab = 0, i.e. b = 0;  [x^2 - y^2]_z = 0 iff a - a b^2 - 2b = 0.")
    # machine check with the engine: membership of y at (a, b) = (1, -1/2) (on the curve) and (1, 1) (off it); x y - 1 at (3, 0) and (3, 1)
    R = RingL()
    P = eng2.PRIME

    def member(h, av, bv, Np=22):
        u, v = R.mono(2, 0), R.mono(-1, 1)
        U = add(u, scal(R.one(), av), -1)
        V = add(v, scal(R.one(), bv), -1)
        Wd = R.mul(U, V)
        Wxd, Wyd = R.dx(Wd), R.dy(Wd)
        E = Echelon(lambda m: (abs(m[0]) + m[1], m))
        for k in range(-Np, Np + 1):
            for l in range(0, Np + 1 - abs(k)):
                m = R.mono(k, l)
                E.insert(add(R.dx(m), R.mul(Wxd, m), -1))
                E.insert(add(R.dy(m), R.mul(Wyd, m), -1))
        return E.contains(h)

    half = inv(2)
    yy = R.mono(0, 1)
    xy1 = add(R.mono(1, 1), R.one(), -1)
    x2y2 = add(R.mono(2, 0), R.mono(0, 2), -1)
    res = (member(yy, 1, (-half) % P), member(yy, 1, 1), member(xy1, 3, 0), member(xy1, 3, 1),
           member(x2y2, 4, half), member(x2y2, 4, 2))
    print("    engine (generators of weight <= 22): y in M at (1, -1/2): %s, at (1, 1): %s;  xy - 1 in M at (3, 0): %s, at (3, 1): %s;  x^2 - y^2 in M at (4, 1/2) [4 - 1 - 1 = 2 != 0]: %s, at (4, 2): %s"
          % res)
    # the curve a - a b^2 - 2 b = 0 at b = 1/2: a (3/4) = 1, a = 4/3
    a43 = 4 * inv(3) % P
    print("    engine: x^2 - y^2 in M at (4/3, 1/2) (on the curve a - a b^2 - 2b = 0): %s" % member(x2y2, a43, half))


def part2():
    X, Y = sp.symbols("X Y")
    bad = 0
    for _ in range(25):
        h = sum(random.randint(-5, 5) * X ** random.randint(0, 4) * Y ** random.randint(0, 4) for _ in range(4))
        H = h
        g = h
        j = 0
        while True:
            j += 1
            g = sp.diff(g, X, Y)
            if g == 0:
                break
            H += g / sp.factorial(j)
        H = sp.Poly(sp.expand(H), X, Y)
        # apply H(E_a, E_b) to 1, E_a = X - d_Y, E_b = Y - d_X

        def Ea(p):
            return sp.expand(X * p - sp.diff(p, Y))

        def Eb(p):
            return sp.expand(Y * p - sp.diff(p, X))
        tot = 0
        for (i, jj), c in H.terms():
            p = sp.Integer(1)
            for _i in range(jj):
                p = Eb(p)
            for _i in range(i):
                p = Ea(p)
            tot += c * p
        if sp.expand(tot - h) != 0:
            bad += 1
    print("(2) standard pair: h = (e^(d_x d_y) h)(E_a, E_b) 1 for 25 random polynomials: %d failures  (so [h]_z = (e^(d_x d_y) h)(a, b) [1]_z)" % bad)


def part3():
    def L(h):
        # (e^(d_x d_y) h)(a, b) for h in C[x, 1/x, y]
        tot, g, j = h, h, 0
        while True:
            j += 1
            g = sp.expand(sp.diff(g, x, y))
            if g == 0:
                break
            tot += g / sp.factorial(j)
        return sp.simplify(tot.subs({x: a, y: b}))
    bad = tot = 0
    for k in range(-5, 6):
        for l in range(0, 5):
            m = x ** k * y ** l
            g1 = sp.expand(sp.diff(m, x) - (y - b) * m)      # T_u m = d_u m - V m
            g2 = sp.expand(sp.diff(m, y) - (x - a) * m)      # T_v m = d_v m - U m
            for g in (g1, g2):
                tot += 1
                if L(g) != 0:
                    bad += 1
    print("(3) open immersion C^* x C: L(h) = (e^(d_x d_y) h)(a, b) kills T_u m, T_v m for |k| <= 5, l <= 4: %d failures of %d;  L(1/x) = %s, L(y/x^2) = %s"
          % (bad, tot, L(1 / x), sp.factor(L(y / x ** 2))))
    print("    at a = 0: (d/dy - x)(x^k y^j) = j x^k y^(j-1) - x^(k+1) y^j for all k in Z gives x^k = 0 and then everything = 0 modulo M_z: S_z = 0 exactly on the line a = 0.")


def part4():
    # double cover, generic centre: h = x^-2 - 1/a is not defined over Q(coefficients) and lies in M_z
    h = x ** -2 - 1 / al ** 2
    print("(4) double cover: h = x^-2 - 1/a: lambda_+ = %s, lambda_- = %s  => [h]_z = 0 at EVERY centre with a != 0, although h != 0 and d_u h != 0."
          % (sp.simplify(lam(h, 1)), sp.simplify(lam(h, -1))))
    print("    So 'if d_u^m h != 0 for all m then [d_u^m h]_z != 0 for all m' (Corollary 5.6) needs h defined over k_0 (m = 0 fails for this h).")


if __name__ == "__main__":
    part1()
    part2()
    part3()
    part4()
