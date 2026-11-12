"""a4_radial_identity.py - audit of Proposition A4 (radial identity) with explicit Stein functionals (sympy, exact).

  R1  standard pair (u, v) = (x, y), centre (a, b) symbolic: the Stein functional is L_z(h) = (exp(d_x d_y) h)(a, b)
      (S24, after Proposition 2.4).  Check L_z((UV)^m h) = L_z((E+1)...(E+m) h), E = (U d_u + V d_v)/2, and
      L_z((U d_u - V d_v) h) = 0, on monomials h = x^i y^j, i, j <= 3, m <= 3.
  R2  double cover (u, v) = (x^2, y/x) on C^* x C (S24, Proposition 6.1), centre (4, 3): two functionals
      lambda_+-(x^k y^l) = l! Res_(x = +-2) x^(k+l+1) exp(b (x^2 - a)) / (x^2 - a)^(l+1).
      Check the same two identities on h = x^k y^l, -3 <= k <= 3, l <= 2, m <= 2.
  R3  [(UV - c)^m] = m! e_m(-c) [1]: the truncated exponentials e_m, e_(m+1) have no common root (resultants), and
      on the double cover at a centre (0, b) the Stein space is zero (S24, Proposition 6.1), so there every UV - c is
      a null element: the remark 'UV - c is never a null element' needs [1]_z != 0.
"""
import sympy as sp

x, y, a, b, c = sp.symbols("x y a b c")


def heat(h, order=12):
    """exp(d_x d_y) h for a polynomial h"""
    res = 0
    term = h
    n = 0
    while term != 0 and n <= order:
        res += term / sp.factorial(n)
        term = sp.diff(term, x, y)
        n += 1
    return sp.expand(res)


def r1():
    U, V = x - a, y - b
    L = lambda h: sp.expand(heat(sp.expand(h)).subs({x: a, y: b}))
    E = lambda h: (U * sp.diff(h, x) + V * sp.diff(h, y)) / 2
    tot = bad = 0
    for i in range(4):
        for j in range(4):
            h = x ** i * y ** j
            for m in range(1, 4):
                rhs = h
                for t in range(m, 0, -1):
                    rhs = sp.expand(E(rhs) + t * rhs)
                lhs = L((U * V) ** m * h)
                tot += 1
                bad += (sp.expand(lhs - L(rhs)) != 0)
            tot += 1
            bad += (L(U * sp.diff(h, x) - V * sp.diff(h, y)) != 0)
    print("R1  standard pair, symbolic centre: %d identities checked, %d failures" % (tot, bad))
    # sanity of the functional itself: it kills (d_x - V) and (d_y - U)
    bad2 = 0
    for i in range(4):
        for j in range(4):
            h = x ** i * y ** j
            bad2 += (L(sp.diff(h, x) - V * h) != 0) + (L(sp.diff(h, y) - U * h) != 0)
    print("    the functional kills (d_u - V) h and (d_v - U) h on 16 monomials: %d failures" % bad2)


def r2():
    A, Bc = 4, 3

    def lam(sign, k, l):
        expr = sp.factorial(l) * x ** (k + l + 1) * sp.exp(Bc * (x ** 2 - A)) / (x ** 2 - A) ** (l + 1)
        return sp.nsimplify(sp.residue(expr, x, sign * 2))

    cache = {}

    def lam_poly(sign, h):
        h = sp.expand(h)
        tot = 0
        for term in sp.Add.make_args(h):
            coeff, mono = term.as_coeff_Mul()
            pw = mono.as_powers_dict()
            k = pw.get(x, 0)
            l = pw.get(y, 0)
            key = (sign, int(k), int(l))
            if key not in cache:
                cache[key] = lam(sign, int(k), int(l))
            tot += coeff * cache[key]
        return sp.nsimplify(tot)

    u, v = x ** 2, y / x
    U, V = u - A, v - Bc
    du = lambda h: sp.diff(h, x) / (2 * x) + y * sp.diff(h, y) / (2 * x ** 2)
    dv = lambda h: x * sp.diff(h, y)
    assert sp.simplify(du(u) - 1) == 0 and sp.simplify(dv(v) - 1) == 0 and sp.simplify(du(v)) == 0 and sp.simplify(dv(u)) == 0
    E = lambda h: (U * du(h) + V * dv(h)) / 2
    # the functionals kill the Stein kernel
    bad0 = tot0 = 0
    for k in range(-2, 3):
        for l in range(0, 3):
            h = x ** k * y ** l
            for sign in (1, -1):
                tot0 += 2
                bad0 += (lam_poly(sign, du(h) - V * h) != 0) + (lam_poly(sign, dv(h) - U * h) != 0)
    print("R2  double cover, centre (4, 3): lambda_+- kill (d_u - V) h and (d_v - U) h: %d checked, %d failures" % (tot0, bad0))
    tot = bad = 0
    for k in range(-3, 4):
        for l in range(0, 3):
            h = x ** k * y ** l
            for m in (1, 2):
                rhs = h
                for t in range(m, 0, -1):
                    rhs = sp.expand(E(rhs) + t * rhs)
                lhs = sp.expand((U * V) ** m * h)
                for sign in (1, -1):
                    tot += 1
                    bad += (sp.simplify(lam_poly(sign, lhs) - lam_poly(sign, rhs)) != 0)
            for sign in (1, -1):
                tot += 1
                bad += (sp.simplify(lam_poly(sign, U * du(h) - V * dv(h))) != 0)
    print("    radial identity [(UV)^m h] = [(E+1)...(E+m) h] (m = 1, 2) and rotation identity: %d checked, %d failures" % (tot, bad))
    print("    lambda_+-(1) =", lam_poly(1, sp.Integer(1)), lam_poly(-1, sp.Integer(1)),
          "; lambda_+-(UV) =", lam_poly(1, U * V), lam_poly(-1, U * V),
          "; lambda_+-((UV)^2) =", lam_poly(1, (U * V) ** 2), lam_poly(-1, (U * V) ** 2))


def r3():
    t = sp.symbols("t")
    e = lambda m: sum(t ** j / sp.factorial(j) for j in range(m + 1))
    res = [sp.resultant(e(m), e(m + 1), t) for m in range(1, 9)]
    print("R3  resultants of consecutive truncated exponentials e_m, e_(m+1), m = 1..8:", res, "-> all nonzero:", all(r != 0 for r in res))
    print("    (e_(m+1) - e_m = t^(m+1)/(m+1)! and e_m(0) = 1: a common root is impossible.)")
    # expansion [(UV - c)^m] = m! e_m(-c) [1] from [(UV)^k] = k! [1]
    ok = True
    for m in range(1, 8):
        lhs = sum(sp.binomial(m, k) * (-c) ** (m - k) * sp.factorial(k) for k in range(m + 1))
        rhs = sp.factorial(m) * sum((-c) ** j / sp.factorial(j) for j in range(m + 1))
        ok = ok and (sp.expand(lhs - rhs) == 0)
    print("    sum_k C(m,k) (-c)^(m-k) k! = m! e_m(-c), m = 1..7:", ok)


if __name__ == "__main__":
    r1()
    r2()
    r3()
