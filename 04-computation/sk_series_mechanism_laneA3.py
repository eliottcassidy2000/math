#!/usr/bin/env python3
"""
Lane A3 -- the k=3 mechanism, plus context / provenance cross-checks, for
    S(k) = sum_{n>=0} C(2n,n) C(4n,2n) / ((kn+1) 64^n).

HEADLINE (cmd `xdeg`):  S(k) = (2k/pi) * area-period of  w^2 = 1 + Im(zeta^k)
over the sector |zeta|<1, |arg zeta| < pi/(2k).   Integrate in x FIRST.
    deg_x ( 1 + Im (x+iy)^k ) = k-1,
so int dx / sqrt(1 + Im zeta^k) is
    k=1 : sqrt of a constant        -> algebraic
    k=2 : sqrt of a linear          -> algebraic
    k=3 : sqrt of a QUADRATIC       -> arcsinh / arcsin      <-- last elementary case
    k=4 : sqrt of a CUBIC           -> incomplete elliptic F  <-- the wall
    k>=5: sqrt of degree >= 4       -> (hyper)elliptic
That is the whole mechanism, and the wall is at k=4 for a completely elementary
reason (Chebyshev / Liouville on ONE variable, not a Schwarz-list or CM argument).

Subcommands
  verify   numeric verification of the reduction chain
  xdeg     THE MECHANISM: x-first integration; k=1,2,3 done in closed form
  k3       the fully explicit k=3 evaluation, step by step, all numeric
  sep      separation-of-variables normal form + the Chebyshev "red herring" table
  thomae   S_5 Thomae orbit of 3F2(1/4,3/4,c;1,1+c;1) and a summability scan
  schwarz  Schwarz triples WITH the Schwarz reduction (corrects a lane-A1 claim)
  known    which nearby 3F2 families ARE classically evaluable
  all      everything
"""
import sys
import mpmath
from mpmath import mp, mpf, sqrt, log, atan, asin, asinh, pi, quad, hyp2f1, hyp3f2, sin, cos

DPS_ID = 80      # identity / closed-form checks
DPS_QUAD = 30    # 2-D quadrature (quadrature-limited anyway); stated everywhere


def F(z):
    return hyp2f1(mpf(1) / 4, mpf(3) / 4, 1, z)


def S_series(k):
    g = mpf(1) / k
    return hyp3f2(mpf(1) / 4, mpf(3) / 4, g, 1, 1 + g, 1)


def S3_closed():
    return (2 * sqrt(3) * log(sqrt(3) + sqrt(2)) - 2 * pi + 6 * atan(sqrt(2))) / pi


KNOWN = {1: lambda: 8 * sqrt(2) / (3 * pi),
         2: lambda: 4 * log(1 + sqrt(2)) / pi,
         3: S3_closed}


def rep(name, a, b=None, tol=None):
    if b is None:
        print(f"    {name:<58s} {mp.nstr(a, 28)}")
        return
    d = abs(a - b)
    t = tol if tol is not None else mpf(10) ** (-mp.dps + 12)
    print(f"    {name:<58s} resid {mp.nstr(d,4):>12s}   {'OK' if d < t else '** MISMATCH **'}")


# ---------------------------------------------------------------- verify ----
def cmd_verify():
    mp.dps = DPS_ID
    print("=" * 78)
    print(f"VERIFY -- reduction chain   (identity checks at dps={DPS_ID})")
    print("=" * 78)

    print("\n[V0] values / known closed forms")
    for k in (1, 2, 3, 4, 5):
        v = S_series(k)
        rep(f"S({k})", v)
        if k in KNOWN:
            rep("   vs closed form", v, KNOWN[k]())
    rep("arctan(sqrt2/5) == pi - 3 arctan(sqrt2)", atan(sqrt(2) / 5), pi - 3 * atan(sqrt(2)))

    print("\n[V1] the hypergeometric ODE is an EXACT derivative:  d/dz[z(1-z)F'] = (3/16)F")
    print("     (holds for 2F1(a,b;c;z) iff a+b = c = 1, i.e. the WHOLE 2F1(lam,1-lam;1;.) family)")
    for z in (mpf(1) / 7, mpf(3) / 5, mpf(9) / 10):
        lhs = mpmath.diff(lambda t: t * (1 - t) * mpmath.diff(F, t), z, h=mpf(10) ** -25)
        rep(f"  z = {mp.nstr(z,4)}", lhs, mpf(3) / 16 * F(z), tol=mpf(10) ** -25)

    print("\n[V2] hence k=1 in one line: S(1) = (16/3) lim_{z->1}(1-z)F'(z) = 16/(3 pi sqrt2)")
    rep("  16/(3 pi sqrt2)", mpf(16) / 3 / (pi * sqrt(2)), S_series(1))

    print("\n[V3] one IBP:  S(k) = 8sqrt2/(3 pi k) + (1-1/k) T(k),  T(k)=int_0^1 2F1(1/4,3/4;2;t^k)dt")
    mp.dps = 40
    for k in (1, 2, 3, 4, 5, 7):
        T = quad(lambda t: hyp2f1(mpf(1) / 4, mpf(3) / 4, 2, t ** k), [0, 1])
        rep(f"  k = {k}", 8 * sqrt(2) / (3 * pi * k) + (1 - mpf(1) / k) * T, S_series(k),
            tol=mpf(10) ** -30)

    print("\n[V4] A2-2 sector form:  S(k) = (2k/pi) int int_{|z|<1,|arg z|<pi/2k} dA / sqrt(1+Im z^k)")
    mp.dps = DPS_QUAD
    for k in (1, 2, 3, 4, 5):
        val = (2 * k / pi) * quad(
            lambda r: quad(lambda th: r / sqrt(1 + r ** k * sin(k * th)),
                           [-pi / (2 * k), 0, pi / (2 * k)]), [0, 1])
        rep(f"  k = {k} (dps={DPS_QUAD})", val, S_series(k), tol=mpf(10) ** -18)

    print("\n[V5] equivalent polar form S(k) = (2/pi) int_{-pi/2}^{pi/2} int_0^1 r dr dpsi/sqrt(1+r^k sin psi)")
    for k in (2, 3, 4):
        val = (2 / pi) * quad(lambda ps: quad(lambda r: r / sqrt(1 + r ** k * sin(ps)), [0, 1]),
                              [-pi / 2, 0, pi / 2])
        rep(f"  k = {k}", val, S_series(k), tol=mpf(10) ** -18)

    print("\n[V6] ODD-k isotrivial rescaling:  S(k) = (2k/pi) int_{-1}^1 m^{k-3} J_k(m) dm/sqrt(1-m^{2k}),")
    print("     J_k(m) = int_0^m u du / sqrt(1+u^k).   (k=3 is the unique k with weight m^0.)")
    for k in (1, 3, 5, 7):
        def Jk(m):
            return quad(lambda u: u / sqrt(1 + u ** k), [0, m])
        val = (2 * k / pi) * quad(lambda m: (m ** (k - 3) if k != 3 else 1) * Jk(m) /
                                  sqrt(1 - m ** (2 * k)), [-1, 0, 1])
        rep(f"  k = {k}", val, S_series(k), tol=mpf(10) ** -15)
    mp.dps = DPS_ID


# ------------------------------------------------------------------ xdeg ----
def cmd_xdeg():
    print("=" * 78)
    print("XDEG -- THE MECHANISM.  Integrate the sector period in x first.")
    print("=" * 78)
    import sympy as sp
    x, y = sp.symbols('x y', real=True)
    print("""
  S(k) = (2k/pi) int int_{Sigma_k} dx dy / sqrt(f_k),   f_k = 1 + Im (x+iy)^k,
         Sigma_k = { x^2+y^2 < 1,  |y| < x tan(pi/2k) }.
  Sigma_k is x-simple: for fixed y, x runs from x0(y)=|y|cot(pi/2k) to x1(y)=sqrt(1-y^2).
  So the inner integral is  int dx / sqrt(f_k(x,y))  with f_k a POLYNOMIAL IN x of degree
""")
    print(f"    {'k':>3s}  {'f_k = 1 + Im (x+iy)^k':<44s} {'deg_x':>5s}  inner int dx/sqrt(f_k)")
    for k in range(1, 8):
        fk = sp.expand(1 + sp.im(sp.expand((x + sp.I * y) ** k)))
        d = int(sp.degree(sp.Poly(fk, x), x))
        kind = {0: "algebraic (constant under the root)",
                1: "algebraic  (2 sqrt(linear))",
                2: "ELEMENTARY: arcsinh / arcsin / log",
                3: "incomplete ELLIPTIC F  <-- the wall",
                }.get(d, ("ELLIPTIC (genus 1)" if d == 4 else f"HYPERELLIPTIC (genus {(d - 1) // 2})"))
        print(f"    {k:3d}  {str(fk):<44s} {d:5d}  {kind}")
    print("""
  deg_x(1 + Im (x+iy)^k) = k-1 exactly (the x^{k-1} y term of the binomial), so
  Liouville/Chebyshev give:  the inner x-integral is elementary  <=>  k-1 <= 2  <=>  k <= 3.
  This is the whole mechanism, and it is one-variable elementary calculus.

  Consequences that match the other lanes exactly:
   * k = 4: inner = incomplete elliptic integral of the FIRST kind -> S(4) is a mean of
     incomplete F's, which is lane A1's independent finding (their A = (1/sqrt2) int F(u|1/2)du).
   * the "weight-1" alphabet {pi, log(alg), arctan(alg)} of S(1),S(2),S(3) is forced: the inner
     integral produces exactly arcsinh/arcsin, and one IBP in y leaves a purely ALGEBRAIC
     integral.
   * the algebraic numbers are read off the two boundary pieces:
       - the rays arg zeta = +-pi/(2k) give cot(pi/2k):  k=2 -> 1, k=3 -> sqrt3   <-- the sqrt3
       - the arc |zeta| = 1 gives 1 + sin(k th), whose square root is
         sqrt2 |cos(k th/2 - pi/4)|                                              <-- the sqrt2
""")
    print("  explicit inner antiderivatives, k <= 3:")
    for k in (1, 2, 3):
        fk = sp.expand(1 + sp.im(sp.expand((x + sp.I * y) ** k)))
        anti = sp.simplify(sp.integrate(1 / sp.sqrt(fk), x))
        print(f"    k={k}:  int dx/sqrt({fk}) = {anti}")
    print("\n  and the k=4 inner integral is genuinely elliptic:")
    f4 = sp.expand(1 + sp.im(sp.expand((x + sp.I * y) ** 4)))
    print(f"    f_4 = {f4}   (cubic in x, leading coeff 4y)")
    print("    w^2 = 4y x^3 - 4y^3 x + 1  ->  Weierstrass  X^3 + aX + b,  a = -y^2, b = 1/(4y)")
    a, b = -y ** 2, sp.Rational(1, 4) / y
    j = sp.simplify(1728 * 4 * a ** 3 / (4 * a ** 3 + 27 * b ** 2))
    print(f"    j(y) = {j}      (non-isotrivial; j -> 1728 nowhere, j=0 at y=0)")


# -------------------------------------------------------------------- k3 ----
def cmd_k3():
    print("=" * 78)
    print("K3 -- the fully explicit k = 3 evaluation")
    print("=" * 78)
    mp.dps = DPS_QUAD
    target = pi * S_series(3) / 6
    print(f"""
  pi S(3)/6 = int int_Sigma dx dy / sqrt(1 + 3x^2 y - y^3),
  Sigma = {{x^2+y^2<1, |y| < x/sqrt3}},  so y in [-1/2,1/2] and
  x from x0 = sqrt3 |y| to x1 = sqrt(1-y^2)  (they MEET at y = +-1/2: the sector corners
  e^{{+-i pi/6}}, which are exactly the zeros of 1 + Im zeta^3 on the boundary).

  target  pi S(3)/6 = {mp.nstr(target, 25)}
""")
    # step 1: raw double integral
    def f3(x, y):
        return 1 + 3 * x ** 2 * y - y ** 3
    v = quad(lambda y: quad(lambda x: 1 / sqrt(f3(x, y)),
                            [sqrt(3) * abs(y), sqrt(1 - y ** 2)]), [-mpf(1) / 2, 0, mpf(1) / 2])
    rep("  [1] raw double integral", v, target, tol=mpf(10) ** -18)

    # step 2: inner integral in closed form
    print("""
  [2] inner integral in closed form.  f_3 = (1-y^3) + 3y x^2, quadratic in x:
        y > 0:  int dx/sqrt(A+Bx^2) = B^{-1/2} asinh(x sqrt(B/A)),  A=1-y^3, B=3y
        y < 0:  int dx/sqrt(A'-B'x^2) = B'^{-1/2} arcsin(x sqrt(B'/A')), A'=1+z^3, B'=3z (y=-z)
      with (all four arguments simplify -- 1-y^3=(1-y)(1+y+y^2) is where sqrt3 will come from):
        x1^2 B/A = 3y(1+y)/(1+y+y^2)        x0^2 B/A = 9y^3/(1-y^3)
        x1^2 B'/A'= 3z(1-z)/(1-z+z^2)       x0^2 B'/A'= 9z^3/(1+z^3)""")

    def inner_pos(y):
        return (asinh(sqrt(3 * y * (1 + y) / (1 + y + y ** 2)))
                - asinh(sqrt(9 * y ** 3 / (1 - y ** 3)))) / sqrt(3 * y)

    def inner_neg(z):
        return (asin(sqrt(3 * z * (1 - z) / (1 - z + z ** 2)))
                - asin(sqrt(9 * z ** 3 / (1 + z ** 3)))) / sqrt(3 * z)

    v2 = quad(inner_pos, [0, mpf(1) / 2]) + quad(inner_neg, [0, mpf(1) / 2])
    rep("  [2] after the closed-form x-integration", v2, target, tol=mpf(10) ** -18)

    # step 3: IBP -> purely algebraic
    print("""
  [3] IBP in y against int (3y)^{-1/2} dy = (2/sqrt3) sqrt y.  The boundary terms at
      y = +-1/2 CANCEL (both arguments equal there: the corner), and vanish at y=0.  What is
      left is PURELY ALGEBRAIC, because
        1 + [3y(1+y)/(1+y+y^2)]  = (1+2y)^2/(1+y+y^2)
        1 + [9y^3/(1-y^3)]       = (1+2y)(1-2y+4y^2)/(1-y^3)
        1 - [3z(1-z)/(1-z+z^2)]  = (1-2z)^2/(1-z+z^2)
        1 - [9z^3/(1+z^3)]       = (1-2z)(1+2z+4z^2)/(1+z^3)
      (squares and cyclotomic factors -- this is the collapse that makes k=3 work).""")
    for expr, val in (("1+g1^2 = (1+2y)^2/(1+y+y^2)", None),):
        pass
    ys = [mpf(1) / 7, mpf(3) / 11, mpf(9) / 20]
    for yv in ys:
        g1sq = 3 * yv * (1 + yv) / (1 + yv + yv ** 2)
        rep(f"  identity check 1+g1^2, y={mp.nstr(yv,4)}", 1 + g1sq,
            (1 + 2 * yv) ** 2 / (1 + yv + yv ** 2), tol=mpf(10) ** -20)
        g2sq = 9 * yv ** 3 / (1 - yv ** 3)
        rep(f"  identity check 1+g2^2, y={mp.nstr(yv,4)}", 1 + g2sq,
            (1 + 2 * yv) * (1 - 2 * yv + 4 * yv ** 2) / (1 - yv ** 3), tol=mpf(10) ** -20)

    import sympy as sp
    t = sp.Symbol('t', positive=True)
    print("\n  [4] symbolic algebraic integrals (sympy), after y = t^2/... substitution")
    y = sp.Symbol('y', positive=True)
    g1 = sp.sqrt(3 * y * (1 + y) / (1 + y + y ** 2))
    g2 = sp.sqrt(9 * y ** 3 / (1 - y ** 3))
    z = sp.Symbol('z', positive=True)
    h1 = sp.sqrt(3 * z * (1 - z) / (1 - z + z ** 2))
    h2 = sp.sqrt(9 * z ** 3 / (1 + z ** 3))
    I1 = sp.simplify(sp.sqrt(y) * sp.diff(g1, y) / sp.sqrt(1 + g1 ** 2))
    I2 = sp.simplify(sp.sqrt(y) * sp.diff(g2, y) / sp.sqrt(1 + g2 ** 2))
    I3 = sp.simplify(sp.sqrt(z) * sp.diff(h1, z) / sp.sqrt(1 - h1 ** 2))
    I4 = sp.simplify(sp.sqrt(z) * sp.diff(h2, z) / sp.sqrt(1 - h2 ** 2))
    for nm, e in (("d(asinh g1)", I1), ("d(asinh g2)", I2), ("d(arcsin h1)", I3), ("d(arcsin h2)", I4)):
        print(f"    sqrt(y)*{nm:<14s} = {sp.simplify(e)}")
    print("""
    -> every one is an ALGEBRAIC function of y with only the radicals sqrt(y),
       sqrt(1+y+y^2), sqrt(1-y^3) = sqrt((1-y)(1+y+y^2)).  Setting y = t^2 rationalises
       sqrt y; the remaining quadratic 1+y+y^2 = 1+t^2+t^4 is where sqrt3 enters
       (disc = -3), and the endpoint y = 1/2 is where sqrt2 enters.""")

    print("""
  [5] the four algebraic integrands simplify completely (4y^3+8y^2+5y+1 = (y+1)(2y+1)^2,
      z^3-2z^2+2z-1 = (z-1)(z^2-z+1), 4z^2-4z+1 = (1-2z)^2):

      A1 = sqrt3 / (2 (1+y+y^2) sqrt(1+y))         A2 = 9y  / (2 (1-y^3) sqrt(1+8y^3))
      B1 = sqrt3 / (2 (1-z+z^2) sqrt(1-z))         B2 = 9z  / (2 (1+z^3) sqrt(1-8z^3))

      pi S(3)/6 = -(2/sqrt3) int_0^{1/2} [ (A1-A2) + (B1-B2) ] dy .""")
    mp.dps = DPS_QUAD

    def A1(y): return sqrt(3) / (2 * (1 + y + y ** 2) * sqrt(1 + y))

    def A2(y): return 9 * y / (2 * (1 - y ** 3) * sqrt(1 + 8 * y ** 3))

    def B1(z): return sqrt(3) / (2 * (1 - z + z ** 2) * sqrt(1 - z))

    def B2(z): return 9 * z / (2 * (1 + z ** 3) * sqrt(1 - 8 * z ** 3))

    half = mpf(1) / 2
    val = -(2 / sqrt(3)) * quad(lambda y: (A1(y) - A2(y)) + (B1(y) - B2(y)), [0, half])
    rep("  [5] purely algebraic form (after IBP)", val, target, tol=mpf(10) ** -18)

    print("""
  [6] where sqrt3 and sqrt2 come from, explicitly.  In int_0^{1/2} A1 dy put 1+y = tau^2:
        int_0^{1/2} A1 dy = sqrt3 int_1^{sqrt(3/2)} d tau / (tau^4 - tau^2 + 1),
        tau^4 - tau^2 + 1 = (tau^2 + sqrt3 tau + 1)(tau^2 - sqrt3 tau + 1)   <-- sqrt3
      and the upper endpoint is tau = sqrt(3/2) = sqrt3/sqrt2                <-- sqrt2 .
      Partial fractions over Q(sqrt3) give exactly the alphabet {log(alg), arctan(alg)}.""")
    lhs = quad(A1, [0, half])
    rhs = sqrt(3) * quad(lambda tau: 1 / (tau ** 4 - tau ** 2 + 1), [1, sqrt(mpf(3) / 2)])
    rep("    int A1 = sqrt3 int d tau/(tau^4-tau^2+1)", lhs, rhs, tol=mpf(10) ** -18)
    import sympy as sp
    tt = sp.Symbol('tau', positive=True)
    anti = sp.simplify(sp.integrate(1 / (tt ** 4 - tt ** 2 + 1), tt))
    print(f"    sympy antiderivative: {anti}")
    print("""
  So  pi S(3) = 6 * (integral of an ALGEBRAIC function over [0,1/2]), which Liouville/Risch
  evaluate in logs and arctans of algebraic numbers -- weight 1, exactly the observed
      pi S(3) = 2 sqrt3 log(sqrt3+sqrt2) - 2 pi + 6 arctan(sqrt2).
  At k = 4 step [2] already fails: the inner x-integral is an incomplete elliptic F, no IBP
  makes it algebraic, and the y-average of incomplete F's is the Eichler-type period that
  lanes A1/A2 identified.""")
    mp.dps = DPS_ID


# ------------------------------------------------------------------- sep ----
def cmd_sep():
    print("=" * 78)
    print("SEP -- separation-of-variables normal form, and why 'inner elementary' is a RED HERRING")
    print("=" * 78)
    mp.dps = DPS_QUAD
    print("""
  From the quadratic transformation (c = a+b+1/2 holds for EVERY k) and Euler,
      S(k) = (2^{1+4/k}/pi) int_0^{sqrt2} int_0^1 D^{-4/k} d rho ds / sqrt(2-s^2),
      D = (1+s) + (1-s) rho^{k/2}.
  The substitution rho = P ((1+s)/(1-s))^{2/k}  (i.e. D = (1+s)(1+P^{k/2})) SEPARATES it:

      S(k) = (2^{1+4/k}/pi) int_0^{sqrt2} (1-s^2)^{-2/k}(2-s^2)^{-1/2}
                            [ int_0^{Lam_k(s)} (1+P^{k/2})^{-4/k} dP ] ds,
      Lam_k(s) = ((1-s)/(1+s))^{2/k}.
  All k-dependence of the REGION is the single curve P = Lam_k(s).
""")
    for k in (2, 4, 6):
        def inner(s):
            Lam = ((1 - s) / (1 + s)) ** (mpf(2) / k)
            return quad(lambda P: (1 + P ** (mpf(k) / 2)) ** (-mpf(4) / k), [0, Lam])
        val = (2 ** (1 + mpf(4) / k) / pi) * quad(
            lambda s: (1 - s ** 2) ** (-mpf(2) / k) * (2 - s ** 2) ** mpf(-0.5) * inner(s),
            [0, 1, sqrt(2)])
        rep(f"  k = {k} (P-chart, dps={DPS_QUAD})", val, S_series(k), tol=mpf(10) ** -14)
    print("  (odd k needs the p-chart rho = p^2 so that D has an integer power of p; same identity)")

    print("""
  CHEBYSHEV on the separated P-integral  int (1+P^{k/2})^{-4/k} dP  (m=0, n=k/2, p=-4/k):
      (m+1)/n = 2/k in Z  <=> k in {1,2};  (m+1)/n+p = -2/k in Z <=> k in {1,2};  p in Z <=> k in {1,2,4}.
  So the inner integral is elementary EXACTLY for k = 1, 2, 4.

  Same computation in the polar chart: int_0^1 r dr / sqrt(1 + sigma r^k)  (m=1,n=k,p=-1/2)
  is elementary EXACTLY for k = 1, 2, 4  as well.

  THEREFORE, rigorously (Chebyshev's theorem, no conjecture):
     * k = 3: BOTH iterated integrals are non-elementary, yet S(3) IS elementary.
     * k = 4: the inner integral IS elementary, yet S(4) is NOT.
  "Inner integrand elementary" is neither necessary nor sufficient.  Only the x-first
  slicing of cmd `xdeg` gives the right criterion.""")
    mp.dps = DPS_ID


# ---------------------------------------------------------------- thomae ----
def cmd_thomae():
    import sympy as sp
    import itertools
    print("=" * 78)
    print("THOMAE -- the S_5 orbit of 3F2(1/4,3/4,c; 1, 1+c; 1) and a summability scan")
    print("=" * 78)
    c = sp.Symbol('c')
    start = ((sp.Rational(1, 4), sp.Rational(3, 4), c), (sp.S(1), 1 + c))

    def thomae(params):
        (a, b, cc), (d, e) = params
        s = sp.simplify(d + e - a - b - cc)
        return ((sp.simplify(d - a), sp.simplify(e - a), s),
                (sp.simplify(s + b), sp.simplify(s + cc)))

    def canon(pr):
        (a, b, cc), (d, e) = pr
        return (tuple(sorted(map(str, (sp.simplify(a), sp.simplify(b), sp.simplify(cc))))),
                tuple(sorted(map(str, (sp.simplify(d), sp.simplify(e))))))

    seen = {canon(start): start}
    frontier = [start]
    while frontier:
        nxt = []
        for pr in frontier:
            (a, b, cc), (d, e) = pr
            for np_ in itertools.permutations((a, b, cc)):
                for dp_ in itertools.permutations((d, e)):
                    img = thomae((np_, dp_))
                    key = canon(img)
                    if key not in seen:
                        seen[key] = img
                        nxt.append(img)
        frontier = nxt
    print(f"  distinct 3F2 in the orbit (mod trivial permutations): {len(seen)}")
    for i, (_, pr) in enumerate(sorted(seen.items(), key=lambda kv: str(kv[0]))):
        (a, b, cc), (d, e) = pr
        print(f"   {i+1:2d}. 3F2({a}, {b}, {cc}; {d}, {e}; 1)")

    # numeric spot-check that Thomae is being applied correctly
    mp.dps = 40
    print("\n  numeric check of the basic Thomae step at c = 1/3 and c = 1/5:")
    from mpmath import gamma
    for kk in (3, 5):
        cv = mpf(1) / kk
        lhs = hyp3f2(mpf(1) / 4, mpf(3) / 4, cv, 1, 1 + cv, 1)
        s = 1 + (1 + cv) - mpf(1) / 4 - mpf(3) / 4 - cv
        pref = (gamma(1) * gamma(1 + cv) * gamma(s)
                / (gamma(cv) * gamma(s + mpf(1) / 4) * gamma(s + mpf(3) / 4)))
        rhs = pref * hyp3f2(1 - cv, 1 - cv + cv, s, s + mpf(1) / 4, s + mpf(3) / 4, 1)
        rhs = pref * hyp3f2(1 - cv, 1, s, s + mpf(1) / 4, s + mpf(3) / 4, 1)
        rep(f"    k={kk}", lhs, rhs, tol=mpf(10) ** -30)
    mp.dps = DPS_ID

    print("\n  --- summability scan over the WHOLE orbit; roots in c reported ---")
    hits = {'gauss(numer=denom)': set(), 'terminating': set(),
            'watson': set(), 'whipple': set(), 'dixon': set()}
    for _, pr in seen.items():
        (a, b, cc), (d, e) = pr
        nums, dens = [a, b, cc], [d, e]
        for X in nums:
            for Y in dens:
                for r in sp.solve(sp.Eq(X, Y), c):
                    hits['gauss(numer=denom)'].add(sp.nsimplify(r))
            for m in range(0, 5):
                for r in sp.solve(sp.Eq(X, -m), c):
                    hits['terminating'].add(sp.nsimplify(r))
        for A, B, C in itertools.permutations(nums):
            for Dn, En in itertools.permutations(dens):
                for so in sp.solve([sp.Eq(Dn, (A + B + 1) / 2), sp.Eq(En, 2 * C)], c, dict=True):
                    if c in so:
                        hits['watson'].add(sp.nsimplify(so[c]))
                for so in sp.solve([sp.Eq(A + B, 1), sp.Eq(Dn + En, 2 * C + 1)], c, dict=True):
                    if c in so:
                        hits['whipple'].add(sp.nsimplify(so[c]))
                for so in sp.solve([sp.Eq(Dn, 1 + A - B), sp.Eq(En, 1 + A - C)], c, dict=True):
                    if c in so:
                        hits['dixon'].add(sp.nsimplify(so[c]))
    for kk, vv in hits.items():
        print(f"    {kk:<20s}: c in {sorted(map(str, vv))}")
    print("""
  S(k) sits at c = 1/k.  The only positive-integer k reachable by ANY classical summation
  anywhere on the orbit is k = 1.  (Independent confirmation of lane A2 section 2a.)""")


# --------------------------------------------------------------- schwarz ----
def cmd_schwarz():
    from fractions import Fraction as Fr
    import itertools
    print("=" * 78)
    print("SCHWARZ -- exponent-difference triples, WITH the Schwarz reduction")
    print("=" * 78)
    print("""
  CAUTION.  'angle sum > 1' on the RAW triple is not the Schwarz criterion.  The triple is
  only defined up to  lam -> +-lam + integer  with an EVEN total shift (these moves do not
  change the projective monodromy).  Membership in Schwarz's list must be tested on the
  reduced triple.
""")

    def reduce_triple(t):
        best = None
        for signs in itertools.product([1, -1], repeat=3):
            for sh in itertools.product(range(-3, 4), repeat=3):
                if sum(sh) % 2:
                    continue
                v = tuple(sorted(abs(Fr(signs[i]) * t[i] + sh[i]) for i in range(3)))
                if best is None or (sum(v), v) < (sum(best), best):
                    best = v
        return best

    for label, mk in (("2F1(1, 4/k; 1+2/k; .)   [A1 normal form integrand]",
                       lambda k: (Fr(1), Fr(4, k), 1 + Fr(2, k))),
                      ("2F1(1/2, 2/k; 1+2/k; .) [preamble (v)]",
                       lambda k: (Fr(1, 2), Fr(2, k), 1 + Fr(2, k)))):
        print(f"  {label}")
        print(f"  {'k':>3s}  {'raw triple':<24s}{'sum':>6s}   {'reduced':<24s}{'sum':>6s}  regime")
        for k in range(1, 13):
            a, b, cp = mk(k)
            t = (abs(1 - cp), abs(a - b), abs(cp - a - b))
            r = reduce_triple(t)
            st, sr = sum(t), sum(r)
            reg = "spherical" if sr > 1 else ("EUCLIDEAN" if sr == 1 else "hyperbolic")
            print(f"  {k:3d}  {str(t):<24s}{str(st):>6s}   {str(r):<24s}{str(sr):>6s}  {reg}")
        print()
    print("""  Verdict: the REDUCED triple has angle sum exactly 1 (or is degenerate) for EVERY k,
  including k = 1,2,3.  So the Schwarz/monodromy classification does NOT separate k<=3 from
  k>=4, and the raw-angle-sum statement in lane A1 section 2.5 is not a valid criterion,
  even though its k<=3 / k>=4 verdict happens to be correct.  The correct separator is the
  elementary one in cmd `xdeg`.""")


# ----------------------------------------------------------------- known ----
def cmd_known():
    from mpmath import gamma
    mp.dps = 40
    print("=" * 78)
    print("KNOWN -- which nearby 3F2 families ARE classically evaluable?")
    print("=" * 78)
    print("""
  (i) the REFLECTED Mellin transform is a Gamma quotient, by GAUSS alone:
        int_0^1 x^{c-1} 2F1(1/4,3/4;1;1-x) dx = (1/c) 3F2(1/4,3/4,1;1,1+c;1)
                                              = (1/c) 2F1(1/4,3/4;1+c;1)
                                              = Gamma(c)^2/(Gamma(c+1/4)Gamma(c+3/4)).
      S(k) is the SAME Mellin transform of the UNreflected f: the free parameter sits in
      the pair (c ; 1+c) instead of the pair (1 ; 1+c), and no Gauss reduction applies.""")
    for cv in (mpf(1), mpf(1) / 2, mpf(1) / 3, mpf(7) / 5):
        lhs = quad(lambda x: x ** (cv - 1) * F(1 - x), [0, 1])
        rhs = gamma(cv) ** 2 / (gamma(cv + mpf(1) / 4) * gamma(cv + mpf(3) / 4))
        rep(f"    reflected, c = {mp.nstr(cv,5)}", lhs, rhs, tol=mpf(10) ** -25)
    print("\n    ... and the UNreflected one is NOT that Gamma quotient:")
    for cv in (mpf(1), mpf(1) / 2, mpf(1) / 3):
        lhs = quad(lambda x: x ** (cv - 1) * F(x), [0, 1])
        rhs = gamma(cv) ** 2 / (gamma(cv + mpf(1) / 4) * gamma(cv + mpf(3) / 4))
        print(f"      c={mp.nstr(cv,5)}: {mp.nstr(lhs,18)} vs {mp.nstr(rhs,18)}"
              f"  {'EQUAL' if abs(lhs-rhs)<mpf(10)**-25 else 'DIFFERENT'}")
    print("""
  (ii) k = 1 for the WHOLE signature family, from the exact-derivative property
       d/dz[z(1-z)y'] = lam(1-lam) y  of 2F1(lam,1-lam;1;z)  (valid iff a+b=c=1):
         S_lam(1) = sin(pi lam)/(pi lam (1-lam)).""")
    for lam in (mpf(1) / 2, mpf(1) / 3, mpf(1) / 4, mpf(1) / 6, mpf(2) / 5):
        v = quad(lambda t: hyp2f1(lam, 1 - lam, 1, t), [0, 1])
        rep(f"    lam = {mp.nstr(lam,5)}", v, mpmath.sin(pi * lam) / (pi * lam * (1 - lam)),
            tol=mpf(10) ** -25)
    print("""
  (iii) Legendre identification: 2F1(1/4,3/4;1;z) = P_{-1/4}(1-2z).  The classical
        Mellin formula int_0^1 x^{s-1} P_nu(2x-1) dx = Gamma(s)^2/(Gamma(s-nu)Gamma(s+nu+1))
        is exactly case (i) (argument 2x-1, i.e. REFLECTED).  There is no classical
        evaluation of int_0^1 x^{s-1} P_nu(1-2x) dx, which is S(k).""")
    for z in (mpf(1) / 5, mpf(7) / 10):
        rep(f"    Legendre, z={mp.nstr(z,4)}", F(z), mpmath.legenp(mpf(-1) / 4, 0, 1 - 2 * z),
            tol=mpf(10) ** -25)
    mp.dps = DPS_ID


def main():
    cmds = sys.argv[1:] or ['all']
    if 'all' in cmds:
        cmds = ['verify', 'xdeg', 'k3', 'sep', 'schwarz', 'known', 'thomae']
    table = {'verify': cmd_verify, 'xdeg': cmd_xdeg, 'k3': cmd_k3, 'sep': cmd_sep,
             'thomae': cmd_thomae, 'schwarz': cmd_schwarz, 'known': cmd_known}
    for c in cmds:
        table[c]()
        print()


if __name__ == '__main__':
    main()
