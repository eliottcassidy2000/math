"""b1b_lemma52_slice_places.py - second audit, Lemma 5.2 on the slice map of README 6.2 (e = 1).

Coordinates r = x w, u (checked in b5 / here): X = C^2_(r,u) minus {om = 0}, om = 1 - r u + 3 r^2 = w,
   v = 2 r^3 - u r^2 + 2 r,  x = r/om,  y = u - 3 r,  d_v = (1/(2 om)) d/dr,  d_u = d/du + (r^2/(2 om)) d/dr.
With tau = 1/u, s = v/u the fibre equation is 2 tau r^3 - r^2 + 2 tau r - s = 0.  Places of C(r, u) over the line at
infinity nu = -deg_(u,v):
   w_1 : r = rho + (rho^2 + 1) tau + ...,  rho^2 = -s     (e = 1, f = 2: residue field C(rho) = C(s)(sqrt(-s)))
   w_2 : r = 1/(2 tau) + ...                              (e = 1, f = 1)
(sum of e f = 3 = generic degree).

Part 1 (exact): Laurent series in tau whose coefficients are Laurent polynomials in rho (at w_1) or polynomials in s
(at w_2), with rational coefficients and explicit precision tracking.  For each test element phi of B:
   - the valuations of phi, d_u phi, d_v phi and the inequalities of Lemma 5.2;
   - the formulas of the proof on the series: d_v = tau D~_2, d_u = -tau D~_1 - s tau D~_2, with D~_1 = tau d/dtau and
     D~_2 = d/ds on the coefficients (d rho/ds = -1/(2 rho)).
Part 2 (numerical, independent of Part 1): along rays (u, v) = t (1, s0), t = 10^6 .. 10^12, on each of the three
sheets, growth exponents of phi, d_u phi, d_v phi computed in the README's coordinates (x, y, w).
"""
import random
from fractions import Fraction as Fr
import sympy as sp

random.seed(521)


# ------------------------------------------------------------------ Laurent polynomials in one variable over Q
class LP(dict):
    @staticmethod
    def const(c):
        c = Fr(c)
        return LP({0: c}) if c else LP()

    @staticmethod
    def var(e=1, c=1):
        return LP({e: Fr(c)})

    def __add__(self, o):
        out = LP(self)
        for e, c in o.items():
            v = out.get(e, 0) + c
            if v:
                out[e] = v
            else:
                out.pop(e, None)
        return out

    def __neg__(self):
        return LP({e: -c for e, c in self.items()})

    def __sub__(self, o):
        return self + (-o)

    def __mul__(self, o):
        out = {}
        for e1, c1 in self.items():
            for e2, c2 in o.items():
                v = out.get(e1 + e2, 0) + c1 * c2
                if v:
                    out[e1 + e2] = v
                else:
                    out.pop(e1 + e2, None)
        return LP(out)

    def scale(self, c):
        c = Fr(c)
        return LP({e: a * c for e, a in self.items()}) if c else LP()

    def deriv(self):
        return LP({e - 1: c * e for e, c in self.items() if e})

    def is_monomial(self):
        return len(self) == 1

    def inv_monomial(self):
        (e, c), = self.items()
        return LP({-e: 1 / c})

    def show(self, name):
        if not self:
            return "0"
        return " + ".join("%s*%s^%d" % (c, name, e) if e else "%s" % c for e, c in sorted(self.items()))


class Ser:
    """sum_j c[j] tau^(val + j) + O(tau^(val + len(c))); exact=True means no error term."""

    def __init__(self, val, c, exact=False):
        c = list(c)
        while c and not c[0]:
            c.pop(0)
            val += 1
        self.val, self.c, self.exact = val, c, exact

    def prec(self):
        return float("inf") if self.exact else self.val + len(self.c)

    def coeff(self, e):
        j = e - self.val
        return self.c[j] if 0 <= j < len(self.c) else LP()

    @staticmethod
    def const(a):
        return Ser(0, [a if isinstance(a, LP) else LP.const(a)], exact=True)

    @staticmethod
    def tau(e=1):
        return Ser(e, [LP.const(1)], exact=True)

    def __add__(self, o):
        p = min(self.prec(), o.prec())
        v = min(self.val, o.val)
        if p == float("inf"):
            top = max(self.val + len(self.c), o.val + len(o.c))
            return Ser(v, [self.coeff(e) + o.coeff(e) for e in range(v, top)], exact=True)
        return Ser(v, [self.coeff(e) + o.coeff(e) for e in range(v, int(p))])

    def __neg__(self):
        return Ser(self.val, [-a for a in self.c], self.exact)

    def __sub__(self, o):
        return self + (-o)

    def __mul__(self, o):
        v = self.val + o.val
        if self.exact and o.exact:
            n = len(self.c) + len(o.c) - 1 if self.c and o.c else 0
            ex = True
        else:
            p = min(self.prec() + o.val, o.prec() + self.val)
            n = int(p) - v
            ex = False
        out = [LP() for _ in range(max(n, 0))]
        for i, a in enumerate(self.c):
            if i >= n:
                break
            for j, b in enumerate(o.c):
                if i + j >= n:
                    break
                out[i + j] = out[i + j] + a * b
        return Ser(v, out, ex)

    def scale(self, c):
        if isinstance(c, LP):
            return Ser(self.val, [a * c for a in self.c], self.exact)
        return Ser(self.val, [a.scale(c) for a in self.c], self.exact)

    def inverse(self, nterms=None):
        assert self.c, "inverse of a series that is zero to the known precision"
        lead = self.c[0]
        assert lead.is_monomial(), "leading coefficient is not a monomial: cannot invert in Laurent polynomials"
        li = lead.inv_monomial()
        a = [cc * li for cc in self.c]                # unit 1 + a1 tau + ...
        n = len(a) if not self.exact else (nterms or 16)
        a = a + [LP()] * (n - len(a))
        b = [LP.const(1)]
        for k in range(1, n):
            acc = LP()
            for j in range(1, k + 1):
                acc = acc + a[j] * b[k - j]
            b.append(-acc)
        return Ser(-self.val, [bb * li for bb in b])

    def valuation(self):
        if self.c:
            return self.val
        return None if not self.exact else "zero"

    def dtau_log(self):                               # tau d/dtau
        return Ser(self.val, [a.scale(self.val + j) for j, a in enumerate(self.c)], self.exact)

    def dcoef(self, f):
        return Ser(self.val, [f(a) for a in self.c], self.exact)


def eval_poly(P, rser, user):
    """P: sympy Poly in (r, u) with rational coefficients."""
    cr, cu = {0: Ser.const(1)}, {0: Ser.const(1)}

    def pw(cache, base, n):
        if n not in cache:
            cache[n] = pw(cache, base, n - 1) * base
        return cache[n]
    tot = Ser(0, [], exact=True)
    for (i, j), c in P.terms():
        tot = tot + (pw(cr, rser, i) * pw(cu, user, j)).scale(Fr(int(c.p), int(c.q)))
    return tot


def eval_rat(expr, rs, us, rser, user):
    num, den = sp.fraction(sp.cancel(expr))
    N = eval_poly(sp.Poly(num, rs, us), rser, user)
    if den.is_number:
        return N.scale(1 / Fr(int(sp.Rational(den).p), int(sp.Rational(den).q)))
    D = eval_poly(sp.Poly(den, rs, us), rser, user)
    return N * D.inverse()


def part1():
    r, u = sp.symbols("r u")
    v = 2 * r ** 3 - u * r ** 2 + 2 * r
    om = 1 - r * u + 3 * r ** 2
    X, Y = r / om, u - 3 * r
    wR = sp.cancel(1 / (1 + X * Y))
    assert sp.cancel(wR - om) == 0 and sp.cancel(Y + 3 * X * wR - u) == 0 and sp.cancel(X * wR + X * wR ** 2 - X ** 3 * wR ** 3 - v) == 0
    Vs = sp.Symbol("V")
    print("coordinates: r = x w gives u = y + 3 x w, v = x w + x w^2 - x^3 w^3 = 2 r^3 - u r^2 + 2 r, w = 1 - r u + 3 r^2   [checked]")
    print("discriminant in r of 2 r^3 - u r^2 + 2 r - V: %s   (the quartic of README 6.2 up to the factor -4)" % sp.factor(sp.discriminant(v - Vs, r)))
    du = lambda h: sp.cancel(sp.diff(h, u) + r ** 2 / (2 * om) * sp.diff(h, r))
    dv = lambda h: sp.cancel(sp.diff(h, r) / (2 * om))
    assert du(u) == 1 and du(v) == 0 and dv(v) == 1 and dv(u) == 0

    NT = 16
    # ---- w_2: r = R / tau, 2 R^3 - R^2 + 2 tau^2 R - s tau^2 = 0; coefficients polynomials in s
    s2 = LP.var(1)
    R = Ser(0, [LP.const(Fr(1, 2))] + [LP()] * (NT - 1))
    for _ in range(8):
        G = (R * R * R).scale(2) - R * R + Ser.tau(2) * (R.scale(2) - Ser.const(s2))
        Gp = (R * R).scale(6) - R.scale(2) + Ser.tau(2).scale(2)
        R = R - G * Gp.inverse()
        R = Ser(R.val, (R.c + [LP()] * NT)[:NT])
    G = (R * R * R).scale(2) - R * R + Ser.tau(2) * (R.scale(2) - Ser.const(s2))
    assert not G.c, "w_2 root not converged"
    r_w2 = Ser.tau(-1) * R
    # ---- w_1: r = rho + ..., 2 tau r^3 - r^2 + 2 tau r + rho^2 = 0; coefficients Laurent polynomials in rho
    rho = LP.var(1)
    r1 = Ser(0, [rho] + [LP()] * (NT - 1))
    for _ in range(8):
        G = Ser.tau(1) * (r1 * r1 * r1).scale(2) - r1 * r1 + Ser.tau(1) * r1.scale(2) + Ser.const(rho * rho)
        Gp = Ser.tau(1) * (r1 * r1).scale(6) - r1.scale(2) + Ser.tau(1).scale(2)
        r1 = r1 - G * Gp.inverse()
        r1 = Ser(r1.val, (r1.c + [LP()] * NT)[:NT])
    G = Ser.tau(1) * (r1 * r1 * r1).scale(2) - r1 * r1 + Ser.tau(1) * r1.scale(2) + Ser.const(rho * rho)
    assert not G.c, "w_1 root not converged"
    u_ser = Ser.tau(-1)
    print("place w_1: r = %s + (%s) tau + (%s) tau^2 + ...,  rho^2 = -s" % (r1.coeff(0).show("rho"), r1.coeff(1).show("rho"), r1.coeff(2).show("rho")))
    print("place w_2: r = (%s)/tau + (%s) tau + (%s) tau^3 + ..." % (r_w2.coeff(-1).show("s"), r_w2.coeff(1).show("s"), r_w2.coeff(3).show("s")))

    base = {"x": X, "y": Y, "w": om, "1/w": 1 / om, "r": r, "x y": X * Y, "x w": X * om, "y w": Y * om, "x^2": X ** 2, "w^2": om ** 2,
            "x^2 y w": X ** 2 * Y * om, "x + y": X + Y, "x y - w": X * Y - om, "y^3 + x": Y ** 3 + X, "u v + x": u * v + X,
            "x/w^2": X / om ** 2, "1 + x y w": 1 + X * Y * om, "y - u": Y - u, "2 r - u": 2 * r - u, "w - 3 r^2": om - 3 * r ** 2,
            "x^3 w^5 - y": X ** 3 * om ** 5 - Y, "u": u, "v": v, "u^2 v - x": u ** 2 * v - X}
    for t in range(12):                      # random elements of B = Q[x, y, w]
        e = sum(random.randint(-3, 3) * X ** random.randint(0, 2) * Y ** random.randint(0, 2) * om ** random.randint(0, 2) for _ in range(3))
        if sp.cancel(e) != 0:
            base["rnd%d" % t] = e
    tests = dict(base)
    for name in list(base):
        tests["d_u(%s)" % name] = du(base[name])
        tests["d_v(%s)" % name] = dv(base[name])

    for place, rser, sLP, ds in (("w_1", r1, -(rho * rho), lambda a: a.deriv() * LP.var(-1, Fr(-1, 2))),
                                 ("w_2", r_w2, s2, lambda a: a.deriv())):
        tot = bad = eq = fbad = fchk = skipped = 0
        rows = []
        for name, h in tests.items():
            if sp.cancel(h) == 0:
                continue
            S0 = eval_rat(h, r, u, rser, u_ser)
            Su = eval_rat(du(h), r, u, rser, u_ser)
            Sv = eval_rat(dv(h), r, u, rser, u_ser)
            w0 = S0.valuation()
            if w0 is None:
                skipped += 1
                continue
            for S1, exact_zero in ((Su, sp.cancel(du(h)) == 0), (Sv, sp.cancel(dv(h)) == 0)):
                tot += 1
                if exact_zero:
                    continue
                w1 = S1.valuation()
                if w1 is None:                     # zero to the known precision: valuation >= precision
                    if S1.prec() < w0 + 1:
                        skipped += 1
                    continue
                if w1 < w0 + 1:
                    bad += 1
                elif w1 == w0 + 1:
                    eq += 1
            # the formulas of the proof on the series
            D1 = S0.dtau_log()
            D2 = S0.dcoef(ds)
            pred_v = Ser.tau(1) * D2
            pred_u = -(Ser.tau(1) * D1) - (Ser.tau(1) * D2).scale(sLP)
            for pred, act in ((pred_u, Su), (pred_v, Sv)):
                top = min(pred.prec(), act.prec())
                top = int(top) if top != float("inf") else max(pred.val + len(pred.c), act.val + len(act.c))
                lo = min(pred.val, act.val)
                fchk += 1
                if any((pred.coeff(e) - act.coeff(e)) for e in range(lo, top)):
                    fbad += 1
            if name in ("x", "y", "w", "r", "1/w"):
                rows.append("%s: %s -> d_u %s, d_v %s" % (name, w0, Su.valuation(), Sv.valuation()))
        print("%s: %d inequalities on %d elements: %d violations, %d equalities, %d undecided at this precision; formulas of the proof on the series: %d mismatches of %d"
              % (place, tot, len(tests), bad, eq, skipped, fbad, fchk))
        print("     valuations  " + ";  ".join(rows))


def part2():
    import mpmath as mp
    mp.mp.dps = 80
    xs, ys = sp.symbols("x y")
    wS = 1 / (1 + xs * ys)
    uS = ys + 3 * xs * wS
    vS = xs * wS + xs * wS ** 2 - xs ** 3 * wS ** 3
    J = sp.cancel(sp.diff(uS, xs) * sp.diff(vS, ys) - sp.diff(uS, ys) * sp.diff(vS, xs))
    assert sp.cancel(J + 2 * wS ** 3) == 0
    duS = lambda h: sp.cancel((sp.diff(h, xs) * sp.diff(vS, ys) - sp.diff(h, ys) * sp.diff(vS, xs)) / J)
    dvS = lambda h: sp.cancel((sp.diff(uS, xs) * sp.diff(h, ys) - sp.diff(uS, ys) * sp.diff(h, xs)) / J)
    elems = {"x": xs, "y": ys, "w": wS, "x y": xs * ys, "x w^2": xs * wS ** 2, "y^2 w": ys ** 2 * wS, "x^2 y - w": xs ** 2 * ys - wS,
             "1/w": 1 / wS, "x y w + y": xs * ys * wS + ys}
    s0 = mp.mpf("0.7") + mp.mpf("0.3") * 1j
    worst = None
    lines = []
    viol = 0
    for name, h in elems.items():
        fs = [sp.lambdify((xs, ys), e, "mpmath") for e in (h, duS(h), dvS(h))]
        expo = {}
        for sheet in range(3):
            vals = []
            for t in (mp.mpf(10) ** 8, mp.mpf(10) ** 12):
                uu, vv = t, s0 * t
                roots = mp.polyroots([2, -uu, 2, -vv], maxsteps=200, extraprec=200)
                # order the sheets: two bounded roots (near +-sqrt(-s0)) and one large root
                roots = sorted(roots, key=lambda z: (abs(z) > 10, mp.re(z)))
                rr = roots[sheet]
                om_ = 1 - rr * uu + 3 * rr ** 2
                xx, yy = rr / om_, uu - 3 * rr
                # the point really is a preimage (README formulas)
                ww = 1 / (1 + xx * yy)
                assert abs(yy + 3 * xx * ww - uu) < mp.mpf(10) ** (-20) * abs(uu) and abs(xx * ww + xx * ww ** 2 - xx ** 3 * ww ** 3 - vv) < mp.mpf(10) ** (-20) * abs(vv)
                vals.append([abs(f(xx, yy)) for f in fs])
            ex = [float(mp.log(vals[1][i] / vals[0][i]) / mp.log(mp.mpf(10) ** 4)) if vals[0][i] != 0 and vals[1][i] != 0 else None for i in range(3)]
            expo[sheet] = ex
            for i in (1, 2):
                if ex[i] is not None and ex[0] is not None:
                    gap = ex[0] - ex[i]              # must be >= 1 (growth exponent = -valuation)
                    if worst is None or gap < worst:
                        worst = gap
                    if gap < 1 - 1e-3:
                        viol += 1
        lines.append("%-10s growth exponents (phi, d_u phi, d_v phi) on the sheets: " % name + " | ".join("(%s)" % ", ".join("%.3f" % e if e is not None else "zero" for e in expo[sh]) for sh in range(3)))
    print("numerical rays (u, v) = t (1, 0.7 + 0.3 i), t = 1e8 -> 1e12, sheets 0, 1 (r bounded: place w_1), 2 (r ~ u/2: place w_2):")
    for l in lines:
        print("   " + l)
    print("   smallest drop of the growth exponent under d_u or d_v: %.4f (Lemma 5.2: >= 1); violations: %d" % (worst, viol))


if __name__ == "__main__":
    part1()
    part2()
