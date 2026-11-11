"""b5_slice_which_generator.py - second audit, README 6.2 "What is proved on the slice" (ii).

The README proves: for AT LEAST ONE of the generators h = x, y, w the iterates d_u^m h never vanish (and the same for
d_v).  Here: for ALL THREE generators and both directions, by the order of pole along the puncture curve.

Coordinates r = x w, u:  X = C^2_(r,u) minus Omega, Omega = {om = 0}, om = 1 - r u + 3 r^2 = w,
   v = 2 r^3 - u r^2 + 2 r,   d_v = (1/(2 om)) d/dr,   d_u = d/du + (r^2/(2 om)) d/dr,
   x = r / om,  y = u - 3 r,  w = om.   Omega is the critical curve of the finite map (r, u) -> (u, v); its image is
the non-properness curve of the slice map.

Lemma (pole growth).  If phi = psi / om^n with n >= 1 and om not dividing psi, then
   d_u phi = [ -n psi r (r u - 2)/2 + om (...) ] / om^(n+2),     d_v phi = [ -n psi (6 r - u)/2 + om (...) ] / om^(n+2),
and r (r u - 2), 6 r - u are not divisible by om; so the order of pole along Omega grows by exactly 2.
Consequently d_u^m x and d_v^m x have a pole of order 2m + 1, d_u^m y, d_v^m y, d_u^m w, d_v^m w a pole of order
2m - 1 (m >= 1): none of them is zero.
"""
import sympy as sp

r, u = sp.symbols("r u")
om = 1 - r * u + 3 * r ** 2
v = 2 * r ** 3 - u * r ** 2 + 2 * r
du = lambda h: sp.cancel(sp.diff(h, u) + r ** 2 / (2 * om) * sp.diff(h, r))
dv = lambda h: sp.cancel(sp.diff(h, r) / (2 * om))
assert du(u) == 1 and du(v) == 0 and dv(v) == 1 and dv(u) == 0


def pole_order(expr):
    """order of pole along om = 0 (negative = order of zero) of a rational function whose denominator is a power of om."""
    num, den = sp.fraction(sp.cancel(expr))
    den = sp.Poly(den, r, u)
    omp = sp.Poly(om, r, u)
    n = 0
    while True:
        qq, rem = sp.div(den, omp)
        if rem.is_zero:
            den = qq
            n += 1
        else:
            break
    assert den.total_degree() == 0, "denominator is not a power of om: %s" % den
    num = sp.Poly(num, r, u)
    z = 0
    while n == 0 and not num.is_zero:
        qq, rem = sp.div(num, omp)
        if rem.is_zero:
            num = qq
            z += 1
        else:
            break
    return n - z


if __name__ == "__main__":
    # the two residues modulo om of the lemma
    usub = (1 + 3 * r ** 2) / r                       # u on Omega
    print("on Omega: r (r u - 2) = %s,   6 r - u = %s   (both nonzero)" % (sp.factor(sp.cancel((r * (r * u - 2)).subs(u, usub))), sp.factor(sp.cancel((6 * r - u).subs(u, usub)))))
    print("d_u om = %s,   d_v om = %s" % (sp.factor(du(om)), sp.factor(dv(om))))
    # symbolic check of the lemma for generic psi, n
    n = sp.symbols("n", positive=True, integer=True)
    psi = sp.Function("psi")(r, u)
    phi = psi / om ** n
    lead_u = sp.simplify((sp.diff(phi, u) + r ** 2 / (2 * om) * sp.diff(phi, r)) * om ** (n + 2) + n * psi * r * (r * u - 2) / 2)
    lead_v = sp.simplify((sp.diff(phi, r) / (2 * om)) * om ** (n + 2) + n * psi * (6 * r - u) / 2)
    print("d_u(psi/om^n) om^(n+2) + n psi r (r u - 2)/2 = %s" % sp.factor(lead_u))
    print("d_v(psi/om^n) om^(n+2) + n psi (6 r - u)/2   = %s" % sp.factor(lead_v))
    gens = {"x": r / om, "y": u - 3 * r, "w": om, "x y": (r / om) * (u - 3 * r), "1/w": 1 / om}
    for name, h in gens.items():
        row_u, row_v = [], []
        a, b = h, h
        for m in range(1, 7):
            a, b = du(a), dv(b)
            row_u.append(pole_order(a) if a != 0 else "ZERO")
            row_v.append(pole_order(b) if b != 0 else "ZERO")
        print("h = %-4s pole order along Omega: h: %2d;  d_u^m h, m = 1..6: %s;  d_v^m h: %s" % (name, pole_order(h), row_u, row_v))
    for name, h in (("u", u), ("v", v), ("u v^2", u * v ** 2)):
        a, cnt = h, 0
        while a != 0 and cnt < 10:
            a = du(a)
            cnt += 1
        print("control: d_u^m (%s) = 0 for m = %d" % (name, cnt))
