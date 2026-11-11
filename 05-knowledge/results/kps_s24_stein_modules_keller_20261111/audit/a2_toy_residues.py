"""a2_toy_residues.py  (independent audit; sympy, exact, symbolic centre (a, b) = (c^2, b))

Proposition 6.1 of the note: on B = C[x, 1/x, y], W = (x^2 - a)(y/x - b), a = c^2 != 0,
    lambda_+-(x^k y^l) = l! Res_(x = +-c) [ x^(k+l+1) e^(b (x^2 - a)) / (x^2 - a)^(l+1) ].
Since the pole has order l+1,  lambda_+(x^k y^l) = (d/dx)^l [ x^(k+l+1) e^(b(x^2-a)) (x + c)^(-l-1) ] at x = c  (and
the mirror formula at x = -c).

Checks:
 1. lambda_+- vanish on (d/dx - W_x)(x^k y^l) and (d/dy - W_y)(x^k y^l) for -7 <= k <= 7, 0 <= l <= 5
 2. lambda_+-(x^(-2) V^m) = (-1)^m m! / (2 a^(m+1)), V = y/x - b, m <= 8
 3. lambda_+-(U^i V^j) = i! [i = j] / 2 for i + j <= 6 (U = x^2 - a); lambda_+-(1) = 1/2, lambda_+-(x) = +- c/2
 4. lambda_+ is one half of the local Gaussian functional at the preimage q_+ = (c, b c) of the centre:
       lambda_+(h) = (1/2) sum_i i! h_ii,   h = sum h_ij U^i V^j the Taylor expansion at q_+   (monomials, |k| <= 5, l <= 4)
    (so the two Stein functionals of the toy are the two local Gaussians; there is no 'extra' functional)
 5. lambda_+-(x^(-1) V^m): odd under the deck involution (lambda_- = - lambda_+), nonzero
"""
import sympy as sp

x, y = sp.symbols("x y")
c, b = sp.symbols("c b", positive=True)
a = c ** 2
W = (x ** 2 - a) * (y / x - b)
Wx, Wy = sp.diff(W, x), sp.diff(W, y)
ok = True


def report(name, cond):
    global ok
    ok = ok and bool(cond)
    print("%-92s %s" % (name, "OK" if cond else "FAIL"))


_cache = {}


def lam_mono(k, l, sign):
    """lambda_sign(x^k y^l), exact."""
    key = (k, l, sign)
    if key not in _cache:
        if sign > 0:
            g = x ** (k + l + 1) * sp.exp(b * (x ** 2 - a)) * (x + c) ** (-(l + 1))
            val = sp.diff(g, x, l).subs(x, c)
        else:
            g = x ** (k + l + 1) * sp.exp(b * (x ** 2 - a)) * (x - c) ** (-(l + 1))
            val = sp.diff(g, x, l).subs(x, -c)
        _cache[key] = sp.simplify(val)
    return _cache[key]


def lam(h, sign):
    """lambda_sign of a Laurent polynomial h in x (polynomial in y) with coefficients in Q(c, b)."""
    h = sp.expand(h)
    tot = 0
    for term in sp.Add.make_args(h):
        coeff, mon = term.as_independent(x, y)
        k = sp.degree(mon * x ** 60, x) - 60 if mon != 1 else 0
        l = sp.degree(mon, y) if mon != 1 else 0
        tot += coeff * lam_mono(int(k), int(l), sign)
    return sp.simplify(tot)


# 1
bad = 0
cnt = 0
for k in range(-7, 8):
    for l in range(0, 6):
        m = x ** k * y ** l
        gx = sp.expand(sp.diff(m, x) - Wx * m)
        gy = sp.expand(sp.diff(m, y) - Wy * m)
        for sgn in (1, -1):
            cnt += 2
            bad += (lam(gx, sgn) != 0) + (lam(gy, sgn) != 0)
report("1. lambda_+- kill (d/dx - W_x) m and (d/dy - W_y) m, m = x^k y^l, |k| <= 7, l <= 5 (%d values)" % cnt, bad == 0)

# 2
V = y / x - b
U = x ** 2 - a
good = True
vals = []
for m in range(0, 9):
    h = sp.expand(x ** (-2) * V ** m)
    target = (-1) ** m * sp.factorial(m) / (2 * a ** (m + 1))
    lp, lm = lam(h, 1), lam(h, -1)
    good = good and sp.simplify(lp - target) == 0 and sp.simplify(lm - target) == 0
    vals.append(str(lp))
report("2. lambda_+-(x^-2 V^m) = (-1)^m m!/(2 a^(m+1)), m <= 8", good)
print("      values of lambda_+ (m = 0..8):", ", ".join(vals))

# 3
good = True
for i in range(0, 7):
    for j in range(0, 7 - i):
        h = sp.expand(U ** i * V ** j)
        target = sp.Rational(sp.factorial(i), 2) if i == j else 0
        good = good and sp.simplify(lam(h, 1) - target) == 0 and sp.simplify(lam(h, -1) - target) == 0
report("3. lambda_+-(U^i V^j) = i! [i = j]/2, i + j <= 6", good)
report("3'. lambda_+-(1) = 1/2, lambda_+-(x) = +- c/2 (independent functionals)",
       lam(sp.Integer(1), 1) == sp.Rational(1, 2) and lam(sp.Integer(1), -1) == sp.Rational(1, 2)
       and sp.simplify(lam(x, 1) - c / 2) == 0 and sp.simplify(lam(x, -1) + c / 2) == 0)

# 4  local Gaussian at q_+: x = sqrt(a + U), y = x (b + V);  h_ii = [U^i V^i] h
Us, Vs = sp.symbols("Us Vs")
good = True
for k in range(-5, 6):
    for l in range(0, 5):
        hloc = (a + Us) ** sp.Rational(k + l, 2) * (b + Vs) ** l
        tot = 0
        for i in range(0, l + 1):
            cV = sp.expand(hloc).coeff(Vs, i) if l > 0 else (hloc if i == 0 else 0)
            # coefficient of Us^i in (a + Us)^((k+l)/2), times the V-coefficient
            cUV = sp.series(cV, Us, 0, i + 1).removeO().coeff(Us, i) if i > 0 else cV.subs(Us, 0)
            tot += sp.factorial(i) * cUV
        tot = sp.simplify(sp.powdenest(sp.simplify(tot), force=True))
        good = good and sp.simplify(lam_mono(k, l, 1) - tot / 2) == 0
report("4. lambda_+ = (1/2) * local Gaussian functional at q_+ on x^k y^l, |k| <= 5, l <= 4", good)

# 5
good = True
for m in range(0, 7):
    h = sp.expand(x ** (-1) * V ** m)
    lp, lm = lam(h, 1), lam(h, -1)
    good = good and sp.simplify(lp + lm) == 0 and lp != 0
report("5. lambda_-(x^-1 V^m) = - lambda_+(x^-1 V^m) != 0, m <= 6", good)

# 6  hand proof that dim S_F <= 2 for a != 0 (and S_F = 0 for a = 0): two explicit elements of M'
Dx = lambda h: sp.diff(h, x) - Wx * h
Dy = lambda h: sp.diff(h, y) - Wy * h
good = True
for i in range(-5, 6):
    for j in range(0, 5):
        # D_y(x^(i-1) y^(j+1)) - D_x(x^i y^j) = 2a x^(i-2) y^(j+1) - (i-j-1) x^(i-1) y^j - 2b x^(i+1) y^j
        lhs = Dy(x ** (i - 1) * y ** (j + 1)) - Dx(x ** i * y ** j)
        rhs = 2 * a * x ** (i - 2) * y ** (j + 1) - (i - j - 1) * x ** (i - 1) * y ** j - 2 * b * x ** (i + 1) * y ** j
        good = good and sp.simplify(lhs - rhs) == 0
    # D_y(x^i) = a x^(i-1) - x^(i+1)
    good = good and sp.simplify(Dy(x ** i) - (a * x ** (i - 1) - x ** (i + 1))) == 0
report("6. reduction identities: y-degree drops by one (a != 0) and x^(i+1) = a x^(i-1) mod M'  => dim S_F <= 2", good)

print("ALL OK" if ok else "SOME CHECK FAILED")
