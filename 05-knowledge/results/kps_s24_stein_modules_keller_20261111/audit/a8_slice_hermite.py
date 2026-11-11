"""a8_slice_hermite.py - Proposition 2.2 (d) and the Hermite classes on the slice (centre (2, 3)), independent engine.

 1. [d_u^m d_v^n h] = [h H_(n,m)(U, V)] for h = w, x, y and 1 <= m + n <= 3.  The dual fields are
        d_u h = (h_x v_y - h_y v_x)/J,  d_v h = (u_x h_y - u_y h_x)/J,  1/J = -q^3/2;
    to stay in exact integer arithmetic the check is done for D_u = 2 d_u, D_v = 2 d_v:
        [D_u^m D_v^n h] = 2^(m+n) [h H_(n,m)(U, V)].
 2. Hermite classes [h H_(n,m)(U, V)], n, m <= 4, for h = w, x, 1 (N = nonzero normal form, 0 = zero).
"""
from math import comb, factorial
from tw import RingC, add, mul, scal, power, deriv, build, codim, maxwt, P1

one = {(0, 0): 1}
x, y, w, q = {(1, 0): 1}, {(-1, 0): 1}, {(0, -1): 1}, {(0, 1): 1}
R = RingC(1, 2, 1)
pw = lambda el, k: power(R, el, k, one)
u = add(y, scal(mul(R, x, w), 3))
v = add(add(mul(R, x, w), mul(R, x, pw(w, 2))), mul(R, pw(x, 3), pw(w, 3)), -1)
J = scal(pw(w, 3), -2)
U = add(u, scal(one, 2), -1)
V = add(v, scal(one, 3), -1)
W = mul(R, U, V)
ux, uy, vx, vy = deriv(R, u, "x"), deriv(R, u, "y"), deriv(R, v, "x"), deriv(R, v, "y")
mq3 = scal(pw(q, 3), -1)


def Du(h):      # 2 d_u h = -q^3 (h_x v_y - h_y v_x)
    return mul(R, mq3, add(mul(R, deriv(R, h, "x"), vy), mul(R, deriv(R, h, "y"), vx), -1))


def Dv(h):      # 2 d_v h = -q^3 (u_x h_y - u_y h_x)
    return mul(R, mq3, add(mul(R, ux, deriv(R, h, "y")), mul(R, uy, deriv(R, h, "x")), -1))


assert Du(u) == scal(one, 2) and Dv(v) == scal(one, 2) and not Du(v) and not Dv(u), "dual fields"


def hermite(n, m):
    out = {}
    for j in range(min(n, m) + 1):
        out = add(out, scal(mul(R, pw(U, n - j), pw(V, m - j)), (-1) ** j * factorial(j) * comb(n, j) * comb(m, j)))
    return out


Np = 70
E = build(R, W, Np, P1)
print("N' = %d: %d pivots; c(20) = %d, c(60) = %d" % (Np, len(E.piv), codim(R, E, 20), codim(R, E, 60)))

tot = bad = skipped = 0
for name, h in (("w", w), ("x", x), ("y", y)):
    for m in range(0, 4):
        for n in range(0, 4 - m):
            if m + n == 0:
                continue
            d = h
            for _ in range(n):
                d = Dv(d)
            for _ in range(m):
                d = Du(d)
            lhs = mul(R, d, J)
            rhs = scal(mul(R, mul(R, h, hermite(n, m)), J), 2 ** (m + n))
            diff = add(lhs, rhs, -1)
            if maxwt(R, diff) > Np - 10:
                skipped += 1
                continue
            tot += 1
            bad += bool(E.nf(diff))
print("1. [d_u^m d_v^n h] = [h H_(n,m)(U,V)], h = w, x, y, 1 <= m + n <= 3: %d checked, %d failures, %d skipped (weight)" % (tot, bad, skipped))

for name, h in (("w", w), ("x", x), ("1", one)):
    print("2. Hermite classes of h = %s (rows n = 0..4, columns m = 0..4):" % name)
    for n in range(5):
        row = []
        for m in range(5):
            el = mul(R, mul(R, h, hermite(n, m)), J)
            row.append("?" if maxwt(R, el) > Np - 10 else ("N" if E.nf(el) else "0"))
        print("     n = %d: %s" % (n, " ".join(row)))
