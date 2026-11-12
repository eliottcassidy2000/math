"""a2_membership_over_Q.py - audit of Theorem A2 (2): explicit membership certificates over Q.

For W = x + x^(l+1) y (l = 1, 2, 3) and a target polynomial T, solve EXACTLY over Q the linear system

    (d/dx - W_x) p + (d/dy - W_y) q = T          (equivalently div(v) - v(W) = T, v = p d/dx + q d/dy)

for polynomials p, q of bounded bidegree, then re-expand the left side and compare with T coefficient by coefficient.
No use of the author's engine (no prime, no truncation relative membership: a solution is a certificate).

Targets: W^m (m >= l), W^m - (1-l)...(m-l) (m < l), x^l W^n - n! x^l.
Controls (must have NO solution in the search box; non-membership itself is proved by the linear form lambda):
1, x^l, x^l W^n.
Also: the certificate predicted by Lemma A1, v = -h xi with xi = x d/dx - l y d/dy, for W h - (xi + delta) h.
"""
from fractions import Fraction
from math import factorial


def add(P, Q, s=1):
    R = dict(P)
    for m, c in Q.items():
        R[m] = R.get(m, 0) + s * c
        if R[m] == 0:
            del R[m]
    return R


def mul(P, Q):
    R = {}
    for (a, b), c in P.items():
        for (a2, b2), c2 in Q.items():
            m = (a + a2, b + b2)
            R[m] = R.get(m, 0) + c * c2
    return {m: c for m, c in R.items() if c != 0}


def scal(P, s):
    return {m: c * s for m, c in P.items() if c * s != 0}


def dx(P):
    return {(k - 1, j): c * k for (k, j), c in P.items() if k > 0}


def dy(P):
    return {(k, j - 1): c * j for (k, j), c in P.items() if j > 0}


def power(P, n):
    R = {(0, 0): Fraction(1)}
    for _ in range(n):
        R = mul(R, P)
    return R


def Wpoly(l):
    return {(1, 0): Fraction(1), (l + 1, 1): Fraction(1)}


def image(p, q, W):
    """(d/dx - W_x) p + (d/dy - W_y) q"""
    return add(add(dx(p), mul(dx(W), p), -1), add(dy(q), mul(dy(W), q), -1))


def solve(T, W, DX, DY):
    """find p, q with monomials x^a y^b, a <= DX, b <= DY, such that image(p, q) = T.  Exact Gauss-Jordan over Q.
    returns (p, q) or None."""
    unknowns = [("p", a, b) for a in range(DX + 1) for b in range(DY + 1)] + \
               [("q", a, b) for a in range(DX + 1) for b in range(DY + 1)]
    cols = []
    rows = set(T.keys())
    for (kind, a, b) in unknowns:
        mono = {(a, b): Fraction(1)}
        im = image(mono, {}, W) if kind == "p" else image({}, mono, W)
        cols.append(im)
        rows |= set(im.keys())
    rows = sorted(rows)
    ridx = {m: i for i, m in enumerate(rows)}
    nr, nc = len(rows), len(unknowns)
    M = [[Fraction(0)] * (nc + 1) for _ in range(nr)]
    for c, im in enumerate(cols):
        for m, v in im.items():
            M[ridx[m]][c] = v
    for m, v in T.items():
        M[ridx[m]][nc] = v
    piv = []
    r = 0
    for c in range(nc):
        pr = None
        for i in range(r, nr):
            if M[i][c] != 0:
                pr = i
                break
        if pr is None:
            continue
        M[r], M[pr] = M[pr], M[r]
        inv = 1 / M[r][c]
        M[r] = [v * inv for v in M[r]]
        for i in range(nr):
            if i != r and M[i][c] != 0:
                f = M[i][c]
                M[i] = [vi - f * vr for vi, vr in zip(M[i], M[r])]
        piv.append(c)
        r += 1
        if r == nr:
            break
    for i in range(r, nr):
        if M[i][nc] != 0:
            return None
    p, q = {}, {}
    for i, c in enumerate(piv):
        v = M[i][nc]
        if v != 0:
            kind, a, b = unknowns[c]
            (p if kind == "p" else q)[(a, b)] = v
    return p, q


def show(P, limit=6):
    items = sorted(P.items())
    s = " + ".join("(%s) x^%d y^%d" % (c, k, j) for (k, j), c in items[:limit])
    if len(items) > limit:
        s += " + ... (%d terms)" % len(items)
    return s if s else "0"


def bidegree(P):
    return (max(k for k, j in P), max(j for k, j in P)) if P else (0, 0)


def main():
    one = {(0, 0): Fraction(1)}
    for l in (1, 2, 3):
        W = Wpoly(l)
        xl = {(l, 0): Fraction(1)}
        print("=== l = %d, W = x + x^%d y ===" % (l, l + 1))
        targets = []
        for m in range(1, 5 if l < 3 else 5):
            c = Fraction(1)
            for t in range(m):
                c *= (1 - l + t)
            Wm = power(W, m)
            targets.append(("W^%d - (%s)" % (m, c) if c != 0 else "W^%d" % m, add(Wm, scal(one, c), -1), True))
        for n in range(1, 4 if l == 3 else 5):
            Wn = power(W, n)
            targets.append(("x^%d W^%d - %d! x^%d" % (l, n, n, l), add(mul(xl, Wn), scal(xl, Fraction(factorial(n))), -1), True))
        targets.append(("1", one, False))
        targets.append(("x^%d" % l, xl, False))
        targets.append(("x^%d W^2" % l, mul(xl, power(W, 2)), False))
        for name, T, expect in targets:
            bx, by = bidegree(T)
            sol = solve(T, W, bx + 1, by + 1)
            if sol is None:
                # control: enlarge once to make the negative answer less dependent on the box
                sol = solve(T, W, bx + 3, by + 3) if not expect else None
            if sol is None:
                print("  %-22s : no solution in the box (deg_x <= %d, deg_y <= %d)%s" % (name, bx + (1 if expect else 3), by + (1 if expect else 3), "" if not expect else "   <-- UNEXPECTED"))
                continue
            p, q = sol
            chk = add(image(p, q, W), T, -1)
            print("  %-22s : v = p d/dx + q d/dy found, bidegrees p %s, q %s, %d + %d terms; div(v) - v(W) - target = %s%s"
                  % (name, bidegree(p), bidegree(q), len(p), len(q), show(chk), "" if expect else "   <-- UNEXPECTED (control)"))
            if l == 1 and name in ("W^1", "x^1 W^1 - 1! x^1"):
                print("      p =", show(p, 12))
                print("      q =", show(q, 12))
        # Lemma A1 certificate: W h - (xi + delta) h = -(div(v) - v(W)) with v = h xi
        bad = tot = 0
        for k in range(0, 6):
            for j in range(0, 4):
                h = {(k, j): Fraction(1)}
                p = scal(mul(h, {(1, 0): Fraction(1)}), -1)           # v = -h xi, xi = x d/dx - l y d/dy
                q = scal(mul(h, {(0, 1): Fraction(1)}), l)
                lhs = image(p, q, W)
                rhs = add(mul(W, h), scal(h, k - l * j + 1 - l), -1)
                tot += 1
                bad += bool(add(lhs, rhs, -1))
        print("  Lemma A1 certificate v = -h xi gives W h - (k - l j + 1 - l) h exactly: %d monomials, %d failures" % (tot, bad))


if __name__ == "__main__":
    main()
