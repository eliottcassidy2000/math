"""b2_lemma53_torsion.py - second audit, Lemma 5.3 (torsion-freeness of B over C[E_a, E_b]) by brute force.

For a polynomial P(a, b) the map h -> P(E_a, E_b) h, E_a = u - d_v, E_b = v - d_u, is restricted to a box V of basis
monomials and its rank is computed modulo the prime 2^61 - 1.  rank = dim V is a proof of injectivity on V over Q
(the matrices have rational entries with denominators powers of 2, and rank mod p <= rank over Q).
Any rank defect would be a torsion candidate and is printed.

Maps:  toy    = double cover (x^2, y/x) on C^* x C               (ring L, x in Z, y >= 0)
       open   = inclusion C^* x C -> C^2                          (ring L)
       torus  = inclusion C^* x C^* -> C^2                        (ring L, both exponents in Z)
       slice  = README 6.2, e = 1, on C^2 minus {1 + xy = 0}      (ring H)
       hyp    = README 6.3, inclusion of C^2 minus {1 + xy = 0}   (ring H)
Also: the commutator [E_a, E_b] = 0 on every box, and (double cover) the exact shift of the monomial valuation
w(P(E_a, E_b) phi) = w(phi) - deg P.
"""
import random
import sys
from fractions import Fraction
import eng2
from eng2 import RingL, RingH, add, scal, lincomb, power, dual_fields, rank_rows, flint_rank, inv

random.seed(53)
P = eng2.PRIME
HALF = inv(2)


def setup(name):
    if name == "toy":
        R = RingL()
        x, y = R.mono(1, 0), R.mono(0, 1)
        u, v = R.mono(2, 0), R.mono(-1, 1)
        Jinv = scal(R.one(), HALF)
        box = lambda N: [(k, l) for k in range(-N, N + 1) for l in range(0, N + 1) if abs(k) + l <= N]
    elif name == "open":
        R = RingL()
        u, v = R.mono(1, 0), R.mono(0, 1)
        Jinv = R.one()
        box = lambda N: [(k, l) for k in range(-N, N + 1) for l in range(0, N + 1) if abs(k) + l <= N]
    elif name == "torus":
        R = RingL()
        u, v = R.mono(1, 0), R.mono(0, 1)
        Jinv = R.one()
        box = lambda N: [(k, l) for k in range(-N, N + 1) for l in range(-N, N + 1) if abs(k) + abs(l) <= N]
    elif name == "slice":
        R = RingH()
        x, y = R.x(), R.y()
        w = R.q(-1)
        u = add(y, scal(R.mul(x, w), 3))
        v = lincomb([(1, R.mul(x, w)), (1, R.mul(x, R.q(-2))), (-1, R.mul(power(R, x, 3), R.q(-3)))])
        Jinv = scal(R.q(3), -HALF)                   # J = -2 w^3
        box = lambda N: [(c, n) for c in range(-N, N + 1) for n in range(-N, N + 1) if abs(c) + 2 * abs(n) <= N]
    elif name == "hyp":
        R = RingH()
        u, v = R.x(), R.y()
        Jinv = R.one()
        box = lambda N: [(c, n) for c in range(-N, N + 1) for n in range(-N, N + 1) if abs(c) + 2 * abs(n) <= N]
    else:
        raise ValueError(name)
    du, dv, J = dual_fields(R, u, v, Jinv)
    Ea = lambda h: add(R.mul(u, h), dv(h), -1)
    Eb = lambda h: add(R.mul(v, h), du(h), -1)
    return R, u, v, du, dv, Ea, Eb, box


def apply_poly(poly, Ea, Eb, h):
    """poly: dict {(i, j): coeff}; returns sum coeff E_a^i E_b^j h."""
    maxi = max(i for i, j in poly)
    maxj = max(j for i, j in poly)
    # table of E_b^j h, then E_a^i of each
    out = {}
    hb = h
    col = []
    for j in range(maxj + 1):
        col.append(hb)
        hb = Eb(hb)
    for j in range(maxj + 1):
        cur = col[j]
        for i in range(maxi + 1):
            c = poly.get((i, j))
            if c:
                out = add(out, scal(cur, c))
            cur = Ea(cur)
    return out


def polys():
    r = lambda: random.randint(-9, 9)
    out = [("E_a", {(1, 0): 1}), ("E_a - 1", {(1, 0): 1, (0, 0): -1}), ("E_a - 2", {(1, 0): 1, (0, 0): -2}),
           ("E_a - 12345", {(1, 0): 1, (0, 0): -12345}),
           ("E_b", {(0, 1): 1}), ("E_b - 3", {(0, 1): 1, (0, 0): -3}), ("E_b + 7/2", {(0, 1): 1, (0, 0): 7 * HALF % P}),
           ("E_a E_b - 1", {(1, 1): 1, (0, 0): -1}), ("E_a E_b", {(1, 1): 1}), ("E_a^2 + E_b", {(2, 0): 1, (0, 1): 1}),
           ("E_a - E_b", {(1, 0): 1, (0, 1): -1}), ("E_a^2 - E_b^2 - 1", {(2, 0): 1, (0, 2): -1, (0, 0): -1}),
           ("E_a^3 - E_b^2", {(3, 0): 1, (0, 2): -1}), ("E_a^2 E_b - E_a + 2", {(2, 1): 1, (1, 0): -1, (0, 0): 2})]
    for t in range(4):
        q = {(i, j): r() for i in range(3) for j in range(3) if i + j <= 2}
        if not any(q[m] for m in q if m[0] + m[1] == 2):
            q[(1, 1)] = 1
        out.append(("random quadratic #%d" % (t + 1), {m: c % P for m, c in q.items() if c % P}))
    return out


def run(name, N):
    R, u, v, du, dv, Ea, Eb, box = setup(name)
    B = box(N)
    # commutation of E_a, E_b on the box
    bad = 0
    for m in B[:: max(1, len(B) // 60)]:
        h = {m: 1}
        bad += bool(add(Ea(Eb(h)), Eb(Ea(h)), -1))
    print("%s: box of %d monomials (N = %d); [E_a, E_b] h = 0 on a sample: %d failures" % (name, len(B), N, bad))
    defects = 0
    for label, poly in polys():
        rows = [apply_poly(poly, Ea, Eb, {m: 1}) for m in B]
        rk = flint_rank(rows) if len(B) > 400 else rank_rows(rows)
        flag = "" if rk == len(B) else "   <== RANK DEFECT %d (torsion candidate)" % (len(B) - rk)
        defects += rk != len(B)
        print("   P = %-22s rank %d of %d%s" % (label, rk, len(B), flag))
        sys.stdout.flush()
    return defects


def toy_valuation():
    """double cover: w(x^k y^l) = -(k + 3 l)/2 is the extension of -deg_(u,v); check w(P(E) phi) = w(phi) - deg P."""
    R, u, v, du, dv, Ea, Eb, box = setup("toy")
    wv = lambda h: min(Fraction(-(k + 3 * l), 2) for (k, l) in h) if h else None
    bad = tot = 0
    for label, poly in polys():
        n = max(i + j for (i, j) in poly)
        for _ in range(40):
            h = {}
            for _t in range(random.randint(1, 4)):
                h = add(h, R.mono(random.randint(-6, 6), random.randint(0, 5), random.randint(1, 50)))
            if not h:
                continue
            g = apply_poly(poly, Ea, Eb, h)
            tot += 1
            if g == {} or wv(g) != wv(h) - n:
                bad += 1
    print("toy: w(P(E_a, E_b) phi) = w(phi) - deg P for the monomial valuation: %d cases, %d failures" % (tot, bad))


if __name__ == "__main__":
    total = 0
    total += run("toy", 30)
    total += run("open", 26)
    total += run("torus", 20)
    total += run("hyp", 26)
    total += run("slice", 22)
    toy_valuation()
    print("TOTAL rank defects: %d" % total)
