#!/usr/bin/env python3
"""Verifier (lane digits): ALL integer roots of an integer polynomial, exact, no factorisation and no floating point.

floors(c) returns a superset of { floor(rho) : rho a real root of the polynomial }, by recursion on the derivative:
between consecutive unit intervals that contain the critical points the polynomial is monotone, so a root there is
located by an integer bisection on exact signs.
"""


def pv(c, x):
    v = 0
    for t in c:
        v = v * x + t
    return v


def deriv(c):
    n = len(c) - 1
    return [c[i] * (n - i) for i in range(n)]


def sgn(x):
    return (x > 0) - (x < 0)


def floors(c):
    c = list(c)
    while len(c) > 1 and c[0] == 0:
        c.pop(0)
    deg = len(c) - 1
    if deg <= 0:
        return []
    if deg == 1:
        return [(-c[1]) // c[0]]
    crit = floors(deriv(c))
    B = 2 + max(abs(x) for x in c[1:]) // abs(c[0])         # Cauchy bound: every root has |rho| < 1 + max|c_i / c_0|
    res = set(crit)
    pts = sorted(set([-B] + [t for t in crit if -B < t < B] + [t + 1 for t in crit if -B < t + 1 < B] + [B]))
    critset = set(crit)
    for u, v in zip(pts, pts[1:]):
        if u in critset and v == u + 1:
            continue                                      # unit interval of a critical point: floor u already in res
        su, sv = sgn(pv(c, u)), sgn(pv(c, v))
        if su == 0:
            res.add(u)
        if sv == 0:
            res.add(v)
        if su * sv < 0:
            lo, hi = u, v                                 # monotone on [u, v]: one sign change
            while hi - lo > 1:
                mid = (lo + hi) // 2
                sm = sgn(pv(c, mid))
                if sm == 0:
                    res.add(mid)
                    lo = hi = mid
                    break
                if sm == su:
                    lo = mid
                else:
                    hi = mid
            res.add(lo)
    return sorted(res)


def int_roots(c):
    out = set()
    for t in floors(c):
        for w in (t, t + 1):
            if pv(c, w) == 0:
                out.add(w)
    return sorted(out)


if __name__ == "__main__":
    import random
    from vlin import int_roots_exact
    rng = random.Random(5)
    bad = 0
    for trial in range(4000):
        deg = rng.choice([2, 3, 4, 5, 6])
        # product of random linear / quadratic integer factors, so that integer roots occur often (also multiple ones)
        c = [1]
        while len(c) - 1 < deg:
            if rng.random() < 0.6 or len(c) - 1 == deg - 1:
                r = rng.choice([rng.randint(-9, 9), rng.randint(-2000, 2000), rng.randint(-10 ** 6, 10 ** 6)])
                f = [1, -r]
            else:
                f = [rng.randint(1, 3), rng.randint(-50, 50), rng.randint(-500, 500)]
            new = [0] * (len(c) + len(f) - 1)
            for i, x in enumerate(c):
                for j, y in enumerate(f):
                    new[i + j] += x * y
            c = new
        a = set(int_roots(c))
        b = set(int_roots_exact(c)) if abs(c[-1]) < 10 ** 14 else None
        direct = None
        if b is not None and a != b:
            bad += 1
    print("int_roots (monotone bisection) = rational-root-theorem roots on 4000 random products:", bad == 0)
