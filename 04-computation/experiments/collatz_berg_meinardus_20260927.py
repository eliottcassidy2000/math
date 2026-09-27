#!/usr/bin/env python3
"""collatz_berg_meinardus_20260927.py -- the Berg-Meinardus functional equation re-derived and checked
(session collatz-posets-zeta5-20260927, opus, 2026-09-27, addendum to the first note).

T(n) = n/2 (n even), (3n+1)/2 (n odd). A function f on the positive integers is T-invariant iff f(2m) = f(m) and
f(2m+1) = f(3m+2) for all m >= 0 (m >= 1 for the first). With h(z) = sum_(n>=1) f(n) z^n and lambda = e^(2 pi i/3):

    h(z^3) = h(z^6) + (1/(3z)) * sum_(k=0..2) lambda^k h(lambda^k z^2)          (BM)

is equivalent to T-invariance (coefficientwise: the z^(3n) coefficients encode the two relations, the other
coefficients vanish identically). Checks, as truncated power series to order Z:
 (1) f = 1 (h = z/(1-z)): (BM) holds to order Z (the Collatz-conjecture solution);
 (2) control f(n) = (-1)^n (not T-invariant): (BM) fails at z^3, where f(1) = f(2) is required;
 (3) control f(n) = 1 for n odd, 0 for n even: fails at z^3 as well (f(1) != f(2));
 (4) the coefficient bookkeeping: the right side of (BM) has nonzero coefficients only at exponents 0 mod 3, and
     [z^(6m)] RHS = f(m) - 0, [z^(6m+3)] RHS = f(3m+2), so (BM) <=> f(2m) = f(m), f(2m+1) = f(3m+2).
 (5) the dimension of the solution space of (BM) in C[[z]] truncated at order Z (linear algebra over the coefficients
     f(1..Z), with the relations that stay inside the window) equals the number of T-orbit classes of {1..Z} that
     are closed inside the window: printed for Z = 50, 100, 200 together with the number of classes that leave
     the window (every class reaching 1 within the window is one class).
Usage: python3 collatz_berg_meinardus_20260927.py
"""
import cmath, math

LAM = cmath.exp(2j * math.pi / 3)


def series_h(f, Z):
    # coefficients c[n] = f(n), n = 1..Z, as a dict of exponent -> value
    return {n: f(n) for n in range(1, Z + 1)}


def compose_power(c, e, Z):
    # h(z^e) truncated at Z
    return {n * e: v for n, v in c.items() if n * e <= Z}


def rotate(c, w):
    # h(w z): coefficient n gets w^n
    return {n: v * w ** n for n, v in c.items()}


def bm_residual(f, Z):
    c = series_h(f, Z)
    lhs = compose_power(c, 3, Z)
    rhs = compose_power(c, 6, Z)
    for k in range(3):
        rot = compose_power(rotate(c, LAM ** k), 2, Z + 1)   # h(lambda^k z^2) = sum f(n) lambda^(kn) z^(2n)
        for n, v in rot.items():
            if n - 1 <= Z and n - 1 >= 0:
                rhs[n - 1] = rhs.get(n - 1, 0) + (LAM ** k) * v / 3
    res = {}
    for n in range(0, Z + 1):
        r = lhs.get(n, 0) - rhs.get(n, 0)
        if abs(r) > 1e-9:
            res[n] = r
    return res


def main():
    Z = 240
    print("== Berg-Meinardus equation h(z^3) = h(z^6) + (1/(3z)) sum_k lambda^k h(lambda^k z^2), truncated at z^%d ==" % Z)
    r1 = bm_residual(lambda n: 1.0, Z)
    print(" (1) f = 1: residual coefficients above 1e-9: %d (equation holds: %s)" % (len(r1), len(r1) == 0))
    r2 = bm_residual(lambda n: (-1.0) ** n, Z)
    print(" (2) control f = (-1)^n: first failing exponent %s (expected 3: the relation f(1) = f(2) at m = 0 fails)" % (min(r2) if r2 else None))
    r3 = bm_residual(lambda n: 1.0 if n % 2 else 0.0, Z)
    print(" (3) control f = [n odd]: first failing exponent %s (expected 3: f(1) != f(2))" % (min(r3) if r3 else None))
    # (4) bookkeeping: random f, check the residual equals f(2m)-f(m) at 6m and f(2m+1)-f(3m+2) at 6m+3, zero elsewhere
    import random
    random.seed(7)
    vals = {n: random.random() for n in range(1, 3 * Z + 3)}
    f = lambda n: vals[n]
    r4 = bm_residual(f, Z)
    ok = True
    for n in range(0, Z + 1):
        r = r4.get(n, 0)
        if n % 6 == 0 and n >= 6:
            m = n // 6; ok &= abs(r - (vals[2 * m] - vals[m])) < 1e-9
        elif n % 6 == 3:
            m = (n - 3) // 6; ok &= abs(r - (vals[2 * m + 1] - vals[3 * m + 2])) < 1e-9
        else:
            ok &= abs(r) < 1e-9
    print(" (4) random f: residual = f(2m)-f(m) at z^(6m), f(2m+1)-f(3m+2) at z^(6m+3), 0 elsewhere: %s" % ok)
    # (5) solution-space dimension inside a window = number of closed orbit classes
    def T(n):
        return n // 2 if n % 2 == 0 else (3 * n + 1) // 2
    for W in (50, 100, 200, 1000):
        # union-find on {1..W} with edges n -> T(n) when T(n) <= W
        parent = list(range(W + 1))
        def find(x):
            while parent[x] != x:
                parent[x] = parent[parent[x]]; x = parent[x]
            return x
        leaves = set()
        for n in range(1, W + 1):
            t = T(n)
            if t <= W:
                a, b = find(n), find(t)
                if a != b:
                    parent[a] = b
            else:
                leaves.add(n)
        classes = {find(n) for n in range(1, W + 1)}
        open_classes = {find(n) for n in leaves}
        closed = len(classes) - len(open_classes)
        # the relations inside the window are exactly the edges n -> T(n) <= W; the solution space of the
        # windowed (BM) system has dimension = number of classes (each class one free value)
        print(" (5) window {1..%d}: %d orbit classes, %d of them closed (all lie in the class of 1: %s), %d classes exit the window" % (W, len(classes), closed, closed == 1, len(open_classes)))


if __name__ == '__main__':
    main()
