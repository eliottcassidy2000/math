"""a3_no_euler_field_example.py - audit of the remark 'What this does not settle: rank-one potentials without an
Euler field, such as x (x^2 y + x + 1)' (README A2, HYP-9668 (Q)).

For W = x^(l+1) y + r(x) (r any polynomial) the residue form of Theorem A2 still works:

    lambda_r(x^k y^j) = j! * [coefficient of x^((l+1)(j+1) - 1 - k) in exp(-r(x))].

Checks (exact rational arithmetic)
  E1  lambda_r kills both families of generators of M_W (monomials k < 14, j < 5) for
      W1 = x^3 y + x^2 + x  (= x (x^2 y + x + 1)),  W2 = x^3 y + x + 1/3 x^2... (a second r),  W3 = x^2 y + x (Broughton)
  E2  one-dimensionality on a box: x^k y^j - (lambda(x^k y^j)/lambda(1)) lies in M_W, with an explicit vector
      field over Q (k <= 4, j <= 2)
  E3  moments lambda(W^m) and lambda(x^l W^n); W1: lambda(W) = lambda(1) and lambda(W^2) = 0, so NO vector field with
      xi(W) = W and constant divergence exists (it would force [W^2] = delta (delta + 1) [1] with delta = 1)
  E4  explicit certificates over Q:  W1^2 in M_W,  W1^4 in M_W,  x^2 W1^n - n! x^2 in M_W (n = 1, 2, 3)
      => M_W is NOT a Mathieu-Zhao subspace for W1 = x (x^2 y + x + 1)  (P = W1^2, h = x^2).
"""
from fractions import Fraction
from math import factorial

from a2_membership_over_Q import add, mul, scal, dx, dy, power, image, solve, show, bidegree


def exp_series(r, order):
    """coefficients of exp(-r(x)) up to x^order; r = dict {deg: coeff} with r(0) = 0 allowed or not"""
    # exp(-r) = exp(-r0) * exp(-(r - r0)); keep r0 = 0 here (all examples have r(0) = 0)
    assert r.get(0, 0) == 0
    s = [Fraction(0)] * (order + 1)
    s[0] = Fraction(1)
    term = [Fraction(0)] * (order + 1)
    term[0] = Fraction(1)
    for n in range(1, order + 1):
        new = [Fraction(0)] * (order + 1)
        for a, ca in enumerate(term):
            if ca == 0:
                continue
            for d, cr in r.items():
                if a + d <= order:
                    new[a + d] += ca * (-cr) / n
        term = new
        s = [u + v for u, v in zip(s, term)]
    return s


class Lam:
    def __init__(self, l, r, order=80):
        self.l = l
        self.e = exp_series(r, order)

    def mono(self, k, j):
        n = (self.l + 1) * (j + 1) - 1 - k
        if n < 0:
            return Fraction(0)
        return factorial(j) * self.e[n]

    def __call__(self, P):
        return sum(c * self.mono(k, j) for (k, j), c in P.items())


def potential(l, r):
    W = {(l + 1, 1): Fraction(1)}
    for d, c in r.items():
        W[(d, 0)] = Fraction(c)
    return W


def main():
    one = {(0, 0): Fraction(1)}
    cases = [("W1 = x^3 y + x^2 + x", 2, {1: Fraction(1), 2: Fraction(1)}),
             ("W2 = x^4 y + x - 2 x^2 + 5 x^3", 3, {1: Fraction(1), 2: Fraction(-2), 3: Fraction(5)}),
             ("W3 = x^2 y + x", 1, {1: Fraction(1)})]
    for name, l, r in cases:
        W = potential(l, r)
        lam = Lam(l, r)
        tot = bad = 0
        for k in range(14):
            for j in range(5):
                h = {(k, j): Fraction(1)}
                for gen in (image(h, {}, W), image({}, h, W)):
                    tot += 1
                    bad += (lam(gen) != 0)
        print("E1  %s: lambda_r kills the generators: %d checked, %d failures" % (name, tot, bad))
        mom = []
        Wm = one
        xl = {(l, 0): Fraction(1)}
        momx = []
        for m in range(0, 8):
            mom.append(lam(Wm))
            momx.append(lam(mul(xl, Wm)))
            Wm = mul(Wm, W)
        print("E3  %s: lambda(W^m), m = 0..7: %s ; lambda(x^%d W^n), n = 0..7: %s"
              % (name, [str(v) for v in mom], l, [str(v) for v in momx]))
        if mom[0] != 0:
            delta = mom[1] / mom[0]
            pred2 = delta * (delta + 1) * mom[0]
            print("    an Euler field of constant divergence would need delta = lambda(W)/lambda(1) = %s and lambda(W^2) = %s; actual lambda(W^2) = %s -> %s"
                  % (delta, pred2, mom[2], "consistent" if pred2 == mom[2] else "IMPOSSIBLE: no such field"))

    name, l, r = cases[0]
    W = potential(l, r)
    lam = Lam(l, r)
    tot = bad = 0
    for k in range(0, 5):
        for j in range(0, 3):
            h = {(k, j): Fraction(1)}
            T = add(h, scal(one, lam(h) / lam(one)), -1)
            sol = solve(T, W, max(k, 3) + 3, j + 2) if T else ({}, {})
            tot += 1
            if sol is None or add(image(sol[0], sol[1], W), T, -1):
                bad += 1
    print("E2  %s: x^k y^j - (lambda(x^k y^j)/lambda(1)) in M_W with an explicit field over Q, k <= 4, j <= 2: %d checked, %d failures" % (name, tot, bad))

    x2 = {(2, 0): Fraction(1)}
    for label, T in [("W1^2", power(W, 2)), ("W1^3", power(W, 3)), ("W1^4", power(W, 4)),
                     ("W1 - 1", add(W, one, -1)),
                     ("x^2 W1 - 1! x^2", add(mul(x2, W), x2, -1)),
                     ("x^2 W1^2 - 2! x^2", add(mul(x2, power(W, 2)), scal(x2, 2), -1)),
                     ("x^2 W1^3 - 3! x^2", add(mul(x2, power(W, 3)), scal(x2, 6), -1)),
                     ("x^2 W1^4 - 4! x^2", add(mul(x2, power(W, 4)), scal(x2, 24), -1))]:
        bx, by = bidegree(T)
        sol = solve(T, W, bx + 1, by + 1)
        if sol is None:
            print("E4  %-20s: no solution in the box" % label)
            continue
        p, q = sol
        chk = add(image(p, q, W), T, -1)
        print("E4  %-20s: explicit v over Q, bidegrees p %s, q %s; div(v) - v(W) - target = %s" % (label, bidegree(p), bidegree(q), show(chk)))
    for label, T in [("1", one), ("x^2", x2), ("x^2 W1^2", mul(x2, power(W, 2)))]:
        bx, by = bidegree(T)
        sol = solve(T, W, bx + 3, by + 3)
        print("E4  control %-12s: %s ; lambda = %s" % (label, "no solution in the box" if sol is None else "SOLUTION FOUND (unexpected)", lam(T)))


if __name__ == "__main__":
    main()
