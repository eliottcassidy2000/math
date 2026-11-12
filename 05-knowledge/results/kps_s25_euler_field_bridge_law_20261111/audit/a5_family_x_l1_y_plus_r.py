"""a5_family_x_l1_y_plus_r.py - audit upgrade U1/U2: the residue form of Theorem A2 settles the whole family

    W = x^(l+1) y + r(x),   l >= 1,  r a polynomial with r(0) = 0      (no Euler field is needed).

Claims checked here in exact arithmetic for several (l, r), including r'(0) = 0 and deg r > l:
  F1  lambda_r(x^k y^j) = j! [x^((l+1)(j+1)-1-k)] exp(-r(x)) kills both families of generators (k < 14, j < 5)
  F2  lambda_r(W^m) = 0 for every m >= l (m <= 9 tested) and lambda_r(x^l W^n) = n! (n <= 9)
  F3  one-dimensionality on a box: every monomial x^k y^j (k <= l + 2, j <= 2) is congruent modulo M_W to a
      multiple of x^l, with an explicit vector field over Q:  x^k y^j - lambda_r(x^k y^j) x^l  lies in M_W
      (lambda_r(x^l) = 1).  The general statement (codimension one for every l and r) is proved by hand in AUDIT.md.
  F4  explicit certificates over Q:  W^l in M_W and x^l W^2 - 2 x^l in M_W.
  F5  the closed formula lambda_r(x^k W^m) = [x^(l-k)] ( m! exp(-r(x)) e_m(r(x)) ), e_m = truncated exponential
      (k <= l + 2, m <= 6); since m! e^(-s) e_m(s) = m! + O(s^(m+1)) this gives F2 for all m.
Hence M_W is not a Mathieu-Zhao subspace for any such W:  P = W^l,  h = x^l.
"""
from fractions import Fraction
from math import factorial

from a2_membership_over_Q import add, mul, scal, power, image, solve, show, bidegree
from a3_no_euler_field_example import Lam, potential


def main():
    one = {(0, 0): Fraction(1)}
    cases = [("x^2 y + x + x^2 + x^3", 1, {1: 1, 2: 1, 3: 1}),
             ("x^3 y + x^2", 2, {2: 1}),
             ("x^3 y + x + x^2  [= x (x^2 y + x + 1)]", 2, {1: 1, 2: 1}),
             ("x^4 y + 2 x - x^4", 3, {1: 2, 4: -1}),
             ("x^3 y - 3 x + x^5", 2, {1: -3, 5: 1})]
    for name, l, r in cases:
        r = {d: Fraction(c) for d, c in r.items()}
        W = potential(l, r)
        lam = Lam(l, r)
        xl = {(l, 0): Fraction(1)}
        bad = tot = 0
        for k in range(14):
            for j in range(5):
                h = {(k, j): Fraction(1)}
                for gen in (image(h, {}, W), image({}, h, W)):
                    tot += 1
                    bad += (lam(gen) != 0)
        mom = []
        momx = []
        Wm = one
        for m in range(0, 10):
            mom.append(lam(Wm))
            momx.append(lam(mul(xl, Wm)))
            Wm = mul(Wm, W)
        f2 = all(v == 0 for v in mom[l:]) and all(momx[n] == factorial(n) for n in range(10))
        print("W = %s (l = %d)" % (name, l))
        print("   F1 generators killed: %d checked, %d failures; lambda(1) = %s, lambda(x^l) = %s" % (tot, bad, lam(one), lam(xl)))
        print("   F2 lambda(W^m), m = 0..9: %s ; zero for m >= l and lambda(x^l W^n) = n!: %s" % ([str(v) for v in mom], f2))
        b3 = t3 = 0
        for k in range(0, l + 3):
            for j in range(0, 3):
                h = {(k, j): Fraction(1)}
                T = add(h, scal(xl, lam(h)), -1)
                t3 += 1
                if not T:
                    continue
                bx, by = bidegree(T)
                sol = solve(T, W, max(bx, l) + l + 2, by + 2)
                if sol is None or add(image(sol[0], sol[1], W), T, -1):
                    b3 += 1
        print("   F3 x^k y^j - lambda(x^k y^j) x^l in M_W with an explicit field, k <= l + 2, j <= 2: %d checked, %d failures" % (t3, b3))
        # F5: lambda_r(x^k W^m) = [x^(l-k)] ( m! exp(-r) e_m(r) ),  e_m the truncated exponential
        b5 = t5 = 0
        order = 40
        rser = [Fraction(0)] * (order + 1)
        for d, c in r.items():
            rser[d] = c

        def smul(a, b):
            out = [Fraction(0)] * (order + 1)
            for i, ai in enumerate(a):
                if ai:
                    for jj, bj in enumerate(b):
                        if bj and i + jj <= order:
                            out[i + jj] += ai * bj
            return out

        for m in range(0, 7):
            em = [Fraction(0)] * (order + 1)
            term = [Fraction(0)] * (order + 1)
            term[0] = Fraction(1)
            for i in range(0, m + 1):
                if i > 0:
                    term = [t / i for t in smul(term, rser)]
                em = [u + v for u, v in zip(em, term)]
            series = [factorial(m) * c for c in smul(lam.e[:order + 1], em)]
            Wm = power(W, m)
            for k in range(0, l + 3):
                lhs = lam(mul({(k, 0): Fraction(1)}, Wm))
                rhs = series[l - k] if l - k >= 0 else Fraction(0)
                t5 += 1
                b5 += (lhs != rhs)
        print("   F5 lambda(x^k W^m) = [x^(l-k)] m! e^(-r) e_m(r), k <= l + 2, m <= 6: %d checked, %d failures" % (t5, b5))
        for label, T in (("W^l", power(W, l)), ("x^l W^2 - 2 x^l", add(mul(xl, power(W, 2)), scal(xl, 2), -1))):
            bx, by = bidegree(T)
            sol = solve(T, W, bx + 2, by + 1)
            ok = sol is not None and not add(image(sol[0], sol[1], W), T, -1)
            print("   F4 %-18s in M_W, explicit vector field over Q: %s" % (label, ok))


if __name__ == "__main__":
    main()
