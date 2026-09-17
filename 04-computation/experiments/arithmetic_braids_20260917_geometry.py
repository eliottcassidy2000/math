"""Exact controls for the 2026-09-17 arithmetic-braids geometry lane.

Reproduce: python 04-computation/experiments/arithmetic_braids_20260917_geometry.py
All checks use integers/Fraction; no numerical orbit or factorization is trusted.
"""
from fractions import Fraction as F
from math import gcd, isqrt


def check(predicate, message):
    if not predicate:
        raise RuntimeError(message)


def trim(p):
    p = list(map(F, p))
    while len(p) > 1 and not p[-1]:
        p.pop()
    return p


def add(p, q):
    out = [F(0)] * max(len(p), len(q))
    for i, v in enumerate(p):
        out[i] += v
    for i, v in enumerate(q):
        out[i] += v
    return trim(out)


def mul(p, q):
    out = [F(0)] * (len(p) + len(q) - 1)
    for i, v in enumerate(p):
        for j, w in enumerate(q):
            out[i+j] += v*w
    return trim(out)


def scale(p, a):
    return trim([a*v for v in p])


def rem(p, q):
    p, q = trim(p), trim(q)
    while len(p) >= len(q) and p != [0]:
        k, a = len(p)-len(q), p[-1]/q[-1]
        for i, v in enumerate(q):
            p[i+k] -= a*v
        p = trim(p)
    return p


def compose(p, q):
    out = [F(0)]
    for v in reversed(p):
        out = add(mul(out, q), [v])
    return out


def marked_cycle(t):
    den = 2*t*(t+1)
    return [(t**3+2*t*t+t+1)/den,
            (t**3-t-1)/den,
            -(t**3+2*t*t+3*t+1)/den]


def main():
    # Independent polynomial-quotient path: no symbolic-algebra dependency.
    X, D, L = [0, 1], [-2, 0, 1], [-F(1, 2), -1]
    q7, q9 = [-1, -2, 1, 1], [1, -3, 0, 1]
    D2, D3 = compose(D, D), compose(D, compose(D, D))
    check(add(D3, scale(X, -1)) == mul(add(D, scale(X, -1)), mul(q7, q9)),
          "Chebyshev period-three factorization")
    for q, c, reduced, multiplier in [
            (q7, -F(7, 4), [1, -1, -1], 1),
            (q9, -F(11, 4), [2, -1, -1], 19)]:
        check(rem(D2, q) == list(map(F, reduced)), "inverse cycle polynomial")
        fL = add(mul(L, L), [c])
        check(rem(add(fL, scale(compose(L, D2), -1)), q) == [0],
              "reversed affine cycle transport")
        # Cubic image under L, obtained by direct composition q(-x-1/2).
        image_cubic = scale(compose(q, L), -1)
        check(-8*image_cubic[0] == multiplier, "cycle multiplier")
    parabolic = [-F(1, 8), -F(9, 4), F(1, 2), 1]
    f = [-F(7, 4), 0, 1]
    check(add(compose(f, compose(f, f)), scale(X, -1)) ==
          mul(add(f, scale(X, -1)), mul(parabolic, parabolic)),
          "parabolic dynatomic square")
    print("polynomial quotient identities: PASS (q7 multiplier 1; q9 multiplier 19)")

    # Rational-cycle chart from THM-4139: all reduced t=p/q in this box.
    cycle_parameters = {F(p, q) for p in range(-12, 13)
                        for q in range(1, 13)} - {F(0), F(-1)}
    cycles = set()
    fixed_parameters = set()
    rho = lambda t: -1/(t+1)
    tau = lambda t: 1/t
    for t in cycle_parameters:
        z = marked_cycle(t)
        c = z[1]-z[0]**2
        check(len(set(z)) == 3, "nondegenerate rational chart")
        check(all(z[(i+1) % 3] == z[i]**2+c for i in range(3)), "chart cycle")
        sigma = sum(z)
        check(c == -(sigma*sigma+sigma+2), "trace parameter")
        cp, sp = c-sigma-F(3, 4), -sigma-F(3, 2)
        reverse = [-z[0]-F(1, 2), -z[2]-F(1, 2), -z[1]-F(1, 2)]
        check(all(reverse[(i+1) % 3] == reverse[i]**2+cp for i in range(3)),
              "universal affine time reversal")
        check(sum(reverse) == sp and cp-sp-F(3, 4) == c, "involution")
        check(t == z[0]+z[1] and t == (z[1]-z[2])/(z[0]-z[1]), "chart recovery")
        check(marked_cycle(rho(t)) == [z[2], z[0], z[1]], "cyclic remarking")
        check(marked_cycle(tau(t)) == [-z[2]-F(1, 2), -z[1]-F(1, 2), -z[0]-F(1, 2)],
              "reciprocity reverses cycle")
        check(rho(rho(rho(t))) == t and tau(tau(t)) == t and
              tau(rho(tau(t))) == rho(rho(t)), "S3 relation")
        eta = 2*sigma+1
        check(eta == (t**3+t*t-2*t-1)/(t*(t+1)), "seventh-root trace quotient")
        check(eta+F(1, 2) == (t-1)*(t+2)*(2*t+1)/(2*t*(t+1)), "fixed locus")
        if cp == c:
            fixed_parameters.add(t)
        cycles.add((c, tuple(sorted(z))))
    check(fixed_parameters == {F(1), F(-2), -F(1, 2)}, "fixed chart orbit")
    print(f"rational chart: {len(cycle_parameters)} markings, {len(cycles)} cycles; reversal PASS")
    print("chart reciprocity/rotation S3: PASS; fixed orbit [-2, -1/2, 1]")

    # Complete Euclid universe, with inherited primitivity/parity filters.
    triple_count = 0
    depth_count = 0
    for m in range(2, 81):
        for n in range(1, m):
            if gcd(m, n) != 1 or (m-n) % 2 != 1:
                continue
            a, b, c = m*m-n*n, 2*m*n, m*m+n*n
            triple_count += 1
            check(c+b == (m+n)**2 and c-b == (m-n)**2, "odd-square half identity")
            short, long = sorted([a, b])
            e, d, ell = F(short*short, c*c), F(long*long, c*c), F(a*b, c*c)
            r, t = e/d, F(short, long)
            check(e+d == 1 and ell*ell == e*d and r == t*t,
                  "unit-diameter semicircle")
            check(ell == t/(1+t*t), "one-dimensional coordinate identity")
            check(F(a*b, c).denominator == c, "primitive altitude denominator")
            initial_b, initial_c, old_a, signed_product = b, c, [], 1
            for k in range(5):
                check(a*a+b*b == c*c and gcd(a, b) == 1, "primitive square orbit")
                check(a % 2 != 0 and b % 2 == 0 and abs(a) > 1, "new prime nontriviality")
                check(all(gcd(a, prev) == 1 for prev in old_a), "pairwise coprime odd legs")
                check(c == initial_c**(2**k), "fixed hypotenuse support")
                check(b == 2**k*initial_b*signed_product, "even leg accumulation")
                ap, bp, cp = a*a-b*b, 2*a*b, c*c
                check(2*F(ap, cp) == (2*F(a, c))**2-2, "Chebyshev semiconjugacy")
                tangent = F(min(abs(a), abs(b)), max(abs(a), abs(b)))
                doubled_tangent = 2*tangent/(1-tangent*tangent)
                check(F(min(abs(ap), abs(bp)), max(abs(ap), abs(bp))) ==
                      min(doubled_tangent, 1/doubled_tangent), "folded angle tent map")
                old_a.append(a)
                signed_product *= a
                a, b, c = ap, bp, cp
                depth_count += 1
    print(f"primitive Euclid universe m<=80: {triple_count} triples; {depth_count} depth checks PASS")

    a, b, c = 3, 4, 5
    orbit = []
    for k in range(4):
        orbit.append((a, b, c))
        a, b, c = a*a-b*b, 2*a*b, c*c
    print(f"signed triple square orbit: {orbit}")

    shapes = {}
    duplicates = []
    max_altitude = F(0)
    max_indices = []
    for k in range(2, 1001):
        a, b, c = k*k-1, 2*k, k*k+1
        check(gcd(gcd(a, b), c) == (1 if k % 2 == 0 else 2), "family primitive content")
        ell = F(a*b, c*c)
        if ell > max_altitude:
            max_altitude, max_indices = ell, [k]
        elif ell == max_altitude:
            max_indices.append(k)
        shape = tuple(sorted([F(a, c), F(b, c)]))
        if shape in shapes:
            duplicates.append((shapes[shape], k))
        shapes[shape] = k
        reciprocal_k = F(k+1, k-1)
        ap, bp, cp = reciprocal_k**2-1, 2*reciprocal_k, reciprocal_k**2+1
        check(tuple(sorted([ap/cp, bp/cp])) == shape, "folded family involution")
    check(max_altitude == F(12, 25) and max_indices == [2, 3], "integer family altitude maximum")
    check(duplicates == [(2, 3)], "unique integer shape duplication")
    print("integer family 2<=k<=1000: maximum 12/25 at k=2,3; sole duplicate (2,3); PASS")

    # Exact third-iterate no-new-prime criterion, without factoring iterates.
    critical = F(0)
    orbit = []
    for k in range(3):
        critical = critical*critical-F(7, 4)
        orbit.append(critical)
    check(orbit == [F(-7, 4), F(21, 16), F(-7, 256)], "critical numerator exception")
    thue_solutions, rational_count = [], 0
    for b in range(2, 2001):
        for a in range(-2*b+1, -b):
            if gcd(a, b) != 1:
                continue
            rational_count += 1
            unit = a**3+2*a*a*b+a*b*b+b**3
            check(gcd(unit, a*(a+b)*b) == 1, "third numerator primitive cofactor")
            if abs(unit) == 1:
                thue_solutions.append((a, b, unit))
    check(thue_solutions == [(-7, 4, 1)], "bounded third-iterate exception census")
    print(f"critical orbit at -7/4: {[str(x) for x in orbit]}")
    print(f"FINITE-EXACT -2<a/b<-1, 2<=b<=2000: {rational_count} parameters; units {thue_solutions}")
    print("All checks passed.")


if __name__ == "__main__":
    main()
