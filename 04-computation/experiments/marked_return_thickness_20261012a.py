"""Exact return-thickness controls; no general Collatz, JC, or HN conclusion."""
from fractions import Fraction
from math import gcd
import sympy as s

CHECKS = 0


def check(ok, label):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise RuntimeError(label)


def zero(expr, label):
    check(s.cancel(expr) == 0, label)


def partitions(n, k):
    if n == 0:
        yield ()
        return
    for prefix in partitions(n-1, k):
        for colour in range(min(k, max(prefix, default=-1)+2)):
            yield prefix+(colour,)


def run():
    x, sigma, eta = s.symbols('x sigma eta')
    c = -(sigma**2 + sigma + 2)

    def cubic(t):
        return x**3-t*x**2-(t*t+2*t+3)*x+t**3+2*t*t+3*t+1

    p = cubic(sigma)
    other = cubic(-sigma-1)
    f = x*x+c
    f2 = s.expand(f*f+c)
    f3 = s.expand(f2*f2+c)
    zero(f3-x-(f-x)*p*other, 'universal third-return factorization')
    zero(s.discriminant(p, x)-(4*sigma**2+6*sigma+9)**2,
         'cycle discriminant')
    zero(s.resultant(p, f-x, x)-(4*sigma**2+6*sigma+9),
         'fixed-point collision resultant')

    # This equality is only modulo p. Its first-order defect is retained.
    ell = -x-s.Rational(1, 2)
    c_rev = c-sigma-s.Rational(3, 4)
    defect = ell**2+c_rev-(-f2-s.Rational(1, 2))
    zero(defect-(x+sigma)*p, 'affine reversal defect')

    mu_sigma = -8*(sigma**3+2*sigma**2+3*sigma+1)
    zero(s.rem(s.diff(f3, x)-mu_sigma, p, x), 'multiplier constant on cycle')
    mu = s.expand(mu_sigma.subs(sigma, (eta-1)/2))
    mu_rev = s.expand(mu.subs(eta, -eta-1))
    q = eta**2+eta+7
    zero(mu-(1-7*eta-eta**2-eta**3), 'multiplier formula')
    zero(mu-1+eta*q, 'parabolic factor')
    zero(mu+mu_rev-(q+2), 'reversal multiplier sum')
    zero(mu-mu_rev+(2*eta+1)*q, 'reversal multiplier difference')
    zero(mu*mu_rev-(-q**3+7*q**2+q+1), 'reversal multiplier product')

    # Actual distinct cycle at each listed trace. Modulus powers retain jets.
    for sig, expected_c, thickness, expected_mu in [
        (s.Rational(-1, 2), s.Rational(-7, 4), 2, s.Integer(1)),
        (s.Rational(-3, 4), s.Rational(-29, 16), 1, s.Rational(35, 8)),
        (s.Integer(-1), s.Integer(-2), 1, s.Integer(8)),
        (s.Integer(0), s.Integer(-2), 1, s.Integer(-8)),
    ]:
        pp = p.subs(sigma, sig)
        cc = c.subs(sigma, sig)
        check(cc == expected_c, 'parameter specialization')
        check(mu_sigma.subs(sigma, sig) == expected_mu, 'multiplier specialization')
        check(s.gcd(pp, x*x+cc-x) == 1, 'not a fixed-point divisor')
        check(s.discriminant(pp, x) != 0, 'distinct cycle points')
        modulus = pp**(thickness+1)
        residue = x
        first_defect = None
        for m in range(1, 17):
            for _ in range(3):
                residue = s.rem(residue*residue+cc, modulus, x)
            defect_m = s.expand(residue-x)
            check(s.rem(defect_m, pp**thickness, x) == 0, 'return at claimed thickness')
            check(defect_m != 0, 'no return at next thickness')
            if first_defect is None:
                first_defect = defect_m
            if thickness == 2:
                zero(defect_m-m*first_defect, 'parabolic defect accumulates linearly')
            else:
                zero(s.rem(s.diff(defect_m, x), pp, x)-(expected_mu**m-1),
                     'first-jet multiplier persists')
        print(f'cycle c={cc}, trace={sig}, multiplier={expected_mu}: '
              f'exact thickness {thickness}, repetitions 1..16 checked')

    # Deliberately separate the critical orbit from the displayed 3-cycle.
    for cc in (Fraction(-7, 4), Fraction(-29, 16)):
        value = Fraction(0)
        seen = set()
        terms = []
        for j in range(1, 9):
            value = value*value+cc
            check(value.denominator == cc.denominator**(2**(j-1)),
                  'exact nonintegral critical denominator')
            check(value not in seen, 'critical orbit has not closed')
            seen.add(value)
            if j <= 3:
                terms.append(str(value))
        print(f'critical orbit c={cc}: {terms}; denominator law checked through 8')
    a, b = -7, 4
    check(a*(a+b)**2+b**3 == 1, '63/64 unit cancellation')
    check(Fraction(a, b)*(Fraction(a, b)+1)**2 == Fraction(-63, 64),
          'critical-orbit bridge to 63')
    check(2**6-1 == 3**2*7, 'Mersenne factorization')
    check(pow(2, 6, 9) == 1 and all(pow(2, j, 9) != 1 for j in range(1, 6)),
          'new order at old prime 3')
    old_product = 1
    stripped = None
    for n in range(1, 7):
        term = 2**n-1
        if n == 6:
            stripped = term
            while gcd(stripped, old_product) > 1:
                stripped //= gcd(stripped, old_product)
        old_product *= term
    check(stripped == 1, 'no new prime at term6, without discarding valuation')
    print('63: no new prime, but ord_9(2)=6; critical-orbit denominator still grows')

    # The next residual of the finite charged grammar is a different question
    # from whether this individual integer has a supplied ROOT certificate.
    source, value = 127, 127
    prefix = (1, 1, 1, 1, 1, 1, 2, 4, 3, 3)
    for a in prefix:
        check(2*(value+1) > source+1, 'no earlier band exit at127')
        z = 3*value+1
        check((z & -z).bit_length()-1 == a, 'exact next valuation at127')
        value = z >> a
    check(source == 2*63+1 and value == 29, 'Mersenne residual endpoint')
    check(2*(value+1) <= source+1, 'first band exit at127')
    check(sum(a % 2 for a in prefix) == 8, 'signed-charge-free actual charge')
    print('grammar residual127=2*63+1: first band word(1^6,2,4,3,3) reaches29, charge8')

    # Complete boundary universes, not claimed realizable geometric gadgets.
    bells = (1, 1, 2, 5, 15, 52, 203)
    for n in range(7):
        rows = tuple(partitions(n, 6))
        check(len(rows) == bells[n], 'small complete partition universe')
        check(len(rows) == len(set(rows)), 'no duplicated boundary state')
        if n >= 2:
            equal = tuple(p for p in rows if p[0] == p[1])
            different = tuple(p for p in rows if p[0] != p[1])
            check(len(equal) == bells[n-1], 'equal ports contract to one port')
            check(len(different) == bells[n]-bells[n-1], 'separating support budget')
    check(tuple(partitions(2, 6)) == ((0, 0), (0, 1)), 'two-port orbit dichotomy')
    check(len(tuple(partitions(7, 6))) == 876, 'palette cutoff excludes seven blocks')
    print('six named ports at palette6: 203 states = 52 equal + 151 separating')
    print('two named ports: only equal/different orbits; no partial separating orbit')
    print(f'PASS: {CHECKS} exact checks; symbolic identities plus finite controls, no universality claim')


if __name__ == '__main__':
    run()
