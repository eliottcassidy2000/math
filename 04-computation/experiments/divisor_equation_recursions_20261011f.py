"""Exact divisor-record recursion and a different reciprocal-equation curve.

No minimum positive solution, full elliptic rank, or plane chromatic conclusion.
Only this file's always-active finite controls are counted.
"""
from dataclasses import dataclass
from fractions import Fraction as F
from itertools import product
from math import gcd, isqrt, lcm, prod

CHECKS = 0
SEED = (F(4), F(4))
POSITIVE = (301361533449900458837600, 49105016933436320224063,
            316629033253501281102807)
AREA_PRIMES = (3, 1237, 10193, 11329, 165713)


def need(ok, message):
    if not ok:
        raise ValueError(message)


def check(ok, message):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise RuntimeError(message)


def natural(n, lower=0):
    need(type(n) is int and n >= lower, 'exact integer in range')


def rational(x):
    need(type(x) in (int, F), 'exact rational')
    return F(x)


def prime(n):
    natural(n, 1)
    return n > 1 and all(n % p for p in range(2, isqrt(n)+1))


def first_primes(count):
    natural(count)
    values, n = [], 2
    while len(values) < count:
        if prime(n):
            values.append(n)
        n += 1
    return tuple(values)


def factor(n):
    natural(n, 1)
    result, p = [], 2
    while p*p <= n:
        exponent = 0
        while n % p == 0:
            n //= p
            exponent += 1
        if exponent:
            result.append((p, exponent))
        p += 1
    if n > 1:
        result.append((n, 1))
    return tuple(result)


def tau(n):
    return prod(a+1 for _, a in factor(n))


def divisor_sieve(limit):
    natural(limit, 1)
    counts = [0]*(limit+1)
    for d in range(1, limit+1):
        for n in range(d, limit+1, d):
            counts[n] += 1
    return tuple(counts)


def record_candidates(limit):
    """All prime-prefix, nonincreasing exponent vectors below an exact cap."""
    natural(limit, 1)
    primes = first_primes(limit.bit_length())
    result = []

    def visit(i, last, n, count, exponents):
        result.append((n, count, exponents))
        power = n
        for a in range(1, last+1):
            power *= primes[i]
            if power > limit:
                break
            visit(i+1, a, power, count*(a+1), exponents+(a,))
    visit(0, limit.bit_length(), 1, 1, ())
    return tuple(sorted(result))


def records(limit):
    result, best = [], 0
    for n, count, exponents in record_candidates(limit):
        if count > best:
            result.append((n, count, exponents))
            best = count
    return tuple(result)


def cubic(a, b, c):
    a, b, c = map(rational, (a, b, c))
    s = a+b
    return s**3+c*s*s-c*c*s-c**3-a*b*(2*s+c)


def equation(a, b, c):
    a, b, c = map(rational, (a, b, c))
    need((a+b)*(a+c)*(b+c) != 0, 'all three original denominators nonzero')
    return a/(b+c)+b/(a+c)-c/(a+b)


def audit_point(point):
    if point is None:
        return None
    need(type(point) is tuple and len(point) == 2, 'finite point pair or infinity')
    x, y = map(rational, point)
    need(y*y == x*(x*x+5*x-32), 'point on the exact new cubic')
    return x, y


def add(left, right):
    left, right = audit_point(left), audit_point(right)
    if left is None:
        return right
    if right is None:
        return left
    x, y = left
    u, v = right
    if x == u and y == -v:
        return None
    slope = (3*x*x+10*x-32)/(2*y) if left == right else (v-y)/(u-x)
    z = slope*slope-5-x-u
    answer = (z, slope*(x-z)-y)
    return audit_point(answer)


def multiply(point, m):
    point = audit_point(point)
    natural(m)
    result = None
    while m:
        if m & 1:
            result = add(result, point)
        point = add(point, point)
        m >>= 1
    return result


def to_curve(triple):
    need(type(triple) is tuple and len(triple) == 3, 'exact triple')
    a, b, c = map(rational, triple)
    need(equation(a, b, c) == 0, 'original equation, not just a cleared cubic')
    denominator = 2*(a+b)+c
    need(denominator != 0, 'finite curve chart')
    return audit_point((4*(a+b+2*c)/denominator, 12*(a-b)/denominator))


def from_curve(point, positive=False):
    point = audit_point(point)
    need(type(positive) is bool and point is not None, 'finite point and exact positivity flag')
    x, y = point
    need(x not in (F(-4), F(0), F(8)), 'original rational denominators survive')
    values = (8-x+y, 8-x-y, 4*x-8)
    scale = lcm(*(v.denominator for v in values))
    values = tuple(int(v*scale) for v in values)
    common = gcd(gcd(abs(values[0]), abs(values[1])), abs(values[2]))
    values = tuple(v//common for v in values)
    if values[0]+values[1] < 0:
        values = tuple(-v for v in values)
    need(equation(*values) == 0 and to_curve(values) == point, 'exact inverse chart')
    if positive:
        need(all(v > 0 for v in values), 'positive solution required')
    return values


def valuation_two(q):
    q = rational(q)
    need(q != 0, 'nonzero valuation argument')
    def order(n):
        n = abs(n)
        return (n & -n).bit_length()-1
    return order(q.numerator)-order(q.denominator)


@dataclass(frozen=True)
class Rotation:
    real: F
    imaginary_coefficient: F
    radicand: int
    source: tuple


def rotation(triple, radicand, square_factor):
    """Compile triangle angle into r+q sqrt(-D), with a supplied checked factorization."""
    to_curve(triple)
    need(all(type(v) is int and v > 0 for v in triple), 'positive integer source triangle')
    natural(radicand, 1)
    natural(square_factor, 1)
    a, b, c = triple
    s = a+b
    ratio = F(2*s+3*c, 2*s+c)
    need(ratio.numerator*ratio.denominator == radicand*square_factor**2,
         'exact square-class witness')
    real = F(a*a+b*b-c*c, 2*a*b)
    q = F(s*s-c*c, 2*a*b)*F(square_factor, ratio.denominator)
    need(q > 0 and real*real+radicand*q*q == 1, 'actual unit vector')
    return Rotation(real, q, radicand, triple)


def rejects(call):
    try:
        call()
    except ValueError:
        return True
    return False


def main():
    global CHECKS
    CHECKS = 0
    limit = 100000
    sieve = divisor_sieve(limit)
    expected, best = [], 0
    for n in range(1, limit+1):
        if sieve[n] > best:
            expected.append((n, sieve[n]))
            best = sieve[n]
    rows = records(limit)
    check(tuple((n, d) for n, d, _ in rows) == tuple(expected), 'complete bounded record recursion')
    for n, count, profile in record_candidates(limit):
        check(count == sieve[n] == tau(n), 'independent divisor counts')
        check(all(a >= b for a, b in zip(profile, profile[1:])), 'ordered exponents')
    first = (1, 2, 4, 6, 12, 24, 36, 48, 60, 120, 180, 240, 360)
    second = (1, 2, 3, 4, 6, 8, 9, 10, 12, 16, 18, 20, 24)
    check(tuple(n for n, _, _ in rows[:13]) == first, 'arguments including omitted360')
    check(tuple(d for _, d, _ in rows[:13]) == second, 'paired record values')
    check(tau(72) == tau(60) and (72, tau(72)) not in expected, 'normal shape not sufficient')
    check(36 % 24 != 0 and 60 % 48 != 0, 'record chain is not divisibility recursion')
    divisor_vertices = 0
    for n, count, profile in rows[:13]:
        box = tuple(product(*(range(a+1) for a in profile)))
        check(len(box) == count, 'divisor exponent box')
        for vector in box:
            for i, a in enumerate(profile):
                if vector[i] < a:
                    other = vector[:i]+(vector[i]+1,)+vector[i+1:]
                    check(sum(vector) % 2 != sum(other) % 2, 'all Hasse edges bipartite')
        if sum(profile) >= 3:
            path = [1]
            for p, a in factor(n):
                for _ in range(a):
                    path.append(path[-1]*p)
            check(all(path[j] % path[i] == 0 for i in range(4) for j in range(i+1, 4)),
                  'full comparability has K4')
        divisor_vertices += len(box)
    algebra_cases = 0
    for a, b, c in product(range(-4, 5), repeat=3):
        den = 2*(a+b)+c
        if den:
            x, y = F(4*(a+b+2*c), den), F(12*(a-b), den)
            check(y*y-x**3-5*x*x+32*x == F(576)*cubic(a,b,c)/den**3, 'polynomial curve identity')
            algebra_cases += 1
        if (a+b)*(a+c)*(b+c):
            check(equation(a,b,c) == cubic(a,b,c)/((a+b)*(a+c)*(b+c)), 'denominator-sensitive cubic')
    current = None
    positives = []
    for m in range(1, 25):
        current = add(current, SEED)
        check(current == multiply(SEED, m), 'two addition schedules')
        triple = from_curve(current)
        check(to_curve(triple) == current and equation(*triple) == 0, 'exact signed solution')
        positive = all(v > 0 for v in triple)
        check(positive == (2 < current[0] < 4), 'exact positivity arc')
        check(from_curve((current[0], -current[1])) == (triple[1], triple[0], triple[2]), 'source swap is curve negation')
        if positive:
            positives.append(m)
            a,b,c = triple
            s = a+b
            check(max(a,b) < c < s, 'actual triangle')
            area16 = (s*s-c*c)*(c*c-(a-b)**2)
            check(F(area16) == (s*s-c*c)**2*F(2*s+3*c, 2*s+c), 'Heron square-class identity')
    check(from_curve(multiply(SEED,8), positive=True) == POSITIVE, 'literal positive8P')
    check(from_curve(multiply(SEED,2)) == (32,-25,-17), 'algebra preserved while positivity fails')
    # Infinite signed subgroup mechanism: doubling a negative binary x-order loses2 each time.
    point = multiply(SEED,3)
    for j in range(6):
        check(valuation_two(point[0]) == -4-2*j, 'unbounded denominator doubling')
        check(2*valuation_two(point[1]) == 3*valuation_two(point[0]), 'dominant cubic valuation')
        point = add(point,point)
    # Bounded height control is not a proof of minimum positive solution.
    hits = []
    for s in range(2, 1000):
        for c in range(s//2+1, s):
            numerator, denominator = (s+c)**2*(s-c), 2*s+c
            if numerator % denominator:
                continue
            ab = numerator//denominator
            disc = s*s-4*ab
            if disc >= 0 and isqrt(disc)**2 == disc and (s-isqrt(disc)) % 2 == 0:
                a = (s-isqrt(disc))//2
                if a:
                    hits.append((a,s-a,c))
    check(not hits, 'only the declared sum<=999 search is empty')
    D = prod(AREA_PRIMES)
    check(all(prime(p) for p in AREA_PRIMES) and len(set(AREA_PRIMES)) == len(AREA_PRIMES),
          'squarefree area radicand')
    rot = rotation(POSITIVE, D, 69739671)
    check(D == 71013535489456671 and D not in (1,3,11,33), 'outside the old biquadratic field')
    check(D % 8 == 7, 'new binary square class')
    check(pow((-3) % 11, 5, 11) == pow((-D) % 11, 5, 11) == 10,
          'retained11-adic obstruction marker')
    # An unrelated rational unit-circle parametrization is not this angle.
    a,b,_ = POSITIVE
    rr, qq = F(a*a-b*b,a*a+b*b), F(2*a*b,a*a+b*b)
    check(rr*rr+qq*qq == 1 and rr != rot.real, 'rational Cayley shadow loses triangle angle')
    for bad in (
        lambda: records(True), lambda: records(0), lambda: tau(1.0),
        lambda: to_curve((1,1,-1)), lambda: to_curve((2,3,4)),
        lambda: from_curve((F(0),F(0))), lambda: from_curve((F(8),F(24))),
        lambda: from_curve((F(-4),F(12))), lambda: from_curve(None),
        lambda: from_curve(multiply(SEED,2),positive=True),
        lambda: add((True,4),SEED), lambda: multiply(SEED,True),
        lambda: rotation(POSITIVE,D,69739672), lambda: from_curve(SEED,positive=1),
    ):
        check(rejects(bad), 'malformed or unguarded source rejected')
    print('divisor_equation_recursions_20261011f')
    print('records through100000:',tuple((n,d) for n,d,_ in rows))
    print('first13arguments:',first)
    print('first13values:',second)
    print('declared divisor vertices:',divisor_vertices,'; rational polynomial controls:',algebra_cases)
    print('curve: y^2=x(x^2+5x-32); seed(4,4); positive multiples1..24:',tuple(positives))
    print('8P positive primitive triple:',POSITIVE)
    print('area squarefree D:',D,'; factors:',AREA_PRIMES)
    print('unit rotation real:',rot.real)
    print('unit rotation sqrt(-D) coefficient:',rot.imaginary_coefficient)
    print('field escape is not screen escape: retained p11 Legendre markers are both-1')
    print('no positive solution with A+B<=999; no global leastness or rank claim')
    print('no six-colour equality or plane chromatic conclusion')
    print('exact checks:',CHECKS)


if __name__ == '__main__':
    main()
