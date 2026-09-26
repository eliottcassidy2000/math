"""Exact arithmetic controls for the crossing/quadric lane; standard library only."""
from collections import defaultdict
from fractions import Fraction as Q
from math import gcd, isqrt
from itertools import product


def require(test, label):
    if not test:
        raise RuntimeError(label)


def prime(n):
    return n >= 2 and all(n % d for d in range(2, isqrt(n) + 1))


def prime_brackets():
    exceptions = []
    pairs = []
    for m in range(1, 101):
        low, high = (2*m-1)**2, (2*m+1)**2
        for n in range(low+1, high+1):
            if 2*n <= high:
                exceptions.append(n)
                if prime(n):
                    pairs.extend((n, j) for j in range(2, high//n+1))
    require(exceptions == [2, 3, 4, 10, 11, 12], 'bracket integer exception list')
    require([n for n in exceptions if prime(n)] == [2, 3, 11], 'prime exception list')
    print('Bracket exceptions', exceptions, 'prime/multiplier pairs', pairs)
    units = [r for r in range(30) if gcd(r, 30) == 1]
    require(units == [1, 7, 11, 13, 17, 19, 23, 29], '30 wheel')
    require(min(n for n in range(2, 100) if gcd(n, 30) == 1 and not prime(n)) == 49,
            'smallest nontrivial wheel false positive')
    require({(r*r-2) % 30 for r in units} == {17, 29}, 'quadratic wheel collapse')
    require(all((r*r-2) % 30 == r for r in (17, 29)), 'two wheel fixed residues')
    print('Wheel units', units, 'Chebyshev image residues', [17, 29], 'hostile11^2-2=119')


def quadratic_controls():
    for c in (0, -1, -2):
        row = []
        for initial in (-1, 0, 1):
            values = [initial]
            for _ in range(4):
                values.append(values[-1]**2+c)
            row.append(values)
        print('Quadratic seed orbits c=', c, row)
    cases = 0
    for denominator in range(2, 18):
        for numerator in range(-2*denominator, 2*denominator+1):
            if gcd(numerator, denominator) != 1:
                continue
            value = Q(numerator, denominator)
            for j in range(1, 5):
                value = value*value-2
                require(value.denominator == denominator**(2**j), 'denominator squares without cancellation')
                cases += 1
    for z in (Q(1, 2), Q(2, 3), Q(3), Q(-2), Q(-3, 4)):
        J = lambda x: x + 1/x
        require(J(z*z) == J(z)**2-2, 'Joukowski semiconjugacy')
    print('Quadratic rational denominator checks', cases)


def rational_chart_controls():
    cases = 0
    for a, b, c, d in product(range(-2, 3), repeat=4):
        if a*d-b*c == 0:
            continue
        # Exclude candidates with a pole at any probe source or its image.
        if any(c*n+d == 0 for n in range(1, 21)):
            continue
        R = lambda n: Q(a*n+b, c*n+d)
        require(any(R(n//2) != R(n)**2-2 for n in range(2,21,2)),
                'no tested nonconstant Mobius chart semiconjugates even steps')
        cases += 1
    for constant in (Q(2), Q(-1)):
        require(constant == constant*constant-2, 'constant positive control')
    def Fiterate(x, depth):
        for _ in range(depth):
            x = x*x-2
        return x
    # A finite unrolled chain works: the obstruction requires a closed finite chart system.
    for j in range(5):
        for x in (Q(1,3), Q(2), Q(5,2)):
            left = Fiterate(2**(j+1)*(x/2), j+1)
            right = Fiterate(2**j*x, j)**2-2
            require(left == right, 'unrolled growing-degree chart control')
    print('Nonconstant Mobius chart hostile census', cases, 'unrolled positive controls15')


MATRICES = {
    'U': ((1,-2,2),(2,-1,2),(2,-2,3)),
    'A': ((1,2,2),(2,1,2),(2,2,3)),
    'D': ((-1,2,2),(-2,1,2),(-2,2,3)),
}
SIGNS = (1, 1, -1)
INVERSES = {key: tuple(tuple(SIGNS[i]*M[j][i]*SIGNS[j] for j in range(3))
                            for i in range(3)) for key, M in MATRICES.items()}


def matvec(M, v):
    return tuple(sum(M[i][j]*v[j] for j in range(3)) for i in range(3))


def triples(cap):
    answer = set()
    for m in range(2, isqrt(cap)+1):
        for n in range(1, m):
            c = m*m+n*n
            if c <= cap and (m-n) % 2 and gcd(m, n) == 1:
                answer.add((m*m-n*n, 2*m*n, c))
    return answer


def double(P):
    a, b, c = P
    return abs(a*a-b*b), 2*a*b, c*c


def inverse_double(P):
    a, b, c = P
    q = isqrt(c)
    require(q*q == c, 'inverse requires square hypotenuse')
    m, n = isqrt((c+a)//2), isqrt((c-a)//2)
    require(m*m == (c+a)//2 and n*n == (c-a)//2 and 2*m*n == b,
            'Euclid inverse')
    odd, even = (m, n) if m % 2 else (n, m)
    return odd, even, q


def address(P):
    reverse = []
    while P != (3, 4, 5):
        candidates = [(key, matvec(M, P)) for key, M in INVERSES.items()]
        candidates = [(key, v) for key, v in candidates if min(v) > 0 and v[2] < P[2]]
        require(len(candidates) == 1, 'unique strictly descending Berggren parent')
        key, P = candidates[0]
        reverse.append(key)
    return ''.join(reversed(reverse))


def pythagorean_controls():
    source = triples(5000)
    for P in source:
        a, b, c = double(P)
        require(a*a+b*b == c*c and gcd(a, b) == 1 and a % 2 and b % 2 == 0,
                'doubling preserves primitive ordered triple')
        require(inverse_double((a,b,c)) == P, 'exact doubling inverse')
        address(P)
    small = triples(500)
    # Independent exhaustive Euclid-parameter universe for target C<=500^2.
    target = {P for P in triples(500**2) if isqrt(P[2])**2 == P[2]}
    require({double(P) for P in small} == target, 'onto every square-hypotenuse target')
    require(double((3,4,5)) == (7,24,25) and address((7,24,25)) == 'UU', 'root double')
    require(double((5,12,13)) == (119,120,169) and address((119,120,169)) == 'AA',
            'doubling fails ancestry preservation')
    print('Primitive triples through5000', len(source), 'square-hypotenuse bijection targets', len(target))
    print('Ancestry hostile: U -> AA, while empty -> UU')
    grouped = defaultdict(int)
    for P in source:
        grouped[P[2]] += 1
    print('Equal-hypotenuse fibres at65,85,1105:', [(c, grouped[c]) for c in (65,85,1105)])


def mobius(n):
    value, p = 1, 2
    while p*p <= n:
        if n % p == 0:
            n //= p
            value = -value
            if n % p == 0:
                return 0
            while n % p == 0:
                n //= p
        p += 1
    return -value if n > 1 else value


def r2_divisors(n):
    return 4*sum(1 if d % 4 == 1 else -1 if d % 4 == 3 else 0
                 for d in range(1, n+1) if n % d == 0)


def theta_controls():
    cap = 2000
    direct = [0]*(cap+1)
    for a in range(-isqrt(cap), isqrt(cap)+1):
        for b in range(-isqrt(cap), isqrt(cap)+1):
            if a*a+b*b <= cap:
                direct[a*a+b*b] += 1
    for n in range(1, cap+1):
        require(direct[n] == r2_divisors(n), 'theta coefficient/divisor convolution')
    grouped = defaultdict(int)
    for P in triples(150):
        grouped[P[2]] += 1
    for c in range(2, 151):
        primitive = sum(mobius(d)*r2_divisors((c//d)**2) for d in range(1,c+1) if c%d == 0)
        require(primitive == 8*grouped[c], 'primitive theta projection')
    require(prime(13) and prime(43) and 13 % 30 == 43 % 30
            and direct[13] == 8 and direct[43] == 0, 'wheel loses sum-of-two-squares prime type')
    one = [0]*301
    one[0] = 1
    for a in range(1, isqrt(300)+1):
        one[a*a] = 2
    three = [sum(one[j]*direct[n-j] for j in range(n+1)) for n in range(301)]
    require(one[2] == 0 and direct[3] == 0 and three[7] == 0, 'one/two/three square support gaps')
    for n in range(1, 301):
        four = sum(direct[j]*direct[n-j] for j in range(n+1))
        jacobi = (8*sum(d for d in range(1,n+1) if n%d == 0) if n%2 else
                  24*sum(d for d in range(1,n+1,2) if n%d == 0))
        require(four == jacobi and four > 0, 'four-square full-support coefficient')
    print('Theta coefficients exact through', cap, 'primitive projection hypotenuses2..150')
    print('Modulo30 hostile: r2(13)=8, r2(43)=0')
    print('Primitive fibre multiplicity for r split primes: 2^(r+2); positivity does not bound v2')
    print('Theta power support hostiles1->2,2->3,3->7; four-square formula exact through300')


if __name__ == '__main__':
    prime_brackets()
    quadratic_controls()
    rational_chart_controls()
    pythagorean_controls()
    theta_controls()
    print('PASS exact finite controls; all-depth statements use the written proofs; no Collatz closure')
