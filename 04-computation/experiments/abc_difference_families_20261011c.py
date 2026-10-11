"""Exact difference families, conditional-ABC premises, and geometric cap.

The checker never assumes or numerically verifies ABC. It verifies primitive
input packets, elementary identities, and finite declared positive/hostile cases.
"""

from collections import deque
from dataclasses import dataclass, replace
from fractions import Fraction
from itertools import combinations
from math import gcd, lcm, prod

import plane_six_reframes_parity_20261011b as parity


def integer(value, lower=0):
    if type(value) is not int or value < lower:
        raise ValueError('exact integer above the stated lower bound required')


def prime(n):
    integer(n)
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    p = 3
    while p * p <= n:
        if n % p == 0:
            return False
        p += 2
    return True


def factors(n):
    integer(n, 1)
    result, p = [], 2
    while p * p <= n:
        if n % p == 0:
            exponent = 0
            while n % p == 0:
                n //= p
                exponent += 1
            result.append((p, exponent))
        p = 3 if p == 2 else p + 2
    if n > 1:
        result.append((n, 1))
    return tuple(result)


def radical(n):
    return prod(p for p, _ in factors(n))


@dataclass(frozen=True)
class ABCPacket:
    a: int
    b: int
    c: int


def audit_abc(packet):
    if type(packet) is not ABCPacket:
        raise ValueError('exact additive packet required')
    for x in (packet.a, packet.b, packet.c):
        integer(x, 1)
    if packet.a + packet.b != packet.c or gcd(packet.a, packet.b) != 1:
        raise ValueError('primitive positive a+b=c required')
    return packet


def packet_radical(packet):
    audit_abc(packet)
    return prod({p for x in (packet.a, packet.b, packet.c) for p, _ in factors(x)})


def centre_packets(m):
    integer(m, 2)
    if m % 2:
        raise ValueError('even centre required for the unscaled gap-two packet')
    return audit_abc(ABCPacket(m - 1, 2, m + 1)), audit_abc(ABCPacket(m*m - 1, 1, m*m))


def fibonacci_lucas(n):
    integer(n)
    f, following = 0, 1
    for _ in range(n):
        f, following = following, f + following
    return f, 2 * following - f


def neighbour_factors(n):
    """Exact ordered factorizations of F_n-2 and F_n+2 for odd n>=9."""
    integer(n, 9)
    if n % 2 == 0:
        raise ValueError('odd index required')
    a, b = (n + 3)//2, (n - 3)//2
    fa, la = fibonacci_lucas(a)
    fb, lb = fibonacci_lucas(b)
    plus_sign_pair, minus_sign_pair = (fa, lb), (fb, la)
    minus, plus = (minus_sign_pair, plus_sign_pair) if b % 2 == 0 else (plus_sign_pair, minus_sign_pair)
    f, _ = fibonacci_lucas(n)
    if min(*minus, *plus) <= 1 or prod(minus) != f - 2 or prod(plus) != f + 2:
        raise ValueError('internal product identity failed')
    return minus, plus


def pell_packet(n):
    integer(n, 1)
    f, ell = fibonacci_lucas(n)
    g = gcd(f, ell)
    small, large = sorted((ell*ell, 5*f*f))
    packet = ABCPacket(small//(g*g), 4//(g*g), large//(g*g))
    return audit_abc(packet)


def sieve_modulus(primes):
    if type(primes) is not tuple or any(type(p) is not int or p <= 2 or not prime(p) for p in primes):
        raise ValueError('tuple of odd primes required')
    if primes != tuple(sorted(set(primes))):
        raise ValueError('distinct sorted sieve primes required')
    return 2 * prod(primes)


def wheel_survives(m, primes):
    integer(m)
    sieve_modulus(primes)
    return m % 2 == 0 and all(m % p not in (1, p - 1) for p in primes)


def crt_pair(a, m, b, n):
    for x in (a, b):
        integer(x)
    for x in (m, n):
        integer(x, 1)
    if gcd(m, n) != 1:
        raise ValueError('coprime CRT moduli required')
    return (a + m * ((b-a) * pow(m, -1, n) % n)) % (m*n), m*n


def composite_survivor_ray(primes):
    """Every value on the returned ray survives S but both neighbours factor."""
    modulus = sieve_modulus(primes)
    q = (primes[-1] + 2) if primes else 3
    while not prime(q):
        q += 2
    r = q + 2
    while not prime(r):
        r += 2
    residue, period = crt_pair(0, modulus, 1, q)
    residue, period = crt_pair(residue, period, r - 1, r)
    # Strictly exceed the two forced prime factors, including boundary cases.
    minimum = max(q + 2, r)
    first = residue + max(0, (minimum - residue + period - 1)//period) * period
    return first, period, q, r


def unit_direction(m, d):
    integer(m, 1)
    integer(d, 1)
    if d >= m:
        raise ValueError('0<d<m required')
    denominator = m*m + d*d
    return Fraction(m*m - d*d, denominator), Fraction(2*m*d, denominator)


def audit_points(points):
    if type(points) is not tuple:
        raise ValueError('point tuple required')
    if any(type(p) is not tuple or len(p) != 2 or any(type(x) is not Fraction for x in p) for p in points):
        raise ValueError('exact rational point coordinates required')
    if len(set(points)) != len(points):
        raise ValueError('distinct points required')


def unit_graph(points):
    audit_points(points)
    return parity.Graph(tuple(range(len(points))),
                        tuple((i, j) for i, j in combinations(range(len(points)), 2)
                              if sum((points[i][a]-points[j][a])**2 for a in (0, 1)) == 1))


def odd_denominator_bit(x):
    if type(x) is not Fraction or x.denominator % 2 == 0:
        raise ValueError('rational with odd reduced denominator required')
    return x.numerator % 2


def rational_unit_colouring(points, graph=None):
    """Component-relative arithmetic colouring; no bipartite-search oracle."""
    audit_points(points)
    if graph is None:
        graph = unit_graph(points)
    parity.audit_graph(graph)
    if graph.vertices != tuple(range(len(points))):
        raise ValueError('complete original vertex set required')
    adjacency = {i: [] for i in graph.vertices}
    for i, j in graph.edges:
        if sum((points[i][a]-points[j][a])**2 for a in (0, 1)) != 1:
            raise ValueError('claimed unit edge has wrong length')
        adjacency[i].append(j)
        adjacency[j].append(i)
    roots = {}
    for root in graph.vertices:
        if root in roots:
            continue
        roots[root] = root
        queue = deque((root,))
        while queue:
            i = queue.popleft()
            for j in adjacency[i]:
                if j not in roots:
                    roots[j] = root
                    queue.append(j)
    colours = tuple(odd_denominator_bit(points[i][0]-points[roots[i]][0]) ^
                    odd_denominator_bit(points[i][1]-points[roots[i]][1]) for i in graph.vertices)
    if any(colours[i] == colours[j] for i, j in graph.edges):
        raise ValueError('arithmetic colour failed')
    return colours


def main():
    checks = 0
    def check(ok, message):
        nonlocal checks
        checks += 1
        if not ok:
            raise RuntimeError(message)
    def reject(call):
        nonlocal checks
        checks += 1
        try:
            call()
        except ValueError:
            return
        raise RuntimeError('invalid input accepted')

    twin_centres = []
    for m in range(2, 2002, 2):
        linear, square = centre_packets(m)
        check(packet_radical(linear) == radical(2*(m*m-1)), 'linear radical')
        check(packet_radical(square) == radical(m)*radical(m*m-1), 'square radical and gcd')
        if prime(m - 1) and prime(m + 1):
            twin_centres.append(m)
            check(packet_radical(linear) == 2*(m*m-1) > linear.c, 'twins are not linear abc exceptions')
            check(packet_radical(square) == radical(m)*(m*m-1) > square.c, 'twins are not square abc exceptions')
            check(m == 4 or m % 6 == 0, 'twin-centre modulus3 guard')
    check([p for p in range(2, 2000) if prime(p) and prime(p+1)] == [2], 'only consecutive integer primes')
    check([p for p in range(2, 2000) if prime(p) and prime(p+2) and prime(p+4)] == [3], 'unique twin overlap')
    check(packet_radical(centre_packets(26)[1]) == 390 < 676, 'square-family abc-quality hostile')
    check(packet_radical(centre_packets(244)[0]) == 210 < 245, 'linear-family abc-quality hostile')
    check(radical(1*2*3) == radical(1*8*9) == 6, 'radical loses multiplicities and height')

    for n in range(1, 501):
        f, ell = fibonacci_lucas(n)
        check(ell*ell-5*f*f == 4*(-1)**n, 'golden signed norm')
        check(gcd(f, ell) == (2 if n % 3 == 0 else 1), 'primitive Pell scale')
        packet = pell_packet(n)
        check(12*packet.c >= 5*f*ell, 'ABC conditional height premise only')
        if n <= 26:
            check(packet_radical(packet) <= 10*radical(f)*radical(ell), 'ABC radical comparison premise')
        if n % 2 and n >= 9:
            minus, plus = neighbour_factors(n)
            check(prod(minus) == f-2 and min(minus) > 1, 'lower neighbour nontrivial product')
            check(prod(plus) == f+2 and min(plus) > 1, 'upper neighbour nontrivial product')
    for d in range(1, 31):
        for k in range(1, 16):
            check(fibonacci_lucas(d*k)[0] % fibonacci_lucas(d)[0] == 0, 'Fibonacci divisibility')
    small_fibonacci_twins = {fibonacci_lucas(n)[0] for n in range(1, 31)
                             if prime(fibonacci_lucas(n)[0]) and
                             (prime(max(0, fibonacci_lucas(n)[0]-2)) or prime(fibonacci_lucas(n)[0]+2))}
    check(small_fibonacci_twins == {3, 5, 13}, 'finite controls for proved Fibonacci twin classification')
    check(fibonacci_lucas(19)[0] == 37*113, 'prime index does not imply Fibonacci prime')
    check(fibonacci_lucas(23)[1] == 139*461, 'prime index does not imply Lucas prime')

    wheels = []
    for primes in ((), (3,), (3, 5), (3, 5, 7), (3, 5, 7, 11)):
        modulus = sieve_modulus(primes)
        count = sum(wheel_survives(m, primes) for m in range(modulus))
        check(count == prod(p-2 for p in primes), 'complete finite-wheel count')
        first, period, q, r = composite_survivor_ray(primes)
        for t in range(64):
            m = first + period*t
            check(wheel_survives(m, primes), 'composite family survives declared wheel')
            check((m-1) % q == 0 and (m+1) % r == 0 and m-1 > q and m+1 > r,
                  'both neighbours are genuinely composite')
        wheels.append((primes, modulus, count, first, period, q, r))

    geometry_count = edge_count = 0
    for m in range(2, 15):
        for d in range(1, m):
            x, y = unit_direction(m, d)
            check(x*x+y*y == 1, 'difference-of-squares unit parameter')
            denominator = lcm(x.denominator, y.denominator)
            a, b = int(x*denominator), int(y*denominator)
            check(denominator % 2 == 1 and (a+b) % 2 == 1, 'primitive unit-step parity')
            tx, ty = Fraction(m, 2*d), Fraction(d, 4*m)
            points = ((tx, ty), (tx+x, ty+y), (tx-y, ty+x), (tx+x-y, ty+y+x))
            graph = unit_graph(points)
            colours = rational_unit_colouring(points, graph)
            check(len(graph.edges) == 4 and sum(colours) == 2, 'rational rotated square')
            check(not parity.compile_boundary(graph, (0,)*4, graph.vertices).odd_cycle,
                  'independent signed-boundary verification')
            geometry_count += 1
            edge_count += len(graph.edges)
    shifted = ((Fraction(1,2), Fraction(0)), (Fraction(3,2), Fraction(0)))
    check(rational_unit_colouring(shifted) == (0, 1), 'component origin is a necessary sidecar')
    reject(lambda: odd_denominator_bit(shifted[0][0]))
    check(rational_unit_colouring(()) == (), 'empty rational graph')
    reject(lambda: centre_packets(True))
    reject(lambda: centre_packets(5))
    reject(lambda: audit_abc(ABCPacket(True, 2, 3)))
    reject(lambda: audit_abc(ABCPacket(2, 4, 6)))
    reject(lambda: neighbour_factors(11.0))
    reject(lambda: neighbour_factors(10))
    reject(lambda: pell_packet(0))
    reject(lambda: sieve_modulus((3, 9)))
    reject(lambda: sieve_modulus((3, 3)))
    reject(lambda: sieve_modulus((3.0, 5)))
    reject(lambda: unit_direction(5, 5))
    reject(lambda: rational_unit_colouring(((0, Fraction(0)),)))
    reject(lambda: rational_unit_colouring((shifted[0], shifted[0])))
    reject(lambda: rational_unit_colouring(((Fraction(0), Fraction(0)), (Fraction(2), Fraction(0))),
                                         parity.Graph((0,1), ((0,1),))))

    print('ABC difference families: exact arithmetic and a geometric transfer boundary')
    print('Twin-centre controls: even m=2..2000;', len(twin_centres), 'twin centres; unique overlapping pair 3-5-7')
    print('Fibonacci/Lucas identities: n=1..500; explicit composite neighbours for every tested odd n>=9')
    print('Proved Fibonacci twin-prime values: 3,5,13; prime-index hostiles F19=37*113, L23=139*461')
    print('Conditional ABC only: log(rad(Fn*Ln))/log(Fn*Ln) tends to1; no primality conclusion')
    for row in wheels:
        print('Wheel / modulus / survivors / both-composite ray(first,period,q,r):', row)
    print(f'Rational rotated-square controls: {geometry_count}; exact unit edges: {edge_count}; all arithmetic two-colourings pass')
    print('ABC-quality hostiles: 675+1=676 radical390; 243+2=245 radical210')
    print(f'Checks: {checks}')
    print('Scope: no abc proof, no twin-prime infinitude, no plane six/seven witness, no new tournament count')


if __name__ == '__main__':
    main()
