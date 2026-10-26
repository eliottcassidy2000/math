#!/usr/bin/env python3
"""Exact finite wheel, self-hit, height and recurrence checks. Stdlib only."""
from fractions import Fraction
from math import gcd, isqrt, prod
import json

checks = 0


def require(test, label):
    global checks
    checks += 1
    if not test:
        raise RuntimeError(label)


def trial_lpf(n):
    for d in range(2, isqrt(n) + 1):
        if n % d == 0:
            return d
    return n


def T(x):
    return x * (x + 2), 2 * (x + 1), (x + 1) ** 2 + 1


def D(v):
    a, b, c = v
    return -a + 2 * b + 2 * c, -2 * a + b + 2 * c, -2 * a + 2 * b + 3 * c


def floor_count(H, modulus, roots):
    return sum((H - r) // modulus - (2 - r) // modulus for r in roots)


def mobius_terms(primes):
    terms = [(1, 1, (0,))]
    for ell in primes:
        local = (0,) if ell == 2 else (0, ell - 2)
        extra = []
        for d, mu, roots in terms:
            inv = pow(d, -1, ell)
            lifts = tuple(sorted(r + d * ((a - r) * inv % ell)
                                 for r in roots for a in local))
            require(len(lifts) == len(set(lifts)) == len(roots) * len(local),
                    'CRT forbidden-root cardinality')
            require(all(r * (r + 2) % (d * ell) == 0 for r in lifts),
                    'all generated forbidden roots are literal roots')
            extra.append((d * ell, -mu, lifts))
        terms += extra
    return terms


def main():
    HMAX = 10000
    lpf = list(range(HMAX + 3))
    for p in range(2, isqrt(HMAX + 2) + 1):
        if lpf[p] == p:
            for n in range(p * p, HMAX + 3, p):
                if lpf[n] == n:
                    lpf[n] = p
    for n in range(2, HMAX + 3):
        require(lpf[n] == trial_lpf(n), 'least-factor sieve vs independent trial division')
    primes = [n for n in range(2, HMAX + 3) if lpf[n] == n]
    cuts = [p for p in primes if p <= 101]
    twins = {x for x in range(3, HMAX + 1) if lpf[x] == x and lpf[x + 2] == x + 2}
    pi2 = [0] * (HMAX + 1)
    for x in range(1, HMAX + 1):
        pi2[x] = pi2[x - 1] + (x in twins)

    # Strict counts by direct gcd, independently compared to least factors.
    strict = {}
    moduli = {}
    for B in cuts:
        M = prod(p for p in primes if p <= B)
        moduli[B] = M
        values = [0] * (HMAX + 1)
        for x in range(3, HMAX + 1):
            good = gcd(x * (x + 2), M) == 1
            require(good == (lpf[x] > B and lpf[x + 2] > B), 'gcd and factor definitions agree')
            values[x] = values[x - 1] + good
        strict[B] = values

    # Exact residue lifts; the complete largest period is 30030.
    M, wheel = 2, [1]
    wheel_rows = []
    test_heights = list(range(3, 201)) + [500, 999, 1000, 5000, 10000]
    for B in [2, 3, 5, 7, 11, 13]:
        if B > 2:
            new = []
            for a in wheel:
                forbidden = {(-a * pow(M, -1, B)) % B,
                             (-(a + 2) * pow(M, -1, B)) % B}
                children = [a + j * M for j in range(B) if j not in forbidden]
                require(len(forbidden) == 2 and len(children) == B - 2, 'two forbidden lifts')
                require(all((v * (v + 2)) % B for v in children), 'all child residues accepted')
                require(all((v - 1) // 2 - (a - 1) // 2 == j * (M // 2)
                            for j in range(B) if j not in forbidden
                            for v in [a + j * M]), 'exact D-depth lift, not a single child')
                new += children
            wheel, M = sorted(new), M * B
        pset = [p for p in primes if p <= B]
        require(wheel == [a for a in range(M) if gcd(a * (a + 2), M) == 1],
                'nested lifts equal full literal wheel')
        require(len(wheel) == prod(p - (1 if p == 2 else 2) for p in pset), 'wheel product count')
        require(M - 1 in wheel, 'coherent negative-one class')
        terms = mobius_terms(pset)
        require(len(terms) == 2 ** len(pset), 'full squarefree divisor universe')
        error_bound = 2 * 3 ** (len(pset) - 1)
        require(sum(len(roots) for d, mu, roots in terms) == error_bound, 'Möbius root error count')
        for H in test_heights:
            literal = strict[B][H]
            wc = floor_count(H, M, wheel)
            mc = sum(mu * floor_count(H, d, roots) for d, mu, roots in terms)
            require(literal == wc == mc, 'gcd, allowed floors, and Möbius floors agree')
            require(abs(Fraction(literal) - Fraction((H - 2) * len(wheel), M))
                    <= min(len(wheel), error_bound), 'exact stated error bound')
        wheel_rows.append({'B': B, 'modulus': M, 'residues': len(wheel),
                           'density': str(Fraction(len(wheel), M)),
                           'mobius_divisor_terms': len(terms), 'root_error_bound': error_bound})

    # Proper-prime-factor exception rule, independently evaluated per prime.
    previous = set(range(3, 201))
    for B in range(2, 102):
        pset = [p for p in primes if p <= B]
        M = prod(pset)
        accepted = {x for x in range(3, 201)
                    if all((x % p or x == p) and ((x + 2) % p or x + 2 == p) for p in pset)}
        rough = {x for x in range(3, 201) if gcd(x * (x + 2), M) == 1}
        patch = {x for x in twins if x <= min(200, B)}
        require(not rough & patch and accepted == rough | patch, 'exact disjoint self-hit patch')
        require(accepted <= previous, 'restored sieve nested')
        previous = accepted
    require(previous == twins & set(range(3, 201)), 'restored sieve stabilizes at finite height')
    require(strict[101][101] == 0 and pi2[101] == 9, 'strict limit loses even genuine twins')

    # Check the complete finite-height recurrence at every H, not a fitted sample.
    for H in range(3, HMAX + 1):
        B = isqrt(H + 2)
        cut = max(p for p in cuts if p <= B)
        require(B < H, 'recurrence decreases height')
        require(pi2[H] == strict[cut][H] + pi2[B], 'square-root recurrence at every height')
    for B in cuts:
        q = next(p for p in primes if p > B)
        for H in range(3, min(HMAX, q * q - 3) + 1):
            require(pi2[H] == strict[B][H] + pi2[min(H, B)], 'next-unsieved-prime height threshold')

    first_false = []
    first_nonsquare = None
    for B in [p for p in cuts if p <= 61]:
        q = next(p for p in primes if p > B)
        bad = [x for x in range(3, 5001) if lpf[x] > B and lpf[x + 2] > B and x not in twins]
        require(bool(bad), 'bounded first-false search has witness')
        x = bad[0]
        require(x >= q * q - 2, 'universal first-false lower bound')
        boundary = q * q - 2
        boundary_survives = lpf[boundary] > B and lpf[boundary + 2] > B
        require(boundary_survives == (lpf[boundary] == boundary), 'boundary survivor iff lower is prime')
        upper_root = isqrt(x + 2)
        square_pattern = (lpf[x] == x and upper_root * upper_root == x + 2
                          and lpf[upper_root] == upper_root)
        if not square_pattern and first_nonsquare is None:
            first_nonsquare = B
        first_false.append({'B': B, 'q': q, 'first_possible': boundary, 'first_false': [x, x + 2],
                            'least_factors': [lpf[x], lpf[x + 2]], 'prime_upper_square_pattern': square_pattern})
    require(first_nonsquare == 61, 'first tested-stage failure of square-only pattern')
    require(first_false[-1]['first_false'] == [4757, 4759], 'exact nonsquare hostile')
    require(4757 == 67 * 71 and trial_lpf(4759) == 4759 and 4487 == 7 * 641,
            'independent hostile factors and prime companion')

    require(T(-1) == (-1, 0, 1), 'coherent negative lift is a degenerate null point')
    for x in range(3, 1000, 2):
        a, b, c = T(x)
        require(a * a + b * b == c * c and gcd(a, b) == 1, 'positive primitive D ray')
        require(D(T(x)) == T(x + 2), 'actual integer branch action')
    require(T(1) == (3, 4, 5), 'excluded unit parameter is nondegenerate root')
    examples = []
    for H in [3, 5, 7, 100, 1000, 10000]:
        B = isqrt(H + 2)
        cut = max(p for p in cuts if p <= B)
        examples.append({'H': H, 'B': B, 'strict': strict[cut][H], 'self_hit_patch': pi2[B], 'twins': pi2[H]})
    print(json.dumps({'status': 'PASS', 'checks': checks,
                      'scope': 'FINITE-EXACT; all-parameter results use README.md proofs',
                      'universe': {'integer_height_max': HMAX, 'least_factor_checked_through': HMAX + 2,
                                   'strict_cutoffs': cuts, 'wheel_cutoffs': [2, 3, 5, 7, 11, 13],
                                   'self_hit_B': [2, 101], 'self_hit_x': [3, 200],
                                   'first_false_prefix_max': 5000, 'first_false_prime_cutoff_max': 61},
                      'wheels': wheel_rows, 'recurrence_examples': examples,
                      'first_false_census': first_false,
                      'first_square_only_failure_cutoff': first_nonsquare,
                      'coherent_degenerate_lift': {'x': -1, 'center': 0, 'triple': T(-1)}},
                     indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
