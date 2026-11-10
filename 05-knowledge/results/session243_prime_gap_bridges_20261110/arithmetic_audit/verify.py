#!/usr/bin/env python3
"""Independent, standard-library arithmetic audit. No producer imports.

The proofs, not these finite controls, establish the all-n implications.
"""
import argparse
from bisect import bisect_left, bisect_right
import json
from math import factorial, gcd, isqrt, prod


def demand(condition, label):
    if not condition:
        raise RuntimeError(label)


def eratosthenes(bound):
    flags = bytearray(b'\1') * (bound + 1)
    flags[:2] = b'\0\0'
    for p in range(2, isqrt(bound) + 1):
        if flags[p]:
            flags[p*p:bound+1:p] = b'\0' * ((bound-p*p)//p + 1)
    return flags


def valuation_factorial(n, p):
    exponent = 0
    while n:
        n //= p
        exponent += n
    return exponent


def unitary_pairs(n, primes):
    """Gray-code walk, retaining both complementary factors at every step."""
    N = factorial(n)
    exponents = {p: valuation_factorial(n, p) for p in primes if p <= n}
    C = N // 4
    A = 1 << (exponents[2] - 2)
    B = C // A
    odd_blocks = [p ** e for p, e in exponents.items() if p != 2]
    demand(prod(odd_blocks) == B, 'Legendre packets incomplete')
    previous = 0
    seen = set()
    hits = []
    for i in range(1 << len(odd_blocks)):
        gray = i ^ (i >> 1)
        if i:
            changed = gray ^ previous
            demand(changed & (changed - 1) == 0, 'Gray transition not one bit')
            block = odd_blocks[changed.bit_length() - 1]
            if gray & changed:
                demand(B % block == 0, 'invalid B to A transfer')
                B //= block
                A *= block
            else:
                demand(A % block == 0, 'invalid A to B transfer')
                A //= block
                B *= block
        demand(A * B == C and gcd(A, B) == 1, 'unitary factor violation')
        demand(A % 2 == 0 and B % 2 == 1, 'canonical even side lost')
        demand(A not in seen, 'duplicate partition')
        seen.add(A)
        if abs(A - B) == 1:
            m = A + B
            demand(m * m == N + 1, 'inverse has wrong quotient')
            demand({(m-1)//2, (m+1)//2} == {A, B}, 'inverse not reversible')
            hits.append((m, min(A, B), max(A, B)))
        previous = gray
    demand(len(seen) == 1 << (len(exponents)-1), 'wrong partition count')
    return sorted(hits), len(seen), exponents


def sign_packet(m, exponents):
    moduli = [1 << (exponents[2]-1)]
    moduli.extend(p ** e for p, e in exponents.items() if p != 2)
    return all((m-1) % q == 0 or (m+1) % q == 0 for q in moduli)


def brown_audit(primes):
    rows = []
    roots = []
    partition_count = 0
    divisor_scan_count = 0
    for n in range(4, 61):
        hits, count, exponents = unitary_pairs(n, primes)
        partition_count += count
        N = factorial(n)
        root = isqrt(N+1)
        expected = [root] if root*root == N+1 else []
        demand([h[0] for h in hits] == expected, 'Brown inverse versus isqrt')
        if n <= 10:
            # An independent all-divisor enumeration, with no prime packets.
            C = N//4
            divisors = [a for a in range(1, isqrt(C)+1)
                        if C % a == 0 and gcd(a, C//a) == 1]
            divisor_scan_count += isqrt(C)
            demand(len(divisors) == count, 'all-divisor partition count')
            direct = sorted((a+C//a, a, C//a) for a in divisors if C//a-a == 1)
            demand(direct == hits, 'all-divisor Brown roots')
        for m, a, b in hits:
            demand(sign_packet(m, exponents), 'true root rejected by sign packet')
            demand(gcd(m, N) == 1, 'Brown root not rough')
            roots.append([n, m, a, b])
        rows.append({'n': n, 'partitions': count, 'roots': expected})
    demand(roots == [[4,5,2,3], [5,11,5,6], [7,71,35,36]], 'finite Brown roots')
    boundary = []
    for n in range(4):
        q = factorial(n)+1
        boundary.append({'n': n, 'factorial_plus_one': q, 'is_square': isqrt(q)**2 == q})
        demand(isqrt(q)**2 != q, 'unexpected n<4 root')
    residues = 0
    residue_rows = []
    for n in range(4, 10):
        exponents = {p: valuation_factorial(n, p) for p in primes if p <= n}
        N = factorial(n)
        good = 0
        # Includes every parity and the complete common period of both tests.
        for m in range(N//2):
            sign = sign_packet(m, exponents)
            demand(sign == ((m*m-1) % N == 0), 'sign iff on complete residues')
            good += sign
            residues += 1
        demand(good == 1 << len(exponents), 'CRT sign-class count')
        residue_rows.append({'n': n, 'period': N//2, 'successful_residues': good})
    signed_controls = 0
    for n in range(4, 61):
        exponents = {p: valuation_factorial(n, p) for p in primes if p <= n}
        for m in range(-128, 129):
            demand(sign_packet(m, exponents) == ((m*m-1) % factorial(n) == 0),
                   'sign iff on negative and nonnegative controls')
            signed_controls += 1
    e4 = {p: valuation_factorial(4, p) for p in primes if p <= 4}
    e5 = {p: valuation_factorial(5, p) for p in primes if p <= 5}
    demand(gcd(7, factorial(5)) == 1 and not sign_packet(7, e5), 'wheel loss control')
    demand(sign_packet(23, e4) and 23*23-1 == 22*factorial(4), 'quotient loss control')
    demand(all((29*29-1) % p == 0 for p in [2,3,5,7]) and
           (29*29-1) % factorial(7) != 0, 'radical versus prime powers')
    demand(all((m*m-1) % 24 == 0 for m in range(24) if gcd(m,24) == 1),
           'n=4 wheel units do satisfy the divisibility equation')
    return {'status': 'PASS', 'n_universe': [4,60], 'partition_count': partition_count,
            'roots': roots, 'n_below_four': boundary, 'all_divisor_scan_n': [4,10],
            'all_divisor_candidates': divisor_scan_count, 'complete_residue_tests': residues,
            'residue_rows': residue_rows, 'signed_integer_tests': signed_controls,
            'hostiles': {'wheel': {'n':5,'m':7}, 'lost_quotient': {'n':4,'m':23,'quotient':22},
                         'lost_exponents': {'n':7,'m':29}}, 'rows': rows}


def gap_audit(primes):
    largest_center = prod(primes[:8])
    flags = eratosthenes(largest_center+1000)
    all_primes = [p for p, flag in enumerate(flags) if flag]
    rows = []
    for n in range(2, 9):
        P = prod(primes[:n])
        lower = all_primes[bisect_left(all_primes, P)-1]
        upper = all_primes[bisect_right(all_primes, P)]
        excluded_lower = all_primes[bisect_left(all_primes, P-1)-1]
        excluded_upper = all_primes[bisect_right(all_primes, P+1)]
        L, F = P-excluded_lower, excluded_upper-P
        l, f = P-lower, upper-P
        lower_bit, upper_bit = bool(flags[P-1]), bool(flags[P+1])
        demand(l == (1 if lower_bit else L) and f == (1 if upper_bit else F),
               'offset-one restoration')
        interior = [q for q in range(excluded_lower+1, excluded_upper) if flags[q]]
        skipped = [q for q in [P-1,P+1] if flags[q]]
        # Stronger than equal-offset case: the property holds for every F,L.
        demand(interior == skipped, 'unexpected interior prime')
        demand(not any(flags[lower+1:upper]), 'crossing pair not consecutive')
        if F == L:
            demand((len(interior) == 0) == (not lower_bit and not upper_bit),
                   'equal-offset consecutiveness iff')
        rows.append({'n':n,'P':P,'F':F,'L':L,'nearest_pair':[lower,upper],
                     'gap':upper-lower,'raw_pair':[excluded_lower,excluded_upper],
                     'interior_primes':interior})
    r6 = rows[4]
    demand(r6 == {'n':6,'P':30030,'F':17,'L':17,'nearest_pair':[30029,30047],
                   'gap':18,'raw_pair':[30013,30047],'interior_primes':[30029]},
           '30030 complete hostile')
    demand(rows[6]['F'] == rows[6]['L'] == 23 and rows[6]['gap'] == 46 and
           not rows[6]['interior_primes'], 'P8 positive consecutive example')
    demand(flags[3] and gcd(6,3) == 3, 'endpoint-equals-divisor exception')
    demand(rows[0]['L'] == 3 and rows[0]['F'] == 5, 'n2 lower exception')
    demand(all(2-d < 2 for d in range(2,20)), 'L1 is undefined')
    demand(30031 == 59*509 and not flags[30031], 'Euclid composite')
    demand(30049 == 151*199 and not flags[30049], 'prime-offset composite')
    obstructions = []
    for d in range(2, 101):
        r = next(p for p in primes if d % p == 0)
        P = 1
        for n, p in enumerate(primes, 1):
            P *= p
            if P % r == 0 and P-d > r:
                break
        demand(P % r == 0 and P-d > r, 'obstruction center not reached')
        # Verify proper factors exactly. This is a certificate even for huge P.
        for c in range(1,6):
            for sign in [-1,1]:
                endpoint = c*P+sign*d
                demand(endpoint % r == 0 and endpoint//r > 1 and r > 1,
                       'fixed multiplier proper factor certificate')
        # Once at this center, every subsequent tested center preserves it.
        for p in primes[n:n+3]:
            P *= p
            demand(P % r == 0 and P-d > r, 'obstruction persistence')
        obstructions.append({'d':d,'prime_divisor':r,'first_obstructing_n':n})
    return {'status':'PASS','sieve_bound':len(flags)-1,'primorial_indices':[2,8],
            'rows':rows,'fixed_offset_universe':[2,100],'fixed_multipliers':[1,5],
            'proper_factor_certificates':99*5*2,'fixed_offset_rows':obstructions,
            'P2_d3_exception': {'lower_endpoint':3,'upper_endpoint':9,'both_prime':False}}


def prime_square_audit(primes):
    flags = eratosthenes(primes[30]**2)
    rows = []
    for p, q in zip(primes[1:30], primes[2:31]):
        d = (q-p)//2
        direct = [x for x in range(p*p+1,q*q) if flags[x]]
        reconstructed = []
        for a in range(p,q,2):
            reconstructed.extend(x for x in range(a*a+1,(a+2)**2) if flags[x])
        demand(direct == reconstructed, 'prime-square decomposition')
        demand(q*q-p*p == 4*d*(p+d), 'prime-square interval length')
        demand(all(not flags[a*a] for a in range(p+2,q,2)), 'shared endpoint prime')
        # Verify the four Oppermann subintervals of each odd bracket align.
        for a in range(p,q,2):
            cuts = [a*a, a*(a+1), (a+1)**2, (a+1)*(a+2), (a+2)**2]
            demand(cuts == sorted(set(cuts)), 'Oppermann cut order')
            demand(all(not flags[x] for x in cuts), 'composite cut point')
        rows.append({'p':p,'q':q,'gap':q-p,'odd_square_brackets':d,'prime_count':len(direct),
                     'Legendre_conditional_minimum':q-p,
                     'Legendre_counting_implication_to_four':d>=2})
    demand([q for q in range(5,9) if flags[q]] == [5,7], 'Brocard p2 exception')
    demand(rows[0]['Legendre_conditional_minimum'] == 2, 'twin Legendre boundary')
    return {'status':'PASS','pair_count':len(rows),'first_pair':[3,5],
            'last_pair':[113,127],'largest_endpoint_prime_index':31,
            'rows':rows,'caution':'Conditional minima are implications, not proofs of either conjecture.'}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--fail', choices=['skip-one','fixed-offset','brown-quotient','sign-exponent'])
    args = parser.parse_args()
    flags = eratosthenes(200)
    primes = [p for p, flag in enumerate(flags) if flag]
    brown = brown_audit(primes)
    gaps = gap_audit(primes)
    squares = prime_square_audit(primes)
    if args.fail == 'skip-one':
        demand(gaps['rows'][4]['gap'] == 34, 'deliberate false consecutive gap')
    if args.fail == 'fixed-offset':
        demand(gaps['P2_d3_exception']['lower_endpoint'] > 3, 'deliberate erased equality exception')
    if args.fail == 'brown-quotient':
        demand(23*23 == factorial(4)+1, 'deliberate sign-to-equality overclaim')
    if args.fail == 'sign-exponent':
        demand((29*29-1) % factorial(7) == 0, 'deliberate radical-to-full-power overclaim')
    print(json.dumps({'status':'PASS','implementation':'independent; Eratosthenes + Legendre valuations + Gray codes',
                      'brown':brown,'prime_gaps':gaps,'prime_square_brocard':squares}, indent=2,sort_keys=True))


if __name__ == '__main__':
    main()
