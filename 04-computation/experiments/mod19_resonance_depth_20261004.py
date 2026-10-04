"""Exact shared ternary depth of prime observers of a six-exponent ray.

Elementary order and CRT proofs are in the companion note. Prime and order
certificates are exact; there is no probable-prime or Collatz orbit oracle.
"""
from collections import Counter
from functools import lru_cache
from math import gcd, isqrt, lcm
import json


def need(ok, message):
    if not ok:
        raise ArithmeticError(message)


@lru_cache(None)
def prime(p):
    if type(p) is not int or p < 2:
        return False
    return all(p % d for d in range(2, isqrt(p)+1))


def factor(n):
    need(type(n) is int and n >= 1, 'positive integer factorization')
    result, d = {}, 2
    while d*d <= n:
        while n % d == 0:
            result[d] = result.get(d, 0)+1
            n //= d
        d += 1
    if n > 1:
        result[n] = result.get(n, 0)+1
    return result


def valuation(n, p):
    need(n != 0, 'nonzero valuation')
    value = 0
    while n % p == 0:
        value += 1
        n //= p
    return value


def certified_order(base, modulus, candidate):
    need(gcd(base, modulus) == 1 and pow(base, candidate, modulus) == 1,
         'candidate modular return')
    divisors = factor(candidate)
    order = candidate
    for q in divisors:
        while order % q == 0 and pow(base, order//q, modulus) == 1:
            order //= q
    need(pow(base, order, modulus) == 1 and all(pow(base, order//q, modulus) != 1
                                              for q in factor(order)), 'exact order certificate')
    return order


@lru_cache(None)
def prime_order(base, p, depth=1):
    need(prime(p) and p not in (2, 3), 'prime observer away from two and three')
    need(type(depth) is int and depth >= 1, 'positive prime-adic precision')
    # Factor a modest p-1 once; the new factor at depth>1 is already the known p.
    candidate, modulus = (p-1)*p**(depth-1), p**depth
    factors = set(factor(p-1)) | ({p} if depth > 1 else set())
    need(pow(base, candidate, modulus) == 1, 'Euler return at the supplied prime')
    order = candidate
    for q in sorted(factors):
        while order % q == 0 and pow(base, order//q, modulus) == 1:
            order //= q
    need(pow(base, order, modulus) == 1 and all(pow(base, order//q, modulus) != 1
                                              for q in factors if order % q == 0), 'all prime shortening exclusions')
    return order


def resonance_depth(p):
    return max(0, valuation(prime_order(2, p), 3)-1)


def compatibility_index(p, ternary_depth, prime_depth, hub_valuation=0):
    need(type(ternary_depth) is int and ternary_depth >= 1, 'positive ternary precision')
    need(type(prime_depth) is int and prime_depth >= 1, 'positive auxiliary precision')
    need(type(hub_valuation) is int and hub_valuation >= 0, 'nonnegative hub valuation')
    need(prime(p) and p not in (2, 3), 'admissible auxiliary prime')
    if prime_depth <= hub_valuation:
        return 1
    return 3**min(ternary_depth-1, resonance_depth(p))


def root_row2_address(block, modulus):
    """(2^(4+6b)-1)/3 modulo modulus, retaining division precision."""
    numerator = (pow(2, 4+6*block, 3*modulus)-1) % (3*modulus)
    need(numerator % 3 == 0, 'exact root-row division')
    return numerator//3


def main():
    primes = [p for p in range(5, 1001) if prime(p)]
    auxiliary = [3511, 5779, 87211]
    universe = primes+auxiliary
    histogram, least, lift_checks = Counter(), {}, 0
    for p in universe:
        d = prime_order(2, p)
        block_order = d//gcd(d, 6)
        depth = resonance_depth(p)
        need(prime_order(64, p) == block_order, 'power-of-element order')
        if p <= 1000:
            histogram[depth] += 1
            least.setdefault(depth, p)
        for k in range(1, 4):
            actual = prime_order(64, p, k)
            need(actual % block_order == 0, 'order projects to base prime')
            quotient = actual//block_order
            while quotient % p == 0:
                quotient //= p
            need(quotient == 1 and valuation(actual, 3) == depth, 'only auxiliary-prime order growth')
            for a in range(1, 7):
                for v in range(4):
                    period = 1 if k <= v else prime_order(64, p, k-v)
                    need(gcd(3**(a-1), period) == compatibility_index(p, a, k, v), 'all retained-depth gcd controls')
                    lift_checks += 1
    need(min(p for p in universe if resonance_depth(p) > 0) == 19, 'smallest ternary-resonant auxiliary prime')
    need((127-1) % 9 == 0 and prime_order(2, 127) == 7 and resonance_depth(127) == 0,
         'unit-group size does not supply a missing base-two order')

    # Exact earlier cyclotomic certificates, independently rerun without imports.
    layers = [(1, {}), (2, {19:1}), (3, {87211:1}),
              (4, {163:1, 135433:1, 272010961:1})]
    layer_report = []
    for r, factors in layers:
        t = 2**(3**(r-1))
        numerator = t*t-t+1
        need(numerator % 3 == 0 and numerator % 9 != 0, 'exactly one exceptional factor three')
        quotient = numerator//3
        recovered = 1
        for p, exponent in factors.items():
            need(prime(p), 'exact trial-division certificate for layer prime')
            recovered *= p**exponent
            need(certified_order(2, p, 2*3**r) == 2*3**r, 'primitive layer order')
        need(recovered == quotient, 'complete inherited layer factorization')
        layer_report.append(dict(layer=r, quotient=quotient, factors=factors,
                                 resonance=None if r == 1 else r-1))

    selected = []
    for p in (19, 73, 109, 127, 163, 3511, 5779, 87211):
        d = prime_order(2, p)
        selected.append(dict(prime=p, order2=d, order2_factors=factor(d),
                             order64=prime_order(64,p), resonance=resonance_depth(p),
                             order64_depths=[prime_order(64,p,k) for k in range(1,4)]))
    need(prime_order(2,3511,2) == prime_order(2,3511), 'a genuine flat first order lift')
    need(prime_order(2,3511,3) == 3511*prime_order(2,3511), 'subsequent lift resumes')
    need(resonance_depth(5779) == resonance_depth(87211) == 2,
         'equal retained ternary depth despite unequal prime clocks')
    need(prime_order(64,5779) == 107*prime_order(64,87211) and pow(2,54,5779) == 2944,
         'resonance quotient does not identify whole phases')

    # Complete joint-address languages on one genuinely completed root ray.
    joint = []
    for p in (5,7,19,73,109,127,163,5779,87211):
        n = prime_order(64,p)
        for a in range(1,5):
            m, common = 3**(a-1), gcd(3**(a-1), n)
            pairs = {(root_row2_address(b,3**a), root_row2_address(b,p))
                     for b in range(lcm(m,n))}
            ternary = {root_row2_address(b,3**a):b for b in range(m)}
            auxiliary_addresses = {root_row2_address(b,p):b for b in range(n)}
            need(len(ternary) == m and len(auxiliary_addresses) == n, 'two individual address bijections')
            predicted = {(x,y) for x,bx in ternary.items() for y,by in auxiliary_addresses.items()
                         if (bx-by) % common == 0}
            need(pairs == predicted and len(pairs)*common == m*n, 'complete joint language and exact compatibility index')
            joint.append(dict(prime=p, ternary_depth=a, candidate_pairs=m*n,
                              compatible_pairs=len(pairs), index=common))
    print(json.dumps(dict(status='PROVED shared-depth and elementary all-level existence; FINITE-EXACT declared controls',
        prime_universe=dict(all_primes_from5_through1000=len(primes), extra_primes=auxiliary,
                            histogram=dict(sorted(histogram.items())), least_prime_per_exact_depth=least),
        lifted_gcd_checks=lift_checks, named_clocks=selected, inherited_binary_layers=layer_report,
        complete_joint_languages=joint,
        boundary='Shared ternary phase is a quotient of a fixed inverse-ray block clock; retain the hub, row, source certificate and original integer identity.'), indent=2))
    print('PASS: exact primality/order certificates and all checks active under -O')


if __name__ == '__main__':
    main()
