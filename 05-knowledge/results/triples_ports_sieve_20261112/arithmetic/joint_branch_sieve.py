"""Exact joint sieve for three formal scalar counts, not graph-return closure.

The all-height proof and scope are in JOINT_BRANCH_SIEVE.md. Every check
uses RuntimeError, including under python -O. Standard library only.
"""
from fractions import Fraction
from math import gcd, isqrt, lcm, prod
from pathlib import Path
import json

from wheel_refinement import need, prime, prime_factors, order2, exact_density

HERE = Path(__file__).resolve().parent
COUNTS = [(99, 23), (297, 71), (297, 73)]
EXPECTED = {
    23: (11, [5, 6]),
    31: (5, [1]),
    47: (23, [22]),
    53: (13, [1]),
    59: (29, [11]),
    71: (35, [18]),
    79: (13, [3]),
    83: (41, [13, 16]),
}


def literal(k, primes):
    return all((a * pow(2, k, p) + b) % p
               for p in primes for a, b in COUNTS)


def compressed(k):
    if k % 12 != 6:
        return False
    t = (k - 6) // 12
    return all(t % m not in residues for m, residues in EXPECTED.values())


def aliquot(n):
    if n == 1:
        return 0
    total = 1
    for d in range(2, isqrt(n) + 1):
        if n % d == 0:
            total += d
            if d * d != n:
                total += n // d
    return total


def main():
    primes = [p for p in range(2, 100) if prime(p)]
    need(len(primes) == 25, 'exact first 25 primes')
    # These three primes already force the single phase k=6 mod12.
    forcing = [r for r in range(12) if literal(r, [5, 7, 13])]
    need(forcing == [6], 'complete forcing phase word')
    need(all((a * pow(2, 1, 2) + b) % 2 for a, b in COUNTS),
         'all three counts odd for positive exponents')
    original = []
    reduced = []
    for p in primes[1:]:
        d = order2(p)
        branch_residues = [[r for r in range(d)
                            if (a * pow(2, r, p) + b) % p == 0]
                           for a, b in COUNTS]
        original.append(dict(prime=p, order=d, forbidden=branch_residues))
        m = d // gcd(d, 12)
        forbidden = [t for t in range(m) if not literal(6 + 12 * t, [p])]
        if forbidden:
            need(p in EXPECTED, 'no undeclared conditional prime')
            need((m, forbidden) == EXPECTED[p], 'complete conditional phase word')
            reduced.append(dict(prime=p, modulus=m, forbidden=forbidden))
        else:
            need(p not in EXPECTED, 'no omitted conditional prime')
    need(len(reduced) == len(EXPECTED), 'all conditional constraints')

    clocks = [35, 11, 23, 13, 29, 41]
    allowed = [27, 9, 22, 11, 28, 39]
    need(all(gcd(a, b) == 1 for i, a in enumerate(clocks)
             for b in clocks[i + 1:]), 'pairwise coprime CRT clocks')
    t_period = prod(clocks)
    period = 12 * t_period
    density_t = prod(Fraction(a, m) for a, m in zip(allowed, clocks))
    density = density_t / 12
    clauses = [dict(modulus=m, residue=r)
               for m, residues in EXPECTED.values() for r in residues]
    independent_density, intersections = exact_density(clauses)
    need(density_t == independent_density, 'factorization vs generalized CRT')
    need(t_period == 136871735 and period == 1642460820, 'period values')
    need(density == Fraction(5346, 136735), 'joint density')
    count = density * period
    need(count.denominator == 1 and count.numerator == 64216152,
         'exact accepted phase count')
    minimality = []
    for q in prime_factors(period):
        shift = period // q
        k = next(k for k in range(1, 100000)
                 if compressed(k) != compressed(k + shift))
        # Test the actual 25-prime predicates, not only the compressed formula.
        need(literal(k, primes) == compressed(k), 'literal minimality left')
        need(literal(k + shift, primes) == compressed(k + shift),
             'literal minimality right')
        minimality.append(dict(prime_factor=q, shift=shift, k=k,
                               at_k=compressed(k), at_shift=compressed(k + shift)))
    for k in range(1, 10001):
        need(literal(k, primes) == compressed(k), 'literal all-prime comparison')

    seed_k = 6
    seed_counts = [a * 2 ** seed_k + b for a, b in COUNTS]
    need(seed_counts == [6359, 19079, 19081], 'seed formal counts')
    need(all(prime(q) for q in seed_counts), 'complete trial division seed primality')
    aliquot_seed = []
    for q in seed_counts:
        need(aliquot(q * q) == q + 1, 'literal proper-divisor sum of prime-square seed')
        aliquot_seed.append(dict(target=q + 1, maximum_source=q * q,
                                 certificate='inherited sharp bound plus prime-square equality'))
    for source in range(2, 5001):
        target = aliquot(source)
        if target < 2:
            continue
        need(source <= (target - 1) ** 2, 'literal sharp aliquot source bound')
        root = isqrt(source)
        need((source == (target - 1) ** 2) == (root * root == source and prime(root)),
             'literal prime-square equality boundary')
    orders = [order2(q) for q in seed_counts]
    composite_step = lcm(period, *orders)
    need(orders == [3179, 9539, 3180], 'seed prime orders')
    need(composite_step == 239978082932247660, 'joint composite progression')
    for t in range(1, 101):
        k = seed_k + composite_step * t
        need(compressed(k) and literal(k, primes), 'joint sieve progression control')
        need(all((a * pow(2, k, q) + b) % q == 0
                 for (a, b), q in zip(COUNTS, seed_counts)),
             'each formal branch has its own fixed prime divisor')

    report = dict(
        status='PASS',
        scope='All positive exponents k, by complete local phase words and CRT. '
              'Three formal scalar counts; no simultaneous tournament-return realization claim.',
        formal_counts=COUNTS, first_25_primes=primes,
        forcing_primes=[5, 7, 13], forcing_modulus=12, forcing_residues=forcing,
        original_prime_clauses=original, conditional_t_clauses=reduced,
        independent_clocks=clocks, allowed_per_clock=allowed,
        t_minimal_period=t_period, k_minimal_period=period,
        accepted_k_phases=count.numerator,
        t_density=[density_t.numerator, density_t.denominator],
        k_density=[density.numerator, density.denominator],
        compatible_CRT_intersections=intersections, minimality_witnesses=minimality,
        direct_exponent_controls=10000,
        all_prime_seed=dict(k=seed_k, counts=seed_counts, primality='complete trial division'),
        aliquot_bridge=dict(seed_maxima=aliquot_seed, bound_source_controls=4999,
                            inherited=True, enumeration_extension=False,
                            simultaneous_composite_consequence='At all three targets count+1, '
                            'A135244 is strictly below count squared along the displayed u>=1 progression.'),
        simultaneous_composites=dict(k0=seed_k, step=composite_step,
                                     fixed_prime_divisors=seed_counts, orders=orders,
                                     t_scope='every integer t>=1', finite_controls=100))
    (HERE / 'joint_branch_sieve.json').write_text(
        json.dumps(report, indent=2, sort_keys=True) + '\n', encoding='utf-8')
    print(json.dumps(dict(status='PASS', minimal_period=period,
                          density=report['k_density'], accepted_phases=count.numerator,
                          prime_seed=seed_counts, composite_step=composite_step)))


if __name__ == '__main__':
    main()
